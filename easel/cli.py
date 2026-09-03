#!/usr/bin/env python3
import argparse
import glob
import logging
import os
import shlex
import shutil
import sys
import textwrap
from ruamel.yaml import YAML, comments
from ruamel.yaml.scalarstring import SingleQuotedScalarString
import subprocess
import signal

## import command line function to set up snakemake environment and run snakemake once DEF file is constructed
from . import preprocess as preprocess_input
from .validate import validate, report

## Logging
logger = logging.getLogger("easel")
logging.basicConfig(level=logging.INFO, format="[%(levelname)-8s] %(message)s")

## Absolute paths on file system
## Workflow files (Snakefile_*, rules/, envs/, prof/) sit alongside the
## package directory
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

## Description for CLI interface
DESCRIPTION = '''\
easel - codon alignment and selection screen over thousands of transcripts/genes.

  easel -toga2 /path/to/genomes/hg38/TOGA2 -asm assemblies.txt \\
         -sb selected.bed -it species_tree.nh --twoBit_path hg38.exons.2bit \\
         -ab srv -dr

Input mode is mutually exclusive: -toga2 extracts alignments from TOGA2
annotations, -free takes FASTA files from other sources. A tree is required unless
--do_alignment_only is set or -ct, --compute_tree is set. Otherwise, supply a species tree with -it, 
or point at precomputed ones with --input_gene_trees.

All selection analyses are off by default.
Available tools are: aBSREL, BUSTED, MEME, Relax, Bayescode.
Enable them with -ab/-bu/-me/-re/-bc.
'''

EPILOG = '''\
Nothing is executed unless -dr (dry run) or -rs (run) is given; without either,
easel only writes DEF.yaml.
'''


class SleasyHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Verbatim description/epilog; option help is dedented and re-wrapped.

    Conventions inside a help string:
      * a blank line starts a new paragraph
      * a paragraph whose first line starts with '*' keeps its hand alignment
    """

    def _split_lines(self, text, width):
        lines = []
        for para in textwrap.dedent(text).strip().split("\n\n"):
            if para.lstrip().startswith("*"):
                for bullet in textwrap.dedent(para).strip().splitlines():
                    lines.extend(
                        textwrap.wrap(bullet, width, subsequent_indent="    ") or [""])
            else:
                lines.extend(textwrap.wrap(" ".join(para.split()), width) or [""])
            lines.append("")
        return lines[:-1] if lines else []

    def _format_action_invocation(self, action):
        ## "-toga2 DIR, --toga2_reference_path DIR" -> "-toga2, --toga2_… DIR"
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        return ", ".join(action.option_strings) + " " + self._format_args(action, default)


def _help_formatter(prog):
    columns = shutil.get_terminal_size((100, 24)).columns
    return SleasyHelpFormatter(prog, max_help_position=6,
                               width=max(60, min(columns - 2, 100)))


ALIGNER_CHOICES = ["prank", "prank_nt", "prank_codon", "macse2", "muscle"]


def selection_parser():
    """Parse CMD args."""
    app = argparse.ArgumentParser(
        prog="easel",
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=_help_formatter,
    )

    ## ---------------------------------------------------------------- input ---
    group_input = app.add_argument_group("Input (exactly one mode required)")
    use_TOGA = group_input.add_mutually_exclusive_group(required=True)

    use_TOGA.add_argument(
        "-toga2", "--toga2_reference_path",
        dest="toga2_directory", metavar="DIR", default=None,
        help="""
      Directory holding the TOGA2 annotation runs to extract alignments from,
      i.e. the parent of the per-query 'vs_<assembly>' directories. The
      reference name is taken from the component above it, so a path ending in
      <reference>/TOGA2 is expected. Requires --twoBit_path.

      Example: /path/to/genomes/hg38/TOGA2
      """)

    use_TOGA.add_argument(
        "-free", "--free_mode",
        dest="free_directory", metavar="DIR", default=None,
        help="""
      Directory of your own per-transcript FASTA files, one file per transcript.
      Combine with --align to align them, or with --do_screen_only if they are
      already aligned.
      """)

    use_TOGA.add_argument(
        "-toga", "--toga_reference_path",
        dest="toga_directory", metavar="DIR", default=None,
        help="""
      TOGA v1 annotations. NOT YET SUPPORTED in this release: the alignment
      extraction rules for TOGA v1 are not part of the repository. Use -toga2,
      or pass already-extracted alignments with -free.
      """)

    group_input.add_argument(
        "-asm", "--assemblies",
        dest="assembly_list", metavar="FILE", default=None,
        help="""
      Single-column list of assemblies to include. Required with -toga2.
      Optional in -free mode, where leaving it out keeps every sequence found.
      Entries may be written either as 'name' or 'vs_name'.
      """)

    group_input.add_argument(
        "-sb", "--selected_bed_file",
        dest="selected_transcripts", metavar="BED", default=None,
        help="""
      BED12 file of the transcripts to screen; column 4 supplies the transcript
      names. Required with -toga2. Not used in -free mode, where every FASTA
      file in the input directory is run.
      """)

    group_input.add_argument(
        "--twoBit_path",
        dest="twoBitPath", metavar="FILE", default=None,
        help="""
      2bit file of the reference exons, passed to TOGA2 alignment extraction.
      Required with -toga2.
      """)

    ## ------------------------------------------------------------ alignment ---
    group_align = app.add_argument_group("Alignment")

    group_align.add_argument(
        "-a", "--align",
        dest="aligner", metavar="NAME", default=None, choices=ALIGNER_CHOICES,
        help="""
      Aligner to use (case sensitive). Leave out only if the input is already
      aligned, and then set --do_screen_only.

      * prank        - TOGA2 input; nucleotide mode, codonified afterwards
      * prank_nt     - free input, nucleotide mode, codonified afterwards
      * prank_codon  - free input, codon mode
      * macse2       - TOGA2 input/free input, codon-aware
      * muscle       - TOGA2 input/free input, nucleotide mode, codonified afterwards

      prank, prank_nt and muscle are codonified against a reference sequence,
      so they need --reference_name in -free mode.
      """)

    ## ---------------------------------------------------------------- trees ---
    group_tree = app.add_argument_group("Trees (choose one)")

    group_tree.add_argument(
        "-it", "--input_tree",
        dest="input_tree", metavar="FILE", default=None,
        help="""
      Newick species tree, or bare topology, used for PRANK and for the
      selection screens. It is pruned per transcript to the species actually
      present in that alignment.
      """)

    group_tree.add_argument(
        "-ct", "--comp_tree",
        dest="tree", action="store_true", default=False,
        help="""
      Compute a gene tree per transcript with IQ-TREE 3 instead of using a
      fixed species tree.
      """)

    ## ISSUE: Precomputed gene trees should be supported with TOGA2
    group_tree.add_argument(
        "--input_gene_trees",
        dest="input_gene_trees", metavar="DIR", default=None,
        help="""
      Directory of precomputed gene trees, one per transcript, named
      <transcript><ext>. Not supported together with -toga2.
      """)

    group_tree.add_argument(
        "-st", "--comp_species_tree",
        dest="species_tree", action="store_true", default=False,
        help="""
      Reconstruct a species tree from the gene trees with ASTRAL-III.
      NOT YET IMPLEMENTED - easel exits if this is set.
      """)

    ## selection
    group_sel = app.add_argument_group("Selection analyses (all off unless enabled)")

    group_sel.add_argument(
        "-ab", "--absrel",
        dest="absrel", metavar="MODES", default=None,
        help="""
      Run aBSREL. Comma-separated list of: std, srv, mh, all. 'std' means
      neither srv nor mh, and cannot be combined with them.

      Example: -ab srv,mh
      """)

    group_sel.add_argument(
        "-bu", "--busted",
        dest="busted", metavar="MODES", default=None,
        help="""
      Run BUSTED. Comma-separated list of: std, srv, mh, error_sink, all.
      error_sink enables BUSTED-E and needs HyPhy >= 2.5.58.

      Example: -bu srv,error_sink
      """)

    group_sel.add_argument(
        "-me", "--meme",
        dest="meme", metavar="MODES", default=None,
        help="""
      Run MEME. Comma-separated list of: std, srv, mh, all.
      """)

    group_sel.add_argument(
        "-re", "--relax",
        dest="relax", metavar="MODES", default=None,
        help="""
      Run RELAX. Comma-separated list of: std, srv, mh, all, plus one integer
      giving the number of replicate runs to average over RELAX's stochasticity
      (default 10). Requires --foreground_list.

      Example: -re srv,10
      """)

    group_sel.add_argument(
        "-bc", "--bayescode",
        dest="bayescode", action="store_true", default=False,
        help="""
      Run BayesCode (mutation-selection per-site omega estimation).
      """)

    ## POTENTIAL ISSUE: We should implement disparate/partial Foreground lists
    group_sel.add_argument(
        "--foreground_list",
        dest="foreground_lst", metavar="FILE", default=None,
        help="""
      Single-column list of foreground species, which are labelled {Foreground}
      in the tree. Required for RELAX. Must be a subset of --assemblies.
      """)

    ##  run control
    group_run = app.add_argument_group("Run control")

    group_run.add_argument(
        "--directory_name",
        dest="directory_name", metavar="NAME", default="snakemake_selection_screen",
        help="""
      Name of the working directory created for this screen.
      (default: %(default)s)
      """)

    group_run.add_argument(
        "-dr", "--dry_run",
        dest="dry_run", action="store_true", default=False,
        help="""
      Write DEF.yaml, then run snakemake --dry-run to show what would be
      executed. Takes precedence over -rs.
      """)

    group_run.add_argument(
        "-rs", "--run_snakemake",
        dest="run_snakemake", action="store_true", default=False,
        help="""
      Write DEF.yaml and launch the pipeline. Without -dr or -rs, easel only
      writes DEF.yaml and stops.
      """)

    group_run.add_argument(
        "-f", "--force_run",
        dest="force_run", action="store_true", default=False,
        help="""
      Overwrite an existing DEF.yaml whose settings differ from the current
      command line. Without -f, easel refuses to touch the directory.
      """)

    group_run.add_argument(
        "--local",
        dest="local_run", action="store_true", default=False,
        help="""
      Run on this machine instead of submitting to SLURM. Without this, the
      bundled prof/config.yaml profile is used, which requires a cluster.
      """)

    group_run.add_argument(
        "--cores",
        dest="cores", metavar="N", type=int, default=4,
        help="""
      Cores to use with --local. (default: %(default)s)
      """)

    group_run.add_argument(
        "--profile",
        dest="profile", metavar="DIR", default=None,
        help="""
      Snakemake profile directory to use instead of the bundled prof/. Point
      this at your own copy to set the SLURM partition and account.
      """)

    group_run.add_argument(
        "--conda_prefix",
        dest="conda_prefix", metavar="DIR", default=None,
        help="""
      Where to create the per-rule conda environments. Defaults to
      .snakemake/conda inside easel's own installation directory, shared
      across every run directory, so environments are built once instead of
      per run.
      """)

    group_run.add_argument(
        "--apptainer_prefix",
        dest="apptainer_prefix", metavar="DIR", default=None,
        help="""
      Where to cache pulled apptainer/singularity images (.sif files).
      Defaults to .snakemake/apptainer inside easel's own installation
      directory, shared across every run directory, so images are pulled
      once instead of per run.
      """)

    group_run.add_argument(
        "--rerun_triggers_mtime",
        dest="rerun_trigger", action="store_true", default=False,
        help="""
      Pass '--rerun-triggers mtime' to snakemake, so that a changed rule does
      not force everything downstream to rerun. (debugging)
      """)

    ## advanced
    group_adv = app.add_argument_group("Advanced")

    group_adv.add_argument(
        "--do_alignment_only",
        dest="doAlignmentOnly", action="store_true", default=False,
        help="""
      Stop after alignment creation; skip trees and all selection analyses.
      """)

    group_adv.add_argument(
        "--do_screen_only",
        dest="doScreenOnly", action="store_true", default=False,
        help="""
      Start at alignment postprocessing. The input must already be aligned,
      so this only makes sense with -free.
      """)

    group_adv.add_argument(
        "--reference_name",
        dest="reference_name", metavar="NAME", default=None,
        help="""
      Name of the reference sequence used for codonification. Required in -free
      mode with a codonifying aligner (prank, prank_nt, muscle).
      """)

    ## ISSUE: TOGA2 will probably be containerized, and then the container image should be used
    ## in all instances
    group_adv.add_argument(
        "--toga2_activate",
        dest="toga2Activate", metavar="FILE", default=None,
        help="""
      Path to an activate script that puts the toga2 executable on PATH. TOGA2
      is not packaged on bioconda, so it has to come from the host. Omit if
      toga2 is already on PATH.
      """)

  ## ISSUE: This relates to the old TOGA1 alignment extraction script
    group_adv.add_argument(
        "--extract_ali_params",
        dest="extract_ali_params", metavar="LIST", default="skip_dups",
        help="""
      Comma-separated flags passed to TOGA2 alignment extraction:

      * skip_dups      - drop non-one2one orthologs
      * allow_one2zero - keep lost transcripts
      * align_entirely - whole-gene instead of exon-by-exon
      * exclude_UL     - drop 'uncertain loss' sequences

      (default: %(default)s)
      """)

    group_adv.add_argument(
        "--max_CDS_length",
        dest="maxCDSLengthReference", metavar="BP", type=int, default=15000,
        help="Skip transcripts whose CDS is longer than this. (default: %(default)s)")

    group_adv.add_argument(
        "--min_CDS_length",
        dest="minCDSLengthReference", metavar="BP", type=int, default=50,
        help="Skip transcripts whose CDS is shorter than this. (default: %(default)s)")

    group_adv.add_argument(
        "--min_num_aligned_species",
        dest="minNumAlignedSpecies", metavar="N", type=int, default=5,
        help="""
      Skip transcripts present in fewer than N species. HyPhy needs at least 3
      taxa, IQ-TREE bootstrapping at least 4. (default: %(default)s)
      """)

    group_adv.add_argument(
        "--do_hmm_cleaning",
        dest="doHMMCleaning", action="store_true", default=False,
        help="""
      Clean alignments with HmmCleaner.pl. HmmCleaner is not on bioconda, so
      it runs from the ghcr.io/hillerlab/hmmcleaner container image instead
      of requiring a host install; this needs apptainer/singularity to be
      available wherever the pipeline actually runs.
      """)

    group_adv.add_argument(
        "--hmm_cleaner_params",
        dest="hmm_cleaner_params", metavar="LIST", default="0.15,0.08,0.15,0.45",
        help="""
      Four HmmCleaner cost values c1,c2,c3,c4. The first two are negated, and
      they must increase: c1 < c2 < 0 < c3 < c4. (default: %(default)s)
      """)

    group_adv.add_argument(
        "--do_manual_cleaning",
        dest="doManualCleaning", action="store_true", default=False,
        help="Clean alignments with the built-in column and row filter.")

    group_adv.add_argument(
        "--manual_cleaner_params",
        dest="manual_cleaner_params", metavar="LIST", default="mc0.6,ms0.3,ml25,m",
        help="""
      Comma-separated filter settings:

      * mc<f> - min fraction of a column that must align for it to stay
      * ms<f> - min fraction of a sequence that must align for it to stay
      * ml<n> - min sequence length in amino acids
      * m     - mask ambiguous and stop codons with NNN

      (default: %(default)s)
      """)

    group_adv.add_argument(
        "--max_failed_fraction",
        dest="max_failed_fraction", metavar="F", type=float, default=0.4,
        help="""
      Stop before the analyses if more than this fraction of transcripts
      produced no alignment at all (phase 1 failures, as opposed to alignments
      that were produced and rejected by validation).

      A few aligner failures no longer take the whole run down: those
      transcripts are excluded and listed in skipped_transcripts.tsv.

      (default: %(default)s)
      """)

    group_adv.add_argument(
        "--do_error_sink_cleaning",
        dest="doErrorCleaning", action="store_true", default=False,
        help="""
      Filter alignments by the empirical Bayes factors from BUSTED-E. Implies
      -bu error_sink.
      """)

    ## resources
    group_res = app.add_argument_group("Resources (per rule)")

    def add_resource(flag, threads, mem_mb, label):
        group_res.add_argument(
            f"--{flag}_threads", dest=f"{flag}_threads", metavar="N",
            type=int, default=threads,
            help=f"Threads for {label}. (default: %(default)s)")
        group_res.add_argument(
            f"--{flag}_mem_mb", dest=f"{flag}_mem_mb", metavar="MB",
            type=int, default=mem_mb,
            help=f"Memory in MB for {label}. (default: %(default)s)")

    add_resource("extract_alignments", 1, 10000, "alignment extraction")
    add_resource("hmmcleaner", 1, 1000, "HmmCleaner")
    add_resource("manualcleaner", 1, 1000, "the manual cleaner")
    add_resource("tree_comp", 4, 15000, "tree computation")
    add_resource("prank", 5, 5000, "PRANK")
    add_resource("absrel", 5, 10000, "aBSREL")
    add_resource("busted", 5, 10000, "BUSTED")
    add_resource("meme", 5, 10000, "MEME")
    add_resource("relax", 10, 10000, "RELAX")

    args = app.parse_args()
    return args

## ---------------------------------------------------------------------------
## Utility helpers
## ---------------------------------------------------------------------------

## Function to reformat lists to be represented as python lists in yaml files
def format_list(in_list):
    formatted = comments.CommentedSeq(in_list)
    formatted.fa.set_flow_style()
    return formatted

## Function to parse Hyphy commands from command line input
## For the future: allow more fine-grained control over srv and MH
def parse_hyphy_modes(in_command, hyphy_program, config):
    if not in_command:
        config["settings"]["selectionSettings"][str(hyphy_program)]["activate"] = False
        logger.info("%s is deactivated for the current selection screen", str(hyphy_program))
        if hyphy_program == "BUSTED":
            config["settings"]["selectionSettings"]["BUSTED"]['error_sink'] = False
        return
    commands = in_command.split(",")
    if "srv" in commands:
        srv_command = " Yes "
    else:
        srv_command = " No "

    if "mh" in commands:
        mh_command = " Double+Triple"
    else:
        mh_command = " None "

    if hyphy_program == "BUSTED":
        if 'error_sink' in commands:
            config["settings"]["selectionSettings"]["BUSTED"]['error_sink'] = True
        else:
            config["settings"]["selectionSettings"]["BUSTED"]['error_sink'] = False

    if commands == ["all"]:
        srv_command = " Yes "
        mh_command = " Double+Triple"
        if hyphy_program == 'BUSTED':
            config["settings"]["selectionSettings"]["BUSTED"]['error_sink'] = True

    config["settings"]["selectionSettings"][str(hyphy_program)]["srv"] = srv_command
    config["settings"]["selectionSettings"][str(hyphy_program)]["multi_hit"] = mh_command

    logger.info("%s is activated with %s enabled", str(hyphy_program), in_command)

    if any(single_command.isdigit() for single_command in commands) and hyphy_program == 'RELAX':
        rounds = [single_command for single_command in commands if single_command.isdigit()]
        rounds = int(''.join(str(integer) for integer in rounds))
        config["settings"]["selectionSettings"][str(hyphy_program)]["relaxRounds"] = rounds
        logger.info("%s is activated with %d rounds", str(hyphy_program), rounds)

## function to parse comma separated lists, which is how parameters are passed to 
## extract codon alignment, hmm cleaner and manual cleaner
def parse_parameter_lists(comma_separated_list, program_name):
    def text_num_split(item):
        for index, letter in enumerate(item, 0):
            if letter.isdigit():
                return item[:index] + ' ' + item[index:]
        return item

    ## ISSUE: Extract should be updated for TOGA2 sequence-alignment
    if program_name == "extract":
        out_string = " ".join(["--" + entry for entry in comma_separated_list.split(",")])
    elif program_name == "hmm":
        out_string = " ".join(["-" + entry  for entry in comma_separated_list.split(",")[0:2]]) + " " + " ".join([entry for entry in comma_separated_list.split(",")[2:]])
    elif program_name == "manual":
        out_string = " ".join(["-" + text_num_split(entry) for entry in comma_separated_list.split(",")])
    return out_string

import re

## ISSUE: Aren't these parameters already processed above?
def parse_manual_cleaner_params(comma_separated_list, param_name):
    entries = [e.strip() for e in comma_separated_list.split(",") if e.strip()]
    for item in entries:
        match = re.match(r"([a-zA-Z]+)([0-9.]*)", item)
        
        if match:
            current_label, current_val = match.groups()
            if current_label == param_name:
              ## In case it's m for mask
                if current_val == "":
                    return True
                try:
                    return float(current_val)
                ## Fallback  (change this at some point)    
                except ValueError:
                    return True 
                    
    ## Parameter not found in the list
    return False


## Input helper functions (from polymeval)
def get_snakefile_path(name="Snakefile"):
    snakefile = os.path.join(BASE_DIR, name)
    return snakefile

def get_cluster_configfile_path():
    """Path to the bundled snakemake profile directory."""
    return os.path.join(BASE_DIR, "prof")


def _ensure_directory(path):
    """Create *path* if it doesn't exist yet.
    """
    os.makedirs(path, exist_ok=True)
    return path


def _log_run_plan(config):
    """State what this run will actually do, before handing over to snakemake.

    On a fresh run there are no verdicts yet, so phase 2's transcript list --
    and with it every analysis rule -- is empty until phase 1 has finished.
    Currently, this is not listing every job in the pipeline (some smaller logging jobs are missing)
    """
    sel = config["settings"]["selectionSettings"]
    aln = config["settings"]["alignmentSettings"]
    cln = config["settings"]["cleaningSettings"]
    tre = config["settings"]["treeSettings"]
    onoff = lambda x: "on" if x else "off"

    mode = ("TOGA2" if config["toga2Mode"]
            else "free" if config["freeMode"] else "TOGA v1")
    analyses = [name for name, key in (("aBSREL", "ABSREL"), ("BUSTED", "BUSTED"),
                                       ("MEME", "MEME"), ("RELAX", "RELAX"),
                                       ("BayesCode", "bayesCode"))
                if sel.get(key, {}).get("activate")]

    logger.info("Run plan")
    logger.info("  mode              : %s, aligner %s", mode, aln["aligner"])
    logger.info("  alignment only    : %s", onoff(aln["doAlignmentOnly"]))
    logger.info("  cleaning          : HMM %s, manual %s",
                onoff(cln["hmmCleaning"]["doHMMCleaning"]),
                onoff(cln["manualCleaning"]["doManualCleaning"]))
    logger.info("  trees             : gene %s, input species tree %s, precomputed %s",
                onoff(tre["computeGeneTrees"]["activate"]),
                onoff(tre["inputSpeciesTree"]["treeFile"]),
                onoff(tre["preCompGeneTrees"]["geneTreePath"]))
    logger.info("  analyses          : %s", ", ".join(analyses) or "none")
    logger.info("  BUSTED error sink : %s, error cleaning %s",
                onoff(sel["BUSTED"]["error_sink"]), onoff(sel["doErrorCleaning"]))
    if analyses:
        logger.info("  NOTE: these run in phase 2, on the transcripts phase 1 "
                    "found usable, so phase 1's job table will not list them.")


## ---------------------------------------------------------------------------
## The two phases
## ---------------------------------------------------------------------------
## A run is two snakemake invocations over the same Snakefile, distinguished by
## a config key that the workflow reads at parse time:
##
##   phase 1  target `validated`, easel_phase=validate
##            align, clean and validate every transcript, with --keep-going.
##            Ends by writing VERDICT_SUMMARY (see finish_validation_phase in
##            rules/common.smk).
##   phase 2  the default target, easel_phase=analyse
##            re-parse; the verdicts are on disk now, so rejected transcripts
##            are filtered out at PARSE time and never enter the DAG.
##

PHASE_VALIDATE = "validate"
PHASE_ANALYSE = "analyse"
VERDICT_SUMMARY = "validation_summary.tsv"


def _phase_cmd(base, target, phase, extra=(), config=()):
    """`base` plus a target, a phase, and any extra flags.

    The target goes BEFORE the options. --config and --quiet each take a
    variable number of values, so a target placed after one of them is read as
    another of its values rather than as a target --config goes last and carries every
    key at once, because argparse keeps only the final --config it sees.
    """
    cmd = list(base)
    if target:
        cmd.insert(1, target)
    cmd += [flag for flag in extra if flag not in cmd]
    return cmd + ["--config", f"easel_phase={phase}", *config]


def _read_verdict_counts(snakemake_dir):
    """-> (counts, total) from phase 1's summary, or (None, 0) if unusable."""
    counts = {"OK": 0, "SKIP": 0, "FAIL": 0}
    total = 0
    try:
        with open(os.path.join(snakemake_dir, VERDICT_SUMMARY),
                  encoding="utf-8") as fh:
            if not fh.readline().startswith("transcript\t"):
                return None, 0
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2 or not fields[0]:
                    continue
                total += 1
                counts[fields[1]] = counts.get(fields[1], 0) + 1
    except OSError:
        return None, 0
    return (counts, total) if total else (None, 0)


def _validation_gate(snakemake_dir, max_failed_fraction):
    """Should phase 2 run? -> (proceed, message, level).

    The point of the two-phase design is that a handful of aligner failures no
    longer takes the whole run down -- those transcripts are simply excluded and
    reported. That tolerance needs a brake, or a run in which 90% of the
    alignments died would quietly "succeed" on the remaining 10%. This is the
    brake.
    """
    counts, total = _read_verdict_counts(snakemake_dir)
    if counts is None:
        return False, (
            f"Phase 1 left no usable {VERDICT_SUMMARY} in {snakemake_dir}, so "
            "there are no verdicts to analyse. Its own error output is above; "
            "nothing was skipped silently."), "critical"

    failed, ok = counts.get("FAIL", 0), counts.get("OK", 0)
    share = failed / total
    tally = (f"{ok} usable, {counts.get('SKIP', 0)} rejected by validation, "
             f"{failed} with no alignment, of {total}")

    if not ok:
        return False, (
            f"No transcript passed validation ({tally}). Nothing to analyse -- "
            f"see skipped_transcripts.tsv and {VERDICT_SUMMARY} for why."
        ), "critical"
    if share > max_failed_fraction:
        return False, (
            f"{failed} of {total} transcripts ({share:.1%}) produced no "
            f"alignment at all, over the --max_failed_fraction of "
            f"{max_failed_fraction:.1%}. Stopping before the analyses rather "
            f"than reporting on {ok} transcripts as if the run were complete. "
            f"Look at codon_alignments/<id>/ for a failing transcript, then "
            f"either fix it or raise --max_failed_fraction to continue."
        ), "critical"
    if failed:
        return True, (
            f"Phase 1: {tally}. The {failed} without an alignment are excluded "
            f"and listed in {VERDICT_SUMMARY}; their tmp/ directories are left "
            f"in place to inspect."), "warning"
    return True, f"Phase 1: {tally}.", "info"


## Transcripts to keep for the two DAG builds that only DESCRIBE the
## workflow. 
GRAPH_SAMPLE = 3


def _log_max_job_estimate(snake_file, snakemake_dir):
    """Log the job counts the workflow would produce if every transcript passed.

    Phase 1's job table covers alignment and cleaning only, and phase 2's
    cannot be known until phase 1 has run. rulegraph_reference requests the full
    analysis output set with assumeAllUsable, so a dry run against it resolves
    completely and yields the upper bound: what would be spawned if every
    alignment were accepted.
    """
    ## --rerun-incomplete because this runs WITHOUT --profile (see above), which
    ## also drops the profile's rerun-incomplete: True. Any output left
    ## half-written by an earlier failed run then makes snakemake refuse with an
    ## IncompleteFilesException. Harmless here: a dry run rewrites nothing.
    ## --config graphSample: resolve_run_inputs (rules/common.smk) truncates
    ## the transcript list, and the per-rule counts below are divided back out
    ## and multiplied by the real total. 
    ## assumeAllUsable: this is an estimate of what a run WOULD produce, so it
    ## must not be filtered by whatever verdicts happen to be on disk from a
    ## previous run. Without it, a rerun in a directory where most transcripts
    ## were rejected would report a far smaller workflow than the next run will
    ## actually have.
    cmd = ["snakemake", "rulegraph_reference",
           "--snakefile", snake_file, "--directory", snakemake_dir,
           "--dry-run", "--rerun-incomplete",
           "--config", f"graphSample={GRAPH_SAMPLE}", "assumeAllUsable=True"]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             check=True, text=True)
    except subprocess.CalledProcessError as exc:
        ## Show snakemake's own message. 
        out = (exc.stdout or "").strip().splitlines()
        reasons = [ln for ln in out
                   if any(k in ln for k in ("Error", "Exception", "error:"))]
        logger.warning("Could not estimate the maximum job count "
                       "(snakemake exited %s).", exc.returncode)
        for line in reasons[:5] or out[:5]:
            logger.warning("  %s", line.rstrip())
        if len(out) > 5:
            logger.warning("  ... last lines:")
            for line in out[-8:]:
                logger.warning("  %s", line.rstrip())
        return
    except OSError as exc:
        logger.warning("Could not estimate the maximum job count: %s", exc)
        return
    lines = res.stdout.splitlines()
    try:
        start = next(i for i, ln in enumerate(lines)
                     if ln.strip().startswith("Job stats:"))
    except StopIteration:
        logger.info("No job-stats table in the dry-run output; skipping estimate")
        return
    m = re.search(r"graphSample=(\d+) of (\d+) transcripts", res.stdout)
    sample, total = (int(m.group(1)), int(m.group(2))) if m else (None, None)

    rows = []
    for ln in lines[start + 1:]:
        if not ln.strip():
            break
        rows.append(ln.rstrip())

    if not sample or not total or sample >= total:
        logger.info("Maximum jobs, if every transcript passes validation:")
        for ln in rows:
            logger.info("  %s", ln)
    else:
        logger.info("Maximum jobs for %s transcripts, if every one passes "
                    "validation (projected from a %s-transcript sample):",
                    f"{total:,}", sample)
        grand = 0
        for ln in rows:
            fields = ln.split()
            ## The table is "job  count", plus a header and a "total" row.
            if len(fields) >= 2 and fields[-1].isdigit():
                name, count = " ".join(fields[:-1]), int(fields[-1])
                if name.lower() == "total":
                    continue
                if count % sample == 0:
                    ## Once per transcript (or a fixed number of times per
                    ## transcript, e.g. RELAX's rounds).
                    projected = count // sample * total
                    grand += projected
                    logger.info("  %-34s %12s   (%s per transcript)",
                                name, f"{projected:,}", count // sample)
                else:
                    ## Once per run: skip_report, and the target rule itself.
                    grand += count
                    logger.info("  %-34s %12s   (once per run)",
                                name, f"{count:,}")
            else:
                logger.info("  %s", ln)
        logger.info("  %-34s %12s", "total", f"{grand:,}")
    logger.info("  (actual counts will be lower by whatever validate_alignment "
                "rejects; rerun this yourself with: snakemake --snakefile %s "
                "--directory %s --dry-run rulegraph_reference)",
                snake_file, snakemake_dir)


def _write_rulegraph(snake_file, snakemake_dir, force=False):
    """Render the rule graph to rulegraph.pdf, if graphviz is available."""
    out_path = os.path.join(snakemake_dir, "rulegraph.pdf")
    ## Regenerate when the workflow or the config is newer than the figure.
    if os.path.exists(out_path) and not force:
        newest = 0.0
        for pattern in ("Snakefile_*", os.path.join("rules", "*.smk"),
                        os.path.join("easel", "*.py")):
            for path in glob.glob(os.path.join(BASE_DIR, pattern)):
                newest = max(newest, os.path.getmtime(path))
        def_yaml = os.path.join(snakemake_dir, "DEF.yaml")
        if os.path.exists(def_yaml):
            newest = max(newest, os.path.getmtime(def_yaml))
        if os.path.getmtime(out_path) >= newest:
            logger.info("rulegraph.pdf is up to date; delete it to force a rebuild")
            return
        logger.info("rulegraph.pdf is older than the workflow or DEF.yaml; regenerating")
    if shutil.which("dot") is None:
        logger.warning(
            "graphviz (`dot`) not found, skipping rulegraph.pdf. "
            "Install it with: conda install -c conda-forge graphviz"
        )
        return
    base = ["snakemake", "--snakefile", snake_file, "--directory", snakemake_dir,
            "--rerun-incomplete", "--config", f"graphSample={GRAPH_SAMPLE}",
            "assumeAllUsable=True"]
    graph_flag = ["--rulegraph"]

    def _graph(cmd):
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              check=True)
    try:
        try:
            dag = _graph(base[:1] + ["rulegraph_reference"] + base[1:] + graph_flag)
        except subprocess.CalledProcessError as exc:
            logger.warning("Could not graph rulegraph_reference (snakemake "
                           "exited %s); falling back to the default target, so "
                           "analysis rules past the checkpoint will be MISSING. "
                           "Its output was:", exc.returncode)
            for line in (exc.stderr or b"").decode(errors="replace").strip().splitlines()[-15:]:
                logger.warning("  %s", line)
            dag = _graph(base + graph_flag)
        except OSError as exc:
            logger.warning("Could not run snakemake for the rulegraph: %s", exc)
            return
        with open(out_path, "wb") as fh:
            subprocess.run(["dot", "-Tpdf"], input=dag.stdout, stdout=fh, check=True)
        logger.info("Rule graph written to %s", out_path)
    except (subprocess.CalledProcessError, OSError) as exc:
        ## A missing figure must never take down the run.
        logger.warning("Could not render rulegraph.pdf: %s", exc)


## Retained, but inert: snakemake emits this only for checkpoint jobs, and the
## workflow has no checkpoints since validate_alignment is a plain rule.
## Kept so that the filter is already in place if one is ever reintroduced, and
## so a reader who finds the message in an old log knows where it came from.

CHECKPOINT_NOTICE = "DAG of jobs will be updated after completion."

## Dropped in the run directory when the user interrupts, and read by the
## workflow's onerror hook (run_was_interrupted in rules/common.smk) so it can
## skip the log bundling and run report. snakemake's onerror cannot tell a
## cancelled run from a failed job on its own.
INTERRUPT_MARKER = ".easel_interrupted"


def _clear_interrupt_marker(cwd):
    """Remove a marker left behind by an earlier run.

    Called before launching and after the process exits, so a session that was
    hard-killed between the two cannot leave a marker that silences the report
    on someone's next, genuinely failing run.
    """
    try:
        os.remove(os.path.join(cwd, INTERRUPT_MARKER))
    except OSError:
        pass


def _is_checkpoint_notice(line):
    return line.strip() == CHECKPOINT_NOTICE


def _invoke_snakemake(cmd, cwd, drop=None):
    """Run one snakemake invocation. Returns (returncode, interrupted).

    *drop*, if given, is a predicate on output lines: snakemake's output is
    relayed through this process and matching lines are left out, with a count
    reported afterwards so nothing is discarded silently.
    """
    ## Forward Ctrl-C to the whole process group so snakemake can cancel its
    ## cluster jobs and clean up. os.setpgrp/os.killpg are POSIX-only.
    popen_kwargs = {}
    if os.name == "posix":
        popen_kwargs["preexec_fn"] = os.setpgrp
    if drop is not None:
        popen_kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, errors="replace", bufsize=1)

    dropped = 0
    interrupted = False
    _clear_interrupt_marker(cwd)
    process = subprocess.Popen(cmd, cwd=cwd, **popen_kwargs)

    def relay():
        """Pass output through, minus the dropped lines. Returns the count."""
        n = 0
        for line in process.stdout:
            if drop(line):
                n += 1
                continue
            sys.stdout.write(line)
        sys.stdout.flush()
        return n

    try:
        if drop is not None:
            dropped = relay()
        process.wait()
    except KeyboardInterrupt:
        interrupted = True
        logger.info("Intercepted Ctrl-C, forwarding to snakemake ...")
        try:
            with open(os.path.join(cwd, INTERRUPT_MARKER), "w") as fh:
                fh.write("interrupted\n")
        except OSError as exc:
            logger.warning("Could not record the interrupt (%s); the run "
                           "report will still be written.", exc)
        try:
            if os.name == "posix":
                os.killpg(process.pid, signal.SIGINT)
            else:
                process.terminate()
        except OSError as exc:
            logger.warning("Could not signal snakemake: %s", exc)
        logger.info("Waiting for snakemake to finish cleaning up ...")
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        if drop is not None:
            ## Keep relaying while it shuts down: the cleanup messages are the
            ## ones worth reading after an interrupt.
            try:
                dropped += relay()
            except (OSError, ValueError):
                pass
        process.wait()
        ## The hook has run by now, so the marker has done its job.
        _clear_interrupt_marker(cwd)

    if dropped:
        logger.info(
            "Left out %d repeated %r line(s), one per checkpoint job. "
            "Snakemake emits it without the metadata its own --quiet filter "
            "keys on, so it cannot be switched off there.",
            dropped, CHECKPOINT_NOTICE)
    return process.returncode, interrupted


## function to run snakemake (adapted from polymeval)
def run_snakemake(snake_file,
                  fasta_path,
                  directory_name="easel_run",
                  dryrun=True,
                  snake_default=False,
                  rerun_triggers=False,
                  working_dir=None,
                  updated_rule=False,
                  local_run=False,
                  cores=4,
                  profile=None,
                  conda_prefix=None,
                  apptainer_prefix=None,
                  max_failed_fraction=0.4):
    """Build and run the snakemake command line(s). Returns the exit code.

    A real run is two invocations -- see the "The two phases" section above. A
    dry run stays a single one: it targets rulegraph_reference, whose plan is
    already the upper bound over both phases, and a second full DAG build would
    double the wait on a large dataset for no extra information.
    """
    working_dir = working_dir or os.getcwd()
    snakemake_dir = _ensure_directory(os.path.join(working_dir, directory_name))

    snake_file = get_snakefile_path(snake_file)
    cmd = ["snakemake", "--snakefile", snake_file, "--directory", snakemake_dir]

    ## Whether to do a local run (no remote execution, head node of an HPC or local machine)
    ## Or a cluster run
    if local_run:
        cmd += ["--cores", str(cores), "--max-threads", str(cores)]
    else:
        cmd += ["--profile", profile or get_cluster_configfile_path()]

    ## Conda prefix must point to easel base directory, otherwise each run will download
    ## all dependencies again, which is unnecessary overload
    cmd += ["--use-conda", "--conda-prefix",
            conda_prefix or os.path.join(BASE_DIR, ".snakemake", "conda")]

    ## Needed for any rule declaring `container:` (e.g. hmm_cleaner/
    ## transfer_cleaner in rules/hmm_cleaning.smk) to actually run in that
    ## container -- without it Snakemake silently falls back to whatever is
    ## on PATH, same as a conda: directive without --use-conda.
    ## Passed unconditionally, same as --use-conda above, so it applies
    ## whether or not the bundled profile (which also sets use-apptainer:
    ## True) is in play.
    ##
    ## --apptainer-prefix mirrors --conda-prefix just above: without it,
    ## Snakemake's default is .snakemake/apptainer under --directory (the run
    ## directory), so every run directory re-pulls the same .sif images that
    ## conda environments already share out of BASE_DIR.
    cmd += ["--use-apptainer", "--apptainer-prefix",
            apptainer_prefix or os.path.join(BASE_DIR, ".snakemake", "apptainer")]

    if rerun_triggers:
        cmd += ["--rerun-triggers", "mtime"]
    if snake_default:
        cmd += ["--rerun-incomplete", "--keep-going"]

    ## Rulegraph, creates PDF for users to review.
    _write_rulegraph(snake_file, snakemake_dir, force=updated_rule)

    if dryrun:
        ## --quiet rules: a dry run over thousands of transcripts otherwise
        ## prints one block per job and buries the summary. The job-stats table
        ## and the reasons still print; only the per-job listing is dropped.

        dry = _phase_cmd(cmd, "rulegraph_reference", PHASE_ANALYSE,
                         config=["assumeAllUsable=True"])
        dry += ["--dry-run", "--quiet", "rules"]
        logger.info("Dry run targets rulegraph_reference: the plan below is the "
                    "UPPER BOUND over both phases (every transcript assumed to "
                    "pass validation).")
        logger.info("Running: %s", " ".join(dry))
        returncode, interrupted = _invoke_snakemake(
            dry, snakemake_dir, drop=_is_checkpoint_notice)
        if returncode == 0:
            logger.info("Dry run finished successfully.")
        else:
            logger.critical("Dry run exited with code %d.", returncode)
        return (returncode or 1) if interrupted else returncode

    ## ---- phase 1: align, clean, validate ---------------------------------
    phase1 = _phase_cmd(cmd, "validated", PHASE_VALIDATE, extra=["--keep-going"])
    logger.info("Phase 1/2 -- align, clean and validate every transcript")
    logger.info("Running: %s", " ".join(phase1))
    returncode, interrupted = _invoke_snakemake(phase1, snakemake_dir)
    if interrupted:
        logger.warning("Interrupted during phase 1. Nothing was analysed; "
                       "rerun the same command to continue where it stopped.")
        return returncode or 1

    ## A non-zero exit code here is EXPECTED whenever some alignments failed,
    ## because of --keep-going. So whether to continue is decided from the
    ## verdicts, not from the exit code.
    proceed, message, level = _validation_gate(snakemake_dir, max_failed_fraction)
    getattr(logger, level)("%s", message)
    if not proceed:
        return returncode or 1

    ## ---- phase 2: trees, screens, reports --------------------------------
    phase2 = _phase_cmd(cmd, None, PHASE_ANALYSE)
    logger.info("Phase 2/2 -- trees, selection screens and reports")
    logger.info("Running: %s", " ".join(phase2))
    returncode, interrupted = _invoke_snakemake(phase2, snakemake_dir)

    if returncode == 0:
        logger.info("Snakemake finished successfully.")
    else:
        logger.critical("Snakemake exited with code %d.", returncode)
    if interrupted:
        return returncode or 1
    return returncode

## ---------------------------------------------------------------------------
## Main function
## ---------------------------------------------------------------------------
def main():

    ## Logged before parsing
    logger.info("Command: %s %s",
                os.path.basename(sys.argv[0]), shlex.join(sys.argv[1:]))

    ## Read command line arguments
    args = selection_parser()


    ## --------------------------------------------------------------------------
    ## Input validation
    ## --------------------------------------------------------------------------

    ## Reject every unusable flag combination and missing path before anything is
    ## created. Collects all problems and reports them in one pass, rather than
    ## one per run. Content checks (FASTA/BED/tree parsing, tree-vs-assembly
    ## coverage, TOGA run directories, 2bit signature) all live here too, so
    ## there is a single pass and a single report -- see easel/validate.py and
    ## easel/formats.py.
    resolved = validate(args)
    report(resolved)

    ## --------------------------------------------------------------------------
    ## Input mode, reference name, assembly list and tree paths
    ## --------------------------------------------------------------------------
    ## All resolved and validated by easel.validate above -- including the
    ## assembly list (formats.read_one_column) and the tree (formats.read_in_tree
    ## / compare_tree_asm). 
    Fasta_path = resolved.fasta_path
    alignment_reference = resolved.alignment_reference
    assemblies = resolved.assemblies
    tree_file = (os.path.abspath(args.input_tree)
                 if resolved.tree_strategy == "input_tree" else None)
    gene_tree_path = (os.path.abspath(args.input_gene_trees)
                      if resolved.tree_strategy == "precomputed" else None)

    ## ISSUE: Check whether this is necessary
    ## Content already validated in validate(); resolve to an absolute path.
    if args.selected_transcripts is not None:
        selected_transcripts_file = os.path.abspath(args.selected_transcripts)
    else:
        selected_transcripts_file = None

    ## --------------------------------------------------------------------------
    ## Build YAML config dictionary
    ## --------------------------------------------------------------------------
    config = {
        "togaMode": resolved.mode == "toga",
        "toga2Mode": resolved.mode == "toga2",
        "freeMode": resolved.mode == "free",
        "fastaPath": SingleQuotedScalarString(Fasta_path),
        "absPath": SingleQuotedScalarString(os.path.join(os.getcwd(), args.directory_name)),
        "speciesLst": format_list(assemblies) if assemblies else None,
        "originalTranscripts": (SingleQuotedScalarString(selected_transcripts_file)
                                if selected_transcripts_file else None),
        "selectedTranscripts": (SingleQuotedScalarString(selected_transcripts_file)
                                if selected_transcripts_file else None),
        "referenceName": SingleQuotedScalarString(alignment_reference) if alignment_reference else None,
        "settings": {
            "alignmentSettings": {
                "doAlignmentOnly": args.doAlignmentOnly,
                "aligner": SingleQuotedScalarString(args.aligner) if args.aligner is not None else None,
                "fromTOGA": {
                  #"includeLostTranscripts": args.includeLostTranscripts, can be removed since it can be adapted in the below step
                  "twoBitPath": SingleQuotedScalarString(args.twoBitPath) if args.twoBitPath else None,
                  "toga2Activate": SingleQuotedScalarString(args.toga2Activate) if args.toga2Activate else None,
                  "extractAliParams": SingleQuotedScalarString(parse_parameter_lists(args.extract_ali_params, program_name="extract")),
                  #"maxCDSLengthReference": args.maxCDSLengthReference  IMPORTANT: Will be part of preprocess script
                },
                #"freeAlignmentPattern": args.TOGA_pattern, ## IMPORTANT: ADAPT THIS. For Precomputed alignments, whether it is mfa or not,
            },
            "cleaningSettings": {
                "manualCleaning": {
                  "doManualCleaning": args.doManualCleaning,
                  "manualCleanerRef": "REFERENCE", ## ISSUE: must be adapted, SingleQuotedScalarString(alignment_reference), 
                  "manualCleanerParams": {
                    "mincodon": parse_manual_cleaner_params(args.manual_cleaner_params,  "mc"),
                    "minseq": parse_manual_cleaner_params(args.manual_cleaner_params, "ms"), 
                    "minaalen": parse_manual_cleaner_params(args.manual_cleaner_params, "ml"), 
                    "mask": parse_manual_cleaner_params(args.manual_cleaner_params, "m"),
                  },
                },
                "hmmCleaning": {
                  "doHMMCleaning": args.doHMMCleaning,
                  "hmmCleanerParams": SingleQuotedScalarString(parse_parameter_lists(args.hmm_cleaner_params, program_name="hmm")),
                },
            },
            "treeSettings": {
              "computeGeneTrees": {
                  "activate": args.tree,
              },
              "inputSpeciesTree": {
                "treeFile": SingleQuotedScalarString(tree_file) if tree_file is not None else None,
              },
              "computeSpeciesTree": {
                "activate": args.species_tree,
              },
              "preCompGeneTrees": {
                "geneTreePath": SingleQuotedScalarString(gene_tree_path) if gene_tree_path is not None else None,
              },
            },
            "selectionSettings": {
                "doScreenOnly": args.doScreenOnly,
                "doErrorCleaning": args.doErrorCleaning,
                "ABSREL": {
                    "activate": True,
                    "srv": " Yes ",
                    "multi_hit": " None "
                },
                "BUSTED": {
                    "activate": True,
                    "srv": " Yes ",
                    "multi_hit": " Double+Triple ",
                    "error_sink": True
                },
                "MEME": {
                    "activate": True,
                    "srv": " Yes ",
                    "multi_hit": " None "
                },
                "RELAX": {
                    "activate": True,
                    "srv": " No ",
                    "multi_hit": " None ",
                    "relaxRounds": 10
                },
                "foregroundLst": SingleQuotedScalarString(args.foreground_lst),
                "bayesCode" : {
                  "activate": args.bayescode,
                },
            },
        },
        "resources": {
            "extractAlignments": {
              "threads": args.extract_alignments_threads,
              "mem_mb": args.extract_alignments_mem_mb
            },
            "hmmCleaner": {
              "threads": args.hmmcleaner_threads,
              "mem_mb": args.hmmcleaner_mem_mb
            },
            "manualCleaner": {
              "threads": args.manualcleaner_threads,
              "mem_mb": args.manualcleaner_mem_mb 
            },
            "tree": {
              "threads": args.tree_comp_threads,
              "mem_mb": args.tree_comp_mem_mb
            },
            "absrel": {
              "threads": args.absrel_threads,
              "mem_mb": args.absrel_mem_mb
            },
            "busted": {
              "threads": args.busted_threads,
              "mem_mb": args.busted_mem_mb
            },
            "meme": {
              "threads": args.meme_threads,
              "mem_mb": args.meme_mem_mb
            },
            "relax": {
              "threads": args.relax_threads,
              "mem_mb": args.relax_mem_mb
            },
            "prank": {
              "threads": args.prank_threads,
              "mem_mb": args.prank_mem_mb
            },
        }
    }

    ## --------------------------------------------------------------------------
    ## Post-config adjustments based on mode flags
    ## --------------------------------------------------------------------------
    ## Every flag combination checked below was already rejected by validate()
    ## if unusable (report() would have exited before this point), so nothing
    ## here re-checks anything -- it only translates already-valid flags into
    ## config values that couldn't be filled in above.


    if resolved.foreground:
        config["settings"]["selectionSettings"]["foregroundLst"] = format_list(resolved.foreground)
    else:
        config["settings"]["selectionSettings"]["foregroundLst"] = None

    parse_hyphy_modes(args.absrel, 'ABSREL', config)
    parse_hyphy_modes(args.busted, 'BUSTED', config)
    parse_hyphy_modes(args.meme,   'MEME',   config)
    parse_hyphy_modes(args.relax,  'RELAX',  config)

    ## --do_error_sink_cleaning implies BUSTED + error_sink.
    if args.doErrorCleaning and (
      (not config["settings"]["selectionSettings"]["BUSTED"]["activate"]) or
      (not config["settings"]["selectionSettings"]["BUSTED"]['error_sink'])
    ):
        config["settings"]["selectionSettings"]["BUSTED"]["activate"] = True
        config["settings"]["selectionSettings"]["BUSTED"]['error_sink'] = True

    ## --do_alignment_only overrides all downstream options.
    if config["settings"]["alignmentSettings"]["doAlignmentOnly"]:
        config["settings"]["selectionSettings"]["BUSTED"]["activate"] = False
        config["settings"]["selectionSettings"]["ABSREL"]["activate"] = False
        config["settings"]["selectionSettings"]["RELAX"]["activate"] = False
        config["settings"]["selectionSettings"]["MEME"]["activate"] = False
        config["settings"]['treeSettings']["computeGeneTrees"]["activate"] = False
        config["settings"]["selectionSettings"]["bayesCode"]["activate"] = False

    ## --------------------------------------------------------------------------
    ## Create DEF.yaml
    ## --------------------------------------------------------------------------

    ## Determine output file name and path
    yaml_name = 'DEF.yaml'

    ## The run directory is needed here regardless of -dr/-rs (DEF.yaml and the
    ## preprocessing outputs below have to go somewhere), and again later inside
    ## run_snakemake() if either flag is given -- _ensure_directory is the one
    ## place that guarantee is implemented, so the two can't drift apart.
    path_for_DEF = os.path.join(os.getcwd(), args.directory_name, yaml_name)
    path_to_snakemake_dir = _ensure_directory(os.path.join(os.getcwd(), args.directory_name))

    ## Load any DEF.yaml already sitting in this run directory once, up front --
    ## used both to decide whether preprocessing needs to rerun (below) and,
    ## further down, whether DEF.yaml itself needs rewriting.
    current_def_file = None
    if os.path.exists(path_for_DEF):
        current_yaml = YAML()
        with open(path_for_DEF, 'r') as file:
            current_def_file = current_yaml.load(file)

    ## preprocess input: Write new bed file of transcripts that are considered for selection screens, or else a link to the fasta files that
    ## are considered for screens/alignments
    ## Exchange input bed for filtered bed
    ## ISSUE: This will not work if the user passes a file called filtered.bed
    ##
    ## Preprocessing's output (filtered.bed, excluded_transcripts.txt) depends
    ## only on this specific subset of settings -- not aligner choice, HyPhy
    ## modes, resources, tree source, etc.
    preprocessing_inputs = {
        "TOGA_mode": config["togaMode"] or config["toga2Mode"],
        "free_mode": config["freeMode"],
        "foreground_list": args.foreground_lst,
        "fasta_path": str(config["fastaPath"]),
        "bed_file": str(config["originalTranscripts"]),
        "assembly_file": args.assembly_list,
        "min_species": args.minNumAlignedSpecies,
        "max_cds": args.maxCDSLengthReference,
        "min_cds": args.minCDSLengthReference,
    }
    ## Content fingerprints, so an edit-in-place (same path, different content --
    ## e.g. a species added to assemblies.txt, or TOGA rerun updating
    ## loss_summary.tsv) is caught too, not just a changed path/parameter.
    ## Kept separate from preprocessing_inputs itself, which is unpacked
    ## straight into run_preprocessing()'s call below.
    preprocessing_fingerprint = dict(preprocessing_inputs, **{
        "assembly_file_fp": preprocess_input.file_fingerprint(args.assembly_list),
        "foreground_list_fp": preprocess_input.file_fingerprint(args.foreground_lst),
        "bed_file_fp": preprocess_input.file_fingerprint(str(config["originalTranscripts"])),
        "fasta_fp": preprocess_input.fingerprint_fasta_path(
            str(config["fastaPath"]), preprocessing_inputs["TOGA_mode"], args.assembly_list),
    })
    config["preprocessingInputs"] = preprocessing_fingerprint

    filtered_bed_path = os.path.join(path_to_snakemake_dir, "filtered.bed")

    if (config["togaMode"] or config["toga2Mode"]):
        config["originalTranscripts"] = SingleQuotedScalarString(SingleQuotedScalarString(selected_transcripts_file))
        config["selectedTranscripts"] = SingleQuotedScalarString(filtered_bed_path)

    config["excludedTranscripts"] = SingleQuotedScalarString(os.path.join(path_to_snakemake_dir, "excluded_transcripts.txt"))

    ## Check whether the rest of DEF.yaml needs to change too. Done here, before
    ## preprocessing actually runs below, so a run that's about to be rejected
    ## (DEF.yaml differs, no -f) exits before paying for preprocessing instead
    ## of after.
    update_DEF = False
    write_def = True

    if current_def_file is not None:
        if current_def_file == config:
            logger.info("This command has been run before with identical settings; "
                        "DEF.yaml is unchanged.")
            write_def = False
        elif not args.force_run:
            logger.critical(
              "DEF.yaml already exists at %s and its content differs from your input commands. "
              "If you want to force a run in this working directory, enable -f/--force_run",
              path_for_DEF
            )
            sys.exit(1)
        else:
            logger.warning(
              "There is an existing snakemake selection directory at %s. "
              "The commands with which it was initialized differ from the currently invoked commands, "
              "but the run will be FORCED with the currently enabled commands.",
              path_for_DEF.split('DEF.yaml')[0]
            )
            update_DEF = True

    preprocessing_outputs = [os.path.join(path_to_snakemake_dir,"excluded_transcripts.txt")]

    if preprocessing_inputs["TOGA_mode"]:
        preprocessing_outputs.append(filtered_bed_path)

    needs_preprocessing = (
        not all(os.path.exists(f) for f in preprocessing_outputs)
        or current_def_file is None
        or current_def_file.get("preprocessingInputs") != preprocessing_fingerprint
    )
  
    if needs_preprocessing:
        preprocess_input.run_preprocessing(
              **preprocessing_inputs,
              output_log= os.path.join(path_to_snakemake_dir, "preprocess_report.log"),
              output_exclusion= os.path.join(path_to_snakemake_dir, "excluded_transcripts.txt"),
              output_filtered_bed = filtered_bed_path,
          )
    else:
        logger.info("Preprocessing inputs unchanged since the last run in this "
                    "directory; reusing the existing filtered.bed.")

    if write_def:
        with open(path_for_DEF, "w") as yaml_file:
            yaml = YAML()
            yaml.boolean_representation = ['False', 'True']
            yaml.default_flow_style = False
            yaml.indent(mapping=2, sequence=4, offset=2)
            yaml.preserve_quotes = True
            yaml.dump(config, yaml_file)
        logger.info("YAML configuration written to %s", yaml_name)

    logger.info("For an overview of all pipeline steps enabled, consult the DAG schematic at %s/rulegraph.pdf", args.directory_name)

    ## --------------------------------------------------------------------------
    ## Run snakemake
    ## --------------------------------------------------------------------------

    _log_run_plan(config)

    if config["settings"]["treeSettings"]["inputSpeciesTree"]["treeFile"]:
        snakefile = "Snakefile_speciesTree"
    else:
        snakefile = "Snakefile_standard"

    ## -dr dry-runs the upper-bound target itself, so snakemake prints the plan. 
    ## On -rs there is no plan printed at all, so log the upper bound before launching.
    if args.run_snakemake and not args.dry_run:
        _log_max_job_estimate(get_snakefile_path(snakefile),
                              os.path.join(os.getcwd(), args.directory_name))

    if args.dry_run or args.run_snakemake:
        returncode = run_snakemake(
            snake_file=snakefile,
            fasta_path=Fasta_path,
            directory_name=args.directory_name,
            dryrun=args.dry_run,
            snake_default=True,
            rerun_triggers=args.rerun_trigger,
            updated_rule=update_DEF,
            local_run=args.local_run,
            cores=args.cores,
            profile=args.profile,
            conda_prefix=args.conda_prefix,
            apptainer_prefix=args.apptainer_prefix,
            max_failed_fraction=args.max_failed_fraction,
        )
        ## Propagate failure. Previously the return code was never inspected, so
        ## `easel -rs` exited 0 even when the whole pipeline failed.
        sys.exit(returncode)
    else:
        logger.info(
          "DEF.yaml written but nothing was executed. Add -dr for a dry run, or "
          "-rs to launch the pipeline."
        )


if __name__ == "__main__":
    main()
