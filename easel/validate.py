#!/usr/bin/env python3
"""easel.validate - flag-combination and path checks.

Two jobs, both done *before* any config is built or any directory is created:

  1. every user-supplied path exists and is the right kind of thing
  2. every flag combination is one the pipeline can actually execute
"""

import logging
import os
import sys
from dataclasses import dataclass, field

from . import formats

logger = logging.getLogger("easel")

## Aligners whose output is fed through codonify_ali, which needs a real
## sequence name to codonify against (rules/codonify_alignment.smk).
CODONIFYING_ALIGNERS = {"prank", "prank_nt", "muscle"}

## Aligners whose output is already codon/in-frame nucleotide, so
## codonify_ali is skipped for them (rules/common.smk's
## get_uncleaned_alignment reads their _ori.fa directly).
CODON_AWARE_ALIGNERS = {"macse2", "prank_codon"}

## Which aligners each input mode has a rule chain for
## (Snakefile_standard include block).
## ISSUE: TOGA2 should be able to also run macse2/prank_codon
ALIGNERS_BY_MODE = {
    "toga":  {"prank", "muscle", "macse2"},
    "toga2": {"prank", "muscle", "macse2"},
    "free":  {"prank_nt", "prank_codon", "macse2", "muscle"},
}

HYPHY_MODES = {
    "absrel": {"std", "srv", "mh", "all"},
    "meme":   {"std", "srv", "mh", "all"},
    "relax":  {"std", "srv", "mh", "all"},
    "busted": {"std", "srv", "mh", "all", "error_sink"},
}

## HyPhy needs >=3 taxa; IQ-TREE bootstrapping wants >=4.
MIN_TAXA_HYPHY = 3
MIN_TAXA_IQTREE = 4


@dataclass
class Resolved:
    """Values derived during validation, so main() need not recompute them."""
    mode: str = ""
    fasta_path: str | None = None
    alignment_reference: str | None = None
    assemblies: list = field(default_factory=list)
    foreground: list = field(default_factory=list)
    tree_strategy: str = "none"          # input_tree | gene_trees | precomputed | none
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    @property
    def ok(self):
        return not self.errors


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def reference_from_toga_path(path, flag, r: Resolved):
    """Derive the reference name from a TOGA/TOGA2 path.

    Expected layout is <...>/<reference>/TOGA2, so the reference is the
    component above the last. Tolerates trailing slashes and relative paths;
    the naive path.split("/")[-2] returns 'TOGA2' for a trailing slash and
    '.' for './TOGA2'.
    """
    parts = [p for p in path.replace("\\", "/").split("/") if p not in ("", ".", "..")]
    if len(parts) < 2:
        r.errors.append(
            f"{flag}: cannot derive a reference name from {path!r}. Expected a "
            f"path ending in <reference>/TOGA2, or pass --reference_name explicitly.")
        return None
    return parts[-2]


def _read_list(path, label, r: Resolved, strip_vs=False):
    """formats.read_one_column, with its InputError turned into a collected error.

    This is the whole point of the split: formats does the parsing and the
    logging, validate decides whether to keep going.
    """
    try:
        return formats.read_one_column(path, label=label, strip_vs=strip_vs)
    except formats.InputError as exc:
        r.errors.append(str(exc))
        return []


def _report_pair(result, r: Resolved):
    """Fold a formats (ok, message) result into the collected errors."""
    ok, message = result
    if not ok:
        r.errors.append(message)
    return ok


def _content_check(func, r: Resolved, *args, **kwargs):
    """Call a formats function that raises InputError, collecting the failure."""
    try:
        func(*args, **kwargs)
        return True
    except formats.InputError as exc:
        r.errors.append(str(exc))
        return False


def _check_paths(args, r: Resolved):
    """Every path the user handed us must exist and be the right kind."""
    files = [
        (args.assembly_list,        "-asm/--assemblies"),
        (args.selected_transcripts, "-sb/--selected_bed_file"),
        (args.input_tree,           "-it/--input_tree"),
        (args.foreground_lst,       "--foreground_list"),
        (args.twoBitPath,           "--twoBit_path"),
        (getattr(args, "toga2Activate", None), "--toga2_activate"),
    ]
    dirs = [
        (args.toga_directory,   "-toga/--toga_reference_path"),
        (args.toga2_directory,  "-toga2/--toga2_reference_path"),
        (args.free_directory,   "-free/--free_mode"),
        (args.input_gene_trees, "--input_gene_trees"),
        (getattr(args, "profile", None), "--profile"),
    ]
    for value, flag in files:
        if not value:
            continue
        if not os.path.exists(value):
            r.errors.append(f"{flag}: no such file: {value}")
        elif os.path.isdir(value):
            r.errors.append(f"{flag}: expected a file, got a directory: {value}")
        elif os.path.getsize(value) == 0:
            r.errors.append(f"{flag}: file is empty: {value}")
    conda_prefix = getattr(args, "conda_prefix", None)
    if conda_prefix and not os.path.isdir(conda_prefix):
        parent = os.path.dirname(os.path.abspath(conda_prefix)) or "."
        if not os.path.isdir(parent):
            r.errors.append(f"--conda_prefix: parent directory does not exist: {parent}")
        elif not os.access(parent, os.W_OK):
            r.errors.append(f"--conda_prefix: parent directory is not writable: {parent}")

    for value, flag in dirs:
        if not value:
            continue
        if not os.path.exists(value):
            r.errors.append(f"{flag}: no such directory: {value}")
        elif not os.path.isdir(value):
            r.errors.append(f"{flag}: expected a directory, got a file: {value}")
        elif not os.listdir(value):
            r.errors.append(f"{flag}: directory is empty: {value}")


def _check_profile(args, r: Resolved):
    """A snakemake profile is a DIRECTORY holding config.yaml/config.v8+.yaml."""
    profile = getattr(args, "profile", None)
    if not profile or not os.path.isdir(profile):
        return                                   # existence handled in _check_paths
    if not any(os.path.isfile(os.path.join(profile, name))
               for name in ("config.yaml", "config.yml", "config.v8+.yaml")):
        r.errors.append(
            f"--profile {profile} contains no config.yaml; a snakemake profile "
            f"is a directory holding one")


def _check_mode_and_reference(args, r: Resolved):
    """Input mode, its mandatory companions, and the reference name."""
    picked = [(args.toga_directory, "toga"), (args.toga2_directory, "toga2"),
              (args.free_directory, "free")]
    picked = [(v, m) for v, m in picked if v]
    if len(picked) != 1:
        r.errors.append(
            "exactly one input mode is required: -toga, -toga2 or -free "
            f"({len(picked)} given)")
        return
    r.fasta_path, r.mode = picked[0]
    r.fasta_path = os.path.abspath(r.fasta_path)

    if r.mode in ("toga", "toga2"):
        flag = "-toga" if r.mode == "toga" else "-toga2"
        r.alignment_reference = reference_from_toga_path(
            picked[0][0], flag, r)
        if args.reference_name:
            r.alignment_reference = args.reference_name   # explicit override wins
        if not args.assembly_list:
            r.errors.append(f"{flag} requires -asm/--assemblies")
        if not args.selected_transcripts:
            r.errors.append(f"{flag} requires -sb/--selected_bed_file")
        if r.mode == "toga2" and not args.twoBitPath:
            r.errors.append("-toga2 requires --twoBit_path (TOGA2 alignment "
                            "extraction passes it to toga2 as -re)")
    else:  ## Free Mode
        r.alignment_reference = args.reference_name
        if args.selected_transcripts:
            r.errors.append(
                "-sb/--selected_bed_file is not supported with -free; easel runs "
                "every FASTA file in the directory.)")
        if r.alignment_reference is None and args.aligner in CODONIFYING_ALIGNERS:
            r.errors.append(
                f"--align {args.aligner} codonifies against a reference sequence, "
                f"but no --reference_name was given. Pass --reference_name NAME, "
                f"where NAME is a sequence header present in the input FASTA files.")
        elif r.alignment_reference is None:
            r.warnings.append(
                f"no --reference_name given; none is required for "
                f"--align {args.aligner}")
        ## Flags that only mean something for -toga2.
        for value, flag in ((args.twoBitPath, "--twoBit_path"),
                            (getattr(args, "toga2Activate", None), "--toga2_activate")):
            if value:
                r.errors.append(f"{flag} has no meaning with -free; it is only "
                                f"used by TOGA2 alignment extraction")
        if args.assembly_list:
            r.assemblies = _read_list(args.assembly_list, "assembly list", r)

    if r.mode in ("toga", "toga2") and args.assembly_list:
        r.assemblies = _read_list(args.assembly_list, "assembly list", r,
                                  strip_vs=True)
        if r.assemblies and r.fasta_path and os.path.isdir(r.fasta_path):
            _content_check(formats.check_toga_run_dirs, r,
                           r.fasta_path, r.assemblies)
        if args.twoBitPath and os.path.isfile(args.twoBitPath):
            _report_pair(formats.check_twobit(args.twoBitPath), r)
    if args.selected_transcripts and os.path.isfile(args.selected_transcripts):
        _report_pair(formats.check_bed_file(args.selected_transcripts), r)
    if r.mode == "free" and r.fasta_path and os.path.isdir(r.fasta_path):
        _report_pair(formats.validate_transcript_directory(r.fasta_path), r)


def _check_aligner(args, r: Resolved):
    if not r.mode:
        return
    allowed = ALIGNERS_BY_MODE[r.mode]
    if args.aligner and args.aligner not in allowed:
        r.errors.append(
            f"--align {args.aligner} has no rule chain in {r.mode} mode; "
            f"choose one of {sorted(allowed)}")
    elif args.aligner and args.aligner not in (CODONIFYING_ALIGNERS | CODON_AWARE_ALIGNERS):
        r.errors.append(
            f"--align {args.aligner} is classified in ALIGNERS_BY_MODE but not "
            f"as codonifying (CODONIFYING_ALIGNERS) or codon-aware "
            f"(CODON_AWARE_ALIGNERS) in easel.validate; add it to one of them.")
    if not args.aligner:
        if r.mode == "free" and not args.doScreenOnly:
            r.errors.append(
                "-free without --align: set an aligner, or --do_screen_only if the "
                "input FASTA files are already aligned")
        if r.mode in ("toga", "toga2"):
            r.errors.append(
                f"-{r.mode} requires --align (TOGA output is unaligned). To screen "
                f"existing alignments, pass them with -free instead.")
    if args.aligner and args.doScreenOnly:
        r.warnings.append(
            f"--align {args.aligner} is ignored because --do_screen_only is set")


def _check_hyphy_modes(args, r: Resolved):
    for value, prog in ((args.absrel, "absrel"), (args.busted, "busted"),
                        (args.meme, "meme"), (args.relax, "relax")):
        if not value:
            continue
        tokens = [t.strip() for t in value.split(",") if t.strip()]
        if not tokens:
            r.errors.append(f"--{prog}: empty mode list")
            continue
        allowed = HYPHY_MODES[prog]
        ints = []
        for t in tokens:
            if t in allowed:
                continue
            if prog == "relax" and t.isdigit():
                ints.append(int(t))
                continue
            r.errors.append(f"--{prog}: unknown mode {t!r} "
                            f"(allowed: {sorted(allowed)}"
                            f"{' plus a round count' if prog == 'relax' else ''})")
        if "std" in tokens and ({"srv", "mh", "all"} & set(tokens)):
            r.errors.append(f"--{prog}: 'std' cannot be combined with srv/mh/all")
        if len(ints) > 1:
            r.errors.append(f"--relax: give at most one round count, got {ints}")
        if ints and ints[0] < 1:
            r.errors.append("--relax: round count must be >= 1 "
                            "(0 gives an empty aggregate)")
    if args.relax and not args.foreground_lst:
        r.errors.append("--relax requires --foreground_list "
                        "(RELAX is run with --test Foreground)")
    if args.doErrorCleaning and not args.busted:
        r.warnings.append("--do_error_sink_cleaning enables BUSTED with error_sink, "
                          "which was not requested explicitly")
    elif args.doErrorCleaning and "error_sink" not in (args.busted or ""):
        r.warnings.append("--do_error_sink_cleaning adds error_sink to --busted")


def _check_trees(args, r: Resolved):
    chosen = [n for n, v in (("input_tree", args.input_tree),
                             ("gene_trees", args.tree),
                             ("precomputed", args.input_gene_trees)) if v]
    if len(chosen) > 1:
        r.errors.append(
            f"choose one tree source, got {len(chosen)}: "
            f"{', '.join(chosen)} (-it / -ct / --input_gene_trees)")
    elif not chosen and not args.doAlignmentOnly:
        r.errors.append(
            "no tree source: pass -it/--input_tree, -ct/--comp_tree, or "
            "--input_gene_trees (or --do_alignment_only to stop after alignment)")
    ## In what case would this be none?        
    r.tree_strategy = chosen[0] if len(chosen) == 1 else "none"

    ## ISSUE: This should be possible (?)
    if args.input_gene_trees and r.mode in ("toga", "toga2"):
        r.errors.append("--input_gene_trees is not supported with -toga/-toga2; "
                        "pass pre-aligned sequences with -free, or use -ct / -it")
    if chosen and args.doAlignmentOnly:
        r.warnings.append(f"--do_alignment_only: the {chosen[0]} tree source "
                          f"will not be used")
    # ## ISSUE: This has to be tested/fixed
    # if args.bayescode and args.foreground_lst and r.tree_strategy == "precomputed":
    #     r.errors.append(
    #         "--bayescode with --foreground_list and --input_gene_trees requests "
    #         "tmp/<id>_tmp.treefile, which only rule compute_tree produces (-ct). "
    #         "Use -ct, or drop --foreground_list for the BayesCode run.")

    if args.input_tree and os.path.isfile(args.input_tree):
        try:
            tree = formats.read_in_tree(args.input_tree)
        except formats.InputError as exc:
            r.errors.append(str(exc))
        else:
            if r.assemblies:
                _report_pair(
                    formats.compare_tree_asm(tree, r.assemblies, r.alignment_reference),
                    r)
            else:
                r.warnings.append(
                    "no assemblies given, but a species input tree was provided")

    if args.input_gene_trees and os.path.isdir(args.input_gene_trees):
        _report_pair(formats.validate_tree_directory(args.input_gene_trees), r)


def _check_stages(args, r: Resolved):
    if args.doAlignmentOnly and args.doScreenOnly:
        r.errors.append("--do_alignment_only and --do_screen_only are mutually "
                        "exclusive")
    if args.doScreenOnly and r.mode in ("toga", "toga2"):
        r.errors.append("--do_screen_only cannot be used with -toga/-toga2; pass "
                        "the existing alignments with -free instead")
    analyses = [n for n, v in (("--absrel", args.absrel), ("--busted", args.busted),
                               ("--meme", args.meme), ("--relax", args.relax),
                               ("--bayescode", args.bayescode)) if v]
    if args.doAlignmentOnly and analyses:
        r.warnings.append("--do_alignment_only: disabling " + ", ".join(analyses))
    if not args.doAlignmentOnly and not analyses:
        r.warnings.append("no analysis enabled; the run will stop after alignment "
                          "postprocessing (did you mean --do_alignment_only?)")
    
    ## ISSUE: Must be implemented still
    if args.species_tree:
        r.errors.append("-st/--comp_species_tree is not implemented yet "
                        "(rules/run_astral.smk does not exist)")


def _check_numbers_and_lists(args, r: Resolved):
    if args.minCDSLengthReference >= args.maxCDSLengthReference:
        r.errors.append(f"--min_CDS_length ({args.minCDSLengthReference}) must be "
                        f"below --max_CDS_length ({args.maxCDSLengthReference})")
    floor = MIN_TAXA_IQTREE if args.tree else MIN_TAXA_HYPHY
    if args.minNumAlignedSpecies < floor:
        r.errors.append(
            f"--min_num_aligned_species must be at least {floor} "
            f"({'IQ-TREE bootstrapping' if args.tree else 'HyPhy'} needs it)")
    for name, val in vars(args).items():
        if (name.endswith(("_threads", "_mem_mb")) and isinstance(val, int)
                and val <= 0):
            r.errors.append(f"--{name} must be positive, got {val}")
    ## Shape (one name per line, no duplicates) is checked whenever the file is
    ## given, regardless of whether -asm was also given: -asm is optional in
    ## free mode, so a foreground list without it is a normal RELAX-only run.
    if args.foreground_lst:
        r.foreground = _read_list(args.foreground_lst, "foreground list", r)
        if args.assembly_list and r.assemblies:
            stray = [s for s in r.foreground if s.replace("vs_", "") not in r.assemblies]
            if stray:
                r.errors.append(
                    f"--foreground_list contains {len(stray)} name(s) absent from "
                    f"--assemblies: {', '.join(stray[:5])}"
                    f"{' …' if len(stray) > 5 else ''}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def validate(args) -> Resolved:
    """Run every check. Order matters only in that mode is resolved first."""
    r = Resolved()
    _check_paths(args, r)
    _check_profile(args, r)
    _check_mode_and_reference(args, r)
    _check_aligner(args, r)
    _check_hyphy_modes(args, r)
    _check_trees(args, r)
    _check_stages(args, r)
    _check_numbers_and_lists(args, r)
    return r


def report(r: Resolved, exit_on_error=True):
    for w in r.warnings:
        logger.warning(w)
    for e in r.errors:
        logger.critical(e)
    if r.errors:
        logger.critical("%d problem(s) found; nothing was written.", len(r.errors))
        if exit_on_error:
            sys.exit(1)
    else:
        logger.info("input OK  mode=%s  reference=%s  tree=%s  assemblies=%d",
                    r.mode, r.alignment_reference, r.tree_strategy, len(r.assemblies))
    return r.ok
