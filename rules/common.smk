## Helpers shared by Snakefile_standard and Snakefile_speciesTree.

## Included with `include: "rules/common.smk"`

import glob
import os
import re
import sys
from pathlib import Path

from easel.validate import CODONIFYING_ALIGNERS
from easel.formats import FASTA_EXTENSIONS, TREE_EXTENSIONS, check_bed_file

#FASTA_EXTENSIONS = ('.fa', '.fasta', '.fna', '.fas')
#TREE_EXTENSIONS = ('.nwk', '.nh', '.newick', '.tree', '.tre', '.nhx', '.treefile')

def _split_extension(filename, extensions):
    """Return (stem, suffix) for a recognised extension, else (None, None).

    Handles a single trailing '.gz'.
    """
    name, gz = filename, False
    if name.endswith('.gz'):
        name, gz = name[:-3], True
    for ext in extensions:
        if name.endswith(ext):
            return name[:-len(ext)], ext + ('.gz' if gz else '')
    return None, None


def read_bed_names(bed_path):
    """Transcript names from column 4 of a BED file, in file order.
    """
    ok, message = check_bed_file(bed_path)
    if not ok:
        sys.exit(f"ERROR: {message}")

    names = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            names.append(line.split("\t")[3])
    return names


def read_transcript_names_from_free(input_dir, excluded_transcripts):
    """Transcript names from the FASTA files in *input_dir*, minus exclusions.

    Returns (names, suffix). *excluded_transcripts* is the path to the exclusion
    list written by preprocess_input.
    """
    exclude = set()
    if excluded_transcripts and os.path.exists(excluded_transcripts):
        with open(excluded_transcripts) as fh:
            for line in fh:
                name = line.strip()
                if name:
                    ## The list may carry file names; compare on the stem too.
                    exclude.add(name)
                    stem, _ = _split_extension(name, FASTA_EXTENSIONS)
                    if stem:
                        exclude.add(stem)

    names, suffixes = [], set()
    for file_path in sorted(Path(input_dir).iterdir()):
        if not file_path.is_file() or file_path.name.startswith("."):
            continue
        stem, suffix = _split_extension(file_path.name, FASTA_EXTENSIONS)
        if stem is None:
            continue
        suffixes.add(suffix)
        if stem not in exclude:
            names.append(stem)

    if not names:
        sys.exit(
            f"ERROR: no usable FASTA files in {input_dir} "
            f"(after removing {len(exclude)} excluded entries)."
        )
    if len(suffixes) > 1:
        sys.exit(
            f"ERROR: {input_dir} mixes FASTA extensions {sorted(suffixes)}. "
            f"The rules build input paths from a single suffix, so please use "
            f"one extension throughout."
        )
    return names, suffixes.pop()


def read_gene_tree_names(input_dir, transcript_list):
    """Match precomputed gene trees to transcripts. Returns (matched, suffix)."""
    wanted = set(transcript_list)
    matches, suffixes, unmatched = [], set(), []

    for file_path in sorted(Path(input_dir).iterdir()):
        if not file_path.is_file() or file_path.name.startswith("."):
            continue
        stem, ext = _split_extension(file_path.name, TREE_EXTENSIONS)
        if stem is None:
            continue
        if stem in wanted:  
            matches.append(stem)
            suffixes.add(file_path.name[len(stem):])
            continue
        ## Only reached when the stem is not itself a transcript, for instance
        ## tree named <transcript>.gene_tree.nh. 
        candidates = [file_path.name[:i] for i in range(1, len(file_path.name))
                      if file_path.name[:i] in wanted]
        if len(candidates) == 1:
            matches.append(candidates[0])
            suffixes.add(file_path.name[len(candidates[0]):])
        elif len(candidates) > 1:
            sys.exit(
                f"ERROR: {file_path.name} matches more than one transcript "
                f"({', '.join(sorted(candidates)[:5])}). Name gene trees "
                f"<transcript><ext> so the mapping is unambiguous."
            )
        else:
            unmatched.append(file_path.name)

    if unmatched:
        shown = ", ".join(unmatched[:5]) + (" …" if len(unmatched) > 5 else "")
        sys.exit(
            f"ERROR: {len(unmatched)} tree file(s) in {input_dir} do not "
            f"correspond to any input transcript: {shown}"
        )
    if len(suffixes) > 1:
        sys.exit(
            f"ERROR: {input_dir} mixes gene-tree suffixes {sorted(suffixes)}; "
            f"the rules build input paths from a single suffix."
        )
    missing = wanted - set(matches)
    if missing:
        shown = ", ".join(sorted(missing)[:5]) + (" …" if len(missing) > 5 else "")
        sys.exit(
            f"ERROR: no gene tree found for {len(missing)} transcript(s): "
            f"{shown}. Provide one tree per transcript."
        )
    return matches, (suffixes.pop() if suffixes else "")


def write_species_list(toga_path, species_list, out_path="species.TOGA.dir.txt"):
    """One TOGA run directory per line, for the alignment-extraction rule."""
    with open(out_path, "w") as fh:
        for entry in species_list:
            name = entry if entry.startswith("vs_") else "vs_" + entry
            fh.write(os.path.join(toga_path, name) + "\n")


def transcript_wildcard_pattern(transcripts):
    """The wildcard_constraints pattern for {transcript_id}: one path segment.

    NOT an alternation of the transcript names, which is what this used to
    return ("|".join(re.escape(t) for t in transcripts)).

    One path segment loses nothing here. Every output template in this
    workflow carries {transcript_id} as its own directory component
    (codon_alignments/{transcript_id}/...), and snakemake requires all
    occurrences of a wildcard within one pattern to take the same value. The
    directory therefore pins the id, and '[^/]+' cannot cross the '/' that
    delimits it.

    *transcripts* is still taken, and still checked, for the one thing that
    would break the argument above: a name containing a path separator.
    """
    bad = [t for t in transcripts if "/" in t or chr(92) in t]
    if bad:
        sys.exit(
            "ERROR: transcript name(s) contain a path separator and cannot be "
            "used as a directory name: " + ", ".join(bad[:5])
        )
    return "[^/]+"


## ---------------------------------------------------------------------------
## Parse-time run setup, shared by both Snakefiles
## ---------------------------------------------------------------------------

## ---------------------------------------------------------------------------
## Transcript list cache
## ---------------------------------------------------------------------------
## In cluster mode snakemake spawns a fresh process per job (spawn_jobs.py) and
## every one of them executes this file top to bottom before running its own
## rule. So whatever resolve_run_inputs() does costs once per JOB, not once per
## run -- a validate_alignment job that reads a single FASTA was re-validating
## the entire BED first.

TRANSCRIPT_CACHE = ".easel_transcripts"

def _source_stamp(*paths):
    """Identity of the inputs the cache was built from.

    For a file, size and mtime. For a directory, mtime only -- POSIX updates it
    when an entry is added, removed or renamed, which is exactly what changes a
    -free transcript list. Counting entries would mean listing the directory,
    i.e. paying the cost this cache exists to avoid.
    """
    parts = []
    for path in paths:
        if not path:
            parts.append("-")
            continue
        try:
            st = os.stat(path)
        except OSError:
            parts.append("missing")
            continue
        parts.append(f"{'d' if os.path.isdir(path) else 'f'}"
                     f"{0 if os.path.isdir(path) else st.st_size}:{st.st_mtime_ns}")
    return "|".join(parts)


def _read_transcript_cache(*sources):
    """-> (transcripts, suffix), or None when absent, stale or unreadable.

    Any failure returns None, which just means doing the work properly. There
    is no failure mode here that should stop a run.
    """
    try:
        with open(TRANSCRIPT_CACHE, encoding="utf-8") as fh:
            if fh.readline().rstrip("\n") != _source_stamp(*sources):
                return None
            suffix = fh.readline().rstrip("\n")
            names = [ln.rstrip("\n") for ln in fh if ln.strip()]
    except OSError:
        return None
    if not names:
        return None
    return names, suffix


def _write_transcript_cache(sources, transcripts, suffix=""):
    tmp = f"{TRANSCRIPT_CACHE}.{os.getpid()}.tmp"
    try:
        with open(tmp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(_source_stamp(*sources) + "\n")
            fh.write((suffix or "") + "\n")
            for name in transcripts:
                fh.write(name + "\n")
        os.replace(tmp, TRANSCRIPT_CACHE)
    except OSError:
        try:
            os.remove(tmp)
        except OSError:
            pass


def resolve_run_inputs():
    """Resolve the transcript list, match precomputed gene trees to it if
    given, write the TOGA/TOGA2 species-directory list, and fill in
    config["fileSuffix"]/config["treeSuffix"] from whatever is actually on
    disk. Returns the transcript list.
    """
    all_species_lst = config["speciesLst"] if config["speciesLst"] else []

    if config.get("fastaPath") and not str(config["fastaPath"]).endswith(("/", "\\")):
        config["fastaPath"] = str(config["fastaPath"]) + "/"

    ## Read from the cache when the sources are unchanged, otherwise do it
    ## properly and leave a cache behind for the jobs this run will spawn.
    if config["selectedTranscripts"] is not None:
        sources = (config["selectedTranscripts"],)
        cached = _read_transcript_cache(*sources)
        if cached is None:
            transcripts = read_bed_names(config["selectedTranscripts"])
            _write_transcript_cache(sources, transcripts)
        else:
            transcripts, _ = cached
    else:
        ## The exclusion list is a source too: it changes which FASTAs count,
        ## without touching the FASTA directory itself.
        sources = (config["fastaPath"], config["excludedTranscripts"])
        cached = _read_transcript_cache(*sources)
        if cached is None:
            transcripts, suffix = read_transcript_names_from_free(
                config["fastaPath"], config["excludedTranscripts"])
            _write_transcript_cache(sources, transcripts, suffix)
        else:
            ## The suffix is cached with the list because it is derived from the
            ## same directory listing and mutates config here -- skipping the
            ## listing without it would leave every input path malformed.
            transcripts, suffix = cached
        config["fileSuffix"] = suffix

    if config["settings"]["treeSettings"]["preCompGeneTrees"]["geneTreePath"]:
        tree_matches, suffix = read_gene_tree_names(
            config["settings"]["treeSettings"]["preCompGeneTrees"]["geneTreePath"], transcripts)
        config["treeSuffix"] = suffix
        if len(tree_matches) != len(transcripts):
            sys.exit(
                "Number of extracted gene trees and transcripts doesn't match. "
                "When providing gene trees, make sure to provide one tree per transcript"
            )

    if config["togaMode"] or config["toga2Mode"]:
        write_species_list(config["fastaPath"], all_species_lst)

    sample = config.get("graphSample")
    if sample:
        total = len(transcripts)
        n = max(1, min(int(sample), total))
        transcripts = transcripts[:n]
        msg = f"easel: graphSample={n} of {total} transcripts"
        print(msg, file=sys.stderr)

    return transcripts

## ---------------------------------------------------------------------------
## Two-phase run: parse-time verdicts instead of a checkpoint
## ---------------------------------------------------------------------------
##
## The Gating is now done with two snakemake invocations, driven by
## easel/cli.py:
##
##   phase 1   --config easel_phase=validate   target: rule validated
##       Align, clean and validate every transcript. No analysis rule is in
##       this DAG at all. Run with --keep-going, so one aligner that dies does
##       not stop the rest
##
##   phase 2   --config easel_phase=analyse    target: rule all
##       Re-parse. The verdicts are on disk now, so the transcript list is
##       filtered HERE, at parse time, and a rejected transcript is simply not
##       in the DAG. Nothing is decided during the run, so nothing re-plans.
PHASE_VALIDATE = "validate"
PHASE_ANALYSE = "analyse"

## One line per transcript: transcript<TAB>verdict<TAB>reason. Written at the
## end of phase 1 (write_verdict_summary), read at parse time in phase 2.
VERDICT_SUMMARY = "validation_summary.tsv"

## OK and SKIP are written by scripts/check_alignment.py. FAIL is never written
## by anything: it is inferred here from the ABSENCE of a verdict.
VERDICT_OK = "OK"
VERDICT_SKIP = "SKIP"
VERDICT_FAIL = "FAIL"


def easel_phase():
    """Which of the two passes this parse belongs to.

    Defaults to "analyse", so a hand-run `snakemake` still targets rule all
    """
    return str(config.get("easel_phase", PHASE_ANALYSE))


def verdict_path(transcript_id):
    return f"codon_alignments/{transcript_id}/validation.txt"


def read_verdict(path):
    """-> (verdict, reason), or (None, reason) when there is no usable verdict.

    Deliberately total. A truncated or unreadable verdict must resolve to "not
    usable", not raise inside an input function -- where the traceback would
    name a rule that has nothing to do with it.
    """
    try:
        with open(path, encoding="utf-8") as fh:
            fields = fh.read().split("\t")
    except OSError:
        return None, ""
    verdict = fields[0].strip()
    reason = fields[1].strip() if len(fields) > 1 else ""
    if verdict not in (VERDICT_OK, VERDICT_SKIP):
        return None, f"unusable verdict {verdict[:40]!r} in {path}"
    return verdict, reason


def write_verdict_summary(transcript_ids):
    """Tally every transcript's verdict into VERDICT_SUMMARY. Returns counts.

    Called from phase 1's onsuccess/onerror hooks, not from a rule. It cannot
    be a rule: under --keep-going a failed aligner means that transcript's
    validation.txt is never produced, so a rule declaring them as inputs could
    never run -- and those are exactly the transcripts this file exists to
    record.

    Rewritten from scratch at the end of every phase 1, never merged with an
    older copy. That is what keeps FAIL from becoming permanent: a transcript
    that died on a bad node is attempted again on the next run instead of being
    excluded forever.
    """
    counts = {VERDICT_OK: 0, VERDICT_SKIP: 0, VERDICT_FAIL: 0}
    tmp = f"{VERDICT_SUMMARY}.{os.getpid()}.tmp"
    with open(tmp, "w", encoding="utf-8", newline="\n") as out:
        out.write("transcript\tverdict\treason\n")
        for tid in transcript_ids:
            verdict, reason = read_verdict(verdict_path(tid))
            if verdict is None:
                verdict = VERDICT_FAIL
                reason = reason or (
                    "no alignment verdict: alignment, cleaning or validation "
                    f"did not complete; see codon_alignments/{tid}/")
            counts[verdict] += 1
            reason = reason.replace("\t", " ").replace("\n", " ")
            out.write(f"{tid}\t{verdict}\t{reason}\n")
    os.replace(tmp, VERDICT_SUMMARY)
    return counts


def _read_verdict_summary():
    """-> {transcript: verdict}, or None when absent, empty or unrecognisable."""
    rows = {}
    try:
        with open(VERDICT_SUMMARY, encoding="utf-8") as fh:
            if not fh.readline().startswith("transcript\t"):
                return None
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[0]:
                    rows[fields[0]] = fields[1]
    except OSError:
        return None
    return rows or None


def verdict_summary_rows():
    """Yield (transcript, verdict, reason) from VERDICT_SUMMARY; nothing if absent.

    skip_report needs the reasons, not just the verdicts, and it needs the FAIL
    rows in particular: a transcript with no alignment has no done.txt either,
    so this file is the only place it appears.
    """
    try:
        with open(VERDICT_SUMMARY, encoding="utf-8") as fh:
            if not fh.readline().startswith("transcript\t"):
                return
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[0]:
                    yield fields[0], fields[1], fields[2] if len(fields) > 2 else ""
    except OSError:
        return


_VERDICT_SETS = None
_USABLE = frozenset()


def resolve_verdicts(transcript_ids):
    """Split the transcript list three ways, at parse time.

    Returns (usable, attempted, failed):

      usable     verdict OK. Trees, selection screens and summary.pdf are
                 requested for these and only these.
      attempted  usable + SKIP: everything that HAS an alignment. Archiving,
                 cleanup and the cleaning report run for all of them -- a
                 rejected alignment still gets tarred and its tmp/ deleted,
                 which is the entire point of the storage work. Filtering these
                 down to `usable` would leave a rejected transcript tmp/ on
                 disk forever.
      failed     no usable verdict at all. These get nothing but a row in
                 skipped_transcripts.tsv, and their tmp/ is deliberately left
                 in place: a failure is the one case you want to inspect.

    Nothing is filtered in phase 1, nor for the describe-only DAG builds
    (assumeAllUsable). Phase 1 in particular MUST attempt every transcript, or
    one that failed on the previous run would never be retried.
    """
    global _VERDICT_SETS, _USABLE
    if _VERDICT_SETS is not None:
        return _VERDICT_SETS

    unfiltered = (list(transcript_ids), list(transcript_ids), [])
    if easel_phase() == PHASE_VALIDATE or config.get("assumeAllUsable"):
        _VERDICT_SETS = unfiltered
    else:
        rows = _read_verdict_summary()
        if rows is None:
            print(f"easel: no {VERDICT_SUMMARY} in this run directory; "
                  "assuming every transcript is usable. Run phase 1 first "
                  f"(snakemake --config easel_phase={PHASE_VALIDATE} validated) "
                  "to gate the analyses on alignment validation.",
                  file=sys.stderr)
            _VERDICT_SETS = unfiltered
        else:
            usable, attempted, failed = [], [], []
            for tid in transcript_ids:
                verdict = rows.get(tid)
                if verdict == VERDICT_OK:
                    usable.append(tid)
                    attempted.append(tid)
                elif verdict == VERDICT_SKIP:
                    attempted.append(tid)
                else:
                    failed.append(tid)
            _VERDICT_SETS = (usable, attempted, failed)
    _USABLE = frozenset(_VERDICT_SETS[0])
    return _VERDICT_SETS


def finish_validation_phase():
    """End of phase 1: write the verdict summary and say what is in it.

    Called from BOTH hooks. onsuccess alone is not enough: phase 1 runs with
    --keep-going, so a run where some aligners died ends in onerror, and those
    are exactly the runs whose FAIL rows matter.
    """
    counts = write_verdict_summary(transcripts)
    total = sum(counts.values())
    print(f"\nValidation of {total} transcript(s): {counts[VERDICT_OK]} usable, "
          f"{counts[VERDICT_SKIP]} rejected by validation, "
          f"{counts[VERDICT_FAIL]} with no alignment.")
    print(f"  verdicts written to {VERDICT_SUMMARY}\n")
    return counts


def is_usable(transcript_id):
    """True when this transcript's alignment passed validation.
    """
    if _VERDICT_SETS is None:
        resolve_verdicts(transcripts)
    return transcript_id in _USABLE

## Helper functions for HyPhy output file generation to rule all. Both
## Snakefiles call these from their own gather_input() with a different
## "keyword" (_g_tree/.g_tree for the standard tree, _i_tree/.i_tree for the
## species tree).

def gather_absrel(keyword, ids=None):
    ## Which transcripts to request for. rule all passes the transcripts
    ## that PASSED validation; rulegraph_reference passes none, meaning
    ## all of them, because it exists to draw the workflow in full.
    ids = transcripts if ids is None else ids
    absrel_input = []
    if config["settings"]["selectionSettings"]["ABSREL"]["activate"]:
        absrel_input.extend(
                expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.absrel_tree{insert}.nh",  transcript_id=ids, insert = keyword),
            ),
        absrel_input.extend(
                expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.absrel_er{insert}.tsv",  transcript_id=ids, insert = keyword),
            ),

    return absrel_input


def gather_meme(keyword, ids=None):
    ## Which transcripts to request for. rule all passes the transcripts
    ## that PASSED validation; rulegraph_reference passes none, meaning
    ## all of them, because it exists to draw the workflow in full.
    ids = transcripts if ids is None else ids
    meme_input = []
    if config["settings"]["selectionSettings"]["MEME"]["activate"]:
        meme_input.extend(
            expand( "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.meme_mle{insert}.tsv", transcript_id=ids, insert = keyword)
            ),
        meme_input.extend(
            expand( "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.meme_er{insert}.tsv", transcript_id=ids, insert = keyword)
            ),
    return meme_input


def gather_relax(keyword, ids=None):
    ## Which transcripts to request for. rule all passes the transcripts
    ## that PASSED validation; rulegraph_reference passes none, meaning
    ## all of them, because it exists to draw the workflow in full.
    ids = transcripts if ids is None else ids
    relax_input = []
    if config["settings"]["selectionSettings"]["RELAX"]["activate"]:
        relax_input.extend(
            expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.relax{insert}.tsv", transcript_id = ids, insert = keyword)),
    return relax_input


def gather_bayescode(keyword, ids=None):
    ## Which transcripts to request for. rule all passes the transcripts
    ## that PASSED validation; rulegraph_reference passes none, meaning
    ## all of them, because it exists to draw the workflow in full.
    ids = transcripts if ids is None else ids
    bayescode_input = []
    if config["settings"]["selectionSettings"]["bayesCode"]["activate"]:
        bayescode_input.extend(
            expand("codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel{insert}.ci0.025.tsv", transcript_id=ids, insert = keyword))
        bayescode_input.extend(
            expand("codon_alignments/{transcript_id}/tmp/{transcript_id}_classic{insert}.ci0.025.tsv", transcript_id=ids, insert = keyword))

    return bayescode_input


def gather_busted(keyword, ids=None):
    ## Which transcripts to request for. rule all passes the transcripts
    ## that PASSED validation; rulegraph_reference passes none, meaning
    ## all of them, because it exists to draw the workflow in full.
    ids = transcripts if ids is None else ids
    busted_input = []
    if config["settings"]["selectionSettings"]["BUSTED"]["activate"] and not config["settings"]["selectionSettings"]["BUSTED"]['error_sink']:
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/{transcript_id}.busted{insert}.json", transcript_id=ids, insert = keyword)
            ),
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_er{insert}.tsv", transcript_id=ids, insert = keyword)
            ),
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_model{insert}.tsv", transcript_id=ids, insert = keyword)
            ),
    elif config["settings"]["selectionSettings"]["BUSTED"]['error_sink'] and not config["settings"]["selectionSettings"]["doErrorCleaning"]:
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/{transcript_id}.busted{insert}.filtered.json",  transcript_id=ids, insert = keyword),
        ),
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_er{insert}.tsv", transcript_id=ids, insert = keyword)
            ),
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_model{insert}.tsv", transcript_id=ids, insert = keyword)
            ),
    elif config["settings"]["selectionSettings"]["BUSTED"]['error_sink'] and config["settings"]["selectionSettings"]["doErrorCleaning"]:
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/tmp/{transcript_id}.busted{insert}.filtered.fa", transcript_id=ids, insert = keyword),
        ),
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_er{insert}.tsv", transcript_id=ids, insert = keyword)
            ),
        busted_input.extend(
            expand("codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_model{insert}.tsv", transcript_id=ids, insert = keyword)
            ),

    return busted_input


## Every file the activated analyses produce, for ONE keyword scheme
## ("_g_tree"/".g_tree" for the standard tree, "_i_tree"/".i_tree" for the
## species tree).
##
## For rule rulegraph_reference ONLY. rule all must never request these: doing
## so is exactly what previously let aBSREL/BUSTED/IQ-TREE run on transcripts
## the validate_alignment gate had rejected.
def analysis_output_files(absrel_keyword, hyphy_keyword):
    return (list(gather_absrel(absrel_keyword))
            + list(gather_meme(hyphy_keyword))
            + list(gather_relax(hyphy_keyword))
            + list(gather_busted(hyphy_keyword))
            + expand("codon_alignments/{transcript_id}/{transcript_id}.summary.pdf",
                     transcript_id=transcripts)
            + expand("codon_alignments/{transcript_id}/validation.txt",
                     transcript_id=transcripts)
            + expand("codon_alignments/{transcript_id}/{transcript_id}.cleaning.txt",
                     transcript_id=transcripts)
            + expand("codon_alignments/{transcript_id}/{transcript_id}.final.fa",
                     transcript_id=transcripts))


## ---------------------------------------------------------------------------
## Input functions shared across multiple rule files
## ---------------------------------------------------------------------------
## Each of these is referenced from a different file than the one that
## produces the alignment it points at (rules/codonify_alignment.smk,
## rules/hmm_cleaning.smk, rules/manual_cleaning.smk, rules/run_absrel.smk,
## rules/run_busted.smk, rules/run_meme.smk, rules/run_relax_cp.smk,
## rules/run_bayescode.smk, rules/run_iqtree.smk, rules/validate_alignment.smk,
## rules/draw_tree.smk, rules/draw_output.smk, rules/precomp_gene_tree.smk),
## so they live here rather than depending on include: order to make them
## visible where needed.

def get_uncleaned_alignment(wildcards):
    """The un-cleaned nucleotide alignment: codonify_ali's output for
    aligners that need codonification, or the aligner's direct output
    otherwise (already codon/in-frame, e.g. macse2/prank_codon -- see
    easel.validate.CODON_AWARE_ALIGNERS)."""
    if config["settings"]["alignmentSettings"]["aligner"] in CODONIFYING_ALIGNERS:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.fa"
    else:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa"


def get_input_for_tree(wildcards):
    """The alignment to hand to the tree/HyPhy/BayesCode steps, once
    whichever cleaning steps are active have run."""
    if config["settings"]["cleaningSettings"]["manualCleaning"]["doManualCleaning"]:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}.manual.fa"
    elif config["settings"]["cleaningSettings"]["hmmCleaning"]["doHMMCleaning"]:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.hmm_cleaned.fa"
    else:
        return get_uncleaned_alignment(wildcards)

## ISSUE: Rewrite this in rules
get_input_for_busted = get_input_for_tree

def get_input_tree_for_hyphy(wildcards):
    """Which tree file the HyPhy/BayesCode/BUSTED steps should use."""
    if config['settings']['treeSettings']['computeGeneTrees']['activate']:
        return "codon_alignments/{transcript_id}/{transcript_id}_iqtree.nh"
    else:
        return "codon_alignments/{transcript_id}/{transcript_id}_pruned_tree.nh"


def get_input_ali_for_hyphy(wildcards):
    """The alignment HyPhy/draw_output should use, after optional
    error-sink-based filtering."""
    if config["settings"]["selectionSettings"]["BUSTED"]['error_sink'] and config["settings"]["selectionSettings"]["doErrorCleaning"]:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}.busted.g_tree.filtered.fa"
    else:
        return get_input_for_busted(wildcards)


## ---------------------------------------------------------------------------
## Per-transcript combined log
## ---------------------------------------------------------------------------

## The pipeline stage each rule's log belongs to, in the order a transcript
## actually moves through them in a normal run.
RULE_LOG_ORDER = [
    ## alignment extraction
    "extract_ali", "extract_ali_macse2",
    "extract_ali_muscle_afa", "extract_ali_muscle_best", "extract_ali_muscle_cc",
    "extract_ali_prank", "rename_prank",
    ## codonification
    "rename_reference", "codonify_ali", "translate_prank_aa",
    ## cleaning
    "convert_to_aa", "hmm_cleaner", "transfer_cleaner", "clean_transfer_cleaner",
    "manual_cleaner",
    ## validation
    "validate_alignment",
    ## tree building
    "compute_tree", "gene_tree_annotation",
    "copy_precomp_gene_tree", "precomp_tree_bayescode",
    "draw_tree", "draw_tree_bayescode",
    ## selection screens
    "absrel_run",
    "busted_run", "error_mask", "clean_error_sink_ali",
    "meme_run",
    "relax_run", "aggregate_relax",
    ## turning HyPhy JSON into the TSVs draw_output reads
    "extract_absrel", "extract_busted", "extract_meme", "extract_relax",
    ## bayescode
    "convert_to_phylip", "mutselomega_w0", "mutselomega_w",
    "read_mutselomega_w0", "read_mutselomega_w",
    ## final output
    "cleaning_report", "draw_output", "keep_final_alignment",
    "transcript_done", "selection_report",
    ## archive-then-delete, twice: tmp/ first, then HyPhy_output/ + *.json
    "tar_compress", "clean_tmp", "tar_results", "clean_results",
]
_RULE_LOG_ORDER_INDEX = {name: i for i, name in enumerate(RULE_LOG_ORDER)}


def _natural_sort_key(text):
    """Splits digit runs out as ints, so e.g. '..._2.log' sorts before
    '..._10.log' -- plain string sort would put '_10' first."""
    return [int(chunk) if chunk.isdigit() else chunk
            for chunk in re.split(r"(\d+)", text)]


def _log_sort_key(path):
    rule_name = os.path.basename(os.path.dirname(path))
    ## Unknown rule names (e.g. a future rule this list hasn't been updated
    ## for) sort after every known one, instead of being silently dropped.
    return (_RULE_LOG_ORDER_INDEX.get(rule_name, len(RULE_LOG_ORDER)),
            rule_name, _natural_sort_key(os.path.basename(path)))


def _index_transcript_logs():
    """Every file under logs/, grouped by the transcript it belongs to.
    """
    index = {}
    if not os.path.isdir("logs"):
        return index
    wanted = set(transcripts)
    for rule_dir in os.scandir("logs"):
        if not rule_dir.is_dir():
            continue
        for entry in os.scandir(rule_dir.path):
            if not entry.name.endswith(".log"):
                continue
            stem = entry.name[:-4]
            ## Logs are named {transcript_id}.log or {transcript_id}_<extra>.log,
            ## where <extra> is a tree keyword, a RELAX round, or both:
            ## _i_tree, _3, _i_tree_3.
            tid = stem if stem in wanted else None
            cut = len(stem)
            while tid is None:
                cut = stem.rfind("_", 0, cut)
                if cut <= 0:
                    break
                if stem[:cut] in wanted:
                    tid = stem[:cut]
            if tid is not None:
                index.setdefault(tid, []).append(entry.path)
    return index


## Written by easel's own SIGINT handler (easel/cli.py) just before it
## forwards the signal on to snakemake.
INTERRUPT_MARKER = ".easel_interrupted"


def run_was_interrupted():
    """True when this run is ending because the user pressed Ctrl-C.
    """
    return os.path.exists(INTERRUPT_MARKER)


def print_interrupt_note():
    """One line instead of a post-mortem. Nothing here is a diagnosis: the run
    was cancelled, the outputs on disk are still valid, and rerunning the same
    command picks up where it stopped."""
    print("\nInterrupted before the run finished. Per-transcript logs were "
          "not bundled\nand no run report was written -- rerun the same "
          "command to continue.\n")


def bundle_transcript_logs():
    """Concatenate every per-rule log for each transcript into one file,
    in pipeline stage order (RULE_LOG_ORDER): codon_alignments/{transcript_id}/{transcript_id}_full.log.

    Called from onsuccess/onerror, not a rule.
    """
    index = _index_transcript_logs()
    for tid in transcripts:
        matches = sorted(index.get(tid, ()), key=_log_sort_key)
        if not matches:
            continue
        out_dir = f"codon_alignments/{tid}"
        os.makedirs(out_dir, exist_ok=True)
        with open(f"{out_dir}/{tid}_full.log", "w") as out:
            for match in matches:
                out.write(f"----- {match} -----\n")
                with open(match) as fh:
                    out.write(fh.read())
                out.write("\n")


def print_run_report():
    """Final tally, printed once the whole run finishes (onsuccess/onerror)

    Reports BOTH the outcome and how far each transcript got. Outcome alone is
    not enough: a transcript with no done.txt may have failed anywhere, or
    simply never been reached because something upstream failed first.

    Milestones are inferred from files, and tmp.tar.gz counts as evidence for
    every earlier one, since tar_compress removes tmp/ on success.
    """
    succeeded = skipped = 0
    reached = {"alignment": 0, "cleaning": 0, "validation": 0, "done": 0}

    for tid in transcripts:
        d = f"codon_alignments/{tid}"
        tmp = f"{d}/tmp"
        archived = os.path.exists(f"{d}/tmp.tar.gz")

        def has(*names):
            return archived or any(os.path.exists(f"{tmp}/{n}") for n in names)

        if has(f"{tid}_ori.fa", f"{tid}_ori.best.fas"):
            reached["alignment"] += 1
        if has(f"{tid}.manual.fa", f"{tid}.masked.hmm_cleaned.fa", f"{tid}.masked.fa"):
            reached["cleaning"] += 1
        if os.path.exists(f"{d}/validation.txt"):
            reached["validation"] += 1

        done_path = f"{d}/done.txt"
        if os.path.exists(done_path):
            reached["done"] += 1
            verdict = Path(done_path).read_text().split("\t")[0].strip()
            if verdict == "OK":
                succeeded += 1
            else:
                skipped += 1

    n = len(transcripts)
    failed = n - succeeded - skipped
    print(f"\n{n} valid input transcript(s): {succeeded} finished successfully, "
          f"{skipped} skipped (see skipped_transcripts.tsv), "
          f"{failed} did not complete.")
    print("  reached: " + ", ".join(f"{k} {v}/{n}" for k, v in reached.items()))
    ## Name the stage the run stalled at instead of leaving it to be inferred.
    order = ["alignment", "cleaning", "validation", "done"]
    drops = []
    for i, stage in enumerate(order):
        prev = n if i == 0 else reached[order[i - 1]]
        if reached[stage] < prev:
            drops.append((prev - reached[stage], stage, prev, reached[stage]))
    if drops:
        lost, stage, prev, got = max(drops)
        print(f"  biggest drop entering '{stage}': {prev} transcripts had the "
              f"previous stage, {got} reached this one ({lost} lost)")
        print(f"  -> look at logs/*{stage[:5]}*/ and the first error in the "
              f"snakemake log")
    print()
