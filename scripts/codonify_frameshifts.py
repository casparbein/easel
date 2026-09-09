#!/usr/bin/env python3

import argparse
import logging
import os
import re
import sys

logger = logging.getLogger(__name__)


class _Formatter(logging.Formatter):
    """Plain text for the report (DEBUG/INFO); "LEVEL: message" for warnings and errors."""

    def format(self, record: logging.LogRecord) -> str:
        if record.levelno <= logging.INFO:
            return record.getMessage()
        return f"{record.levelname}: {record.getMessage()}"


def configure_logging(level: int = logging.WARNING, log_file=None,
                      console: bool = True) -> None:
    """Send the unit-by-unit report to stdout and warnings/errors to stderr.

    Keeping the two apart preserves the streams the script used before it was
    switched from print() to logging: `script.py -q sp1 > report.txt` still
    captures the report, while warnings stay visible on the terminal.

    `log_file` adds a handler that appends everything (report and warnings) to
    that file, which is how a Snakemake rule's `log:` directive is honoured -
    the same thing `>> {log} 2>&1` did for the shell version of the rule. When
    a log file is used the console handlers are normally dropped (console=False)
    so the output is not duplicated into Snakemake's own job output.
    """
    logger.setLevel(level)
    logger.propagate = False

    handlers = []
    if console:
        report = logging.StreamHandler(sys.stdout)
        report.setLevel(logging.DEBUG)
        report.addFilter(lambda record: record.levelno <= logging.INFO)
        report.setFormatter(_Formatter())

        problems = logging.StreamHandler(sys.stderr)
        problems.setLevel(logging.WARNING)
        problems.setFormatter(_Formatter())
        handlers += [report, problems]

    if log_file:
        to_file = logging.FileHandler(str(log_file), mode="a", encoding="utf-8")
        to_file.setLevel(logging.DEBUG)
        to_file.setFormatter(_Formatter())
        handlers.append(to_file)

    logger.handlers[:] = handlers


def _verbosity(verbose: bool, query) -> int:
    """DEBUG reports every sequence, INFO only the -q one, WARNING neither."""
    if verbose:
        return logging.DEBUG
    if query:
        return logging.INFO
    return logging.WARNING


configure_logging()


DESCRIPTION = '''\
Script to turn a prank nucleotide alignment into codons, based on reference triplets.
Removes columns that are masked N or - entirely, masks all nucleotides in codons of query
species if they are not part of full codons (for example AG-). Also masks stop codons
introduced through frameshift removal.

Output:
-> Masked alignment in fasta format
-> file identifying frameshift insertions and deletions based on reference codons.
   Two column format with triplet position in the reference in one column and the frameshift
   type in the second column
-> file that contains a premasked version of the alignment, where no columns are removed yet
   and reference based triplets are separated by "***".
-> TSV table with per-unit before/after masking for all queries (.units.tsv)

Input:
-> Prank/Muscle nucleotide alignment (fasta format) (-i)
-> reference name in alignment (-r)
-> output base name (-o)
-> masking regime (-m N|all): only mask frameshifts, or mask N downstream in-frame codons
-> upstream masking (-u N): additionally mask the N reference codons immediately preceding
   an affected codon (0 = off, the default; 1 = the codon right before the frameshift)
-> masking of suspect C-terminus regions (-c CODON_NUM): any frameshift within the last
   CODON_NUM reference codons triggers masking of the entire remaining C-terminal zone,
   even if the reading frame recovers
-> query (-q): gives verbose output of what codons were masked why, for a given sequence
-> --verbose: like query, but for all sequences

Snakemake:
The script doubles as a Snakemake `script:` target. When Snakemake injects its
`snakemake` object the command line is bypassed and the run is configured from
the rule instead:
    input:  ali            (or the first input file)
    params: reference      (required)
            out_head       (optional; otherwise inferred from the output paths)
            mask_downstream / mask_upstream / cterm_mask / query / verbose (optional)
    output: codon_ali      -> <out_head>.masked.fa
            premask_ali    -> <out_head>.premasked.fa
            frameshifts    -> <out_head>.frameshifts.txt
            units          -> <out_head>.units.tsv
    log:    everything the script reports is appended to log[0]
'''

# REF: single source of truth for the character classes that were previously
#      re-created as set literals ~12 times inside the hot loops.
GAP = "-"
BASES = frozenset("ACGT")
IUPAC_AMBIGUITY = frozenset("NRYSWKMBDHV")
NUCLEOTIDE_CHARS = BASES | IUPAC_AMBIGUITY
STOP_CODONS = frozenset(("TAA", "TGA", "TAG"))

# The four files every run produces, keyed by the output name used in the
# Snakemake rule. Single source of truth for the writers below and for the
# check that a rule's `output:` paths agree with the chosen base name.
OUTPUT_SUFFIXES = {
    "codon_ali": ".masked.fa",
    "premask_ali": ".premasked.fa",
    "frameshifts": ".frameshifts.txt",
    "units": ".units.tsv",
}


class CodonifyError(Exception):
    """A problem with the input or the requested output that the user has to fix."""


# ----------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------
def count_bases(seq) -> int:
    """Count nucleotide positions (ACGT + IUPAC ambiguity codes) in a sequence."""
    return sum(1 for c in seq if c in NUCLEOTIDE_CHARS)


def compress_gaps(seq) -> str:
    """Collapse runs of 4+ dashes into -(N) notation for display."""
    s = "".join(seq)
    return re.sub(r"-{4,}", lambda m: f"-({len(m.group())})", s)


cg = compress_gaps


def mask_unit(unit: list) -> list:
    """Return a copy of a unit with every non-gap position replaced by N."""
    return ["N" if c != GAP else c for c in unit]


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def parse_mask_downstream(value: str):
    """argparse type for --mask-downstream: 'all' or a non-negative integer."""
    if value == "all":
        return "all"
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"must be 'all' or a non-negative integer, got '{value}'"
        )
    if n < 0:
        raise argparse.ArgumentTypeError(f"must be >= 0, got {n}")
    return n


def non_negative_int(value: str) -> int:
    """argparse type for --mask-upstream / --cterm-mask."""
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"must be a non-negative integer, got '{value}'")
    if n < 0:
        raise argparse.ArgumentTypeError(f"must be >= 0, got {n}")
    return n


def argument_parser() -> argparse.Namespace:
    """Parse CMD args."""
    app = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(
            prog, max_help_position=6, indent_increment=2
        ),
    )
    app.add_argument(
        "-i", "--input_fasta",
        action="store", dest="input", type=str, required=True,      
        help="Input fasta file created by a nucleotide aligner like prank",
    )
    app.add_argument(
        "-r", "--reference",
        action="store", dest="reference", type=str, required=True,  
        help='Name of the reference. The ">" prefix is added automatically if missing.',
    )
    app.add_argument(
        "-o", "--output_fasta",
        action="store", dest="output", type=str, required=True,
        help='Base name of the output files (e.g. out).\n'
             'Files created: out.masked.fa, out.frameshifts.txt, out.premasked.fa, out.units.tsv',
    )
    app.add_argument(
        "-q", "--query",
        action="store", dest="query", type=str, default=None,
        help='Name of a specific query for detailed unit-by-unit logging.\n'
             'The ">" prefix is added automatically if missing.',
    )
    app.add_argument(
        "-v", "--verbose",
        action="store_true", dest="verbose",
        help="Enable verbose/debug output (all queries)",
    )
    app.add_argument(
        "-m", "--mask-downstream",
        action="store", dest="mask_downstream",
        type=parse_mask_downstream, default="all", 
        metavar="N|all",
        help='How many in-frame codons to mask downstream of an unresolved frameshift.\n'
             '"all" (default): mask until frameshift is resolved.\n'
             '0: only mask the frameshift codon itself, no downstream masking.\n'
             'N (integer): mask at most N downstream in-frame codons.',
    )
    app.add_argument(
        "-u", "--mask-upstream",
        action="store", dest="mask_upstream",
        type=non_negative_int, default=0,
        metavar="N",
        help='How many reference codons immediately UPSTREAM of an affected codon to mask.\n'
             '0 (default): mask the affected codon and downstream only.\n'
             '1: also mask the codon immediately preceding the frameshift.\n'
             'N: mask the N reference codons preceding it. Reference-insertion units\n'
             '   inside that window are masked too, so the masked block stays contiguous.',
    )
    app.add_argument(
        "-c", "--cterm-mask",
        action="store", dest="cterm_mask",
        type=non_negative_int, default=None,
        metavar="N",
        help='C-terminal masking: define the last N reference codons as the C-terminal zone.\n'
             'Any frameshift occurring within the zone sets a permanent latch that masks\n'
             'all remaining C-terminal units, regardless of -m and even if the frame recovers.\n'
             'Example: -c 100',
    )

    if len(sys.argv) == 1:
        app.print_help(sys.stderr)
        sys.exit(1)

    return app.parse_args()


# ----------------------------------------------------------------------------
# IO
# ----------------------------------------------------------------------------
def read_fasta(fasta_path: str) -> dict:
    """Read a fasta file. Returns dict {header: sequence_string}, sequences upper-cased."""
    chunks = {}                    
    header = None
    saw_lowercase = False

    with open(fasta_path, "r") as fasta:
        for line_no, line in enumerate(fasta, 1):
            fasta_line = line.strip()
            if not fasta_line:
                continue
            if fasta_line.startswith(">"):
                if fasta_line in chunks:
                    raise ValueError(
                        f"duplicate header '{fasta_line}' at line {line_no} in {fasta_path}"
                    )
                header = fasta_line
                chunks[header] = []
            else:
                if header is None:
                    raise ValueError(
                        f"sequence data before the first '>' header at line {line_no} in {fasta_path}"
                    )
                if not saw_lowercase and fasta_line != fasta_line.upper():
                    saw_lowercase = True
                chunks[header].append(fasta_line.upper())
    if saw_lowercase:
        logger.warning("Input contains lower-case sequence characters; they were upper-cased.")
    if not chunks:
        raise ValueError(f"no sequences found in {fasta_path}")
    return {head: "".join(parts) for head, parts in chunks.items()}


def validate_alignment(fasta_dict: dict, reference: str) -> None:
    """Fail fast on ragged input; every downstream index assumes equal lengths.""" 
    ref_len = len(fasta_dict[reference])
    ragged = [(h, len(s)) for h, s in fasta_dict.items() if len(s) != ref_len]
    if ragged:
        shown = ", ".join(f"{h} ({n} cols)" for h, n in ragged[:5])
        more = f" (+{len(ragged) - 5} more)" if len(ragged) > 5 else ""
        raise ValueError(
            f"input is not an alignment: {len(ragged)} sequence(s) differ in length from the "
            f"reference ({ref_len} cols): {shown}{more}"
        )
    unexpected = {c for s in fasta_dict.values() for c in s} - NUCLEOTIDE_CHARS - {GAP}
    if unexpected:
        logger.warning(
            "Unexpected characters in the alignment: %s. They are neither masked nor counted.",
            ", ".join(sorted(repr(c) for c in unexpected)),
        )


def define_units(fasta_dict: dict, reference: str) -> list:
    """Define codon units based on the reference sequence.

    Returns a list of [start, end] inclusive column index pairs. A unit is either a
    reference codon (3 reference bases plus any gap columns interleaved with them)
    or a run of reference-only gap columns (an insertion relative to the reference).
    """
    fasta_string = fasta_dict[reference]
    last = len(fasta_string) - 1

    units = []
    current_nts = 0
    in_gap_run = False
    start_nt_index = 0
    start_gap_index = 0

    for i, char in enumerate(fasta_string):
        if char != GAP:
            if current_nts == 0:
                start_nt_index = i
                current_nts = 1
            elif current_nts != 2:
                current_nts += 1
            else:
                units.append([start_nt_index, i])
                current_nts = 0
        else:
            if current_nts != 0:
                continue                  # gap inside an incomplete codon: swallowed by that codon
            if not in_gap_run:
                start_gap_index = i
                in_gap_run = True
            if i == last or fasta_string[i + 1] != GAP:
                units.append([start_gap_index, i])
                in_gap_run = False

    if current_nts != 0:                       
        units.append([start_nt_index, last])   
        logger.warning(
            "Reference ends with an incomplete codon (%d base(s)); columns %d-%d were kept "
            "as a trailing partial unit and will be masked by the sanity check.",
            current_nts, start_nt_index + 1, last + 1,
        )
    return units


def write_premasked_fasta(fasta_dict: dict, units: list, reference: str, out_base: str) -> None:
    """Write the premasked fasta file with ***-delimited units."""
    with open(out_base + OUTPUT_SUFFIXES["premask_ali"], "w") as out:
        ref_string = fasta_dict[reference]
        out.write(f"{reference}\n")
        out.write("***".join(ref_string[start:end + 1] for start, end in units) + "\n")
        for header, fasta_string in fasta_dict.items():      # REF: items() instead of keys() + re-lookup
            if header == reference:
                continue
            out.write(f"{header}\n")
            out.write("***".join(fasta_string[start:end + 1] for start, end in units) + "\n")


# ----------------------------------------------------------------------------
# Masking helpers
# ----------------------------------------------------------------------------
def reference_codon_flags(reference_substrings: list) -> list:
    """True for every unit that is a real reference codon (i.e. not an all-gap unit).

    Uses "contains a non-gap character" rather than "contains ACGT" so the answer is
    invariant to reference masking (stop codons / ambiguities become N, never gaps).
    """
    return [any(c != GAP for c in sub) for sub in reference_substrings]


def compute_cterm_after(codon_flags: list) -> list:
    """For each unit i, the number of reference codons at positions j > i.

    Unit i is in the C-terminal zone when cterm_after[i] < cterm_mask.
    """
    cterm_after = [0] * len(codon_flags)
    count = 0
    for j in range(len(codon_flags) - 1, -1, -1):
        cterm_after[j] = count
        if codon_flags[j]:
            count += 1
    return cterm_after


def upstream_span_start(codon_flags: list, i: int, n_codons: int) -> int:
    """First unit index of the window covering the n_codons reference codons before unit i."""  # NEW(F)
    start = i
    found = 0
    j = i - 1
    while j >= 0 and found < n_codons:
        if codon_flags[j]:
            found += 1
        start = j
        j -= 1
    return start


def _mask_action(cterm_masked: bool, downstream_masked: bool,
                 downstream_count: int, mask_downstream) -> str:
    """Build a consistent action string for the per-unit log."""
    lim = str(mask_downstream)
    if cterm_masked and downstream_masked:
        return f"masked (C-term + downstream {downstream_count}/{lim})"
    if cterm_masked:
        return "masked (C-term latch)"
    if downstream_masked:
        return f"masked (downstream {downstream_count}/{lim})"
    return f"SKIP (limit {lim} reached)"


def _should_mask_downstream(mask_downstream, downstream_good_masked: int) -> bool:
    """REF: this predicate was duplicated verbatim in cases 2b, 4a and 5."""
    return mask_downstream == "all" or (
        isinstance(mask_downstream, int) and downstream_good_masked < mask_downstream
    )


# ----------------------------------------------------------------------------
# Core
# ----------------------------------------------------------------------------
def reformat_queries(fasta_dict: dict, unit_list: list, reference: str,
                     debug_query=None, mask_downstream="all",
                     mask_upstream: int = 0, cterm_mask=None):
    """Core masking logic.

    Returns: (reformatted_dict, all_frameshifts, unit_rows)
      unit_rows: list of dicts {unit, reference, query, before, after}
                 containing pre- and post-masking unit strings for all queries.
    """
    all_frameshifts = {}
    reformatted_dict = {}
    unit_rows = []

    reference_fasta_string = fasta_dict[reference]
    reference_substrings = [
        list(reference_fasta_string[start:end + 1]) for start, end in unit_list
    ]

    # Mask reference stop codons and IUPAC ambiguous nucleotides
    ref_frameshifts = {}
    for i, sub in enumerate(reference_substrings):
        letters_only = "".join(c for c in sub if c in BASES)
        if len(letters_only) == 3 and letters_only in STOP_CODONS:
            logger.debug(f"Reference stop codon '{letters_only}' at unit {i+1}, masking to N")
            reference_substrings[i] = ["N" if c in BASES else c for c in sub]
            ref_frameshifts[i + 1] = "STOP_CODON"
    for i, sub in enumerate(reference_substrings):
        if i + 1 in ref_frameshifts:
            continue
        letters_only = "".join(c for c in sub if c != GAP)
        if len(letters_only) == 3 and any(c in IUPAC_AMBIGUITY for c in letters_only):
            logger.debug(f"Reference ambiguous nucleotide(s) '{letters_only}' at unit {i+1}, masking to N")
            reference_substrings[i] = mask_unit(sub)
            ref_frameshifts[i + 1] = "AMBIGUOUS"

    all_frameshifts["__reference__"] = ref_frameshifts
    # REF: keep every sequence as a list of unit-lists (the reference used to be a list of
    #      strings, which forced final_ali_filtering to cope with two different types).
    reformatted_dict[reference] = [list(entry) for entry in reference_substrings]

    codon_flags = reference_codon_flags(reference_substrings)
    cterm_after = compute_cterm_after(codon_flags) if cterm_mask is not None else None

    # REF: hoisted out of the per-query loop - it does not depend on the query
    first_cterm_unit = None
    if cterm_mask is not None:
        first_cterm_unit = next(
            (i for i in range(len(reference_substrings)) if cterm_after[i] < cterm_mask), None
        )

    n_units = len(reference_substrings)

    for header, fasta_string in fasta_dict.items():
        if header == reference:
            continue

        is_debug_query = (header == debug_query) or logger.isEnabledFor(logging.DEBUG)
        if is_debug_query:
            logger.info(f"\n{'='*72}")
            logger.info(f"  DETAILED CODONIFY LOG FOR: {header}")
            logger.info(f"{'='*72}")
            if cterm_mask is not None:
                logger.info(f"  C-terminal zone (-c {cterm_mask}): starts at unit "
                            f"{first_cterm_unit + 1 if first_cterm_unit is not None else 'n/a'} "
                            f"(units with < {cterm_mask} ref codons remaining)")
            if mask_upstream:
                logger.info(f"  Upstream masking (-u {mask_upstream}): {mask_upstream} reference codon(s) "
                            f"before every affected codon are masked as well")

        query_substrings = [
            list(fasta_string[start:end + 1]) for start, end in unit_list
        ]
        query_adapted = [sublist[:] for sublist in query_substrings]

        frameshift_counter = 0
        # One-way latch: set to True the first time a frameshift fires within the
        # C-terminal zone. Never reset - ensures all subsequent C-term units are
        # masked even after the reading frame recovers (counter returns to 0).
        frameshift_in_cterm = False
        downstream_good_masked = 0
        upstream_masked_units = set()
        frameshifts = {}

        def register_frameshift(idx: int, n_bases: int, in_cterm: bool) -> str:
            """Advance the frameshift state and retroactively mask upstream codons.

            Returns a suffix for the per-unit debug line (empty when nothing was masked).
            """
            nonlocal frameshift_counter, frameshift_in_cterm, downstream_good_masked
            if not frameshift_in_cterm:
                frameshift_in_cterm = in_cterm
            frameshift_counter += n_bases
            if frameshift_counter % 3 == 0:
                frameshift_counter = 0
                downstream_good_masked = 0 
            if mask_upstream <= 0 or idx == 0:
                return ""
            start = upstream_span_start(codon_flags, idx, mask_upstream)
            changed = []
            for k in range(start, idx):
                upstream_masked_units.add(k + 1)
                masked = mask_unit(query_adapted[k])
                if masked != query_adapted[k]:
                    query_adapted[k] = masked
                    changed.append(k + 1)
            if not changed:
                return ""
            span = f"{changed[0]}-{changed[-1]}" if len(changed) > 1 else str(changed[0])
            return f" | upstream -u {mask_upstream}: masked unit(s) {span}"

        for i in range(n_units):
            ref_unit = reference_substrings[i]
            qry_unit = query_substrings[i]
            ref_str = "".join(ref_unit)
            qry_str = "".join(qry_unit)
            ref_all_gaps = ref_str == GAP * len(ref_unit)
            in_cterm = cterm_mask is not None and cterm_after[i] < cterm_mask
            cterm_info = ""
            if cterm_mask is not None:
                cterm_info = (f" [cterm: {cterm_after[i]}<{cterm_mask}]"
                              if in_cterm
                              else f" [cterm: {cterm_after[i]}>={cterm_mask},out]")

            # ------------------------------------------------------------------
            # Case 2/3: the reference unit is nothing but gaps (insertion in the query).
            # ------------------------------------------------------------------
            if ref_all_gaps:
                if len(ref_unit) % 3 == 0:
                    # Case 2a: query has IUPAC or mixed bases+gaps -> mask non-gaps
                    if any(c in IUPAC_AMBIGUITY for c in qry_str) or GAP in qry_str:
                        query_adapted[i] = mask_unit(qry_unit)
                        up_note = ""
                        # Only count as a frameshift if the query has both ACGT and gaps
                        if any(b in BASES for b in qry_str) and GAP in qry_str:
                            up_note = register_frameshift(i, count_bases(qry_unit), in_cterm)
                        if is_debug_query:
                            logger.info(f"  Unit {i+1:>4d} | Case 2a (ref gap, mod3) | "
                                        f"ref: {cg(ref_str)} | "
                                        f"query: {cg(qry_str)} -> {cg(query_adapted[i])} | "
                                        f"FS counter: {frameshift_counter} | masked (partial){up_note}")

                    # Case 2b: clean query but active frameshift -> downstream masking
                    elif frameshift_counter % 3 != 0:
                        cterm_should_mask = frameshift_in_cterm and in_cterm
                        downstream_should_mask = _should_mask_downstream(
                            mask_downstream, downstream_good_masked)
                        if cterm_should_mask or downstream_should_mask:
                            query_adapted[i] = mask_unit(qry_unit)
                            downstream_good_masked += 1
                        action = _mask_action(cterm_should_mask, downstream_should_mask,
                                              downstream_good_masked, mask_downstream)
                        if is_debug_query:
                            logger.info(f"  Unit {i+1:>4d} | Case 2b (ref gap, mod3) | "
                                        f"ref: {cg(ref_str)} | "
                                        f"query: {cg(qry_str)} -> {cg(query_adapted[i])} | "
                                        f"FS counter: {frameshift_counter} | {action}{cterm_info}")

                    # Case 2c: in-frame, keep as-is
                    else:
                        query_adapted[i] = qry_unit[:]       # FIX(B12): was an alias of qry_unit
                        if is_debug_query:
                            logger.info(f"  Unit {i+1:>4d} | Case 2c (ref gap, mod3) | "
                                        f"ref: {cg(ref_str)} | "
                                        f"query: {cg(qry_str)} -> {cg(query_adapted[i])} | "
                                        f"FS counter: {frameshift_counter} | in-frame (kept)")

                # ------------------------------------------------------------------
                # Case 3: ref is all-gaps, not mod 3 -> frameshift deletion in query
                # ------------------------------------------------------------------
                else:
                    frameshifts[i + 1] = "DEL"
                    query_adapted[i] = mask_unit(qry_unit)
                    up_note = ""
                    if ref_unit != qry_unit and any(b in BASES for b in qry_str):
                        up_note = register_frameshift(i, count_bases(qry_unit), in_cterm)
                    if is_debug_query:
                        logger.info(f"  Unit {i+1:>4d} | Case 3 (FS deletion)    | "
                                    f"ref: {cg(ref_str)} | "
                                    f"query: {cg(qry_str)} -> {cg(query_adapted[i])} | "
                                    f"FS counter: {frameshift_counter} | DEL{up_note}")

            # ------------------------------------------------------------------
            # Case 1: reference unit is exactly 3 slots and query has a gap or IUPAC
            # ------------------------------------------------------------------
            elif len(ref_unit) == 3 and (
                GAP in qry_unit or any(c in IUPAC_AMBIGUITY for c in qry_unit)
            ):
                # Case 1a: all-gap codon, or IUPAC only (no real gap mixed with bases)
                if qry_str == GAP * 3 or GAP not in qry_str:
                    query_adapted[i] = (qry_unit[:] if qry_str == GAP * 3      # FIX(B12)
                                        else ["N"] * len(qry_unit))            # FIX(B13): was len(ref_unit)
                    if is_debug_query:
                        logger.info(f"  Unit {i+1:>4d} | Case 1a (gap/IUPAC)     | "
                                    f"ref: {cg(ref_str)} | "
                                    f"query: {cg(qry_str)} -> {cg(query_adapted[i])} | "
                                    f"FS counter: {frameshift_counter}")

                # Case 1b: partial gap mixed with bases -> frameshift
                else:
                    current_frameshift = count_bases(qry_unit)
                    query_adapted[i] = mask_unit(qry_unit)
                    frameshifts[i + 1] = "INS"
                    up_note = register_frameshift(i, current_frameshift, in_cterm)
                    if is_debug_query:
                        logger.info(f"  Unit {i+1:>4d} | Case 1b (frameshift)    | "
                                    f"ref: {cg(ref_str)} | "
                                    f"query: {cg(qry_str)} -> {cg(query_adapted[i])} | "
                                    f"FS +{current_frameshift}, counter: {frameshift_counter} | INS{up_note}")

            # ------------------------------------------------------------------
            # Case 4: ref has mixed bases and gaps (partial codon with insertions)
            # ------------------------------------------------------------------
            elif len(ref_unit) != 3 and GAP in ref_unit:
                frameshifts[i + 1] = "DEL"
                # positions where reference and query agree on "gap vs. base"
                equal_counter = sum(1 for r, q in zip(ref_unit, qry_unit)
                                    if (r == GAP) == (q == GAP))
                query_adapted[i] = mask_unit(qry_unit)
                up_note = ""

                # Sub-case 4a: query positions align with reference
                if equal_counter == len(ref_unit):
                    if frameshift_counter % 3 == 0:
                        query_adapted[i] = qry_unit[:]         # FIX(B12)
                        frameshift_counter = 0
                        downstream_good_masked = 0
                        action = "in-frame (unmasked)"
                    else:
                        cterm_should_mask = frameshift_in_cterm and in_cterm
                        downstream_should_mask = _should_mask_downstream(
                            mask_downstream, downstream_good_masked)
                        if cterm_should_mask or downstream_should_mask:
                            downstream_good_masked += 1
                        else:
                            query_adapted[i] = qry_unit[:]     # FIX(B12)
                        action = _mask_action(cterm_should_mask, downstream_should_mask,
                                              downstream_good_masked, mask_downstream)

                # Sub-case 4b: query bases are in-frame (mod 3) but don't align with reference
                elif count_bases(qry_unit) % 3 == 0 and count_bases(qry_unit) > 0:
                    action = "masked, bases in frame"

                # Sub-case 4c: actual frameshift
                else:
                    up_note = register_frameshift(i, count_bases(qry_unit), in_cterm)
                    action = "DEL"

                if is_debug_query:
                    logger.info(f"  Unit {i+1:>4d} | Case 4 (mixed ref+gap)  | "
                                f"ref: {cg(ref_str)} | "
                                f"query: {cg(qry_str)} -> {cg(query_adapted[i])} | "
                                f"FS counter: {frameshift_counter} | {action}{cterm_info}{up_note}")

            # ------------------------------------------------------------------
            # Case 5: normal in-frame codon inside a frameshifted block
            # ------------------------------------------------------------------
            elif (
                len(ref_unit) == 3
                and GAP not in qry_unit
                and GAP not in ref_unit
                and frameshift_counter % 3 != 0
            ):
                cterm_should_mask = frameshift_in_cterm and in_cterm
                downstream_should_mask = _should_mask_downstream(
                    mask_downstream, downstream_good_masked)
                if cterm_should_mask or downstream_should_mask:
                    query_adapted[i] = mask_unit(qry_unit)
                    downstream_good_masked += 1
                action = _mask_action(cterm_should_mask, downstream_should_mask,
                                      downstream_good_masked, mask_downstream)
                if is_debug_query:
                    logger.info(f"  Unit {i+1:>4d} | Case 5 (shifted frame)  | "
                                f"ref: {cg(ref_str)} | "
                                f"query: {cg(qry_str)} -> {cg(query_adapted[i])} | "
                                f"FS counter: {frameshift_counter} | {action}{cterm_info}")

            # ------------------------------------------------------------------
            # Post-case: C-terminus latch override
            # Applies after any case when frameshift_in_cterm is set and we are
            # still inside the C-terminal zone - catches units left unmasked when
            # the frameshift counter drops back to 0 (frame recovery within C-term).
            # ------------------------------------------------------------------
            if frameshift_in_cterm and in_cterm:
                pre_latch = "".join(query_adapted[i])
                query_adapted[i] = mask_unit(query_adapted[i])
                if is_debug_query and "".join(query_adapted[i]) != pre_latch:
                    logger.info(f"  Unit {i+1:>4d} | C-term latch (override) | "
                                f"ref: {cg(ref_str)} | "
                                f"{cg(pre_latch)} -> {cg(query_adapted[i])} | "
                                f"masked (C-term latch, FS counter={frameshift_counter}){cterm_info}")

        # Post-masking: scan all units for in-frame stop codons
        for i, unit in enumerate(query_adapted):
            base_positions = [j for j, c in enumerate(unit) if c in BASES]
            if not base_positions or len(base_positions) % 3 != 0:
                continue
            before = "".join(unit)
            stops_found = []                     # REF: dropped the redundant `has_stop` flag
            for k in range(0, len(base_positions), 3):
                triplet_pos = base_positions[k:k + 3]
                triplet = "".join(unit[j] for j in triplet_pos)
                if triplet in STOP_CODONS:
                    stops_found.append(triplet)
                    for j in triplet_pos:
                        unit[j] = "N"
            if stops_found:
                frameshifts[i + 1] = "STOP_CODON"
                if is_debug_query:
                    logger.info(f"  Unit {i+1:>4d} | Stop codon masking      | "
                                f"{cg(before)} -> {cg(unit)} | "
                                f"codons: {stops_found}")

        # Sanity check: every unit must have 0 or a multiple of 3 raw ACGT bases.
        # N and other IUPAC codes are already masked and are fine - only check ACGT.
        for i, unit in enumerate(query_adapted):
            base_count = sum(1 for c in unit if c in BASES)
            if base_count > 0 and base_count % 3 != 0:
                before = "".join(unit)
                query_adapted[i] = ["N" if c in NUCLEOTIDE_CHARS else c for c in unit]
                msg = (f"Sanity fix for {header} at unit {i+1}: "
                       f"base count {base_count} not divisible by 3. "
                       f"Masked: {cg(before)} -> {cg(query_adapted[i])}")
                logger.warning(msg)

        if is_debug_query:
            logger.info(f"\n  {'-'*68}")
            logger.info("  PREMASKED UNITS (before masking) vs MASKED UNITS (after masking):")
            logger.info(f"  {'-'*68}")
            logger.info(f"  {'Unit':>6s}  {'Reference':^20s}  {'Query (before)':^20s}  {'Query (after)':^20s}")
            logger.info(f"  {'-'*68}")
            for i in range(n_units):
                ref_s = cg(reference_substrings[i])
                pre_s = cg(query_substrings[i])
                post_s = cg(query_adapted[i])
                changed = " *" if query_substrings[i] != query_adapted[i] else ""
                logger.info(f"  {i+1:>6d}  {ref_s:^20s}  {pre_s:^20s}  {post_s:^20s}{changed}")
            logger.info(f"  {'-'*68}")
            logger.info("  (* = unit was modified during masking)")
            logger.info(f"  Frameshifts detected: {frameshifts if frameshifts else 'none'}")
            if mask_upstream:
                logger.info(f"  Units masked by -u {mask_upstream} (upstream of an affected codon): "
                            f"{sorted(upstream_masked_units) if upstream_masked_units else 'none'}")
            logger.info(f"{'='*72}\n")

        # Collect per-unit rows for TSV output (pre-filtered, pre-column-removal)
        for i in range(n_units):
            unit_rows.append({
                "unit": i + 1,
                "reference": "".join(reference_substrings[i]),
                "query": header,
                "before": "".join(query_substrings[i]),
                "after": "".join(query_adapted[i]),
            })

        reformatted_dict[header] = query_adapted
        all_frameshifts[header] = frameshifts

    return reformatted_dict, all_frameshifts, unit_rows


def write_unit_tsv(unit_rows: list, out_base: str) -> None:
    """Write the per-unit before/after masking table for all queries as a TSV.

    Columns: unit (1-based), reference unit string, query header,
             query unit before masking, query unit after masking, modified (yes/no).
    Note: unit strings are pre-column-filtering (the .masked.fa output may be shorter)
    and long gap runs are shown compressed as -(N).
    """
    with open(out_base + OUTPUT_SUFFIXES["units"], "w") as out:
        out.write("unit\treference\tquery\tbefore\tafter\tmodified\n")
        for row in unit_rows:
            modified = "yes" if row["before"] != row["after"] else "no"
            out.write(
                f"{row['unit']}\t{cg(row['reference'])}\t{row['query']}\t"
                f"{cg(row['before'])}\t{cg(row['after'])}\t{modified}\n"
            )


def final_ali_filtering(final_filtering: dict, all_frameshifts: dict, out_base: str) -> None:
    """Drop alignment columns without a single ACGT, write .masked.fa and .frameshifts.txt."""
    headers = list(final_filtering)
    # REF: flatten once instead of building a defaultdict of per-unit column lists and then
    #      calling list.index() per sequence (O(n^2), and silently wrong for duplicate headers).
    flat = [[c for unit in final_filtering[h] for c in unit] for h in headers]
    kept = [col for col in zip(*flat) if any(c in BASES for c in col)]

    with open(out_base + OUTPUT_SUFFIXES["codon_ali"], "w") as out:
        for idx, header in enumerate(headers):
            out.write(header + "\n")
            out.write("".join(col[idx] for col in kept) + "\n")

    merged_frameshifts = {}
    for header_fs in all_frameshifts.values():
        for pos, fs_type in header_fs.items():
            merged_frameshifts.setdefault(pos, fs_type)

    with open(out_base + OUTPUT_SUFFIXES["frameshifts"], "w") as out:
        for fs in sorted(merged_frameshifts):
            out.write(f"{fs}\t{merged_frameshifts[fs]}\n")


def run(input_fasta: str, output_base: str, reference: str,
        query=None, mask_downstream="all", mask_upstream: int = 0,
        cterm_mask=None) -> None:
    """Read, codonify and write all four output files.

    Independent of how the arguments arrived (command line or Snakemake rule).
    Raises CodonifyError for anything the user has to fix.
    """
    if not reference.startswith(">"):
        reference = ">" + reference

    debug_query = query
    if debug_query is not None and not debug_query.startswith(">"):
        debug_query = ">" + debug_query

    try:
        fasta_dict = read_fasta(input_fasta)
    except OSError as exc:
        raise CodonifyError(f"cannot read input FASTA: {exc}") from exc
    except ValueError as exc:
        raise CodonifyError(str(exc)) from exc

    if reference not in fasta_dict:
        raise CodonifyError(f"reference '{reference}' not found in input FASTA. "
                            f"Available headers: {list(fasta_dict)}")

    if debug_query is not None and debug_query not in fasta_dict:
        raise CodonifyError(f"query '{debug_query}' not found in input FASTA. "
                            f"Available headers: {list(fasta_dict)}")

    try:
        validate_alignment(fasta_dict, reference)
    except ValueError as exc:
        raise CodonifyError(str(exc)) from exc

    units = define_units(fasta_dict, reference)
    if not units:
        raise CodonifyError(f"reference '{reference}' contains no codon units.")

    try:
        write_premasked_fasta(fasta_dict, units, reference, output_base)

        query_premask, frameshifts, unit_rows = reformat_queries(
            fasta_dict, units, reference,
            debug_query=debug_query,
            mask_downstream=mask_downstream,
            mask_upstream=mask_upstream,
            cterm_mask=cterm_mask,
        )

        write_unit_tsv(unit_rows, output_base)
        final_ali_filtering(query_premask, frameshifts, output_base)
    except OSError as exc:
        raise CodonifyError(f"cannot write output files: {exc}") from exc


def main() -> None:
    """Command line entry point."""
    args = argument_parser()
    configure_logging(_verbosity(args.verbose, args.query))
    try:
        run(args.input, args.output, args.reference,
            query=args.query,
            mask_downstream=args.mask_downstream,
            mask_upstream=args.mask_upstream,
            cterm_mask=args.cterm_mask)
    except CodonifyError as exc:
        logger.error("%s", exc)
        sys.exit(1)


# ----------------------------------------------------------------------------
# Snakemake `script:` support
# ----------------------------------------------------------------------------
def _rule_value(container, key, default=None):
    """Read `key` from a Snakemake input/output/params object.

    Those are Namedlist instances: named entries are attributes, and `.get()`
    exists on newer Snakemake versions. Missing entries yield `default`.
    """
    if container is None:
        return default
    value = getattr(container, key, None)
    if value is None:
        getter = getattr(container, "get", None)
        if callable(getter):
            try:
                value = getter(key, None)
            except TypeError:
                value = None
    return default if value is None else value


def _rule_input(smk) -> str:
    """The alignment to codonify: input.ali, or the only/first input file."""
    container = getattr(smk, "input", None)
    value = _rule_value(container, "ali")
    if value is None and container is not None and len(container):
        value = container[0]
    if value is None:
        raise CodonifyError("the rule has no input alignment "
                            "(expected `input: ali=...` or a single input file)")
    return str(value)


def _rule_output_base(smk) -> str:
    """params.out_head, or the base name implied by the declared outputs."""
    base = _rule_value(getattr(smk, "params", None), "out_head")
    if base is not None:
        return str(base)

    output = getattr(smk, "output", None)
    for key, suffix in OUTPUT_SUFFIXES.items():
        value = _rule_value(output, key)
        if value is not None and str(value).endswith(suffix):
            return str(value)[: -len(suffix)]
    for value in (output or []):
        for suffix in OUTPUT_SUFFIXES.values():
            if str(value).endswith(suffix):
                return str(value)[: -len(suffix)]

    raise CodonifyError(
        "cannot work out the output base name: set `params: out_head=...`, or declare "
        "an output path ending in one of " + ", ".join(sorted(OUTPUT_SUFFIXES.values()))
    )


def _check_declared_outputs(smk, output_base: str) -> None:
    """Fail early if the rule declares an output this script will never write.

    Snakemake would otherwise run the whole job and then kill it with a generic
    "missing output files" message.
    """
    output = getattr(smk, "output", None)
    if output is None:
        return
    written = [output_base + suffix for suffix in OUTPUT_SUFFIXES.values()]
    # compare normalised, but report the paths as the user wrote them
    written_norm = {os.path.normpath(path) for path in written}
    declared_norm = {os.path.normpath(str(path)) for path in output}
    stray = [str(path) for path in output if os.path.normpath(str(path)) not in written_norm]
    if stray:
        raise CodonifyError(
            "declared output(s) " + ", ".join(stray) + " will never be written. "
            "With out_head '" + output_base + "' this rule writes " +
            ", ".join(sorted(written)) + ". Fix params.out_head or the output: paths."
        )
    undeclared = sorted(path for path in written
                        if os.path.normpath(path) not in declared_norm)
    if undeclared:
        logger.info("Files written but not declared as rule outputs: %s "
                    "(Snakemake will not track or clean them)", ", ".join(undeclared))


def _rule_masking(smk):
    """Masking options from params, validated exactly like the command line flags."""
    params = getattr(smk, "params", None)
    try:
        mask_downstream = parse_mask_downstream(str(_rule_value(params, "mask_downstream", "all")))
        mask_upstream = non_negative_int(str(_rule_value(params, "mask_upstream", 0)))
        cterm = _rule_value(params, "cterm_mask")
        cterm_mask = None if cterm is None else non_negative_int(str(cterm))
    except argparse.ArgumentTypeError as exc:
        raise CodonifyError(f"invalid masking parameter in the rule's params: {exc}") from exc
    return mask_downstream, mask_upstream, cterm_mask


def main_snakemake(smk) -> None:
    """Entry point used when Snakemake runs this file through a `script:` directive."""
    params = getattr(smk, "params", None)
    log = getattr(smk, "log", None)
    log_file = str(log[0]) if log is not None and len(log) else None

    query = _rule_value(params, "query")
    verbose = bool(_rule_value(params, "verbose", False))
    configure_logging(_verbosity(verbose, query), log_file=log_file,
                      console=log_file is None)

    try:
        input_fasta = _rule_input(smk)
        reference = _rule_value(params, "reference")
        if reference is None:
            raise CodonifyError("the rule is missing `params: reference=...` "
                                "(the name of the reference sequence in the alignment)")
        output_base = _rule_output_base(smk)
        _check_declared_outputs(smk, output_base)
        mask_downstream, mask_upstream, cterm_mask = _rule_masking(smk)

        logger.info("codonify: %s -> %s.* (reference %s, -m %s, -u %s, -c %s)",
                    input_fasta, output_base, reference,
                    mask_downstream, mask_upstream, cterm_mask)

        run(input_fasta, output_base, str(reference),
            query=None if query is None else str(query),
            mask_downstream=mask_downstream,
            mask_upstream=mask_upstream,
            cterm_mask=cterm_mask)
    except CodonifyError as exc:
        # Log it for the rule's log file, then re-raise so Snakemake fails the job.
        logger.error("%s", exc)
        raise


# `snakemake` only exists in the globals when a `script:` directive runs this file.
if "snakemake" in globals():
    main_snakemake(globals()["snakemake"])
elif __name__ == "__main__":
    main()
