# Filters a codon-based multiple sequence alignment in FASTA format in two steps:
# 1) Column filtering: remove codon columns where fewer than `--mincodon` fraction
#    of sequences have a valid (non-gap, non-ambiguous) codon.
# 2) Row filtering: remove sequences where fewer than `--minseq` fraction of the
#    remaining codon columns are valid, or whose retained alignment is shorter
#    than `--minaalen` amino acids.
# Optionally masks bad codons and stop codons in the output with NNN.

import logging
import math
import sys
import re
import numpy as np
import pyfastx

__author__ = "Ekaterina Osipova, 2020, adapted by Bernhard Bein, 2026"

## Logging
log = logging.getLogger(__name__)

## A "good" codon: exactly three unambiguous ATGC bases (case-insensitive)
GOOD_CODON_PATTERN = re.compile(r"^[ATGCatgc]{3}$")
## Compared upper-cased. GOOD_CODON_PATTERN accepts mixed case, so "Tag" and
## "tAA" passed is_good_codon and were never masked.
STOP_CODONS = frozenset({"TAG", "TGA", "TAA"})

## Read/Write input
def read_fasta(fasta_file: str) -> dict[str, str]:
    """
    Read a FASTA alignment into an ordered dict {name: sequence}.
    Validates that all sequences have the same length (required for an MSA).
    """
    fasta_dict: dict[str, str] = {}
    for name, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fasta_dict[name] = seq

    lengths = {len(s) for s in fasta_dict.values()}
    if len(lengths) != 1:
        raise ValueError(
            f"Sequences in {fasta_file!r} do not all have the same length "
            f"(found lengths: {sorted(lengths)}). Not a valid alignment."
        )

    aln_len = lengths.pop()
    n_seqs = len(fasta_dict)
    n_codons = aln_len // 3

    if aln_len % 3 != 0:
        log.warning(
            "Alignment length %d is not a multiple of 3. "
            "The last %d nucleotide(s) will be ignored.",
            aln_len,
            aln_len % 3,
        )

    log.info(
        "Read %d sequences, alignment length %d nt (%d codon columns).",
        n_seqs,
        aln_len,
        n_codons,
    )
    return fasta_dict


def write_alignment(
    fasta_dict: dict[str, str], names_to_keep: list[str], out_file: str | None = None
) -> None:
    """Write FASTA records for *names_to_keep* to *out_file* (or stdout if None)."""
    if out_file:
        with open(out_file, "w") as f:
            for name in names_to_keep:
                f.write(f">{name}\n")
                f.write(f"{fasta_dict[name]}\n")
    else:
        for name in names_to_keep:
            print(f">{name}")
            print(fasta_dict[name])

## Small helper functions
def is_good_codon(codon: str) -> bool:
    """Return True if *codon* is a valid, non-ambiguous, non-gap triplet."""
    return bool(GOOD_CODON_PATTERN.match(codon))


def seq_to_codons(seq: str) -> list[str]:
    """Split *seq* into a list of codon strings (triplets)."""
    return [seq[i : i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]


def codon_goodness_vector(seq: str) -> list[bool]:
    """Return a per-codon boolean list (True = good codon)."""
    return [is_good_codon(c) for c in seq_to_codons(seq)]


## Column Filtering 
def select_good_columns(
    fasta_dict: dict[str, str], min_good_fraction: float
) -> list[int]:
    """
    Return the list of codon-column indices where at least
    *min_good_fraction* of sequences carry a valid codon.

    Logs every dropped column with its actual fraction.
    """
    names = list(fasta_dict)
    n_seqs = len(names)

    # Build boolean matrix: rows = sequences, columns = codon positions
    bool_matrix = np.array(
        [codon_goodness_vector(fasta_dict[name]) for name in names],
        dtype=bool,
    )
    n_cols = bool_matrix.shape[1]
    ## ceil, not round: round(0.8 * 3) == 2, so an '80% coverage' filter
    ## admitted 67%; banker's rounding also made round(0.5 * 1) == 0.
    threshold = math.ceil(min_good_fraction * n_seqs)

    good_cols: list[int] = []
    dropped_cols: list[tuple[int, float]] = []

    for col_idx in range(n_cols):
        n_good = bool_matrix[:, col_idx].sum()
        fraction = n_good / n_seqs
        if n_good >= threshold:
            good_cols.append(col_idx)
        else:
            dropped_cols.append((col_idx, fraction))

    ## Summary
    log.info(
        "Column filtering (min good fraction = %.2f): "
        "keeping %d / %d codon columns, dropping %d.",
        min_good_fraction,
        len(good_cols),
        n_cols,
        len(dropped_cols),
    )

    if dropped_cols:
        for col_idx, frac in dropped_cols:
            log.debug(
                "  Dropped codon column %d: %.1f%% good codons (threshold %.1f%%).",
                col_idx,
                frac * 100,
                min_good_fraction * 100,
            )
        ## Log a compact summary at INFO level instead of one line per column
        dropped_indices = [str(c) for c, _ in dropped_cols]
        log.info("  Dropped codon column indices: %s", ", ".join(dropped_indices))

    return good_cols


## Apply column filtering + masking codons (masking premature stop codons probably obsolete due to codonification doing this already)
def apply_column_filter(
    fasta_dict: dict[str, str],
    good_columns: list[int],
    mask: bool,
) -> dict[str, str]:
    """
    Retain only *good_columns* in each sequence.
    If *mask* is True, replace each remaining bad codon and stop codon with NNN.
    """
    filtered: dict[str, str] = {}
    for name, seq in fasta_dict.items():
        codons = seq_to_codons(seq)
        kept = [codons[i] for i in good_columns]

        if mask:
            masked: list[str] = []
            for codon in kept:
                if codon.upper() in STOP_CODONS:
                    masked.append("NNN")
                elif not is_good_codon(codon):
                    masked.append("NNN")
                else:
                    masked.append(codon)
            kept = masked

        filtered[name] = "".join(kept)

    return filtered


## Filter out rows (seqs) that are below threshold in full codons
def select_good_sequences(
    fasta_dict: dict[str, str],
    min_seq_fraction: float,
    min_aa_len: int,
) -> list[str]:
    """
    Return names of sequences that pass both criteria:
      1. At least *min_seq_fraction* of codon columns are valid (good codons).
      2. Total number of codon columns >= *min_aa_len*.

    Logs every dropped sequence with the reason.
    """
    good_names: list[str] = []
    dropped: list[tuple[str, str]] = []

    for name, seq in fasta_dict.items():
        goodness = codon_goodness_vector(seq)
        n_codons = len(goodness)
        n_good = sum(goodness)
        fraction = n_good / n_codons if n_codons > 0 else 0.0
        threshold = math.ceil(min_seq_fraction * n_codons)

        fails_fraction = n_good < threshold
        fails_length = n_codons < min_aa_len

        if fails_fraction or fails_length:
            reasons = []
            if fails_fraction:
                reasons.append(
                    f"only {fraction:.1%} good codons "
                    f"({n_good}/{n_codons}), threshold {min_seq_fraction:.1%}"
                )
            if fails_length:
                reasons.append(
                    f"sequence too short: {n_codons} codons "
                    f"(< {min_aa_len} AA threshold)"
                )
            dropped.append((name, "; ".join(reasons)))
        else:
            good_names.append(name)

    log.info(
        "Row filtering (min seq fraction = %.2f, min AA len = %d): "
        "keeping %d / %d sequences, dropping %d.",
        min_seq_fraction,
        min_aa_len,
        len(good_names),
        len(fasta_dict),
        len(dropped),
    )
    for name, reason in dropped:
        log.info("  Dropped sequence %r: %s.", name, reason)

    return good_names


## CLI args (masked for snakemake executrion)
def main() -> None:

    in_alignment = snakemake.input[0]
    out_alignment = snakemake.output[0]
    ## Named, not positional: the rule supplies these as named params, so
    ## reordering that block would silently swap a threshold with the flag.
    mincodon = snakemake.params.mincodon
    minseq   = snakemake.params.minseq
    minaalen = snakemake.params.minaalen
    mask     = snakemake.params.mask  #'-mc 0.6 -ms 0.3 -ml 25 -m'
    logfile = snakemake.log[0]

    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="[%(levelname)s] %(message)s",
    )

    log.info("=== manual_filter_msa ===")
    log.info("Input alignment : %s", in_alignment)
    log.info("Min codon frac  : %.2f", mincodon)
    log.info("Min seq frac    : %.2f", minseq)
    log.info("Min AA length   : %d", minaalen)
    log.info("Mask bad codons : %s", mask)

    ## Step 0: Read alignment
    fasta_dict = read_fasta(in_alignment)

    ## Step 1: Column filtering
    good_columns = select_good_columns(fasta_dict, mincodon)

    ## Step 2: Apply column filter (and optionally mask)
    column_filtered = apply_column_filter(fasta_dict, good_columns, mask)

    ## Step 3: Row filtering
    good_sequences = select_good_sequences(column_filtered, minseq, minaalen)

    ## Output
    log.info(
        "Writing final alignment: %d sequences x %d codon columns to %s.",
        len(good_sequences),
        len(good_columns),
        out_alignment,
    )
    write_alignment(column_filtered, good_sequences, out_alignment)


if __name__ == "__main__":
    main()
