#!/usr/bin/env python3
# Masks the trailing (terminal) stop codon of every sequence in a codon-based
# multiple sequence alignment in FASTA format.

import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _seqio import iter_fasta  # noqa: E402

__author__ = "Bernhard Bein, 2026"

## Logging
log = logging.getLogger(__name__)
LOG_LEVEL = logging.INFO

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
## Characters that mean "this sequence has nothing here". MACSE's frameshift
## marker `!` is deliberately NOT a gap: it is real, if broken, data.
GAP_CHARS = frozenset("-.~?")


## Read/Write input
def read_alignment(fasta_file: str) -> dict[str, str]:
    """
    Read a FASTA alignment into an ordered dict {name: sequence}.
    Validates that all sequences have the same length (required for an MSA).

    Names are truncated at the first whitespace. .gz input is handled by _seqio.
    """
    fasta_dict: dict[str, str] = {}
    for name, seq in iter_fasta(fasta_file, full_header=False):
        if name in fasta_dict:
            raise ValueError(f"{fasta_file!r}: duplicate FASTA header {name!r}")
        fasta_dict[name] = seq

    if not fasta_dict:
        raise ValueError(f"No sequences found in {fasta_file!r}.")

    lengths = {len(s) for s in fasta_dict.values()}
    if len(lengths) != 1:
        raise ValueError(
            f"Sequences in {fasta_file!r} do not all have the same length "
            f"(found lengths: {sorted(lengths)}). Not a valid alignment."
        )

    aln_len = lengths.pop()
    if aln_len % 3 != 0:
        log.warning(
            "Alignment length %d is not a multiple of 3. The trailing %d "
            "nucleotide(s) lie outside any codon and are copied through unchanged.",
            aln_len,
            aln_len % 3,
        )

    log.info(
        "Read %d sequences, alignment length %d nt (%d codon columns).",
        len(fasta_dict),
        aln_len,
        aln_len // 3,
    )
    return fasta_dict


def write_alignment(fasta_dict: dict[str, str], out_file: str | None = None) -> None:
    """Write every record of *fasta_dict* to *out_file* (or stdout if None)."""
    if out_file is None:
        for name, seq in fasta_dict.items():
            print(f">{name}")
            print(seq)
        return
    with open(out_file, "w") as f:
        for name, seq in fasta_dict.items():
            f.write(f">{name}\n{seq}\n")


## Small helper functions
def is_gap_codon(codon: str) -> bool:
    """Return True if *codon* carries no data at all (only gap characters)."""
    return not (set(codon) - GAP_CHARS)


def terminal_codon_index(codons: list[str]) -> int | None:
    """
    Return the index of the last codon that is not entirely gaps, or None.

    Scanning back from the 3-prime end rather than taking the last codon column
    is what makes this correct for an alignment with ragged ends: a sequence
    that stops early still gets *its own* last codon inspected, not the run of
    trailing gaps that pads it out to the full alignment width.
    """
    for i in range(len(codons) - 1, -1, -1):
        if not is_gap_codon(codons[i]):
            return i
    return None


## The single transformation this script performs
def mask_terminal_stop(seq: str, mask_with: str = "NNN") -> tuple[str, int | None]:
    """
    Return (sequence, masked_codon_index).

    If the sequence's terminal codon is a stop codon it is replaced by
    *mask_with* and its codon index is returned; otherwise the sequence comes
    back unchanged with None. The sequence length never changes, and trailing
    1-2 nt that do not form a full codon are copied through verbatim.
    """
    n_full = len(seq) - len(seq) % 3
    codons = [seq[i: i + 3] for i in range(0, n_full, 3)]
    tail = seq[n_full:]

    idx = terminal_codon_index(codons)
    if idx is None or codons[idx].upper() not in STOP_CODONS:
        return seq, None

    codons[idx] = mask_with
    return "".join(codons) + tail, idx


def mask_alignment(fasta_dict: dict[str, str], mask_with: str = "NNN") -> dict[str, str]:
    """Apply mask_terminal_stop to every sequence and log what was changed."""
    masked: dict[str, str] = {}
    n_masked = 0

    for name, seq in fasta_dict.items():
        new_seq, idx = mask_terminal_stop(seq, mask_with)
        masked[name] = new_seq
        if idx is None:
            log.debug("  %s: terminal codon is not a stop codon, left unchanged.", name)
            continue
        n_masked += 1
        log.debug(
            "  %s: masked terminal stop codon %s at codon column %d -> %s.",
            name,
            seq[idx * 3: idx * 3 + 3].upper(),
            idx,
            mask_with,
        )

    log.info(
        "Masked a terminal stop codon in %d / %d sequences (mask = %s).",
        n_masked,
        len(fasta_dict),
        mask_with,
    )
    return masked


## Entry points
def run(in_alignment: str, out_alignment: str | None, mask_with: str = "NNN") -> None:
    log.info("=== mask_terminal_stops ===")
    log.info("Input alignment : %s", in_alignment)
    log.info("Mask with       : %s", mask_with)

    fasta_dict = read_alignment(in_alignment)
    masked = mask_alignment(fasta_dict, mask_with)

    log.info(
        "Writing alignment: %d sequences x %d nt to %s.",
        len(masked),
        len(next(iter(masked.values()))),
        out_alignment or "<stdout>",
    )
    write_alignment(masked, out_alignment)


## CLI args (bypassed for snakemake execution)
def main_snakemake(smk) -> None:
    """Configured from the rule; `mask_with` is optional and defaults to NNN."""
    logging.basicConfig(
        filename=smk.log[0],
        filemode="w",
        level=LOG_LEVEL,
        format="[%(levelname)s] %(message)s",
    )
    run(smk.input[0], smk.output[0], getattr(smk.params, "mask_with", "NNN"))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Mask the trailing (terminal) stop codon of every sequence "
                    "in a codon alignment. Nothing else is changed."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="input codon alignment (FASTA, optionally .gz)")
    parser.add_argument("-o", "--output", default=None,
                        help="output FASTA (default: stdout)")
    parser.add_argument("--mask-with", default="NNN", metavar="CODON",
                        help="replacement for the terminal stop codon "
                             "(default: NNN; use --- to gap it out instead)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="log one line per sequence")
    args = parser.parse_args()

    if len(args.mask_with) != 3:
        parser.error(f"--mask-with must be exactly 3 characters, got {args.mask_with!r}")

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.DEBUG if args.verbose else LOG_LEVEL,
        format="[%(levelname)s] %(message)s",
    )
    try:
        run(args.input, args.output, args.mask_with)
    except (OSError, ValueError) as exc:
        log.error("%s", exc)
        sys.exit(1)


if __name__ == "__main__":
    if "snakemake" in globals():
        main_snakemake(globals()["snakemake"])
    else:
        main()
