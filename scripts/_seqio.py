"""Minimal FASTA I/O, shared by the pipeline scripts.

Deliberately stdlib-only and dependency-free. Snakemake puts a `script:` directive's own directory on
sys.path, so a sibling script can do `from _seqio import read_fasta`.
"""

import gzip
import io
import os

FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fas")


def _open_text(path):
    """Open a plain or gzipped file in text mode."""
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def iter_fasta(path, full_header=True):
    """Yield (header, sequence) in file order.

    full_header=True keeps the whole header line, which is what the tree code
    needs. full_header=False truncates at the first whitespace, matching
    pyfastx's notion of a sequence name.
    """
    header, chunks = None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:] if full_header else line[1:].split(None, 1)[0]
                chunks = []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def read_fasta(path, full_header=True):
    """Ordered {header: sequence}. Raises on a duplicate header."""
    out = {}
    for header, seq in iter_fasta(path, full_header=full_header):
        if header in out:
            raise ValueError(f"{path}: duplicate FASTA header {header!r}")
        out[header] = seq
    return out


def read_names(path, full_header=True):
    """Just the headers, in file order."""
    return [h for h, _ in iter_fasta(path, full_header=full_header)]


def alignment_stats(path):
    """(n_sequences, set_of_lengths, n_all_gap_columns) for an alignment."""
    lengths, seqs = set(), []
    for _, seq in iter_fasta(path):
        lengths.add(len(seq))
        seqs.append(seq)
    if len(lengths) != 1 or not seqs:
        return len(seqs), lengths, None
    width = lengths.copy().pop()
    empty_cols = sum(
        1 for i in range(width)
        if all(s[i] in "-Nn?*" for s in seqs)
    )
    return len(seqs), lengths, empty_cols


def is_fasta(path):
    """True if the first non-blank line starts with '>'. Handles .gz."""
    try:
        with _open_text(path) as fh:
            for line in fh:
                if line.strip():
                    return line.lstrip().startswith(">")
    except (OSError, UnicodeDecodeError, EOFError):
        return False
    return False


def has_fasta_extension(filename):
    """True for .fa/.fasta/.fna/.fas, optionally .gz-suffixed."""
    name = filename[:-3] if filename.endswith(".gz") else filename
    return name.endswith(FASTA_EXTENSIONS)


def write_fasta(path, records, width=None):
    """Write {header: seq} or [(header, seq)] to *path*."""
    items = records.items() if hasattr(records, "items") else records
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for header, seq in items:
            fh.write(f">{header}\n")
            if width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")
