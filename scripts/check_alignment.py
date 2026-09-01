#!/usr/bin/env python3
"""Gate a transcript's alignment before anything expensive depends on it.
Writes a single line to snakemake.output.status:
    OK
    SKIP<TAB><reason>
and always exits 0. 
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _seqio import iter_fasta  # noqa: E402

GAP_CHARS = set("-.?nN*")

def check(path, min_taxa=4, foreground=(), min_ungapped_fraction=0.1):
    """Return (verdict, reason). verdict is 'OK' or 'SKIP'."""
    if not os.path.exists(path):
        return "SKIP", f"alignment missing: {path}"
    if os.path.getsize(path) == 0:
        return "SKIP", "alignment file is empty (0 bytes)"

    try:
        records = list(iter_fasta(path))
    except (OSError, UnicodeDecodeError) as exc:
        return "SKIP", f"alignment unreadable: {exc}"

    if not records:
        return "SKIP", "no FASTA records"

    names = [h.split()[0] if h.split() else h for h, _ in records]
    seqs = [s for _, s in records]

    if len(set(names)) != len(names):
        dup = sorted({n for n in names if names.count(n) > 1})
        return "SKIP", f"duplicate sequence names: {', '.join(dup[:5])}"

    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        return "SKIP", (f"not aligned: {len(lengths)} distinct sequence lengths "
                        f"({min(lengths)}-{max(lengths)})")

    width = lengths.pop()
    if width == 0:
        return "SKIP", "alignment width is 0"
    if width % 3 != 0:
        return "SKIP", f"alignment width {width} is not a multiple of 3"

    ## Sequences that are effectively all gaps carry no information and make
    ## downstream tree estimation and HyPhy fits fail.
    informative = [
        n for n, s in zip(names, seqs)
        if sum(1 for c in s if c not in GAP_CHARS) >= max(3, width * min_ungapped_fraction)
    ]
    if len(informative) < min_taxa:
        return "SKIP", (f"only {len(informative)} informative sequence(s) of "
                        f"{len(seqs)}, need {min_taxa}")

    if all(all(c in GAP_CHARS for c in col) for col in zip(*seqs)):
        return "SKIP", "every alignment column is gaps or N"

    ## RELAX is run with --test Foreground and hard-fails if no branch carries
    ## the label, so a transcript that lost all its foreground taxa during
    ## cleaning has to be skipped rather than attempted.
    if foreground:
        present = [f for f in foreground
                   if any(f == n or n.startswith(f) for n in informative)]
        if not present:
            return "SKIP", "no foreground species left after cleaning"

    return "OK", ""


def main():
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        path = smk.input.ali
        min_taxa = int(smk.params.min_taxa)
        foreground = list(smk.params.foreground or [])
        out_path = smk.output.status
        log_path = smk.log[0] if smk.log else None
    ## CLI for debugging
    else:
        if len(sys.argv) < 2:
            sys.exit(__doc__)
        path = sys.argv[1]
        min_taxa = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        foreground = sys.argv[3:]
        out_path = None
        log_path = None

    ## Any unexpected exception becomes a SKIP verdict rather than a failed
    ## job. check() already returns SKIP for every malformed input it
    ## anticipates; this covers the ones it does not. 
    try:
        verdict, reason = check(path, min_taxa=min_taxa, foreground=foreground)
    except Exception as exc:                                  # noqa: BLE001
        verdict = "SKIP"
        reason = f"validation error ({type(exc).__name__}): {exc}"
        print(f"unexpected error validating {path}: {exc!r}", file=sys.stderr)
    line = verdict if verdict == "OK" else f"{verdict}\t{reason}"

    message = (f"{os.path.basename(path)}: {verdict}"
               + (f" - {reason}" if reason else ""))
    if log_path:
        with open(log_path, "w") as fh:
            fh.write(message + "\n")
    print(message, file=sys.stderr)

    if out_path:
        with open(out_path, "w", newline="\n") as fh:
            fh.write(line + "\n")
    ## Always exit 0: a skip is a recorded outcome, not a rule failure.

main()
