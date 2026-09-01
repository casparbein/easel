#!/usr/bin/env python3
"""Report what actually differs between the cleaning stages of one transcript.

Answers the three things that can make a stage look 'empty' in
plot_cleaning.R: mismatched dimensions, mismatched sequence names, or a mask
character the plot script does not recognise.

    python3 diag_cleaning.py [DIR] [ID]
"""
import os
import sys
from collections import Counter

## Same chain as plot_cleaning.R: .premasked.fa is a codonify internal and is
## deliberately not a stage here.
STAGES = [
    ("ori",         "_ori.fa"),
    ("ren",         "_ren.fa"),
    ("masked",      ".masked.fa"),
    ("hmm_cleaned", ".masked.hmm_cleaned.fa"),
    ("manual",      ".manual.fa"),
]
## Longest first so ".masked.fa" never shadows ".masked.hmm_cleaned.fa".
SUFFIXES = [".masked.hmm_cleaned.fa", ".masked.fa", ".manual.fa",
            "_ren.fa", "_ori.fa"]


def read_fasta(path):
    """-> dict name -> seq, names whitespace-trimmed like the R script does."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0] if line[1:].split() else line[1:], []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def infer_id(d):
    files = os.listdir(d)
    for sfx in SUFFIXES:
        hits = [f for f in files if f.endswith(sfx)]
        if hits:
            return hits[0][: len(hits[0]) - len(sfx)]
    return None


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "."
    tid = sys.argv[2] if len(sys.argv) > 2 else infer_id(d)
    if not tid:
        sys.exit(f"no alignment files found in {d}")
    print(f"transcript: {tid}\ndirectory : {os.path.abspath(d)}\n")

    loaded = {}
    for label, sfx in STAGES:
        p = os.path.join(d, tid + sfx)
        if not os.path.exists(p):
            print(f"{label:12} MISSING            ({os.path.basename(p)})")
            continue
        if os.path.getsize(p) == 0:
            print(f"{label:12} EMPTY FILE         ({os.path.basename(p)})")
            continue
        s = read_fasta(p)
        widths = {len(v) for v in s.values()}
        loaded[label] = s
        print(f"{label:12} {len(s):3d} seqs x {sorted(widths)} cols"
              f"{'   <-- NOT ALIGNED' if len(widths) != 1 else ''}")
        alpha = Counter("".join(s.values()))
        top = ", ".join(f"{c!r}:{n}" for c, n in alpha.most_common(12))
        print(f"{'':12} alphabet: {top}")
    print()

    order = [l for l, _ in STAGES if l in loaded]
    for a, b in zip(order, order[1:]):
        A, B = loaded[a], loaded[b]
        print(f"--- {a} -> {b} " + "-" * 40)
        na, nb = set(A), set(B)
        if na != nb:
            print(f"  names differ: {len(na - nb)} only in {a}, {len(nb - na)} only in {b}")
            for x in list(na - nb)[:3]:
                print(f"    only in {a}: {x!r}")
            for x in list(nb - na)[:3]:
                print(f"    only in {b}: {x!r}")
        shared = sorted(na & nb)
        print(f"  shared sequences: {len(shared)}")
        wa = {len(A[k]) for k in shared}
        wb = {len(B[k]) for k in shared}
        if wa == wb:
            note = ""
        elif max(wb) < max(wa):
            ## plot_cleaning.R maps these back with a two-pointer walk and
            ## still compares the surviving columns for masking.
            note = "   <-- width SHRANK: read as column removal"
        else:
            note = "   <-- width GREW: cells cannot be mapped onto the reference grid"
        print(f"  widths: {a}={sorted(wa)}  {b}={sorted(wb)}{note}")
        if wa == wb and shared:
            trans = Counter()
            for k in shared:
                for x, y in zip(A[k], B[k]):
                    if x != y:
                        trans[(x, y)] += 1
            total = sum(trans.values())
            print(f"  cells changed: {total}")
            if total:
                for (x, y), n in trans.most_common(10):
                    print(f"    {x!r} -> {y!r}  x{n}")
            else:
                print("    (identical content -- nothing for this stage to show)")
        print()


main()
