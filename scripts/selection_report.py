#!/usr/bin/env python3
"""One plain-text answer per transcript: what did the selection screens find?

Reads only what it is given, so it adapts to whichever screens were enabled --
a missing argument means that screen was off and its section is omitted. The
point is that the headline result stays readable on disk after HyPhy_output/
and the JSONs have been archived away.

    selection_report.py --transcript ENST1 --out ENST1.selection.txt \
        [--status validation.txt] \
        [--absrel-tsv F] [--busted-json F] [--meme-tsv F] [--relax-tsv F] \
        [--pvalue 0.05]
"""
import argparse
import csv
import json
import os
import statistics
import sys

WIDTH = 66

def write_report(path, lines):
    """LF endings, so the file reads the same from the cluster and from Windows."""
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

## From validate
def read_verdict(path):
    """validate_alignment's verdict: (verdict, reason), or (None, None).
    """
    if not path or not os.path.exists(path):
        return None, None
    try:
        fields = open(path, encoding="utf-8").read().strip().split("\t")
    except OSError as exc:
        print(f"warning: cannot read {path}: {exc}", file=sys.stderr)
        return None, None
    if not fields or not fields[0]:
        return None, None
    return fields[0].strip(), (fields[1].strip() if len(fields) > 1 else "")


def read_tsv(path):
    """-> (fieldnames, rows). ([], []) if unusable, with a note on stderr."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return [], []
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            return (r.fieldnames or []), list(r)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        print(f"warning: cannot read {path}: {exc}", file=sys.stderr)
        return [], []


## Why this?
def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def find_column(fields, *wanted, exclude=()):
    """First column containing all of *wanted and none of *exclude*.

    `exclude` because, for example "Branch Uncorrected P-value" contains
    "corrected" as a substring of "uncorrected", so a plain match for
    ("corrected", "p-value") silently picked the UNcorrected column and
    under-reported the branches under selection.
    """
    for f in fields or []:
        low = f.lower()
        if any(x.lower() in low for x in exclude):
            continue
        if all(w.lower() in low for w in wanted):
            return f
    return None


## ── aBSREL: which branches are under selection ──────────────────────────────
def absrel_section(path, alpha):
    fields, rows = read_tsv(path)
    if not rows:
        return ["aBSREL", "  no usable output at " + str(path)]
    col_p = (find_column(fields, "corrected", "p-value", exclude=("uncorrected",))
             or find_column(fields, "p-value", exclude=("uncorrected",)))
    if not col_p:
        return ["aBSREL", f"  no p-value column in {os.path.basename(path)}"
                          f" (columns: {', '.join(fields[:6])})"]
    ## The table is per branch per site, so one p-value repeats many times.
    per_branch = {}
    for row in rows:
        b = row.get("Branch")
        p = num(row.get(col_p))
        if b and p is not None:
            per_branch[b] = p
    if not per_branch:
        return ["aBSREL", "  no branch p-values found"]
    hits = sorted((p, b) for b, p in per_branch.items() if p <= alpha)
    out = ["aBSREL",
           f"  branches tested        : {len(per_branch)}",
           f"  under selection (p<={alpha}): {len(hits)}"]
    if hits:
        out.append("  " + ", ".join(f"{b} (p={p:.3g})" for p, b in hits[:12]))
        if len(hits) > 12:
            out.append(f"  ... and {len(hits) - 12} more")
    else:
        out.append("  -> no branch under selection")
    return out


## ── BUSTED: was selection detected at all ───────────────────────────────────
def busted_section(path, alpha):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ["BUSTED", f"  no usable output at {path}"]
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError) as exc:
        return ["BUSTED", f"  cannot read {os.path.basename(path)}: {exc}"]
    tr = data.get("test results") or {}
    p = num(tr.get("p-value"))
    lrt = num(tr.get("LRT"))
    out = ["BUSTED"]
    if p is None:
        out.append("  no 'test results' p-value in the JSON")
    else:
        verdict = "SELECTION DETECTED" if p <= alpha else "no evidence of selection"
        out.append(f"  p-value                : {p:.3g}"
                   + (f"   (LRT {lrt:.4g})" if lrt is not None else ""))
        out.append(f"  -> {verdict} (p<={alpha})")
    ## Report the omega classes too.
    try:
        dist = data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
        pos = [(num(d.get("omega")), num(d.get("proportion")))
               for d in dist.values()]
        pos = [(o, w) for o, w in pos if o is not None and w is not None and o > 1]
        if pos:
            out.append("  omega>1 classes        : "
                       + ", ".join(f"omega={o:.3g} (weight {w:.3g})" for o, w in pos))
    except (KeyError, TypeError, AttributeError):
        pass
    return out


## ── MEME: which sites are under episodic selection ──────────────────────────
def meme_section(path, alpha):
    fields, rows = read_tsv(path)
    if not rows:
        return ["MEME", f"  no usable output at {path}"]
    ## MEME's column names come from the JSON, so match rather than assume.
    col_p = find_column(fields, "p-value")
    if not col_p:
        return ["MEME", f"  no p-value column in {os.path.basename(path)}"
                        f" (columns: {', '.join(fields[:8])})"]
    #col_site = (find_column(fields, "codon") or find_column(fields, "site")
    #            or fields[0])
    hits = []
    for i, row in enumerate(rows, start=1):
        p = num(row.get(col_p))
        if p is not None and p <= alpha:
            site = i
            hits.append((site, p))
    out = ["MEME",
           f"  sites tested           : {len(rows)}",
           f"  under selection (p<={alpha}): {len(hits)}"]
    if hits:
        out.append("  sites: " + ", ".join(str(s) for s, _ in hits[:25]))
        if len(hits) > 25:
            out.append(f"  ... and {len(hits) - 25} more")
    else:
        out.append("  -> no site under episodic selection")
    return out


## ── RELAX: relaxed or intensified, and the spread of k ──────────────────────
def relax_section(path, alpha):
    fields, rows = read_tsv(path)
    if not rows:
        return ["RELAX", f"  no usable output at {path}"]
    col_k = find_column(fields, "k") or "k"
    col_p = find_column(fields, "p-value")
    fg = {}
    for row in rows:
        b = row.get("Branch")
        p = num(row.get(col_p)) if col_p else None
        k = num(row.get(col_k))
        if not b or p is None:
            continue
        fg.setdefault(b, {"k": [], "p": []})
        if k is not None:
            fg[b]["k"].append(k)
        fg[b]["p"].append(p)
    if not fg:
        return ["RELAX",
                "  no foreground branch carried a p-value; RELAX reports the",
                "  test only for branches labelled 'Test'"]
    out = ["RELAX", f"  foreground branches    : {len(fg)}"]
    ## k > 1 is intensified selection, k < 1 relaxed (RELAX's own convention).
    for b, v in sorted(fg.items()):
        ks, ps = v["k"], v["p"]
        p_med = statistics.median(ps) if ps else None
        if not ks:
            out.append(f"  {b}: no k values")
            continue
        mean_k = statistics.fmean(ks)
        direction = "RELAXED" if mean_k < 1 else "INTENSIFIED"
        sig = ("significant" if p_med is not None and p_med <= alpha
               else "not significant")
        out.append(f"  {b}: {direction} ({sig}, median p="
                   f"{'n/a' if p_med is None else format(p_med, '.3g')})")
        out.append(f"      k over {len(ks)} run(s): mean {mean_k:.4g}, "
                   f"min {min(ks):.4g}, max {max(ks):.4g}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--transcript", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--status", help="validate_alignment's validation.txt")
    ap.add_argument("--absrel-tsv")
    ap.add_argument("--busted-json")
    ap.add_argument("--meme-tsv")
    ap.add_argument("--relax-tsv")
    ap.add_argument("--pvalue", type=float, default=0.05)
    args = ap.parse_args()

    lines = [f"Selection screen summary: {args.transcript}", "=" * WIDTH, ""]

    verdict, reason = read_verdict(args.status)
    if verdict is not None and verdict != "OK":
        lines.append("No screen was run: the alignment validation checkpoint")
        lines.append(f"rejected this transcript (verdict {verdict}).")
        lines.append(f"  reason: {reason or 'not recorded'}")
        lines.append("")
        lines.append("-" * WIDTH)
        lines.append("See the cleaning report for how much cleaning removed.")
        write_report(args.out, lines)
        print(f"wrote {args.out} (transcript {verdict})")
        return

    sections = [
        (args.absrel_tsv, absrel_section),
        (args.busted_json, busted_section),
        (args.meme_tsv, meme_section),
        (args.relax_tsv, relax_section),
    ]
    ran = 0
    for path, fn in sections:
        if not path:
            continue          # screen not enabled for this run
        ran += 1
        lines.extend(fn(path, args.pvalue))
        lines.append("")

    if not ran:
        lines.append("No selection screen was enabled for this run.")
        lines.append("")
    lines.append("-" * WIDTH)
    lines.append(f"significance threshold: p <= {args.pvalue}")
    lines.append("Full HyPhy tables and JSONs are in results.tar.gz.")

    write_report(args.out, lines)
    print(f"wrote {args.out}")


main()
