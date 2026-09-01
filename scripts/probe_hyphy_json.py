#!/usr/bin/env python3
"""probe_hyphy_json.py -- map the likelihood/evidence structure of a HyPhy JSON.

Run it on a real aBSREL (or BUSTED / MEME) JSON to find every array that could
carry site-level likelihoods, and print its shape. Use the reported paths to
wire up the ER extraction, instead of guessing at key names.

    python probe_hyphy_json.py codon_alignments/<id>/<id>.g_tree.absrel.json
"""
import json
import sys

KEYWORDS = ("log likelihood", "loglikelihood", "site", "evidence", "ratio",
            "posterior", "rate distribution", "lrt", "p-value")


def shape(v, depth=0):
    """Describe a nested list/dict shape compactly, e.g. list[1] -> list[412] -> float."""
    if isinstance(v, list):
        if not v:
            return "list[0]"
        return f"list[{len(v)}] -> {shape(v[0], depth + 1)}"
    if isinstance(v, dict):
        if not v:
            return "dict{}"
        k = next(iter(v))
        return f"dict{{{len(v)} keys, e.g. {k!r}}} -> {shape(v[k], depth + 1)}"
    if isinstance(v, bool):
        return "bool"
    if isinstance(v, (int, float)):
        return type(v).__name__
    if isinstance(v, str):
        return f"str({v[:40]!r}{'...' if len(v) > 40 else ''})"
    return type(v).__name__


def interesting(path):
    low = path.lower()
    return any(k in low for k in KEYWORDS)


def walk(node, path="", out=None, max_depth=5):
    if out is None:
        out = []
    if path.count("/") > max_depth:
        return out
    if isinstance(node, dict):
        for k, v in node.items():
            p = f"{path}/{k}"
            if interesting(p):
                out.append((p, shape(v)))
            walk(v, p, out, max_depth)
    elif isinstance(node, list) and node and isinstance(node[0], (dict, list)):
        walk(node[0], f"{path}[0]", out, max_depth)
    return out


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    with open(path) as f:
        data = json.load(f)

    print(f"=== {path} ===\n")
    print("-- top-level keys --")
    for k in data:
        print(f"   {k!r}: {shape(data[k])}")

    print("\n-- paths matching likelihood / evidence / site keywords --")
    seen = set()
    for p, s in walk(data):
        if p not in seen:
            seen.add(p)
            print(f"   {p}\n       {s}")

    ba = data.get("branch attributes", {}).get("0", {})
    if ba:
        name = next(iter(ba))
        print(f"\n-- per-branch keys (branch {name!r}, {len(ba)} branches total) --")
        for k, v in ba[name].items():
            print(f"   {k!r}: {shape(v)}")

    n_sites = data.get("input", {}).get("number of sites")
    print(f"\n-- input/number of sites: {n_sites} --")
    print("   Any array whose innermost length equals that is a per-site vector.")
    print("   Two such arrays (one full model, one constrained) = site-level ERs.")


if __name__ == "__main__":
    main()
