#!/usr/bin/env python3

"""
extract_hyphy.py – Combined extractor for HyPhy JSON result files.

Sub-commands
------------
  busted   Extract site-level evidence ratios/log-likelihoods, per-branch/site
           posteriors/EBF, and the MG94 fitted tree from a BUSTED JSON.
           -> CHECK what exactly is extracted!!
  meme     Extract per-branch/site posteriors/EBF AND the MLE table from a MEME JSON
           (two output TSVs).
           -> CHECK what exactly is extracted!!
  absrel   Extract per-branch/site posteriors/EBF and a fitted tree from an aBSREL JSON.
           -> CHECK what exactly is extracted!!
  relax    Extract per-branch k and p-Values of one to several relax runs
           -> CHECK what exactly is extracted!!

Usage examples
--------------
  python extract_hyphy.py busted results.busted.json
  python extract_hyphy.py meme   results.meme.json
  python extract_hyphy.py absrel results.absrel.json
  python extract_hyphy.py relax results.relax1.json,results.relax2.json

Output file names are derived from the input JSON name unless overridden with options.
Run any sub-command with --help for details.
"""

import json
import csv
import glob
import os
import re
import sys
import argparse
import math


## Clamp for posterior probabilities at 0 or 1, so the odds stay finite.
EPS = 1e-10

# ---------------------------------------------------------------------------
# Tree-parser helpers
# ---------------------------------------------------------------------------

def parse_newick_proper(newick_str):
    """Parse a Newick string and return a child->parent mapping."""
    newick_str = re.sub(r':[0-9\.eE+-]+', '', newick_str).rstrip(';')
    tokens = re.findall(r'\(|\)|,|[^(),]+', newick_str)

    parent_map = {}
    stack = []
    internal_node_counter = [0]

    class Node:
        def __init__(self, name=None):
            self.name = name
            self.children = []

    root = Node("root")
    current = root

    for token in tokens:
        if token == '(':
            new_node = Node()
            current.children.append(new_node)
            stack.append(current)
            current = new_node
        elif token == ',':
            new_node = Node()
            current = stack[-1]
            current.children.append(new_node)
            current = new_node
        elif token == ')':
            current = stack.pop()
        else:
            current.name = token

    def traverse(node, parent_name):
        if not node.name:
            node.name = f"Node_{internal_node_counter[0]}"
            internal_node_counter[0] += 1
        if parent_name:
            parent_map[node.name] = parent_name
        for child in node.children:
            traverse(child, node.name)

    traverse(root.children[0] if len(root.children) == 1 else root, None)
    return parent_map


def propagate_codons_for_site(parent_map, site_subs):
    """Propagate codon substitutions down the tree for a single site."""
    children_map = {}
    for child, parent in parent_map.items():
        children_map.setdefault(parent, []).append(child)

    all_nodes = set(parent_map.keys()) | set(parent_map.values())
    roots = all_nodes - set(parent_map.keys())
    if not roots:
        return {}
    root = list(roots)[0]

    node_to_codon = {}

    def dfs(node, current_codon):
        node_codon = site_subs.get(node, current_codon)
        node_to_codon[node] = node_codon
        for child in children_map.get(node, []):
            dfs(child, node_codon)

    root_codon = site_subs.get(root, site_subs.get('root', '---'))
    dfs(root, root_codon)
    return node_to_codon


# ---------------------------------------------------------------------------
# Evidence ratios and empirical Bayes factors
#
#   ER (Evidence Ratio) - a LIKELIHOOD ratio between two *models* for one site or one branch:
#         exp(lnL_alt - lnL_null), equivalently exp(LRT/2).
#         Dimension: site (BUSTED, MEME) or branch x site (aBSREL).
#
#   EBF (Empirical Bayes Factor) - a POSTERIOR-odds / PRIOR-odds ratio for the omega>1 rate class
#         *within a single* selection-allowing model.
#         Dimension: branch x site.
#
# ---------------------------------------------------------------------------

def evidence_ratio(ll_alt, ll_null):
    """ER = exp(lnL_alt - lnL_null), evaluated as a difference of logs.
    """
    if ll_alt is None or ll_null is None:
        return None
    try:
        return math.exp(float(ll_alt) - float(ll_null))
    except (TypeError, ValueError, OverflowError):
        return None

# def er_from_lrt(lrt):
#     """ER implied by a likelihood-ratio statistic, since LRT = 2*(lnL_a - lnL_0)."""
#     if lrt is None:
#         return None
#     try:
#         return math.exp(float(lrt) / 2.0)
#     except (TypeError, ValueError, OverflowError):
#         return None


def empirical_bayes_factor(post_prob, prior_prob):
    """EBF = posterior odds / prior odds for the omega>1 class.

    Returns (ebf, log10_ebf), or (None, None) when the EBF is undefined:
      * prior == 0  -> the fitted model has no omega>1 class here, so there is
                       nothing to shift belief toward. That is NOT "EBF = 0".
      * prior == 1  -> degenerate.
      * prior outside [0, 1] -> the caller passed something that is not a
                       probability.
    """
    if post_prob is None or prior_prob is None:
        return None, None
    try:
        prior_prob = float(prior_prob)
        post_prob = float(post_prob)
    except (TypeError, ValueError):
        return None, None
    if not (0.0 < prior_prob < 1.0):
        return None, None
    p = max(min(post_prob, 1.0 - EPS), EPS)
    ebf = (p / (1.0 - p)) / (prior_prob / (1.0 - prior_prob))
    return ebf, (math.log10(ebf) if ebf > 0 else None)


def _blank(value):
    """Empty cell for a missing value, so None never prints as the text 'None'."""
    return "" if value is None else value


def _unwrap(arr):
    """HyPhy wraps per-site vectors in a one-element outer list."""
    if isinstance(arr, list) and arr and isinstance(arr[0], list):
        return arr[0]
    return arr


def _at(arr, i):
    """arr[i], unwrapping a one-element list, or None if out of range."""
    if arr is None or i >= len(arr):
        return None
    v = arr[i]
    return v[0] if isinstance(v, list) else v


## Check whether this function is needed
def resolve_site_ll(data):
    """Locate the site log-likelihood vectors, whatever layout this JSON uses.

    Returns (global_vectors, per_branch_vectors):
    "Site Log Likelihood": {              <- aBSREL
          "unconstrained": [[...]],             the alternative, global
          "tested": {"<branch>": [[...]]}       each TESTED branch's own null
      }

      branch attributes -> 0 -> <branch> -> "Site Log Likelihood"
                                                  some other HyPhy layouts
    Key names differ between HyPhy versions, so this discovers rather than
    hardcodes.
    """
    glob, per_branch = {}, {}
    for key, value in (data.get("Site Log Likelihood") or {}).items():
        if isinstance(value, list):
            glob[key] = _unwrap(value)
        elif isinstance(value, dict):
            ## Keyed by branch name, one vector each.
            for branch, vec in value.items():
                if isinstance(vec, list):
                    per_branch.setdefault(branch, {})[key] = _unwrap(vec)

    for branch, attrs in (data.get("branch attributes", {}).get("0", {}) or {}).items():
        sll = attrs.get("Site Log Likelihood")
        if isinstance(sll, dict):
            per_branch.setdefault(branch, {}).update(
                {k: _unwrap(v) for k, v in sll.items() if isinstance(v, list)})
        elif isinstance(sll, list):
            per_branch.setdefault(branch, {})["constrained"] = _unwrap(sll)
    return glob, per_branch


def absrel_site_er(data):
    """Return er(branch, site_idx): full adaptive model vs that branch's null.

    aBSREL's omega>1-vs-omega<=1 contrast is defined per branch, so its site ERs
    are branch x site, and the null has to be that branch's own vector.
    """
    glob, per_branch = resolve_site_ll(data)
    alt = glob.get("unconstrained") or glob.get("Full adaptive model")

    def er(branch, site_idx):
        b = per_branch.get(branch) or {}
        null = (b.get("tested")                 # aBSREL, this branch's own null
                or b.get("constrained")
                or b.get("Baseline MG94xREV")
                or next(iter(b.values()), None))
        if null is None:
            ## No branch-specific null. aBSREL emits one only for the branches
            ## it actually tested, so for the rest the contrast is undefined and
            ## the cell is blank.
            return None
        return evidence_ratio(_at(alt, site_idx), _at(null, site_idx))
    return er


def error_sink_enabled(data):
    """True when the HyPhy run was given --error-sink.

    Read from the JSON: HyPhy records its own settings at
    analysis -> settings. BUSTED v4.7 writes "error-sink": 1 there.
    """
    settings = (data.get("analysis") or {}).get("settings") or {}
    flag = settings.get("error-sink")
    if flag is None:
        return False
    return str(flag).strip().lower() not in ("0", "no", "false", "")

def busted_error_sink_index(data):
    """Index of the error-sink omega class, or None when the run had no sink.
    analysis -> settings -> "error-sink" says whether --error-sink was on.

    Which class is the sink is then determined by its omega, which the model
    fixes unrealistically high so that alignment error collects there instead
    of inflating the genuine positive-selection class. 

    Returns (index, note).
    """
    if not error_sink_enabled(data):
        return None, None

    try:
        dist = (data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"])
        omegas = {int(k): float(v["omega"]) for k, v in dist.items()}
    except (KeyError, TypeError, ValueError) as exc:
        return None, (f"error-sink is on but the omega classes are unreadable "
                      f"({exc}); the sink cannot be separated")
    if not omegas:
        return None, "error-sink is on but there are no omega classes"

    sink = max(omegas, key=lambda k: omegas[k])
    note = None
    if sink != 0:
        note = (f"error-sink is on and the largest omega is class {sink} "
                f"(omega {omegas[sink]:.4g}), not class 0 as expected; "
                f"treating class {sink} as the sink")
    return sink, note

def codon_map_by_site(parent_map, substitutions, num_sites):
    """{site_idx: {node: codon}}, computed once per site.
    """
    return {
        s: propagate_codons_for_site(parent_map, substitutions.get(str(s)) or {})
        for s in range(num_sites)
    }


# ===========================================================================
# BUSTED – site-level
# ===========================================================================

def extract_busted_to_csv(json_path: str, csv_path: str):
    """Extract site-level evidence ratios and log-likelihoods from a BUSTED JSON.
    In current (2026, HyPhy v.2.5.101) BUSTED json file, Evidence Ratios can be read 
    out from the file directly, but computation based on site LogL exp(LogL-alt - LogL-null)
    yields the same results.
    Note: Posterior mean of synonymous rate is syn site posterior1 * syn rate1 +
    syn site posterior2 * syn rate2 + syn site posterior3 * syn rate3
    """
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_path}: {e}", file=sys.stderr)
        sys.exit(1)

    for key in ["Evidence Ratios", "Site Log Likelihood"]:
        if key not in data:
            print(f"Error: Missing required key '{key}' in {json_path}", file=sys.stderr)
            sys.exit(1)

    if ["Synonymous site-posteriors"] in data:
        syn_posteriors = data["Synonymous site-posteriors"]

    ev_ratios      = data["Evidence Ratios"]
    site_logl      = data["Site Log Likelihood"]
    

    def _unwrap(arr):
        """Return the flat list, unwrapping a one-element outer list if present."""
        if arr and isinstance(arr[0], list):
            return arr[0]
        return arr

    def _get_optional(mapping, key):
        """Return the unwrapped array for key, or None if key is absent."""
        if key not in mapping:
            return None
        return _unwrap(mapping[key])

    constrained_er      = _get_optional(ev_ratios, "constrained")
    optimized_null_er   = _get_optional(ev_ratios, "optimized null")
    logl_constrained    = _get_optional(site_logl, "constrained")
    logl_optimized_null = _get_optional(site_logl, "optimized null")
    logl_unconstrained  = _get_optional(site_logl, "unconstrained")

    if constrained_er is None:
        print(f"Note: 'constrained' model not found in {os.path.basename(json_path)}; "
              "those columns will be empty.", file=sys.stderr)

    # Determine number of sites from the first available array
    num_sites = next(
        (len(a) for a in [constrained_er, optimized_null_er, logl_unconstrained] if a is not None),
        0,
    )
    num_syn_classes = len(syn_posteriors)

    if num_sites == 0:
        print(f"Error: Could not determine number of sites from {json_path}", file=sys.stderr)
        sys.exit(1)

    headers = [
        "Codon",
        "Constrained Evidence Ratio",
        "Optimized Null Evidence Ratio",
        "Site LogL Constrained",
        "Site LogL Optimized Null",
        "Site LogL Unconstrained"
    ]

    if num_syn_classes:
        for i in range(num_syn_classes):
            headers.append(f"Synonymous Site Posterior {i + 1}")

    def _val(arr, i):
        """Safely retrieve site i from arr; return '' if arr is None."""
        if arr is None:
            return ""
        v = arr[i]
        return v[0] if isinstance(v, list) else v

    try:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(headers)
            for i in range(num_sites):
                row = [
                    i + 1,
                    _val(constrained_er,      i),
                    _val(optimized_null_er,   i),
                    _val(logl_constrained,    i),
                    _val(logl_optimized_null, i),
                    _val(logl_unconstrained,  i),
                ]
                if num_syn_classes:
                    for j in range(num_syn_classes):
                        syn_val = syn_posteriors[j][i]
                        row.append(syn_val[0] if isinstance(syn_val, list) else syn_val)
                writer.writerow(row)
        print(f"Successfully extracted BUSTED data from {os.path.basename(json_path)} to {os.path.basename(csv_path)}")
    except Exception as e:
        print(f"Error writing to {csv_path}: {e}", file=sys.stderr)


# ===========================================================================
# BUSTED – branch / site level
# ===========================================================================

def extract_busted_branches_to_csv(json_path: str, csv_path: str):
    """Per-branch per-site posteriors, EBF and the site-level ER from a BUSTED JSON."""

    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_path}: {e}", file=sys.stderr)
        sys.exit(1)

    ## EBF prior: BUSTED's mixture weights are GLOBAL, so a global prior is
    ## the right prior here.
    rate_dists = None
    try:
        rate_dists = data["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]
    except KeyError as e:
        print(f"Warning: missing rate distributions in {json_path}: {e}. "
              f"EBF columns will be empty.", file=sys.stderr)

    sink_index, sink_note = busted_error_sink_index(data)
    if sink_note:
        print(f"Warning: {sink_note} ({os.path.basename(json_path)})",file=sys.stderr)

    omega_gt_1_indices = []
    global_prior = None
    sink_prior = None
    if isinstance(rate_dists, dict):
        prior = 0.0
        for idx_str, props in rate_dists.items():
            try:
                omega = float(props["omega"])
                proportion = float(props["proportion"])
            except (KeyError, TypeError, ValueError) as e:
                print(f"Warning: unreadable rate class {idx_str} in {json_path}: {e}",
                      file=sys.stderr)
                #continue
            if sink_index is not None and int(idx_str) == sink_index:
                sink_prior = proportion
                sink_omega = omega
                continue
            if omega > 1.0:
                omega_gt_1_indices.append(int(idx_str))
                prior += proportion
        if 0.0 < prior < 1.0:
            global_prior = prior
        else:
            print(f"Note: prior probability of omega>1 is {prior!r} in "
                f"{os.path.basename(json_path)}; EBF columns will be empty.",
                file=sys.stderr)
        if sink_index is not None and sink_prior is not None:
            print(f"Note: error-sink class {sink_index} with value {sink_omega} carried proportion "
                  f"{sink_prior!r} in {os.path.basename(json_path)} and is "
                  f"excluded from the omega>1 prior, which is {prior!r} rather "
                  f"than {prior + sink_prior!r}. Its posterior is reported in "
                  f"its own column.", file=sys.stderr)
        elif sink_index is not None and not sink_prior:
            print(f"Warning: error-sink class {sink_index} with value {sink_omega} was identified in "
                  f"{os.path.basename(json_path)} but no proportion was read "
                  f"for it; the omega>1 prior may include the sink.",
                  file=sys.stderr)

    try:
        branch_attrs = data["branch attributes"]["0"]
        substitutions = data.get("substitutions", {}).get("0", {})
        tree_string = data["input"]["trees"]["0"]
        global_num_sites = data["input"].get("number of sites", 0)
    except KeyError as e:
        print(f"Error: missing branch attributes or tree in {json_path}: {e}",
              file=sys.stderr)
        sys.exit(1)

    ## site-level ER: HyPhy computes it, so prefer its value and only fall
    ## back to recomputing from the two site log-likelihood vectors.
    ## This is already recorded in the previous function, but I'll leave it in for now
    site_er_arr = _unwrap((data.get("Evidence Ratios") or {}).get("optimized null"))
    site_ll = data.get("Site Log Likelihood") or {}
    ll_unconstrained = _unwrap(site_ll.get("unconstrained"))
    ll_constrained = _unwrap(site_ll.get("optimized null")) ## called constraint in the sense that omega =< 1

    def busted_site_er(i):
        v = _at(site_er_arr, i)
        if v is not None:
            return v
        return evidence_ratio(_at(ll_unconstrained, i), _at(ll_constrained, i))

    parent_map = parse_newick_proper(tree_string)
    n_sites_max = max(
        [global_num_sites]
        + [len(a.get("Posterior prob omega class by site", [[]])[0])
           for a in branch_attrs.values()
           if a.get("Posterior prob omega class by site")]
        or [0]
    )
    codons = codon_map_by_site(parent_map, substitutions, n_sites_max)

    headers = [
        "Branch", "Codon Position", "Codon",
        "Posterior Prob (omega>1)", "Prior Prob (omega>1)", "EBF", "log10 EBF",
        "Site ER", "Site ER Model",
        "Posterior Prob (error sink)", "Prior Prob (error sink)",
    ]

    try:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(headers)

            for branch_name, attrs in branch_attrs.items():
                pps_by_site = attrs.get("Posterior prob omega class by site", [])
                if pps_by_site:
                    num_sites = len(pps_by_site[0])
                elif global_num_sites:
                    num_sites = global_num_sites
                else:
                    continue

                for site_idx in range(num_sites):
                    codon_str = codons.get(site_idx, {}).get(branch_name, "")

                    post_prob = None
                    if pps_by_site and global_prior is not None:
                        post_prob = 0.0
                        for idx in omega_gt_1_indices:
                            if idx < len(pps_by_site) and site_idx < len(pps_by_site[idx]):
                                val = pps_by_site[idx][site_idx]
                                post_prob += val[0] if isinstance(val, list) else val

                    sink_post = None
                    if (pps_by_site and sink_index is not None
                            and sink_index < len(pps_by_site)
                            and site_idx < len(pps_by_site[sink_index])):
                        val = pps_by_site[sink_index][site_idx]
                        sink_post = val[0] if isinstance(val, list) else val

                    ebf, log_ebf = empirical_bayes_factor(post_prob, global_prior)
                    writer.writerow([
                        branch_name, site_idx + 1, codon_str,
                        _blank(post_prob), _blank(global_prior),
                        _blank(ebf), _blank(log_ebf),
                        _blank(busted_site_er(site_idx)),
                        "unconstrained vs constrained",
                        _blank(sink_post), _blank(sink_prior)
                    ])

        print(f"Successfully extracted BUSTED branch/site data from "
              f"{os.path.basename(json_path)} to {os.path.basename(csv_path)}")
    except Exception as e:
        print(f"Error writing to {csv_path}: {e}", file=sys.stderr)
        sys.exit(1)


# ===========================================================================
# MEME – branch / site level
# ===========================================================================

def extract_meme_branches_to_csv(json_path: str, csv_path: str):
    """Extract data such as alpha, beta+ and beta- values from a MEME JSON."""

    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_path}: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        branch_attrs = data["branch attributes"]["0"]
        substitutions = data.get("substitutions", {}).get("0", {})
        tree_string = data["input"]["trees"]["0"]
    except KeyError as e:
        print(f"Error: missing branch attributes or tree in {json_path}: {e}",
              file=sys.stderr)
        sys.exit(1)

    col_index, mle_rows = _meme_mle_table(data)
    i_pplus = next((col_index[k] for k in
                    ("p+", "Prop.  omega>1", "Prop. omega>1", "p+ (omega>1), p<sup>+</sup>")
                    if k in col_index), None)
    i_lrt = col_index.get("LRT")
    if i_pplus is None:
        print(f"Warning: no per-site omega>1 mixture weight column in "
              f"{os.path.basename(json_path)} (MLE headers: "
              f"{sorted(col_index)[:8]}). EBF columns will be empty.",
              file=sys.stderr)

    parent_map = parse_newick_proper(tree_string)
    num_sites_global = len(mle_rows)
    codons = codon_map_by_site(parent_map, substitutions, num_sites_global)

    headers = [
        "Branch", "Codon Position", "Codon",
        "Posterior Prob (beta (-))", "Prior Prob (beta(+))", "EBF", "log10 EBF"
    ]

    try:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(headers)

            for branch_name, attrs in branch_attrs.items():
                pps_by_site = attrs.get("Posterior prob omega class by site", [])
                if not pps_by_site or len(pps_by_site) < 2:
                    continue                    # MEME needs the <=1 and >1 classes

                num_sites = len(pps_by_site[1])
                for site_idx in range(num_sites):
                    codon_str = codons.get(site_idx, {}).get(branch_name, "")
                    post_prob = _at(pps_by_site[1], site_idx)   # index 1 = omega>1

                    site_prior = None
                    if i_pplus is not None and site_idx < len(mle_rows):
                        row = mle_rows[site_idx]
                        if i_pplus < len(row):
                            site_prior = row[i_pplus]

                    ## Removed for now
                    # site_er = None
                    # if i_lrt is not None and site_idx < len(mle_rows):
                    #     row = mle_rows[site_idx]
                    #     if i_lrt < len(row):
                    #         site_er = er_from_lrt(row[i_lrt])

                    ## Note: EBF should work but has not been tested since HyPhy vision currently does not display MEME's EBF
                    ebf, log_ebf = empirical_bayes_factor(post_prob, site_prior)
                    writer.writerow([
                        branch_name, site_idx + 1, codon_str,
                        _blank(post_prob), _blank(site_prior),
                        _blank(ebf), _blank(log_ebf)
                    ])

        print(f"Successfully extracted MEME branch/site data from "
              f"{os.path.basename(json_path)} to {os.path.basename(csv_path)}")
    except Exception as e:
        print(f"Error writing to {csv_path}: {e}", file=sys.stderr)
        sys.exit(1)


def _meme_mle_table(data):
    """({header_name: column_index}, rows) from the MLE table."""
    mle = data.get("MLE", {})
    headers = [h[0] if isinstance(h, (list, tuple)) else str(h)
               for h in mle.get("headers", [])]
    content = mle.get("content", {})
    rows = content.get("0", []) if isinstance(content, dict) else content
    return {name: i for i, name in enumerate(headers)}, (rows or [])


# ===========================================================================
# aBSREL – branch / site level + fitted tree
# ===========================================================================

def extract_absrel_branches_to_csv(json_path: str, csv_path: str, tree_path: str,
                                   p_cutoff: float = 0.05):
    """Per-branch per-site posteriors, EBF, site ER and the fitted tree."""

    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_path}: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        branch_attrs = data["branch attributes"]["0"]
        substitutions = data.get("substitutions", {}).get("0", {})
        tree_string = data["input"]["trees"]["0"]
    except KeyError as e:
        print(f"Error: missing keys in {json_path}: {e}", file=sys.stderr)
        sys.exit(1)

    parent_map = parse_newick_proper(tree_string)
    site_er = absrel_site_er(data)

    def safe_p(val):
        try:
            return float(val) if val is not None else 1.0
        except (TypeError, ValueError):
            return 1.0

    selected_branches = [
        br for br, attrs in branch_attrs.items()
        if safe_p(attrs.get("Corrected P-value", 1.0)) <= p_cutoff
    ]
    if selected_branches:
        print(f"Branches under selection (corrected p <= {p_cutoff}): "
              f"{len(selected_branches)}")
    else:
        print(f"No branches under selection in {os.path.basename(json_path)}; "
              f"output files are still written.")

    ## fitted tree, with {Selected} on the significant branches
    new_tree = tree_string
    for branch_name, attrs in branch_attrs.items():
        baseline_len = attrs.get("Baseline MG94xREV")
        if baseline_len is None:
            continue
        node_label = (f"{branch_name}{{Selected}}"
                      if branch_name in selected_branches else branch_name)
        pattern = re.compile(
            rf"\b({re.escape(branch_name)})(?![a-zA-Z0-9_])(:[0-9\.eE+-]+)?")
        new_tree = pattern.sub(f"{node_label}:{baseline_len}", new_tree)

    try:
        with open(tree_path, "w", encoding="utf-8") as f:
            f.write(new_tree.rstrip().rstrip(";") + ";\n")
        print(f"Fitted tree written to {os.path.basename(tree_path)}")
    except Exception as e:
        print(f"Error writing tree to {tree_path}: {e}", file=sys.stderr)
        sys.exit(1)

    ## site count: prefer the posteriors, then the input metadata. Using
    ## len(substitutions) alone produced a header-only tsv whenever HyPhy
    ## was run without substitution reporting.
    num_sites = 0
    for attrs in branch_attrs.values():
        post = attrs.get("posterior") or []
        if post:
            num_sites = max(num_sites, max(len(p) for p in post))
    num_sites = num_sites or data.get("input", {}).get("number of sites", 0) or len(substitutions)
    if not num_sites:
        print(f"Error: cannot determine the number of sites in {json_path}",
              file=sys.stderr)
        sys.exit(1)

    codons = codon_map_by_site(parent_map, substitutions, num_sites)

    headers = [
        "Branch", "Codon Position", "Codon",
        "Posterior Prob (omega>1)", "Prior Prob (omega>1)", "EBF", "log10 EBF",
        "Site ER", "Site ER Model", "Branch LRT",
        "Branch Uncorrected P-value", "Branch Corrected P-value",
    ]

    try:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(headers)

            for branch_name, attrs in branch_attrs.items():
                ## EBF prior: this branch's OWN rate distribution ---
                prior = 0.0
                omega_gt_1_indices = []
                for idx, rate_prop in enumerate(attrs.get("Rate Distributions") or []):
                    if len(rate_prop) == 2:
                        try:
                            rate, prop = float(rate_prop[0]), float(rate_prop[1])
                        except (TypeError, ValueError):
                            continue
                        ## This implemnetation is probably wrong but harmless, as aBSREL will only allow one omega >1, unless error-sink is used
                        if rate > 1.0:
                            prior += prop
                            omega_gt_1_indices.append(idx)
                
                ## A branch with no omega>1 class has no EBF -- not an EBF of 0.
                branch_prior = prior if 0.0 < prior < 1.0 else None


                ## This must be extracted from the Site Log Likelihood at the beginning
                branch_lrt = attrs.get("LRT")
                #branch_er = 
                posteriors = attrs.get("posterior") or []
                er_log_sum, er_seen = 0.0, 0

                for site_idx in range(num_sites):
                    codon_str = codons.get(site_idx, {}).get(branch_name, "")

                    post_prob = None
                    if posteriors and branch_prior is not None:
                        post_prob = 0.0
                        for idx in omega_gt_1_indices:
                            if idx < len(posteriors) and site_idx < len(posteriors[idx]):
                                val = posteriors[idx][site_idx]
                                post_prob += val[0] if isinstance(val, list) else val

                    ebf, log_ebf = empirical_bayes_factor(post_prob, branch_prior)
                    er = site_er(branch_name, site_idx)
                    if er and er > 0:
                        er_log_sum += math.log(er)
                        er_seen += 1

                    writer.writerow([
                        branch_name, site_idx + 1, codon_str,
                        _blank(post_prob), _blank(branch_prior),
                        _blank(ebf), _blank(log_ebf),
                        _blank(er), "full adaptive vs branch null",
                        _blank(branch_lrt),
                        _blank(attrs.get("Uncorrected P-value")),
                        _blank(attrs.get("Corrected P-value")),
                    ])

                ## Probably obsolete
                if er_seen == num_sites and branch_lrt is not None:
                    if not math.isclose(2.0 * er_log_sum, float(branch_lrt),
                                        rel_tol=1e-6, abs_tol=1e-6):
                        print(f"Warning: site ERs for {branch_name} do not "
                              f"reconstruct its LRT (2*sum(log ER)="
                              f"{2 * er_log_sum:.6f}, LRT={branch_lrt}). The "
                              f"resolved null likelihood vector is probably not "
                              f"the right one.", file=sys.stderr)

        print(f"Successfully extracted aBSREL branch/site data to "
              f"{os.path.basename(csv_path)}")
    except Exception as e:
        print(f"Error writing to {csv_path}: {e}", file=sys.stderr)
        sys.exit(1)


# ===========================================================================
# BUSTED – MG94 fitted tree
# ===========================================================================

def extract_busted_mg94_tree(json_path: str, out_path: str):
    """Write the input Newick with branch lengths replaced by MG94xREV values."""

    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_path}: {e}", file=sys.stderr)
        sys.exit(1)

    try:
        newick      = data["input"]["trees"]["0"]
        branch_attrs = data["branch attributes"]["0"]
    except KeyError as e:
        print(f"Error: Missing expected key in {json_path}: {e}", file=sys.stderr)
        sys.exit(1)

    MG94_KEY = "MG94xREV with separate rates for branch sets"

    # Build name -> MG94 length mapping
    mg94_lengths = {
        name: attrs[MG94_KEY]
        for name, attrs in branch_attrs.items()
        if MG94_KEY in attrs
    }

    if not mg94_lengths:
        print(f"Warning: No MG94 branch lengths found in {json_path}; writing input tree unchanged.",
              file=sys.stderr)
        mg94_tree = newick
    else:
        # Regex: match node_name:number and substitute the length when name is known
        def _replacer(m):
            name = m.group(1)
            if name in mg94_lengths:
                return f"{name}:{mg94_lengths[name]:.10g}"
            return m.group(0)  # root / unknown nodes: keep original

        # Node names in HyPhy Newick: letters, digits, underscores, dots, hyphens
        mg94_tree = re.sub(r'([\w.\-]+):([\d.eE+\-]+)', _replacer, newick)

    try:
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(mg94_tree.rstrip().rstrip(";") + ";\n")
        print(f"MG94 tree written to {os.path.basename(out_path)}")
        n_replaced = sum(1 for name in mg94_lengths if name in newick)
        print(f"  Branch lengths replaced: {n_replaced} / {len(mg94_lengths)} known branches")
    except Exception as e:
        print(f"Error writing tree to {out_path}: {e}", file=sys.stderr)


# ===========================================================================
# MLE – generic MLE table extraction
# ===========================================================================

def extract_mle_to_csv(json_path: str, csv_path: str):
    """Extract the MLE content table from a HyPhy JSON.
    The different ML estimators are all dumped in a table-like object in the beginning 
    of the json file, from where they have to be extracted.
    Data looks like {MLE : {content: {0}}, {headers}}
    """
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_path}: {e}", file=sys.stderr)
        sys.exit(1)

    if "MLE" not in data:
        print(f"Error: 'MLE' key not found in {json_path}", file=sys.stderr)
        sys.exit(1)

    mle_data = data["MLE"]
    
    headers = [h[0] if isinstance(h, (list, tuple)) else str(h)
               for h in mle_data.get("headers", [])]
    content = mle_data.get("content", {})
    if isinstance(content, dict):
        if "0" not in content:
            print(f"Error: MLE content in {json_path} has no '0' partition "
                  f"(found {sorted(content)[:5]}).", file=sys.stderr)
            sys.exit(1)
        rows = content["0"]
    else:
        rows = content

    try:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter='\t')
            if headers:
                writer.writerow(headers)
            writer.writerows(rows)
        print(f"Successfully extracted MLE data from {os.path.basename(json_path)} to {os.path.basename(csv_path)}")
    except Exception as e:
        print(f"Error writing to {csv_path}: {e}", file=sys.stderr)

# ===========================================================================
# RELAX – branch level
# ===========================================================================

def extract_relax_branches_to_csv(json_paths: list, csv_path: str):
    """Extract per-branch k and p-values from one or more RELAX JSONs."""
    all_rows = []
    
    for json_path in json_paths:
        try:
            with open(json_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception as e:
            print(f"Error reading {json_path}: {e}", file=sys.stderr)
            continue
            
        try:
            branch_attrs = data.get("branch attributes", {}).get("0", {})
            tested = data.get("tested", {}).get("0", {})
            test_results = data.get("test results", {})
            p_val = test_results.get("p-value", "")
        except AttributeError as e:
            print(f"Error: Unexpected structure in {json_path}: {e}", file=sys.stderr)
            continue

        for branch_name, attrs in branch_attrs.items():
            k_val = attrs.get("k (general descriptive)", "")
            group = tested.get(branch_name, "")
            
            # If it's the foreground ('Test') branch, we capture its specific test p-value and its k
            branch_p_val = p_val if group == "Test" else ""
            all_rows.append([branch_name, k_val, branch_p_val])

    if not all_rows:
        print("Error: No branch data extracted from any RELAX JSONs.", file=sys.stderr)
        sys.exit(1)

    headers = ["Branch", "k", "p-value"]
    try:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(headers)
            writer.writerows(all_rows)
        print(f"Successfully extracted RELAX branch data to {os.path.basename(csv_path)}")
    except Exception as e:
        print(f"Error writing to {csv_path}: {e}", file=sys.stderr)

# ===========================================================================
# CLI
# ===========================================================================

def _find_json(pattern: str, label: str) -> str:
    """Glob for JSON files matching *pattern* and return the first hit."""
    matches = sorted(glob.glob(f"*{pattern}*.json"))
    if not matches:
        print(f"No {label} JSON files found in the current directory. Please provide one as an argument.", file=sys.stderr)
        sys.exit(1)
    json_file = matches[0]
    print(f"No file specified. Using: {json_file}")
    return json_file


def _stem(json_path: str) -> str:
    """Return the path without the .json extension."""
    return json_path[:-5] if json_path.endswith(".json") else json_path


def main():
    parser = argparse.ArgumentParser(
        prog="extract_hyphy.py",
        description="Extract TSV data from HyPhy JSON result files.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- busted ---
    p_busted = subparsers.add_parser(
        "busted",
        help="BUSTED: extract site-level stats, per-branch/site posteriors, and MG94 tree",
    )
    p_busted.add_argument("json", nargs="?", help="Path to BUSTED JSON (auto-detected if omitted)")
    p_busted.add_argument("-o", "--output",          help="Site-level TSV path (default: <stem>_extracted.tsv)")
    p_busted.add_argument("-b", "--branches-output", help="Branch/site TSV path (default: <stem>_branches_per_site.tsv)")
    p_busted.add_argument("-t", "--tree",             help="MG94 tree path (default: <stem>_mg94.nwk)")

    # --- meme ---
    p_meme = subparsers.add_parser(
        "meme",
        help="MEME: extract per-branch/site posteriors AND MLE table (two TSVs)",
    )
    p_meme.add_argument("json", nargs="?", help="Path to MEME JSON (auto-detected if omitted)")
    p_meme.add_argument("-o", "--output",      help="Branch/site TSV path (default: <stem>_branches_per_site.tsv)")
    p_meme.add_argument("-m", "--mle-output",  help="MLE TSV path (default: <stem>_mle.tsv)")

    # --- absrel ---
    p_absrel = subparsers.add_parser(
        "absrel",
        help="aBSREL: extract per-branch/site posteriors and fitted tree",
    )
    p_absrel.add_argument("json", nargs="?", help="Path to aBSREL JSON (auto-detected if omitted)")
    p_absrel.add_argument("-o", "--output", help="Branch/site TSV path (default: <stem>_selected_branches.tsv)")
    p_absrel.add_argument("-t", "--tree",   help="Fitted-tree path (default: <stem>_fitted_tree.nwk)")
    p_absrel.add_argument("-p", "--p-cutoff", type=float, default=0.05,
                          help="Corrected p-value below which a branch counts as "
                               "selected, for the {Selected} tree labels. "
                               "(default: %(default)s)")

    # --- relax ---
    p_relax = subparsers.add_parser(
        "relax",
        help="RELAX: extract branch-level k and foreground p-values from multiple JSONs",
    )
    p_relax.add_argument("jsons", nargs="+", help="Paths to one or more RELAX JSONs (supports wildcards)")
    p_relax.add_argument("-o", "--output", help="Output TSV path (default: <stem of first json>_relax.tsv)")

    args = parser.parse_args()

    # ---- dispatch ----
    if args.command == "busted":
        json_file     = args.json or _find_json("busted", "BUSTED")
        site_csv      = args.output          or _stem(json_file) + "_extracted.tsv"
        branches_csv  = args.branches_output or _stem(json_file) + "_branches_per_site.tsv"
        tree_out      = args.tree            or _stem(json_file) + "_mg94.nwk"
        extract_busted_to_csv(json_file, site_csv)
        extract_busted_branches_to_csv(json_file, branches_csv)
        extract_busted_mg94_tree(json_file, tree_out)

    elif args.command == "meme":
        json_file    = args.json or _find_json("meme", "MEME")
        branches_csv = args.output     or _stem(json_file) + "_branches_per_site.tsv"
        mle_csv      = args.mle_output or _stem(json_file) + "_mle.tsv"
        extract_meme_branches_to_csv(json_file, branches_csv)
        extract_mle_to_csv(json_file, mle_csv)

    elif args.command == "absrel":
        json_file = args.json or _find_json("absrel", "aBSREL")
        csv_file  = args.output or _stem(json_file) + "_selected_branches.tsv"
        tree_file = args.tree   or _stem(json_file) + "_fitted_tree.nwk"
        extract_absrel_branches_to_csv(json_file, csv_file, tree_file, args.p_cutoff)

    elif args.command == "relax":
        expanded_jsons = []
        for j_path in args.jsons:
            if "*" in j_path or "?" in j_path:
                expanded_jsons.extend(glob.glob(j_path))
            else:
                expanded_jsons.append(j_path)
        
        if not expanded_jsons:
            print("Error: No RELAX JSON files found matching the provided arguments.", file=sys.stderr)
            sys.exit(1)
            
        csv_file = args.output or _stem(expanded_jsons[0]) + "_relax.tsv"
        extract_relax_branches_to_csv(expanded_jsons, csv_file)


if __name__ == "__main__":
    main()
