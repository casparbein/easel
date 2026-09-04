#!/usr/bin/env python3

from ete4 import Tree
import argparse
import sys
import logging

DESCRIPTION = '''\
Helper script for tree manipulation of phylogenetic trees in newick format.
Trees can be pruned (leaves are then removed), ancestral branches labelled
(outermost children leaves combined by an underscore, for example: A_B),
and Fore/Background labels can be added to nodes in the tree,
optionally including ancestral nodes as well (if all descendants are also labelled as Fore/Background).
Additionally, if ancestors are named, bootstrap support values are removed.
'''

__author__ = "Bernhard Bein, 2025."

## Logging
log = logging.getLogger(__name__)

## CLI arguments (not used inside snakemake)
def tree_parser():
    """Parse CMD args."""
    app = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(
            prog, max_help_position=6, indent_increment=2
        ),
    )

    app.add_argument(
        "-t",
        "--input_tree",
        action="store",
        dest="in_tree",
        type=str,
        help="path to input tree. Must be in newick format.\n      ",
    )

    app.add_argument(
        "-p",
        "--prune",
        action="store",
        dest="to_prune",
        default=None,
        type=str,
        help=(
            "comma-separated string of branches to prune off the tree.\n"
            "Must be tips, as internal branches are automatically removed.\n"
            "Example: -p tip1,tip2,tip3\n      "
        ),
    )

    app.add_argument(
        "-k",
        "--keep",
        action="store",
        dest="to_keep",
        type=str,
        help=(
            "comma-separated string of branches to keep in tree.\n"
            "Must be tips. Will be employed after pruning,\n"
            "everything that was pruned cannot be rescued.\n"
            "Example: -k tip1,tip2,tip3\n      "
        ),
    )

    app.add_argument(
        "-l",
        "--label",
        action="append",
        nargs="+",
        dest="to_label",
        help=(
            "comma-separated string of branches that should be labelled as Foreground (or anything other).\n"
            "This is necessary for a RELAX screen (and for assigning different omega rate classes in\n"
            "selection simulations). If -a/--ancestor is enabled, ancestral branches\n"
            "which leaves are all Foreground will also be labelled as Foreground.\n"
            "Example: -l tip1,tip2,tip3\n      "
        ),
    )

    app.add_argument(
        "-ls",
        "--label_string",
        action="append",
        dest="label_string",
        help=(
            "What the list of labelled species should be labelled,\n"
            "for instance 'Foreground'.\n"
            "Default: Foreground\n      "
        ),
    )

    app.add_argument(
        "-la",
        "--label_ancestor",
        action="store_true",
        help=(
            "If invoked, will add Foreground labels to all internal branches\n"
            "if all their leaves are listed in -l\n      "
        ),
    )

    app.add_argument(
        "-a",
        "--ancestor",
        action="store_true",
        help=(
            "If invoked, will label all internal branches\n"
            "by leaves from their left/right subtrees.\n"
            "Example naming: tip1_tip2.\n      "
        ),
    )

    app.add_argument(
        "-z",
        "--no_branchlen",
        action="store_true",
        help=(
            "If invoked, tree will be output as a dendrogram,\n"
            "with branchlengths removed.\n      "
        ),
    )

    app.add_argument(
        "-o",
        "--out_tree",
        action="store",
        dest="out_tree",
        type=str,
        help="Path to output tree.\nExample: -o output_tree.nh\n      ",
    )

    args = app.parse_args()
    return args


## Read newick tree
def read_in_tree(tree_path):
    in_tree = Tree(open(tree_path), parser=1)
    return in_tree


## Read labels from a comma-separated string
def read_labels(label_input):
    labels = label_input.split(",")
    return labels


## Read labels from a list file (not yet implemented)
def read_labels_from_list(label_input_path):
    label_list = []
    with open(label_input_path, "r") as labels:
        for line in labels:
            label = line.strip().split("\n")
            label_list.append(label)
    return label_list


## Return tree that is pruned by tips
def prune_tips(tree, not_prune_labels):
    tree.prune(not_prune_labels, preserve_branch_length=True)
    return tree


## Name internal ancestors in tree
def name_ancestors(tree):
    for node in tree.traverse("levelorder"):
        if not node.is_leaf:
            left_leaf = list(node.children[0].leaves())[0].name
            right_leaf = list(node.children[-1].leaves())[0].name
            node.name = f"{left_leaf}_{right_leaf}"
    return tree

## Name Foreground branches in tree
def label_nodes(tree, tip_labels, label_string, ancestor=True):
    found = False
    label = label_string
    for leaf in tree.leaves():
        if leaf.name in tip_labels:
            found = True
            leaf.name = f"{leaf.name}{{{label}}}"
    if not found:
        log.warning("No leaf with the name/s {} was found".format(tip_labels))

    ## Build the set of labelled leaf names for ancestor detection
    tip_labels_foreground = {f"{tip}{{{label}}}" for tip in tip_labels}

    ## Label ancestors if all children are included in tip_labels
    if ancestor:
        for node in tree.traverse("postorder"):
            if node.is_leaf:
                continue
            descendant_labels = {
                leaf.name
                for leaf in node.leaves()
                if leaf.name in tip_labels_foreground
            }
            if len(descendant_labels) == len(list(node.leaves())):
                node.name = f"{node.name}{{{label}}}"
    return tree


def main():
    #args = tree_parser()

    input_tree = snakemake.input[0]
    output_tree = snakemake.output[0]
    branches_to_prune = None
    branches_to_keep = getattr(snakemake.params, "keep", False)
    branches_to_label = getattr(snakemake.params, "label_nodes", False)
    branches_to_label_nested = [branches_to_label]
    ancestor = True
    zerolen = False
    label_ancestor = True
    label_string = ["Foreground"]
    logfile = snakemake.log[0]

    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="[%(levelname)s] %(message)s",
        force=True
    )

    log.info(f"Branches to keep: {branches_to_keep}")
    log.info(f"Branches to label: {branches_to_label}")

    ## test that labels do not overlap if there are several lists passed
    if branches_to_label:
        ## Each entry in branches_to_label is a list with one comma-separated string. This is for when more than Foreground should be labelled
        all_tips = [
            tip
            for group in branches_to_label_nested
            for entry in group
            for tip in entry.split(",")
        ]
        if len(set(all_tips)) != len(all_tips):
            log.critical("Clade labels overlap: each tip must belong to only one clade.")
            sys.exit(1)

    ## set Foreground if label is not defined
    if branches_to_label and not label_string:
        label_string = ["Foreground"]

    ## read newick tree
    in_tree = read_in_tree(input_tree)

    ## name ancestors
    if ancestor:
        log.info("All internal branches will be labelled.")
        in_tree = name_ancestors(in_tree)

    ## prune tree based on labels
    if branches_to_prune is not None:
        log.info(
            "The following branches will be pruned from the tree: {}".format(
                branches_to_prune
            )
        )

        ## create prune and label lists
        branches_to_prune_list = read_labels(branches_to_prune)

        ## Which branches to keep
        branches_not_to_prune = []
        found = False
        for leaf in in_tree.leaves():
            if leaf.name in branches_to_prune_list:
                found = True
                continue
            else:
                branches_not_to_prune.append(leaf.name)

        ## Prune tree
        in_tree = prune_tips(in_tree, branches_not_to_prune)

        ## Notify user if nothing was pruned
        if not found:
            log.info(
                "No tip was pruned as none of {} was found in the tree".format(
                    branches_to_prune_list
                )
            )

    if branches_to_keep:
        log.info(
            "The following branches will be kept in tree: {}".format(branches_to_keep)
        )

        ## Create list
        branches_to_keep_list = read_labels(branches_to_keep)

        branches_to_keep_clean = [
            leaf.name
            for leaf in in_tree.leaves()
            if leaf.name in branches_to_keep_list
        ]

        ## prune
        in_tree = prune_tips(in_tree, branches_to_keep_clean)

    ## label branches and leaves
    ## Currently only implemented for "Foreground", Middle- and Background will be implemented at some point
    if branches_to_label:
        ## Works also for single lists, as they will be wrapped in another list
        original_groups = list(branches_to_label_nested)
        for i, subclade in enumerate(original_groups):
            subclade_string = ",".join(subclade)
            ## Fall back to "Foreground" if fewer label strings than label groups
            lstr = label_string[i] if i < len(label_string) else "Foreground"
            log.info(
                "The following tips will be labelled as '{}': {}".format(
                    lstr, subclade_string
                )
            )
            log.info(
                "Dynamic ancestor labelling is turned on: {}".format(
                    label_ancestor
                )
            )

            tip_list = read_labels(subclade_string)
            in_tree = label_nodes(in_tree, tip_list, lstr, ancestor=label_ancestor)

    if zerolen:
        log.info("Zero branchlength invoked: All branchlengths will be set to 0")
        for node in in_tree.traverse():
            node.dist = 0.0

    in_tree.write(parser=1, outfile=output_tree)


if __name__ == "__main__":
    main()
