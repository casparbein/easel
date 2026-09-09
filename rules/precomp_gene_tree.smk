import os

## Directory of precomputed gene trees, hoisted so the rule bodies below stay
## readable and parse on Python < 3.12 (nested same-quote f-strings are 3.12+).
GENE_TREE_DIR = config["settings"]["treeSettings"]["preCompGeneTrees"]["geneTreePath"]


## Copy gene trees
rule copy_precomp_gene_tree:
    input:
        in_tree = os.path.join(GENE_TREE_DIR, "{transcript_id}" + config["treeSuffix"]),
        ali=get_input_for_busted,
    output:
        out_tree = "codon_alignments/{transcript_id}/{transcript_id}_pruned_tree.nh"
    params:
        label_nodes = config["settings"]["selectionSettings"]["foregroundLst"] if config["settings"]["selectionSettings"]["foregroundLst"] else None,
    localrule: True,
    resources:
        runtime = "5m",
        mem_mb = 500
    log:
        "logs/copy_precomp_gene_tree/{transcript_id}.log"
    #conda:
    #    "../envs/newick_tree_manipulator.yaml" ## This is probably contained in snakemake environment
    script:
        "../scripts/newick_tree_manipulator.py"

## For Bayescode
rule precomp_tree_bayescode:
    input: 
        in_tree = os.path.join(GENE_TREE_DIR, "{transcript_id}" + config["treeSuffix"]),
        ali=get_input_for_busted,
    output:
        ans_tree = "codon_alignments/{transcript_id}/tmp/{transcript_id}_pruned_tree_bayescode.nh"
    params:
        label_nodes = None,
    resources:
        runtime = "5m",
        mem_mb = 500
    threads: config["resources"]["extractAlignments"]["threads"]
    localrule: True,
    log:
        "logs/precomp_tree_bayescode/{transcript_id}.log"
    #conda:
    #    "../envs/newick_tree_manipulator.yaml"
    script:
        "../scripts/newick_tree_manipulator.py"