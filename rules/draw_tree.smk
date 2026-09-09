## Prune a tree, keeping all species that are present in input alignment and removing all else (+reference removal, which is probably obsolete and will be discontinued)
## Adapt to directly run script from snakemake
rule draw_tree:
    input: 
        in_tree = config['settings']['treeSettings']['inputSpeciesTree']["treeFile"],
        ali=get_input_for_busted,
    output:
        ans_tree = "codon_alignments/{transcript_id}/{transcript_id}_pruned_tree.nh"
    params:
        label_nodes = config["settings"]["selectionSettings"]["foregroundLst"] if config["settings"]["selectionSettings"]["foregroundLst"] else None, #foreground_string, "None",
    resources:
        runtime = "5m",
        mem_mb = 500
    threads: config["resources"]["extractAlignments"]["threads"]
    localrule: True,
    log:
        "logs/draw_tree/{transcript_id}.log"
    # conda:
    #     "../envs/newick_tree_manipulator.yaml"
    script:
        "../scripts/newick_tree_manipulator.py"

rule draw_tree_bayescode:
    input: 
        in_tree = config['settings']['treeSettings']['inputSpeciesTree']["treeFile"],
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
        "logs/draw_tree_bayescode/{transcript_id}.log"
    conda:
        "../envs/newick_tree_manipulator.yaml"
    script:
        "../scripts/newick_tree_manipulator.py"