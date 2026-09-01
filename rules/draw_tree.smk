## The tree including all species that one wants to analyse has to be created at the beginning as one of the input files.
## Here, the tree is pruned so only leaves that have representation in the alignment are left
def extract_names_from_fasta(fasta):
    tmp_list = []
    with open(fasta, 'r') as fas:
        for line in fas:
            if line.startswith('>'):
                tmp_list.append(line.strip('>').strip('\n'))
    if len(tmp_list) > 0:
        return(",".join(tmp_list))
    else:
        return ""
               
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
        keep = lambda wildcards, input: extract_names_from_fasta(input.ali)
    resources:
        runtime = "10m",
        mem_mb = config["resources"]["extractAlignments"]["mem_mb"]
    threads: config["resources"]["extractAlignments"]["threads"]
    group: "tree_prune"
    log:
        "logs/draw_tree/{transcript_id}.log"
    conda:
        "../envs/newick_tree_manipulator.yaml"
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
        keep = lambda wildcards, input: extract_names_from_fasta(input.ali)
    resources:
        runtime = "10m",
        mem_mb = config["resources"]["extractAlignments"]["mem_mb"]
    threads: config["resources"]["extractAlignments"]["threads"]
    group: "tree_prune"
    log:
        "logs/draw_tree_bayescode/{transcript_id}.log"
    conda:
        "../envs/newick_tree_manipulator.yaml"
    script:
        "../scripts/newick_tree_manipulator.py"