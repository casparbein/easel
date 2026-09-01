## Directory of precomputed gene trees, hoisted so the rule bodies below stay
## readable and parse on Python < 3.12 (nested same-quote f-strings are 3.12+).
GENE_TREE_DIR = config["settings"]["treeSettings"]["preCompGeneTrees"]["geneTreePath"]

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

## Copy gene trees
rule copy_precomp_gene_tree:
    input:
        in_tree = f"{GENE_TREE_DIR}{{transcript_id}}{config['treeSuffix']}",
        ali=get_input_for_busted,
    output:
        out_tree = "codon_alignments/{transcript_id}/{transcript_id}_pruned_tree.nh"
    params:
        label_nodes = config["settings"]["selectionSettings"]["foregroundLst"] if config["settings"]["selectionSettings"]["foregroundLst"] else None,
        keep = lambda wildcards, input: extract_names_from_fasta(input.ali),
    group: "tree_prep"
    resources:
        runtime = "10m",
    log:
        "logs/copy_precomp_gene_tree/{transcript_id}.log"
    conda:
        "../envs/newick_tree_manipulator.yaml" ## This is probably contained in snakemake environment
    script:
        "../scripts/newick_tree_manipulator.py"

## For Bayescode
rule precomp_tree_bayescode:
    input: 
        in_tree = f"{GENE_TREE_DIR}{{transcript_id}}{config['treeSuffix']}",
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
    group: "tree_prep"
    log:
        "logs/precomp_tree_bayescode/{transcript_id}.log"
    conda:
        "../envs/newick_tree_manipulator.yaml"
    script:
        "../scripts/newick_tree_manipulator.py"