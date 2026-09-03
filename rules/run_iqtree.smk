## Gene tree reconstruction with IQ tree2.
rule compute_tree:
    input:
        ali = get_input_for_tree
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_tmp.treefile"
    params:
        bootstrap = 1000,
        prefix = "codon_alignments/{transcript_id}/tmp/{transcript_id}_tmp",
        #reference = config["referenceName"] ## Will be outgroup taxon in IQTree. Does this make sense??, in command, it was -o {params.reference}
    threads:
        config["resources"]["tree"]["threads"]
    resources:
        runtime = "8h",
        mem_mb = config["resources"]["tree"]["mem_mb"]
    group: "tree_gene"
    log:
        "logs/compute_tree/{transcript_id}.log"
    conda:
        "../envs/iqtree.yaml"
    shell:
        """
        iqtree3 \
        -s {input.ali} \
        --seqtype CODON \
        --prefix {params.prefix} \
        -B {params.bootstrap} \
        --mem {resources.mem_mb}M \
        -T {threads} \
        -blmin 0.001 \
        -mset MG,MGK \
        -nstop 50 \
        >> {log} 2>&1
        """

## Tree annotation for downstream HyPhy runs
## Adapt to directly run script from snakemake
## ISSUE: logging does not work properly yet
rule gene_tree_annotation:
    input:
        in_tree_boot = "codon_alignments/{transcript_id}/tmp/{transcript_id}_tmp.treefile"
    output:
        ans_tree = "codon_alignments/{transcript_id}/{transcript_id}_iqtree.nh"
    params:
        label_nodes = config["settings"]["selectionSettings"]["foregroundLst"] if config["settings"]["selectionSettings"]["foregroundLst"] else None,
    group: "tree_gene"
    resources:
        runtime = "10m",
    log:
        "logs/gene_tree_annotation/{transcript_id}.log"
    conda:
        "../envs/newick_tree_manipulator.yaml"
    script:
        "../scripts/newick_tree_manipulator.py"


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

rule draw_tree_bayescode_ct:
    input: 
        in_tree = "codon_alignments/{transcript_id}/tmp/{transcript_id}_tmp.treefile",
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


