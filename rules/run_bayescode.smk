## Run bayescode, as per https://www.pnas.org/doi/10.1073/pnas.2214977120#sec-1

def get_input_tree_bayescode(wildcards):
    if config["settings"]["selectionSettings"]["foregroundLst"]:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}_tmp.treefile"
    elif config["settings"]["treeSettings"]["preCompGeneTrees"]["geneTreePath"]:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}_pruned_tree_bayescode.nh"
    else:
        return "codon_alignments/{transcript_id}/{transcript_id}_iqtree.nh"

rule convert_to_phylip:
    input:
        get_input_for_busted
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}.phy"
    group: "bayescode_stack"
    log:
        "logs/convert_to_phylip/{transcript_id}.log"
    threads: 
        1
    resources: 
        runtime = "10m",
        mem_mb = 5000
    conda:
        "../envs/bayescode.yaml"
    shell:
        """
        fasta_to_ali.py \
        -i {input} \
        -o {output} \
        >> {log} 2>&1
        """

rule mutselomega_w0:
    input:
        ali= "codon_alignments/{transcript_id}/tmp/{transcript_id}.phy",
        tree = get_input_tree_bayescode
    output:
        temp("codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.trace"),
        temp("codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.chain"),
        temp("codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.param"),
        temp("codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.run"),
    params:
        out_name = "codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree",
    group: "bayescode_stack"
    log:
        "logs/mutselomega_w0/{transcript_id}.log"
    threads:
        1
    resources:
        runtime = "8h",
        mem_mb = 5000
    conda:
        "../envs/bayescode.yaml"
    shell:
        """
        mutselomega \
        --ncat 30 \
        -a {input.ali} \
        -t {input.tree} \
        --until 2000 \
        {params.out_name} \
        >> {log} 2>&1
        """

rule mutselomega_w:
    input:
        ali= "codon_alignments/{transcript_id}/tmp/{transcript_id}.phy",
        tree = get_input_tree_bayescode
    output:
        temp("codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.trace"),
        temp("codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.chain"),
        temp("codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.param"),
        temp("codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.run"),
    params:
        out_name = "codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree",
    group: "bayescode_stack"
    log:
        "logs/mutselomega_w/{transcript_id}.log"
    threads:
        1
    resources:
        runtime = "8h",
        mem_mb = 5000
    conda:
        "../envs/bayescode.yaml"
    shell:
        """
        mutselomega \
        --freeomega \
        --omegancat 30 \
        --flatfitness \
        -a {input.ali} \
        -t {input.tree} \
        --until 2000 \
        {params.out_name} \
        >> {log} 2>&1
        """

rule read_mutselomega_w:
    input:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.trace",
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.chain",
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.param",
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.run",
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.ci0.025.tsv"
    params:
        in_name = "codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree",
    group: "bayescode_stack"
    log:
        "logs/read_mutselomega_w/{transcript_id}.log"
    threads:
        1
    resources:
        runtime = "10m",
        mem_mb = 5000
    conda:
        "../envs/bayescode.yaml"
    shell:
        """
        readmutselomega \
        --every 1 \
        --until 2000 \
        --burnin 1000 \
        --confidence_interval 0.025 \
        --omega \
        {params.in_name} \
        >> {log} 2>&1
        """

rule read_mutselomega_w0:
    input:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.trace",
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.chain",
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.param",
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.run",
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.ci0.025.tsv"
    params:
        in_name = "codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree",
    group: "bayescode_stack"
    log:
        "logs/read_mutselomega_w0/{transcript_id}.log"
    threads:
        1
    resources:
        runtime = "10m",
        mem_mb = 5000
    conda:
        "../envs/bayescode.yaml"
    shell:
        """
        readmutselomega \
        --every 1 \
        --until 2000 \
        --burnin 1000 \
        --confidence_interval 0.025 \
        --omega_0 \
        {params.in_name} \
        >> {log} 2>&1
        """