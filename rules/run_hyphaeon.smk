rule hyphaeon_meme:
    input:
        ali  = get_input_ali_for_hyphy,
        tree = get_input_tree_for_hyphy
    output:
        json = "codon_alignments/{transcript_id}/hyphaeon/{transcript_id}.meme.json",
        csv = "codon_alignments/{transcript_id}/hyphaeon/{transcript_id}.meme.csv"
    localrule: True
    resources:
        runtime = "5min",
        mem_mb  = 1000
    log: "logs/hyphaeon_meme/{transcript_id}.log"
    conda: "../envs/hyphaeon.yaml"
    shell:
        """
        hyphaeon meme \
        -a {input.ali} \
        -t {input.tree} \
        -o {output.json} \
        -c {output.csv}
        >> {log} 2>&1
        """

rule hyphaeon_attribute:
    input:
        ali  = get_input_ali_for_hyphy,
        tree = get_input_tree_for_hyphy
    output:
        json = "codon_alignments/{transcript_id}/hyphaeon/{transcript_id}.attr.json",
    localrule: True
    resources:
        runtime = "5min",
        mem_mb  = 1000
    log: "logs/hyphaeon_attribute/{transcript_id}.log"
    conda: "../envs/hyphaeon.yaml"
    shell:
        """
        hyphaeon meme \
        -a {input.ali} \
        -t {input.tree} \
        --attribute \
        -o {output.json} 
        >> {log} 2>&1
        """

rule hyphaeon_dms:
    input:
        ali  = get_input_ali_for_hyphy,
        tree = get_input_tree_for_hyphy
    output:
        json = "codon_alignments/{transcript_id}/hyphaeon/{transcript_id}.dms.json",
        csv = "codon_alignments/{transcript_id}/hyphaeon/{transcript_id}.dms.csv"
    localrule: True
    resources:
        runtime = "5min",
        mem_mb  = 1000
    log: "logs/hyphaeon_dms/{transcript_id}.log"
    conda: "../envs/hyphaeon.yaml"
    shell:
        """
        hyphaeon dms \
        -a {input.ali} \
        -t {input.tree} \
        -o {output.json} \
        -c {output.csv}
        >> {log} 2>&1
        """

rule hyphaeon_epistasis:
    input:
        ali  = get_input_ali_for_hyphy,
        tree = get_input_tree_for_hyphy
    output:
        json = "codon_alignments/{transcript_id}/hyphaeon/{transcript_id}.ep.json",
        csv = "codon_alignments/{transcript_id}/hyphaeon/{transcript_id}.ep.csv",
        graph = "codon_alignments/{transcript_id}/hyphaeon/{transcript_id}.ep.graphml",
    localrule: True
    resources:
        runtime = "5min",
        mem_mb  = 1000
    log: "logs/hyphaeon_epistasis/{transcript_id}.log"
    conda: "../envs/hyphaeon.yaml"
    shell:
        """
        hyphaeon epistasis \
        -a {input.ali} \
        -t {input.tree} \
        --min-sim 0.30 \
        --n-permutations 10000 \
        --max-perm-p 0.05 \
        -o {output.json} \
        -c {output.csv} \
        --graphml {output.graph}
        >> {log} 2>&1
        """

