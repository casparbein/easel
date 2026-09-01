## Align with Muscle 5
rule extract_ali_muscle_afa:
    input:
        in_fasta =  f"{config['fastaPath']}{{transcript_id}}{config['fileSuffix']}"
    output:
        out_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.afa",
    threads:
        10
    resources:
        runtime = "2h",
        mem_mb = 20000
    group: "align_clean"
    log:
        "logs/muscle_afa/{transcript_id}.log"
    conda:
        "../envs/muscle.yaml"
    shell:
        """
        muscle \
        -align {input.in_fasta} \
        -diversified \
        -output {output.out_fasta} \
        >> {log} 2>&1
        """

## Get best alignment (MaxCC)
rule extract_ali_muscle_best:
    input:
        in_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.afa",
    output:
        out_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa",
    group: "align_clean"
    resources:
        runtime = "10m",
    log:
        "logs/muscle_best/{transcript_id}.log"
    conda:
        "../envs/muscle.yaml"
    shell:
        """
        muscle \
        -maxcc {input.in_fasta} \
        -output {output.out_fasta} \
        >> {log} 2>&1
        """

## Get column confidence scores
rule extract_ali_muscle_cc:
    input:
        in_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa" 
    output:
        out_conf = "codon_alignments/{transcript_id}/tmp/{transcript_id}_conf.fa"
    params:
        in_ensemble = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.afa"
    group: "align_clean"
    resources:
        runtime = "10m",
    log:
        "logs/muscle_cc/{transcript_id}.log"
    conda:
        "../envs/muscle.yaml"
    shell:
        """
        muscle \
        -letterconf {params.in_ensemble} \
        -ref {input.in_fasta} \
        -output {output.out_conf} \
        >> {log} 2>&1
        """
