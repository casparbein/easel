## Align with Muscle 5
rule extract_ali_muscle_afa:
    input:
        in_fasta =  f"{config["fastaPath"]}{{transcript_id}}{config["fileSuffix"]}"
    output:
        out_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.afa",
    threads:
        10
    resources:
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
        -output {output.out_fasta}
        """

## Get best alignment (MaxCC)
rule extract_ali_muscle_best:
    input:
        in_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.afa",
    output:
        out_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa",
    group: "align_clean"
    log:
        "logs/muscle_best/{transcript_id}.log"
    conda:
        "../envs/muscle.yaml"
    shell:
        """
        muscle \
        -maxcc {input.in_fasta} \
        -output {output.out_fasta}
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
    log:
        "logs/muscle_cc/{transcript_id}.log"
    conda:
        "../envs/muscle.yaml"
    shell:
        """
        muscle \
        -letterconf {params.in_ensemble} \
        -ref {input.in_fasta} \
        -output {output.out_conf}
        """

rule filter_muscle_confidence:
    input:
        in_ali = "codon_alignments/{transcript_id}/{transcript_id}_ori.fa",
        conf_scores = "codon_alignments/{transcript_id}/{transcript_id}_conf.fa"
    output:
        filtered_ali = "codon_alignments/{transcript_id}/{transcript_id}_muscle_filtered.fa"
    group: "align_clean"
    params:
        reference = config.get("referenceName") or "",
    log: "logs/filter_muscle_confidence/{transcript_id}.log"
    script: 
        "../scripts/filter_muscle_conf.py"