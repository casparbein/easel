## Align with MACSE2
rule extract_ali_macse2:
    input:
        in_fasta = f"{config['fastaPath']}{{transcript_id}}{config['fileSuffix']}"
    output:
        out_fasta_nt = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori_raw.fa",
        out_fasta_aa = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori_aa.fa"
    threads:
        5
    resources:
        runtime = "4h",
        mem_mb = 20000
    conda:
        "../envs/macse.yaml"
    group: "align_clean"
    log:
        "logs/macse2/{transcript_id}.log"
    shell:
        """
        macse \
        -prog alignSequences \
        -seq {input.in_fasta} \
        -out_AA {output.out_fasta_aa} \
        -out_NT {output.out_fasta_nt} \
        >> {log} 2>&1
        """

## Clean away stop codons
rule clean_macse_stop:
    input:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori_raw.fa"
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa",
    params:
        mincodon = 0,
        minseq =  0,
        minaalen = 0,
        mask = True,
    threads: 1,
    resources:
        runtime = "5m"
    group: "align_clean"
    log:
        "logs/clean_macse_stop/{transcript_id}.log"
    conda:
        "../envs/manual_cleaner.yaml"
    script:
        "../scripts/manual_filter_msa.py"
