## Align with MACSE2
rule extract_ali_macse2:
    input:
        in_fasta = f"{config['fastaPath']}{{transcript_id}}{config['fileSuffix']}"
    output:
        out_fasta_nt = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa",
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
