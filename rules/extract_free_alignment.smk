## Here, existing alignments are copied over to the working directory, in free path mode
# rule extract_ali:
#     input:
#         in_fasta = f"{config['fastaPath']}{{transcript_id}}{config['fileSuffix']}"
#     output:
#         out_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa",
#     localrule: True,
#     threads: 1,
#     resources:
#         runtime = "5m",
#     log:
#         "logs/extract_ali/{transcript_id}.log"
#     shell:
#         "cp {input.in_fasta} {output.out_fasta} 2>> {log}"

rule extract_ali:
    input:
        f"{config['fastaPath']}{{transcript_id}}{config['fileSuffix']}"
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
        "logs/extract_ali/{transcript_id}.log"
    conda:
        "../envs/manual_cleaner.yaml"
    script:
        "../scripts/manual_filter_msa.py"