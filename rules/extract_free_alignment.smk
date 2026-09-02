## Here, existing alignments are copied over to the working directory, in free path mode
rule extract_ali:
    input:
        in_fasta = f"{config['fastaPath']}{{transcript_id}}{config['fileSuffix']}"
    output:
        out_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa",
    localrule: True,
    threads: 1,
    resources:
        runtime = "5m",
    log:
        "logs/extract_ali/{transcript_id}.log"
    shell:
        "cp {input.in_fasta} {output.out_fasta} 2>> {log}"