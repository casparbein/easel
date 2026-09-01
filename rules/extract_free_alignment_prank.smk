## Align with PRANK
## Add input tree if available
rule extract_ali_prank:
    input:
        in_fasta = f"{config['fastaPath']}{{transcript_id}}{config['fileSuffix']}"
    output:
        out_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.best.fas",
    params:
        codon_mode = " -codon " if config["settings"]["alignmentSettings"]["aligner"] == "prank_codon" else " ",
        tree = f'-t={config["settings"]["treeSettings"]["inputSpeciesTree"]["treeFile"]} -prunetree -prunedata -once' if config["settings"]["treeSettings"]["inputSpeciesTree"]["treeFile"] else " -iterate=10 ",
        out_head = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori",
        seed = config["settings"]["alignmentSettings"].get("prankBaseSeed", 12345),
        timeout = "5h",
    conda:
        "../envs/prank.yaml"
    threads:
        config["resources"]["prank"]["threads"]
    resources:
        mem_mb = config["resources"]["prank"]["mem_mb"],
        runtime = "1h",
    group: "align_clean"
    log:
        "logs/extract_ali_prank/{transcript_id}.log"
    shell:
        """
        prank \
        -d={input.in_fasta} \
        -o={params.out_head} \
        {params.codon_mode} \
        {params.tree} \
        -seed={params.seed} \
        -F \
        -f=fasta \
        < /dev/null >> {log} 2>&1 
        """


rule rename_prank:
    input:
        in_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.best.fas"
    output:
        out_fasta = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa"
    group: "align_clean"
    resources:
        runtime = "5m",
    log:
        "logs/rename_prank/{transcript_id}.log",
    shell:
        "cp {input.in_fasta} {output.out_fasta} 2>> {log}"


rule translate_prank_aa:
    input:
        get_uncleaned_alignment,
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori_aa.fa",
    threads:
        config["resources"]["hmmCleaner"]["threads"]
    resources:
        mem_mb = config["resources"]["hmmCleaner"]["mem_mb"],
        runtime = "10m"
    group: "align_clean"
    log:
        "logs/translate_prank_aa/{transcript_id}.log"
    script:
        "../scripts/transform_aa.py"
