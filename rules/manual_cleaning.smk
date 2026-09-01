## Manual cleaning: remove alignment columns that are mostly gaps or highly
## dissimilar. Split out of alignment_postprocessing.smk.

## Input function for manual_cleaner:
def get_input_for_manual(wildcards):
    if config["settings"]["cleaningSettings"]["hmmCleaning"]["doHMMCleaning"]:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.hmm_cleaned.fa"
    else:
        return get_uncleaned_alignment(wildcards)

## This is the manual cleaning step removing alignment columns that are mostly gaps or highly dissimilar
rule manual_cleaner:
    input:
        get_input_for_manual
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}.manual.fa",
    params:
        mincodon = config["settings"]["cleaningSettings"]["manualCleaning"]["manualCleanerParams"]["mincodon"],
        minseq =  config["settings"]["cleaningSettings"]["manualCleaning"]["manualCleanerParams"]["minseq"],
        minaalen = config["settings"]["cleaningSettings"]["manualCleaning"]["manualCleanerParams"]["minaalen"],
        mask = config["settings"]["cleaningSettings"]["manualCleaning"]["manualCleanerParams"]["mask"],
    threads:
        config["resources"]["manualCleaner"]["threads"]
    resources:
        mem_mb = config["resources"]["manualCleaner"]["mem_mb"],
        runtime = "15m"
    group: "align_clean"
    log:
        "logs/manual_cleaner/{transcript_id}.log"
    conda:
        "../envs/manual_cleaner.yaml"
    script:
        "../scripts/manual_filter_msa.py"