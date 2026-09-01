## HmmCleaner stage: translate to amino acids, run HmmCleaner.pl, transfer
## its cleaning decisions back onto the nucleotide alignment.
## Split out of alignment_postprocessing.smk; get_uncleaned_alignment lives
## in rules/common.smk since rules/manual_cleaning.smk needs it too.

## Convert nucleotide alignment to amino acids for HMM cleaner. For a
## codon-aware aligner (macse2/prank_codon) that is the aligner's direct
## output.
rule convert_to_aa:
    input:
        get_uncleaned_alignment,
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.aa.fa",
    threads:
        config["resources"]["hmmCleaner"]["threads"]
    resources:
        mem_mb = config["resources"]["hmmCleaner"]["mem_mb"],
        runtime = "5m"
    group: "align_clean"
    log:
        "logs/convert_to_aa/{transcript_id}.log"
    script:
        "../scripts/transform_aa.py"


## This is cleaning step where parts of alignment ROWS (stretches from one input sequence) are removed based on HMM predictions
## Note that HMM cleaner works on amino acid sequences only, so the results have to be transferred to the original nucleotide alignment
rule hmm_cleaner:
    input:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.aa.fa",
    output:
        aa_ali = "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.aa_hmm.fasta",
        aa_log = "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.aa_hmm.log",
    params:
        cost = config["settings"]["cleaningSettings"]["hmmCleaning"]["hmmCleanerParams"]
    threads:
        config["resources"]["hmmCleaner"]["threads"]
    resources:
        mem_mb = config["resources"]["hmmCleaner"]["mem_mb"],
        runtime = "20m"
    ## HmmCleaner is not on bioconda
    container:
        "docker://ghcr.io/hillerlab/hmmcleaner:latest"
    group: "align_clean"
    log:
        "logs/hmm_cleaner/{transcript_id}.log"
    shell:
        """
        HmmCleaner.pl \
        -costs {params.cost} \
        {input} \
        >> {log} 2>&1
        """

## Transfer HMM cleaner results back onto the nucleotide alignment.
rule transfer_cleaner:
    input:
        log="codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.aa_hmm.log",
        ali=get_uncleaned_alignment,
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked_cleaned.ali",
    threads:
        config["resources"]["hmmCleaner"]["threads"]
    resources:
        mem_mb = config["resources"]["hmmCleaner"]["mem_mb"],
        runtime = "20m"
    container:
        "docker://ghcr.io/hillerlab/hmmcleaner:latest"
    group: "align_clean"
    log:
        "logs/transfer_cleaner/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        ALI="{input.ali}"
        PRODUCED="${{ALI%.fa}}_cleaned.ali"
        transferCleaner.pl \
        -log={input.log} \
        "$ALI" \
        -delchar "N" \
        >> {log} 2>&1
        if [ ! -f "$PRODUCED" ]; then
            echo "transferCleaner.pl exited 0 but did not write $PRODUCED; directory contains:" | tee -a {log} >&2
            ls -1 "$(dirname "$ALI")" | tee -a {log} >&2
            exit 1
        fi
        if [ "$PRODUCED" != "{output}" ]; then
            mv "$PRODUCED" "{output}"
        fi
        """

## Transfer Cleaner creates ## in the header of the alignment file, must be removed
rule clean_transfer_cleaner:
    input:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked_cleaned.ali",
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.hmm_cleaned.fa",
    threads:
        config["resources"]["hmmCleaner"]["threads"]
    resources:
        mem_mb = config["resources"]["hmmCleaner"]["mem_mb"],
        runtime = "5m"
    group: "align_clean"
    log:
        "logs/clean_transfer_cleaner/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        grep -v "#" {input} > {output} 2>> {log} || true
        if ! grep -q '^>' {output}; then
            echo "no sequences left after HMM-cleaner transfer for {input}" | tee -a {log} >&2
            exit 1
        fi
        """
