## Renaming and codonification -- turns the raw aligner output into the
## masked, codon-numbered alignment the rest of the pipeline works from.
## Split out of alignment_postprocessing.smk; see rules/hmm_cleaning.smk and
## rules/manual_cleaning.smk for the stages downstream of this one.

## Input function for codonifying alignments
def get_input_codonify(wildcards):
    ## TOGA2 and free mode both reduced to the same muscle-vs-other test; the
    ## third (implicit TOGA v1) case ignored the aligner and always returned
    ## _ori.fa, which -toga does too by having no working rule chain at all.
    if config["settings"]["alignmentSettings"]["aligner"] == "muscle" and config["freeMode"]:
        ## This is a temporary adjustment for Muscle testing, not yet implemented
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}_muscle_filtered.fa"
    else:
        return "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa"


## Temporary rule: Rename REFERNECE to config["referenceName"] because ref in species tree will likely not be named REFERENCE and
## otherwise TOGA2 alignment extraction will not work.
rule rename_reference:
    input:
        get_input_codonify,
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}_ren.fa"
    params:
        ref_name = config["referenceName"],
    resources:
        runtime = "5m",
    group: "align_clean"
    log:
        "logs/rename_reference/{transcript_id}.log"
    shell:
        """
        cat {input} | sed 's/REFERENCE/{params.ref_name}/g' > {output} 2>> {log}
        """

## Codonify should be here as it should also be applied to Free alignments
rule codonify_ali:
    input:
        ali="codon_alignments/{transcript_id}/tmp/{transcript_id}_ren.fa",
    output:
        codon_ali = "codon_alignments/{transcript_id}/tmp/{transcript_id}.masked.fa",
        premask_ali = "codon_alignments/{transcript_id}/tmp/{transcript_id}.premasked.fa",
        frameshifts = "codon_alignments/{transcript_id}/tmp/{transcript_id}.frameshifts.txt"
    params:
        reference = config["referenceName"],
        out_head = "codon_alignments/{transcript_id}/tmp/{transcript_id}"
    threads:
        config["resources"]["hmmCleaner"]["threads"]
    resources:
        mem_mb = config["resources"]["hmmCleaner"]["mem_mb"],
        runtime = "15m"
    ## ISSUE: Update to newest codonify script
    ## ISSUE: codonify is currently local
    group: "align_clean"
    log:
        "logs/codonify_ali/{transcript_id}.log"
    shell:
        """
        codonify_prank_alignment.py \
        -i {input} \
        -o {params.out_head} \
        -r {params.reference} \
        >> {log} 2>&1
        """
