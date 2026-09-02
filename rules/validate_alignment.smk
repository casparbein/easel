## Gate each transcript after alignment and cleaning, before anything expensive
## depends on it.
##
## Empty or degenerate alignments were previously discovered by whichever
## downstream tool choked on them first -- HyPhy, IQ-TREE or the R plotter --
## which produced an error that named the wrong rule.
##
## check_alignment.py treats every problem it can see -- missing, empty,
## unreadable, unaligned, too few informative sequences, no foreground left,
## and any unexpected error in the checker itself -- as a SKIP verdict and
## exits 0, so a degenerate alignment produces a verdict rather than a failed
## job.
##
## This was a `checkpoint`, because rule all's file list is fixed when the DAG
## is built and only a checkpoint can change what is required mid-run. It is a
## plain rule now: the decision is read at PARSE time instead, from the verdicts
## phase 1 left on disk. See the "Two-phase run" section of rules/common.smk for
## why -- in short, a checkpoint made snakemake re-plan the entire DAG once per
## transcript, from inside the lock that blocks job submission.

from pathlib import Path

rule validate_alignment:
    input:
        ali = get_input_for_tree
    output:
        status = "codon_alignments/{transcript_id}/validation.txt"
    params:
        ## HyPhy needs >= 3 taxa; IQ-TREE bootstrapping wants >= 4.
        min_taxa = max(
            3 if not config["settings"]["treeSettings"]["computeGeneTrees"]["activate"] else 4,
            0,
        ),
        foreground = config["settings"]["selectionSettings"]["foregroundLst"] or [],
    threads: 1,
    resources:
        runtime = "10m",
    group: "align_clean"
    log:
        "logs/validate_alignment/{transcript_id}.log"
    script:
        "../scripts/check_alignment.py"


def alignment_only_target(transcript_id):
    """Terminal file of the alignment branch, mirroring gather_input()."""
    clean = config["settings"]["cleaningSettings"]
    aligner = config["settings"]["alignmentSettings"]["aligner"]
    base = f"codon_alignments/{transcript_id}/tmp/{transcript_id}"
    if clean["manualCleaning"]["doManualCleaning"]:
        return base + ".manual.fa"
    if clean["hmmCleaning"]["doHMMCleaning"]:
        return base + ".masked.hmm_cleaned.fa"
    if aligner in ("prank", "prank_nt", "muscle"):
        return base + ".masked.fa"
    return base + "_ori.fa"


def transcript_targets(wildcards):
    """Everything this transcript still owes, or nothing if it was skipped.
    """
    tid = wildcards.transcript_id
    if not is_usable(tid):
        return []

    extra = []
    if config["settings"]["alignmentSettings"]["aligner"] in ("prank_nt", "prank_codon"):
        extra.append(f"codon_alignments/{tid}/tmp/{tid}_ori_aa.fa")
    if config["settings"]["treeSettings"]["computeGeneTrees"]["activate"]:
        extra.append(f"codon_alignments/{tid}/{tid}_iqtree.nh")
    if config["settings"]["treeSettings"]["preCompGeneTrees"]["geneTreePath"] is not None:
        extra.append(f"codon_alignments/{tid}/{tid}_pruned_tree.nh")

    if config["settings"]["alignmentSettings"]["doAlignmentOnly"]:
        ## Alignment-only runs have no summary.pdf, so done.txt waits for the
        ## terminal file of the alignment branch instead.
        return [alignment_only_target(tid)] + extra
    return [f"codon_alignments/{tid}/{tid}.summary.pdf"] + extra


## How much the cleaning stages removed, for every transcript -- passed or
## rejected. A rejected alignment is exactly the case where you want to see it,
## which is why this keys off the cleaned alignment rather than the verdict: it
## runs for everything with an alignment (rule all requests it for `attempted`,
## i.e. OK plus SKIP), not only for what passed.
rule cleaning_report:
    input:
        ali = get_input_for_tree,
    output:
        plot   = "codon_alignments/{transcript_id}/{transcript_id}.cleaning.pdf",
        report = "codon_alignments/{transcript_id}/{transcript_id}.cleaning.txt",
    params:
        tmp_dir = "codon_alignments/{transcript_id}/tmp",
        script  = f"{workflow.basedir}/scripts/plot_cleaning.R",
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = "10m",
    group: "align_clean"
    log:
        "logs/cleaning_report/{transcript_id}.log"
    conda:
        "../envs/draw_output.yaml"
    shell:
        r"""
        if ! Rscript "{params.script}" \
              --dir "{params.tmp_dir}" \
              --id "{wildcards.transcript_id}" \
              --out "{output.plot}" \
              --report "{output.report}" \
              >> {log} 2>&1; then
            echo "cleaning report failed for {wildcards.transcript_id}; see {log}" | tee -a {log} >&2
            printf 'cleaning report could not be generated for %s\nsee %s\n' \
                "{wildcards.transcript_id}" "{log}" > "{output.report}"
            : > "{output.plot}"
        fi
        """


rule transcript_done:
    """Per-transcript completion marker: 'OK' or 'SKIP<TAB>reason'.

    rule all asks for these instead of for summary.pdf directly, so a skipped
    transcript resolves to zero downstream jobs and the workflow still finishes
    green.
    """
    input:
        status = "codon_alignments/{transcript_id}/validation.txt",
        targets = transcript_targets,
    output:
        "codon_alignments/{transcript_id}/done.txt"
    group: "finalize"
    resources:
        runtime = "5m",
    log:
        "logs/transcript_done/{transcript_id}.log"
    shell:
        "cp {input.status} {output} 2>> {log}"


rule skip_report:
    """One table of every transcript that did not make it, and why.

    Two sources:

      SKIP  the alignment was produced and rejected. done.txt exists, and
            carries the verdict -- read from the marker itself.
      FAIL  the alignment was never produced. Nothing downstream of it ran, so
            there is no done.txt at all; the row exists only in phase 1's
            verdict summary. Leaving these out would make a run that lost 200
            transcripts to aligner timeouts look like a clean one.
    """
    input:
        markers = expand("codon_alignments/{transcript_id}/done.txt",
                         transcript_id=attempted),
    output:
        "skipped_transcripts.tsv"
    resources:
        runtime = "15m",
    log:
        "logs/skip_report.log"
    run:
        rows = []
        for marker in input.markers:
            fields = Path(marker).read_text().strip().split("\t")
            if fields[0] != "OK":
                rows.append((Path(marker).parent.name, fields[0],
                             fields[1] if len(fields) > 1 else ""))
        rows.extend((tid, verdict, reason)
                    for tid, verdict, reason in verdict_summary_rows()
                    if verdict == VERDICT_FAIL)
        ## Input order, so the table reads the same way the BED does.
        order = {tid: i for i, tid in enumerate(transcripts)}
        rows.sort(key=lambda row: order.get(row[0], len(order)))

        with open(output[0], "w", newline="\n") as out, open(log[0], "w") as logfh:
            out.write("transcript\tverdict\treason\n")
            for tid, verdict, reason in rows:
                out.write(f"{tid}\t{verdict}\t{reason}\n")
            n_skip = sum(1 for row in rows if row[1] == VERDICT_SKIP)
            n_fail = sum(1 for row in rows if row[1] == VERDICT_FAIL)
            summary = (f"{len(rows)} of {len(transcripts)} transcript(s) produced "
                       f"no analysis: {n_skip} rejected by validation, {n_fail} "
                       f"did not align; see {output[0]}")
            print(summary)
            logfh.write(summary + "\n")
