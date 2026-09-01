## Keeping a run's storage footprint small, without making the headline
## results unreadable.
##
## Three things happen here:
##   1. the final cleaned alignment is copied OUT of tmp/ before tmp/ is
##      archived, so the one file people actually want is not inside a tarball
##   2. a plain-text selection summary is written per transcript
##   3. HyPhy_output/ and the HyPhy JSONs -- by far the largest outputs -- are
##      archived and removed
##
## Archive-then-delete is split into two rules for the same reason
## tar_compress/clean_tmp are (rules/draw_output.smk): the archive must be a
## committed output before anything is destroyed, so a discarded archive is
## always recoverable from the originals.


## The final alignment is whatever the screens actually ran on -- .manual.fa,
## .masked.hmm_cleaned.fa, .masked.fa or _ori.fa depending on the aligner and
## which cleaning steps are on. get_input_for_tree (rules/common.smk) already
## encodes that choice, so this does not re-derive it.
rule keep_final_alignment:
    input:
        ali = get_input_for_tree,
    output:
        ali = "codon_alignments/{transcript_id}/{transcript_id}.final.fa",
    group: "align_clean"
    resources:
        runtime = "5m",
    log:
        "logs/keep_final_alignment/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        cp "{input.ali}" "{output.ali}" 2>> {log}
        echo "kept $(basename {input.ali}) as {output.ali}" >> {log}
        """

def _selection_report_args(wildcards):
    sel = config["settings"]["selectionSettings"]
    tid = wildcards.transcript_id
    base = f"codon_alignments/{tid}"
    args = []
    if sel["ABSREL"]["activate"]:
        args.append(f"--absrel-tsv {base}/HyPhy_output/{tid}.absrel_er{ABSREL_KEYWORD}.tsv")
    if sel["BUSTED"]["activate"]:
        args.append(f"--busted-json {base}/{tid}.busted{HYPHY_KEYWORD}.json")
    if sel["MEME"]["activate"]:
        args.append(f"--meme-tsv {base}/HyPhy_output/{tid}.meme_mle{HYPHY_KEYWORD}.tsv")
    if sel["RELAX"]["activate"]:
        args.append(f"--relax-tsv {base}/HyPhy_output/{tid}.relax{HYPHY_KEYWORD}.tsv")
    return " ".join(args)


rule selection_report:
    input:
        done = "codon_alignments/{transcript_id}/done.txt",
        status = "codon_alignments/{transcript_id}/validation.txt",
    output:
        report = "codon_alignments/{transcript_id}/{transcript_id}.selection.txt",
    group: "finalize"
    params:
        script = f"{workflow.basedir}/scripts/selection_report.py",
        python = sys.executable,
        args = _selection_report_args,
        ## Hand-editable in DEF.yaml; there is no CLI flag for it.
        pvalue = config["settings"]["selectionSettings"].get("reportPvalue", 0.05),
    resources:
        runtime = "10m",
    log:
        "logs/selection_report/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        if ! "{params.python}" "{params.script}" \
              --transcript "{wildcards.transcript_id}" \
              --out "{output.report}" \
              --status "{input.status}" \
              --pvalue {params.pvalue} \
              {params.args} \
              >> {log} 2>&1; then
            echo "selection report failed for {wildcards.transcript_id}; see {log}" | tee -a {log} >&2
            printf 'selection report could not be generated for %s\nsee %s\n' \
                "{wildcards.transcript_id}" "{log}" > "{output.report}"
        fi
        """


## HyPhy_output/ and the *.json files are the bulk of a run's disk usage. They
## are archived only after selection_report and draw_output have read them.
rule tar_results:
    input:
        report = "codon_alignments/{transcript_id}/{transcript_id}.selection.txt",
        summary = "codon_alignments/{transcript_id}/done.txt",
    output:
        archive = "codon_alignments/{transcript_id}/results.tar.gz",
    params:
        parent = "codon_alignments/{transcript_id}",
    resources:
        runtime = "15m",
    log:
        "logs/tar_results/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        PARENT="{params.parent}"
        LIST=$(mktemp)
        trap 'rm -f "$LIST"' EXIT
        if [ -d "$PARENT/HyPhy_output" ]; then
            echo "HyPhy_output" >> "$LIST"
        fi
        find "$PARENT" -maxdepth 1 -name '*.json' -printf '%f\n' >> "$LIST"
        if [ ! -s "$LIST" ]; then
            echo "nothing to archive for {wildcards.transcript_id}" >> {log}
            tar -czf {output.archive} -T /dev/null
            exit 0
        fi
        echo "archiving:" >> {log}
        sed 's/^/  /' "$LIST" >> {log}
        tar -C "$PARENT" -czf {output.archive}.part -T "$LIST" 2>> {log}
        tar -tzf {output.archive}.part > /dev/null
        mv {output.archive}.part {output.archive}
        """


## Deletion is its own job, and re-verifies the archive first: this is the last
## moment the originals still exist.
rule clean_results:
    input:
        archive = "codon_alignments/{transcript_id}/results.tar.gz",
    output:
        marker = "codon_alignments/{transcript_id}/results.cleaned",
    params:
        parent = "codon_alignments/{transcript_id}",
    resources:
        runtime = "10m",
    log:
        "logs/clean_results/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        ARCHIVE="{input.archive}"
        verified=no
        for attempt in 1 2 3 4 5; do
            if tar -tzf "$ARCHIVE" > /dev/null 2>> {log}; then
                verified=yes
                break
            fi
            if [ "$attempt" -lt 5 ]; then
                echo "attempt $attempt: $ARCHIVE not readable yet, waiting 10s" >> {log}
                sleep 10
            fi
        done
        if [ "$verified" != yes ]; then
            echo "archive still unreadable after 5 attempts; NOT deleting anything" | tee -a {log} >&2
            echo "originals are intact. Remove the bad archive and rerun:" | tee -a {log} >&2
            echo "  rm $ARCHIVE" | tee -a {log} >&2
            exit 1
        fi
        rm -rf "{params.parent}/HyPhy_output" 2>> {log}
        find "{params.parent}" -maxdepth 1 -name '*.json' -delete 2>> {log}
        echo "removed HyPhy_output/ and *.json (archive verified)" >> {log}
        date -u +"%Y-%m-%dT%H:%M:%SZ" > {output.marker}
        """
