def get_activated_inputs(wildcards, step):
    inputs = []

    mapping = {"BUSTED": ["codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_er.g_tree.tsv", 
                            "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_tree.g_tree.nh"],
               "MEME":  "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.meme_mle.g_tree.tsv",
               "ABSREL": ["codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.absrel_er_g_tree.tsv",
               "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.absrel_tree_g_tree.nh" ],
               "bayesCode": ["codon_alignments/{transcript_id}/tmp/{transcript_id}_classic.g_tree.ci0.025.tsv", 
                            "codon_alignments/{transcript_id}/tmp/{transcript_id}_mutsel.g_tree.ci0.025.tsv"],
               "RELAX": "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.relax.g_tree.tsv"}
    
    for feature, path in mapping.items():
        if feature == step:
            #print(config["settings"]["selectionSettings"].get(feature))
            if config["settings"]["selectionSettings"].get(feature)["activate"]:
                inputs.append(path)
            
    return inputs

## Extract output tables from HyPhy Screens
rule draw_output:
    input:
        asbrel = get_activated_inputs("{transcript_id}", "ABSREL"),
        busted = get_activated_inputs("{transcript_id}", "BUSTED"),
        meme_mle = get_activated_inputs("{transcript_id}", "MEME"),
        bayescode = get_activated_inputs("{transcript_id}", "bayesCode"),
        alignment = get_input_ali_for_hyphy,
        relax_k = get_activated_inputs("{transcript_id}", "RELAX")
    output:
        summary_plot = "codon_alignments/{transcript_id}/{transcript_id}.summary.pdf",
        overview = "codon_alignments/{transcript_id}/{transcript_id}.overview.pdf"
    group: "finalize"
    log:
        "logs/draw_output/{transcript_id}.log"
    threads:
        1
    resources:
        runtime = "30m",
        mem_mb = 5000
    conda:
        "../envs/draw_output.yaml"
    script:  
        "../scripts/draw_output.R"

rule tar_compress:
    input:
        done = "codon_alignments/{transcript_id}/done.txt",
        cleaning = (lambda w: f"codon_alignments/{w.transcript_id}/{w.transcript_id}.cleaning.txt"
                    if WANT_CLEANING_REPORT else []),
    output:
        archive = "codon_alignments/{transcript_id}/tmp.tar.gz",
    params:
        parent = "codon_alignments/{transcript_id}",
    resources:
        runtime = "15m",
    log:
        "logs/tar_compress/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        if [ ! -d "{params.parent}/tmp" ]; then
            echo "no tmp/ directory for {wildcards.transcript_id}" >> {log}
            tar -czf {output.archive} -T /dev/null
            exit 0
        fi
        tar -C {params.parent} -czf {output.archive}.part tmp 2>> {log}
        tar -tzf {output.archive}.part > /dev/null
        mv {output.archive}.part {output.archive}
        echo "archived {params.parent}/tmp" >> {log}
        """

## Deleting tmp/ is its own rule, and its own job, so the destructive step runs
## only once the archive above is a committed output.
rule clean_tmp:
    input:
        archive = "codon_alignments/{transcript_id}/tmp.tar.gz",
        final = "codon_alignments/{transcript_id}/{transcript_id}.final.fa",
    output:
        marker = "codon_alignments/{transcript_id}/tmp.cleaned",
    params:
        parent = "codon_alignments/{transcript_id}",
    resources:
        runtime = "10m",
    log:
        "logs/clean_tmp/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        ARCHIVE="{input.archive}"
        TMP="{params.parent}/tmp"
        echo "verifying $ARCHIVE" >> {log}
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
            echo "archive still unreadable after 5 attempts; NOT deleting $TMP" | tee -a {log} >&2
            echo "tmp/ is intact, so nothing is lost. Remove the bad archive and" | tee -a {log} >&2
            echo "rerun to rebuild it:  rm $ARCHIVE" | tee -a {log} >&2
            exit 1
        fi
        if [ -d "$TMP" ]; then
            rm -rf "$TMP" 2>> {log}
            echo "removed $TMP (archive verified)" >> {log}
        else
            echo "$TMP already gone; nothing to remove" >> {log}
        fi
        date -u +"%Y-%m-%dT%H:%M:%SZ" > {output.marker}
        """