## Perform a busted run, error-sink can be turned on and off in the config file
## output names might have to be changed
rule busted_run:
    input:
        ali= get_input_for_busted,
        tree = get_input_tree_for_hyphy
    output:
        "codon_alignments/{transcript_id}/{transcript_id}.busted.g_tree.json"
    params:
        error_sink = "Yes" if config["settings"]["selectionSettings"]["BUSTED"]["error_sink"] else "No",
        mh = config["settings"]["selectionSettings"]["BUSTED"]["multi_hit"],
        srv = config["settings"]["selectionSettings"]["BUSTED"]["srv"],
        foreground = " --branches Foreground " if config["settings"]["selectionSettings"]["foregroundLst"] and not config["settings"]["selectionSettings"]["doErrorCleaning"] else " "
    group: "busted_stack"
    log:
        "logs/busted_run/{transcript_id}.log"
    threads:
        config["resources"]["busted"]["threads"]
    resources:
        runtime = "4h",
        mem_mb = config["resources"]["busted"]["mem_mb"]
    conda:
        "../envs/hyphy.yaml"
    shell:
        """
        set -euo pipefail
        hyphy \
        busted \
        ENV='TOLERATE_NUMERICAL_ERRORS=1;' \
        --alignment {input.ali} \
        --tree {input.tree} \
        --error-sink {params.error_sink} \
        --srv {params.srv} \
        --multiple-hits {params.mh} \
        {params.foreground} \
        --output {output} \
        >> {log} 2>&1
        """

## Create error-corrected alignment from busted run
rule error_mask:
    input:
        json = "codon_alignments/{transcript_id}/{transcript_id}.busted.g_tree.json"
    output:
        json = "codon_alignments/{transcript_id}/{transcript_id}.busted.g_tree.filtered.json",
        ali = temp("codon_alignments/{transcript_id}/{transcript_id}.busted.g_tree.filtered.tmp.fa")
    threads:
        config["resources"]["busted"]["threads"]
    resources:
        runtime = "30m",
        mem_mb = config["resources"]["busted"]["mem_mb"]
    group: "busted_stack"
    log:
        "logs/error_mask/{transcript_id}.log"
    conda:
        "../envs/hyphy.yaml"
    shell:
        """
        hyphy \
        error-filter \
        ENV='TOLERATE_NUMERICAL_ERRORS=1;' \
        --json {input.json} \
        --output {output.ali} \
        --output-json {output.json} \
        >> {log} 2>&1
        """

rule clean_error_sink_ali:
    input:
        ali = "codon_alignments/{transcript_id}/{transcript_id}.busted.g_tree.filtered.tmp.fa"
    output:
        ali = "codon_alignments/{transcript_id}/tmp/{transcript_id}.busted.g_tree.filtered.fa"
    group: "busted_stack"
    resources:
        runtime = "5m",
    log:
        "logs/clean_error_sink_ali/{transcript_id}.log"
    shell:
        r"""
        set -euo pipefail
        grep -v '(' {input.ali} > {output.ali} 2>> {log} || true
        """

## TO-DO: properly implement BUSTED error sink module