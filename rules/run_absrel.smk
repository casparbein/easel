## Perform an aBSREL run, parameters are defined in the config file
rule absrel_run:
    input:
        ali=get_input_ali_for_hyphy,
        tree=get_input_tree_for_hyphy
    output:
        "codon_alignments/{transcript_id}/{transcript_id}.g_tree.absrel.json" 
    threads:
        config["resources"]["absrel"]["threads"]
    resources:
        runtime = "4h",
        mem_mb = config["resources"]["absrel"]["mem_mb"]
    params:
        mh = config["settings"]["selectionSettings"]["ABSREL"]["multi_hit"],
        srv = config["settings"]["selectionSettings"]["ABSREL"]["srv"],
        foreground = " --branches Foreground " if config["settings"]["selectionSettings"]["foregroundLst"] else " "
    group: "absrel_stack"
    log:
        "logs/absrel_run/{transcript_id}.log",
    conda:
        "../envs/hyphy.yaml"
    shell:
        """
        set -euo pipefail
        hyphy \
        absrel \
        ENV='TOLERATE_NUMERICAL_ERRORS=1;' \
        --multiple-hits {params.mh} \
        --srv {params.srv} \
        --alignment {input.ali} \
        --tree {input.tree} \
        {params.foreground} \
        --output {output} \
        >> {log} 2>&1
        if [ ! -s {output} ]; then
            echo "absrel wrote no output for {wildcards.transcript_id}; see {log}" | tee -a {log} >&2
            exit 1
        fi
        if ! grep -q '"input"' {output}; then
            echo "absrel output for {wildcards.transcript_id} is not a usable HyPhy JSON; see {log}" | tee -a {log} >&2
            exit 1
        fi
        """