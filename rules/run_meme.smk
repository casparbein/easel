## Perform an MEME run, parameters are defined in the config file
rule meme_run:
    input:
        ali=get_input_ali_for_hyphy,
        tree =get_input_tree_for_hyphy
    output:
        "codon_alignments/{transcript_id}/{transcript_id}.meme.g_tree.json" 
    threads:
        config["resources"]["meme"]["threads"]
    resources:
        runtime = "4h",
        mem_mb = config["resources"]["meme"]["mem_mb"]
    params:
        mh = config["settings"]["selectionSettings"]["MEME"]["multi_hit"],
        srv = config["settings"]["selectionSettings"]["MEME"]["srv"],
        foreground = " --branches Foreground " if config["settings"]["selectionSettings"]["foregroundLst"] else " "
    group: "meme_stack"
    log:
        "logs/meme_run/{transcript_id}.log"
    conda:
        "../envs/hyphy.yaml"
    shell:
        """
        set -euo pipefail
        hyphy \
        meme \
        ENV='TOLERATE_NUMERICAL_ERRORS=1;' \
        --multiple-hits {params.mh} \
        --srv {params.srv} \
        --alignment {input.ali} \
        --tree {input.tree} \
        {params.foreground} \
        --output {output} \
        >> {log} 2>&1
        """