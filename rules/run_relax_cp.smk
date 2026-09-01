## Perform Relax run, MH and SRV can be set in the config file. Runs of 10 batches should be performed for consistency reasons, but that can be adapted later
rule relax_run: 
    input:
        ali = get_input_ali_for_hyphy,
        tree = get_input_tree_for_hyphy
    output:
        "codon_alignments/{transcript_id}/tmp/{transcript_id}.g_tree.relax_{it}.json" 
    threads:
        config["resources"]["relax"]["threads"]
    group: "relax_stack"
    resources:
        runtime = "4h",
        mem_mb = config["resources"]["relax"]["mem_mb"]
    params:
        mh = config["settings"]["selectionSettings"]["RELAX"]["multi_hit"],
        srv = config["settings"]["selectionSettings"]["RELAX"]["srv"]
    log:
        "logs/relax_run/{transcript_id}_{it}.log"
    conda:
        "../envs/hyphy.yaml"
    shell:
        """
        set -euo pipefail
        hyphy \
        relax \
        ENV='TOLERATE_NUMERICAL_ERRORS=1;' \
        --alignment {input.ali} \
        --tree {input.tree}  \
        --srv {params.srv} \
        --multiple-hits {params.mh} \
        --test Foreground \
        --output {output} \
        >> {log} 2>&1
        """

## aggregate relax runs:
rule aggregate_relax:
    input:
        expand("codon_alignments/{{transcript_id}}/tmp/{{transcript_id}}.g_tree.relax_{it}.json",
        it=range(config["settings"]["selectionSettings"]["RELAX"]["relaxRounds"]))
    output:
        "codon_alignments/{transcript_id}/{transcript_id}.g_tree.relax.aggregated.json"
    group: "relax_stack"
    resources:
        runtime = "5m",
    log:
        "logs/aggregate_relax/{transcript_id}.log"
    shell:
        "cat {input} > {output} 2>> {log}"
