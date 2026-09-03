## Extract output tables from HyPhy Screens
rule extract_busted:
    input:
        json ="codon_alignments/{transcript_id}/{transcript_id}.busted.g_tree.json",
    output:
        busted_er = "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_er.g_tree.tsv",
        busted_model_fit = "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_model.g_tree.tsv",
        busted_tree = "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.busted_tree.g_tree.nh",
    params:
        script = f"{workflow.basedir}/scripts/extract_hyphy.py",
    group: "busted_stack"
    log:
        "logs/extract_busted/{transcript_id}.log"
    threads:
        config["resources"]["busted"]["threads"]
    resources:
        runtime = "15m",
        mem_mb = config["resources"]["busted"]["mem_mb"]
    shell:
        """
        {params.script} \
        busted \
        -o {output.busted_model_fit} \
        -b {output.busted_er} \
        -t {output.busted_tree} \
        {input.json} >> {log} 2>&1
        """

rule extract_absrel:
    input:
        json ="codon_alignments/{transcript_id}/{transcript_id}.g_tree.absrel.json",
    output:
        absrel_tree = "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.absrel_tree_g_tree.nh",
        absrel_er = "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.absrel_er_g_tree.tsv",
    params:
        script = f"{workflow.basedir}/scripts/extract_hyphy.py",
    group: "absrel_stack"
    log:
        "logs/extract_absrel/{transcript_id}.log"
    threads:
        config["resources"]["busted"]["threads"]
    resources:
        runtime = "15m",
        mem_mb = config["resources"]["busted"]["mem_mb"]
    shell:
        """
        {params.script} \
        absrel \
        -o {output.absrel_er} \
        -t {output.absrel_tree} \
        {input.json} >> {log} 2>&1
        """

rule extract_meme:
    input:
        json ="codon_alignments/{transcript_id}/{transcript_id}.meme.g_tree.json",
    output:
        meme_mle = "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.meme_mle.g_tree.tsv",
        meme_er = "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.meme_er.g_tree.tsv",
    params:
        script = f"{workflow.basedir}/scripts/extract_hyphy.py",
    group: "meme_stack"
    log:
        "logs/extract_meme/{transcript_id}.log"
    threads:
        config["resources"]["busted"]["threads"]
    resources:
        runtime = "15m",
        mem_mb = config["resources"]["busted"]["mem_mb"]
    shell:
        """
        {params.script} \
        meme \
        -o {output.meme_er} \
        -m {output.meme_mle} \
        {input.json} >> {log} 2>&1
        """

rule extract_relax:
    input:
        json=expand("codon_alignments/{{transcript_id}}/tmp/{{transcript_id}}.g_tree.relax_{it}.json",
        it=range(config["settings"]["selectionSettings"]["RELAX"]["relaxRounds"]))
    output:
        relax = "codon_alignments/{transcript_id}/HyPhy_output/{transcript_id}.relax.g_tree.tsv"
    params:
        script = f"{workflow.basedir}/scripts/extract_hyphy.py",
    group: "relax_post"
    log:
        "logs/extract_relax/{transcript_id}.log"
    threads: 1,
    resources:
        runtime = "15m",
        mem_mb = config["resources"]["busted"]["mem_mb"]
    shell:
        """
        {params.script} \
        relax \
        -o {output.relax} \
        {input.json} >> {log} 2>&1
        """