## Here, codon alignment for all specified transcripts are created.
rule extract_ali:
    input:
        spec_names = "species.TOGA.dir.txt",
    output:
        codon_ali = "codon_alignments/{transcript_id}/tmp/{transcript_id}_ori.fa",
    params:
        transcript_id = "{transcript_id}",
        reference = config["referenceName"],
        additional =  "-t {}".format(config["settings"]["treeSettings"]["inputSpeciesTree"]["treeFile"]) if config["settings"]["alignmentSettings"]["aligner"] == "prank" and config["settings"]["treeSettings"]["inputSpeciesTree"]["treeFile"] is not None else " ",
        aligner = config["settings"]["alignmentSettings"]["aligner"],
        seed = config["settings"]["alignmentSettings"].get("prankBaseSeed", 12345),
        twoBit = config["settings"]["alignmentSettings"]["fromTOGA"]["twoBitPath"],
        ## ISSUE: TOGA container not yet available
        activate = config["settings"]["alignmentSettings"]["fromTOGA"].get("toga2Activate") or "",
        tmp_dir = "codon_alignments/{transcript_id}/tmp"
    threads:
        config["resources"]["extractAlignments"]["threads"],
    resources:
        mem_mb = config["resources"]["extractAlignments"]["mem_mb"],
        runtime = "10h",
    group: "align_clean"
    log:
        "logs/extract_ali/{transcript_id}.log"
    shell:
        """
        set -euo pipefail
        toga2 sequence-alignment \
              -v \
              -a {params.aligner} \
              -re {params.twoBit} \
              -asref {params.reference} \
              -l FI,I,PI \
              --input_dirs {input.spec_names} \
              --transcript_id {params.transcript_id} \
              {params.additional} \
              --seed {params.seed} \
              --tmp_dir {params.tmp_dir} \
              -o {output.codon_ali} \
              >> {log} 2>&1
        """
