# eaSel

Run selection screens with the [HyPhy](https://hyphy.org) suite over thousands
of coding sequences from hundreds of species. eaSel (EA**S**y **S**ELection) is a 
snakemake pipeline that extracts or ingests per-transcript alignments, cleans them,
optionally builds trees, and runs aBSREL, BUSTED, MEME, RELAX and BayesCode.

> **Status:** pre-release. 
Feel free to open a GitHub issue whenever something does not work properly

## Install

```bash
git clone https://github.com/casparbein/eaSel.git
cd easel
conda env create -f environment.yaml
conda activate easel
pip install -e .
```

`easel --help` should now work.


Everything else comes from the per-rule environments in `envs/`.

## Quick start

```bash
# TOGA2 annotations, fixed species tree, aBSREL with synonymous rate variation
easel -toga2 /path/to/genomes/hg38/TOGA2 \
       -asm assemblies.txt \
       -sb selected.bed \
       --twoBit_path hg38.exons.2bit \
       -it species_tree.nh \
       -a prank \
       -ab srv \
       -rs

# your own FASTA files, per-gene trees, BUSTED with error-sink
easel -free my_transcripts/ -a prank_nt --reference_name hg38 \
       -ct -bu srv,error_sink -rs
```

Nothing runs unless you pass `-dr` (dry run) or `-rs`. Without either, easel
writes `DEF.yaml` and stops.
