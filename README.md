# easel

easel (EAsy SELection) is a [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline for large-scale screens of 
positive (episodic diversifying) and relaxed selection. Easel uses tools from the [HyPhy](https://hyphy.org) suite and scales to hundreds of species and thousands of 
transcripts.

> **Status:** pre-release. 
Feel free to open a GitHub issue whenever something does not work properly

## Requirements
While easel has a `--local` mode, it's best used on an HPC. \
Currently, easel requires the slurm scheduling system and mamba or conda (see "Install").

## Install

```bash
git clone https://github.com/casparbein/eaSel.git
cd easel
conda env create -f environment.yaml
conda activate easel
pip install -e .
```

`easel --help` should now work.

All other necessary tools will be automatically provided by snakemake through containers and conda environments.\
To run easel on your own HPC with slurm, you only have to change the name of the partition in `prof/config.yaml`:

```yaml
executor: slurm
jobs: 1000
slurm-delete-logfiles-older-than: 0
slurm-keep-successful-logs: True
slurm-no-requeue: True
use-conda: True
use-apptainer: True
restart-times: 1
rerun-incomplete: True
printshellcmds: True
scheduler: greedy
latency-wait: 120
default-resources:
    slurm_partition: "batch" #<- change here for your own partition 
    runtime: 24h #<- change according to what your partition allows
    mem_mb: 1000 #<- hange according to what your partition allows
```

> **Untested for other HPC environments:**
Alternatively, you can pass your own `my_config.yaml` file to easel, but it has to be configured for slurm:
>```bash
>easel -free my_transcripts/ \
>       --config my_config.yaml \
>       -a prank_codon
>       -ct \
>       -bu srv \
>       -rs
>```

## Pipeline overview
<img width="2704" height="2850" alt="easel_pipeline" src="https://github.com/user-attachments/assets/02964856-502d-4e4f-b156-239570e8381a" />

## Quick start

Easel has two modes: `-toga2` and `-free`. It can be seamlessly run with [TOGA2](https://github.com/hillerlab/TOGA2) output like this
(with TOGA2 runs living in /path/to/genomes/hg38/TOGA2):
```bash
# TOGA2 annotations, fixed species tree, aBSREL with synonymous rate variation
easel -toga2 /path/to/genomes/hg38/TOGA2 \
       -asm assemblies.txt \
       -sb selected.bed \
       --twoBit_path hg38.exons.2bit \
       --reference_name hg38
       -it species_tree.nh \
       -a prank \
       -ab srv \
       -rs
```

'Free' means that you can pass your own unaligned or pre-aligned fasta files,
living in my_transcripts/. 
```bash
# your own FASTA files, per-gene trees, BUSTED with error-sink
easel -free my_transcripts/ \
       -a prank_nt \
       --reference_name hg38 \
       -ct \
       -bu srv,error_sink \
       -rs
```

Nothing runs unless you pass `-dr` (dry run) or `-rs`. Without either, easel
writes `DEF.yaml` and stops.

## Test data
TBC

## Command Line Interface
TBC

## Tools used by easel
TBC
