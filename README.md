# easel

easel (EAsy SELection) is a [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline for large-scale screens of 
positive (episodic diversifying) and relaxed selection. Easel uses tools from the [HyPhy](https://hyphy.org) suite and [BayesCode](https://github.com/ThibaultLatrille/bayescode) and scales to hundreds of species and thousands of transcripts.

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

> **Under development:** 
> Easel has two modes: `-toga2` and `-free`. It can be seamlessly run with [TOGA2](https://github.com/hillerlab/TOGA2) output like this
> (with TOGA2 runs living in /path/to/genomes/hg38/TOGA2):
>```bash
># TOGA2 annotations, fixed species tree, aBSREL with synonymous rate variation
>easel -toga2 /path/to/genomes/hg38/TOGA2 \
>       -asm assemblies.txt \
>       -sb selected.bed \
>       --twoBit_path hg38.exons.2bit \
>       --reference_name hg38
>       -it species_tree.nh \
>       -a prank \
>       -ab srv \
>       -rs
>```

In Free Mode, you can pass your own unaligned or pre-aligned fasta files,
living in my_transcripts/. 
```bash
# your own FASTA files, nucleotide FASTA alignment (computed with Prank), per-gene trees (computed with IQtree), BUSTED with error-sink
easel -free my_transcripts/ \
       -a prank_nt \
       --reference_name hg38 \
       -ct \
       -bu srv,error_sink \
       -dr
```
Nothing runs unless you pass `-dr` (dry run) or `-rs`. Without either, easel
writes `DEF.yaml` and stops.

It is recommended to first perform a 'dry run' (`-dr` option in easel) to get an idea of how many jobs the pipeline will perform. Depending on number of input sequences, this can take some time (see Known Issues). When actually running the screen (`-rs`), it might again take some time for the pipeline to schedule all necessary jobs.

Since snakemake manages job submissions, easel is mostly intended to run on an HPC login node, through a **screen session** (Given that runs can take between minutes and several hours). Set it up like this:

```bash
screen -S easel_run
```

detaching from a running screen can be done with Ctrl + A + D, and reattaching works like this: 

```bash
screen -r easel_run ## if a screen with the name easel_run exists
```

> **Untested for remote submission:** 
Potentially, easel jobs can also be submitted to another node through `sbatch`, although there might arise conflicts for where snakemake will submit easel-spawned jobs, and sometimes remote nodes do not have internet access, leading to failed conda downloads and other cryptic errors.


## Test data
TBC

## Command Line Interface
TBC

## Tools used by easel

Easel automates alignment, gene-tree inference and selection screen using the latest state-of-the-art tools. If you use easel, please cite these tools accordingly.
A list of all relevant articles is available [here](https://github.com/casparbein/easel/docs/citations.md).

Orthology inference:
- [TOGA2](https://github.com/hillerlab/TOGA2)
<br/>


Alignment:
- [Prank](https://github.com/ariloytynoja/prank-msa) (Codon and Nucleotide modes)
- [Macse2](https://www.agap-ge2pop.org/macsee-pipelines/) (Nucleotide)
- [Muscle5](https://github.com/rcedgar/muscle) (Codon)
<br/>


Gene tree inference:
- [IQtree3](https://github.com/iqtree/iqtree3)
<br/>


Selection screens:
- [HyPhy](https://github.com/veg/hyphy)
  - [BUSTED](https://hyphy.org/methods/busted/)
  - [aBSREL](https://hyphy.org/methods/absrel/)
  - [MEME](https://hyphy.org/methods/meme/)
  - [RELAX](https://hyphy.org/methods/relax/)
- [BayesCode](https://github.com/ThibaultLatrille/bayescode)
<br/>


## Known issues
- easel's [TOGA2](https://github.com/hillerlab/TOGA2) mode is still under development and not yet functional.
- Depending on the number of input sequences, DAG creation can take up to ~20-30 Minutes (>5000 input sequences), for both dry runs and actual runs. 

