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
git clone https://github.com/casparbein/easel.git
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

Alternatively, you can pass your own `my_config.yaml` file to easel, but it has to be configured for slurm:
```bash
easel -free my_transcripts/ \
       --profile my_config.yaml \
       -a prank_codon
       -ct \
       -bu srv \
       -rs
```

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

### TOGA2
TOGA2 runs with different references are available [here](https://genome.senckenberg.de/download/TOGA2). \
For an easel run with TOGA2 annotations as input for the alignment, we need:

- TOGA2 result directories (one per assembly analyzed) 
- a reference exon 2bit file (fetch here)
- a bed file with transcripts we want to analyse (fetch here)
- a list of assemblies that we want to analyse (based on the TOGA2 result directory names)

For a quick easel test, we'll download 6 TOGA annotations based on the human (hg38) reference and save them into toga2_runs in our working directory.
We then also save the human exon 2bit file (used for exon-by-exon alignment) and a bed file with reference annotations (in our test file, we have 20),
and we create a list of assemblies we want to analyse (in this case all test assemblies we downloaded).

```bash
## Create a directory where to save the input files
mkdir easel_test_run_toga

## Create a directory where to download TOGA2 result directories into (structure is reference/TOGA2)
mkdir -p hg38/TOGA2
cd hg38/TOGA2/
wget -r -nH --cut-dirs=4 --no-parent -R 'index.html*' https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Acinonyx_jubatus__cheetah__HLaciJub2__GCA_003709585.1/
wget -r -nH --cut-dirs=4 --no-parent -R 'index.html*' https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Bison_bonasus__European_bison__HLbisBon1__GCA_963879515.1/
wget -r -nH --cut-dirs=4 --no-parent -R 'index.html*' https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Canis_lupus_familiaris__dog__canFam4__GCA_011100685.1/
wget -r -nH --cut-dirs=4 --no-parent -R 'index.html*' https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Choloepus_didactylus__southern_two-toed_sloth__HLchoDid2__GCF_015220235.1/
wget -r -nH --cut-dirs=4 --no-parent -R 'index.html*' https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Equus_caballus__Horse__HLequCaba5__GCA_052818215.1/
wget -r -nH --cut-dirs=4 --no-parent -R 'index.html*' https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Gorilla_gorilla__western_gorilla__HLgorGor7__GCA_029281585.3/
rm robots.txt

## Create assembly list from directory names
ls > ../../easel_test_run_toga/assembly_list.txt
cd ../..

## Download bed file () and 2bit file
cd easel_test_run_toga
wget https://raw.githubusercontent.com/casparbein/easel/refs/heads/main/test_input/toga2/hg38.test.bed
wget https://raw.githubusercontent.com/casparbein/easel/refs/heads/main/test_input/toga2/hg38.test.exons.2bit
cd ..

## Run easel on TOGA2 test data:
## Aligner: Prank (Nucleotide); Tree: Compute Gene Trees; Selection Screen: BUSTED+SRV
easel \
-toga2 hg38/TOGA2/ \
-a prank \
-ct \
-bu srv \
-twoBit easel_test_run_toga/hg38.exon.2bit \
-asm easel_test_run_toga/assembly_list.txt \
-sb easel_test_run_toga/hg38.toga.transcripts.bed \
-rs
```

### Free (any coding FASTA files)
In "Free" mode, transcript/gene fasta files from other sources (such as OrthoFinder, Annevo, etc.) can be supplied to easel. \
In this case, only fasta files are necessary, which can be downloaded from the [test data directory](https://github.com/casparbein/easel/tree). \
This data comes from [Enard et al. (2016)](https://elifesciences.org/articles/12469), accessible [here](https://datadryad.org/dataset/doi:10.5061/dryad.fs756).

```bash
## Make a directory for where to store the fasta files
mkdir easel_test_run_free
cd easel_test_run_free

## Download fasta files 
wget https://raw.githubusercontent.com/casparbein/easel/refs/heads/main/test_input/free/ENSG00000262505.fa
wget https://raw.githubusercontent.com/casparbein/easel/refs/heads/main/test_input/free/ENSG00000263042.fa
wget https://raw.githubusercontent.com/casparbein/easel/refs/heads/main/test_input/free/ENSG00000266412.fa
cd ..

## Run Test run
## Aligner: Prank (Nucleotide); Tree: Compute Gene Trees; Selection Screen: BUSTED+SRV
easel -free easel_test_run_free/ \
-ct \
-a prank_codon \
-bu srv \
-rs
```

## Command Line Interface
Exactly one input mode is required (`-toga2`, `-free`). Nothing runs unless `-dr` or `-rs` is given.

### Input (exactly one mode required)

| Flag | Value | Default | Description |
|---|---|---|---|
| `-toga2`, `--toga2_reference_path` | `DIR` | — | Directory holding the TOGA2 annotation runs to extract alignments from. The reference name is taken from the component above it, so a path ending in &lt;reference&gt;/TOGA2 is expected. Requires `--twoBit_path`. |
| `-free`, `--free_mode` | `DIR` | — | Directory of your own per-transcript FASTA files, one file per transcript. Combine with `--align` to align them, or with `--do_screen_only` if they are already aligned. |
| `-toga`, `--toga_reference_path` | `DIR` | — | TOGA v1 annotations. NOT YET SUPPORTED in this release: the alignment extraction rules for TOGA v1 are not part of the repository. Use `-toga2`, or pass already-extracted alignments with `-free`. |
| `-asm`, `--assemblies` | `FILE` | — | Single-column list of assemblies to include. Required with `-toga2`. Optional in `-free` mode, where leaving it out keeps every sequence found. Entries may be written either as 'name' or 'vs_name'. |
| `-sb`, `--selected_bed_file` | `BED` | — | BED12 file of the transcripts to screen; column 4 supplies the transcript names. Required with `-toga2`. Not used in `-free` mode, where every FASTA file in the input directory is run. |
| `--twoBit_path` | `FILE` | — | 2bit file of the reference exons, passed to TOGA2 alignment extraction. Required with `-toga2`. |

### Alignment

| Flag | Value | Default | Description |
|---|---|---|---|
| `-a`, `--align` | `NAME` | — | Aligner to use (case sensitive). Leave out only if the input is already aligned, and then set `--do_screen_only`. Options: <br>**prank** — TOGA2 input; nucleotide mode, codonified afterwards; <br> **prank_nt** — free input, nucleotide mode, codonified afterwards; <br> **prank_codon** — free input, codon mode; <br> **macse2** —  free input, codon-aware; muscle — free input, nucleotide mode, codonified afterwards. |

### Trees (choose one)

| Flag | Value | Default | Description |
|---|---|---|---|
| `-it`, `--input_tree` | `FILE` | — | Newick species tree, or bare topology, used for PRANK and for the selection screens. It is pruned per transcript to the species actually present in that alignment. |
| `-ct`, `--comp_tree` | — | flag | Compute a gene tree per transcript with IQ-TREE 3 instead of using a fixed species tree. |
| `--input_gene_trees` | `DIR` | — | Directory of precomputed gene trees, one per transcript, named &lt;transcript&gt;&lt;ext&gt;. Not supported together with `-toga2`. |
| `-st`, `--comp_species_tree` | — | flag | Reconstruct a species tree from the gene trees with ASTRAL-III. NOT YET IMPLEMENTED - easel exits if this is set. |

### Selection analyses (all off unless enabled)

| Flag | Value | Default | Description |
|---|---|---|---|
| `-ab`, `--absrel` | `MODES` | — | Run aBSREL. Comma-separated list of: std, srv, mh, all. 'std' means neither srv nor mh, and cannot be combined with them. |
| `-bu`, `--busted` | `MODES` | — | Run BUSTED. Comma-separated list of: std, srv, mh, error_sink, all. error_sink enables BUSTED-E and needs HyPhy &gt;= 2.5.58. |
| `-me`, `--meme` | `MODES` | — | Run MEME. Comma-separated list of: std, srv, mh, all. |
| `-re`, `--relax` | `MODES` | — | Run RELAX. Comma-separated list of: std, srv, mh, all, plus one integer giving the number of replicate runs to average over RELAX's stochasticity (default 10). Requires `--foreground_list`. |
| `-bc`, `--bayescode` | — | flag | Run BayesCode (mutation-selection per-site omega estimation). |
| `--foreground_list` | `FILE` | — | Single-column list of foreground species, which are labelled {Foreground} in the tree. Required for RELAX. Must be a subset of `--assemblies`. |

### Run control

| Flag | Value | Default | Description |
|---|---|---|---|
| `--directory_name` | `NAME` | `snakemake_selection_screen` | Name of the working directory created for this screen. |
| `-dr`, `--dry_run` | — | flag | Write DEF.yaml, then run snakemake `--dry-run` to show what would be executed. Takes precedence over `-rs`. |
| `-rs`, `--run_snakemake` | — | flag | Write DEF.yaml and launch the pipeline. Without `-dr` or `-rs`, easel only writes DEF.yaml and stops. |
| `-f`, `--force_run` | — | flag | Overwrite an existing DEF.yaml whose settings differ from the current command line. Without `-f`, easel refuses to touch the directory. |
| `--local` | — | flag | Run on this machine instead of submitting to SLURM. Without this, the bundled prof/config.yaml profile is used, which requires a cluster. |
| `--cores` | `N` | `4` | Cores to use with `--local`. |
| `--profile` | `DIR` | — | Snakemake profile directory to use instead of the bundled prof/. Point this at your own copy to set the SLURM partition and account. |
| `--conda_prefix` | `DIR` | — | Where to create the per-rule conda environments. Defaults to .snakemake/conda inside easel's own installation directory, shared across every run directory, so environments are built once instead of per run. |
| `--apptainer_prefix` | `DIR` | — | Where to cache pulled apptainer/singularity images (.sif files). Defaults to .snakemake/apptainer inside easel's own installation directory, shared across every run directory, so images are pulled once instead of per run. |
| `--rerun_triggers_mtime` | — | flag | Pass '`--rerun-triggers` mtime' to snakemake, so that a changed rule does not force everything downstream to rerun. |

### Advanced

| Flag | Value | Default | Description |
|---|---|---|---|
| `--do_alignment_only` | — | flag | Stop after alignment creation; skip trees and all selection analyses. |
| `--do_screen_only` | — | flag | Start at alignment postprocessing. The input must already be aligned, so this only makes sense with `-free`. |
| `--reference_name` | `NAME` | — | Name of the reference sequence used for codonification. Required in `-free` mode with a codonifying aligner (prank, prank_nt, muscle). |
| `--toga2_activate` | `FILE` | — | Path to an activate script that puts the toga2 executable on PATH. TOGA2 is not packaged on bioconda, so it is supplied through a container in easel. CURRENTLY NOT IMPLEMENTED. |
| `--extract_ali_params` | `LIST` | `skip_dups` | Comma-separated flags passed to TOGA1 alignment extraction: Options: skip_dups — drop non-one2one orthologs; allow_one2zero — keep lost transcripts; align_entirely — whole-gene instead of exon-by-exon; exclude_UL — drop 'uncertain loss' sequences. CURRENTLY NOT IMPLEMENTED. |
| `--max_CDS_length` | `BP` | `15000` | Skip transcripts whose CDS is longer than this. |
| `--min_CDS_length` | `BP` | `50` | Skip transcripts whose CDS is shorter than this. |
| `--min_num_aligned_species` | `N` | `5` | Skip transcripts present in fewer than N species. HyPhy needs at least 3 taxa, IQ-TREE bootstrapping at least 4. |
| `--do_hmm_cleaning` | — | flag | Clean alignments with HmmCleaner.pl. HmmCleaner is not on bioconda, so it runs from the ghcr.io/hillerlab/hmmcleaner container image ; this needs apptainer/singularity to be available wherever the pipeline actually runs. |
| `--hmm_cleaner_params` | `LIST` | `0.15,0.08,0.15,0.45` | Four HmmCleaner cost values c1,c2,c3,c4. The first two are negated, and they must increase: c1 &lt; c2 &lt; 0 &lt; c3 &lt; c4. |
| `--do_manual_cleaning` | — | flag | Clean alignments with the built-in column and row filter. |
| `--manual_cleaner_params` | `LIST` | `mc0.6,ms0.3,ml25,m` | Comma-separated filter settings: Options: mc&lt;f&gt; — min fraction of a column that must align for it to stay; ms&lt;f&gt; — min fraction of a sequence that must align for it to stay; ml&lt;n&gt; — min sequence length in amino acids; m — mask ambiguous and stop codons with NNN. |
| `--max_failed_fraction` | `F` | `0.4` | Stop before the analyses if more than this fraction of transcripts produced no alignment at all (phase 1 failures, as opposed to alignments that were produced and rejected by validation). |
| `--do_error_sink_cleaning` | — | flag | Filter alignments by the empirical Bayes factors from BUSTED-E. Implies `-bu` error_sink. |

### Resources (per rule)

| Flag | Value | Default | Description |
|---|---|---|---|
| `--extract_alignments_threads` | `N` | `1` | Threads for alignment extraction. |
| `--extract_alignments_mem_mb` | `MB` | `10000` | Memory in MB for alignment extraction. |
| `--hmmcleaner_threads` | `N` | `1` | Threads for HmmCleaner. |
| `--hmmcleaner_mem_mb` | `MB` | `1000` | Memory in MB for HmmCleaner. |
| `--manualcleaner_threads` | `N` | `1` | Threads for manual cleaner. |
| `--manualcleaner_mem_mb` | `MB` | `1000` | Memory in MB for manual cleaner. |
| `--tree_comp_threads` | `N` | `4` | Threads for tree computation. |
| `--tree_comp_mem_mb` | `MB` | `15000` | Memory in MB for tree computation. |
| `--prank_threads` | `N` | `5` | Threads for PRANK. |
| `--prank_mem_mb` | `MB` | `5000` | Memory in MB for PRANK. |
| `--absrel_threads` | `N` | `5` | Threads for aBSREL. |
| `--absrel_mem_mb` | `MB` | `10000` | Memory in MB for aBSREL. |
| `--busted_threads` | `N` | `5` | Threads for BUSTED. |
| `--busted_mem_mb` | `MB` | `10000` | Memory in MB for BUSTED. |
| `--meme_threads` | `N` | `5` | Threads for MEME. |
| `--meme_mem_mb` | `MB` | `10000` | Memory in MB for MEME. |
| `--relax_threads` | `N` | `10` | Threads for RELAX. |
| `--relax_mem_mb` | `MB` | `10000` | Memory in MB for RELAX. |



## Tools used by easel

Easel automates alignment, gene-tree inference and selection screen using the latest state-of-the-art tools. If you use easel, please cite these tools accordingly.
A list of all relevant articles is available [here](https://github.com/casparbein/easel/docs/citations.md).

Orthology inference:
- [TOGA2](https://github.com/hillerlab/TOGA2)

Alignment:
- [Prank](https://github.com/ariloytynoja/prank-msa) (Codon and Nucleotide modes)
- [Macse2](https://www.agap-ge2pop.org/macsee-pipelines/) (Codon, Free mode only)
- [Muscle5](https://github.com/rcedgar/muscle) (Nucleotide, Free mode only)

Gene tree inference:
- [IQtree3](https://github.com/iqtree/iqtree3)

Selection screens:
- [HyPhy](https://github.com/veg/hyphy)
  - [BUSTED](https://hyphy.org/methods/busted/)
  - [aBSREL](https://hyphy.org/methods/absrel/)
  - [MEME](https://hyphy.org/methods/meme/)
  - [RELAX](https://hyphy.org/methods/relax/)
- [BayesCode](https://github.com/ThibaultLatrille/bayescode)

## Known issues
- Depending on the number of input sequences, DAG creation can take up to ~20-30 Minutes (>5000 input sequences), mostly for dry-runs 

