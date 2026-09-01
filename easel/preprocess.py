#!/usr/bin/env python3
"""
preprocess_input.py

Pre-filters transcripts / genes before launching a easel selection screen.

Can be used:
  (a) as a standalone CLI tool, or
  (b) imported and called programmatically from cli.py via
      run_preprocessing().

TOGA / TOGA2 mode
-----------------
For every entry in the BED file the script:
  1. Looks up  <fasta_path>/<assembly>/loss_summary.tsv  for provided assemblies in
  <assembly_lst>
  2. In each loss_summary.tsv it searches for a row whose first column is
     exactly "TRANSCRIPT" and whose second column matches the BED entry name.
     The third column of that row is then reported (one line per assembly).
  3. Counts the total number of assemblies where the transcript was found as
     "TRANSCRIPT" across the whole dataset.
  4. Derives the CDS length of the transcript from the BED12 block sizes
     (column 10: comma-separated list of exon lengths).
  5. Flags and logs any transcript that fails one of the three quality filters:
       - fewer than <min_species> assemblies found    (default: 5)
       - CDS length > <max_cds>  base pairs           (default: 15 000)
       - CDS length < <min_cds>  base pairs           (default: 50)

Free mode
---------
For every FASTA file found directly inside the free path the script:
  1. Counts the number of sequences (entries) in the file.
  2. Finds the maximum sequence length across all entries.
  3. Applies the same three filters as above, with "max sequence length"
     in place of "CDS length".

Output
------
  - A log file (--output_log, default: preprocess_report.log) summarising
    all findings.
  - A tab-separated report table (<output_log stem>_table.tsv).
  - A plain-text exclusion list (--output_exclusion, default:
    excluded_transcripts.txt) with the names of transcripts / FASTA files
    that failed at least one filter.
"""

import argparse
import hashlib
import logging
import os
import sys

## ---------------------------------------------------------------------------
## Logging
## ---------------------------------------------------------------------------
logger = logging.getLogger("preprocess")


## ---------------------------------------------------------------------------
## Argument parsing (standalone CLI use)
## ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    app = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    mode_group = app.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--fasta_path",
        dest="fasta_path",
        metavar="PATH",
        help=(
            "Root directory of TOGA / TOGA2 / Free annotation runs.\n"
            "Each sub-directory is expected to be one assembly and to contain\n"
            "a loss_summary.tsv file."
        ),
    )
    # mode_group.add_argument(
    #     "--free_path",
    #     dest="free_path",
    #     metavar="PATH",
    #     help=(
    #         "Directory containing user-supplied FASTA files (free mode).\n"
    #         "Each .fa / .fasta / .fna file is treated as one gene/locus."
    #     ),
    # )

    app.add_argument(
        "--bed_file",
        dest="bed_file",
        metavar="FILE",
        default=None,
        help=(
            "BED12 file of transcripts to consider (TOGA / TOGA2 mode).\n"
            "Can be omitted for free mode"
        ),
    )

    app.add_argument(
        "--assembly_list",
        dest="assembly_list",
        metavar="ASM",
        default=None,
        help=(
            "Single column list of assemblies to be included in screen (TOGA / TOGA2 mode)"
        ),
    )

    app.add_argument(
        "--min_species",
        dest="min_species",
        type=int,
        default=5,
        metavar="N",
        help="Minimum number of assemblies in which a transcript must be\n"
             "found with status TRANSCRIPT. Default: 5",
    )
    app.add_argument(
        "--max_cds",
        dest="max_cds",
        type=int,
        default=15000,
        metavar="BP",
        help="Maximum CDS / sequence length in base pairs. Default: 15000",
    )
    app.add_argument(
        "--min_cds",
        dest="min_cds",
        type=int,
        default=50,
        metavar="BP",
        help="Minimum CDS / sequence length in base pairs. Default: 50",
    )
    app.add_argument(
        "--output_log",
        dest="output_log",
        default="preprocess_report.log",
        metavar="FILE",
        help="Path for the detailed report log. Default: preprocess_report.log",
    )
    app.add_argument(
        "--output_exclusion",
        dest="output_exclusion",
        default="excluded_transcripts.txt",
        metavar="FILE",
        help=(
            "Path for the exclusion list (one entry per line).\n"
            "Default: excluded_transcripts.txt"
        ),
    )
    app.add_argument(
        "--filtered_out_bed",
        dest="output_bed",
        default="filtered.bed",
        metavar="OBED",
        help=(
            "Path for filtered bed file (TOGA/TOGA2 mode).\n"
            "Default: filtered.bed"
        ),
    )
    return app


## ---------------------------------------------------------------------------
## BED12 helpers
## ---------------------------------------------------------------------------
def parse_bed12(bed_path: str) -> list[dict]:
    """Parse a BED12 file and return a list of transcript dicts.

    Each dict contains:
      name     : transcript / gene name (column 3)
      cds_len  : total CDS length in bp (sum of block sizes, column 10)
    """
    transcripts = []
    with open(bed_path, "r") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                logger.warning(
                    "BED line %d has fewer than 12 fields – skipped: %s",
                    lineno, line[:80],
                )
                continue
            name = fields[3]
            try:
                ## from Bogdan's old script
                chrom_start = int(fields[1])
                thick_start = int(fields[6])
                thick_end = int(fields[7])
                block_count = int(fields[9])
                block_sizes = [int(x) for x in fields[10].split(",") if x]
                block_starts = [int(x) for x in fields[11].split(",") if x]
                block_ends = [block_starts[i] + block_sizes[i] for i in range(block_count)]
                block_abs_starts = [block_starts[i] + chrom_start for i in range(block_count)]
                block_abs_ends = [block_ends[i] + chrom_start for i in range(block_count)]
                block_new_starts, block_new_ends = [], []

                for block_num in range(block_count):
                    # go block-by-block
                    block_start = block_abs_starts[block_num]
                    block_end = block_abs_ends[block_num]

                    # skip the block if it is entirely UTR
                    if block_end <= thick_start:
                        continue
                    elif block_start >= thick_end:
                        continue

                    block_new_start = block_start if block_start >= thick_start else thick_start
                    block_new_end = block_end if block_end <= thick_end else thick_end
                    block_new_starts.append(block_new_start - thick_start)
                    block_new_ends.append(block_new_end - thick_start)

                cds_len = sum(block_new_ends[i] - block_new_starts[i] for i in range(len(block_new_starts)))

            except ValueError:
                logger.warning(
                    "Could not parse block sizes for %s (line %d) – "
                    "falling back to end-start length.",
                    name, lineno,
                )
                cds_len = int(fields[2]) - int(fields[1])

            transcripts.append({"name": name, "cds_len": cds_len})

    logger.info("Parsed %d transcript entries from %s", len(transcripts), bed_path)
    return transcripts



## ---------------------------------------------------------------------------
## TOGA / TOGA2 helpers
## ---------------------------------------------------------------------------
LOSS_SUMMARY_FILENAME = "loss_summary.tsv"
ORTHOLOGY_FILENAME = "orthology_classification.tsv"
TRANSCRIPT_STATUS = "TRANSCRIPT"

def _get_assembly_list(path_to_assembly):
    """Assembly names as vs_<name> run-directory names.

    Delegates the parsing to formats.read_one_column so blanks, comments and
    stray whitespace are handled the same way here as in the CLI.
    """
    from .formats import read_one_column
    return ["vs_" + name for name in
            read_one_column(path_to_assembly, label="assembly list", strip_vs=True)]


def read_list(path_to_list):
    """One-column list, via the shared reader."""
    from .formats import read_one_column
    return read_one_column(path_to_list, label="list")


def _iter_assembly_dirs(fasta_path: str, assembly_list: list) -> list[str]:
    """Return sorted list of first-level subdirectory names under fasta_path."""
    dirs = [a for a in assembly_list if os.path.isdir(os.path.join(fasta_path, a))]
    logger.info("Found %d assembly directories under %s", len(dirs), fasta_path)
    return dirs


def _load_loss_summary(tsv_path: str) -> dict[str, str]:
    """Parse loss_summary.tsv and return {transcript_name: col2_value}
    for all rows where column 0 == 'TRANSCRIPT'.
    """
    hits: dict[str, str] = {}
    if not os.path.isfile(tsv_path):
        return hits
    with open(tsv_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            status, name, value = fields[0], fields[1], fields[2]
            if status == TRANSCRIPT_STATUS:
                hits[name] = value
    return hits

def _load_orthology_classification(tsv_path: str) -> dict[str, str]:
    """Parse orthology_classification.tsv and return {transcript_name: orthology}
    for all rows where column 2 == transcript_name.
    """
    orthos: dict[str, str] = {}
    if not os.path.isfile(tsv_path):
        return orthos
    with open(tsv_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("t_gene"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            transcript_id, orthology_class = fields[1], fields[4]
            orthos[transcript_id] = orthology_class
    return orthos


def run_toga_mode(
    fasta_path: str,
    bed_path: str | None,
    foreground_list: str | None,
    assembly_path : str,
    min_species: int,
    max_cds: int,
    min_cds: int,
) -> tuple[list[str], list[dict]]:
    """Main logic for TOGA / TOGA2 mode.

    Returns
    -------
    excluded : list[str]
        Transcript names that failed at least one filter.
    report_rows : list[dict]
        One dict per transcript with full details for the report.
    """
    assembly_list  = _get_assembly_list(assembly_path)
    assembly_dirs = _iter_assembly_dirs(fasta_path, assembly_list)
    if not assembly_dirs:
        logger.critical("No assembly subdirectories found under %s", fasta_path)
        sys.exit(1)

    ## --- Determine transcript set ---
    if bed_path:
        transcripts = parse_bed12(bed_path)
        print(len(transcripts))
    else:
        logger.critical(
            "No bed_file provided. Please provide a bed file if TOGA transcripts should be run."
            )
        sys.exit(1)

    if foreground_list:
        foreground_list = _get_assembly_list(foreground_list)
        logger.info("Foreground list was provided, any transcript has to at least contain one foreground assembly")

    ## --- Pre-load all loss_summary.tsv files ---
    logger.info("Loading loss_summary.tsv and orthology classification files from %d assemblies…", len(assembly_dirs))
    asm_hits: dict[str, dict[str, str]] = {}
    asm_orthos: dict[str, dict[str, str]] = {}
    missing_tsv: list[str] = []
    missing_ortho: list[str] = []
    for asm in assembly_dirs:

        loss_path = os.path.join(fasta_path, asm, LOSS_SUMMARY_FILENAME)
        ortho_path = os.path.join(fasta_path, asm, ORTHOLOGY_FILENAME)
        #hits = _load_loss_summary(loss_path)
        #orthos = _load_orthology_classification(ortho_path)
        
        if not os.path.isfile(loss_path):
            missing_tsv.append(asm)
        elif not os.path.isfile(ortho_path):
            missing_ortho.append(asm)
        else:
            asm_hits[asm] = _load_loss_summary(loss_path)
            asm_orthos[asm] = _load_orthology_classification(ortho_path)

    if missing_tsv or missing_ortho:
        logger.warning(
            "%d/%d assembly directories have no loss_summary.tsv and/or no orthology classification file: %s%s",
            len(missing_tsv),
            len(missing_ortho),
            ", ".join(missing_tsv[:5]),
            " …" if len(missing_tsv) > 5 else "",
        )

    ## --- Evaluate each transcript ---
    excluded: list[str] = []
    report_rows: list[dict] = []
    #mask_for_orthology: list[str] = []

    for tx in transcripts:
        name = tx["name"]
        cds_len = tx["cds_len"]

        found_in_loss: dict[str, str] = {}
        found_in_ortho: dict[str, str] = {}
        for asm in asm_hits:
            status = asm_hits[asm].get(name)
            if status is None:
                continue
            if status in ("FI", "I", "PI"):
                found_in_loss[asm] = status
            else:
                ## UL / M / L / N / PG: no intact copy, so the orthology call
                ## for this assembly is not trusted.
                continue

            orthology = asm_orthos.get(asm, {}).get(name)
            if orthology in ("many2one", "one2one"):
                found_in_ortho[asm] = orthology

        species_count = len(set(found_in_loss.keys()).intersection(set(found_in_ortho.keys())))
        species_data = set(found_in_loss.keys()).intersection(set(found_in_ortho.keys()))

        failures: list[str] = []
        if species_count < min_species:
            failures.append(
                f"too few assemblies with Intact copies and/or right orthology: {species_count} < {min_species}"
            )
        if cds_len is not None and cds_len > max_cds:
            failures.append(f"CDS too long: {cds_len} bp > {max_cds} bp")
        if cds_len is not None and cds_len < min_cds:
            failures.append(f"CDS too short: {cds_len} bp < {min_cds} bp")
        if foreground_list and not any(species in species_data for species in foreground_list):
            failures.append("No foreground species were found in alignment")

        row = {
            "name": name,
            "cds_len": cds_len,
            "species_count": species_count,
            "found_in": set(found_in_loss.keys()).intersection(set(found_in_ortho.keys())),
            "failures": failures,
        }
        report_rows.append(row)

        if failures:
            logger.warning(
                "EXCLUDE %s  |  species=%d  |  cds=%s bp  |  reason: %s",
                name,
                species_count,
                str(cds_len) if cds_len is not None else "n/a",
                "; ".join(failures),
            )
            excluded.append(name)
        else:
            logger.info(
                "OK       %s  |  species=%d  |  cds=%s bp",
                name,
                species_count,
                str(cds_len) if cds_len is not None else "n/a",
            )

    return excluded, report_rows


## ---------------------------------------------------------------------------
## Free-mode helpers
## ---------------------------------------------------------------------------
FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".fas"}


def _iter_fasta_files(free_path: str) -> list[str]:
    """Return sorted list of FASTA file paths directly inside free_path."""
    files = []
    try:
        for entry in sorted(os.listdir(free_path)):
            full = os.path.join(free_path, entry)
            if os.path.isfile(full) and os.path.splitext(entry)[1].lower() in FASTA_EXTENSIONS:
                files.append(full)
    except OSError as exc:
        logger.critical("Cannot list free_path directory: %s", exc)
        sys.exit(1)
    logger.info("Found %d FASTA files under %s", len(files), free_path)
    return files


def _read_fasta_stats(fasta_path: str) -> tuple[int, int]:
    """Count sequences and find the maximum sequence length in a FASTA file.

    Returns (n_sequences, max_length).  Alignment gaps ('-') are excluded from
    length counts.
    """
    n_seqs = 0
    max_len = 0
    current_len = 0
    seq_names = []
    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                seq_names.append(line)
                if n_seqs > 0:
                    max_len = max(max_len, current_len)
                n_seqs += 1
                current_len = 0
            else:
                current_len += len(line.replace("-", "").replace("N", ""))
        if n_seqs > 0:
            max_len = max(max_len, current_len)
    return n_seqs, max_len, seq_names


## To-Do: enable removal of certain assemblies/queries in fasta file
def run_free_mode(
    fasta_path: str,
    min_species: int,
    foreground_list: str | None,
    max_cds: int,
    min_cds: int,
    assembly_path:  str | None = None,
    precomp_gene_tree_path: str | None = None,
) -> tuple[list[str], list[dict]]:
    """Main logic for free mode.

    Returns
    -------
    excluded : list[str]
        FASTA file basenames that failed at least one filter.
    report_rows : list[dict]
    """
    fasta_files = _iter_fasta_files(fasta_path)
    if not fasta_files:
        logger.critical("No FASTA files found under %s", fasta_path)
        sys.exit(1)

    excluded: list[str] = []
    report_rows: list[dict] = []

    if foreground_list:
        foreground_list = read_list(foreground_list)
        logger.info("Foreground list was provided, any transcript has to at least contain one foreground assembly")

    for fpath in fasta_files:
        name = os.path.basename(fpath)
        n_seqs, max_len, seq_names = _read_fasta_stats(fpath)

        failures: list[str] = []
        if n_seqs < min_species:
            failures.append(f"too few sequences (species): {n_seqs} < {min_species}")
        if max_len > max_cds:
            failures.append(f"max sequence length too long: {max_len} bp > {max_cds} bp")
        if max_len < min_cds:
            failures.append(f"max sequence length too short: {max_len} bp < {min_cds} bp")
        ## Was inverted (it flagged a failure when foreground species WERE
        ## present) and could never match anyway, because seq_names holds raw
        ## header lines including the leading '>'.
        if foreground_list:
            ids = {n.lstrip(">").split()[0] if n.lstrip(">").split() else ""
                   for n in seq_names}
            if not (ids & set(foreground_list)):
                failures.append("No foreground species were found in alignment")
        
        row = {
            "name": name,
            "n_sequences": n_seqs,
            "max_length": max_len,
            "failures": failures,
        }
        report_rows.append(row)

        if failures:
            logger.warning(
                "EXCLUDE %s  |  n_seqs=%d  |  max_len=%d bp  |  reason: %s",
                name, n_seqs, max_len, "; ".join(failures),
            )
            excluded.append(name)
        else:
            logger.info(
                "OK       %s  |  n_seqs=%d  |  max_len=%d bp",
                name, n_seqs, max_len,
            )

    return excluded, report_rows


## ---------------------------------------------------------------------------
## Report writing
## ---------------------------------------------------------------------------
def write_report(report_rows: list[dict], output_path: str, mode: str) -> None:
    """Write a tab-separated report of all evaluated entries to output_path."""
    with open(output_path, "w") as fh:
        if mode == "toga":
            fh.write("transcript\tcds_len_bp\tspecies_count\tstatus\treasons\tfound_in\n")
            for row in report_rows:
                status = "EXCLUDED" if row["failures"] else "PASS"
                fh.write(
                    "\t".join([
                        row["name"],
                        str(row["cds_len"]) if row["cds_len"] is not None else "n/a",
                        str(row["species_count"]),
                        status,
                        "; ".join(row["failures"]) if row["failures"] else "-",
                        ",".join(sorted(row["found_in"])),
                    ]) + "\n"
                )
        else:
            fh.write("fasta_file\tn_sequences\tmax_len_bp\tstatus\treasons\n")
            for row in report_rows:
                status = "EXCLUDED" if row["failures"] else "PASS"
                fh.write(
                    "\t".join([
                        row["name"],
                        str(row["n_sequences"]),
                        str(row["max_length"]),
                        status,
                        "; ".join(row["failures"]) if row["failures"] else "-",
                    ]) + "\n"
                )
    logger.info("Full report written to %s", output_path)


def write_exclusion_list(excluded: list[str], output_path: str) -> None:
    """Write the exclusion list (one entry per line)."""
    with open(output_path, "w") as fh:
        for name in excluded:
            fh.write(name + "\n")
    logger.info(
        "Exclusion list with %d entr%s written to %s",
        len(excluded),
        "ies" if len(excluded) != 1 else "y",
        output_path,
    )

def write_filtered_bed(excluded: list[str], input_bed: str, output_path: str) -> None:
    out_bed_list = []
    with open(input_bed, "r") as ib:
        for line in ib:
            fields = line.rstrip().split('\t')
            if fields[3] in excluded:
                continue
            else:
                out_bed_list.append(fields)
    
    with open(output_path, "w") as ob:
        for entry in out_bed_list:
            ob.write("\t".join(entry) + "\n")
    
    logger.info(
        "Filtered bed file written to %s",
        output_path,
    )


## ---------------------------------------------------------------------------
## Staleness fingerprinting
## ---------------------------------------------------------------------------
## Lets a caller (easel.cli) detect that an input's *content* changed even
## though its path didn't -- a path/parameter comparison alone would miss an
## edit-in-place, e.g. someone adding a species to assemblies.txt without
## renaming it.

def file_fingerprint(path: str | None) -> str | None:
    """sha256 of a file's bytes, or None if path is falsy/missing.

    Cheap enough for the small single-column/BED files (assembly list,
    foreground list, bed file). See fingerprint_fasta_path for fasta_path
    itself, where hashing content would be nearly as expensive as
    preprocessing.
    """
    if not path or not os.path.isfile(path):
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint_fasta_path(
    fasta_path: str, toga_mode: bool, assembly_file: str | None = None
) -> list:
    """[relpath, mtime_ns, size] for exactly the files preprocessing would
    read under fasta_path: the per-assembly loss_summary.tsv /
    orthology_classification.tsv in TOGA/TOGA2 mode, or the FASTA files
    themselves in free mode.

    mtime+size, not a content hash. This gets written into DEF.yaml and compared against
    the copy reloaded from it on the next run to decide whether preprocessing
    is still up to date.
    """
    if toga_mode:
        assemblies = _get_assembly_list(assembly_file) if assembly_file else []
        dirs = _iter_assembly_dirs(fasta_path, assemblies)
        files = [os.path.join(fasta_path, d, name)
                 for d in dirs for name in (LOSS_SUMMARY_FILENAME, ORTHOLOGY_FILENAME)]
    else:
        files = _iter_fasta_files(fasta_path)

    stats = []
    for f in files:
        try:
            st = os.stat(f)
            stats.append([os.path.relpath(f, fasta_path), st.st_mtime_ns, st.st_size])
        except OSError:
            stats.append([os.path.relpath(f, fasta_path), None, None])
    return sorted(stats)


## ---------------------------------------------------------------------------
## Public API — called programmatically from cli.py
## ---------------------------------------------------------------------------
def run_preprocessing(
    TOGA_mode: bool = True,
    free_mode: bool = False, 
    fasta_path: str | None = None,
    foreground_list: str | None = None,
    bed_file: str | None = None,
    assembly_file: str | None = None,
    min_species: int = 5,
    max_cds: int = 15000,
    min_cds: int = 50,
    output_log: str = "preprocess_report.log",
    output_exclusion: str = "excluded_transcripts.txt",
    output_filtered_bed: str = "filtered.bed",
) -> list[str]:
    """Run the preprocessing pipeline and return the excluded transcript list.

    Parameters
    ----------
    TOGA_mode:
    free_mode:
    bed_file : str, optional
        BED12 file of transcripts to consider (TOGA / TOGA2 mode only).
    assembly_file : str, optional
        Single-column assembly list of assemblies to be searched (TOGA / TOGA2 mandatory, free mode optional).
    min_species : int
        Minimum number of assemblies required per transcript. Default 5.
    max_cds : int
        Maximum CDS / sequence length in bp. Default 15 000.
    min_cds : int
        Minimum CDS / sequence length in bp. Default 50.
    output_log : str
        Path for the log / report file. Default 'preprocess_report.log'.
    output_exclusion : str
        Path for the exclusion list. Default 'excluded_transcripts.txt'.

    Returns
    -------
    list[str]
        Names of transcripts / FASTA files that failed at least one filter.
        Empty list when everything passes.
    """
    if not TOGA_mode and not free_mode:
        logger.critical(
            "run_preprocessing() requires either TOGA input or free input to be set."
        )
        sys.exit(1)
    if TOGA_mode and free_mode:
        logger.critical(
            "run_preprocessing() received both TOGA mode and free mode directive — "
            "only one mode can be active at a time."
        )
        sys.exit(1)

    ## Add a file handler so the log is also persisted to disk
    file_handler = logging.FileHandler(output_log, mode="w")
    file_handler.setFormatter(logging.Formatter("[%(levelname)-8s] %(message)s"))
    logger.addHandler(file_handler)

    logger.info("preprocess_input started")
    logger.info(
        "Filters: min_species=%d  max_cds=%d bp  min_cds=%d bp",
        min_species, max_cds, min_cds,
    )

    if TOGA_mode:
        logger.info("Mode: TOGA / TOGA2  (fasta_path: %s)", fasta_path)
        excluded, report_rows = run_toga_mode(
            fasta_path=fasta_path,
            bed_path=bed_file,
            foreground_list=foreground_list,
            assembly_path=assembly_file,
            min_species=min_species,
            max_cds=max_cds,
            min_cds=min_cds,
        )
        mode = "toga"
    else:
        logger.info("Mode: free  (free_path: %s)", fasta_path)
        excluded, report_rows = run_free_mode(
            fasta_path=fasta_path,
            foreground_list=foreground_list,
            min_species=min_species,
            assembly_path=assembly_file,
            max_cds=max_cds,
            min_cds=min_cds,
        )
        mode = "free"

    total = len(report_rows)
    n_excluded = len(excluded)
    logger.info(
        "Summary: %d entries evaluated — %d PASS, %d EXCLUDED",
        total, total - n_excluded, n_excluded,
    )
    if n_excluded:
        logger.warning("Excluded entries: %s", ", ".join(excluded))
    if total == n_excluded:
        logger.critical("All entries were removed from bed file. Adjust your filtering or check the input files")
        sys.exit(1)

    table_path = output_log.replace(".log", "_table.tsv")
    write_report(report_rows, table_path, mode)
    write_exclusion_list(excluded, output_exclusion)
    if mode == "toga":
        write_filtered_bed(excluded, bed_file, output_filtered_bed)

    ## Remove file handler so it does not bleed into callers
    logger.removeHandler(file_handler)
    file_handler.close()

    return excluded


## ---------------------------------------------------------------------------
## Standalone CLI entry point
## ---------------------------------------------------------------------------
def main() -> None:
    ## Ensure the module logger has a console handler when used standalone
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="[%(levelname)-8s] %(message)s")

    parser = build_parser()
    args = parser.parse_args()

    ## Basic path validation
    if args.fasta_path and not os.path.isdir(args.fasta_path):
        logger.critical(
            "--fasta_path does not exist or is not a directory: %s", args.fasta_path
        )
        sys.exit(1)
    if args.bed_file and not os.path.isfile(args.bed_file):
        logger.critical("--bed_file not found: %s", args.bed_file)
        sys.exit(1)
    # if args.free_path and not os.path.isdir(args.free_path):
    #     logger.critical(
    #         "--free_path does not exist or is not a directory: %s", args.free_path
    #     )
    #     sys.exit(1)
    # if args.free_path and args.bed_file:
    #     logger.warning(
    #         "--bed_file is ignored in free mode "
    #         "(sequences are read from the FASTA files directly)"
    #     )

    run_preprocessing(
        fasta_path=args.fasta_path,
        bed_file=args.bed_file,
        assembly_file=args.assembly_list,
        min_species=args.min_species,
        max_cds=args.max_cds,
        min_cds=args.min_cds,
        output_log=args.output_log,
        output_exclusion=args.output_exclusion,
        output_filtered_bed=args.output_bed,
    )


if __name__ == "__main__":
    main()
