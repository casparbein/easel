"""Validation of the biological input files themselves.

Functions a caller can reasonably recover from return ``(ok, message)``.
Functions whose failure means the run cannot proceed log CRITICAL and raise
:class:`InputError`; library code should not call sys.exit itself.
"""

import gzip
from collections import Counter
import logging
import os
import struct
from pathlib import Path

logger = logging.getLogger("easel.formats")

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".fas"}
TREE_EXTENSIONS = {".nwk", ".nh", ".newick", ".tree", ".tre", ".nhx", ".treefile"}

## UCSC 2bit signature 0x1A412743, stored in the writer's byte order, so a
## reader has to accept both.
TWOBIT_MAGIC = (0x1A412743, 0x43274119)


class InputError(Exception):
    """Raised when an input file cannot be used. Already logged when raised."""


def _tree_class():
    """Import ete4 only when a tree actually has to be parsed.
    """
    try:
        from ete4 import Tree
    except ImportError as exc:                       # pragma: no cover
        _fail("ete4 is required to read Newick trees but is not installed: %s", exc)
    return Tree


def _fail(message, *args):
    logger.critical(message, *args)
    raise InputError(message % args if args else message)


## ---------------------------------------------------------------------------
## One-column lists (assemblies, foreground species)
## ---------------------------------------------------------------------------

def read_one_column(path, label="list", strip_vs=False, unique=True):
    """Parse a one-column text file into a list of names.
    Blank lines and '#' comments are skipped, surrounding whitespace is
    stripped, and a line with more than one field is an error rather than
    something each caller splits differently.
    """
    names = []
    try:
        with open(path, "r", encoding="utf-8") as handle:
            for lineno, raw in enumerate(handle, 1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) > 1:
                    _fail("%s %s:%d has %d fields, expected one name per line: %r",
                          label, path, lineno, len(fields), line)
                name = fields[0]
                if strip_vs and name.startswith("vs_"):
                    name = name[len("vs_"):]
                names.append(name)
    except OSError as exc:
        _fail("cannot read %s %s: %s", label, path, exc)

    if not names:
        _fail("%s %s contains no entries", label, path)

    if unique:
        counts = Counter(names)
        duplicates = sorted(n for n, c in counts.items() if c > 1)
        if duplicates:
            shown = ", ".join(duplicates[:5]) + (" ..." if len(duplicates) > 5 else "")
            _fail("%s %s contains duplicate entries: %s", label, path, shown)

    logger.info("Read %d entr%s from %s", len(names),
                "y" if len(names) == 1 else "ies", path)
    return names


def check_toga_run_dirs(toga_path, assemblies):
    """Every assembly must have a vs_<name> run directory under *toga_path*.
    """
    missing = [a for a in assemblies
               if not os.path.isdir(os.path.join(toga_path, "vs_" + a))]
    if missing:
        shown = ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else "")
        _fail("%d of %d assemblies have no run directory under %s: %s",
              len(missing), len(assemblies), toga_path, shown)
    logger.info("All %d assemblies have a run directory under %s",
                len(assemblies), toga_path)
    return True


## ---------------------------------------------------------------------------
## Trees
## ---------------------------------------------------------------------------

def read_in_tree(tree_path, min_tips=3):
    """Parse a Newick file and return the ete4 Tree."""
    if not os.path.isfile(tree_path):
        _fail("tree file not found: %s", tree_path)
    try:
        with open(tree_path, "r", encoding="utf-8") as handle:
            tree = _tree_class()(handle, parser=1)
    except Exception as exc:
        _fail("%s is not readable as Newick (ete4 parser=1): %s", tree_path, exc)

    tips = [leaf.name for leaf in tree.leaves()]
    unnamed = sum(1 for name in tips if not name)
    if unnamed:
        _fail("%s has %d unnamed tip(s); every leaf needs a name to be matched "
              "against the alignments", tree_path, unnamed)
    duplicates = sorted(n for n, c in Counter(tips).items() if c > 1)
    if duplicates:
        shown = ", ".join(duplicates[:5]) + (" ..." if len(duplicates) > 5 else "")
        _fail("%s has duplicate tip labels: %s. HyPhy keys its branch attributes "
              "by name, so duplicates are ambiguous.", tree_path, shown)
    if len(tips) < min_tips:
        _fail("%s has %d tips; at least %d are needed",
              tree_path, len(tips), min_tips)

    logger.info("Tree %s is valid Newick with %d tips", tree_path, len(tips))
    return tree


def compare_tree_asm(in_tree, assembly_list, reference=None):
    """Check that the tree covers every assembly. Returns (ok, message).
    """
    expected = set(assembly_list)
    if reference:
        expected = expected | {reference}
    tips = {leaf.name for leaf in in_tree.leaves()}

    extra = sorted(tips - expected)
    if extra:
        logger.warning(
            "%d tip(s) in the input tree are not in the assembly list and will "
            "be pruned per transcript: %s%s",
            len(extra), ", ".join(extra[:5]), " ..." if len(extra) > 5 else "")

    absent = sorted(expected - tips)
    if absent:
        shown = ", ".join(absent[:5]) + (" ..." if len(absent) > 5 else "")
        return False, (f"{len(absent)} assembl{'y' if len(absent) == 1 else 'ies'} "
                       f"absent from the input tree: {shown}. Every assembly in "
                       f"the screen must be a tip in the tree.")
    return True, f"Input tree covers all {len(expected)} expected taxa"

## This extends tree format checks for precomputed gene trees
def validate_tree_directory(input_dir, extensions=None):
    """Every tree file in *input_dir* must parse. Returns (ok, message)."""
    extensions = extensions or TREE_EXTENSIONS
    files = [p for p in sorted(Path(input_dir).iterdir())
             if p.is_file() and not p.name.startswith(".")]
    if not files:
        return False, f"Directory '{input_dir}' contains no files."

    checked, skipped = 0, []
    for path in files:
        name = path.name[:-3] if path.name.endswith(".gz") else path.name
        if not any(name.endswith(ext) for ext in extensions):
            skipped.append(path.name)
            continue
        try:
            ## Maybe use read_in_tree and compare_tree_asm here (?)
            with open(path, "r", encoding="utf-8") as handle:
                _tree_class()(handle, parser=1)
        except Exception as exc:
            return False, f"{path} is not readable as Newick: {exc}"
        checked += 1

    if skipped:
        logger.warning("Ignored %d non-tree file(s) in %s: %s%s", len(skipped),
                       input_dir, ", ".join(skipped[:5]),
                       " ..." if len(skipped) > 5 else "")
    if not checked:
        return False, f"No tree files with a recognised extension in '{input_dir}'."
    return True, f"All {checked} tree files in '{input_dir}' parsed successfully."


## ---------------------------------------------------------------------------
## BED
## ---------------------------------------------------------------------------

def check_bed_file(bed_path, min_fields=12):
    """Validate a BED12 file. Returns (ok, message).
    """
    names, problems, n_records = [], [], 0
    try:
        with open(bed_path, "r", encoding="utf-8") as handle:
            for lineno, raw in enumerate(handle, 1):
                line = raw.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split("\t")
                n_records += 1
                if len(fields) < min_fields:
                    problems.append(f"line {lineno}: {len(fields)} fields, "
                                    f"expected at least {min_fields}")
                    continue
                name = fields[3]
                names.append(name)
                if "/" in name or "\\" in name or name != name.strip():
                    problems.append(f"line {lineno}: transcript name {name!r} "
                                    f"cannot be a directory component")
                try:
                    chrom_start, chrom_end = int(fields[1]), int(fields[2])
                    thick_start, thick_end = int(fields[6]), int(fields[7])
                    block_count = int(fields[9])
                except ValueError:
                    problems.append(f"line {lineno}: non-numeric coordinate field")
                    continue
                if chrom_start >= chrom_end:
                    problems.append(f"line {lineno}: chromStart >= chromEnd")
                if thick_start > thick_end:
                    problems.append(f"line {lineno}: thickStart > thickEnd")
                if fields[5] not in ("+", "-", "."):
                    problems.append(f"line {lineno}: strand is {fields[5]!r}")
                sizes = [x for x in fields[10].split(",") if x]
                starts = [x for x in fields[11].split(",") if x]
                if not (len(sizes) == len(starts) == block_count):
                    problems.append(
                        f"line {lineno}: blockCount is {block_count} but there "
                        f"are {len(sizes)} sizes and {len(starts)} starts")
                if len(problems) >= 20:
                    problems.append("... further problems not listed")
                    break
    except OSError as exc:
        return False, f"cannot read {bed_path}: {exc}"

    if not n_records:
        return False, f"{bed_path} contains no BED records."

    duplicates = sorted(n for n, c in Counter(names).items() if c > 1)
    if duplicates:
        shown = ", ".join(duplicates[:5]) + (" ..." if len(duplicates) > 5 else "")
        problems.insert(0, f"duplicate transcript name(s): {shown}. Each name "
                           f"becomes an output directory, there would be naming inconsistencies in snakemake.")
    if problems:
        return False, (f"{bed_path} is not a usable BED12 file:\n    "
                       + "\n    ".join(problems))

    logger.info("BED file %s validated: %d records", bed_path, n_records)
    return True, f"{bed_path} validated: {n_records} records"


## ---------------------------------------------------------------------------
## FASTA
## ---------------------------------------------------------------------------

def is_fasta(file_path):
    """True if the first non-blank line starts with '>'. Handles .gz."""
    try:
        if str(file_path).endswith(".gz"):
            handle = gzip.open(file_path, "rt", encoding="utf-8")
        else:
            handle = open(file_path, "r", encoding="utf-8")
        with handle:
            for line in handle:
                if line.strip():
                    return line.lstrip().startswith(">")
    except (OSError, UnicodeDecodeError, EOFError) as exc:
        logger.warning("Cannot read %s: %s", file_path, exc)
        return False
    return False


def validate_transcript_directory(input_dir, extensions=None):
    """Every FASTA in *input_dir* must look like FASTA. Returns (ok, message).
    Non-FASTA files are skipped with a warning rather than being fatal.
    """
    extensions = extensions or FASTA_EXTENSIONS
    files = [p for p in sorted(Path(input_dir).iterdir())
             if p.is_file() and not p.name.startswith(".")]
    if not files:
        return False, f"Directory '{input_dir}' contains no files."

    checked, skipped, bad = 0, [], []
    for path in files:
        name = path.name[:-3] if path.name.endswith(".gz") else path.name
        if not any(name.endswith(ext) for ext in extensions):
            skipped.append(path.name)
            continue
        if not is_fasta(path):
            bad.append(path.name)
            continue
        checked += 1

    if skipped:
        logger.warning("Ignored %d non-FASTA file(s) in %s: %s%s", len(skipped),
                       input_dir, ", ".join(skipped[:5]),
                       " ..." if len(skipped) > 5 else "")
    if bad:
        shown = ", ".join(bad[:5]) + (" ..." if len(bad) > 5 else "")
        return False, (f"{len(bad)} file(s) in '{input_dir}' have a FASTA "
                       f"extension but do not start with '>': {shown}")
    if not checked:
        return False, (f"No FASTA files with a recognised extension "
                       f"({', '.join(sorted(extensions))}) in '{input_dir}'.")
    return True, f"All {checked} FASTA files in '{input_dir}' validated successfully."


## ---------------------------------------------------------------------------
## 2bit
## ---------------------------------------------------------------------------

def check_twobit(path):
    """Verify a UCSC .2bit file by its magic number. Returns (ok, message).
    """
    try:
        with open(path, "rb") as handle:
            head = handle.read(4)
    except OSError as exc:
        return False, f"cannot read --twoBit_path {path}: {exc}"
    if len(head) < 4:
        return False, f"--twoBit_path {path} is too short to be a 2bit file"
    little = struct.unpack("<I", head)[0]
    big = struct.unpack(">I", head)[0]
    if little not in TWOBIT_MAGIC and big not in TWOBIT_MAGIC:
        return False, (f"--twoBit_path {path} is not a 2bit file "
                       f"(signature 0x{little:08X}, expected 0x1A412743)")
    #logger.info("2bit file %s has a valid signature", path)
    return True, f"{path} is a valid 2bit file"
