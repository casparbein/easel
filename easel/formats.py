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

def read_fasta_headers(path):
    """Sequence identifiers in a FASTA file: each header up to the first space.

    The whole file is read -- headers are interleaved with sequence, not
    grouped at the top -- so this is the expensive half of the checks below.
    """
    names = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                parts = line[1:].split()
                if parts:
                    names.append(parts[0])
    return names


def _strip_extension(name, extensions):
    """*name* without a trailing (optionally .gz'd) extension from the set."""
    if name.endswith(".gz"):
        name = name[:-3]
    for ext in sorted(extensions, key=len, reverse=True):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def fasta_stems(fasta_dir, extensions=None):
    """{transcript stem: path} for the FASTA files in a directory."""
    extensions = extensions or FASTA_EXTENSIONS
    out = {}
    for path in sorted(Path(fasta_dir).iterdir()):
        if not path.is_file() or path.name.startswith("."):
            continue
        name = path.name[:-3] if path.name.endswith(".gz") else path.name
        if any(name.endswith(ext) for ext in extensions):
            out[_strip_extension(path.name, extensions)] = path
    return out


def compare_tree_headers(tips, headers, tree_label="the tree", reference=None):
    """Check that every FASTA header is a tip in the tree. Returns (ok, message).

    The direction is deliberately asymmetric:

      headers - tips   FATAL. The pipeline prunes each tree to its alignment's
                       taxa (prune_tips in scripts/newick_tree_manipulator.py),
                       and ete4 raises on a label that is not a tip -- so the
                       run dies in phase 2 with a message naming neither the
                       header nor the tree.
      tips - headers   NORMAL, not reported. Pruning a tree that covers more
                       species than one transcript is the entire point;
                       compare_tree_asm() already warns where that matters.
    """
    tips = set(tips)
    headers = set(headers)
    if reference:
        ## Renamed to REFERENCE inside the pipeline, but a normal tip in the
        ## input tree, so it is covered either way.
        headers = headers - {reference}
    absent = sorted(headers - tips)
    if not absent:
        return True, "{} covers all {} sequence header(s)".format(
            tree_label, len(headers))
    shown = ", ".join(absent[:5])
    if len(absent) > 5:
        shown += " ..."
    return False, ("{} sequence header(s) absent from {}: {}".format(
        len(absent), tree_label, shown))


def check_headers_against_tree(fasta_dir, tips, tree_label="the input tree",
                               reference=None, extensions=None, max_files=None):
    """Every header in every alignment must be a tip in one shared tree.

    Reports which alignments failed rather than a bare count: a naming
    mismatch shows up as every alignment failing on the same names, and that
    is a different fix from one transcript carrying one stray taxon.
    """
    stems = fasta_stems(fasta_dir, extensions)
    if not stems:
        return True, "no FASTA files under {} to check against {}".format(
            fasta_dir, tree_label)
    paths = list(stems.values())
    if max_files:
        paths = paths[:max_files]

    tips = set(tips)
    offenders, all_missing = [], set()
    for path in paths:
        headers = set(read_fasta_headers(path))
        if reference:
            headers = headers - {reference}
        missing = headers - tips
        if missing:
            offenders.append((path.name, sorted(missing)))
            all_missing |= missing

    if not offenders:
        return True, "{} covers the headers of all {} alignment(s)".format(
            tree_label, len(paths))
    return False, _coverage_message(offenders, all_missing, len(paths), tree_label)


def _coverage_message(offenders, all_missing, n_files, tree_label):
    """Shared wording for the two coverage checks."""
    shown = ", ".join(sorted(all_missing)[:6])
    if len(all_missing) > 6:
        shown += " ..."
    detail = "; ".join("{} ({})".format(name, ", ".join(miss[:3]))
                       for name, miss in offenders[:3])
    if len(offenders) > 3:
        detail += "; and {} more".format(len(offenders) - 3)
    systematic = ""
    if len(offenders) == n_files and n_files > 1:
        systematic = (" Every alignment is affected, which points at a naming "
                      "convention rather than at individual transcripts.")
    return ("{} of {} alignment(s) contain sequence headers absent from {}. "
            "Missing tip(s): {}.{} First offenders: {}. Every header must be a "
            "tip in the tree its alignment is screened against, or the run "
            "fails when that tree is pruned to the alignment.".format(
                len(offenders), n_files, tree_label, shown, systematic, detail))


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


def _match_stem(tree_name, tree_extensions, wanted):
    """Which alignment stem this tree file belongs to, or None.

    Prefers the exact stem, then the longest prefix of the file name that is an
    alignment stem -- so <id>.nh and <id>_iqtree.nh both pair with <id>. Tests
    prefixes of the tree name against the stem set rather than scanning the
    stems, which is O(name length) instead of O(alignments); the same trick,
    for the same reason, as read_gene_tree_names() in rules/common.smk.
    """
    stem = _strip_extension(tree_name, tree_extensions)
    if stem in wanted:
        return stem
    for i in range(len(tree_name) - 1, 0, -1):
        if tree_name[:i] in wanted:
            return tree_name[:i]
    return None


## This extends tree format checks for precomputed gene trees
def validate_tree_directory(input_dir, extensions=None, fasta_dir=None,
                               reference=None, fasta_extensions=None):
    """Every tree file in *input_dir* must parse. Returns (ok, message)."""
    extensions = extensions or TREE_EXTENSIONS
    files = [p for p in sorted(Path(input_dir).iterdir())
             if p.is_file() and not p.name.startswith(".")]
    if not files:
        return False, f"Directory '{input_dir}' contains no files."

    wanted = fasta_stems(fasta_dir, fasta_extensions) if fasta_dir else {}
    paired = set()
    offenders, all_missing = [], set()

    checked, skipped = 0, []
    for path in files:
        name = path.name[:-3] if path.name.endswith(".gz") else path.name
        if not any(name.endswith(ext) for ext in extensions):
            skipped.append(path.name)
            continue
        try:
            ## Maybe use read_in_tree and compare_tree_asm here (?)
            with open(path, "r", encoding="utf-8") as handle:
                tmp_tree = _tree_class()(handle, parser=1)
        except Exception as exc:
            return False, f"{path} is not readable as Newick: {exc}"
        checked += 1
        
        if wanted:
            stem = _match_stem(path.name, extensions, wanted)
            if stem is None:
                continue
            paired.add(stem)
            tips = {leaf.name for leaf in tmp_tree.leaves()}
            headers = set(read_fasta_headers(wanted[stem]))
            if reference:
                headers = headers - {reference}
            missing = headers - tips
            if missing:
                offenders.append((wanted[stem].name, sorted(missing)))
                all_missing |= missing

    if skipped:
        logger.warning("Ignored %d non-tree file(s) in %s: %s%s", len(skipped),
                       input_dir, ", ".join(skipped[:5]),
                       " ..." if len(skipped) > 5 else "")
    if not checked:
        return False, f"No tree files with a recognised extension in '{input_dir}'."
    
    problems = []
    if wanted:
        unpaired = sorted(set(wanted) - paired)
        if unpaired:
            shown = ", ".join(unpaired[:5]) + (" ..." if len(unpaired) > 5 else "")
            problems.append(
                f"{len(unpaired)} of {len(wanted)} alignment(s) have no matching "
                f"tree in '{input_dir}': {shown}")
        if offenders:
            problems.append(_coverage_message(
                offenders, all_missing, len(wanted), "their own gene tree"))
    if problems:
        return False, ". ".join(problems) + "."

    covered = (f"; all {len(wanted)} alignment(s) are covered by their own tree"
               if wanted else "")
    return True, (f"All {checked} tree files in '{input_dir}' parsed "
                  f"successfully{covered}.")



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
