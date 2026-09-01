#!/usr/bin/env Rscript
## Visualise what each cleaning stage did to one transcript's alignment, and
## write a high-level cleaning report for the user to investigate
##
## The pipeline chain, earliest first:
##   _ori.fa                  original fasta file         <- (reference grid)
##   _ren.fa                  REFERENCE renamed only      <- reference grid
##   .masked.fa               codonify                    masks CELLS
##   .masked.hmm_cleaned.fa   HmmCleaner + transferCleaner masks CELLS
##   .manual.fa               manual cleaner              drops ROWS and COLS
##
## _ori.fa is accepted only as a fallback grid when _ren.fa is absent (the two
## differ in the reference header, not in the matrix). 
##
## Codonify and HmmCleaner mask individual matrix entries; they do not drop
## sequences, and only codonify can occasionally change the column count when
## removing columns that consist of gaps and Ns only.
## Only manual cleaner routinely drops whole rows and
## columns. The report states all three numbers -- cells masked, columns
## removed, sequences removed -- for every stage
##
## Everything is placed on the _ren.fa/_ori.fa grid. Column indices are carried
## forward through each transition, so a stage that both drops columns and
## masks cells (codonify sometimes) still reports and places both.
##
## Standalone (no arguments needed inside a transcript's tmp/):
##   Rscript plot_cleaning.R
##   Rscript plot_cleaning.R --dir codon_alignments/ENST1/tmp --id ENST1
## As a snakemake script: reads snakemake@input / @output instead.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(Biostrings)
})

## ── Stage chain ──────────────────────────────────────────────────────────────
## Order of alignment based on pipeline
CHAIN <- c(ori         = "_ori.fa",
           ren         = "_ren.fa",
           masked      = ".masked.fa",
           hmm_cleaned = ".masked.hmm_cleaned.fa",
           manual      = ".manual.fa")

## What each transition means, for the report.
TRANSITION_LABEL <- c(
  "ori->ren"              = "Reference renaming",
  "ori->masked"           = "Codonify",
  "ren->masked"           = "Codonify",
  "masked->hmm_cleaned"   = "HmmCleaner (via transferCleaner)",
  "ren->hmm_cleaned"      = "HmmCleaner (via transferCleaner)",
  "ori->hmm_cleaner"      = "HmmCleaner (via transferCleaner)",
  "hmm_cleaned->manual"   = "Manual cleaner",
  "masked->manual"        = "Manual cleaner",
  "ren->manual"           = "Manual cleaner",
  "ori->manual"           = "Manual cleaner"
)

## Longest suffix first, so ".masked.hmm_cleaned.fa" is never read as
## ".masked.fa" and "_ori.fa" never shadows "_ren.fa".
SUFFIXES_BY_LENGTH <- CHAIN[order(-nchar(CHAIN))]

## ── Arguments ────────────────────────────────────────────────────────────────
## Accepts both "--key value" and "--key=value".
parse_cli <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (startsWith(a, "--")) {
      if (grepl("=", a, fixed = TRUE)) {
        out[[sub("^--([^=]+)=.*$", "\\1", a)]] <- sub("^--[^=]+=", "", a)
      } else if (i < length(args) && !startsWith(args[i + 1], "--")) {
        out[[sub("^--", "", a)]] <- args[i + 1]
        i <- i + 1L
      } else {
        out[[sub("^--", "", a)]] <- TRUE
      }
    }
    i <- i + 1L
  }
  out
}

## Recover the transcript id from whatever alignment files are present, so
## running this inside a transcript's tmp/ needs no --id at all.
infer_id <- function(dir) {
  files <- list.files(dir)
  for (sfx in SUFFIXES_BY_LENGTH) {
    hit <- files[endsWith(files, sfx)]
    if (length(hit) > 0) {
      ids <- unique(substr(hit, 1, nchar(hit) - nchar(sfx)))
      if (length(ids) > 1) warning("several transcripts in ", dir, "; using ", ids[1])
      return(ids[1])
    }
  }
  NULL
}

## For a snakemake run
if (exists("snakemake")) {
  opt <- list(dir = NULL, id = snakemake@wildcards$transcript_id,
              out = snakemake@output$plot, report = snakemake@output$report)
  paths <- as.list(snakemake@input)[names(CHAIN)]
} else {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  dir <- cli$dir %||% "."
  id  <- cli$id  %||% infer_id(dir)
  if (is.null(id)) {
    stop("could not work out the transcript id: no file in '", dir,
         "' ends in one of ", paste(CHAIN, collapse = ", "),
         ".\n  Pass --id explicitly, or point --dir at a transcript's tmp/ .",
         call. = FALSE)
  }
  if (is.null(cli$id)) message("inferred transcript id: ", id)
  opt   <- list(dir = dir, id = id,
                out    = cli$out    %||% paste0(id, ".cleaning.pdf"),
                report = cli$report %||% paste0(id, ".cleaning.txt"))
  paths <- map(CHAIN, ~ file.path(dir, paste0(id, .x)))
  ## Any stage can be overridden with --<stage> /path/to/file
  for (nm in names(CHAIN)) if (!is.null(cli[[nm]])) paths[[nm]] <- cli[[nm]]
}

## ── Reading ──────────────────────────────────────────────────────────────────
read_matrix <- function(path) {
  if (is.null(path) || length(path) != 1 || is.na(path)) return(NULL)
  if (!file.exists(path) || file.size(path) == 0) return(NULL)
  seqs <- readBStringSet(path)
  if (length(seqs) == 0) return(NULL)
  if (length(unique(width(seqs))) != 1) {
    warning(basename(path), " is not aligned; skipping")
    return(NULL)
  }
  m <- do.call(rbind, strsplit(as.character(seqs), "", fixed = TRUE))
  rownames(m) <- sub("\\s.*$", "", names(seqs))
  m
}

## HmmCleaner/transferCleaner rewrite gaps as '*' (the bulk
## of them) but can also produce 'N' as well (through transferCleaner's -delChar). 
MASK_CHARS <- c("N", "n", "-", ".", "?", "X", "x", "*")


keep_shape <- function(x, res) {
  if (!is.null(dim(x))) { dim(res) <- dim(x); dimnames(res) <- dimnames(x) }
  res
}
is_mask <- function(x) keep_shape(x, x %in% MASK_CHARS)
is_gap  <- function(x) keep_shape(x, x %in% c("-", "."))

## A cell was masked if it held a real nucleotide before and a masking character
## after. The tools differ in which character they write -- codonify and
## HmmCleaner/transferCleaner use N and Manual cleaner use a gap -- so rather
## than hardcoding one per stage any masking character counts, 
## and the report names the ones actually seen for each stage. This might lead 
## to cases were actual maskings and gaps are confounded in the report
## but normally it is clear where grid cells/columns/rows were removed.
##
## Case is deliberately NOT treated as masking: none of these tools soft-mask
## by lowercasing, and genomic sequence routinely carries lowercase for
## soft-masked repeats, so reading a case change as masking invents findings.
##
## Two other non-events, both deliberate. A gap becoming an N is not counted:
## there was no nucleotide there to mask, and counting it would inflate
## codonify's number with cells that never held data. Nor is one masking
## character becoming another (N -> gap), which is a later stage restating an
## earlier stage's decision rather than masking anything new (N -> gap should
## not happen in this pipeline, this is a safeguard).
masked_cells <- function(before, after) {
  res <- !is_mask(before) & is_mask(after)
  dim(res)      <- dim(before)
  dimnames(res) <- dimnames(before)
  res
}

stages <- map(paths, read_matrix)
present <- names(CHAIN)[!map_lgl(stages[names(CHAIN)], is.null)]
if (length(present) == 0) {
  tried <- unlist(paths)
  stop("no readable alignment found for '", opt$id, "'.\n  looked for:\n",
       paste0("    ", names(tried), ": ", tried,
              ifelse(file.exists(tried), "  [exists but unreadable/empty]", "  [missing]"),
              collapse = "\n"), call. = FALSE)
}

## ── Reference grid: the raw/renamed alignment, per the pipeline's own start ──
grid_name <- if ("ren" %in% present) "ren" else if ("ori" %in% present) "ori" else present[1]
grid      <- stages[[grid_name]]
if (!grid_name %in% c("ren", "ori")) {
  warning("neither _ren.fa nor _ori.fa found; falling back to '", grid_name,
          "' as the reference grid")
}
message("reference grid: ", grid_name, " (", nrow(grid), " seqs x ", ncol(grid), " cols)")

## ── Walk the chain, one transition at a time ─────────────────────────────────
## colmap[[stage]][j] = which grid column stage-column j corresponds to, or NA
## once a stage's columns can no longer be traced back (a width increase, 
## which should not happen).
status <- matrix("retained", nrow(grid), ncol(grid), dimnames = dimnames(grid))
status[is_gap(grid)] <- "gap"

colmap <- list()
colmap[[grid_name]] <- seq_len(ncol(grid))
summaries <- list()

## Columns are only ever deleted, never reordered, so a two-pointer walk
## recovers which ones went. Comparison ignores positions the stage masked
## itself, otherwise a masked column looks like a deleted one.
map_removed_columns <- function(P, C) {
  ## Tolerate the target being masked ANYWHERE, not only where the source was
  ## an unmasked base.A mask in the target is never evidence that a column was dropped.
  cols_match <- function(a, b) all(a == b | is_mask(b))
  removed <- integer(0); i <- 1L; j <- 1L
  while (i <= ncol(P) && j <= ncol(C)) {
    if (cols_match(P[, i], C[, j])) { i <- i + 1L; j <- j + 1L }
    else { removed <- c(removed, i); i <- i + 1L }
  }
  if (j <= ncol(C)) {
    return(list(removed = integer(0), ok = FALSE,
                why = sprintf("ran out of source columns at %d/%d with %d of %d target columns unmatched",
                              i - 1L, ncol(P), ncol(C) - j + 1L, ncol(C))))
  }
  if (i <= ncol(P)) removed <- c(removed, seq(i, ncol(P)))
  list(removed = removed, ok = TRUE, why = NULL)
}

## Compare A's columns `acols` against all of B's, classify the differences as
## masking, and paint them onto `status` at grid columns `gridcols`. Shared by
## the equal-width and the column-dropping branches, since a stage can both
## drop columns and mask cells in the ones it keeps.
classify_masking <- function(A, B, shared, acols, gridcols, to_stage, status) {
  Asub <- A[shared, acols, drop = FALSE]
  Bsub <- B[shared, , drop = FALSE]
  if (ncol(Asub) != ncol(Bsub)) {
    return(list(status = status, n = 0L, per_seq = integer(0),
                note = sprintf("could not line up columns for masking (%d vs %d)",
                               ncol(Asub), ncol(Bsub))))
  }
  hits    <- masked_cells(Asub, Bsub)
  changed <- Asub != Bsub

  ## Cells that changed without being masked, split by kind. The big one is a
  ## masking character being rewritten as another: HmmCleaner cannot handle
  ## '-' and rewrites gaps as *. Those cells never held a nucleotide, so
  ## it must not count as masking and must
  ## not be reported as if something were wrong.
  other <- changed & !hits
  note  <- NULL
  if (sum(other) > 0) {
    ob <- Asub[other]; oa <- Bsub[other]
    mb <- is_mask(ob); mo <- is_mask(oa)
    k_rewrite  <- mb & mo                                    # N <-> gap: format only
    k_unmasked <- mb & !mo                                   # N/gap -> base: masking LOST (should not happen)
    k_case     <- !mb & !mo & toupper(ob) == toupper(oa)     # A -> a 
    k_subst    <- !mb & !mo & toupper(ob) != toupper(oa)     # A -> G: real base change (should not happen)
    ex <- function(sel)
      paste(head(unique(sprintf("%s->%s", ob[sel], oa[sel])), 4), collapse = ", ")
    parts <- character(0)
    if (any(k_rewrite))
      parts <- c(parts, sprintf("%d gap/N rewrites (ignored, not masking)", sum(k_rewrite)))
    if (any(k_unmasked))
      parts <- c(parts, sprintf("%d cells UNMASKED, earlier masking lost (%s)",
                                sum(k_unmasked), ex(k_unmasked)))
    if (any(k_case))
      parts <- c(parts, sprintf("%d case changes (%s)", sum(k_case), ex(k_case)))
    if (any(k_subst))
      parts <- c(parts, sprintf("%d base substitutions (%s)", sum(k_subst), ex(k_subst)))
    note <- paste(parts, collapse = "; ")
    if (sum(hits) == 0 && any(k_case | k_subst | k_unmasked))
      note <- paste("no masking detected;", note)
  }
  ## Which character this stage actually masks with, counted from the data.
  chars <- if (sum(hits) > 0) {
    tb <- sort(table(Bsub[hits]), decreasing = TRUE)
    paste(sprintf("%s (%d)", names(tb), as.integer(tb)), collapse = ", ")
  } else NA_character_
  label <- switch(to_stage, masked = "codonify_masked",
                            hmm_cleaned = "hmm_masked",
                            manual = "manual_masked", NA_character_)
  if (!is.na(label) && sum(hits) > 0)
    status <- paint(status, shared, gridcols, hits, label)
  list(status = status, n = sum(hits), chars = chars,
       per_seq = setNames(rowSums(hits), shared), note = note)
}

paint <- function(status, rows, gridcols, hits, label) {
  if (is.null(dim(hits)) || nrow(hits) != length(rows) ||
      ncol(hits) != length(gridcols)) {
    stop("paint(): hits is not a ", length(rows), " x ", length(gridcols),
         " matrix (got ", paste(dim(hits) %||% length(hits), collapse = "x"),
         ")", call. = FALSE)
  }
  for (jj in seq_along(gridcols)) {
    g <- gridcols[jj]
    if (is.na(g)) next
    sel <- rows[hits[, jj]]
    if (length(sel)) status[sel, g] <- label
  }
  status
}

## Start the walk AT the reference grid, not at the earliest file on disk.
chain <- present[seq(match(grid_name, present), length(present))]
for (k in seq_len(length(chain) - 1L)) {
  a <- chain[k]; b <- chain[k + 1L]
  A <- stages[[a]]; B <- stages[[b]]
  key <- paste0(a, "->", b)
  lab <- unname(TRANSITION_LABEL[key])
  if (is.na(lab)) lab <- paste(a, "->", b)

  rows_removed <- setdiff(rownames(A), rownames(B))
  shared       <- intersect(rownames(A), rownames(B))
  s <- list(from = a, to = b, label = lab,
            rows_removed = rows_removed, cols_removed = integer(0),
            cells_masked = 0L, per_seq = integer(0), chars = NA_character_,
            note = NULL)

  map_a <- colmap[[a]]
  if (is.null(map_a) || length(map_a) != ncol(A)) map_a <- rep(NA_integer_, ncol(A))

  if (length(shared) == 0) {
    s$note <- sprintf("no shared sequence names (%d vs %d)", nrow(A), nrow(B))
    colmap[[b]] <- rep(NA_integer_, ncol(B))

  } else if (ncol(A) == ncol(B)) {
    res <- classify_masking(A, B, shared, seq_len(ncol(A)), map_a, b, status)
    s$cells_masked <- res$n; s$per_seq <- res$per_seq
    s$chars <- res$chars; s$note <- res$note
    status <- res$status
    colmap[[b]] <- map_a

  } else if (ncol(B) < ncol(A)) {
    ## Columns disappeared. A stage can do both -- codonify may drop a column
    ## while fixing a frame AND mask cells in the columns it keeps -- so the
    ## surviving columns are still compared cell by cell.
    mp <- map_removed_columns(A[shared, , drop = FALSE], B[shared, , drop = FALSE])
    if (!mp$ok) {
      s$note <- paste("could not trace which columns were removed; omitted from",
                      "the figure --", mp$why)
      colmap[[b]] <- rep(NA_integer_, ncol(B))
    } else {
      s$cols_removed <- mp$removed
      kept <- setdiff(seq_len(ncol(A)), mp$removed)
      ## Masking first, then paint the dropped columns over the top: a dropped
      ## column outranks a mask, since the column is simply gone.
      res <- classify_masking(A, B, shared, kept, map_a[kept], b, status)
      s$cells_masked <- res$n; s$per_seq <- res$per_seq
    s$chars <- res$chars; s$note <- res$note
      status <- res$status
      gcols <- map_a[mp$removed]
      gcols <- gcols[!is.na(gcols)]
      if (length(gcols)) status[, gcols] <- "dropped_col"
      colmap[[b]] <- map_a[kept]
    }

  } else {
    ## Width grew (should not happen)
    s$note <- sprintf("alignment width grew %d -> %d; later stages cannot be mapped onto the %s grid",
                      ncol(A), ncol(B), grid_name)
    colmap[[b]] <- rep(NA_integer_, ncol(B))
  }

  if (length(rows_removed)) {
    keep <- intersect(rows_removed, rownames(status))
    if (length(keep)) status[keep, ] <- "dropped_row"
  }
  summaries[[key]] <- s
}

## Consistency check: the figure and the report are built from different
## objects (status vs the per-stage summaries), so they could disagree (should not happen)
reported_masked <- sum(map_dbl(summaries, ~ .x$cells_masked))
painted_masked  <- sum(status %in% c("codonify_masked", "hmm_masked", "manual_masked"))
if (reported_masked > 0 && painted_masked == 0) {
  warning("the report counts ", reported_masked, " masked cells but none were ",
          "placed on the grid -- the figure will show no masking. This is a bug ",
          "in the script, not a property of your data.")
} else if (painted_masked < reported_masked) {
  message(sprintf("note: %d of %d masked cells are hidden by a later stage's ",
                  reported_masked - painted_masked, reported_masked),
          "dropped row/column")
}

## ── Plot ─────────────────────────────────────────────────────────────────────
lvls <- c("retained", "gap", "codonify_masked", "hmm_masked", "manual_masked",
          "dropped_col", "dropped_row")
labs <- c(retained = "retained", gap = "gap (original)",
          codonify_masked = "cell masked by codonify",
          hmm_masked      = "cell masked by HmmCleaner",
          manual_masked   = "cell masked by manual cleaner (stop/unsafe codon)",
          dropped_col     = "column dropped (codonify/manual cleaner)",
          dropped_row     = "sequence dropped (manual cleaner)")
pal <- c(retained = "#cfd8e3", gap = "#f7f7f7", codonify_masked = "#f0a30a",
         hmm_masked = "#d1495b", manual_masked = "#8e44ad",
         dropped_col = "#5b5b5b", dropped_row = "#1f1f1f")

affected     <- !status %in% c("retained", "gap")
dim(affected) <- dim(status)
dimnames(affected) <- dimnames(status)
removed_frac <- rowMeans(affected)
row_levels   <- names(sort(removed_frac, decreasing = TRUE))

long <- tibble(
  seq    = factor(rep(rownames(status), times = ncol(status)), levels = rev(row_levels)),
  col    = rep(seq_len(ncol(status)), each = nrow(status)),
  status = factor(as.vector(status), levels = lvls)
)

main <- ggplot(long, aes(col, seq, fill = status)) +
    geom_raster() +
    scale_fill_manual(values = pal, labels = labs, name = NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = sprintf("alignment column (nucleotides)"), y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), legend.position = "bottom",
          axis.text.y = element_text(size = 6))

col_bar <- long %>%
  group_by(col) %>%
  summarise(frac = mean(!status %in% c("retained", "gap")), .groups = "drop") %>%
  ggplot(aes(col, frac)) +
  geom_col(width = 1, fill = "darkgrey") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(y = "per column", x = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor   = element_blank(),                          
        panel.grid.major.x = element_blank(),                          
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
        axis.ticks.y       = element_line(color = "grey50", linewidth = 0.3),
        axis.text.x        = element_blank())

row_bar <- tibble(seq = factor(names(removed_frac), levels = rev(row_levels)),
                  frac = removed_frac) %>%
  ggplot(aes(frac, seq)) +
  geom_col(fill = "darkgrey") +
  scale_x_continuous() +
  labs(x = "per sequence", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.y = element_blank())

combined <- (col_bar + plot_spacer() + main + row_bar) +
  plot_layout(ncol = 2, widths = c(1, 0.18), heights = c(0.18, 1)) +
  plot_annotation(title = paste("Alignment cleaning:", opt$id),
                  subtitle = paste("stages:", paste(present, collapse = " -> ")))

h <- max(3.5, min(20, 1.5 + nrow(grid) * 0.12))
ggsave(opt$out, combined, width = 11, height = h, limitsize = FALSE)
message("wrote ", opt$out)

## ── Report ───────────────────────────────────────────────────────────────────
pct   <- function(n, d) if (d > 0) sprintf("%.2f%%", 100 * n / d) else "n/a"
cells <- length(status)

lines <- c(
  sprintf("Cleaning report: %s", opt$id),
  strrep("=", 66), "",
  "Dimensions by stage (sequences x columns)",
  map_chr(present, ~ sprintf("  %-12s %4d x %6d%s", .x,
                             nrow(stages[[.x]]), ncol(stages[[.x]]),
                             if (.x == grid_name) "   <- reference grid" else "")),
  "",
  sprintf("Reference grid: %s, %d cells", grid_name, cells),
  "",
  "Per stage. Codonify can drop individual columns after masking (if only N or - are retained in that column),"
  "HmmCleaner masks individual nucleotides; only",
  "Manual cleaner routinely drops whole sequences or columns; masked cells are either not allowed (AC-) or stop codons (TGA).",
  strrep("-", 66)
)

for (s in summaries) {
  if (s$label == "Reference renaming") next
  n_seq_aff <- sum(s$per_seq > 0)
  top <- if (n_seq_aff > 0) {
    o <- sort(s$per_seq[s$per_seq > 0], decreasing = TRUE)
    paste(sprintf("%s (%d)", names(head(o, 3)), head(o, 3)), collapse = ", ")
  } else "-"
  lines <- c(lines,
    sprintf("%s   [%s -> %s]", s$label, s$from, s$to),
    sprintf("  cells masked      : %d (%s of grid)", s$cells_masked,
            pct(s$cells_masked, cells)),
    sprintf("  masking character : %s",
            if (is.na(s$chars)) "-" else s$chars),
    sprintf("  columns removed   : %d%s", length(s$cols_removed),
            if (length(s$cols_removed)) sprintf(" (%s of grid)",
                pct(length(s$cols_removed), ncol(grid))) else ""),
    sprintf("  sequences removed : %d%s", length(s$rows_removed),
            if (length(s$rows_removed))
              paste0(" (", paste(head(s$rows_removed, 5), collapse = ", "),
                     if (length(s$rows_removed) > 5) ", ..." else "", ")") else ""),
    sprintf("  sequences affected: %d of %d", n_seq_aff, nrow(stages[[s$from]])),
    sprintf("  most masked       : %s", top))
  if (!is.null(s$note)) lines <- c(lines, sprintf("  NOTE              : %s", s$note))
  lines <- c(lines, "")
}

final         <- stages[[tail(present, 1)]]
masked_labels <- c("codonify_masked", "hmm_masked", "manual_masked")
events        <- sum(map_dbl(summaries, ~ .x$cells_masked))
visible       <- sum(status %in% masked_labels)
dropped_cells <- sum(status %in% c("dropped_col", "dropped_row"))
cols_dropped  <- sum(map_dbl(summaries, ~ length(.x$cols_removed)))
rows_dropped  <- length(unique(unlist(lapply(summaries, function(z) z$rows_removed))))

lines <- c(lines, strrep("-", 66),
  sprintf("Overall: %d x %d (%s)  ->  %d x %d (%s)",
          nrow(grid), ncol(grid), grid_name,
          nrow(final), ncol(final), tail(present, 1)),
  "",
  "The per-stage counts above are masking EVENTS. The figure shows the LAST",
  "thing that happened to each cell, so its totals are lower: a cell masked",
  "by codonify and masked again by HmmCleaner is two events but one cell, and",
  "a masked cell inside a column that a later stage drops is shown as",
  "dropped rather than masked.",
  "",
  sprintf("  masking events, summed over stages : %d", events),
  sprintf("  cells still shown as masked        : %d (%s of grid)",
          visible, pct(visible, cells)),
  sprintf("  masked, then re-masked or dropped  : %d", events - visible),
  "",
  sprintf("  columns dropped                    : %d of %d (%s)",
          cols_dropped, ncol(grid), pct(cols_dropped, ncol(grid))),
  sprintf("  sequences dropped                  : %d of %d", rows_dropped, nrow(grid)),
  sprintf("  grid cells in dropped columns/rows : %d (%s of grid)",
          dropped_cells, pct(dropped_cells, cells)))

if (rows_dropped == 0 && dropped_cells == cols_dropped * nrow(grid)) {
  lines <- c(lines, sprintf("      = %d columns x %d sequences, consistent",
                            cols_dropped, nrow(grid)))
}

## Independent check on the column bookkeeping.
lines <- c(lines, "",
  sprintf("  column arithmetic: %d - %d dropped = %d, final width %d  [%s]",
          ncol(grid), cols_dropped, ncol(grid) - cols_dropped, ncol(final),
          if (ncol(grid) - cols_dropped == ncol(final)) "consistent"
          else "MISMATCH -- column mapping is wrong, treat positions as unreliable"))

writeLines(unlist(lines), opt$report)
message("wrote ", opt$report)
