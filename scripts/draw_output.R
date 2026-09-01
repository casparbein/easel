library(tidyverse)
library(patchwork)

## select_sites() falls back to sample(); fix the seed so figures reproduce.
set.seed(1)
library(Biostrings)
library(ggtree)
library(tidytree)
library(scales)

#is.waive <- function(x) inherits(x, "waiver")

## ── Color scheme ──────────────────────────────────────────────────────────────
clustal_colors <- c(
  "A" = "#80a0f0", "R" = "#f01505", "N" = "#00ff00", "D" = "#c048c0",
  "C" = "#f08080", "Q" = "#00ff00", "E" = "#c048c0", "G" = "#f09048",
  "H" = "#15a4a4", "I" = "#80a0f0", "L" = "#80a0f0", "K" = "#f01505",
  "M" = "#80a0f0", "F" = "#80a0f0", "P" = "#ffff00", "S" = "#00ff00",
  "T" = "#00ff00", "W" = "#80a0f0", "Y" = "#15a4a4", "V" = "#80a0f0",
  "-" = "white",   "X" = "grey"
)

## ── Page geometry ─────────────────────────────────────────────────────────────
PER_TAXON_IN    <- 0.16   # page height per alignment row
PER_CODON_IN    <- 0.34   # page width per codon column in the zoom panel
PER_SITE_IN     <- 0.012  # page width per codon column in the overview (no text)
CODON_TEXT_SIZE <- 2.4    # ggplot mm; ~6.8 pt, fits 3 chars in PER_CODON_IN
MIN_MAIN_W      <- 15
MIN_OVER_W      <- 25
MAX_PAGE_IN     <- 200    # R's pdf device refuses anything larger

## ── Helpers ───────────────────────────────────────────────────────────────────

# Returns FALSE when a snakemake input slot is absent/empty, otherwise the value
safe_input <- function(x) {
  if (is.null(x) || length(x) == 0) return(FALSE)
  vals <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(vals) == 0) FALSE else vals
}

# Parse a DNAStringSet into a long codon tibble
get_codon_data <- function(seqs) {
  seq_mat       <- as.matrix(seqs)
  codon_indices <- seq(1, ncol(seq_mat) - 2, by = 3)
  map_df(names(seqs), function(id) {
    s <- as.character(seqs[[id]])
    tibble(
      sequence_id = id,
      codon_pos   = seq_along(codon_indices),
      dna_pos     = codon_indices,
      codon       = substring(s, codon_indices, codon_indices + 2)
    )
  })
}

# Build a ggtree plot from a file path; returns list(plot, ordered_labels, scale_data)
build_tree_plot <- function(tree_path) {
  tree_obj <- read.tree(tree_path)
  ttib     <- as_tibble(tree_obj) %>%
    mutate(is_selected = grepl("\\{Selected\\}", label))
  p <- ggtree(as.treedata(ttib), aes(color = is_selected, size = is_selected)) +
    scale_color_manual(values = c("TRUE" = "red",  "FALSE" = "black")) +
    scale_size_manual (values = c("TRUE" = 1.2,    "FALSE" = 0.5)) +
    #theme(legend.position = "none") +
    guides(color = "none", size = "none") +
    geom_tiplab(aes(label = ""), align = TRUE)

  ordered <- p$data %>%
    arrange(desc(y)) %>%
    filter(!is.na(label), label != "") %>%
    mutate(label = str_remove(label, "\\{Selected\\}")) %>%
    pull(label)

  p$data <- p$data %>%
    arrange(y) %>%
    mutate(
      ## ape sets label = NA for unnamed internal nodes; NA != "" is NA, and
      ## cumsum propagates it through every later row, emptying the panel.
      labeled_index = cumsum(!is.na(label) & label != ""),
      raw_pos       = row_number(),
      y             = if_else(!is.na(label) & label != "", as.numeric(labeled_index), raw_pos - 0.5)
    )

  scale_data <- p$data %>%
    filter(!is.na(label), label != "") %>%
    mutate(label = str_remove(label, "\\{Selected\\}"))

  list(plot = p, ordered_labels = ordered, scale_data = scale_data)
}

# Returns a sorted numeric vector of up to n_top codon positions selected by value_col.
# Falls back to first n_top positions when value_col is absent/all-NA.
select_sites <- function(df, value_col = NULL, n_top = 20, n_random = 10) {
  all_pos <- sort(unique(df[["Codon Position"]]))

  use_values <- !is.null(value_col) &&
    value_col %in% names(df) &&
    any(!is.na(suppressWarnings(as.numeric(df[[value_col]]))))

  if (!use_values) return(head(all_pos, n_top))

  summary <- df %>%
    mutate(.val = suppressWarnings(as.numeric(.data[[value_col]]))) %>%
    group_by(`Codon Position`) %>%
    summarise(max_val = max(.val, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_val))

  top <- summary %>% slice_head(n = n_top) %>% pull(`Codon Position`)

  if (length(top) < n_top) {
    rest  <- setdiff(all_pos, top)
    extra <- if (length(rest) > 0) sample(rest, min(n_random, length(rest))) else integer(0)
    top   <- sort(unique(c(top, extra)))
  }
  sort(top)
}

### ── CLI Arguments ──────────────────────────────────────────────────────────────
#option_list <- list(
#  make_option(c("-a", "--alignment"), type = "character", default = NULL,
#              help = "Path to sequence alignment (required)"),
#  make_option(c("-o", "--output"), type = "character", default = "summary_plot.pdf",
#              help = "Path to main output PDF [default %default]"),
#  make_option(c("-t", "--absrel_tree"), type = "character", default = NULL,
#              help = "Path to aBSREL fitted tree (.nwk)"),
#  make_option(c("-r", "--absrel_er"), type = "character", default = NULL,
#              help = "Path to aBSREL branch-site data (.tsv)"),
#  make_option(c("-b", "--busted_er"), type = "character", default = NULL,
#              help = "Path to BUSTED branch-site data (.tsv)"),
#  make_option(c("-m", "--busted_tree"), type = "character", default = NULL,
#              help = "Path to BUSTED MG94 fitted tree (.nwk)"),
#  make_option(c("-e", "--meme"), type = "character", default = NULL,
#              help = "Path to MEME MLE data (.csv)"),
#  make_option(c("-c", "--classic"), type = "character", default = NULL,
#              help = "Path to Bayescode classic omega data"),
#  make_option(c("-u", "--mutsel"), type = "character", default = NULL,
#              help = "Path to Bayescode mutsel data")
#)
#opt_parser <- OptionParser(option_list = option_list)
#opt <- parse_args(opt_parser)
#
#if (is.null(opt$alignment)) {
#  print_help(opt_parser)
#  stop("Missing required argument: --alignment", call. = FALSE)
#}

dna_seqs_path <- snakemake@input[["alignment"]] #opt$alignment
output_pdf    <- snakemake@output[["summary_plot"]] #opt$output
overview_pdf  <- snakemake@output[["overview"]] # sub("\\.pdf$", "_overview.pdf", output_pdf)

absrel_tree_raw  <- safe_input(snakemake@input[["asbrel"]][2])
absrel_er_raw    <- safe_input(snakemake@input[["asbrel"]][1])
busted_er_raw    <- safe_input(snakemake@input[["busted"]][1])
busted_tree_raw  <- safe_input(snakemake@input[["busted"]][2])
relax_k          <- safe_input(snakemake@input[["relax_k"]])
meme_raw         <- safe_input(snakemake@input[["meme_mle"]])
classic <- snakemake@input[["bayescode"]][1]
mutsel <- snakemake@input[["bayescode"]][2]

# Bayescode requires both files
bayescode_raw <- FALSE
if (!is.null(classic) && !is.null(mutsel)) {
  bayescode_raw <- c(classic, mutsel)
}

has_absrel      <- !isFALSE(absrel_tree_raw)
has_absrel_er   <- !isFALSE(absrel_er_raw)
has_busted      <- !isFALSE(busted_er_raw)
has_busted_tree <- !isFALSE(busted_tree_raw)
has_meme        <- !isFALSE(meme_raw)
has_bayescode   <- !isFALSE(bayescode_raw) && length(bayescode_raw) >= 2
has_tree        <- has_absrel || has_busted_tree
has_relax       <- !isFALSE(relax_k)

message(
  "Input status — aBSREL: ", has_absrel, "  aBSREL-er: ", has_absrel_er,
  "  BUSTED: ", has_busted, "  BUSTED-tree: ", has_busted_tree, 
  "  MEME: ", has_meme, "  Bayescode: ", has_bayescode, "  RELAX: ", has_relax
)

## ── Alignment (always) ────────────────────────────────────────────────────────
dna_seqs <- readDNAStringSet(dna_seqs_path)
plot_data <- get_codon_data(dna_seqs)
all_seq_ids <- names(dna_seqs)

## ── Tree loading ──────────────────────────────────────────────────────────────
mg_tree_plot <- NULL
ordered_labels <- all_seq_ids # default Y ordering
tree_scale_y <- NULL

if (has_absrel) {
  message("Loading aBSREL tree...")
  tr <- build_tree_plot(absrel_tree_raw[1])
  mg_tree_plot <- tr$plot
  ordered_labels <- tr$ordered_labels
  tree_scale_y <- tr$scale_data
} else if (has_busted_tree) {
  message("Loading BUSTED MG94 tree...")
  tr <- build_tree_plot(busted_tree_raw[1])
  mg_tree_plot <- tr$plot
  ordered_labels <- tr$ordered_labels
  tree_scale_y <- tr$scale_data
}

## ── Branch-site data ──────────────────────────────────────────────────────────
busted_df <- NULL
absrel_df <- NULL
branch_msa <- NULL
has_ebf <- FALSE
selected_sites <- NULL

## extract_hyphy.py now separates the two quantities: EBF is the
## posterior/prior odds ratio, "Site ER" the likelihood ratio between models.
EBF_COL   <- "EBF"
LOGEBF_COL <- "log10 EBF"
ER_COL    <- "Site ER"

process_hyphy_branch_data <- function(filepath) {
  df <- read_delim(filepath, col_names = TRUE, show_col_types = FALSE)
  ## An undefined EBF is written as an empty cell now, not as 0.0
  has_stats <- EBF_COL %in% names(df) &&
    any(!is.na(suppressWarnings(as.numeric(df[[EBF_COL]]))))
  merged_msa <- df %>%
    full_join(plot_data, by = c("Codon Position" = "codon_pos", "Branch" = "sequence_id")) %>%
    mutate(
      Codon      = ifelse(is.na(codon), Codon, codon),
      aa_residue = GENETIC_CODE[Codon],
      aa_residue = ifelse(is.na(aa_residue), "X", aa_residue),
      Branch     = factor(Branch, levels = rev(ordered_labels))
    )
    
  list(df = df, merged = merged_msa, has_stats = has_stats)
}

if (has_busted) {
  message("Loading BUSTED data...")
  res <- process_hyphy_branch_data(busted_er_raw[1])
  busted_df  <- res$df
  branch_msa <- res$merged
  has_ebf    <- res$has_stats
  
  selected_sites <- if (has_ebf) {
    select_sites(busted_df, EBF_COL, n_top = 20)
  } else {
    NULL#head(sort(unique(busted_df[["Codon Position"]])), 20)
  }
} else if (has_absrel_er) {
  message("Loading aBSREL branch data...")
  res <- process_hyphy_branch_data(absrel_er_raw[1])
  absrel_df  <- res$df
  branch_msa <- res$merged
  has_ebf    <- res$has_stats
  
  selected_sites <- if (has_ebf) {
    select_sites(absrel_df, EBF_COL, n_top = 20)
  } else {
    NULL#head(sort(unique(absrel_df[["Codon Position"]])), 20)
  }
}

## ── MEME ──────────────────────────────────────────────────────────────────────
meme_df <- NULL

if (has_meme) {
  message("Loading MEME data...")
  meme_df <- read_delim(meme_raw[1], col_names = TRUE, show_col_types = FALSE) %>%
    mutate(`Codon Position` = row_number())

  if (is.null(selected_sites)) {
    selected_sites <- sort(
      meme_df %>% arrange(`p-value`) %>% slice_head(n = 20) %>% pull(`Codon Position`)
    )
  }
}

## ── Bayescode ─────────────────────────────────────────────────────────────────
both_bc <- NULL

if (has_bayescode) {
  message("Loading Bayescode data...")
  read_bc <- function(path, run_label) {
    read_delim(path, skip = 2, col_names = FALSE, show_col_types = FALSE) %>%
      dplyr::rename(site = X1, CI_0025 = X2, omega = X3, CI_0975 = X4) %>%
      mutate(run = run_label, site = as.integer(site))
  }

  classical_df <- read_bc(bayescode_raw[1], "classic")
  mutsel_df <- read_bc(bayescode_raw[2], "mutsel")
  both_bc <- bind_rows(classical_df, mutsel_df)

  if (is.null(selected_sites)) {
    selected_sites <- sort(
      classical_df %>%
        mutate(dev = abs(omega - 1)) %>%
        arrange(desc(dev)) %>%
        slice_head(n = 20) %>%
        pull(site)
    )
  }
}

## Fallback site selection
if (is.null(selected_sites)) {
  selected_sites <- head(sort(unique(plot_data$codon_pos)), 20)
}

message("Selected sites (first 5): ", paste(head(selected_sites, 5), collapse = ", "))

## ── RELAX ─────────────────────────────────────────────────────────────────────
relax_plot <- NULL
if (has_relax) {
  message("Loading RELAX data...")
  relax_data <- read_delim(relax_k, col_names = TRUE, show_col_types = FALSE)
  
  relax_data <- relax_data %>%
    filter(!is.na(`p-value`))
  
  if (nrow(relax_data) > 0) {
    # Ensure branch labels match the tree order if possible
    known_branches <- intersect(relax_data$Branch, ordered_labels)
    other_branches <- setdiff(relax_data$Branch, ordered_labels)
    relax_data <- relax_data %>%
      mutate(Branch = factor(Branch, levels = c(known_branches, other_branches)))
      
    relax_plot <- ggplot(relax_data, aes(x = Branch, y = k, color = `p-value`)) +
      geom_boxplot(outlier.shape = NA) +
      ## height = 0: without it the k values themselves are displaced vertically,
    ## so the plotted relaxation coefficients are not the ones in the table.
    geom_jitter(width = 0.2, height = 0) +
      expand_limits(y = 0) +
      theme_bw() +
      #scale_color_gradient2(low = "red", mid = "orange", high = "yellow", midpoint = 0.05, limits = c(0,1)) +
      scale_colour_gradientn(
        colours = c("red", "orange", "yellow", "grey"),
        values = c(0, 0.05, 0.1, 1),
        limits = c(0,1)
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.title.x = element_blank()
      ) +
      ylab("Relaxation (k)") +
      guides(
    #fill = guide_legend(title.position = "top", title.hjust = 0.5),
    color = guide_colorbar(title.position = "top", title.hjust = 0.5)
  )
  }
}

## ── Base MSA tibble for plotting ──────────────────────────────────────────────
msa_data <- if (!is.null(branch_msa)) {
  branch_msa
} else {
  plot_data %>%
    ## dplyr:: qualified: library(Biostrings) attaches S4Vectors, which masks
    ## dplyr::rename. This is the branch taken on any MEME-only or BayesCode-only
    ## run, so the unqualified call would have failed exactly there.
    dplyr::rename(`Codon Position` = codon_pos, Branch = sequence_id, Codon = codon) %>%
    mutate(
      aa_residue = GENETIC_CODE[Codon],
      aa_residue = ifelse(is.na(aa_residue), "X", aa_residue),
      Branch     = factor(Branch, levels = rev(ordered_labels))
    )
}

## ── Overview alignment (all positions, numerical x, no text, no EBF) ─────────
alignment_overview <- ggplot(
  msa_data,
  aes(x = as.numeric(`Codon Position`), y = Branch, fill = aa_residue)
) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(values = clustal_colors, guide = "none") +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  theme_minimal() +
  theme(
    axis.text.x     = element_text(size = 6),
    axis.text.y     = element_text(size = 5),
    panel.grid      = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Codon Position (all)", y = NULL)

## ── Detail alignment zoom (selected sites, factor x) ─────────────────────────
zoom_data <- msa_data %>% filter(`Codon Position` %in% selected_sites)

make_zoom <- function(dat, with_ebf) {
  p <- ggplot(dat, aes(x = as.factor(`Codon Position`), y = Branch, fill = aa_residue))

  if (with_ebf) {
    p <- p +
      geom_tile(
        aes(color = suppressWarnings(as.numeric(`log10 EBF`))),
        width = 0.9, height = 0.9
      ) +
      scale_color_gradient(
        low = "grey80", high = "black",
        limits = c(1, NA), 
        ## NA (transparent), not the string "NA", which is not a colour.
        na.value = NA,
        ## This column is an empirical Bayes factor, not an evidence ratio.
        name = "log10 EBF"
      )
  } else {
    p <- p + geom_tile(width = 0.95, height = 0.95)
  }

  p +
    geom_text(aes(label = Codon), size = CODON_TEXT_SIZE, fontface = "bold", color = "black") +
    scale_fill_manual(values = clustal_colors) +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5,ncol = 6),
            color = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(labels = setNames(as.character(tree_scale_y$label),
                                       as.character(tree_scale_y$label))) +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      panel.grid      = element_blank(),
      legend.position = "right",
      axis.text.y     = element_text(size = 5)
    ) +
    labs(x = "Codon Position", y = NULL, fill = "Amino Acid")
}

alignment_zoom <- make_zoom(zoom_data, with_ebf = has_ebf)

## ── Optional panels ───────────────────────────────────────────────────────────
mutsel_plot <- NULL
mutsel_overview <- NULL
if (!is.null(both_bc)) {
  # Main zoom plot
  mutsel_plot <- ggplot(
    both_bc %>% filter(site %in% selected_sites),
    aes(factor(site, levels = as.character(sort(selected_sites))),
      omega,
      color = run, group = run
    )
  ) +
    geom_point() +
    geom_errorbar(aes(ymin = CI_0025, ymax = CI_0975), width = 0.2, alpha = 0.6) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    scale_color_manual(
      name = "Model",
      labels = c("classic" = "Classic Omega", "mutsel" = "Mut.-Sel. Baseline"),
      values = c("darkorange", "gold")
    ) +
    ylab("Omega") +
        guides(
    #fill = guide_legend(title.position = "top", title.hjust = 0.5),
    color = guide_legend(title.position = "top", title.hjust = 0.5)
  )
    
  # Overview plot (all sites)
  mutsel_overview <- ggplot(
    both_bc,
    aes(as.numeric(site), omega, color = run, group = run)
  ) +
    geom_point(alpha = 0.5, size = 1) +
    geom_errorbar(aes(ymin = CI_0025, ymax = CI_0975), width = 0.2, alpha = 0.3) +
    theme_bw() +
    scale_x_continuous(expand = c(0.005, 0.005)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    scale_color_manual(
      name = "Model",
      labels = c("classic" = "Classic Omega", "mutsel" = "Mut.-Sel. Baseline"),
      values = c("darkorange", "gold")
    ) +
    ylab("Omega")
}

meme_plot <- NULL
meme_overview <- NULL
if (!is.null(meme_df)) {
  # Main zoom plot
  meme_plot <- ggplot(
    meme_df %>% filter(`Codon Position` %in% selected_sites),
    aes(factor(`Codon Position`, levels = as.character(sort(selected_sites))),
      -log10(`p-value`),
     color = `# branches under selection`
    )
  ) +
    geom_point() +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    scale_color_gradient(name = "Branches under selection", low = "blue", high = "red", 
                          breaks = breaks_width(1)) +
    #scale_colour_gradient2(name = "Branches under selection") +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    ylab("-Log10(p-value)") +
    expand_limits(y = 0) +
    guides(
    #fill = guide_legend(title.position = "top", title.hjust = 0.5),
    color = guide_colorbar(title.position = "top", title.hjust = 0.5)
  )
    
  # Overview plot (all sites)
  meme_overview <- ggplot(
    meme_df,
    aes(as.numeric(`Codon Position`), -log10(`p-value`), color = `# branches under selection`)
  ) +
    geom_point(alpha = 0.6, size = 1) +
    theme_bw() +
    scale_x_continuous(expand = c(0.005, 0.005)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    scale_color_gradient(name = "Branches under selection", low = "blue", high = "red",
                          breaks = breaks_width(1)) +
    #scale_colour_gradient2(name = "Branches under selection") +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    ylab("-Log10(p-value)")
}

## ── Layout composition (Patchwork internal design) ──────────────────────────

# Build main plot layout
main_plots  <- list()
main_design <- character()
main_heights <- numeric()

# 1) Top panels
#if (!is.null(relax_plot)) {
#  main_plots <- c(main_plots, list(R = relax_plot))
#  main_design <- c(main_design, if (has_tree) "##RRRRR" else "RRRRRRR")
#  main_heights <- c(main_heights, 0.8)
#}
#
#if (!is.null(mutsel_plot)) {
#  main_plots <- c(main_plots, list(M = mutsel_plot))
#  main_design <- c(main_design, if (has_tree) "##MMMMM" else "MMMMMMM")
#  main_heights <- c(main_heights, 0.4)
#}
#
#if (!is.null(meme_plot)) {
#  main_plots <- c(main_plots, list(E = meme_plot))
#  main_design <- c(main_design, if (has_tree) "##EEEEE" else "EEEEEEE")
#  main_heights <- c(main_heights, 0.4)
#}
#
## 2) Main tree and zoom panels
#if (!is.null(mg_tree_plot)) {
#  main_plots <- c(main_plots, list(T = mg_tree_plot, Z = alignment_zoom))
#  main_design <- c(main_design, "TTZZZZZ")
#} else {
#  main_plots <- c(main_plots, list(Z = alignment_zoom))
#  main_design <- c(main_design, "ZZZZZZZ")
#}
#main_heights <- c(main_heights, 3)
#
#output_plot <- wrap_plots(main_plots, design = paste(main_design, collapse = "\n")) +
#  plot_layout(heights = main_heights, guides = "collect") &
#  theme(legend.position = "bottom")

# 1) Initialize containers
main_plots <- list()
main_design <- c()
main_heights <- c()

# 2) Identify which plots are present
has_relax  <- !is.null(relax_plot)
has_mutsel <- !is.null(mutsel_plot)
has_meme   <- !is.null(meme_plot)
has_side   <- has_mutsel | has_meme # Check if there's anything to put next to Relax

# 3) Handle Relax Plot (R)
if (has_relax) {
  main_plots <- c(main_plots, list(R = relax_plot))
  
  # If NO side plots exist, give Relax its own full row now
  if (!has_side) {
    main_design <- c(main_design, if (has_tree) "##RRRRR" else "RRRRRRR")
    main_heights <- c(main_heights, 0.8)
  }
}

# 4) Handle Mutsel Plot (M)
if (has_mutsel) {
  main_plots <- c(main_plots, list(M = mutsel_plot))
  # If Relax exists, it takes the first 2 slots ("RR"), otherwise uses empty space "##"
  row_str <- if (has_relax) "RRMMMMM" else (if (has_tree) "##MMMMM" else "MMMMMMM")
  main_design <- c(main_design, row_str)
  main_heights <- c(main_heights, 0.4)
}

# 5) Handle Meme Plot (E)
if (has_meme) {
  main_plots <- c(main_plots, list(E = meme_plot))
  # If Relax exists, it continues in the first 2 slots ("RR")
  row_str <- if (has_relax) "RREEEEE" else (if (has_tree) "##EEEEE" else "EEEEEEE")
  main_design <- c(main_design, row_str)
  main_heights <- c(main_heights, 0.4)
}

# 6) Main tree and zoom panels
if (!is.null(mg_tree_plot)) {
  main_plots <- c(main_plots, list(T = mg_tree_plot, Z = alignment_zoom))
  main_design <- c(main_design, "TTZZZZZ")
} else {
  main_plots <- c(main_plots, list(Z = alignment_zoom))
  main_design <- c(main_design, "ZZZZZZZ")
}
main_heights <- c(main_heights, 3)

# 7) Assemble using patchwork
output_plot <- wrap_plots(main_plots, design = paste(main_design, collapse = "\n")) +
  plot_layout(heights = main_heights, guides = "collect") &
  theme(
      legend.position = "bottom",
      legend.box = "horizontal",       # Place different legend blocks side-by-side
      legend.direction = "horizontal", # Flow items horizontally
      legend.margin = margin(t = 0),    # Remove top margin to save space
      legend.spacing.x = unit(0.2, "cm"), 
      legend.text = element_text(size = 8)
    )
#theme(legend.position = "bottom")

# Build overview layout
over_plots   <- list()
over_design  <- character()
over_heights <- numeric()

if (!is.null(mutsel_overview)) {
  over_plots <- c(over_plots, list(M = mutsel_overview))
  over_design <- c(over_design, "M")
  over_heights <- c(over_heights, 1)
}

if (!is.null(meme_overview)) {
  over_plots <- c(over_plots, list(E = meme_overview))
  over_design <- c(over_design, "E")
  over_heights <- c(over_heights, 1)
}

over_plots <- c(over_plots, list(A = alignment_overview))
over_design <- c(over_design, "A")
over_heights <- c(over_heights, 3)

overview_plot <- wrap_plots(over_plots, design = paste(over_design, collapse = "\n")) +
  plot_layout(heights = over_heights, guides = "collect") &
  theme(legend.position = "bottom")

## ── Save ──────────────────────────────────────────────────────────────────────
n_taxa <- n_distinct(msa_data$Branch)
n_zoom <- n_distinct(zoom_data[["Codon Position"]])
n_all  <- n_distinct(msa_data[["Codon Position"]])

clamp <- function(x) min(MAX_PAGE_IN, x)

zoom_w_frac <- if (!is.null(mg_tree_plot)) 5 / 7 else 1
main_h_frac <- 3 / sum(main_heights)
plot_w <- clamp(max(MIN_MAIN_W, (n_zoom * PER_CODON_IN) / zoom_w_frac))
plot_h <- clamp(max(7 + 2 * (length(main_heights) - 1),
                    (n_taxa * PER_TAXON_IN) / main_h_frac))

over_h_frac <- 3 / sum(over_heights)
over_w <- clamp(max(MIN_OVER_W, n_all * PER_SITE_IN))
over_h <- clamp(max(5 + 2 * (length(over_heights) - 1),
                    (n_taxa * PER_TAXON_IN) / over_h_frac))

if (max(plot_w, plot_h, over_w, over_h) >= MAX_PAGE_IN)
  message("Note: page size hit the ", MAX_PAGE_IN,
          " in device limit; this alignment may still look crowded")

## limitsize = FALSE: ggsave refuses anything over 50 in by default, which
## these pages legitimately exceed for large alignments.
ggsave(output_pdf, plot = output_plot, device = "pdf",
       width = plot_w, height = plot_h, units = "in", limitsize = FALSE)
message(sprintf("Main detail plot saved to: %s  (%.1f x %.1f in; %d taxa x %d sites)",
                output_pdf, plot_w, plot_h, n_taxa, n_zoom))

ggsave(overview_pdf, plot = overview_plot, device = "pdf",
       width = over_w, height = over_h, units = "in", limitsize = FALSE)
message(sprintf("Overview plot saved to: %s  (%.1f x %.1f in; %d taxa x %d sites)",
                overview_pdf, over_w, over_h, n_taxa, n_all))
