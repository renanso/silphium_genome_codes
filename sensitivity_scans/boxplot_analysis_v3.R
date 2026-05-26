# Box plot analysis for GWAS subgroups results
# Creates box plots with data points for MAF, N, and p-values

#Setup
args <- commandArgs(trailingOnly = FALSE)
script_file <- args[grep("--file=", args)]
if (length(script_file) > 0) {
  script_dir <- dirname(sub("--file=", "", script_file[1]))
  setwd(script_dir)
}
if (!dir.exists("plots")) dir.create("plots")

#Read Data

maf_data <- read.csv("all_maf.csv", na.strings ="NA")
n_data <- read.csv("all_n.csv", na.strings ="NA")
pvals_data <- read.csv("all_pval.csv", na.strings ="NA")

#maf_data <- read.csv("capitulum_maf.csv", na.strings ="NA")
#n_data <- read.csv("capitulum_n.csv", na.strings ="NA")
#pvals_data <- read.csv("capitulum_pval.csv", na.strings ="NA")

#Transform p-values
pvals_data_log <- -log10(pvals_data[, grep("p.value", names(pvals_data))])

#  Helper for significance labels 
get_sig_label <- function(p_val) {
  if (p_val < 0.001) return("***")
  if (p_val < 0.01) return("**")
  if (p_val < 0.05) return("*")
  return("ns")
}

#  Pairwise standard t-test (pooled variance) 
pairwise_ttest <- function(data_list) {
  n_groups <- length(data_list)
  pval_matrix <- matrix(NA, nrow = n_groups, ncol = n_groups)
  
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      x <- na.omit(data_list[[i]])
      y <- na.omit(data_list[[j]])
      pval_matrix[i, j] <- t.test(x, y, var.equal = TRUE)$p.value  # use regular t-test
    }
  }
  
  return(pval_matrix)
}

#  Function to plot significance brackets 
plot_significance <- function(pvals_matrix, y_max, positions, cex_text = 1.4, spacing_factor = 1.6) {
  n_groups <- length(positions)
  n_comparisons <- n_groups * (n_groups - 1) / 2
  available_space <- y_max * 0.35  # Allow more headroom for spacing
  bracket_height <- available_space / n_comparisons * 1.2  # Increase vertical distance between brackets
  y_pos <- y_max * 0.95
  
  for (i in 1:(n_groups - 1)) {
    for (j in (i + 1):n_groups) {
      if (!is.na(pvals_matrix[i, j])) {
        p_val <- pvals_matrix[i, j]
        sig_label <- get_sig_label(p_val)
        
        x1 <- positions[i]
        x2 <- positions[j]
        y_offset <- bracket_height * 0.25
        
        segments(x1, y_pos, x1, y_pos - y_offset, lwd = 1.2)
        segments(x2, y_pos, x2, y_pos - y_offset, lwd = 1.2)
        segments(x1, y_pos, x2, y_pos, lwd = 1.2)
        text((x1 + x2) / 2, y_pos + y_offset * 0.8, sig_label,
             cex = cex_text, adj = c(0.5, 0), font = 2)
        
        y_pos <- y_pos - bracket_height * spacing_factor  # more vertical space between brackets
      }
    }
  }
}

#  Compact boxplot function with larger fonts 
plot_boxplot <- function(data_subset, title, ylab, add_line = FALSE, y_factor = 1.4, ylim_override = NULL) {
  y_max <- max(data_subset, na.rm = TRUE)
  if (is.null(ylim_override)) {
    ylim <- c(0, y_max * y_factor)
  } else {
    ylim <- ylim_override
  }
  
  bp <- boxplot(data_subset,
                main = title,
                ylab = ylab,
                xlab = "Population Group",
                ylim = ylim,
                #col = c("lightblue", "lightgreen", "lightcoral", "lightyellow"),
                names = c("Full", "S. integ.", "Wild", "Breeding"),
                pch = NA,
                cex.lab = 1.3,
                cex.axis = 1.2)
  
  # Jittered points
  for (i in 1:ncol(data_subset)) {
    x_pos <- jitter(rep(i, nrow(data_subset)), amount = 0.12)
    y_vals <- data_subset[, i]
    points(x_pos, y_vals, pch = 16, cex = 0.7, col = rgb(0, 0, 0, 0.4))
  }
  
  if (add_line) abline(h = 7, lty = 2, col = "red", lwd = 2)
  grid(nx = NA, ny = NULL, lty = "dotted")
  
  return(list(data = data_subset, max = y_max, ylim = ylim))
}

#  Pairwise tests (standard t-tests) 
maf_pvals <- pairwise_ttest(list(maf_data$MAF_full_population,
                                  maf_data$MAF_S.integrifolium_genets,
                                  maf_data$MAF_wild_genets,
                                  maf_data$MAF_breeding_genets))

n_pvals <- pairwise_ttest(list(n_data$N_full_population,
                                n_data$N_S.integrifolium_genets,
                                n_data$N_wild_genets,
                                n_data$N_breeding_genets))

pvals_log_pvals <- pairwise_ttest(list(pvals_data_log$p.value_full_population,
                                        pvals_data_log$p.value_S.integrifolium_genets,
                                        pvals_data_log$p.value_wild_genets,
                                        pvals_data_log$p.value_breeding_genets))

#  Data selection 
maf_cols <- maf_data[, c("MAF_full_population", "MAF_S.integrifolium_genets",
                         "MAF_wild_genets", "MAF_breeding_genets")]
n_cols <- n_data[, c("N_full_population", "N_S.integrifolium_genets",
                     "N_wild_genets", "N_breeding_genets")]
pvals_cols <- pvals_data_log[, c("p.value_full_population", "p.value_S.integrifolium_genets",
                                 "p.value_wild_genets", "p.value_breeding_genets")]

#  Create updated plots 
#png("plots/maf_boxplot.png", width = 600, height = 750, res = 120)
png("plots/all_maf_boxplot.png", width = 600, height = 750, res = 120)
result <- plot_boxplot(maf_cols, "Minor Allele Frequency (MAF)", "MAF")
plot_significance(maf_pvals, result$max * 1.2, c(1, 2, 3, 4), cex_text = 1.4, spacing_factor = 1)
dev.off()

#png("plots/n_boxplot.png", width = 600, height = 750, res = 120)
png("plots/all_n_boxplot.png", width = 600, height = 750, res = 120)
result <- plot_boxplot(n_cols, "Sample Size (N)", "Number of Samples")
plot_significance(n_pvals, result$max * 1.4, c(1, 2, 3, 4), cex_text = 1.4, spacing_factor = 1)
dev.off()

pval_ymax <- max(pvals_cols, na.rm = TRUE)
pval_ylim <- c(0, min(pval_ymax * 1.2, 40))  # cap at 12 for compact plots, adjust as needed

#png("plots/pvals_boxplot.png", width = 600, height = 750, res = 120)
png("plots/all_pvals_boxplot.png", width = 600, height = 750, res = 120)
result <- plot_boxplot(pvals_cols,
                       "P-values (log10 scale)",
                       "-log10(p-value)",
                       add_line = TRUE,
                       ylim_override = pval_ylim)
plot_significance(pvals_log_pvals, result$max * 1.1, c(1, 2, 3, 4))
dev.off()