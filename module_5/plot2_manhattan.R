# plot_manhattan_qqman.R
# Manhattan plot from PLINK logistic output using qqman

## ---------- 1. Read GWAS results ----------

result_file <- "additive.AGE.SEX.analysis.assoc.logistic"  # change if needed

gwas <- read.table(result_file, header = TRUE)

# Keep only additive test, autosomes, non-missing P
gwas_add <- subset(gwas,
                   TEST == "ADD" &
                   !is.na(P) &
                   CHR >= 1 & CHR <= 22)

if (nrow(gwas_add) == 0) {
  stop("No ADD rows found in ", result_file)
}

## ---------- 2. Ensure qqman is available ----------

if (!requireNamespace("qqman", quietly = TRUE)) {
  message("qqman not found, trying to install from CRAN...")
  tryCatch(
    {
      install.packages("qqman", repos = "https://cloud.r-project.org")
    },
    error = function(e) {
      message("CRAN install failed, trying GitHub (stephenturner/qqman)...")
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes", repos = "https://cloud.r-project.org")
      }
      remotes::install_github("stephenturner/qqman")
    }
  )
}

library(qqman)

## ---------- 3. Manhattan plot ----------

pdf("additive.AGE.SEX.manhattan.qqman.pdf")

manhattan(gwas_add,
          chr = "CHR",
          bp  = "BP",
          snp = "SNP",
          p   = "P",
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5),
          main = "Additive model (AGE + SEX adjusted)")

dev.off()

cat("Manhattan plot written to additive.AGE.SEX.manhattan.qqman.pdf\n")

