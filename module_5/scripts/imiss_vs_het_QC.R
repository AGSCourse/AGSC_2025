### Load PLINK .imiss and .het output files
imiss <- read.table("raw-GWA-data.imiss", header = TRUE)
het   <- read.table("raw-GWA-data.het",   header = TRUE)

### Merge on individual IDs
dat <- merge(imiss, het, by = c("FID", "IID"))

### Calculate missingness rate per individual
dat$miss_rate <- dat$N_MISS / dat$N_GENO

### Calculate observed heterozygosity rate
dat$het_rate <- (dat$N.NM. - dat$O.HOM.) / dat$N.NM.

### Determine heterozygosity thresholds (mean ± 3 SD)
mean_het <- mean(dat$het_rate, na.rm = TRUE)
sd_het   <- sd(dat$het_rate, na.rm = TRUE)
low_het  <- mean_het - 3 * sd_het
high_het <- mean_het + 3 * sd_het

### Missingness threshold (from instructions)
miss_thresh <- 0.03

### Identify individuals failing imiss/het QC
fail <- subset(
  dat,
  miss_rate >= miss_thresh | het_rate < low_het | het_rate > high_het
)

### Print number of failing individuals
cat("Number of failing individuals:", nrow(fail), "\n")

### Write FID IID list for PLINK
write.table(
  fail[, c("FID", "IID")],
  file = "fail-imisshet-qc.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("Output written to fail-imisshet-qc.txt\n")
