# Analysis of IOVS mouse lens data (Supplemental Table S01):
# Khan, Shahid Y., et al. "Proteome Profiling of Developing Murine Lens Through Mass Spectrometry." 
# Investigative Ophthalmology & Visual Science 59.1 (2018): 100-107.

# load libraries
library(tidyverse)
library(limma)
library(edgeR)
library(sva)

# read the Supplemental 01 file (saved as a CSV export from XLSX file)
data_start <- read_csv("iovs-58-13-55_s01.csv")

# filter out proteins not seen in all three runs
data_no_na <- na.omit(data_start)

# fix the column headers
col_headers <- colnames(data_no_na)
col_headers <- str_replace(col_headers, " {2,3}", " ")
col_headers <- str_replace(col_headers, "Reporter ion intensities ", "")
colnames(data_no_na) <- col_headers

# save the annotation columns (gene symbol and protein accession) for later and remove from data frame
annotate_df <- data_no_na[1:2]
data_raw <- data_no_na[3:20]
row.names(data_raw) <- annotate_df$`Protein Accession No.`

# separate the TMT data by experiment
exp1_raw <- data_raw[c(1:6)]
exp2_raw <- data_raw[c(7:12)]
exp3_raw <- data_raw[c(13:18)]

# figure out the global scaling value
target <- mean(c(colSums(exp1_raw), colSums(exp2_raw), colSums(exp3_raw)))

# do the sample loading normalization before the IRS normalization
# there is a different correction factor for each column
# seems like a loop could be used here somehow...
norm_facs <- target / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3_raw)
exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")

# make a pre-IRS data frame after sample loading norms
data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl)

# make working frame with row sums from each frame
irs <- tibble(rowSums(exp1_sl), rowSums(exp2_sl), rowSums(exp3_sl))
colnames(irs) <- c("sum1", "sum2", "sum3")

# get the geometric average intensity for each protein
irs$average <- apply(irs, 1, function(x) exp(mean(log(x))))

# compute the scaling factor vectors
irs$fac1 <- irs$average / irs$sum1
irs$fac2 <- irs$average / irs$sum2
irs$fac3 <- irs$average / irs$sum3

# make new data frame with normalized data
data_irs <- exp1_sl * irs$fac1
data_irs <- cbind(data_irs, exp2_sl * irs$fac2)
data_irs <- cbind(data_irs, exp3_sl * irs$fac3)

# some experimental design setup
mygroups <- c(rep.int(c(1, 2, 3, 4, 5, 6), 3))
design <- read.csv("design.csv")
mygroups <- design$group
mod <- model.matrix(~ group, data = design)
batch <- design$batch

# run ComBat as alernative to IRS
# NOTE: sample loading is probably better to do before feeding the batch corection
data_combat <- ComBat(dat = data_sl, batch = batch, mod = mod, par.prior = TRUE, mean.only = TRUE)
par(mfrow = c(1, 1))

# explore the different normalizations with boxplots
# we have raw data (un-norm), sample loading, IRS, library size, and TMM
col_vec1 <- c(rep_len("red", 18), rep_len("green", 18))
boxplot(cbind(log2(data_raw), log2(data_sl)), main = "Raw data (red) vs SL data (green)", col = col_vec1)
# # cpm function does the library size scaling (the same relative scalings as SL)
# boxplot(cbind(log2(data_sl), cpm(data_raw, log = TRUE)), main = "SL data (red) vs cpm() (green)\noriginal scale 10,000 times larger than cpm", col = col_vec1)
# # cpm divides column by column sum scalar then multiplies by 1e+06 scalar
# # cpm scales to library size of 1,000,000 in each sample, this resets any prior multiplicative factors
# # adjust SL data to one million per sample
# temp <- (1e+06/target) * data_sl
# boxplot(cbind(log2(temp), cpm(data_sl, log = TRUE)), main = "scaled SL data (red) vs cpm() (green)", col = col_vec1)
# # double check the data are the same
# summary(temp[1:5])
# summary(cpm(data_sl)[, 1:5])
# finally, does IRS norm appear different with this visualization?
boxplot(cbind(log2(data_sl), log2(data_irs)), main = "SL data (red) vs IRS data (green)", col = col_vec1)
boxplot(cbind(log2(data_sl), log2(data_combat)), main = "SL data (red) vs ComBat data (green)", col = col_vec1)
boxplot(cbind(log2(data_irs), log2(data_combat)), main = "IRS data (red) vs ComBat data (green)", col = col_vec1)
# interestingly, the cpm library scaling effectively removes any prior scaling
# conclusions: sample loading and library size corrections are similar (overall scale differs), 
#   cpm function log option is a log2, cpm() of either Raw or SL data will be the same

# see what TMM does
# raw data
raw_tmm <- calcNormFactors(data_raw)
data_raw_tmm <- sweep(data_raw, 2, raw_tmm, FUN = "/") # this is raw data after TMM on original scale
boxplot(cbind(log2(data_raw), log2(data_raw_tmm)), main = "Raw data (red) vs TMM of raw (green)", col = col_vec1)
# SL data
sl_tmm <- calcNormFactors(data_sl)
data_sl_tmm <- sweep(data_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
boxplot(cbind(log2(data_sl), log2(data_sl_tmm)), main = "SL data (red) vs TMM of SL (green)", col = col_vec1)
# IRS data
irs_tmm <- calcNormFactors(data_irs)
data_irs_tmm <- sweep(data_irs, 2, irs_tmm, FUN = "/") # this is data after SL, IRS, and TMM on original scale
boxplot(cbind(log2(data_irs), log2(data_irs_tmm)), main = "IRS data (red) vs TMM of IRS (green)", col = col_vec1)
# ComBat data
combat_tmm <- calcNormFactors(data_combat)
data_combat_tmm <- sweep(data_combat, 2, combat_tmm, FUN = "/") # this is data after SL, IRS, and TMM on original scale
boxplot(cbind(log2(data_combat), log2(data_combat_tmm)), main = "ComBat data (red) vs TMM of ComBat (green)", col = col_vec1)

# density plots are also informative (compare w/o and w/ TMM)
par(mfrow = c(2, 1))
plotDensities(log2(data_raw), main = "raw data")
plotDensities(log2(data_raw_tmm), main = "Raw data after TMM scaling factors")
plotDensities(log2(data_sl), main = "SL data")
plotDensities(log2(data_sl_tmm), main = "SL data after TMM scaling factors")
plotDensities(log2(data_irs), main = "IRS data")
plotDensities(log2(data_irs_tmm), main = "IRS data after TMM scaling factors")
plotDensities(log2(data_combat), main = "ComBat data (x-axis scale may differ)")
plotDensities(log2(data_combat_tmm), main = "ComBat data After TMM scaling factors")
par(mfrow = c(1, 1))

# look at some clustering
# Correct outcome is that data should group by developmental time points
# Grouping by TMT experiment means that correction of measurement factors has failed 
col_vec2 <- c("red", "orange", "black", "green", "blue", "violet")
plotMDS(log2(data_raw), col = col_vec2, main = 'Raw data (log2)') # groups by TMT exp, not as tight
plotMDS(cpm(data_raw, log = TRUE), col = col_vec2, main = 'cpm(w/ log) of raw data') # cpm does library size correction
plotMDS(DGEList(counts = data_raw, group = mygroups), col = col_vec2, main = 'DGEList of raw data (same as cpm())')
plotMDS(log2(data_sl), col = col_vec2, main = 'SL data (log2) [SL and library scaling are similar]') # SL just like library correction
plotMDS(DGEList(counts = data_sl, group = mygroups), col = col_vec2, main = 'DGEList of SL data') # changes scale not differences
plotMDS(y_sl_tmm, col = col_vec2, main = 'DGEList of TMM-normalized SL data\n(best non-IRS treatment)')
# explore deeper dimensions to see when biology pops out
plotMDS(y_sl_tmm, dim.plot = c(2, 3), col = col_vec2, main = 'Vertical is still by TMT Exp., Horizontal is by time')
plotMDS(y_sl_tmm, dim.plot = c(3, 4), col = col_vec2, main = 'TMT Exp batches are finally gone')
# see how well IRS is at batch removal
plotMDS(DGEList(counts = data_irs, group = mygroups), col = col_vec2, main = 'DGEList of IRS data')
plotMDS(y_irs_tmm, col = col_vec2, main = 'DGEList of TMM-normalized IRS data\nData finally groups by biological factors')
# plotMDS, when given a DGEList, looks like a cpm transformation is done



# compare E15+E18 to P6+P9 with edgeR and limma
# before IRS data (SL normed only)
group = as.factor(rep(c("early", "early", "middle", "middle", "late", "late"), 3))
group <- factor(group, levels(group)[c(1, 3, 2)]) # set the factor order
tmt <- as.factor(rep(c("exp1", "exp2", "exp3"), c(6, 6, 6))) # possible batch factors
yb <- DGEList(counts = SLNorm, group = group)
yb$samples$tmt <- tmt
row.names(yb) <- annotate_df$`Protein Accession No.`
cpm_yb <- cpm(yb)
lcpm_yb <- cpm(yb, log = TRUE)





# previous code *********************************************
yb <- calcNormFactors(yb)
yb <- estimateDisp(yb)
plotBCV(yb, main = "Before IRS")
etb <- exactTest(yb, pair = c(1, 3))
summary(decideTestsDGE(etb))
ttb <- topTags(etb, n = 5000, sort.by = "none")
ttb <- ttb$table
ttb$candidate <- "no"
ttb[which(ttb$FDR <= 0.10 & ttb$FDR > 0.05), dim(ttb)[2]] <- "low"
ttb[which(ttb$FDR <= 0.05 & ttb$FDR > 0.01), dim(ttb)[2]] <- "med"
ttb[which(ttb$FDR <= 0.01), dim(ttb)[2]] <- "high"
ttb$candidate <- as.factor(ttb$candidate)
DE_before <- data.frame(rowMeans(SLNorm[c(1, 2, 7, 8, 13, 14)]), rowMeans(SLNorm[c(5, 6, 11, 12, 17, 18)]), ttb$candidate)
colnames(DE_before) <- c("Embryo", "Late", "candidate")
ggplot(ttb, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(ttb$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("SLNorm only p-value distribution (about 300 candidates)")

# after IRS
ya <- DGEList(counts = IRSNorm, group = group)
ya <- calcNormFactors(ya)
ya <- estimateDisp(ya)
plotBCV(ya, main = "After IRS")
eta <- exactTest(ya, pair = c(1, 3))
summary(decideTestsDGE(eta))
tta <- topTags(eta, n = 5000, sort.by = "none")
tta <- tta$table
tta$candidate <- "no"
tta[which(tta$FDR <= 0.10 & tta$FDR > 0.05), dim(tta)[2]] <- "low"
tta[which(tta$FDR <= 0.05 & tta$FDR > 0.01), dim(tta)[2]] <- "med"
tta[which(tta$FDR <= 0.01), dim(tta)[2]] <- "high"
tta$candidate <- as.factor(tta$candidate)
DE_after <- data.frame(rowMeans(IRSNorm[c(1, 2, 7, 8, 13, 14)]), rowMeans(IRSNorm[c(5, 6, 11, 12, 17, 18)]), tta$candidate)
colnames(DE_after) <- c("Embryo", "Late", "candidate")
ggplot(tta, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tta$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("IRSNorm p-value distribution (about 1800 candidates)")

# make the combined candidate corelation plot
ggplot(DE_before, aes(x = Embryo, y = Late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Before IRS embryo vs late") + 
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down

ggplot(DE_after, aes(x = Embryo, y = Late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("After IRS, embryo vs late") + 
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate corelation plots
ggplot(DE_before, aes(x = Embryo, y = Late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("Before, separated by candidate")

ggplot(DE_after, aes(x = Embryo, y = Late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("After, separated by candidate")
