# Analysis of IOVS mouse lens data (Supplemental Table S01):
# Khan, Shahid Y., et al. "Proteome Profiling of Developing Murine Lens Through Mass Spectrometry." 
# Investigative Ophthalmology & Visual Science 59.1 (2018): 100-107.

# load libraries
library(tidyverse)
library(limma)
library(edgeR)
library(sva)
library(reshape2)

# create some functions
subset_proteins <- function(df, proteins, labels, value_col) {
  df_subset <- subset(df, df$protein %in% proteins)
  df_subset$proteins <- labels
  keycol <- "time"
  valuecol <- value_col
  gathercols <- colnames(df)[1:6]
  return(gather_(df_subset, keycol, valuecol, gathercols))
}

plot_proteins <- function(df_list, annotate_df) {
  ave_df <- df_list[[1]]
  sd_df <- df_list[[2]]
  
  # lists of lens proteins to extract and plot
  alphas <- c("P24622", "P23927")
  alpha_labels <- c("alphaA", "alphaB")
  betas <- c("P02525", "Q9JJV1", "Q9JJV0","Q9WVJ5",	"P62696", "Q9JJU9", "O35486")	
  beta_labels <- c("betaA1", "betaA2", "betaA4", "betaB1", "betaB2", "betaB3", "betaS")
  gammas <- c("P04345", "P04344", "Q61597", "P04342", "Q03740", "Q9CXV3", "Q8VHL5")
  gamma_labels <- c("gammaA", "gammaB", "gammaC", "gammaD", "gammaE", "gammaF", "gammaN")
  ave_df$protein <- annotate_df$`Protein Accession No.`
  sd_df$protein <- annotate_df$`Protein Accession No.`
  # plot alphas
  ave_df_wide <- subset_proteins(ave_df, alphas, alpha_labels, "intensity")
  sd_df_wide <- subset_proteins(sd_df, alphas, alpha_labels, "sd")
  ave_df_wide$sd <- sd_df_wide$sd
  alpha_plot <- ggplot(ave_df_wide, aes(time, intensity, fill=proteins)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin = intensity-sd, ymax = intensity+sd), position = position_dodge(0.9), width = 0.4) +
    xlab("Developmental Time") + ylab("Total Protein Intensity") +
    ggtitle("Alpha crystallin expression")
  print(alpha_plot)
  # plot betas
  ave_df_wide <- subset_proteins(ave_df, betas, beta_labels, "intensity")
  sd_df_wide <- subset_proteins(sd_df, betas, beta_labels, "sd")
  ave_df_wide$sd <- sd_df_wide$sd
  beta_plot <- ggplot(ave_df_wide, aes(time, intensity, fill=proteins)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin = intensity-sd, ymax = intensity+sd), position = position_dodge(0.9), width = 0.4) +
    xlab("Developmental Time") + ylab("Total Protein Intensity") +
  ggtitle("Beta crystallin expression") 
  print(beta_plot)
  # plot gammas
  ave_df_wide <- subset_proteins(ave_df, gammas, gamma_labels, "intensity")
  sd_df_wide <- subset_proteins(sd_df, gammas, gamma_labels, "sd")
  ave_df_wide$sd <- sd_df_wide$sd
  gamma_plot <- ggplot(ave_df_wide, aes(time, intensity, fill=proteins)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin = intensity-sd, ymax = intensity+sd), position = position_dodge(0.9), width = 0.4) +
    xlab("Developmental Time") + ylab("Total Protein Intensity") +
  ggtitle("Gamma crystallin expression")
  print(gamma_plot)
  return(NULL)
}

# computes CVs per time point
make_CVs <- function(df) {
  # separate by time points
  E15 <- df[c(1, 7, 13)]
  E18 <- df[c(2, 8, 14)]
  P0 <- df[c(3, 9, 15)]
  P3 <- df[c(4, 10, 16)]
  P6 <- df[c(5, 11, 17)]
  P9 <- df[c(6, 12, 18)]
  
  E15$ave <- rowMeans(E15)
  E15$sd <- apply(E15[1:3], 1, sd)
  E15$cv <- 100 * E15$sd / E15$ave
  E18$ave <- rowMeans(E18)
  E18$sd <- apply(E18[1:3], 1, sd)
  E18$cv <- 100 * E18$sd / E18$ave
  P0$ave <- rowMeans(P0)
  P0$sd <- apply(P0[1:3], 1, sd)
  P0$cv <- 100 * P0$sd / P0$ave
  P3$ave <- rowMeans(P3)
  P3$sd <- apply(P3[1:3], 1, sd)
  P3$cv <- 100 * P3$sd / P3$ave
  P6$ave <- rowMeans(P6)
  P6$sd <- apply(P6[1:3], 1, sd)
  P6$cv <- 100 * P6$sd / P6$ave
  P9$ave <- rowMeans(P9)
  P9$sd <- apply(P9[1:3], 1, sd)
  P9$cv <- 100 * P9$sd / P9$ave
  
  ave_df <- data.frame(E15$ave, E18$ave, P0$ave, P3$ave, P6$ave, P9$ave)
  sd_df <- data.frame(E15$sd, E18$sd, P0$sd, P3$sd, P6$sd, P9$sd)
  cv_df <- data.frame(E15$cv, E18$cv, P0$cv, P3$cv, P6$cv, P9$cv)
  return(list(ave_df, sd_df, cv_df))
}

# kind of hardcoded test where data frame can change
run_edgeR <- function(df, group, comparison, norm) {
  df <- df[apply(df, 1, function(x) all(x>=0)), ] 
  y <- DGEList(counts = df, group = group)
  if(norm) y <- calcNormFactors(y)
  df <- sweep(df, 2, y$samples$norm.factors, FUN = "/") # scale by TMM factors (if any)
  y <- estimateDisp(y)
  plotBCV(y, main = comparison)
  et <- exactTest(y, pair = c("early", "late"))
  print(summary(decideTestsDGE(et)))
  tt <- topTags(et, n = 5000, sort.by = "none")
  tt <- tt$table
  tt$candidate <- "no"
  tt[which(tt$FDR <= 0.10 & tt$FDR > 0.05), dim(tt)[2]] <- "low"
  tt[which(tt$FDR <= 0.05 & tt$FDR > 0.01), dim(tt)[2]] <- "med"
  tt[which(tt$FDR <= 0.01), dim(tt)[2]] <- "high"
  tt$candidate <- as.factor(tt$candidate)
  DE <- data.frame(rowMeans(df[c(1, 2, 7, 8, 13, 14)]), rowMeans(df[c(5, 6, 11, 12, 17, 18)]), tt$candidate)
  colnames(DE) <- c("Embryo", "Late", "candidate")
  ggplot(tt, aes(PValue)) + 
    geom_histogram(bins = 100, fill = "white", color = "black") + 
    geom_hline(yintercept = mean(hist(tt$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
    ggtitle("p-value distribution")
  return(DE)
}

make_candidate_plot <- function(df, title) {
  ggplot(df, aes(x = Embryo, y = Late)) +
    geom_point(aes(color = candidate, shape = candidate)) +
    scale_y_log10() +
    scale_x_log10() +
    ggtitle(title) + 
    geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
    geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down
}

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
design <- read.csv("design.csv")
design$group <- factor(design$group, levels = c("early", "middle", "late"))
mod <- model.matrix(~ group, data = design)
batch <- design$batch

# run ComBat as alernative to IRS
# NOTE: sample loading is probably better to do before feeding the batch corection
#data_combat <- ComBat(dat = data_sl, batch = batch, mod = mod, par.prior = TRUE, mean.only = TRUE)
data_combat <- ComBat(dat = data_sl, batch = batch, mod = mod, par.prior = TRUE)
par(mfrow = c(1, 1)) # any plotting in the ComBat call leaves plots as 2x2

# explore the different normalizations with boxplots
# we have raw data (un-norm), sample loading, IRS, library size, and TMM
col_vec1 <- c(rep_len("red", 18), rep_len("green", 18))
boxplot(cbind(log2(data_raw), log2(data_sl)), notch = TRUE, main = "Raw data (red) vs SL data (green)", col = col_vec1)
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
boxplot(cbind(log2(data_sl), log2(data_irs)), notch = TRUE, main = "SL data (red) vs IRS data (green)", col = col_vec1)
boxplot(cbind(log2(data_sl), log2(data_combat)), notch = TRUE, main = "SL data (red) vs ComBat data (green)", col = col_vec1)
boxplot(cbind(log2(data_irs), log2(data_combat)), notch = TRUE, main = "IRS data (red) vs ComBat data (green)", col = col_vec1)
# interestingly, the cpm library scaling effectively removes any prior scaling
# conclusions: sample loading and library size corrections are similar (overall scale differs), 
#   cpm function log option is a log2, cpm() of either Raw or SL data will be the same

# see exactly what TMM does
# raw data
raw_tmm <- calcNormFactors(data_raw)
data_raw_tmm <- sweep(data_raw, 2, raw_tmm, FUN = "/") # this is raw data after TMM on original scale
boxplot(cbind(log2(data_raw), log2(data_raw_tmm)), notch = TRUE, main = "Raw data (red) vs TMM of raw (green)", col = col_vec1)
# SL data
sl_tmm <- calcNormFactors(data_sl)
data_sl_tmm <- sweep(data_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale
boxplot(cbind(log2(data_sl), log2(data_sl_tmm)), notch = TRUE, main = "SL data (red) vs TMM of SL (green)", col = col_vec1)
# IRS data
irs_tmm <- calcNormFactors(data_irs)
data_irs_tmm <- sweep(data_irs, 2, irs_tmm, FUN = "/") # this is data after SL, IRS, and TMM on original scale
boxplot(cbind(log2(data_irs), log2(data_irs_tmm)), notch = TRUE, main = "IRS data (red) vs TMM of IRS (green)", col = col_vec1)
# ComBat data
combat_tmm <- calcNormFactors(data_combat)
data_combat_tmm <- sweep(data_combat, 2, combat_tmm, FUN = "/") # this is data after SL, IRS, and TMM on original scale
boxplot(cbind(log2(data_combat), log2(data_combat_tmm)), notch = TRUE, main = "ComBat data (red) vs TMM of ComBat (green)", col = col_vec1)

# density plots are also informative (compare w/o and w/ TMM)
plotDensities(log2(data_raw), main = "raw data")
plotDensities(log2(data_raw_tmm), main = "Raw data after TMM scaling factors")
plotDensities(log2(data_sl), main = "SL data")
plotDensities(log2(data_sl_tmm), main = "SL data after TMM scaling factors")
plotDensities(log2(data_irs), main = "IRS data")
plotDensities(log2(data_irs_tmm), main = "IRS data after TMM scaling factors")
plotDensities(log2(data_combat), main = "ComBat data (x-axis scale may differ)")
plotDensities(log2(data_combat_tmm), main = "ComBat data After TMM scaling factors")

# look at some clustering
# Correct outcome is that data should group by developmental time points
# Grouping by TMT experiment means that there are effectively "batch" effects
col_vec2 <- c("red", "orange", "black", "green", "blue", "violet")
plotMDS(log2(data_raw), col = col_vec2, main = 'Raw data (log2)') # groups by TMT exp, not as tight per cluster
plotMDS(log2(data_sl), col = col_vec2, main = 'SL data (log2)') # SL just like library correction
plotMDS(log2(data_sl_tmm), col = col_vec2, main = 'TMM-normalized SL data\n(best non-IRS treatment)')
# explore deeper dimensions to see when biology pops out
plotMDS(log2(data_sl_tmm), dim.plot = c(2, 3), col = col_vec2, 
        main = 'SL/TMM: Vertical is still by TMT Exp., Horizontal is by time\n(dimensions 2 and 3)')
plotMDS(log2(data_sl_tmm), dim.plot = c(3, 4), col = col_vec2, 
        main = 'SL/TMM: TMT Exp batches are finally gone\n(dimensions 3 and 4)')
# see how well is IRS at batch removal?
plotMDS(log2(data_irs_tmm), col = col_vec2, main = 'TMM-normalized IRS data\nData finally groups by biological factors')
plotMDS(log2(data_combat_tmm), col = col_vec2, main = 'TMM-normalized Combat data')
# plotMDS, when given a DGEList, looks like a cpm transformation is done

# get CVs and averages
list_raw <- make_CVs(data_raw_tmm)
list_sl <- make_CVs(data_sl_tmm)
list_irs <- make_CVs(data_irs_tmm)
list_combat <- make_CVs(data_combat_tmm)

# compare CV distributions
par(mfrow = c(2, 2))
boxplot(list_raw[[3]], notch = TRUE, main = "RAW/TMM CVs", ylim = c(0, 200))
boxplot(list_sl[[3]], notch = TRUE, main = "SL/TMM CVs", ylim = c(0, 200))
boxplot(list_irs[[3]], notch = TRUE, main = "IRS/TMM CVs", ylim = c(0, 150))
boxplot(list_combat[[3]], notch = TRUE, main = "ComBat/TMM CVs", ylim = c(0, 150))
par(mfrow = c(1, 1))

(raw_med_cv <- mean(apply(list_raw[[3]], 2, median)))
(sl_med_cv <- mean(apply(list_sl[[3]], 2, median)))
(irs_med_cv <- mean(apply(list_irs[[3]], 2, median)))
(combat_med_cv <- mean(apply(list_combat[[3]], 2, median)))

# do a couple of multi-panel to see effects of batch correction
par(mfrow = c(2,2))
boxplot(list_raw[[3]], notch = TRUE, main = "RAW/TMM CVs (Ave = 62%)")
boxplot(list_sl[[3]], notch = TRUE, main = "SL/TMM CVs (Ave = 55%)")
plotMDS(log2(data_raw_tmm), col = col_vec2, main = 'Raw data with TMM')
plotMDS(log2(data_sl_tmm), col = col_vec2, main = 'TMM-normalized SL data\n(best non-IRS treatment)')

boxplot(list_irs[[3]], notch = TRUE, main = "IRS/TMM CVs (Ave = 13%)")
boxplot(list_combat[[3]], notch = TRUE, main = "ComBat/TMM CVs (Ave = 18%)", ylim = c(0, 150))
plotMDS(log2(data_irs_tmm), col = col_vec2, main = "TMM-normalized IRS data")
plotMDS(log2(data_combat_tmm), col = col_vec2, main = "TMM-normalized Combat data")
par(mfrow = c(1, 1))

# run several edgeR analyses
de_raw <- run_edgeR(data_raw, design$group, "RAW data, no TMM", FALSE)
de_raw_tmm <- run_edgeR(data_raw, design$group, "RAW data after TMM", TRUE)
de_sl <- run_edgeR(data_sl, design$group, "SL data, no TMM", FALSE)
de_sl_tmm <- run_edgeR(data_sl, design$group, "SL data after TMM", TRUE)
de_irs <- run_edgeR(data_irs, design$group, "IRS data, no TMM", FALSE)
de_irs_tmm <- run_edgeR(data_irs, design$group, "IRS data after TMM", TRUE)
# seems ComBat puts in some negative numbers that edgeR does not like
de_combat <- run_edgeR(data_combat, design$group, "ComBat data, no TMM", FALSE)
de_combat_tmm <- run_edgeR(data_combat, design$group, "ComBat data after TMM", TRUE)

# make the combined candidate corelation plot
make_candidate_plot(de_raw, "RAW data, no TMM")
make_candidate_plot(de_raw_tmm, "RAW data after TMM")
make_candidate_plot(de_sl, "SL data, no TMM")
make_candidate_plot(de_sl_tmm, "SL data after TMM")
make_candidate_plot(de_irs, "IRS data, no TMM")
make_candidate_plot(de_irs_tmm, "IRS data after TMM")
make_candidate_plot(de_combat, "ComBat data, no TMM")
make_candidate_plot(de_combat_tmm, "ComBat data after TMM")

# make separate corelation plots
ggplot(de_irs_tmm, aes(x = Embryo, y = Late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("IRS w/TMM, separated by candidate")

ggplot(de_combat_tmm, aes(x = Embryo, y = Late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("ComBat w/TMM, separated by candidate")
