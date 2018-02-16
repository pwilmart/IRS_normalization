# Analysis of IOVS mouse lens data (Supplemental Table S01):
# Khan, Shahid Y., et al. "Proteome Profiling of Developing Murine Lens Through Mass Spectrometry."
# Investigative Ophthalmology & Visual Science 59.1 (2018): 100-107.
# load libraries
library(tidyverse) 
library(limma) 
library(edgeR) 
library(sva)

#read the Supplemental 01 file (saved as a CSV export from XLSX file)
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

# check what the final data frame looks like
head(data_raw)

# let's see what the starting data look like
boxplot(log2(data_raw), col = rep(c("red", "green", "blue"), each = 6), 
        notch = TRUE, main = "RAW data: Exp1 (red), Exp2 (green), Exp3 (blue)")

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_raw), col = rep(c("red", "green", "blue"), 6), main = "Raw data")

# check column totals to see how similar they are
format(round(colSums(data_raw), digits = 0), big.mark = ",")

# we will do some normalization within each TMT experiment, so separate the data
# separate the TMT data by experiment
exp1_raw <- data_raw[c(1:6)]
exp2_raw <- data_raw[c(7:12)]
exp3_raw <- data_raw[c(13:18)]

# first basic normalization is to adjust each TMT experiment to equal signal per channel

# figure out the global scaling value
target <- mean(c(colSums(exp1_raw), colSums(exp2_raw), colSums(exp3_raw)))

# do the sample loading (SL) normalization before the IRS normalization
# there is a different correction factor for each column
norm_facs <- target / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3_raw)
exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")

# make a pre-IRS data frame after sample loading norms
data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl)

# see what the SL normalized data look like
boxplot(log2(data_sl), col = rep(c("red", "green", "blue"), each = 6), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nExp1 (red), Exp2 (green), Exp3 (blue)")

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_sl), col = rep(c("red", "green", "blue"), 6), main = "SL normalization")

# check the columnn totals
format(round(colSums(data_sl), digits = 0), big.mark = ",")

# do TMM on raw data - we use it later for CV analysis
raw_tmm <- calcNormFactors(data_raw)
data_raw_tmm <- sweep(data_raw, 2, raw_tmm, FUN = "/") # this is data after SL and TMM on original scale

# see exactly what TMM does with SL data
sl_tmm <- calcNormFactors(data_sl)
data_sl_tmm <- sweep(data_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale

boxplot(log2(data_sl_tmm), notch = TRUE, col = rep(c("red", "green", "blue"), each = 6), 
        main = "TMM normalization of SL data\nExp1 (red), Exp2 (green), Exp3 (blue)")

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_sl_tmm), col = rep(c("red", "green", "blue"), 6), main = "SL/TMM normalization")

# check column totals after TMM
format(round(colSums(data_sl_tmm), digits = 0), big.mark = ",")

# see how things cluster after we have gotten the boxplots and desity plots looking nice
plotMDS(log2(data_sl_tmm), col = rep(c("red", "green", "blue"), each = 6), 
        main = "SL/TMM clusters by TMT exeriment")

# data is not clustering correctly, we need to do Internal Reference Scaling (IRS)
# make working frame with row sums from each frame
irs <- tibble(rowSums(exp1_sl), rowSums(exp2_sl), rowSums(exp3_sl))
colnames(irs) <- c("sum1", "sum2", "sum3")

# get the geometric average intensity for each protein (arithmetic average works fine, too)
irs$average <- apply(irs, 1, function(x) exp(mean(log(x))))

# compute the scaling factor vectors
irs$fac1 <- irs$average / irs$sum1
irs$fac2 <- irs$average / irs$sum2
irs$fac3 <- irs$average / irs$sum3

# make new data frame with normalized data
data_irs <- exp1_sl * irs$fac1
data_irs <- cbind(data_irs, exp2_sl * irs$fac2)
data_irs <- cbind(data_irs, exp3_sl * irs$fac3)

# see what the data look like
boxplot(log2(data_irs), col = rep(c("red", "green", "blue"), each = 6), notch = TRUE, 
        main = "Internal Reference Scaling (IRS) normalized data: \nExp1 (red), Exp2 (green), Exp3 (blue)")

# can also look at density plots (like a distribution histogram)    
plotDensities(log2(data_irs), col = rep(c("red", "green", "blue"), 6), main = "IRS data")

# see what column totals look like
format(round(colSums(data_irs), digits = 0), big.mark = ",")

# we can add TMM normalization after IRS
irs_tmm <- calcNormFactors(data_irs)
data_irs_tmm <- sweep(data_irs, 2, irs_tmm, FUN = "/") # this is data after SL, IRS, and TMM on original scale

boxplot(log2(data_irs_tmm), notch = TRUE, col = rep(c("red", "green", "blue"), each = 6), 
        main = "TMM normalization of IRS data\nExp1 (red), Exp2 (green), Exp3 (blue)")

# can also look at density plots (like a distribution histogram)
plotDensities(log2(data_irs_tmm), col = rep(c("red", "green", "blue"), 6), main = "IRS/TMM data")

# and check column totals
format(round(colSums(data_irs_tmm), digits = 0), big.mark = ",")

# see how things cluster after IRS and TMM
col_vec <- c("red", "orange", "black", "green", "blue", "violet")
plotMDS(log2(data_irs_tmm), col = col_vec, main = "IRS/TMM clusters are grouped by time")

# is IRS just another "batch" correction?
# some experimental design setup
design <- read.csv("design.csv")
design$group <- factor(design$group, levels = c("early", "middle", "late"))
mod <- model.matrix(~ group, data = design)
batch <- design$batch

# run ComBat as alternative to IRS
# NOTE: sample loading is probably better to do before the batch corection
data_combat <- ComBat(dat = data_sl, batch = batch, mod = mod, par.prior = TRUE)
par(mfrow = c(1, 1)) # any plotting in the ComBat call leaves plots as 2x2

# ComBat introduces some negative corrected counts
data_combat <-  data_combat[apply(data_combat, 1, function(x) all(x > 0)), ] 

# usual views
boxplot(log2(data_combat), notch = TRUE, col = rep(c("red", "green", "blue"), each = 6), 
        main = "ComBat batch correction of SL data\nExp1 (red), Exp2 (green), Exp3 (blue)")

plotDensities(log2(data_combat), col = rep(c("red", "green", "blue"), 6), main = "ComBat data")

# check column totals
format(round(colSums(data_combat), digits = 0), big.mark = ",")

# add TMM normalization step and check distributions
combat_tmm <- calcNormFactors(data_combat)
data_combat_tmm <- sweep(data_combat, 2, combat_tmm, FUN = "/")

boxplot(log2(data_combat_tmm), notch = TRUE, col = rep(c("red", "green", "blue"), each = 6), 
        main = "ComBat batch corrected with TMM\nExp1 (red), Exp2 (green), Exp3 (blue)")

plotDensities(log2(data_combat_tmm), col = rep(c("red", "green", "blue"), 6), main = "ComBat/TMM data")

# check column totals and clustering
format(round(colSums(data_combat_tmm), digits = 0), big.mark = ",")
plotMDS(log2(data_combat_tmm), col = col_vec, main = "TMM-normalized Combat data")

# function computes CVs per time point
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

# CVs probe individual protein behaviour more than intensity distributions
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

# print out the average median CVs
(raw_med_cv <- round(mean(apply(list_raw[[3]], 2, median)), 2))
(sl_med_cv <- round(mean(apply(list_sl[[3]], 2, median)), 2))
(irs_med_cv <- round(mean(apply(list_irs[[3]], 2, median)),2))
(combat_med_cv <- round(mean(apply(list_combat[[3]], 2, median)), 2))

# do a couple of multi-panel to see effects of batch correction
par(mfrow = c(2,2))
plotMDS(log2(data_raw_tmm), col = col_vec, main = 'Raw data with TMM')
plotMDS(log2(data_sl_tmm), col = col_vec, main = 'TMM-normalized SL data\n(best non-IRS treatment)')
boxplot(list_raw[[3]], notch = TRUE, main = "RAW/TMM CVs (Ave = 62%)")
boxplot(list_sl[[3]], notch = TRUE, main = "SL/TMM CVs (Ave = 55%)")

# now after batch corrections
par(mfrow = c(2, 2))
plotMDS(log2(data_irs_tmm), col = col_vec, main = "TMM-normalized IRS data")
plotMDS(log2(data_combat_tmm), col = col_vec, main = "TMM-normalized Combat data")
boxplot(list_irs[[3]], notch = TRUE, main = "IRS/TMM CVs (Ave = 13%)")
boxplot(list_combat[[3]], notch = TRUE, main = "ComBat/TMM CVs (Ave = 18%)", ylim = c(0, 150))

library(psych)
# lets compare the combination of SL and TMM normalizations to SL/IRS/TMM 
# again using the idea that replicates of the same time point should be similar
pairs.panels(log2(data_sl_tmm[c(1, 7, 13)]), lm = TRUE, main = "SL/TMM E15")
pairs.panels(log2(data_irs_tmm[c(1, 7, 13)]), lm = TRUE, main = "SL/IRS/TMM E15")
pairs.panels(log2(data_sl_tmm[c(2, 8, 14)]), lm = TRUE, main = "SL/TMM E18")
pairs.panels(log2(data_irs_tmm[c(2, 8, 14)]), lm = TRUE, main = "SL/IRS/TMM E18")
pairs.panels(log2(data_sl_tmm[c(3, 9, 15)]), lm = TRUE, main = "SL/TMM P0")
pairs.panels(log2(data_irs_tmm[c(3, 9, 15)]), lm = TRUE, main = "SL/IRS/TMM P0")
pairs.panels(log2(data_sl_tmm[c(4, 10, 16)]), lm = TRUE, main = "SL/TMM P3")
pairs.panels(log2(data_irs_tmm[c(4, 10, 16)]), lm = TRUE, main = "SL/IRS/TMM P3")
pairs.panels(log2(data_sl_tmm[c(5, 11, 17)]), lm = TRUE, main = "SL/TMM P6")
pairs.panels(log2(data_irs_tmm[c(5, 11, 17)]), lm = TRUE, main = "SL/IRS/TMM P6")
pairs.panels(log2(data_sl_tmm[c(6, 12, 18)]), lm = TRUE, main = "SL/TMM P9")
pairs.panels(log2(data_irs_tmm[c(6, 12, 18)]), lm = TRUE, main = "SL/IRS/TMM P9")

# at this point, normalized data could be written to CSV files for saving