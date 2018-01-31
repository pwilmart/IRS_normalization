# Analysis of IOVS mouse lens data (Supplemental Table S01)

# load libraries
library(tidyverse)
library(edgeR)
library(ggExtra)
library(car)

# read the Supplemental 01 file (saved as a CSV file)
data_start <- read_csv("iovs-58-13-55_s01.csv")

# filter out proteins not seen in all three runs
data_no_na <- na.omit(data_start)

# fix the column headers
col_headers <- colnames(data_no_na)
col_headers <- str_replace(col_headers, " {2,3}", " ")
col_headers <- str_replace(col_headers, "Reporter ion intensities ", "")
colnames(data_no_na) <- col_headers

# save the annotation columns for later and remove from data frame
annotate_df <- data_no_na[1:2]
data_no_na <- data_no_na[3:20]

# separate the TMT data by experiment
exp1 <- data_no_na[c(1:6)]
exp2 <- data_no_na[c(7:12)]
exp3 <- data_no_na[c(13:18)]

# figure out the global scaling value
target <- mean(c(colSums(exp1), colSums(exp2), colSums(exp3)))

# do the sample loading normalization before the IRS normalization
# there is a different correction factor for each column
norm_facs <- target / colSums(exp1)
exp1 <- sweep(exp1, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2)
exp2 <- sweep(exp2, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3)
exp3 <- sweep(exp3, 2, norm_facs, FUN = "*")

# make a pre-IRS data frame after sample loading norms
SLNorm <- cbind(exp1, exp2, exp3)

# make working frame with row sums from each frame
irs <- tibble(rowSums(exp1), rowSums(exp2), rowSums(exp3))
colnames(irs) <- c("sum1", "sum2", "sum3")
# get the average intensity for each protein
irs$average <- rowMeans(irs)
# compute the scaling factor vectors
irs$fac1 <- irs$average / irs$sum1
irs$fac2 <- irs$average / irs$sum2
irs$fac3 <- irs$average / irs$sum3

# make new data frame with normalized data
IRSNorm <- exp1 * irs$fac1
IRSNorm <- cbind(IRSNorm, exp2 * irs$fac2)
IRSNorm <- cbind(IRSNorm, exp3 * irs$fac3)

# look at the clustering
g <- c(rep.int(c(1, 2, 3, 4, 5, 6), 3))
col_vec2 <- c("red", "orange", "black", "green", "blue", "violet")
temp = DGEList(counts = SLNorm, group = g)
plotMDS(temp$counts, col = col_vec2, main = "Sample Loading Normed intensities - no IRS\nUnderlying data for DGEList object")
plotMDS(temp, col = col_vec2, main = "Sample Loading (SL) Norm Only")
temp = DGEList(counts = IRSNorm, group = g)
plotMDS(temp, col = col_vec2, main = 'Internal Reference Scaling (IRS) Norm')
temp <- NULL

# look at some data summaries
col_vec1 <- c(rep_len("red", 18), rep_len("green", 18))
boxplot(cbind(log10(SLNorm), log10(IRSNorm)), main = "Before (red) and after (green) IRS", col = col_vec1)

# check CV distributions before and after
E15b <- SLNorm[c(1, 7, 13)]
E18b <- SLNorm[c(2, 8, 14)]
P0b <- SLNorm[c(3, 9, 15)]
P3b <- SLNorm[c(4, 10, 16)]
P6b <- SLNorm[c(5, 11, 17)]
P9b <- SLNorm[c(6, 12, 18)]

E15b$ave <- rowMeans(E15b)
E15b$sd <- apply(E15b[1:3], 1, sd)
E15b$cv <- 100 * E15b$sd / E15b$ave

E18b$ave <- rowMeans(E18b)
E18b$sd <- apply(E18b[1:3], 1, sd)
E18b$cv <- 100 * E18b$sd / E18b$ave

P0b$ave <- rowMeans(P0b)
P0b$sd <- apply(P0b[1:3], 1, sd)
P0b$cv <- 100 * P0b$sd / P0b$ave

P3b$ave <- rowMeans(P3b)
P3b$sd <- apply(P3b[1:3], 1, sd)
P3b$cv <- 100 * P3b$sd / P3b$ave

P6b$ave <- rowMeans(P6b)
P6b$sd <- apply(P6b[1:3], 1, sd)
P6b$cv <- 100 * P6b$sd / P6b$ave

P9b$ave <- rowMeans(P9b)
P9b$sd <- apply(P9b[1:3], 1, sd)
P9b$cv <- 100 * P9b$sd / P9b$ave

cv_df <- data.frame(E15b$cv, E18b$cv, P0b$cv, P3b$cv, P6b$cv, P9b$cv)

E15a <- IRSNorm[c(1, 7, 13)]
E18a <- IRSNorm[c(2, 8, 14)]
P0a <- IRSNorm[c(3, 9, 15)]
P3a <- IRSNorm[c(4, 10, 16)]
P6a <- IRSNorm[c(5, 11, 17)]
P9a <- IRSNorm[c(6, 12, 18)]

E15a$ave <- rowMeans(E15a)
E15a$sd <- apply(E15a[1:3], 1, sd)
E15a$cv <- 100 * E15a$sd / E15a$ave

E18a$ave <- rowMeans(E18a)
E18a$sd <- apply(E18a[1:3], 1, sd)
E18a$cv <- 100 * E18a$sd / E18a$ave

P0a$ave <- rowMeans(P0a)
P0a$sd <- apply(P0a[1:3], 1, sd)
P0a$cv <- 100 * P0a$sd / P0a$ave

P3a$ave <- rowMeans(P3a)
P3a$sd <- apply(P3a[1:3], 1, sd)
P3a$cv <- 100 * P3a$sd / P3a$ave

P6a$ave <- rowMeans(P6a)
P6a$sd <- apply(P6a[1:3], 1, sd)
P6a$cv <- 100 * P6a$sd / P6a$ave

P9a$ave <- rowMeans(P9a)
P9a$sd <- apply(P9a[1:3], 1, sd)
P9a$cv <- 100 * P9a$sd / P9a$ave

cv_df <- cbind(cv_df, data.frame(E15a$cv, E18a$cv, P0a$cv, P3a$cv, P6a$cv, P9a $cv))
col_vec3 <- c(rep_len("red", 6), rep_len("green", 6))
boxplot(cv_df, notch = TRUE, col = col_vec3, main = "CV distributions before (red) and after (green)")

# compute average median CVs for before and after
(before <- mean(apply(cv_df[1:6], 2, median)))
(after <- mean(apply(cv_df[7:12], 2, median)))

# compare E15+E18 to P6+P9 with edgeR 
# before IRS
yb <- DGEList(counts = SLNorm, group = rep(c(1, 1, 2, 2, 3, 3), 3))
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
ya <- DGEList(counts = IRSNorm, group = rep(c(1, 1, 2, 2, 3, 3), 3))
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

# add marginal distrubution histograms to basic correlation plot (good starting point)
ggplot()
corr_plot <- ggplot(DE_before, aes(x = log10(Embryo), y = log10(Late))) +
  geom_point() + ggtitle("Before IRS")
ggMarginal(corr_plot, type = "histogram")

ggplot()
corr_plot <- ggplot(DE_after, aes(x = log10(Embryo), y = log10(Late))) +
  geom_point() + ggtitle("After IRS")
ggMarginal(corr_plot, type = "histogram")

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

# compare replicates at same time point between TMT experiments
# E15, Before IRS
scatterplotMatrix(log10(E15b[1:3]))
# E15, After IRS
scatterplotMatrix(log10(E15a[1:3]))
# E18, Before IRS
scatterplotMatrix(log10(E18b[1:3]))
# E18, After IRS
scatterplotMatrix(log10(E18a[1:3]))
# P0, Before IRS
scatterplotMatrix(log10(P0b[1:3]))
# P0, After IRS
scatterplotMatrix(log10(P0a[1:3]))
# P3, Before IRS
scatterplotMatrix(log10(P3b[1:3]))
# P3, After IRS
scatterplotMatrix(log10(P3a[1:3]))
# P6, Before IRS
scatterplotMatrix(log10(P6b[1:3]))
# P6, After IRS
scatterplotMatrix(log10(P6a[1:3]))
# P9, Before IRS
scatterplotMatrix(log10(P9b[1:3]))
# P9, After IRS
scatterplotMatrix(log10(P9a[1:3]))