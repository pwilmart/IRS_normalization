# Analysis of IOVS mouse lens data (Supplemental Table S01):
# Khan, Shahid Y., et al. "Proteome Profiling of Developing Murine Lens Through Mass Spectrometry." 
# Investigative Ophthalmology & Visual Science 59.1 (2018): 100-107.

# load libraries
library(tidyverse) # modern R packages for big data analysis
library(limma) # edgeR will load this if we do not
library(edgeR)

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

# figure out the global target scaling value
target <- mean(c(colSums(exp1_raw), colSums(exp2_raw), colSums(exp3_raw)))

# do the sample loading normalizations (scale to the target value)
norm_facs <- target / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs, FUN = "*")
norm_facs <- target / colSums(exp3_raw)
exp3_sl <- sweep(exp3_raw, 2, norm_facs, FUN = "*")
data_sl <- cbind(exp1_sl, exp2_sl, exp3_sl)

# make data frame with row means from each TMT experiment
refs <- tibble(rowMeans(exp1_sl), rowMeans(exp2_sl), rowMeans(exp3_sl))
colnames(refs) <- c("ave1", "ave2", "ave3")

# make new data frame with ratio values (time points divided by the "reference")
data_ratio <- exp1_sl / refs$ave1
data_ratio <- cbind(data_ratio, exp2_sl / refs$ave2)
data_ratio <- cbind(data_ratio, exp3_sl / refs$ave3)
row.names(data_ratio) <- annotate_df$`Protein Accession No.`

# let's see what box plots look like
boxplot(data_ratio, col = rep(c("red", "green", "blue"), each = 6), main = "Ratio to 'average' standard")
boxplot(log2(data_ratio), col = rep(c("red", "green", "blue"), each = 6), 
        main = "Log2 ratio to 'average' standard")
plotDensities(log2(data_ratio), col = rep(c("red", "green", "blue"), each = 6), main = "Ratios")

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

list_ratio <- make_CVs(data_ratio)

# compare CV distributions
boxplot(list_ratio[[3]], notch = TRUE, main = "Ratios CVs (14.3% ave)", ylim = c(0, 150))

round(mean(apply(list_ratio[[3]], 2, median)), 2)

plotMDS(log2(data_ratio), col = rep(c("red", "green", "blue"), each = 6),
        main = "log2 of ratios to 'standard'")
plotMDS(log2(data_ratio), col = c("red", "orange", "black", "green", "blue", "violet"),
        main = "log2 of ratios to 'standard'")

# collect the 6 early versus the 6 late columns
early <- data_ratio[c(1, 2, 7, 8, 13, 14)]
late <- data_ratio[c(5, 6, 11, 12, 17, 18)]

# get the data into a separate data frame so we can apply the t-test to each row
# to be safe, take the log2 of the ratios (maybe this helps for ratios less than one?)
pair <- data.frame(log2(early), log2(late))

# add average ratio columns (non-logged ratios), fold-change column, and row names
pair$ave_early <- rowMeans(early)
pair$ave_late  <- rowMeans(late)
pair$logFC <- log2(pair$ave_late / pair$ave_early)
row.names(pair) <- annotate_df$`Protein Accession No.`

# apply the t-test (from:https://stackoverflow.com/questions/25496693/
# how-to-do-t-test-two-samples-to-rows-of-a-dataframe-using-apply-family-functio)
# data is set up to do the two-sample t.test (a built-in) on each row
t.result <- apply(pair, 1, function(x) t.test(x[1:6], x[7:12]))
# extract the p-value column from the t-test thingy 
pair$p_value <- unlist(lapply(t.result, function(x) x$p.value))
# do a Benjamini-Hochberg multiple testing correction
pair$fdr <- p.adjust(pair$p_value, method = "BH")
head(pair)

# add a DE candidate status column
pair$candidate <- "no"
pair[which(pair$fdr <= 0.10 & pair$fdr > 0.05), dim(pair)[2]] <- "low"
pair[which(pair$fdr <= 0.05 & pair$fdr > 0.01), dim(pair)[2]] <- "med"
pair[which(pair$fdr <= 0.01), dim(pair)[2]] <- "high"
pair$candidate <- as.factor(pair$candidate)

# count up, down and the rest (FDR less than 0.05)
all <- dim(pair)[1]
up <- dim(pair[(pair$fdr <= 0.05) & (pair$logFC > 0.0), ])[1]
down <- dim(pair[(pair$fdr <= 0.05) & (pair$logFC <= 0.0), ])[1]
up 
all - up - down
down

# what does the test p-value distribution look like?
ggplot(pair, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(pair$p_value, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("log2 ratios t-test: p-value distribution")

# make the combined candidate corelation plot
de_ratio <- data.frame(pair$ave_early, pair$ave_late, pair$candidate)
colnames(de_ratio) <- c("early", "late", "candidate")
ggplot(de_ratio, aes(x = early, y = late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  ggtitle("t-test ratios: early vs late")

# make separate corelation plots
ggplot(de_ratio, aes(x = early, y = late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  facet_wrap(~ candidate) +
  ggtitle("t-test ratios: separated by candidate")

# make a volcano plot
vc <- data.frame(-1*pair$logFC, -1*log10(pair$fdr), pair$candidate)
colnames(vc) <- c("logFC", "fdr", "candidate")
ggplot(vc, aes(x = logFC, y = fdr)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ggtitle("t-test ratios: Volcano Plot")

# create a basic design matrix
group = as.factor(rep(c("early", "early", "middle", "middle", "late", "late"), 3))
group <- factor(group, levels(group)[c(1, 3, 2)]) # set the factor order
design <- model.matrix(~ 0 + group)
colnames(design) <- c("early", "middle", "late")
design

# the contrast is where the early versus late sub-comparison is selected
contrast <- makeContrasts(late-early, levels = design)
contrast

# do the linear model fitting
data_limma <- log2(data_ratio)
fit <- lmFit(data_limma, design)

# get the fit for the contrast of interest
fit2 <- contrasts.fit(fit, contrast)

# do the empirical Bayes moderation of the test statistic
fit2 <- eBayes(fit2, trend = TRUE)

# grab the information in topTable so we can get the data to plot candidates
# the coef parameter has to do with the contrast of interest
# specify no sorting of results and a number that is longer than the data table
tt_limma <- topTable(fit2, coef = 1, sort.by = "none", number = 10000)

# let's see what columns we have
head(tt_limma)
summary(decideTests(fit2))

# add average ratio columns (not logged), fold-change, and row names
pair <- data_limma
pair$ave_early <- rowMeans(early)
pair$ave_late  <- rowMeans(late)
pair$logFC <- log2(pair$ave_late / pair$ave_early)
row.names(pair) <- annotate_df$`Protein Accession No.`

# get the p-values and FDRs from the top tests
pair$p_value <- tt_limma$P.Value
pair$fdr <- tt_limma$adj.P.Val

# see if we have what we want
# head(pair)

# add a DE candidate status column
pair$candidate <- "no"
pair[which(pair$fdr <= 0.10 & pair$fdr > 0.05), dim(pair)[2]] <- "low"
pair[which(pair$fdr <= 0.05 & pair$fdr > 0.01), dim(pair)[2]] <- "med"
pair[which(pair$fdr <= 0.01), dim(pair)[2]] <- "high"
pair$candidate <- as.factor(pair$candidate)

# count up, down and the rest (FDR less than 0.05)
all <- dim(pair)[1]
up <- dim(pair[(pair$fdr <= 0.05) & (pair$logFC > 0.0), ])[1]
down <- dim(pair[(pair$fdr <= 0.05) & (pair$logFC <= 0.0), ])[1]
up 
all - up - down
down

# what does the test p-value distribution look like?
ggplot(pair, aes(p_value)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(pair$p_value, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("log2 ratios t-test: p-value distribution")

# make the combined candidate corelation plot
de_ratio <- data.frame(pair$ave_early, pair$ave_late, pair$candidate)
colnames(de_ratio) <- c("early", "late", "candidate")
ggplot(de_ratio, aes(x = early, y = late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  ggtitle("limma log ratios: early vs late")

# make separate corelation plots
ggplot(de_ratio, aes(x = early, y = late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  facet_wrap(~ candidate) +
  ggtitle("limma log ratios: separated by candidate")

# make a volcano plot
vc <- data.frame(-1*pair$logFC, -1*log10(pair$fdr), pair$candidate)
colnames(vc) <- c("logFC", "fdr", "candidate")
ggplot(vc, aes(x = logFC, y = fdr)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ggtitle("limma log ratios: Volcano Plot")

# collect up some results to write out
final <- cbind(annotate_df, data_raw, data_sl, refs, data_ratio, pair)
write.csv(final, file = "final_part3.csv")

library(psych)
# lets compare the combination of SL and TMM normalizations to SL/IRS/TMM 
# again using the idea that replicates of the same time point should similar
pairs.panels(log2(data_ratio[c(1, 7, 13)]), lm = TRUE, main = "SL/Ratio E15")
pairs.panels(log2(data_ratio[c(2, 8, 14)]), lm = TRUE, main = "SL/Ratio E18")
pairs.panels(log2(data_ratio[c(3, 9, 15)]), lm = TRUE, main = "SL/Ratio P0")
pairs.panels(log2(data_ratio[c(4, 10, 16)]), lm = TRUE, main = "SL/Ratio P3")
pairs.panels(log2(data_ratio[c(5, 11, 17)]), lm = TRUE, main = "SL/Ratio P6")
pairs.panels(log2(data_ratio[c(6, 12, 18)]), lm = TRUE, main = "SL/Ratio P9")

# create some functions
subset_proteins <- function(df, proteins, labels, value_col) {  
  df_subset <- subset(df, df$protein %in% proteins)
  df_subset <- df_subset[match(proteins, df_subset$protein), ]
  df_subset$proteins <- as.factor(labels)
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
    xlab("Developmental Time") + ylab("Protein Intensity Ratio") +
    ggtitle("Alpha crystallin expression")
  print(alpha_plot)
  # plot betas
  ave_df_wide <- subset_proteins(ave_df, betas, beta_labels, "intensity")
  sd_df_wide <- subset_proteins(sd_df, betas, beta_labels, "sd")
  ave_df_wide$sd <- sd_df_wide$sd
  beta_plot <- ggplot(ave_df_wide, aes(time, intensity, fill=proteins)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin = intensity-sd, ymax = intensity+sd), position = position_dodge(0.9), width = 0.4) +
    xlab("Developmental Time") + ylab("Protein Intensity Ratio") +
    ggtitle("Beta crystallin expression") 
  print(beta_plot)
  # plot gammas
  ave_df_wide <- subset_proteins(ave_df, gammas, gamma_labels, "intensity")
  sd_df_wide <- subset_proteins(sd_df, gammas, gamma_labels, "sd")
  ave_df_wide$sd <- sd_df_wide$sd
  gamma_plot <- ggplot(ave_df_wide, aes(time, intensity, fill=proteins)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin = intensity-sd, ymax = intensity+sd), position = position_dodge(0.9), width = 0.4) +
    xlab("Developmental Time") + ylab("Protein Intensity Ratio") +
    ggtitle("Gamma crystallin expression")
  print(gamma_plot)
  return(NULL)
}

plot_proteins(list_ratio, annotate_df)
