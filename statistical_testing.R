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

# set up the sample mapping
group <- rep(c("early", "early", "middle", "middle", "late", "late"), 3)

# make group into factors and set the order
group <- factor(group, levels = c("early", "middle", "late"))
str(group)

# create a DGEList object with our data
y_sl <- DGEList(counts = data_sl, group = group)

# y_sl is a list: y_sl$counts is the data, and y_sl$samples has interesting content
y_sl$samples

# we will skip running TMM (using the calcNormFactors function)
# we need to estimate the dispersion terms (global and local)
y_sl <- estimateDisp(y_sl)
plotBCV(y_sl, main = "Biological variation SL only (no IRS)")

# the exact test object has columns like fold-change, CPM, and p-values
et_sl <- exactTest(y_sl, pair = c("early", "late"))
summary(decideTestsDGE(et_sl)) # this counts up, down, and unchanged genes (here it is proteins)

# the topTags function adds the BH FDR values to an exactTest data frame. Make sure we do not change the row order!
tt_sl <- topTags(et_sl, n = Inf, sort.by = "none")
tt_sl <- tt_sl$table # tt_sl is a list. We just need the data frame table

# add the default value as a new column
tt_sl$candidate <- "no"
tt_sl[which(tt_sl$FDR <= 0.10 & tt_sl$FDR > 0.05), dim(tt_sl)[2]] <- "low"
tt_sl[which(tt_sl$FDR <= 0.05 & tt_sl$FDR > 0.01), dim(tt_sl)[2]] <- "med"
tt_sl[which(tt_sl$FDR <= 0.01), dim(tt_sl)[2]] <- "high"
tt_sl$candidate <- as.factor(tt_sl$candidate)

# what does tt_sl look like?
head(tt_sl)

# what does the test p-value distribution look like?
ggplot(tt_sl, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_sl$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("SLNorm only p-value distribution")

# for plotting results, we will use the average intensities for the 6 early and the 6 late samples
early <- c(1, 2, 7, 8, 13, 14)
late <- c(5, 6, 11, 12, 17, 18)
de_sl <- data.frame(rowMeans(data_sl[early]), rowMeans(data_sl[late]), tt_sl$candidate)
colnames(de_sl) <- c("early", "late", "candidate")
head(de_sl)
volcano_sl <- data.frame(log2(rowMeans(data_sl[early])/rowMeans(data_sl[late])), 
                         log10(tt_sl$FDR)*(-1), tt_sl$candidate)
colnames(volcano_sl) <- c("FoldChange", "FDR", "candidate")
head(volcano_sl)

# start with MA plot
library(scales)
temp <- data.frame(log2((de_sl$early + de_sl$late)/2), log2(de_sl$late/de_sl$early), de_sl$candidate)
colnames(temp) <- c("Ave", "FC", "candidate")
ggplot(temp, aes(x = Ave, y = FC)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_continuous("FC (late / early)") +
  scale_x_continuous("Ave_intensity") +
  ggtitle("Before IRS early vs late (MA plot)") + 
  geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
  geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate MA plots
ggplot(temp, aes(x = Ave, y = FC)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_continuous("FC (late / early)") +
  scale_x_continuous("Ave_intensity") +
  geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
  geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("Before IRS, separated by candidate (MA plots)")

# make the combined candidate corelation plot
ggplot(de_sl, aes(x = early, y = late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Before IRS early vs late") + 
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate corelation plots
ggplot(de_sl, aes(x = early, y = late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("Before IRS, separated by candidate")

# make a volcano plot
ggplot(volcano_sl, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ggtitle("Before IRS Volcano Plot")

# create a DGEList object with the IRS data
y_irs <- DGEList(counts = data_irs, group = group)

# y_irs is a list: y_irs$counts is the data, and y_irs$samples has interesting content
y_irs$samples

# we will skip running TMM (using the calcNormFactors function)
# we need to estimate the dispersion terms (global and local)
y_irs <- estimateDisp(y_irs)
plotBCV(y_irs, main = "Biological variation with IRS")

# the exact test object has columns like fold-change, CPM, and p-values
et_irs <- exactTest(y_irs, pair = c("early", "late"))
summary(decideTestsDGE(et_irs)) # this counts up, down, and unchanged genes

# the topTags function adds the BH FDR values to an exactTest data frame. Make sure we do not change the row order!
tt_irs <- topTags(et_irs, n = Inf, sort.by = "none")
tt_irs <- tt_irs$table # tt_sl is a list. We just need the data frame table

# add the default value as a new column
tt_irs$candidate <- "no"
tt_irs[which(tt_irs$FDR <= 0.10 & tt_irs$FDR > 0.05), dim(tt_irs)[2]] <- "low"
tt_irs[which(tt_irs$FDR <= 0.05 & tt_irs$FDR > 0.01), dim(tt_irs)[2]] <- "med"
tt_irs[which(tt_irs$FDR <= 0.01), dim(tt_irs)[2]] <- "high"
tt_irs$candidate <- as.factor(tt_irs$candidate)

# what does tt_sl look like?
head(tt_irs)

# what does the test p-value distribution look like?
ggplot(tt_irs, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_sl$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("IRS data p-value distribution")

# for plotting results, we will use the average intensities for the 6 early and the 6 late samples
de_irs <- data.frame(rowMeans(data_irs[early]), rowMeans(data_irs[late]), tt_irs$candidate)
colnames(de_irs) <- c("early", "late", "candidate")
head(de_irs)
volcano_irs <- data.frame(-1*tt_irs$logFC, -1*log10(tt_irs$FDR), tt_irs$candidate)
colnames(volcano_irs) <- c("FoldChange", "FDR", "candidate")
head(volcano_sl)

# start with MA plot
temp <- data.frame(log2((de_irs$early + de_irs$late)/2), log2(de_irs$late/de_irs$early), de_irs$candidate)
colnames(temp) <- c("Ave", "FC", "candidate")
ggplot(temp, aes(x = Ave, y = FC)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_continuous("FC (late / early)") +
  scale_x_continuous("Ave_intensity") +
  ggtitle("After IRS early vs late (MA plot)") + 
  geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
  geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate MA plots
ggplot(temp, aes(x = Ave, y = FC)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_continuous("FC (late / early)") +
  scale_x_continuous("Ave_intensity") +
  geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
  geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("After IRS, separated by candidate (MA plots)")

# make the combined candidate corelation plot
ggplot(de_irs, aes(x = early, y = late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("After IRS early vs late") + 
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate corelation plots
ggplot(de_irs, aes(x = early, y = late)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("IRS, separated by candidate")

# make a volcano plot
ggplot(volcano_irs, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ggtitle("After IRS Volcano Plot")

library(ggExtra)
# add marginal distrubution histograms to basic correlation plot (good starting point)
ggplot()
corr_plot <- ggplot(de_sl, aes(x = log10(early), y = log10(late))) +
  geom_point() + ggtitle("Before IRS")
ggMarginal(corr_plot, type = "histogram")

ggplot()
corr_plot <- ggplot(de_irs, aes(x = log10(early), y = log10(late))) +
  geom_point() + ggtitle("After IRS")
ggMarginal(corr_plot, type = "histogram")

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

# get the averages and SDs for the SL data and plot the crystallins
sl_list <- make_CVs(data_sl)
plot_proteins(sl_list, annotate_df)

# get the averages and SDs for the IRS data and plot the crystallins
irs_list <- make_CVs(data_irs)
plot_proteins(irs_list, annotate_df)

