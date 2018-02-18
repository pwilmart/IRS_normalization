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
group <- rep(c("E15", "E18", "P0", "P3", "P6", "P9"), 3)

# make group into factors and set the order
group <- factor(group, levels = c("E15", "E18", "P0", "P3", "P6", "P9"))

# create a DGEList object with our data
y_sl <- DGEList(counts = data_sl, group = group)

# we need to run TMM norm and estimate the dispersion
y_sl <- calcNormFactors(y_sl)
y_sl <- estimateDisp(y_sl)
y_sl$samples
plotBCV(y_sl, main = "Biological variation SL only (no IRS)")

# we need to transform the SL data with TMM factors for later plotting
sl_tmm <- calcNormFactors(data_sl)
data_sl_tmm <- sweep(data_sl, 2, sl_tmm, FUN = "/") # this is data after SL and TMM on original scale

# the exact test object has columns like fold-change, CPM, and p-values
et_sl <- exactTest(y_sl, pair = c("P0", "P3"))
summary(decideTestsDGE(et_sl)) # this counts up, down, and unchanged genes (here it is proteins)

# the topTags function adds the BH FDR values to an exactTest data frame. Make sure we do not change the row order!
tt_sl <- topTags(et_sl, n = Inf, sort.by = "none")
tt_sl <- tt_sl$table # tt_sl is a list. We just need the data frame table

# add the default value as a new column
tt_sl$candidate <- "no"
tt_sl[which(tt_sl$FDR <= 0.10 & tt_sl$FDR > 0.05), dim(tt_sl)[2]] <- "low"
tt_sl[which(tt_sl$FDR <= 0.05 & tt_sl$FDR > 0.01), dim(tt_sl)[2]] <- "med"
tt_sl[which(tt_sl$FDR <= 0.01), dim(tt_sl)[2]] <- "high"
tt_sl$candidate <- factor(tt_sl$candidate, levels = c("high", "med", "low", "no"))

# what does tt_sl look like?
head(tt_sl)

# what does the test p-value distribution look like?
ggplot(tt_sl, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_sl$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("P0 vs P3 (SLNorm only) p-value distribution")

# for plotting results, we will use the average intensities for the 3 P0 and the 3 P3 samples
P0 <- c(3, 9, 15)
P3 <- c(4, 10, 16)
de_sl <- data.frame(rowMeans(data_sl_tmm[P0]), rowMeans(data_sl_tmm[P3]), tt_sl$candidate)
colnames(de_sl) <- c("P0", "P3", "candidate")
volcano_sl <- data.frame(log2(rowMeans(data_sl[P0])/rowMeans(data_sl[P3])), log10(tt_sl$FDR)*(-1), tt_sl$candidate)
colnames(volcano_sl) <- c("FoldChange", "FDR", "candidate")

# start with MA plot
library(scales)
temp <- data.frame(log2((de_sl$P0 + de_sl$P3)/2), log2(de_sl$P3/de_sl$P0), de_sl$candidate)
colnames(temp) <- c("Ave", "FC", "candidate")
ggplot(temp, aes(x = Ave, y = FC)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_continuous("FC (P3 / P0)") +
  scale_x_continuous("Ave_intensity") +
  ggtitle("Without IRS P0 vs P3 (MA plot)") + 
  geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
  geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate MA plots
ggplot(temp, aes(x = Ave, y = FC)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_continuous("FC (P3 / P0)") +
  scale_x_continuous("Ave_intensity") +
  geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
  geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("Without IRS, separated by candidate (MA plots)")

# make the combined candidate corelation plot
ggplot(de_sl, aes(x = P0, y = P3)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("Without IRS P0 vs P3") + 
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate corelation plots
ggplot(de_sl, aes(x = P0, y = P3)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("Without IRS, separated by candidate")

# make a volcano plot
ggplot(volcano_sl, aes(x = FoldChange, y = FDR)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  xlab("Fold-Change (Log2)") +
  ylab("-Log10 FDR") +
  ggtitle("Without IRS Volcano Plot")

# create a DGEList object with the IRS data
y_irs <- DGEList(counts = data_irs, group = group)

# we need to normalize and estimate the dispersion terms (global and local)
y_irs <- calcNormFactors(y_irs)
y_irs <- estimateDisp(y_irs)
y_irs$samples
plotBCV(y_irs, main = "Biological variation with IRS")

# we need to transform the IRS data with TMM factors for later plotting
irs_tmm <- calcNormFactors(data_sl)
data_irs_tmm <- sweep(data_irs, 2, irs_tmm, FUN = "/") # this is data after SL and TMM on original scale

# the exact test object has columns like fold-change, CPM, and p-values
et_irs <- exactTest(y_irs, pair = c("P0", "P3"))
summary(decideTestsDGE(et_irs)) # this counts up, down, and unchanged genes

# the topTags function adds the BH FDR values to an exactTest data frame. Make sure we do not change the row order!
tt_irs <- topTags(et_irs, n = Inf, sort.by = "none")
tt_irs <- tt_irs$table # tt_sl is a list. We just need the data frame table

# add the default value as a new column
tt_irs$candidate <- "no"
tt_irs[which(tt_irs$FDR <= 0.10 & tt_irs$FDR > 0.05), dim(tt_irs)[2]] <- "low"
tt_irs[which(tt_irs$FDR <= 0.05 & tt_irs$FDR > 0.01), dim(tt_irs)[2]] <- "med"
tt_irs[which(tt_irs$FDR <= 0.01), dim(tt_irs)[2]] <- "high"
tt_irs$candidate <- factor(tt_irs$candidate, levels = c("high", "med", "low", "no"))

# what does tt_sl look like?
head(tt_irs)

# what does the test p-value distribution look like?
ggplot(tt_irs, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(tt_irs$PValue, breaks = 100, plot = FALSE)$counts[26:100])) +
  ggtitle("P0 vs P3 (after IRS) p-value distribution")

# for plotting results, we will use the average intensities for the P0 and the P3 samples
de_irs <- data.frame(rowMeans(data_irs_tmm[P0]), rowMeans(data_irs_tmm[P3]), tt_irs$candidate)
colnames(de_irs) <- c("P0", "P3", "candidate")
volcano_irs <- data.frame(-1*tt_irs$logFC, -1*log10(tt_irs$FDR), tt_irs$candidate)
colnames(volcano_irs) <- c("FoldChange", "FDR", "candidate")

# start with MA plot
temp <- data.frame(log2((de_irs$P0 + de_irs$P3)/2), log2(de_irs$P3/de_irs$P0), de_irs$candidate)
colnames(temp) <- c("Ave", "FC", "candidate")
ggplot(temp, aes(x = Ave, y = FC)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_continuous("FC (P3 / P0)") +
  scale_x_continuous("Ave_intensity") +
  ggtitle("After IRS P0 vs P3 (MA plot)") + 
  geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
  geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate MA plots
ggplot(temp, aes(x = Ave, y = FC)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_continuous("FC (P3 / P0)") +
  scale_x_continuous("Ave_intensity") +
  geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
  geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") + # 2-fold down
  facet_wrap(~ candidate) +
  ggtitle("After IRS, separated by candidate (MA plots)")

# make the combined candidate corelation plot
ggplot(de_irs, aes(x = P0, y = P3)) +
  geom_point(aes(color = candidate, shape = candidate)) +
  scale_y_log10() +
  scale_x_log10() +
  ggtitle("After IRS P0 vs P3") + 
  geom_abline(intercept = 0.0, slope = 1.0, color = "black") + # one-to-one line
  geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted") + # 2-fold up
  geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted") # 2-fold down

# make separate corelation plots
ggplot(de_irs, aes(x = P0, y = P3)) +
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
corr_plot <- ggplot(de_sl, aes(x = log10(P0), y = log10(P3))) +
  geom_point() + ggtitle("Before IRS")
ggMarginal(corr_plot, type = "histogram")

ggplot()
corr_plot <- ggplot(de_irs, aes(x = log10(P0), y = log10(P3))) +
  geom_point() + ggtitle("After IRS")
ggMarginal(corr_plot, type = "histogram")
