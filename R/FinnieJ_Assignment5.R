#### Assignment 5: Analysis of Normalization Methods in Gene Expression Studies of Rat Spermatocytes ----
#
# Author: Jessica Finnie
#
# GitHub: https://github.com/finniej/6210_Assignment5
#
# Date: December 6, 2024
#
# Course: BINF*6210 Software Tools
#
# Instructor: Dr. Karl Cottenie
#
# Purpose: To investigate the impacts of normalizing RNAseq data prior to differential expression analysis.
#
# Attributions: Data sourced from Tian et al. (2020).
#
####

# Load Required Libraries & Setting Working Directory ----

# Core Libraries
library(tidyverse)
library(viridis)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(reshape2)

# RNA-seq Analysis
library(edgeR)
library(limma)
library(org.Rn.eg.db)
library(biomaRt)

# Utilities & Plotting
library(R.utils)
library(gplots)
library(DESeq2)
library(RColorBrewer)
library(gridExtra)

# Working Directory: Please set your working directory to the filepath of the submission folder.
# setwd(dir = "yourfilepath/FinnieJ_Assignment5")

# Data Acquisition ----
# Data for this analysis was sourced from Tian et al.'s (2020) study of the impacts of stress on spermatocyte gene expression in Rattus norvegicus (Rats).The data was retrieved from the Gene Expression Omnibus on November 26, 2024 @ 7:10PM. Do not uncomment this section unless necessary.

# Set working directory to data folder.
# setwd(dir = "./data/")

# Acquire the data, and unpack the .tar file.
# url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144485&format=file"
# utils::download.file(url, destfile="GSE144485_RAW.tar", mode="wb")
# utils::untar("GSE144485_RAW.tar", exdir = ".")
# files <- c("GSM4289045_c1.txt",      # Control Treatments 1-3
#            "GSM4289046_c2.txt",
#            "GSM4289047_c3.txt",
#            "GSM4289048_str31.txt",   # 3 Day Stress Treatment 1-3
#            "GSM4289049_str32.txt",
#            "GSM4289050_str33.txt",
#            "GSM4289051_str141.txt",  # 14 Day Stress Treatment 1-3
#            "GSM4289052_str142.txt",
#            "GSM4289053_str143.txt",
#            "GSM4289054_str211.txt",  # 21 Day Stress Treatment 1-3
#            "GSM4289055_str212.txt",
#            "GSM4289056_str213.txt")
# output_dir <- "./decompressed_files/"
# if (!dir.exists(output_dir)) dir.create(output_dir)
#
# for(i in paste(files, ".gz", sep="")) {
#   output_path <- paste0(output_dir, basename(i))  # Define the output path
#   R.utils::gunzip(i, destname = output_path, overwrite = TRUE)
# }
#
# Return to main working directory.
# setwd(dir = "..")

# Code Part 1 - Exploratory Analysis ----
# _ Loading Data -----
# Pull data from the data folder in the FinnieJ_Assignment5 directory.
# setwd(dir = "AddRelativeFilePathIfNeeded!")
decompressed_files <- list.files("./data/decompressed_files/", full.names = TRUE)

# Take a quick check of the files to confirm that they have been loaded correctly. The files contain Target IDs, lengths, effective lengths, estimated counts and transcripts per million for each transcript ID. Each file corresponds to one of the four treatments - Control, and three stress treatments (3, 14 and 21 days respectively.)
read.delim(decompressed_files[1], nrow = 5)

# _ Initial Processing of Transcript Data ----
# Creating an Estimated Counts object for exploring differential expression. Transcipts per million (tpm) will not be used, as the data will be manipulated prior to calculating tpm in this analysis.
dge_est_counts <- readDGE(decompressed_files, columns = c(1, 4))

# Checking the structure of the data by examining the first few rows of each object, and ensuring that are are no NAs in the the data.
head(dge_est_counts$counts) # for estimated counts
sum(is.na(dge_est_counts$counts)) # checking for NAs
head(dge_est_counts$samples) # Check sample names and metadata

# _ Manipulating Data (Assigning treatment groups, editing row and column names) ----
# For clarity and ease of analysis, each treatment will recieve a group delimiter.
group <- factor(c(
  "Control", "Control", "Control", "3Days", "3Days", "3Days",
  "14Days", "14Days", "14Days", "21Days", "21Days", "21Days"
))
dge_est_counts$samples$group <- group # Assign groups
dge_est_counts$samples$group <- factor(dge_est_counts$samples$group) # Factor groups

# Removing the Gene Expression Omnibus sample IDs, for simplicity in downstream analysis.
colnames(dge_est_counts) # Check the column names to find out how to change them
samplenames <- substring(colnames(dge_est_counts$counts), first = 38, last = nchar(colnames(dge_est_counts$counts))) # Remove the first 38 characters, assign to vector
colnames(dge_est_counts$counts) <- samplenames # Assign the cleaned names to the counts
rownames(dge_est_counts$samples) <- samplenames # Assign the cleaned names to the samples

# Check the updated names
head(colnames(dge_est_counts$counts))
head(rownames(dge_est_counts$samples))

# Edit the Transcript IDs to remove the decimal and trailing digits. Tong et al. recommend grouping transcript data, and biomaRt will be used as a reference database and it doesn't use decimal places in identifiers.
cleaned_target_ids <- gsub("\\..*", "", rownames(dge_est_counts$counts))
rownames(dge_est_counts$counts) <- cleaned_target_ids # Update the transcript IDs
head(rownames(dge_est_counts$counts)) # Check the updated IDs

# Keeping only counts with at least three instances of expression, to remove any unexpressed outliers.
keep <- rowSums(dge_est_counts$counts > 1) >= 3
dge_est_counts <- dge_est_counts[keep, ]

# Lane or batch effects will not be considered in this study, as the lanes/batches were not made available by Tong et al. (2021). In fact, no GEO data sets examined in an extensive search contained lane or batch information in the metadata.

# _ Simple plots of Raw Counts ----
# As exploratory analysis, a variety of plots will be constructed to examine the spread of the data. This is necessary to check for further outliers and to perform any needed quality control.

# Plot 1: Bar Plot of Total Reads with Variance
# Calculate total reads per sample
total_reads <- colSums(dge_est_counts$counts)

# Combine total reads with group information
reads_df <- data.frame(
  Sample = colnames(dge_est_counts$counts),
  TotalReads = total_reads,
  Group = dge_est_counts$samples$group
)
# Reorder the groups
reads_df <- reads_df %>%
  mutate(Group = factor(Group, levels = c("Control", "3Days", "14Days", "21Days")))

# Summarize the variance for each group
variance_summary <- reads_df %>%
  group_by(Group) %>%
  summarize(
    MeanReads = mean(TotalReads),
    Variance = var(TotalReads)  )

# Adding the mean reads and variance to the data frame by Group (Treatment).
reads_df <- reads_df %>%
  left_join(variance_summary, by = "Group")

# Calculate error bar limits for each sample
reads_df <- reads_df %>%
  mutate(
    Lower = TotalReads - sqrt(Variance),
    Upper = TotalReads + sqrt(Variance)  )

# Create the bar plot
ggplot(reads_df, aes(x = Sample, y = TotalReads, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, position = position_dodge(0.7)) +
  labs(
    title = "Total Reads per Sample with Variance",
    x = "Sample",
    y = "Total Reads"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Shows that there is little variance in total reads in control groups, with slightly more variance in the treatment groups. However, the variance is not extreme, suggesting that there are no outliers.

# Plot 2: PCA of Counts to Establish Relationships between Treatments & Variance
# Perform PCA
pca <- prcomp(t(dge_est_counts$counts))

# Storing PCA results and metadata to a dataframe for plotting
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Sample = colnames(dge_est_counts$counts),
  Group = dge_est_counts$samples$group)

# Calculate variance explained to make the PCA plot more informative
variance_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

# Create the PCA plot, with sample labels and colours. Variance explained is also included for each Principal Component.
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 3, box.padding = 0.5, max.overlaps = Inf) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "PCA of Gene Expression Data",
    x = paste0("PC1 (", variance_explained[1], "% Variance)"),
    y = paste0("PC2 (", variance_explained[2], "% Variance)")
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)  )

# Shows decent grouping among treatments. "str212.txt" (21 Day Stress) shows more variance than expected. Once relevant genes are examined, I will keep this in mind to look for any effects in that subset of data.

# _ Retrieving Gene Annotations ----
# Initial analysis is complete. No extreme outliers detected in the transcript reads. Now, I must convert the Transcript IDs to Gene IDs, in order to obtain the Gene Symbol, from Ensembl.

# Connect to Ensembl for Rats (Rattus norvegicus). Add 'host' as needed, as a mirror may not be needed to connect. Only uncomment if needed, Ensembl has been unreliable.
# mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")#, host = "https://useast.ensembl.org")
#
# Save the mart object to a file (e.g., "mart_ensembl.rda") Do not uncomment unless needed.
# save(mart, file = "./data/mart_ensembl.rda")
#
# Load the mart object from file, as Ensembl is slow/unreliable.
load("./data/mart_ensembl.rda")

# _ Mapping Transcript -> Gene ID -> Gene Symbol
# Exploring the mart object for ensembl information
dim(listAttributes(mart)) # Huge!

# Standardize transcript IDs to upper case prior to filtering.
rownames(dge_est_counts$counts) <- toupper(rownames(dge_est_counts$counts))

# Based on Tian et al.'s research, I will only be examining the following genes. The authors found the expression of these genes to be significantly impacted by stress. I want to see if normalization methods impact the significance.
target_genes <- c("Stra8", "Sycp3", "Piwil1", "Tnp1")

# # Retrieve transcript information for these genes. Only uncomment as needed.
# gene_four_annotations <- getBM(
#   attributes = c("ensembl_transcript_id", "ensembl_gene_id",
#                  "external_gene_name", "description"),
#   filters = "external_gene_name",  # Filter by gene names
#   values = target_genes,           # Target gene names
#   mart = mart
# )
#
# # Print the results
# print(gene_four_annotations)
# save(gene_four_annotations, file = "./data/four_gene_annotations.csv")

# Load the Gene annotation file.
load("./data/four_gene_annotations.csv")

# Extract transcript IDs of interest
transcript_ids <- gene_four_annotations$ensembl_transcript_id

# Filter the count matrix for the desired transcript IDs
filtered_counts <- dge_est_counts$counts[rownames(dge_est_counts$counts) %in% transcript_ids, ]

# Update the DGEList object with the filtered counts
dge_est_counts_filtered <- dge_est_counts
dge_est_counts_filtered$counts <- filtered_counts
dge_est_counts_filtered$samples <- dge_est_counts$samples

# Check the dimensions of the filtered counts
dim(dge_est_counts_filtered$counts)

# Check the rownames to ensure only the desired transcripts are present
rownames(dge_est_counts_filtered$counts)

# Create a named vector for mapping transcript IDs to gene IDs
transcript_to_gene <- setNames(
  gene_four_annotations$external_gene_name,
  gene_four_annotations$ensembl_transcript_id)

# Update row names of the counts matrix to gene IDs
rownames(dge_est_counts_filtered$counts) <- transcript_to_gene[rownames(dge_est_counts_filtered$counts)]

# Check the updated row names
dge_est_counts_filtered$counts

# Recalculate library sizes, to ensure downstream analysis is accurate.
dge_est_counts_filtered$samples$lib.size <- colSums(dge_est_counts_filtered$counts)

# A quick PCA to examine the filtered data!
pca <- prcomp(t(dge_est_counts_filtered$counts))
plot(pca$x[, 1:2], col = factor(dge_est_counts_filtered$samples$group), main = "PCA of Gene Expression Data")

# Looks good! Now that the data has been filtered for the four genes of interest, normalization methods can be explored.

# Code Part 2: Main Analysis ----
# In this section, normalization methods will be applied. Their initial effects will be examined, and then downstream limma analysis will be performed. Then the impacts of normalization will be established for the gene expression.

# Applying Normalization Methods!
# Control (No Normalization)
dge_none <- calcNormFactors(dge_est_counts_filtered, method = "none")
cpm_none <- cpm(dge_none, log = TRUE)

# TMM Normalized & CPM: accounts for library size and composition
dge_tmm <- calcNormFactors(dge_est_counts_filtered, method = "TMM")
cpm_tmm <- cpm(dge_tmm, log = TRUE)

# RLE Normalized CPM: tries to make the median of log ratios zero.
dge_rle <- calcNormFactors(dge_est_counts_filtered, method = "RLE")
cpm_rle <- cpm(dge_rle, log = TRUE)

# Upper Quartile Normalized & CPM: scales based on upper quantile gene counts
dge_uq <- calcNormFactors(dge_est_counts_filtered, method = "upperquartile")
cpm_uq <- cpm(dge_uq, log = TRUE)

# Combine all the CPM values into a single matrix for boxplot visualization
normalized_data <- data.frame(
  NONE = as.vector(cpm_none),
  TMM = as.vector(cpm_tmm),
  RLE = as.vector(cpm_rle),
  UQ = as.vector(cpm_uq))

# Create a label for the samples and genes
samples <- rep(colnames(dge_est_counts_filtered$counts), each = nrow(dge_est_counts_filtered$counts))
genes <- rep(rownames(dge_est_counts_filtered$counts), times = ncol(dge_est_counts_filtered$counts))

# Combine the labels into the data
normalized_data$Samples <- samples
normalized_data$Genes <- genes

# Reshape data for ggplot
normalized_long <- melt(normalized_data,
  id.vars = c("Samples", "Genes"),
  variable.name = "Normalization", value.name = "Log2_CPM")

# Plot 3: Density plot for normalized data to examine changes in spread.
ggplot(normalized_long, aes(x = Log2_CPM, color = Normalization)) +
  geom_density(linewidth = 0.8) +
  labs(
    title = "Density Plot of Log2 CPM After Normalization",
    subtitle = "Variance between Normalization Methods",
    x = "Log2 CPM",
    y = "Density"
  ) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))

# RLE shows the greatest variance in the spread of counts per million after normalization.

# Calculating Coefficient of Variation: Lower value = better normalization.
cv_none <- sd(cpm_none) / mean(cpm_none)
cv_tmm <- sd(cpm_tmm) / mean(cpm_tmm)
cv_rle <- sd(cpm_rle) / mean(cpm_rle)
cv_uq <- sd(cpm_uq) / mean(cpm_uq)

data.frame(
  Method = c("NONE", "TMM", "RLE", "UQ"),
  CV = c(cv_none, cv_tmm, cv_rle, cv_uq))
# Lowest CV is TMM, which was expected.

# Combine the data into a single data frame for PCA analysis.
normalized_data_long <- data.frame(
  gene = rep(rownames(dge_est_counts_filtered$counts), times = ncol(dge_est_counts_filtered$counts)),
  sample = rep(colnames(dge_est_counts_filtered$counts), each = nrow(dge_est_counts_filtered$counts)),
  TMM = as.vector(cpm_tmm),
  RLE = as.vector(cpm_rle),
  UQ = as.vector(cpm_uq),
  NONE = as.vector(cpm_none))

# Reshape the data to long format for the PCA.
normalized_data_long <- normalized_data_long %>%
  gather(key = "method", value = "value", NONE, TMM, RLE, UQ)

# Check the structure of the reshaped data
head(normalized_data_long)    # Looks good!

# Plot 4: Comparing the means of the normalization methods for each gene.
ggplot(normalized_data_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = list(c("TMM", "RLE"), c("TMM", "UQ"), c("RLE", "UQ"), c("NONE", "TMM", c("NONE", "UQ", c("NONE", "RLE")))),
    method = "wilcox.test") +
  labs(
    title = "Comparison of Normalization Methods",
    y = "Log2 CPM",
    x = "Normalization Method") +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")

# Boxes represent the log2 of the variance in counts per million. The veritcal bars and scattered dots represent the spread of the variance. The horizontal bars represnt the wilcox test comparison between the means of the methods. The normalization methods were successfully compared to each other using this wilcox test. The control test with no normalization could not be compared to all methods due to limitations of the parameters. This plot shows that the methods are similar, with the RLE method again showing greatest difference in normalization (see Wilcox Test values).

# Combine normalized CPMs into a single data frame with methods as labels
normalized_combined <- data.frame(
  NONE = t(cpm_none),
  TMM = t(cpm_tmm),
  RLE = t(cpm_rle),
  UQ = t(cpm_uq))
rownames(normalized_combined) <- colnames(cpm_none) # Add sample names

# Reshape data to a long format for better handling
normalized_long <- pivot_longer(
  as.data.frame(normalized_combined),
  cols = everything(),
  names_to = "Normalization",
  values_to = "Log2_CPM")

# Split by normalization method and perform PCA
pca_results <- lapply(list(cpm_none, cpm_tmm, cpm_rle, cpm_uq), function(cpm_matrix) {
  prcomp(t(cpm_matrix), scale. = TRUE)})

# Extract PCA results and bind them together
pca_data <- do.call(rbind, lapply(1:length(pca_results), function(i) {
  data.frame(
    PC1 = pca_results[[i]]$x[, 1],
    PC2 = pca_results[[i]]$x[, 2],
    Samples = rownames(pca_results[[i]]$x),
    Normalization = c("NONE", "TMM", "RLE", "UQ")[i])}))

# Plot 5: Plot PCA results with ggplot
ggplot(pca_data, aes(x = PC1, y = PC2, color = Normalization)) +
  geom_point(size = 3) +
  facet_wrap(~Normalization, scales = "free") +
  theme_minimal() +
  labs(title = "PCA of Normalized Data",
    x = "Principal Component 1",
    y = "Principal Component 2",
    color = "Normalization" ) +
  geom_text_repel(size = 3, box.padding = 0.5, max.overlaps = Inf, aes(label = Samples))

# PCA comparison shows that there are many changes in the spread of the samples after normalization. The most notable change is in the effect on the "str212.txt" treatment. Both the TMM and UQ show the best grouping of the three 21 day treatments, removing the outlier effect of "str212.txt". These normalization methods have also produced the best clustering of all of the treatment groups.

# Continuing to Differential Expression Analysis ----
# Define the factor for the four treatment conditions
group <- factor(c("Control", "Control", "Control", "Three", "Three", "Three", "Fourteen", "Fourteen", "Fourteen", "TwentyOne", "TwentyOne", "TwentyOne"))
levels(group) <- c("Control", "Three", "Fourteen", "TwentyOne")

# Create a design matrix for the four treatments
design <- model.matrix(~ 0 + group) # Set comparisons between all levels
colnames(design) <- levels(group) # Set column names to the treatment levels

# Visually check the design matrix
design

# Define contrasts for the treatment groups vs. control, in order to detect changes in gene expression.
contrasts <- makeContrasts(
  `3Days_vs_Control` = `Three` - `Control`,
  `14Days_vs_Control` = `Fourteen` - `Control`,
  `21Days_vs_Control` = `TwentyOne` - `Control`,
  levels = design)

# Fit the a linear model to the CPM data, as a function!
fit_linear_model <- function(cpm_data, design, contrasts) {
  # Fit the linear model
  fit <- lmFit(cpm_data, design)
  fit <- eBayes(fit)
  
  # Apply contrasts
  contrast_results <- contrasts.fit(fit, contrasts)
  contrast_results <- eBayes(contrast_results)
  
  return(contrast_results)}

# Apply the function to each normalization method
none_results <- fit_linear_model(cpm_none, design, contrasts)
tmm_results <- fit_linear_model(cpm_tmm, design, contrasts)
rle_results <- fit_linear_model(cpm_rle, design, contrasts)
uq_results <- fit_linear_model(cpm_uq, design, contrasts)

# Differential expression has now been performed for each normalization method, using the limma package.

# Differential Gene Expression Results ----
# The data will be processed and visualized to compare between treatments.
# This function extracts the significant results for the normalization methods in each treatment comparison.
extract_dge_results <- function(contrast_coef, results_list, significance_threshold = 0.05) {
    results_data <- lapply(names(results_list), function(method) {
    result_table <- topTable(results_list[[method]], coef = contrast_coef, adjust = "fdr", number = Inf)
    result_df <- as.data.frame(result_table)
    result_df$Gene_ID <- rownames(result_df)  # Adds Gene ID column
    result_df$Method <- method                # Adds normalization method column
    return(result_df)})
  
  # Combines all results into a single data frame
  combined_results <- do.call(rbind, results_data)
  
  # Filters for significant genes based on the specified threshold
  significant_genes <- combined_results[combined_results$P.Value < significance_threshold, ]
  
  return(list(combined_results = combined_results, significant_genes = significant_genes))}

# Make a list for result objects
results_list <- list(
  NONE = none_results,
  TMM = tmm_results,
  RLE = rle_results,
  UQ = uq_results)

# Extract results for "3Days_vs_Control", and save them for plotting.
results_3days <- extract_dge_results("3Days_vs_Control", results_list)
combined_3days_results <- results_3days$combined_results
significant_3days_genes <- results_3days$significant_genes

# Extract results for "14Days_vs_Control", and save them for plotting.
results_14days <- extract_dge_results("14Days_vs_Control", results_list)
combined_14days_results <- results_14days$combined_results
significant_14days_genes <- results_14days$significant_genes

# Extract results for "21Days_vs_Control", and save them for plotting.
results_21days <- extract_dge_results("21Days_vs_Control", results_list)
combined_21days_results <- results_21days$combined_results
significant_21days_genes <- results_21days$significant_genes

# Add a new "Treatment" column to each combined result dataset
combined_3days_results$Treatment <- "3 Days vs. Control"
combined_14days_results$Treatment <- "14 Days vs. Control"
combined_21days_results$Treatment <- "21 Days vs. Control"

# Combine all results into one dataset
all_combined_results <- rbind(
  combined_3days_results,
  combined_14days_results,
  combined_21days_results)

# Ensure Treatment is a factor with the correct order
all_combined_results$Treatment <- factor(
  all_combined_results$Treatment,
  levels = c("3 Days vs. Control", "14 Days vs. Control", "21 Days vs. Control"))

# Add the same Treatment column to the significant genes datasets
significant_3days_genes$Treatment <- "3 Days vs. Control"
significant_14days_genes$Treatment <- "14 Days vs. Control"
significant_21days_genes$Treatment <- "21 Days vs. Control"

# Combine all significant genes into one dataset
all_significant_genes <- rbind(
  significant_3days_genes,
  significant_14days_genes,
  significant_21days_genes)

# Ensure the Treatment column has the same factor levels
all_significant_genes$Treatment <- factor(
  all_significant_genes$Treatment,
  levels = c("3 Days vs. Control", "14 Days vs. Control", "21 Days vs. Control"))

# Plot 6: Comparison across treatments of gene expression changes, with normalization methods.
volcano_facet_plot <- ggplot(all_combined_results, aes(x = logFC, y = -log10(P.Value), color = Method)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(data = all_significant_genes, aes(label = Gene_ID),
    size = 4, color = "black", box.padding = 0.5, max.overlaps = 20) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Volcano Plots for Normalization Methods across Treatment Methods",
    x = "Log2 Fold Change", y = "-Log10 p-value") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") +
  facet_wrap(~Treatment, ncol = 3, scales = "fixed")

# Print the plot
volcano_facet_plot

# The final plot shows that the different normalization methods had significant impacts on the detection of signficiance in changes in gene expression. For example, Tian et al.(2021) found that in the 21 day treatment, all four genes' expression was signficantly altered. In my analysis, the Upper Quantile method was the only normalization method to replicate their results. Interestingly, no normalization method out performed TMM and RLE in terms of replicating results.

