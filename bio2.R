# Libraries
library(clusterProfiler)
library(DESeq2)
library(survival)
library(tidyverse)
library(survminer)

##-----Step 2-----------------------------------------------------------------------------------
# Untar folder
untar("C:/Users/willi/Documents/School/Bio/brca_tcga_pan_can_atlas_2018.tar.gz")

# Set working directory
setwd("C:/Users/willi/Documents/School/Bio/brca_tcga_pan_can_atlas_2018")

##-----Steps 3, 4 & 5---------------------------------------------------------------------------
# Read in each dataset
data_rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")
data_patient = read.delim("data_clinical_patient.txt")
data_aberrations = read.delim("data_cna.txt")

# Get index of ERBB2 gene
erbb2_indx = which(data_aberrations[,1] == 'ERBB2')

##-----Step 6-----------------------------------------------------------------------------------
# Remove duplicates from RNA sequencing dataset
keep = !duplicated(data_rnaseq[,1])
data_rnaseq = data_rnaseq[keep,]

# Set row labels as hugo value of genes
rownames(data_rnaseq)  = data_rnaseq[,1]

# Remove info rows and get list of patients
aberrations_patients = as.matrix(data_aberrations[,-c(1,2)])
aberrations_list = colnames(aberrations_patients)

rnaseq_patients = as.matrix(data_rnaseq[,-c(1,2)])
rnaseq_list = colnames(rnaseq_patients)

# Create intersection patient list
intersect_list = intersect(aberrations_list, rnaseq_list)


# Filter all datasets, keeping only the patients present in both
aberrations_filtered <- data_aberrations[, names(data_aberrations) %in% intersect_list]
rnaseq_filtered <- data_rnaseq[, names(data_rnaseq) %in% intersect_list]


##-----Step 7-----------------------------------------------------------------------------------
# Create metadata from these columns
metadata = matrix(0, dim(aberrations_filtered)[2],1)

for (i in 1:dim(aberrations_filtered)[2]){
  
  if(aberrations_filtered[erbb2_indx, i] > 0){
    metadata[i,1] = 1
  }
  else {
    metadata[i,1] = 0
  }
}
metadata[is.na(metadata)] = 0

# Label metadata column
colnames(metadata) = "Ampified_ERBB2"

# Make histogram showing breakdown on HER2+ amplified vs not
hist(as.numeric(metadata), main="HER2 Amplified - Proportion", xlab = "Not Amplified v Amplified")

##-----Step 8-----------------------------------------------------------------------------------

# Run Deseq
dds <- DESeqDataSetFromMatrix(countData = round(rnaseq_filtered),
                              colData = metadata,
                              design = ~ Ampified_ERBB2)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


# Normalize
dds <- DESeq(dds)

# Get Results
res <- results(dds)

##-----Step 9-----------------------------------------------------------------------------------
# Order by Log2Fold Change
resOrdered <- res[order(res$padj),]

# Significantly Differentially Expressed
signif = which(res$padj<0.05)
results_signif = res[signif,]

# Get top differentially expressed genes
results_signif_ordered <- results_signif[order(results_signif$log2FoldChange),]

# Top up & down regulated genes in HER2+ amplified
up_regulated_genes = head(rev(results_signif_ordered), 10)
down_regulated_genes = head(results_signif_ordered, 10)

##-----Step 10-----------------------------------------------------------------------------------

# For Pathway Enrichment we need Entrez IDs
up_regulated_signif = results_signif[results_signif$log2FoldChange > 0,]
signif_rna = data_rnaseq[data_rnaseq$Hugo_Symbol %in% row.names(up_regulated_signif) ,]
signif_entrez_ids = signif_rna$Entrez_Gene_Id

# Pathway Enrichment
up_paths = enrichKEGG(gene = signif_entrez_ids, organism = 'hsa', pAdjustMethod = 'BH', pvalueCutoff = 0.05)
head(up_paths)


# Show percentage of genes kegg actually picked up on
bitr_kegg(signif_entrez_ids, fromType = "kegg", toType = "Module", organism = "hsa")

# Plot enrichment
barplot(up_paths)

##-----Steps 11 & 12---------------------------------------------------------------------------------

# Transform the data to visualize
rld <- vst(dds, blind=FALSE)

# Do Principal Components Analysis
pc = prcomp(assay(rld), scale. = FALSE)

# See how much variance PCs account for
head(summary(pc)$importance[2,])

# Get PC1 & PC2
pc1 <- pc$rotation[,1]
pc2 <- pc$rotation[,2]
pc_comb <- data.frame(pc1, pc2)

# Plot first 2 principle components

ggplot(pc_comb, aes(x = pc1, y = pc2, color = factor(metadata))) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2",
       color = "Her2+ Amplified?")

####################extra portion#######################################
##-----Step 13-----------------------------------------------------------------------------------

# Perform clustering on the first two principal components
kmeans_result <- kmeans(pc_comb, centers = 4)

# Separate data by clusters
cluster1 <- rnaseq_filtered[kmeans_result$cluster == 1]
cluster2 <- rnaseq_filtered[kmeans_result$cluster == 2]
cluster3 <- rnaseq_filtered[kmeans_result$cluster == 3]
cluster4 <- rnaseq_filtered[kmeans_result$cluster == 4]

# Generate PCA for each cluster
pc_cl1 = prcomp(cluster1, scale. = FALSE)
pc_cl2 = prcomp(cluster2, scale. = FALSE)
pc_cl3 = prcomp(cluster3, scale. = FALSE)
pc_cl4 = prcomp(cluster4, scale. = FALSE)

# Extract PCs from PCAs for each cluster
pc1_cl1 <- pc_cl1$rotation[,1]
pc2_cl1 <- pc_cl1$rotation[,2]
pc_comb_cl1 <- data.frame(pc1_cl1, pc2_cl1)

pc1_cl2 <- pc_cl2$rotation[,1]
pc2_cl2 <- pc_cl2$rotation[,2]
pc_comb_cl2 <- data.frame(pc1_cl2, pc2_cl2)

pc1_cl3 <- pc_cl3$rotation[,1]
pc2_cl3 <- pc_cl3$rotation[,2]
pc_comb_cl3 <- data.frame(pc1_cl3, pc2_cl3)

pc1_cl4 <- pc_cl4$rotation[,1]
pc2_cl4 <- pc_cl4$rotation[,2]
pc_comb_cl4 <- data.frame(pc1_cl4, pc2_cl4)

## Plot the results
# Cluster 1,2,3 & 4
ggplot(pc_comb_cl1, aes(x = pc1_cl1, y = pc2_cl1, color = "red")) +
  geom_point() +
  scale_color_manual(values = "red") +
  theme(legend.position="none") +
  labs(title = "Cluster 1",
       x = "PC1 for Cluster 1",
       y = "PC2 for Cluster 1")

ggplot(pc_comb_cl2, aes(x = pc1_cl2, y = pc2_cl2, color="blue")) +
  geom_point() +
  scale_color_manual(values = "blue") +
  theme(legend.position="none") +
  labs(title = "Cluster 2",
       x = "PC1 for Cluster 2",
       y = "PC2 for Cluster 2")

ggplot(pc_comb_cl3, aes(x = pc1_cl3, y = pc2_cl3, color="green")) +
  geom_point() +
  scale_color_manual(values = "green") +
  theme(legend.position="none") +
  labs(title = "Cluster 3",
       x = "PC1 for Cluster 3",
       y = "PC2 for Cluster 3")

ggplot(pc_comb_cl4, aes(x = pc1_cl4, y = pc2_cl4, color="orange")) +
  geom_point() +
  scale_color_manual(values = "orange") +
  theme(legend.position="none") +
  labs(title = "Cluster 4",
       x = "PC1 for Cluster 4",
       y = "PC2 for Cluster 4")


##-----Step 14 & 15---------------------------------------------------------------------------------
## Create Survival model
rld_df_transpose = data.frame(t(assay(rld)))

# Change patient naming to match other data
data_patient$X.Patient.Identifier <- gsub("-", ".", data_patient$X.Patient.Identifier)
data_patient$X.Patient.Identifier <- paste0(data_patient$X.Patient.Identifier, ".01")

# Filter all names not in VST values dataframe
patient_survival <- data_patient[data_patient$X.Patient.Identifier %in% row.names(rld_df_transpose), ]

# Change survival status to number
patient_survival$Overall.Survival.Status <- ifelse(as.character(patient_survival$Overall.Survival.Status) 
                                                   == "0:LIVING", 0, 1)
# Create survival object
surv <- Surv(time = as.numeric(patient_survival$Overall.Survival..Months), event=patient_survival$Overall.Survival.Status)

# Apply cox regression to each gene and place in list
cox_regression_list <- lapply(rld_df_transpose[, -1], function(gene_val) {coxph(surv ~ gene_val, data = rld_df_transpose)})

# Extract coefficients from list and place in separate vector 
cox_coefficient_list <- sapply(cox_regression_list, function(gene_val) {gene_val$coefficients[1]})

# Tidying and sort by lowest coefficients
names(cox_coefficient_list) <- gsub(".gene_val", "", names(cox_coefficient_list))
survival_genes <- sort(cox_coefficient_list, decreasing = FALSE)[1:5]

#Print top 5 genes for indicating survival with their coefficients
print(survival_genes)

# Extra Chart to map survival rate for amplified group
fit <- survfit(surv ~metadata,data = patient_survival)
ggsurvplot(fit, data = patient_survival)
