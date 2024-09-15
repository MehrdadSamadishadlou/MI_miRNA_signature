# Load required libraries
library(GEOquery)
library(Biobase)
library(pheatmap)
library(ggplot2)
library(limma)
library(gplots)
library(stringr)

#### Reading and preparing data ####
# Retrieve dataset from GEO
Geneset <- getGEO("GSE61741", GSEMatrix =TRUE, AnnotGPL=TRUE)
geneset <- Geneset[[1]]

# Extract expression data and apply log2 transformation
exp <- exprs(geneset)
exp <- log2(exp + 1)

# Extract and clean up disease information
disease <- geneset@phenoData@data[["characteristics_ch1"]]
disease <- str_remove(disease, "disease: ")
disease <- str_replace_all(disease, " ", "_")
disease <- str_replace_all(disease, "-", "_")

# Add disease information to expression data
exp <- rbind(exp, disease)
Type <- as.character(exp["disease", ])
exp <- exp[-849,]

# Convert expression values to numeric
for (col in names(exp)) {
  exp[[col]] <- as.numeric(exp[[col]])
}

# Load color palette
library(rcartocolor)
prism_colors <- carto_pal(name = "ArmyRose")


# Separate normal and myocardial_infarction samples
mi <- exp[, disease == "myocardial_infarction"]
normal <- exp[, disease == "normal"]
mi_normal2 <- cbind(mi, normal)
mi_normal <- apply(mi_normal2,  FUN = as.numeric, MARGIN = 2)
rownames(mi_normal) <- rownames(mi_normal2)

#### PCA ####
# Perform PCA
pc <- prcomp(mi_normal)

# Save PCA plot
pdf("Results/MI-LLI.pdf")
plot(pc)
dev.off()

# Prepare data for PCA visualization
Disease <- c(rep("myocardial_infarction", 62), rep("normal", 94))
pcr <- data.frame(pc$r[, 1:3] , Group = Disease)

# Create and save PCA plot
pdf("Results/PCA_mi_LLI.pdf", width = 20, height = 13)
ggplot(pcr, aes(PC1 , PC2 , color= Disease)) + 
  geom_point(size = 7, alpha = 0.5) +
  ggtitle("") +
  scale_color_manual(breaks = levels(as.factor(Disease)),
                     values = c("red", "darkgreen")) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=25, face="bold"),
        legend.title = element_text(size = 18, face="bold"),
        legend.text = element_text(size = 14),
        plot.title = element_text(color="black", size=28, face="bold.italic", hjust = 0.5))
dev.off()

#### Differential Expression Analysis with Limma ####
# Prepare design matrix
sample_status <- factor(Disease)
Design <- model.matrix(~ sample_status + 0, as.data.frame(mi_normal))
colnames(Design) <- levels(sample_status)
colSums(Design) #For checking correctness of the procedure.

# Fit linear model
fit <- lmFit(mi_normal, Design)
cont.matrix <- makeContrasts(myocardial_infarction-normal, levels=Design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

# Get top differentially expressed genes
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("adj.P.Val", "logFC"))
write.csv(tT, "Results/mi_normal_statistical_data.csv", row.names=T, quote = F)

# Filter and save up-regulated genes in MI and normal samples
MI.up <- subset(tT, logFC>1 & adj.P.Val<0.05)
write.csv(MI.up, file="Results/MI_up.csv", quote = F, row.names = T)
normal.up <- subset(tT, logFC < -1 & adj.P.Val<0.05)
write.csv(normal.up, file="Results/normal_up.csv", quote = F, row.names = T)

# Combine and save all differentially expressed genes
DEGs <- rbind(MI.up,normal.up)
write.csv(DEGs, file="Results/DEGs_LLI.csv", quote = F, row.names = T)

# Extracting DEMs expression
dems_exp <- mi_normal[rownames(DEGs), ]
dems_exp <- rbind(dems_exp, c(rep("MI", 62), rep("Healthy", 94)))
write.csv(dems_exp, file="Results/DEMs_Expression.csv", quote = F, row.names = T)


#### LASSO ####
library(glmnet)

# Prepare data for LASS)
disease <- dems_exp[101, ]
disease <- data.frame(disease)
y <- model.matrix(~disease - 1, data = disease)
X <- t(apply(dems_exp[-101,], c(1, 2), as.numeric))

# Perform cross-validated LASSO
cvfit <- cv.glmnet(X, y, alpha=1, family="binomial", weights = NULL)

# Plot and save cross-validation results
pdf("lass_cv_fit_class.pdf", width = 6, height = 4.5, family = "Times")
plot(cvfit)
dev.off()

# Get optimal lambda and fit LASSO model
lambda_min <- cvfit$lambda.min
fit <- glmnet(X, y, family="binomial", alpha=1, weights = NULL, lambda = lambda_min)

# Extract and save LASSO-selected features
lasso_mirs <- fit[["beta"]]@Dimnames[[1]][(fit[["beta"]]@i + 1)]
write.table(lasso_mirs, "lasso_selected_class.txt", quote = F)
lass_mir_exp <- mi_normal[c(lasso_mirs, "disease"), ]
write.csv(lass_mir_exp, "lasso_mirs_exp.csv", quote = FALSE)


#### Volcano Plot ####
library(tidyverse) 
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations


# Prepare data for volcano plot
df <- tT %>%
  mutate(diffexp=ifelse(logFC>1 & adj.P.Val<0.05, "Up-regulated", "Not significant")) %>%
  mutate(diffexp=ifelse(logFC < -1 & adj.P.Val<0.05, "Down-regulated", diffexp)) %>% 
  mutate(label=ifelse(`...1` %in% lasso_mirs$V1, `...1`, NA))

# Create and save volcano plot
pdf("Results/volcano.pdf", width = 17, height = 13, family = "Times")
ggplot(data = df, aes(x = logFC, y = -log10(adj.P.Val), col = diffexp, label=label)) +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +
  geom_point(size = 5, alpha = 0.6) +
  scale_color_manual(values = c("orange", "grey", "darkolivegreen")) +
  labs(color = "", x = expression(bold("log"[2]*"FC")), y = expression(bold("-log"[10]*" adj. p-value"))) +
  geom_text_repel(data          = subset(df, logFC < -1),
                  nudge_x       = -2.5 - subset(df, logFC < -1)$logFC,
                 segment.size  = 0.5,
                 segment.color = "grey35",
                 direction     = "y",
                 hjust         = 1,
                 size          = 10,
                 color         = "darkorange4",
                 box.padding   = 1) +
  geom_text_repel(data          = subset(df, logFC > 1),
                  nudge_x       = 3 - subset(df, logFC > 1)$logFC,
                 segment.size  = 0.5,
                 segment.color = "grey35",
                 direction     = "y",
                 hjust         = 1,
                 size          = 10,
                 color         = "darkgreen",
                 box.padding   = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 33),
        legend.text = element_text(size = 31),
  )
dev.off()


#### Venn Diagram ####

# Load required libraries
library(ggvenn)
library(ggplot2)
library(stringr)

# Read LASSO-selected miRNAs
lasso <- rbind(lasso_mirs, data.frame(V1=seq(1, 219)))
row.names(lasso) <- lasso$V1
lasso$V2 <- lasso$V1
colnames(lasso) <- c("value", "V1")

# Mark LASSO-selected miRNAs
lasso[!str_detect(lasso$V1, "hsa"), "V1"] = FALSE
lasso[str_detect(lasso$V1, "hsa"), "V1"] = TRUE

# Read GPL5175 miRNAs
gpl <- read.table("GPL5175mirs.txt", check.names = F)

# Combine LASSO and GPL5175 data
df <- cbind(lasso, gpl)
colnames(df) <- c("value", "LASSO", "GPL5175")

# Mark GPL5175 miRNAs
df[!str_detect(df$GPL5175, "no"), "GPL5175"] = TRUE
df[str_detect(df$GPL5175, "no"), "GPL5175"] = FALSE

# Convert to logical type and rename columns
df$LASSO <- as.logical(df$LASSO)
df$GPL5175 <- as.logical(df$GPL5175)
colnames(df) <- c("value", "LASSO-selected miRNAs", "GPL5175 miRNAs")

# Find common miRNAs between LASSO and GPL5175
common_mirs <- df[df$`GPL5175 miRNAs` == TRUE, ]
common_mirs <- common_mirs[common_mirs$`LASSO-selected miRNAs` == TRUE, ]
common_mirs <- common_mirs$value

# Save common miRNAs
write.csv(common_mirs, "common_mir_10.csv", quote = F)

# Extract expression data for common miRNAs
common_10_exp <- dems_exp[c(common_mirs, "disease"), ]
write.csv(common_10_exp, "common_mirs_exp.csv", quote = F)

# Format common miRNA names for display
common_mirs <- str_remove_all(common_mirs, "hsa-miR-")
common_mirs <- c(common_mirs, "miRNAs:")

# Create and save Venn diagram
pdf("venn.pdf", width = 15, height = 15, fonts = "Times")
ggplot(df, aes(A=`LASSO-selected miRNAs`, B=`GPL5175 miRNAs`)) + 
  geom_venn(show_percentage = FALSE,
            fill_color = c("#0073C2FF", "#CD534CFF"),
            fill_alpha = 0.7,
            set_name_size = 11,
            text_size = 0 ) +
  # Add text annotations for common miRNAs
  annotate("text", x = 0, y = seq(-0.52, 0.48, length.out = 11), label = common_mirs, size = 8) +
  # Add count annotations
  annotate("text", x = -0.8, y = 0, label = "11", size = 14) +
  annotate("text", x = 0.8, y = 0, label = "219", size = 14) +
  # Customize theme
  theme_bw() +
  theme(text=element_text(family="Times"), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank(), axis.title = element_blank())
dev.off()


#### Extracting LASSO-selected miRs expression from test set####

# Download and process GSE29532 dataset
Geneset <- getGEO("GSE29532", GSEMatrix =TRUE, AnnotGPL=TRUE)
geneset <- Geneset[[1]]

# Extract sample information
MI_samples <- geneset@phenoData@data[["geo_accession"]][geneset@phenoData@data[["time point:ch1"]] == "Time 1"]
cnt_sample <- geneset@phenoData@data[["geo_accession"]][geneset@phenoData@data[["time point:ch1"]] == "Control"]
samples <- c(MI_samples, cnt_sample)

# Load additional libraries for microarray processing
library(oligo)
library(affy)
library(pd.huex.1.0.st.v2)
library(huex.1.0.st.v2frmavecs)
library(frma)
data("huex.1.0.st.v2frmavecs")

# Process raw CEL files
setwd("RAW")
celFiles = list.celfiles()
affyRaw = read.celfiles(celFiles)
eset <- frma(affyRaw, input.vecs=huex.1.0.st.v2frmavecs, 
             target = "full", summarize="robust_weighted_average", 
             background = "rma", normalize = "quantile")

# Extract expression data
exp <- exprs(eset)
colnames(exp) <- samples


# Extract miRNAs of interest
ids <- as.character(lasso$gpl5175probeID)
exp_mirs <- as.matrix(exp[ids,])
row.names(exp_mirs) <- mirs_name$x

# Add sample status information
exp_mirs <- rbind(exp_mirs, status = c(rep("MI", 8), rep("Healthy", 6)))

# Save miRNA expression data
write.table(exp_mirs, "GSE29532_mirs_expression.txt", quote = F, sep = "\t")

