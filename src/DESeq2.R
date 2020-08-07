#!/usr/bin/env R

#RNA-seq analysis using DESeq2 and a counts table from featureCounts
#Emory EICC 2020-08-06
#Sean McKenzie, modified from script by Jessica Randal

#Loads required packages and sets enviroment variables
pacman::p_load("readr", "dplyr", "knitr", "DESeq2", "vsn", "ggplot2", 
               "pheatmap", "EnhancedVolcano", "apeglm", "tinytex","ggbeeswarm")
theme_set(theme_minimal())

files <- list(
  counts = file.path("../input_data/feature_counts.tab"),
  sampledata = file.path("../input_data/Barbian_metadata.csv"),
  
  comp = file.path("../output/DEresults.csv"),
  meannormcounts = file.path("../output/meannormcounts.csv"),
  vstcounts = file.path("../output/vstcounts.csv"),
  
  pca = file.path("../output/graphs/PCA.png"),
  pca_sex = file.path("../output/graphs/PCA_sex.png"),
  pca_animal = file.path("../output/graphs/PCA_animal.png"),
  countplot = file.path("../output/graphs/countsplot.png"),
  heat = file.path("../output/graphs/heatmap.png"),
  volc = file.path("../output/graphs/volcanoplot.png"),
  go = file.path("../output/GOresults.csv")
  
)


#Reads in count data and sets column names (column name settings will likely need to be adjusted ba)
countdata <- read.table(files$counts, sep = "\t", header = TRUE, row.names="Geneid")
sample_cols <- colnames(countdata)[grep("Aligned.bam",colnames(countdata))]
countdata <- countdata[,sample_cols]
sample_cols <- sub("_trim.*","",sub(".*star.","",sample_cols))
colnames(countdata) <- sample_cols

sampledata <- read.table(files$sampledata,sep=",",header=TRUE,row.names = "SampleID")
sampledata$Prenatal_exposure <- as.factor(sampledata$Prenatal_exposure)
sampledata$animalID <- rownames(sampledata)

countdata <- countdata[,rownames(sampledata)]

stopifnot(rownames(sampledata) %in% colnames(countdata))

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sampledata,
                              design = ~ Prenatal_exposure)
dds$Prenatal_exposure <- relevel(dds$Prenatal_exposure, ref = "Control")
#### exploratory PCA ###
# perform variance stabilizing transformation for PCA and Heatmaps
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# calculate the variance for each gene
rv <- rowVars(assay(vsd))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]

# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

# the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

# assembly the data for each plot
d_group <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], 
                      group=sampledata$Prenatal_exposure, 
                      name=colnames(vsd))

# look at pca plot overall
(pca1 <- ggplot(d_group, aes_string(x="PC1", y="PC2", color = "group")) + 
    geom_point(size = 3) + 
    ggtitle("Principal Components Analysis (by group)") +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    scale_color_manual(values=c("#2c7bb6", "#b2182b")))

pca1_data = plotPCA(vsd, intgroup=c("Prenatal_exposure"), returnData = TRUE)

#export PCA
pca_all <- png(files$pca, width=600, height=500)
ggplot(d_group, aes_string(x="PC1", y="PC2", color = "group")) + 
  geom_point(size = 3) + 
  ggtitle("Principal Components Analysis (by group)") +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  scale_color_manual(values=c("#2c7bb6", "#b2182b"))
dev.off()

### save transformed data as df for heatmaps ###
vsd_df <-as.data.frame(assay(vsd)) %>%
  mutate(ID = as.factor(row.names(vsd)))

# export vst transformed count data 
write_excel_csv(vsd_df, files$vstcounts)

### analysis ###
dds <- DESeq(dds)
res <- results(dds, name="Prenatal_exposure_Butyrate_vs_Control",
               alpha=0.05)
summary(res)

LFC_dds <-as.data.frame(lfcShrink(dds, 
                                  coef="Prenatal_exposure_Butyrate_vs_Control", 
                                  type = "apeglm", svalue = TRUE)) %>%
  arrange(log2FoldChange)

res_df <- as.data.frame(res) %>%
  mutate(ID = as.factor(row.names(res))) %>%
  arrange(log2FoldChange) %>%
  mutate(svals = LFC_dds$svalue) %>%
  arrange(padj)
row.names(res_df) <- res_df$ID


topGene <- rownames(res)[which.min(res$pvalue)]

# see response in aggregate
cp_agg <- plotCounts(dds, gene = topGene, intgroup=c("Prenatal_exposure", "animalID"), 
                     returnData = TRUE)
ggplot(cp_agg, aes(x = Prenatal_exposure, y = count, fill = Prenatal_exposure)) +
  scale_y_log10() + 
  geom_dotplot(binaxis='y', stackdir='center') +
  ggtitle(paste("Normalized counts of expression of",topGene)) +
  scale_fill_manual(values=c("#2c7bb6", "#b2182b"))

#export
cp1 <- png(files$countplot, width=600, height=400)
ggplot(cp_agg, aes(x = Prenatal_exposure, y = count, fill = Prenatal_exposure, color = Prenatal_exposure)) +
  scale_y_log10() + 
  geom_dotplot(binaxis='y', stackdir='center') +
  ggtitle(ggtitle(paste("Normalized counts of expression of",topGene))) +
  scale_fill_manual(values=c( "#b2182b","#2c7bb6")) +
  scale_color_manual(values=c( "#b2182b","#2c7bb6"))
dev.off()

# extract and export mean normalized counts
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))%>%
  mutate (ID = row.names(counts(dds)))
write_excel_csv(norm_counts, files$meannormcounts)

full_results <- res_df

full_results <- full_results %>%
  dplyr::select(-ends_with("y")) %>%
  dplyr::select("ID", "baseMean", "log2FoldChange", "lfcSE", 
                "stat", "pvalue", "padj", "svals") %>%
  arrange(padj)

# export results
write_excel_csv(full_results, files$comp)

### heatmap ###

# select only top n genes
res_heat <- as.data.frame(res)
select <- order(res_heat$padj, decreasing = FALSE)[1:20]

#select vars of interest
heat <-as.data.frame(colData(dds)[c("Prenatal_exposure")])

#specify that reference level is on the left
callback <- function(hc, mat){
  sv <- svd(mat)$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

#specify color palette of choice
RdBu11 <-c("#313695", "#4575b4", "#74add1", "#abd9e9", 
           "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", 
           "#f46d43", "#d73027", "#a50026")

#plot
(ht = pheatmap(assay(vsd)[select,],
               cluster_rows = FALSE,
               show_rownames = TRUE,
               cluster_cols = TRUE,
               annotation_col = heat,
               color = RdBu11,
               clustering_callback = callback, 
               width = 1))

#export
ggsave(files$heat, plot=ht, dpi=600)

# color palette of choice
RdBu4 <- c("#92c5de", "#0571b0", "#f4a582", "#ca0020")

#volcanoplot
(volc<- EnhancedVolcano(res_df,
                        lab = res_df$ID,
                        x = "log2FoldChange",
                        y = "padj",
                        xlim = c(-12, 12),
                        title= NULL,
                        col = RdBu4,
                        subtitle = "Log2 fold-change vs adj. p-value",
                        pCutoff = 0.01,
                        legendPosition = "bottom",
                        ))

#export
volc <- png(files$volc, width=700, height=600)
EnhancedVolcano(res_df,
                lab = res_df$ID,
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-12, 12),
                title= NULL,
                col = RdBu4,
                subtitle = "Log2 fold-change vs adj. p-value",
                pCutoff = 0.01,
                legendPosition = "bottom",
                legendLabels = c("NS", "Log2 fold-change", "adj P-value",
                         "adj P-value & Log2 fold-change"))
dev.off()
