#!/usr/bin/env R

#RNA-seq analysis using DESeq2 and a counts table from featureCounts
#Emory EICC 2020-08-06
#Sean McKenzie, modified from script by Jessica Randall

#Loads required packages and sets enviroment variables
pacman::p_load("here", "janitor", "tidyverse", "assertr", 
               "DESeq2", "vsn", "pheatmap", "EnhancedVolcano", 
               "apeglm")

# links to documentation for here, janitor, and tidyverse 
#   dplyr, ggplot2 and readr are all in the tidyverse package

theme_set(theme_minimal())

files <- list(
  counts = here("RJO12354/rna_seq/input_data/feature_counts.tab"),
  sampledata = here("RJO12354/rna_seq/input_data/Barbian_metadata.csv"),
  
  comp = here("RJO12354/rna_seq/output/DEresults.csv"),
  meannormcounts = here("RJO12354/rna_seq/output/meannormcounts.csv"),
  vstcounts = here("RJO12354/rna_seq/output/vstcounts.csv"),
  
  pca = here("RJO12354/rna_seq/output/graphs/PCA.png"),
  pca_sex = here("RJO12354/rna_seq/output/graphs/PCA_sex.png"),
  pca_animal = here("RJO12354/rna_seq/output/graphs/PCA_animal.png"),
  countplot = here("RJO12354/rna_seq/output/graphs/countsplot.png"),
  heat = here("RJO12354/rna_seq/output/graphs/heatmap.png"),
  volc = here("RJO12354/rna_seq/output/graphs/volcanoplot.png"),
  go = here("RJO12354/rna_seq/output/GOresults.csv")
  
)

#Reads in count data and set column names

countdata <- as.data.frame(read.table(files$counts, sep = "\t", header = TRUE, row.names="Geneid")) %>%
  clean_names() %>%
  select(-c(chr, start, end, strand, length))

# change column names to be just sample ids
names(countdata) <- c("rjo10", "rjo11", "rjo12", "rjo13", "rjo14", "rjo15", 
                      "rjo16", "rjo17", "rjo18", "rjo19", "rjo02", "rjo03", 
                      "rjo05", "rjo6", "rj07", "rjo8", "rjo09")

sample_cols <- names(countdata)

# there's a nifty function called gsub that works similarly to sub too, 
# I usually use it within a mutate statement to make a new variable  in a df
# sample_cols <- sub("_trim.*","",sub(".*star.","",sample_cols))

sampledata <- as.data.frame(read_csv(files$sampledata, col_names =TRUE)) %>%
  clean_names()

# we need to create sampledata before we can give it rownames

row.names(sampledata) <- c("rjo02", "rjo03", "rjo05", "rjo6", "rj07", "rjo8",
                           "rjo09", "rjo10", "rjo11", "rjo12", "rjo13", "rjo14",
                           "rjo15", "rjo16", "rjo17", "rjo18", "rjo19")
sampledata <- sampledata %>%
  clean_names() %>%
  mutate(prenat_exp = as.factor(prenatal_exposure),
         animal_id = as.character(rownames(sampledata)))

row.names(sampledata) <- sampledata$animal_id

# put count data in same order as sample data
countdata <- countdata[,rownames(sampledata)]

# changed stopifnot to a verify statement within assertr to stop the program 
# if this is not true

sampledata <- sampledata %>% 
  verify(rownames(sampledata) %in% colnames(countdata))

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sampledata,
                              design = ~ prenat_exp)

dds$prenat_exp <- relevel(dds$prenat_exp, ref = "Control")

#### exploratory PCA ###
# perform variance stabilizing transformation for PCA and Heatmaps
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# love that you did your own PCA, it's much more fun to customize it than
# use DESeq's function

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
                      group=sampledata$prenat_exp, 
                      name=colnames(vsd))

# look at pca plot overall
(pca1 <- ggplot(d_group, aes_string(x="PC1", y="PC2", color = "group")) + 
    geom_point(size = 3) + 
    ggtitle("Principal Components Analysis (by group)") +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    scale_color_manual(values=c("#2c7bb6", "#b2182b")))

pca1_data = plotPCA(vsd, intgroup=c("prenat_exp"), returnData = TRUE)

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

#FIXME: add unit test

# export vst transformed count data 
write_excel_csv(vsd_df, files$vstcounts)

### analysis ###
# renames dds to dds_obj just because it makes me uneasy to use the same name
# for a new object, easy to use the wrong one by accident 

dds_obj <- DESeq(dds)

res <- results(dds_obj, name="prenat_exp_Butyrate_vs_Control", alpha=0.05)

summary(res)

# results with FDR < 0.05 #

# out of 20626 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 11, 0.053%
# LFC < 0 (down)     : 26, 0.13%
# outliers [1]       : 0, 0%
# low counts [2]     : 3189, 15%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# results with FDR < 0.1 #

# out of 20626 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 26, 0.13%
# LFC < 0 (down)     : 39, 0.19%
# outliers [1]       : 0, 0%
# low counts [2]     : 3585, 17%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

LFC_dds <-as.data.frame(lfcShrink(dds_obj, 
                                  coef="prenat_exp_Butyrate_vs_Control", 
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
cp_agg <- plotCounts(dds_obj, gene = topGene, intgroup=c("prenat_exp", "animal_id"), 
                     returnData = TRUE)
cp_agg %>%
  ggplot(aes(x = prenat_exp, y = count, fill = prenat_exp)) +
  scale_y_log10() + 
  geom_dotplot(binaxis='y', stackdir='center') +
  ggtitle(paste("Normalized counts of expression of",topGene)) +
  scale_fill_manual(values=c("#2c7bb6", "#b2182b"))

#export
cp1 <- png(files$countplot, width=600, height=400)
ggplot(cp_agg, aes(x = prenat_exp, y = count, fill = prenat_exp, color = prenat_exp)) +
  scale_y_log10() + 
  geom_dotplot(binaxis='y', stackdir='center') +
  ggtitle(ggtitle(paste("Normalized counts of expression of",topGene))) +
  scale_fill_manual(values=c( "#b2182b","#2c7bb6")) +
  scale_color_manual(values=c( "#b2182b","#2c7bb6"))
dev.off()

# extract and export mean normalized counts
norm_counts <- as.data.frame(counts(dds_obj, normalized = TRUE)) %>%
  mutate(ID = row.names(counts(dds_obj)))

# added unit tests to check that the output had the expected dimensions before
# we export

norm_counts <- norm_counts %>%
  verify(ncol(norm_counts) == 18 & nrow(norm_counts) == 24421) %>%
  write_excel_csv(files$meannormcounts)

# arrange full results file

full_results <- res_df  %>%
  dplyr::select(-ends_with("y")) %>%
  dplyr::select("ID", "baseMean", "log2FoldChange", "lfcSE", 
                "stat", "pvalue", "padj", "svals") %>%
  arrange(padj)

# added unit tests to check that the output had the expected dimensions before
# we export

full_results <- full_results %>%
  verify(ncol(full_results) == 8 & nrow(full_results) == 24421) %>%
  write_excel_csv(files$comp)

### heatmap ###

# select only top n genes
res_heat <- as.data.frame(res)
select <- order(res_heat$padj, decreasing = FALSE)[1:20]

#select vars of interest
heat <-as.data.frame(colData(dds_obj)[c("prenat_exp")])

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
                        xlim = c(-6, 6),
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
                xlim = c(-6, 6),
                title= NULL,
                col = RdBu4,
                subtitle = "Log2 fold-change vs adj. p-value",
                pCutoff = 0.01,
                legendPosition = "bottom",
                legendLabels = c("NS", "Log2 fold-change", "adj P-value",
                         "adj P-value & Log2 fold-change"))
dev.off()

