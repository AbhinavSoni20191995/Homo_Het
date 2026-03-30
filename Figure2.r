#Bulk rna seq PBS wt
#Deseq2 processing----
#Remove duplicates as DEseq2 does not accept non duplicated genes
duplicated(bfx1419_various_countTables_de_expl_copy$Gene_Symbol,  )
sum(duplicated(bfx1419_various_countTables_de_expl_copy$Gene_Symbol,  )) #we can see how many are duplicates
bfx1419_various_countTables_de_expl_duplicate<-bfx1419_various_countTables_de_expl_copy[duplicated(bfx1419_various_countTables_de_expl_copy$Gene_Symbol,  ), ]#details about duplicate
bfx1419_various_countTables_de_expl_copy<-bfx1419_various_countTables_de_expl_copy[!(duplicated(bfx1419_various_countTables_de_expl_copy$Gene_Symbol,  )), ]#remove duplicate column from count table

bfx1419_various_countTables_de_expl_copy<-as.data.frame(bfx1419_various_countTables_de_expl_copy)
row.names(bfx1419_various_countTables_de_expl_copy)<-bfx1419_various_countTables_de_expl_copy$Gene_Symbol
bfx1419_various_countTables_de_expl_final <- bfx1419_various_countTables_de_expl_copy[,-1]

#Upload the meta data file
Meta_combined$Condition<-as.factor(Meta_combined$Condition)

#All the conditions together, then we can also run separately
dds_hao <- DESeqDataSetFromMatrix(countData = bfx1419_various_countTables_de_expl_final,
                                  colData = Meta_combined,
                                  design = ~ Condition)

dds_hao<-DESeq(dds_hao)
rld_hao <-rlog(dds_hao, blind = FALSE)
rld_hao_df<-as.data.frame(assay(rld_hao))
rld_hao_df_wt_lps<-rld_hao_df %>% dplyr::select(PBS_WT_1,PBS_WT_2,PBS_WT_3,PBS_WT_4,LPS_WT_1,LPS_WT_2,LPS_WT_3,LPS_WT_4)

mod_mat <- model.matrix(design(dds_hao), colData(dds_hao))        
mod_mat #see how it looks like
LPS_E269R <- colMeans(mod_mat[dds_hao$Condition == "LPS_E269R", ])
LPS_WT <- colMeans(mod_mat[dds_hao$Condition == "LPS_WT", ])
PBS_E269R <- colMeans(mod_mat[dds_hao$Condition == "PBS_E269R", ])
PBS_WT <- colMeans(mod_mat[dds_hao$Condition == "PBS_WT", ])
?model.matrix
#Run differential analysis comparing different condition
LPS_wt_vs_PBS_wt <- results(dds_hao, contrast = LPS_WT - PBS_WT, alpha = 0.05)
summary(LPS_wt_vs_PBS_wt)
LPS_wt_vs_PBS_wt.df<-as.data.frame(LPS_wt_vs_PBS_wt)
write.csv(LPS_wt_vs_PBS_wt.df, file= "LPS_wt_vs_PBS_wt_05.csv")

#make violin plot for all differentially expressed genes in pbs vs lps in wildtype----
library(ggplot2)
library(EnhancedVolcano)
library(BuenColors)
Labels<-c('Fos','Mafb')
keyvals_pbs_lps <- ifelse(
  (LPS_wt_vs_PBS_wt.df$log2FoldChange< 0 & LPS_wt_vs_PBS_wt.df$padj < 0.05), 'royalblue',
  ifelse((LPS_wt_vs_PBS_wt.df$log2FoldChange > 0 & LPS_wt_vs_PBS_wt.df$padj< 0.05), 'chocolate1',
         'azure2'))
keyvals_pbs_lps[is.na(keyvals_pbs_lps)] <- 'azure'
names(keyvals_pbs_lps)[keyvals_pbs_lps == 'royalblue'] <- 'Decreasing'
names(keyvals_pbs_lps)[keyvals_pbs_lps == 'azure2'] <- 'None'
names(keyvals_pbs_lps)[keyvals_pbs_lps == 'chocolate1'] <- 'Increasing'

EnhancedVolcano(LPS_wt_vs_PBS_wt.df,
                lab= rownames(LPS_wt_vs_PBS_wt.df),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Custom colour over-ride',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0,
                pointSize = 2,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = FALSE,
                colAlpha = 3/5,
                legendPosition = 'right',
                drawConnectors = FALSE,
                colConnectors = 'black',
                colCustom = keyvals_pbs_lps)

#Gene ontology for all LPS driven genes----
library(BuenColors)
library(org.Mm.eg.db)
library(clusterProfiler)

LPS_PBS_comparecluster_data<-LPS_wt_vs_PBS_wt.df %>% dplyr::select(padj,log2FoldChange)
LPS_PBS_comparecluster_data <- LPS_PBS_comparecluster_data %>% filter(padj<0.05)
LPS_PBS_comparecluster_data <- LPS_PBS_comparecluster_data %>% mutate(Gene = rownames(LPS_PBS_comparecluster_data)) 
LPS_PBS_comparecluster_data <- LPS_PBS_comparecluster_data %>% mutate(ENTREZ = mapIds(org.Mm.eg.db,LPS_PBS_comparecluster_data$Gene, "ENTREZID", "SYMBOL")) 
LPS_PBS_comparecluster_data$Direction <- "upregulated"
LPS_PBS_comparecluster_data$Direction[LPS_PBS_comparecluster_data$log2FoldChange < 0] <- "downregulated"
lookup <- c(FC = "log2FoldChange", Entrez = "ENTREZ")
LPS_PBS_comparecluster_data<-rename(LPS_PBS_comparecluster_data, all_of(lookup))

LPS_PBS_comparecluster_data <- compareCluster(Entrez~Direction, data=LPS_PBS_comparecluster_data, 
                                                  fun=enrichGO, OrgDb = org.Mm.eg.db, ont="BP")

dotplot(LPS_PBS_comparecluster_data, x="Direction",showCategory = 10,size = "zScore", color = "p.adjust")

#Now genes in PBS vs LPS which are decreasing or increasing common with Mafb driven clusters in figure1----
#Trying different method to find which Mafb regulated genes are going down

LPS_wt_vs_PBS_wt.df_DECREASE<-LPS_wt_vs_PBS_wt.df %>% dplyr:: filter(log2FoldChange<0 & padj<0.05) #4887
LPS_wt_vs_PBS_wt.df_INCREASE<-LPS_wt_vs_PBS_wt.df %>% dplyr:: filter(log2FoldChange>0 & padj<0.05) #5037

LPS_wt_vs_PBS_wt.df_DECREASE_cluster1 <- intersect(Cluster1_genes,rownames(LPS_wt_vs_PBS_wt.df_DECREASE)) #out of 284 cluster 1 genes 77 are decreasing
LPS_wt_vs_PBS_wt.df_INCREASE_cluster1 <- intersect(Cluster1_genes,rownames(LPS_wt_vs_PBS_wt.df_INCREASE)) #out of 284 cluster 1 genes 155 are increasing 

LPS_wt_vs_PBS_wt.df_DECREASE_cluster2 <- intersect(Cluster2_genes,rownames(LPS_wt_vs_PBS_wt.df_DECREASE)) #out of 302 cluster 1 genes 92 are decreasing
LPS_wt_vs_PBS_wt.df_INCREASE_cluster2 <- intersect(Cluster2_genes,rownames(LPS_wt_vs_PBS_wt.df_INCREASE)) #out of 302 cluster 1 genes 98 are increasing 

#include this information into the file 
library(dplyr)
LPS_wt_vs_PBS_wt.df <- LPS_wt_vs_PBS_wt.df %>%   mutate(
  Mafb_DKO_BMMs_Cluster_Figure1 = dplyr::case_when(
    rownames(.) %in% Cluster1_genes ~ "Cluster_1",
    rownames(.) %in% Cluster2_genes ~ "Cluster_2",
    TRUE                             ~ "Not affected"
  )
)
  
write.csv(LPS_wt_vs_PBS_wt.df, file="LPS_wt_vs_PBS_wt.csv")
library(dplyr)
library(ggplot2)

fraction_table <- data.frame(
  Cluster = c("Cluster1","Cluster1","Cluster2","Cluster2"),
  Direction = c("Decrease","Increase","Decrease","Increase"),
  Overlap = c(length(LPS_wt_vs_PBS_wt.df_DECREASE_cluster1),
              length(LPS_wt_vs_PBS_wt.df_INCREASE_cluster1),
              length(LPS_wt_vs_PBS_wt.df_DECREASE_cluster2),
              length(LPS_wt_vs_PBS_wt.df_INCREASE_cluster2)),
  Total = c(length(Cluster1_genes),
            length(Cluster1_genes),
            length(Cluster2_genes),
            length(Cluster2_genes))
)

fraction_table$Fraction <- fraction_table$Overlap / fraction_table$Total

ggplot(fraction_table, aes(x = Cluster, y = Fraction, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(aes(label = paste0(Overlap,"/",Total)),
            position = position_dodge(width = 0.6),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Increase" = "#FC9272", "Decrease" = "#3182BD")) +
  ylab("Fraction of cluster genes") +
  theme_classic(base_size = 14)




#How are these genes different Gene ontology for the above cluster genes----
#Gene ontology for LPS vs PBS to show general known LPS driven genes are changing 


#Gene ontology to see major categories or function for which Mafb regulated genes go down
library(org.Mm.eg.db)
LPS_wt_vs_PBS_wt.df_DECREASE_cluster1_EntrezID <- mapIds(org.Mm.eg.db, LPS_wt_vs_PBS_wt.df_DECREASE_cluster1, "ENTREZID", "SYMBOL")
LPS_wt_vs_PBS_wt.df_INCREASE_cluster1_EntrezID <- mapIds(org.Mm.eg.db, LPS_wt_vs_PBS_wt.df_INCREASE_cluster1, "ENTREZID", "SYMBOL")
LPS_wt_vs_PBS_wt.df_DECREASE_cluster2_EntrezID <- mapIds(org.Mm.eg.db, LPS_wt_vs_PBS_wt.df_DECREASE_cluster2, "ENTREZID", "SYMBOL")
LPS_wt_vs_PBS_wt.df_INCREASE_cluster2_EntrezID <- mapIds(org.Mm.eg.db, LPS_wt_vs_PBS_wt.df_INCREASE_cluster2, "ENTREZID", "SYMBOL")

MakeGo_mouse<- function(geneid)
  
{
  enrichGO(gene = geneid,
           keyType = "ENTREZID",
           OrgDb         = org.Mm.eg.db,
           ont           = "BP",
           pAdjustMethod = "fdr",
           pvalueCutoff  = 0.01,
           qvalueCutoff  = 0.01,
           readable      = TRUE)
  
}


#EnrichGO for GO enrichment analysis
GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster1_EntrezID <- MakeGo_mouse(LPS_wt_vs_PBS_wt.df_DECREASE_cluster1_EntrezID)
GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster1_EntrezID <- MakeGo_mouse(LPS_wt_vs_PBS_wt.df_INCREASE_cluster1_EntrezID)
GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster2_EntrezID <- MakeGo_mouse(LPS_wt_vs_PBS_wt.df_DECREASE_cluster2_EntrezID)
GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster2_EntrezID <- MakeGo_mouse(LPS_wt_vs_PBS_wt.df_INCREASE_cluster2_EntrezID)

GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster1_EntrezID.df<-as.data.frame(GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster1_EntrezID)
GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster1_EntrezID.df<-as.data.frame(GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster1_EntrezID)
GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster2_EntrezID.df<-as.data.frame(GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster2_EntrezID)
GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster2_EntrezID.df<-as.data.frame(GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster2_EntrezID)


write.csv(GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster1_EntrezID.df, file="GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster1_EntrezID.csv")
write.csv(GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster1_EntrezID.df, file="GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster1_EntrezID.csv")
write.csv(GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster2_EntrezID.df, file="GO_LPS_wt_vs_PBS_wt.df_DECREASE_cluster2_EntrezID.csv")
write.csv(GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster2_EntrezID.df, file="GO_LPS_wt_vs_PBS_wt.df_INCREASE_cluster2_EntrezID.csv")

#Dotplot
#I think its better to use compare cluster 
library(BuenColors)
LPS_mafb_comparecluster_data<-LPS_wt_vs_PBS_wt.df %>% dplyr::select(log2FoldChange,Mafb_DKO_BMMs_Cluster_Figure1) %>% filter(Mafb_DKO_BMMs_Cluster_Figure1 %in% c("Cluster_1","Cluster_2"))
LPS_mafb_comparecluster_data <- LPS_mafb_comparecluster_data %>% mutate(Gene = rownames(LPS_mafb_comparecluster_data)) 
LPS_mafb_comparecluster_data <- LPS_mafb_comparecluster_data %>% mutate(ENTREZ = mapIds(org.Mm.eg.db,LPS_mafb_comparecluster_data$Gene, "ENTREZID", "SYMBOL")) 
LPS_mafb_comparecluster_data$Direction <- "upregulated"
LPS_mafb_comparecluster_data$Direction[LPS_mafb_comparecluster_data$log2FoldChange < 0] <- "downregulated"
lookup <- c(FC = "log2FoldChange", Entrez = "ENTREZ")
LPS_mafb_comparecluster_data<-rename(LPS_mafb_comparecluster_data, all_of(lookup))

LPS_Mafb_cluster_compare <- compareCluster(Entrez~Direction+Mafb_DKO_BMMs_Cluster_Figure1, data=LPS_mafb_comparecluster_data, 
                              fun=enrichGO, OrgDb = org.Mm.eg.db, ont="BP")

head(LPS_Mafb_cluster_compare)
dotplot(LPS_Mafb_cluster_compare, x="Direction",showCategory = 6,size = "zScore", color = "p.adjust") + facet_grid(~Mafb_DKO_BMMs_Cluster_Figure1)

#use individual tables as supplementary table



#Make a corresponding heatplot for decreasing genes overlapping with cluster 1, highlight C1qs----
rld_hao_df_wt_lps_cluster1_decrease <- rld_hao_df_wt_lps[rownames(rld_hao_df_wt_lps) %in% LPS_wt_vs_PBS_wt.df_DECREASE_cluster1, ]
rld_hao_df_wt_lps_cluster1_decrease<-t(scale(t(rld_hao_df_wt_lps_cluster1_decrease)))

#Annotation row
conditions <- rep(c("PBS", "LPS"), each = 4)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(rld_hao_df_wt_lps_cluster1_decrease)
col_ha_c1_decrease_lps_wt <- HeatmapAnnotation(df = annotation_col)


cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)


ht <- Heatmap(
  rld_hao_df_wt_lps_cluster1_decrease,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,             # hide all labels
  show_column_names = FALSE,
  col = col_fun,
  top_annotation = col_ha_c1_decrease_lps_wt
)

draw(ht)


#Now make combined heatmap for all Mafb regulated genes
rld_hao_df_wt_lps_cluster_mafb <- rld_hao_df_wt_lps[rownames(rld_hao_df_wt_lps) %in% LPS_mafb_comparecluster_data$Gene, ]
rld_hao_df_wt_lps_cluster_mafb<-t(scale(t(rld_hao_df_wt_lps_cluster_mafb)))
library(ComplexHeatmap)
library(colorRamp2)
#Annotation row
conditions <- rep(c("PBS", "LPS"), each = 4)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(rld_hao_df_wt_lps_cluster_mafb)
col_ha_cluster_lps_wt <- HeatmapAnnotation(df = annotation_col)

cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)

meta_rows <- LPS_mafb_comparecluster_data[
  match(
    rownames(rld_hao_df_wt_lps_cluster_mafb),
    LPS_mafb_comparecluster_data$Gene
  ),
]

row_split_factor <- paste(
  meta_rows$Mafb_DKO_BMMs_Cluster_Figure1,
  meta_rows$Direction,
  sep = "_"
)

# enforce order of panels
row_split_factor <- factor(
  row_split_factor,
  levels = c(
    "Cluster_1_downregulated",
    "Cluster_1_upregulated",
    "Cluster_2_downregulated",
    "Cluster_2_upregulated"
  )
)

combined_row_ha <- rowAnnotation(
  Cluster_Mafb_DKO = meta_rows$Mafb_DKO_BMMs_Cluster_Figure1,
  Direction        = meta_rows$Direction,
  show_annotation_name = TRUE
)


ht <- Heatmap(
  rld_hao_df_wt_lps_cluster_mafb,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = row_split_factor,  
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = col_fun,
  left_annotation = combined_row_ha,
  top_annotation = col_ha_cluster_lps_wt
)

draw(ht)


#===============================================
#Archived----
#Make GSEA for cluster 1 genes in PBS vs Wt condition for all cluster 1 genes
termtogene_cluster1_Fig1 <- data.frame(
  gs_name   = rep("Mafb_activated_transcripts", length(Cluster1_genes)),
  gs_symbol = Cluster1_genes
)

termtogene_cluster1_Fig1<-na.omit(termtogene_cluster1_Fig1)

#Make rnk dataset from Wt vs LPS peritoneum bulk rna 
# Pull out just the columns corresponding to gene symbols and LogFC
t_wt_vs_LPS_peritoneum_bulk <- as.data.frame(cbind(rownames(LPS_wt_vs_PBS_wt.df), LPS_wt_vs_PBS_wt.df$log2FoldChange)) 
t_wt_vs_LPS_peritoneum_bulk<- t_wt_vs_LPS_peritoneum_bulk[!duplicated(t_wt_vs_LPS_peritoneum_bulk$V1), ]
t_wt_vs_LPS_peritoneum_bulk <- t_wt_vs_LPS_peritoneum_bulk[!duplicated(t_wt_vs_LPS_peritoneum_bulk$V2), ]
t_wt_vs_LPS_peritoneum_bulk <- na.omit(t_wt_vs_LPS_peritoneum_bulk)

t_wt_vs_LPS_peritoneum_bulk_rnk <- as.numeric(t_wt_vs_LPS_peritoneum_bulk$V2)
names(t_wt_vs_LPS_peritoneum_bulk_rnk)<-as.character(t_wt_vs_LPS_peritoneum_bulk$V1)
t_wt_vs_LPS_peritoneum_bulk_rnk <- sort(t_wt_vs_LPS_peritoneum_bulk_rnk, decreasing = TRUE)


GSEA.results.Cluster1.wtvsLPS_bulk<- GSEA(t_wt_vs_LPS_peritoneum_bulk_rnk, TERM2GENE=termtogene_cluster1_Fig1, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.Cluster1.wtvsLPS_bulk<-as.tibble(GSEA.results.Cluster1.wtvsLPS_bulk)
gseaplot2(GSEA.results.Cluster1.wtvsLPS_bulk, geneSetID = 1,color= "#99000D")

leading_cluster1_WTVSLPS_bulk<-strsplit(GSEA.results.Cluster1.wtvsLPS_bulk@result[["core_enrichment"]],"/")
leading_cluster1_WTVSLPS_bulk<-leading_cluster1_WTVSLPS_bulk[[1]]

#Make GSEA for KANG DE genes enriched in Mafb overexpressed in DE
#Here expectation is that these signature must go down, as it happens in inflammation, if it is immunosuppresive
termtogene_KANG_leading_edge_Fig1 <- data.frame(
  gs_name   = rep("KANG_leading_edge_Fig1", length(leading_KANG_DE_10hr)),
  gs_symbol = leading_KANG_DE_10hr
)

termtogene_KANG_leading_edge_Fig1 <-na.omit(termtogene_KANG_leading_edge_Fig1)

GSEA.results.KANG_Depressed_onlyleadingedgeFig1.wtvsLPS_bulk<- GSEA(t_wt_vs_LPS_peritoneum_bulk_rnk, TERM2GENE=termtogene_KANG_leading_edge_Fig1, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.KANG_Depressed.wtvsLPS_bulk<-as.tibble(GSEA.results.KANG_Depressed_onlyleadingedgeFig1.wtvsLPS_bulk)
gseaplot2(GSEA.results.KANG_Depressed_onlyleadingedgeFig1.wtvsLPS_bulk, geneSetID = 1,color= "#99000D")
#its not significant
#which suggest that IFNgamma based immunosuppresive signature is not getting downregulated when treated with LPS
#Mafb regulated LPS based immunosuppresive signature is as I found below




#Heatmap for inspecting trends of of cluster 1 genes in Wt vs LPS condition

rld_hao_df_wt_lps <- as.data.frame(rld_hao_df_wt_lps)
rld_hao_df_wt_lps_leading_edge_cluster1 <- rld_hao_df_wt_lps[rownames(rld_hao_df_wt_lps) %in% leading_cluster1_WTVSLPS_bulk, ]

#Make a corresponding heatplot
rld_hao_df_wt_lps_leading_edge_cluster1<-t(scale(t(rld_hao_df_wt_lps_leading_edge_cluster1)))

#Annotation row
conditions <- rep(c("PBS", "LPS"), each = 4)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(rld_hao_df_wt_lps_leading_edge_cluster1)
col_ha_c1_lps_wt <- HeatmapAnnotation(df = annotation_col)


cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)


ht <- Heatmap(
  rld_hao_df_wt_lps_leading_edge_cluster1,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,             # hide all labels
  show_column_names = FALSE,
  col = col_fun,
  top_annotation = col_ha_c1_lps_wt
)

draw(ht)
#looks mostly increased.
#=======================================

#See dynamics of downregulated genes and molecular regulators using DecoupleR in single cell dataset----

library(Seurat)
library(BuenColors)
library(Nebulosa)

 #drawing some initial graphs----
Integrated_sc<-readRDS("All_Javi_seurat_CCA_processed_singleR_replicatecombined.RDS")
Unintegrated_sc<-readRDS("All_Javi_seurat_processed_singleR_replicatecombined.RDS")


a1<-DimPlot(Unintegrated_sc,
        group.by = "SingleR.labels", 
        cols = c("firebrick", "forestgreen", "darkolivegreen3", "deepskyblue3", "gray", "navy", "orchid", "darkorchid4","royalblue", "goldenrod", "azure4", "black", "darkgoldenrod4", "olivedrab1", "tan2","azure3", "forestgreen", "salmon"))

a2<-DimPlot(
  Unintegrated_sc,
  reduction = "umap",
  group.by = c("sampleid"),
  combine = FALSE,
  cols = c('#FCBBA1','#FB6A4A','#CB181D','#67000D'))

a1+a2

DimPlot(Integrated_sc,
            group.by = "SingleR.labels", 
            reduction="umap.cca",
            cols = c("firebrick", "forestgreen", "darkolivegreen3", "deepskyblue3", "gray", "navy", "orchid", "darkorchid4","royalblue", "goldenrod", "azure4", "black", "darkgoldenrod4", "olivedrab1", "tan2","azure3", "forestgreen", "salmon"))



#What are the proportion of resident to non-resident macrophages with time 

 #We subset the data to just Macrophage----
Macrophages_unitegrated <- subset(Unintegrated_sc, subset = SingleR.labels %in% c("Macrophages"))
Macrophages_integrated <- subset(Integrated_sc, subset = SingleR.labels %in% c("Macrophages"))

Macrophages_integrated@meta.data$SingleR.labels.fine <-Macrophages_unitegrated@meta.data$SingleR.labels.fine
Macrophages_integrated@meta.data$SingleR.labels.TRM <-
  ifelse(
    Macrophages_integrated@meta.data$SingleR.labels.fine %in%
      c("Macrophages (MF.II-480HI)", "Macrophages (MF.II+480LO)"),
    Macrophages_integrated@meta.data$SingleR.labels.fine,
    "Others"
  )

Macrophages_unitegrated@meta.data$SingleR.labels.TRM <-
  ifelse(
    Macrophages_unitegrated@meta.data$SingleR.labels.fine %in%
      c("Macrophages (MF.II-480HI)", "Macrophages (MF.II+480LO)"),
    Macrophages_unitegrated@meta.data$SingleR.labels.fine,
    "Others"
  )

p <- plot_density(
  Macrophages_unintegrated,
  c("Mafb"),
  pal = "magma"
)


#calculate fraction of TRM overtime

DimPlot(Macrophages_integrated, reduction = "umap.cca", group.by = "SingleR.labels.TRM", cols = jdb_palette("Moonrise3"), split.by = "samplid_combined")

library(dplyr)
library(dplyr)

prop_trm <- Macrophages_integrated@meta.data %>%
  dplyr::count(samplid_combined, SingleR.labels.TRM) %>%
  dplyr::group_by(samplid_combined) %>%
  dplyr::mutate(frac = n / sum(n))

library(ggplot2)

ggplot(prop_trm,
       aes(x = samplid_combined,
           y = frac,
           fill = SingleR.labels.TRM)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = jdb_palette("Moonrise3")) +
  ylab("Fraction of cells") +
  xlab("Time") +
  theme_classic()


library(BuenColors)
library(Seurat)
library(BuenColors)
library(ggplot2)


 #Heatmap for macrophaes regulated genes in bulk, this forms the basis of some dynamic changes in Mafb regulated genes----
DoHeatmap(
  Macrophages_unitegrated,
  features = LPS_wt_vs_PBS_wt.df_DECREASE_cluster1,
  group.by = "sampleid"
) +
  scale_fill_gradientn(
    colors = jdb_palette("solar_flare"),
    limits = c(-2, 2),
    oob = scales::squish
  )


genes_use <- intersect(
  LPS_wt_vs_PBS_wt.df_DECREASE_cluster2,
  rownames(Macrophages_unitegrated)
)

mat <- GetAssayData(
  Macrophages_unitegrated,
  assay = "RNA",
  layer  = "data"
)[genes_use, ]

library(Matrix)
meta <- Macrophages_unitegrated@meta.data

mat_avg <- sapply(
  unique(meta$sampleid),
  function(s) {
    cells <- rownames(meta)[meta$sampleid == s]
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  }
)

mat_avg <- as.matrix(mat_avg)

mat_z <- t(scale(t(mat_avg)))
mat_z[is.na(mat_z)] <- 0
sample_order <- unique(Macrophages_unitegrated@meta.data$sampleid)
mat_z <- mat_z[, sample_order, drop = FALSE]

library(ComplexHeatmap)
library(circlize)
library(BuenColors)
set.seed(123)
Heatmap(
  mat_z,
  name = "Z-score",
  cluster_rows = TRUE,        # 
  cluster_columns = FALSE,    # usually better for sample structure
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_km=4,
  col = jdb_palette("solar_flare"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9)
)

?heatmap

#Table for all type of cells

#Done SIMILAR for increasing genes 
#Very interesting pattern appear, the decreasing genes are more sensitive to LPS treatment than increasing genes?
 #Visualize average score in feature plot and vlnplots----
genes_use <- intersect(
  LPS_wt_vs_PBS_wt.df_INCREASE_cluster2,
  rownames(Macrophages_integrated)
)
length(genes_use)

expr <- GetAssayData(
  Macrophages_integrated,
  assay = "RNA",
  layer  = "data"
)

mean.exp_Cluster2_increasing_LPSvsPBS <- Matrix::colMeans(
  expr[genes_use, , drop = FALSE],
  na.rm = TRUE
)

Macrophages_integrated@meta.data$Cluster2_increasing_Mafb_DKO <-
  mean.exp_Cluster2_increasing_LPSvsPBS
library(dplyr)
Macrophages_integrated@meta.data <- Macrophages_integrated@meta.data %>%
  mutate(
    samplid_combined = case_when(
      grepl("^LPS_24H", sampleid) ~ "LPS_24H",
      grepl("^LPS_4H",  sampleid) ~ "LPS_4H",
      grepl("^LPS_1H",  sampleid) ~ "LPS_1H",
      grepl("^NT_0H",   sampleid) ~ "NT_0H",
      TRUE ~ NA_character_
    )
  )

Macrophages_integrated@meta.data$samplid_combined <- factor(
  Macrophages_integrated@meta.data$samplid_combined,
  levels = c("NT_0H", "LPS_1H", "LPS_4H", "LPS_24H")
)

FeaturePlot(
  object = Macrophages_integrated,
  features = "Cluster2_decreasing_Mafb_DKO",
  reduction = "umap.cca",
  split.by = "samplid_combined",
  min.cutoff = "q10",
  max.cutoff = "q90"
) &
  scale_colour_gradientn(
    colours = jdb_palette("solar_extra"),
  )

median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}

p4<-VlnPlot(object = Macrophages_integrated,features="Cluster2_increasing_Mafb_DKO",group.by = "samplid_combined",pt.size = 0) +
  stat_summary(fun.y = median.stat, geom='point', size = 10, colour = "black",shape=95) 
ylims <- range(
  Macrophages_integrated@meta.data$Cluster1_increasing_Mafb_DKO,
  Macrophages_integrated@meta.data$Cluster2_increasing_Mafb_DKO,
  na.rm = TRUE
)

p3 <- p3 + coord_cartesian(ylim = ylims)
p4 <- p4 + coord_cartesian(ylim = ylims)

library(patchwork)

p3 | p4

#Make a ggplot showing average gene expression per condition for 

#DecoupleR based networks----
library(decoupleR)
library(SeuratObject)
library(Seurat)
library(SeuratWrappers)
net <- decoupleR::get_collectri(organism = 'mouse', 
                                split_complexes = FALSE)

saveRDS(Macrophages_unitegrated, file="Macrophages_unitegrated.RDS")
saveRDS(Macrophages_integrated, file="Macrophages_integrated.RDS")

bulk <- AggregateExpression(Macrophages_unitegrated, group.by = "sampleid", return.seurat = TRUE)
Cells(bulk)

#take out the matrix of normalized count
mat_bulk_macrophages <- GetAssayData(
  bulk,
  assay = "RNA",
  layer  = "data"
)

mat_bulk_macrophages <- as.matrix(mat_bulk_macrophages)

mat_sc_macrophages <- GetAssayData(
  Macrophages_unitegrated,
  assay = "RNA",
  layer  = "data"
)


 #SINGLE CELL----
genes_use <- Maf_significant_genes 
mat_sc_macrophages_sub <- mat_sc_macrophages[
  rownames(mat_sc_macrophages) %in% genes_use,
  ,
  drop = FALSE
]

mat_sc_macrophages_sub <- as.matrix(mat_sc_macrophages_sub)


#run ulm
net_sub<- net %>%
  dplyr::filter(target %in% Maf_significant_genes)

# Run ulm just for Mafb regulated genes 
acts <- decoupleR::run_ulm(mat = mat_sc_macrophages_sub, 
                           net = net_sub, 
                           .source = 'source', 
                           .target = 'target',
                           .mor='mor', 
                           minsize = 10)

Macrophages_unitegrated[['tfsulm_Mafb_significant_genes_DKO_BMM']] <- acts %>%
  tidyr::pivot_wider(id_cols = 'source', 
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = Macrophages_unitegrated) <- "tfsulm_Mafb_significant_genes_DKO_BMM"

# Scale the data
Macrophages_unitegrated <- Seurat::ScaleData(Macrophages_unitegrated)
Macrophages_unitegrated@assays$tfsulm@data <- Macrophages_unitegrated@assays$tfsulm@scale.data

rownames(Macrophages_unitegrated)

Seurat::FeaturePlot(Macrophages_unitegrated, features = c("Jun"), pt.size = 1) & scale_colour_gradient2(
  low = jdb_palette("solar_flare")[1],
  mid = "white",
  high = jdb_palette("solar_flare")[length(jdb_palette("solar_flare"))],
  midpoint = 0
)

VlnPlot(object = Macrophages_unitegrated,features="Spi1",group.by = "samplid_combined",pt.size = 0) +
  stat_summary(fun.y = median.stat, geom='point', size = 10, colour = "black",shape=95) 

#top 25 tfs

n_tfs <- 25
Macrophages_unitegrated@meta.data <- Macrophages_unitegrated@meta.data %>%
  mutate(
    samplid_combined = case_when(
      grepl("^LPS_24H", sampleid) ~ "LPS_24H",
      grepl("^LPS_4H",  sampleid) ~ "LPS_4H",
      grepl("^LPS_1H",  sampleid) ~ "LPS_1H",
      grepl("^NT_0H",   sampleid) ~ "NT_0H",
      TRUE ~ NA_character_
    )
  )
unique(Macrophages_unitegrated$samplid_combined)
Macrophages_unitegrated@meta.data$samplid_combined <- factor(
  Macrophages_unitegrated@meta.data$samplid_combined,
  levels = c("NT_0H", "LPS_1H", "LPS_4H", "LPS_24H")
)

# Extract activities from object as a long dataframe
df <- t(as.matrix(Macrophages_unitegrated@assays$tfsulm_Mafb_significant_genes_DKO_BMM@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(Condition = Macrophages_unitegrated@meta.data$samplid_combined) %>%
  tidyr::pivot_longer(cols = -Condition, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(Condition, source) %>%
  dplyr::summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(mean)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  dplyr::filter(source %in% tfs) %>%
  tidyr::pivot_wider(id_cols = 'Condition', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames('Condition') %>%
  as.matrix()

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))


# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20,
                   cluster_rows = FALSE) 

#Take out top 20 from each condition and see common 
library(dplyr)

top20_per_condition <- df %>%
  group_by(Condition) %>%
  arrange(desc(mean)) %>%
  slice_head(n = 20) %>%
  ungroup()

tf_list <- top20_per_condition %>%
  group_by(Condition) %>%
  summarise(tfs = list(unique(source)), .groups = "drop") %>%
  tibble::deframe()


common_all <- Reduce(intersect, tf_list)
common_all

combn(names(tf_list), 2, function(x) {
  intersect(tf_list[[x[1]]], tf_list[[x[2]]])
}, simplify = FALSE)

library(ggVennDiagram)

ggVennDiagram(
  tf_list
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_void()


#Also do for negative 

# Top 20 lowest mean per condition
top20_neg <- df %>%
  group_by(Condition) %>%
  arrange(mean, .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

# Make TF list per condition
tf_list_neg <- top20_neg %>%
  group_by(Condition) %>%
  summarise(tfs = list(unique(source)), .groups = "drop") %>%
  tibble::deframe()

# Common TFs across all conditions
common_neg_all <- Reduce(intersect, tf_list_neg)

common_neg_all

combn(names(tf_list_neg), 2, function(x) {
  intersect(tf_list_neg[[x[1]]], tf_list[[x[2]]])
}, simplify = FALSE)

library(ggVennDiagram)

ggVennDiagram(
  tf_list_neg
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_void()


#Network plot ----
library(dplyr)
"Fosl2","Fos","Jun","Junb","Cebp"
tfs_to_plot <- c("Fos","Fosl1","Fosl2")  # change as needed
genes_to_plot<-Cluster1_genes
net_sub <- net %>%
  mutate(
    source_key = source,
    target_key = target
  ) %>%
  filter(source_key %in% tfs_to_plot) %>%
  filter(target_key %in% genes_to_plot)

library(dplyr)
library(tidyr)

tf_activity_mat <- df %>%
  tidyr::pivot_wider(
    id_cols = source,
    names_from = Condition,
    values_from = mean
  ) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()


library(igraph)
library(tidygraph)
library(ggraph)
library(scales)
library(ggnewscale)

make_tf_graph_sc <- function(net, logfc_mat, tf_mat, timepoint) {
  
  logfc_df <- tibble(
    gene_key = rownames(logfc_mat),
    logFC = logfc_mat[, timepoint]
  )
  
  tf_df <- tibble(
    source_key = rownames(tf_mat),
    tf_score = tf_mat[, timepoint]
  )
  
  df <- net %>%
    left_join(logfc_df, by = c("target_key" = "gene_key")) %>%
    filter(!is.na(logFC))
  
  edges <- df %>%
    transmute(
      from = source,
      to = target
    )
  
  nodes <- bind_rows(
    
    df %>%
      distinct(source, source_key) %>%
      transmute(
        name = source,
        name_key = source_key,
        is_tf = TRUE
      ),
    
    df %>%
      distinct(target, target_key) %>%
      transmute(
        name = target,
        name_key = target_key,
        is_tf = FALSE
      )
    
  ) %>%
    distinct(name, .keep_all = TRUE) %>%
    left_join(logfc_df, by = c("name_key" = "gene_key")) %>%
    left_join(tf_df, by = c("name_key" = "source_key")) %>%
    mutate(
      size = ifelse(is_tf, 10, 5)
    )
  
  graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
}

plot_tf_network_sc <- function(g, layout, title) {
  
  g_tbl <- as_tbl_graph(g)
  
  ggraph(g_tbl, layout = layout) +
    
    geom_edge_link(alpha = 0.4,
                   arrow = arrow(length = unit(2.5, "mm"))) +
    
    # Gene nodes
    geom_node_point(
      data = function(x) dplyr::filter(x, !is_tf),
      aes(fill = logFC),
      shape = 21,
      size = 6,
      color = "black"
    ) +
    
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-3, 3),
      name = "log2FC"
    ) +
    
    ggnewscale::new_scale_fill() +
    
    # TF nodes
    geom_node_point(
      data = function(x) dplyr::filter(x, is_tf),
      aes(fill = tf_score),
      shape = 22,
      size = 10,
      color = "black"
    ) +
    
    scale_fill_viridis_c(
      option = "plasma",
      limits = c(-3, 3),
      name = "TF activity"
    ) +
    
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    
    theme_void() +
    ggtitle(title)
}

timepoints <- c("LPS_1H", "LPS_4H", "LPS_24H")

plots <- lapply(timepoints, function(tp) {
  
  g <- make_tf_graph_sc(
    net = net_sub,
    logfc_mat = logfc_mat,
    tf_mat = tf_activity_mat,
    timepoint = tp
  )
  
  set.seed(123)
  layout_fixed <- layout_with_fr(g)
  
  plot_tf_network_sc(g, layout_fixed, tp)
})

library(patchwork)
wrap_plots(plots)


library(dplyr)

# Use normalized data
mat <- GetAssayData(
  Macrophages_unitegrated,
  assay = "RNA",
  layer = "data"
)

mat <- as.matrix(mat)

meta <- Macrophages_unitegrated@meta.data
head(meta)
# Ensure correct order
meta$samplid_combined <- factor(
  meta$samplid_combined,
  levels = c("NT_0H", "LPS_1H", "LPS_4H", "LPS_24H")
)
unique(meta$samplid_combined)
# Split cells by condition
cells_by_condition <- split(colnames(mat), meta$samplid_combined)

# Compute mean per gene per condition
mean_expr <- sapply(cells_by_condition, function(cells) {
  rowMeans(mat[, cells, drop = FALSE])
})

mean_expr <- as.matrix(mean_expr)

pseudo <- 1e-6

logfc_mat <- sweep(
  log2(mean_expr + pseudo),
  1,
  log2(mean_expr[, "NT_0H"] + pseudo),
  FUN = "-"
)

# Remove NT column (it will be zero)
logfc_mat <- logfc_mat[, colnames(logfc_mat) != "NT_0H"]
head(logfc_mat)


logfc_mat <- logfc_mat[
  rownames(logfc_mat) %in% genes_to_plot,
  ,
  drop = FALSE
]






