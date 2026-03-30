#Figure 3

#Here we introduce the concept of Homo het 
#Then in the mice , we show that those genes are already affected , specifically the downregulated ones

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
rld_hao_df_pbs_wt_het<-rld_hao_df %>% dplyr::select(PBS_WT_1,PBS_WT_2,PBS_WT_3,PBS_WT_4,PBS_E269R_1,PBS_E269R_2,PBS_E269R_3,PBS_E269R_4)

mod_mat <- model.matrix(design(dds_hao), colData(dds_hao))        
mod_mat #see how it looks like
LPS_E269R <- colMeans(mod_mat[dds_hao$Condition == "LPS_E269R", ])
LPS_WT <- colMeans(mod_mat[dds_hao$Condition == "LPS_WT", ])
PBS_E269R <- colMeans(mod_mat[dds_hao$Condition == "PBS_E269R", ])
PBS_WT <- colMeans(mod_mat[dds_hao$Condition == "PBS_WT", ])
?model.matrix
#Run differential analysis comparing different condition
library(DESeq2)
PBS_e269r_vs_PBS_wt <- results(dds_hao, contrast = PBS_E269R - PBS_WT, alpha = 0.05)
LPS_e269r_vs_PBS_wt <- results(dds_hao, contrast = LPS_E269R - PBS_WT, alpha = 0.05)
LPS_e269r_vs_LPS_wt <- results(dds_hao, contrast = LPS_E269R - LPS_WT, alpha = 0.05)
LPS_wt_vs_PBS_wt <- results(dds_hao, contrast = LPS_WT - PBS_WT, alpha = 0.05)

summary(PBS_e269r_vs_PBS_wt)
PBS_e269r_vs_PBS_wt.df<-as.data.frame(PBS_e269r_vs_PBS_wt)
write.csv(PBS_e269r_vs_PBS_wt.df, file= "PBS_e269r_vs_PBS_wt.csv")

LPS_e269r_vs_PBS_wt.df<-as.data.frame(LPS_e269r_vs_PBS_wt)
write.csv(LPS_e269r_vs_PBS_wt.df, file= "LPS_e269r_vs_PBS_wt.csv")

LPS_e269r_vs_LPS_wt.df<-as.data.frame(LPS_e269r_vs_LPS_wt)
write.csv(LPS_e269r_vs_LPS_wt.df, file= "LPS_e269r_vs_LPS_wt.csv")

LPS_wt_vs_PBS_wt.df<-as.data.frame(LPS_wt_vs_PBS_wt)
write.csv(LPS_wt_vs_PBS_wt.df, file= "LPS_wt_vs_PBS_wt.csv")

#Make volcano plot, PCA and common plots as filler
library(ggplot2)
library(dplyr)
library(patchwork)

plot_volcano <- function(df, title, ymax_cut = 10){
  
  df <- df %>%
    mutate(Regulation = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    neglog10 = -log10(padj),
    neglog10_cap = ifelse(neglog10 > ymax_cut, ymax_cut, neglog10))
  
  # Count genes
  up_n <- sum(df$Regulation == "Up", na.rm = TRUE)
  down_n <- sum(df$Regulation == "Down", na.rm = TRUE)
  
  ggplot(df, aes(log2FoldChange, neglog10_cap, color = Regulation)) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_manual(values = c(
      Up = "red",
      Down = "blue",
      NS = "grey80"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    
    annotate("text",
             x = max(df$log2FoldChange, na.rm = TRUE),
             y = ymax_cut * 0.5,
             label = paste0("Up: ", up_n),
             hjust = 1, color = "red", size = 4) +
    
    annotate("text",
             x = min(df$log2FoldChange, na.rm = TRUE),
             y = ymax_cut * 0.5,
             label = paste0("Down: ", down_n),
             hjust = 0, color = "blue", size = 4) +
    
    coord_cartesian(ylim = c(0, ymax_cut)) +
    
    theme_classic() +
    labs(title = title,
         x = "log2FC",
         y = "-log10(adj P)")
}

p1 <- plot_volcano(PBS_e269r_vs_PBS_wt.df, "PBS_e269r_vs_PBS_wt")
p2 <- plot_volcano(LPS_e269r_vs_LPS_wt.df, "LPS_e269r_vs_LPS_wt")
p3 <- plot_volcano(LPS_wt_vs_PBS_wt.df, "LPS_wt_vs_PBS_wt.df")

final_plot <- (p1 | p2 | p3 )

#PCA with selected set of genes----
LPS_wt_vs_PBS_wt.df.significant<-filter(LPS_wt_vs_PBS_wt.df,padj<0.05)
LPS_wt_vs_PBS_wt.df_cluster1<-intersect(rownames(LPS_wt_vs_PBS_wt.df.significant),Cluster1_genes)
Cluster1_genes_PCA<-intersect(LPS_wt_vs_PBS_wt.df_cluster1, rownames(rld_hao_df))
rld_cluster1 <- rld_hao_df[Cluster1_genes_PCA, ]
pca_cluster1 <- prcomp(t(rld_cluster1), scale. = TRUE)
pca_df <- as.data.frame(pca_cluster1$x)
pca_df$Condition <- Meta_combined$Condition

library(ggplot2)
percentVar <- pca_cluster1$sdev^2 / sum(pca_cluster1$sdev^2)
percentVar <- round(percentVar * 100, 1)
ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    title = "PCA using Cluster2 genes ",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  )

#plotPCA(rld_hao[Cluster1_genes_PCA, ], intgroup = "Condition")

#Now genes in PBS vs LPS which are decreasing or increasing common with Mafb driven clusters in figure1----
#Trying different method to find which Mafb regulated genes are going down

#PBS_e269r_vs_PBS_wt
PBS_e269r_vs_PBS_wt.df_DECREASE<-PBS_e269r_vs_PBS_wt.df %>% dplyr:: filter(log2FoldChange<0 & padj<0.05) #209
PBS_e269r_vs_PBS_wt.df_INCREASE<-PBS_e269r_vs_PBS_wt.df %>% dplyr:: filter(log2FoldChange>0 & padj<0.05) #141

PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1 <- intersect(Cluster1_genes,rownames(PBS_e269r_vs_PBS_wt.df_DECREASE)) #out of 284 cluster 1 genes 24 are decreasing,
PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1_wtvsLPS <- intersect(LPS_wt_vs_PBS_wt.df_DECREASE_cluster1,PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1) #6 are already decreasing which were decreasing in LPS

PBS_e269r_vs_PBS_wt.df_INCREASE_cluster1 <- intersect(Cluster1_genes,rownames(PBS_e269r_vs_PBS_wt.df_INCREASE)) #out of 284 cluster 1 genes 0 are decreasing  
PBS_e269r_vs_PBS_wt.df_INCREASE_cluster1_wtvsLPS <- intersect(LPS_wt_vs_PBS_wt.df_INCREASE_cluster1,PBS_e269r_vs_PBS_wt.df_INCREASE_cluster1) #6 are already decreasing which were decreasing in LPS

PBS_e269r_vs_PBS_wt.df_DECREASE_cluster2 <- intersect(Cluster2_genes,rownames(PBS_e269r_vs_PBS_wt.df_DECREASE)) #out of 302 cluster 1 genes 4 are decreasing
PBS_e269r_vs_PBS_wt.df_DECREASE_cluster2_wtvsLPS <- intersect(LPS_wt_vs_PBS_wt.df_DECREASE_cluster2,PBS_e269r_vs_PBS_wt.df_DECREASE_cluster2) #1 are already decreasing which were decreasing in LPS

PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2 <- intersect(Cluster2_genes,rownames(PBS_e269r_vs_PBS_wt.df_INCREASE)) #out of 302 cluster 1 genes 15 are increasing 
PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2_wtvsLPS <- intersect(LPS_wt_vs_PBS_wt.df_INCREASE_cluster1,PBS_e269r_vs_PBS_wt.df_INCREASE_cluster1) #0 are already decreasing which were decreasing in LPS

#LPS_e269r_vs_PBS_wt
LPS_e269r_vs_PBS_wt.df_DECREASE<-LPS_e269r_vs_PBS_wt.df %>% dplyr:: filter(log2FoldChange<0 & padj<0.05) #5134
LPS_e269r_vs_PBS_wt.df_INCREASE<-LPS_e269r_vs_PBS_wt.df %>% dplyr:: filter(log2FoldChange>0 & padj<0.05) #5326

LPS_e269r_vs_PBS_wt.df_DECREASE_cluster1 <- intersect(Cluster1_genes,rownames(LPS_e269r_vs_PBS_wt.df_DECREASE)) #out of 284 cluster 1 genes 73 are decreasing,
LPS_e269r_vs_PBS_wt.df_DECREASE_cluster1_wtvsLPS <- intersect(LPS_wt_vs_PBS_wt.df_DECREASE_cluster1,LPS_e269r_vs_PBS_wt.df_DECREASE_cluster1) #67 are same as in LPS condition (compare fold change)

LPS_e269r_vs_PBS_wt.df_INCREASE_cluster1 <- intersect(Cluster1_genes,rownames(LPS_e269r_vs_PBS_wt.df_INCREASE)) #out of 284 cluster 1 genes 123 are increasing  
LPS_e269r_vs_PBS_wt.df_INCREASE_cluster1_wtvsLPS <- intersect(LPS_wt_vs_PBS_wt.df_INCREASE_cluster1,LPS_e269r_vs_PBS_wt.df_INCREASE_cluster1) #107 are same as LPS

LPS_e269r_vs_PBS_wt.df_DECREASE_cluster2 <- intersect(Cluster2_genes,rownames(LPS_e269r_vs_PBS_wt.df_DECREASE)) #out of 302 cluster 2 genes 85 are decreasing
LPS_e269r_vs_PBS_wt.df_DECREASE_cluster2_wtvsLPS <- intersect(LPS_wt_vs_PBS_wt.df_DECREASE_cluster2,LPS_e269r_vs_PBS_wt.df_DECREASE_cluster2) #73 are decreasing which were decreasing in LPS

LPS_e269r_vs_PBS_wt.df_INCREASE_cluster2 <- intersect(Cluster2_genes,rownames(LPS_e269r_vs_PBS_wt.df_INCREASE)) #out of 302 cluster 2 genes 112 are increasing
LPS_e269r_vs_PBS_wt.df_INCREASE_cluster2_wtvsLPS <- intersect(LPS_wt_vs_PBS_wt.df_INCREASE_cluster2,LPS_e269r_vs_PBS_wt.df_INCREASE_cluster2) #96 are increasing which were increasing in LPS


#LPS_e269r_vs_LPS_wt
LPS_e269r_vs_LPS_wt.df_DECREASE<-LPS_e269r_vs_LPS_wt.df %>% dplyr:: filter(log2FoldChange<0 & padj<0.05) #5134
LPS_e269r_vs_LPS_wt.df_INCREASE<-LPS_e269r_vs_LPS_wt.df %>% dplyr:: filter(log2FoldChange>0 & padj<0.05) #5326

LPS_e269r_vs_LPS_wt.df_DECREASE_cluster1 <- intersect(Cluster1_genes,rownames(LPS_e269r_vs_LPS_wt.df_DECREASE)) #out of 284 cluster 1 genes 12 are DIFFERENTLY decreasing,

LPS_e269r_vs_LPS_wt.df_INCREASE_cluster1 <- intersect(Cluster1_genes,rownames(LPS_e269r_vs_LPS_wt.df_INCREASE)) #out of 284 cluster 1 genes 14 are differently increasing  

LPS_e269r_vs_LPS_wt.df_DECREASE_cluster2 <- intersect(Cluster2_genes,rownames(LPS_e269r_vs_LPS_wt.df_DECREASE)) #out of 302 cluster 2 genes 1 are differently decreasing decreasing

LPS_e269r_vs_LPS_wt.df_INCREASE_cluster2 <- intersect(Cluster2_genes,rownames(LPS_e269r_vs_LPS_wt.df_INCREASE)) #out of 302 cluster 2 genes 28 are increasing



# Pull out just the columns corresponding to gene symbols and LogFC
t_PBS_vs_PBS_het <- as.data.frame(cbind(rownames(PBS_e269r_vs_PBS_wt.df), PBS_e269r_vs_PBS_wt.df$log2FoldChange)) 
t_PBS_vs_PBS_het <- t_PBS_vs_PBS_het[!duplicated(t_PBS_vs_PBS_het$V1), ]
t_PBS_vs_PBS_het <- t_PBS_vs_PBS_het[!duplicated(t_PBS_vs_PBS_het$V2), ]
t_PBS_vs_PBS_het <- na.omit(t_PBS_vs_PBS_het)

t_PBS_vs_PBS_het_rnk <- as.numeric(t_PBS_vs_PBS_het$V2)
names(t_PBS_vs_PBS_het_rnk)<-as.character(t_PBS_vs_PBS_het$V1)
t_PBS_vs_PBS_het_rnk <- sort(t_PBS_vs_PBS_het_rnk, decreasing = TRUE)

t_LPS_vs_LPS_het <- as.data.frame(cbind(rownames(LPS_e269r_vs_LPS_wt.df), LPS_e269r_vs_LPS_wt.df$log2FoldChange)) 
t_LPS_vs_LPS_het <- t_LPS_vs_LPS_het[!duplicated(t_LPS_vs_LPS_het$V1), ]
t_LPS_vs_LPS_het <- t_LPS_vs_LPS_het[!duplicated(t_LPS_vs_LPS_het$V2), ]
t_LPS_vs_LPS_het<- na.omit(t_LPS_vs_LPS_het)

t_LPS_vs_LPS_het_rnk <- as.numeric(t_LPS_vs_LPS_het$V2)
names(t_LPS_vs_LPS_het_rnk)<-as.character(t_LPS_vs_LPS_het$V1)
t_LPS_vs_LPS_het_rnk <- sort(t_LPS_vs_LPS_het_rnk, decreasing = TRUE)



#Show by GSEA for those genesets----
gene_sets <- list(
  LPS_wt_vs_PBS_wt_DECREASE_cluster1 = na.omit(LPS_wt_vs_PBS_wt.df_DECREASE_cluster1),
  LPS_wt_vs_PBS_wt_INCREASE_cluster1 = na.omit(LPS_wt_vs_PBS_wt.df_INCREASE_cluster1),
  LPS_wt_vs_PBS_wt_DECREASE_cluster2 = na.omit(LPS_wt_vs_PBS_wt.df_DECREASE_cluster2),
  LPS_wt_vs_PBS_wt_INCREASE_cluster2 = na.omit(LPS_wt_vs_PBS_wt.df_INCREASE_cluster2)
)

gene_sets <- list(
  cluster1 = na.omit(Cluster1_genes),
  cluster2 = na.omit(Cluster2_genes)
)



# convert to TERM2GENE
termtogene <- stack(gene_sets)
colnames(termtogene) <- c("gs_symbol","gs_name")
termtogene <- termtogene[,c("gs_name","gs_symbol")]

# 2. Run GSEA (6hr and 10hr)
gsea_PBSvsPBS_het <- GSEA(t_PBS_vs_PBS_het_rnk,
                 TERM2GENE = termtogene,
                 verbose = FALSE,
                 eps = 0,
                 pvalueCutoff = 1,
                 by="fgsea",
                 minGSSize = 5)

gsea_LPSvsLPShet <- GSEA(t_LPS_vs_LPS_het_rnk,
                  TERM2GENE = termtogene,
                  verbose = FALSE,
                  eps = 0,
                  pvalueCutoff = 1,
                  by="fgsea",
                  minGSSize = 5)
# 3. Plot all GSEA curves

p1 <- gseaplot2(gsea_PBSvsPBS_het,
               geneSetID = c("LPS_wt_vs_PBS_wt_DECREASE_cluster1","LPS_wt_vs_PBS_wt_INCREASE_cluster1","LPS_wt_vs_PBS_wt_DECREASE_cluster2","LPS_wt_vs_PBS_wt_INCREASE_cluster2"),
               color=c("#756bb1","#efedf5","#fec44f","#fff7bc"),
                title="PBS vs PBS-het")
                
p2 <- gseaplot2(gsea_LPSvsLPShet,
                geneSetID = c("LPS_wt_vs_PBS_wt_DECREASE_cluster1","LPS_wt_vs_PBS_wt_INCREASE_cluster1","LPS_wt_vs_PBS_wt_DECREASE_cluster2","LPS_wt_vs_PBS_wt_INCREASE_cluster2"),
                color=c("#756bb1","#efedf5","#fec44f","#fff7bc"),
                title="LPS vs LPS-het")

p3 <- gseaplot2(gsea_PBSvsPBS_het,
                geneSetID = c("cluster1","cluster2"),
                color=c("#756bb1","#fec44f"),
                title="PBS vs PBS-het")

p4 <- gseaplot2(gsea_LPSvsLPShet,
                geneSetID = c("cluster1","cluster2"),
                color=c("#756bb1","#fec44f"),
                title="LPS vs LPS-het")


#Making Venn diagram, to see effects which are LPS specific or Mafb-Fos specific---- 

library(ggVennDiagram)
#gene common in het mouse in PBS with Maf DKO overexpression
PBS_e269r_vs_PBS_wt.df_DECREASE_genes<-rownames(PBS_e269r_vs_PBS_wt.df_DECREASE)
PBS_e269r_vs_PBS_wt.df_INCREASE_genes<-rownames(PBS_e269r_vs_PBS_wt.df_INCREASE)
genes <- list(
 Set1 = PBS_e269r_vs_PBS_wt.df_DECREASE_genes,
 Set2 = Cluster1_genes,
 Set3 = PBS_e269r_vs_PBS_wt.df_INCREASE_genes,
 Set4 = Cluster2_genes
)

ggVennDiagram(genes, label_alpha = 0) +
  scale_fill_gradient(low="white", high="darkred")

#gene common in het mouse in LBS with Maf DKO overexpression
LPS_e269r_vs_PBS_wt.df_DECREASE_genes<-rownames(LPS_e269r_vs_PBS_wt.df_DECREASE)
LPS_e269r_vs_PBS_wt.df_INCREASE_genes<-rownames(LPS_e269r_vs_PBS_wt.df_INCREASE)

LPS_wt_vs_PBS_wt.df_DECREASE_genes<-rownames(LPS_wt_vs_PBS_wt.df_DECREASE)
LPS_wt_vs_PBS_wt.df_INCREASE_genes<-rownames(LPS_wt_vs_PBS_wt.df_INCREASE)
genes <- list(
  Set1 = LPS_e269r_vs_PBS_wt.df_DECREASE_genes,
  Set2 = Cluster1_genes,
  Set3 = LPS_e269r_vs_PBS_wt.df_INCREASE_genes,
  Set4 = Cluster2_genes
)

ggVennDiagram(genes, label_alpha = 0) +
  scale_fill_gradient(low="white", high="darkred")

#gene common in het mouse and wt mouse in LBS
genes <- list(
  Set1 = LPS_wt_vs_PBS_wt.df_DECREASE_genes,
  Set2 =LPS_e269r_vs_PBS_wt.df_DECREASE_genes,
  Set3 = LPS_wt_vs_PBS_wt.df_INCREASE_genes,
  Set4 = LPS_e269r_vs_PBS_wt.df_INCREASE_genes
)

ggVennDiagram(genes, label_alpha = 0) +
  scale_fill_gradient(low="white", high="darkred")

#common between lps het wt and Cluster 12
LPS_e269r_vs_LPS_wt.df_DECREASE_genes<-rownames(LPS_e269r_vs_LPS_wt.df_DECREASE)
LPS_e269r_vs_LPS_wt.df_INCREASE_genes<-rownames(LPS_e269r_vs_LPS_wt.df_INCREASE)
genes <- list(
  Set1 = LPS_e269r_vs_LPS_wt.df_DECREASE_genes,
  Set2 = Cluster1_genes,
  Set3 = LPS_e269r_vs_LPS_wt.df_INCREASE_genes,
  Set4 = Cluster2_genes
)

ggVennDiagram(genes, label_alpha = 0) +
  scale_fill_gradient(low="white", high="darkred")


#Fractions

library(dplyr)
library(ggplot2)

fraction_table <- data.frame(
  Cluster = c("Cluster1","Cluster1","Cluster2","Cluster2"),
  Direction = c("Decrease","Increase","Decrease","Increase"),
  Overlap = c(length(LPS_e269r_vs_LPS_wt.df_DECREASE_cluster1),
              length(LPS_e269r_vs_LPS_wt.df_INCREASE_cluster1),
              length(LPS_e269r_vs_LPS_wt.df_DECREASE_cluster2),
              length(LPS_e269r_vs_LPS_wt.df_INCREASE_cluster2)),
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

#Fraction of LPS
fraction_table <- data.frame(
  Cluster = c("Cluster1","Cluster1","Cluster2","Cluster2"),
  Direction = c("Decrease","Increase","Decrease","Increase"),
  Overlap = c(length(PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1),
              length(PBS_e269r_vs_PBS_wt.df_INCREASE_cluster1),
              length(PBS_e269r_vs_PBS_wt.df_DECREASE_cluster2),
              length(PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2)),
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



#include this information into the file----
library(dplyr)
PBS_e269r_vs_PBS_wt.df <- PBS_e269r_vs_PBS_wt.df %>%   mutate(
  Mafb_DKO_BMMs_Cluster_Figure1 = dplyr::case_when(
    rownames(.) %in% Cluster1_genes ~ "Cluster_1",
    rownames(.) %in% Cluster2_genes ~ "Cluster_2",
    TRUE                             ~ "Not affected"
  )
)

PBS_e269r_vs_PBS_wt.df <- PBS_e269r_vs_PBS_wt.df %>%   mutate(
  PBS_vs_LPS = dplyr::case_when(
    rownames(.) %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE) ~ "Decreasing_PBSvsLPS",
    rownames(.) %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE) ~ "Increasing_PBSvsLPS",
    TRUE                             ~ "Not affected"
  )
)

write.csv(LPS_wt_vs_PBS_wt.df, file="PBS_e269r_vs_PBS_wt.csv")


library(dplyr)
LPS_e269r_vs_PBS_wt.df <- LPS_e269r_vs_PBS_wt.df %>%   mutate(
  Mafb_DKO_BMMs_Cluster_Figure1 = dplyr::case_when(
    rownames(.) %in% Cluster1_genes ~ "Cluster_1",
    rownames(.) %in% Cluster2_genes ~ "Cluster_2",
    TRUE                             ~ "Not affected"
  )
)

LPS_e269r_vs_PBS_wt.df <- LPS_e269r_vs_PBS_wt.df %>%   mutate(
  PBS_vs_LPS = dplyr::case_when(
    rownames(.) %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE) ~ "Decreasing_PBSvsLPS",
    rownames(.) %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE) ~ "Increasing_PBSvsLPS",
    TRUE                             ~ "Not affected"
  )
)

write.csv(LPS_wt_vs_PBS_wt.df, file="LPS_e269r_vs_PBS_wt.csv")

library(dplyr)
LPS_e269r_vs_LPS_wt.df <- LPS_e269r_vs_LPS_wt.df %>%   mutate(
  Mafb_DKO_BMMs_Cluster_Figure1 = dplyr::case_when(
    rownames(.) %in% Cluster1_genes ~ "Cluster_1",
    rownames(.) %in% Cluster2_genes ~ "Cluster_2",
    TRUE                             ~ "Not affected"
  )
)

LPS_e269r_vs_LPS_wt.df <- LPS_e269r_vs_LPS_wt.df %>%   mutate(
  PBS_vs_LPS = dplyr::case_when(
    rownames(.) %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE) ~ "Decreasing_PBSvsLPS",
    rownames(.) %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE) ~ "Increasing_PBSvsLPS",
    TRUE                             ~ "Not affected"
  )
)

write.csv(LPS_e269r_vs_LPS_wt.df, file="LPS_e269r_vs_LPS_wt.csv")

#Now we will visualize log fold change of Mafb specific genes in all the set 

library(dplyr)
library(tidyr)
library(ggplot2)

# rld_df: genes x samples (rlog counts), rownames = genes
# meta: sample metadata with columns: sample, condition
# genes_X: vector of genes

df <- rld_hao_df[rownames(rld_hao_df) %in% LPS_e269r_vs_PBS_wt.df_DECREASE_cluster1, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "Samples", values_to = "rld") %>%
  left_join(Meta_combined, by = "Samples") %>%
  group_by(gene, Condition) %>%
  summarise(mean_rld = mean(rld), sd_rld=sd(rld),.groups = "drop")

  df$Condition<-factor(df$Condition,levels = c("PBS_WT", "PBS_E269R", "LPS_WT", "LPS_E269R"))
ggplot(df, aes(x = Condition, y = mean_rld,fill=Condition)) +
  geom_boxplot()
  


#Heatmap for different sets----
library(ComplexHeatmap)
 #first we just show what happens in PBS----
#heatmap for all the decreasing genes in LPS situation]
rld_hao_df_pbs_wt_het<-rld_hao_df[ ,c("PBS_WT_1","PBS_WT_2","PBS_WT_3","PBS_WT_4","PBS_E269R_1","PBS_E269R_2","PBS_E269R_3","PBS_E269R_4")]
rld_hao_df_lps_cluster1_decrease <- rld_hao_df_pbs_wt_het[rownames(rld_hao_df_pbs_wt_het) %in% PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2, ]
rld_hao_df_lps_cluster1_decrease <- t(scale(t(rld_hao_df_lps_cluster1_decrease)))

#Annotation row
conditions <- rep(c("PBS","PBS_E269R"), each = 4)
annotation_col <- data.frame(
  Condition = factor(conditions),
  row.names = colnames(rld_hao_df_lps_cluster1_decrease)
)
col_ha_c1_decrease_lps_wt <- HeatmapAnnotation(df = annotation_col)


library(colorRamp2)
cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)


ht <- Heatmap(
  rld_hao_df_lps_cluster1_decrease,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  # row_km = 4,
  show_row_names = TRUE,             # hide all labels
  show_column_names = FALSE,
  col = col_fun,
  top_annotation = col_ha_c1_decrease_lps_wt
)

draw(ht)

 #first we just show what happens in lps----
rld_hao_df_lps_wt_het<-rld_hao_df[ ,c("LPS_WT_1","LPS_WT_2","LPS_WT_3","LPS_WT_4","LPS_E269R_1","LPS_E269R_2","LPS_E269R_3","LPS_E269R_4")]
rld_hao_df_lps_cluster1_decrease <- rld_hao_df_lps_wt_het[rownames(rld_hao_df_lps_wt_het) %in% LPS_e269r_vs_LPS_wt.df_DECREASE_cluster2, ]
rld_hao_df_lps_cluster1_decrease <- t(scale(t(rld_hao_df_lps_cluster1_decrease)))

#Annotation row
conditions <- rep(c("LPS","LPS_E269R"), each = 4)
annotation_col <- data.frame(
  Condition = factor(conditions),
  row.names = colnames(rld_hao_df_lps_cluster1_decrease)
)
col_ha_c1_decrease_lps_wt <- HeatmapAnnotation(df = annotation_col)


library(colorRamp2)
cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)


ht <- Heatmap(
  rld_hao_df_lps_cluster1_decrease,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  # row_km = 4,
  show_row_names = TRUE,             # hide all labels
  show_column_names = FALSE,
  col = col_fun,
  top_annotation = col_ha_c1_decrease_lps_wt
)

draw(ht)

 #Now for combined condition and we see the genes which were downregulating in LPS, how do they look like in het in LPS----

#heatmap for all the decreasing genes in LPS situation
rld_hao_df<-rld_hao_df[ ,c("PBS_WT_1","PBS_WT_2","PBS_WT_3","PBS_WT_4","PBS_E269R_1","PBS_E269R_2","PBS_E269R_3","PBS_E269R_4",
                                          "LPS_WT_1","LPS_WT_2","LPS_WT_3","LPS_WT_4","LPS_E269R_1","LPS_E269R_2","LPS_E269R_3","LPS_E269R_4")]

rld_hao_df_lps_cluster1_decrease <- rld_hao_df[rownames(rld_hao_df) %in% LPS_wt_vs_PBS_wt.df_DECREASE_cluster2, ]
rld_hao_df_lps_cluster1_decrease <- t(scale(t(rld_hao_df_lps_cluster1_decrease)))


#Annotation row
conditions <- rep(c("PBS","PBS_E269R","LPS","LPS_E269R"), each = 4)
annotation_col <- data.frame(
  Condition = factor(conditions),
  row.names = colnames(rld_hao_df_lps_cluster1_decrease)
)
col_ha_c1_decrease_lps_wt <- HeatmapAnnotation(df = annotation_col)


library(colorRamp2)
cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)


ht <- Heatmap(
  rld_hao_df_lps_cluster1_decrease,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_km = 4,
  show_row_names = TRUE,             # hide all labels
  show_column_names = FALSE,
  col = col_fun,
  top_annotation = col_ha_c1_decrease_lps_wt
)

draw(ht)


#Comparative gene ontology----

#I think its better to use compare cluster 
#we want to see the property of these genes which are favored by Maf-Maf complex or Maf-Fos complex

library(BuenColors)
library(org.Mm.eg.db)
library(clusterProfiler)

PBS_mafb_maf_het_comparecluster_data<-PBS_e269r_vs_PBS_wt.df %>% dplyr::select(padj,log2FoldChange,Mafb_DKO_BMMs_Cluster_Figure1) %>% filter(Mafb_DKO_BMMs_Cluster_Figure1 %in% c("Cluster_1","Cluster_2"))
PBS_mafb_maf_het_comparecluster_data <- PBS_mafb_maf_het_comparecluster_data %>% filter(padj<0.05)
PBS_mafb_maf_het_comparecluster_data <- PBS_mafb_maf_het_comparecluster_data %>% mutate(Gene = rownames(PBS_mafb_maf_het_comparecluster_data)) 
PBS_mafb_maf_het_comparecluster_data <- PBS_mafb_maf_het_comparecluster_data %>% mutate(ENTREZ = mapIds(org.Mm.eg.db,PBS_mafb_maf_het_comparecluster_data$Gene, "ENTREZID", "SYMBOL")) 
PBS_mafb_maf_het_comparecluster_data$Direction <- "upregulated"
PBS_mafb_maf_het_comparecluster_data$Direction[PBS_mafb_maf_het_comparecluster_data$log2FoldChange < 0] <- "downregulated"
lookup <- c(FC = "log2FoldChange", Entrez = "ENTREZ")
PBS_mafb_maf_het_comparecluster_data<-rename(PBS_mafb_maf_het_comparecluster_data, all_of(lookup))

PBS_mafb_maf_het_comparecluster <- compareCluster(Entrez~Direction+Mafb_DKO_BMMs_Cluster_Figure1, data=PBS_mafb_maf_het_comparecluster_data, 
                                           fun=enrichGO, OrgDb = org.Mm.eg.db, ont="BP")

head(PBS_mafb_maf_het_comparecluster )
dotplot(PBS_mafb_maf_het_comparecluster, x="Direction",showCategory = 6,size = "zScore", color = "p.adjust") + facet_grid(~Mafb_DKO_BMMs_Cluster_Figure1)

#lps
LPS_mafb_maf_het_comparecluster_data<-LPS_e269r_vs_LPS_wt.df %>% dplyr::select(padj,log2FoldChange,Mafb_DKO_BMMs_Cluster_Figure1) %>% filter(Mafb_DKO_BMMs_Cluster_Figure1 %in% c("Cluster_1","Cluster_2"))
LPS_mafb_maf_het_comparecluster_data <- LPS_mafb_maf_het_comparecluster_data %>% filter(padj<0.05)
LPS_mafb_maf_het_comparecluster_data <- LPS_mafb_maf_het_comparecluster_data %>% mutate(Gene = rownames(LPS_mafb_maf_het_comparecluster_data)) 
LPS_mafb_maf_het_comparecluster_data <- LPS_mafb_maf_het_comparecluster_data %>% mutate(ENTREZ = mapIds(org.Mm.eg.db,LPS_mafb_maf_het_comparecluster_data$Gene, "ENTREZID", "SYMBOL")) 
LPS_mafb_maf_het_comparecluster_data$Direction <- "upregulated"
LPS_mafb_maf_het_comparecluster_data$Direction[LPS_mafb_maf_het_comparecluster_data$log2FoldChange < 0] <- "downregulated"
lookup <- c(FC = "log2FoldChange", Entrez = "ENTREZ")
LPS_mafb_maf_het_comparecluster_data<-rename(LPS_mafb_maf_het_comparecluster_data, all_of(lookup))

LPS_mafb_maf_het_comparecluster_data <- compareCluster(Entrez~Direction+Mafb_DKO_BMMs_Cluster_Figure1, data=LPS_mafb_maf_het_comparecluster_data, 
                                                  fun=enrichGO, OrgDb = org.Mm.eg.db, ont="BP")

head(PBS_mafb_maf_het_comparecluster )
dotplot(LPS_mafb_maf_het_comparecluster_data, x="Direction",showCategory = 6,size = "zScore", color = "p.adjust") + facet_grid(~Mafb_DKO_BMMs_Cluster_Figure1)

#Comparison of Wt vs het signature and similarity to known immunosuppresive signature----
#Make microarray dataset for homo het as a model and comparison with invivo

library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(BuenColors)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(org.Mm.eg.db)
library(stringr)

#Analysing microarray in R----- 
library(limma)
library(Biobase)
library(affycoretools)
#Useful links:
#https://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf
#https://www.bioconductor.org/packages/release/bioc/vignettes/affy/inst/doc/builtinMethods.pdf

Microarray<-as.data.frame(Microarray)
Microarray_meta<-as.data.frame(Microarray_meta)
row.names(Microarray_meta)<-Microarray_meta$Sample
Microarray_meta<-Microarray_meta[,-1]
Microarray_meta$Condition<-as.factor(Microarray_meta$Condition)

f_microarray <- factor(Microarray_meta$Condition, levels=c(
  "DKO_FOS_NoDOX", "DKO_MafB_NoDOX", "DKO_Mafb_FOS_NoDOX", "DKO_E269R_NoDOX", "DKO_E269R_FOS_NoDOX",
  "DKO_FOS_DOX_6", "DKO_MafB_DOX_6", "DKO_Mafb_FOS_DOX_6", "DKO_E269R_DOX_6", "DKO_E269R_FOS_DOX_6",
  "DKO_FOS_DOX_10", "DKO_MafB_DOX_10", "DKO_Mafb_FOS_DOX_10", "DKO_E269R_DOX_10", "DKO_E269R_FOS_DOX_10"
))
design <- model.matrix(~ 0+ f_microarray)
colnames(design)<-c(
  "DKO_FOS_NoDOX", "DKO_MafB_NoDOX", "DKO_Mafb_FOS_NoDOX", "DKO_E269R_NoDOX", "DKO_E269R_FOS_NoDOX",
  "DKO_FOS_DOX_6", "DKO_MafB_DOX_6", "DKO_Mafb_FOS_DOX_6", "DKO_E269R_DOX_6", "DKO_E269R_FOS_DOX_6",
  "DKO_FOS_DOX_10", "DKO_MafB_DOX_10", "DKO_Mafb_FOS_DOX_10", "DKO_E269R_DOX_10", "DKO_E269R_FOS_DOX_10"
)


contrast.Matrix_microarray <- makeContrasts(
  FOS_DOX_10_vs_NoDOX = DKO_FOS_DOX_10 - DKO_FOS_NoDOX,
  FOS_DOX_6_vs_NoDOX = DKO_FOS_DOX_6 - DKO_FOS_NoDOX,
  MafB_DOX_10_vs_NoDOX = DKO_MafB_DOX_10 - DKO_MafB_NoDOX,
  MafB_DOX_6_vs_NoDOX = DKO_MafB_DOX_6 - DKO_MafB_NoDOX,
  Mafb_FOS_DOX_10_vs_NoDOX = DKO_Mafb_FOS_DOX_10 - DKO_Mafb_FOS_NoDOX,
  Mafb_FOS_DOX_6_vs_NoDOX = DKO_Mafb_FOS_DOX_6 - DKO_Mafb_FOS_NoDOX,
  E269R_DOX_10_vs_NoDOX = DKO_E269R_DOX_10 - DKO_E269R_NoDOX,
  E269R_DOX_6_vs_NoDOX = DKO_E269R_DOX_6 - DKO_E269R_NoDOX,
  E269R_FOS_DOX_10_vs_NoDOX = DKO_E269R_FOS_DOX_10 - DKO_E269R_FOS_NoDOX,
  E269R_FOS_DOX_6_vs_NoDOX = DKO_E269R_FOS_DOX_6 - DKO_E269R_FOS_NoDOX,
  levels = design
)

expression_set_microarray <- ExpressionSet(assayData = as.matrix(Microarray),
                                           design = design)
library(pd.mogene.1.0.st.v1)
#expression_set_microarray_annotated<-annotateEset(expression_set_microarray , mogene10sttranscriptcluster.db)
expression_set_microarray_annotated<-annotateEset(expression_set_microarray , pd.mogene.1.0.st.v1)

?annotateEset
?annotate
fit_microarray_annotated <- lmFit(expression_set_microarray_annotated, design)
fit_microarray_annotated <- eBayes(fit_microarray_annotated)
colnames(fit_microarray_annotated)
summary(fit_microarray_annotated)

fit2_microarray_annotated <- contrasts.fit(fit_microarray_annotated, contrast.Matrix_microarray)
fit2_microarray_annotated <- eBayes(fit2_microarray_annotated)
results_Microarray_annotated_all <- decideTests(fit2_microarray_annotated)

plot(fit2_microarray_annotated$df.residual)

onlyFos_10hr_annotated<-topTable(fit2_microarray_annotated, coef=1, adjust="BH", number = 35556)
onlyFos_6hr_annotated<-topTable(fit2_microarray_annotated, coef=2, adjust="BH", number = 35556)
Mafb_10hr_annotated<-topTable(fit2_microarray_annotated, coef=3, adjust="BH", number = 35556)
Mafb_6hr_annotated<-topTable(fit2_microarray_annotated, coef=4, adjust="BH", number = 35556)
Mafb_Fos_10hr_annotated<-topTable(fit2_microarray_annotated, coef=5, adjust="BH", number = 35556)
Mafb_Fos_6hr_annotated<-topTable(fit2_microarray_annotated, coef=6, adjust="BH", number = 35556)
E269R_10hr_annotated<-topTable(fit2_microarray_annotated, coef=7, adjust="BH", number = 35556)
E269R_6hr_annotated<-topTable(fit2_microarray_annotated, coef=8, adjust="BH", number = 35556)
E269R_Fos_annotated_10hr<-topTable(fit2_microarray_annotated, coef=9, adjust="BH", number = 35556)
E269R_Fos_annotated_6hr<-topTable(fit2_microarray_annotated, coef=10, adjust="BH", number = 35556)

?topTable
#Significant 
library(dplyr)
onlyFos_10hr_annotated_significant<-filter(onlyFos_10hr_annotated, onlyFos_10hr_annotated$adj.P.Val<0.05
                                           & abs(onlyFos_10hr_annotated$logFC)>0.0)
onlyFos_6hr_annotated_significant<-filter(onlyFos_6hr_annotated, onlyFos_6hr_annotated$adj.P.Val<0.05
                                          & abs(onlyFos_6hr_annotated$logFC)>0.0)
Mafb_10hr_annotated_significant<-filter(Mafb_10hr_annotated, Mafb_10hr_annotated$adj.P.Val<0.05
                                        & abs(Mafb_10hr_annotated$logFC)>0.0)
Mafb_6hr_annotated_significant<-filter(Mafb_6hr_annotated, Mafb_6hr_annotated$adj.P.Val<0.05
                                       & abs(Mafb_6hr_annotated$logFC)>0.0)
Mafb_Fos_10hr_annotated_significant<-filter(Mafb_Fos_10hr_annotated, Mafb_Fos_10hr_annotated$adj.P.Val<0.05
                                            & abs(Mafb_Fos_10hr_annotated$logFC)>0.0)
Mafb_Fos_6hr_annotated_significant<-filter(Mafb_Fos_6hr_annotated, Mafb_Fos_6hr_annotated$adj.P.Val<0.05
                                           & abs(Mafb_Fos_6hr_annotated$logFC)>0.0)
E269R_10hr_annotated_significant<-filter(E269R_10hr_annotated, E269R_10hr_annotated$adj.P.Val<0.05
                                         & abs(E269R_10hr_annotated$logFC)>0.0)
E269R_6hr_annotated_significant<-filter(E269R_6hr_annotated, E269R_6hr_annotated$adj.P.Val<0.05
                                        & abs(E269R_6hr_annotated$logFC)>0.0)
E269R_Fos_annotated_10hr_significant<-filter(E269R_Fos_annotated_10hr, E269R_Fos_annotated_10hr$adj.P.Val<0.05
                                             & abs(E269R_Fos_annotated_10hr$logFC)>0.0)
E269R_Fos_annotated_6hr_significant<-filter(E269R_Fos_annotated_6hr, E269R_Fos_annotated_6hr$adj.P.Val<0.05
                                            & abs(E269R_Fos_annotated_6hr$logFC)>0.0)


library(ggplot2)
library(dplyr)
library(patchwork)


plot_volcano <- function(df, title, ymax_cut = 10){
  
  df <- df %>%
    mutate(Regulation = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Up",
      adj.P.Val < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    neglog10 = -log10(adj.P.Val),
    neglog10_cap = ifelse(neglog10 > ymax_cut, ymax_cut, neglog10))
  
  # Count genes
  up_n <- sum(df$Regulation == "Up", na.rm = TRUE)
  down_n <- sum(df$Regulation == "Down", na.rm = TRUE)
  
  ggplot(df, aes(logFC, neglog10_cap, color = Regulation)) +
    geom_point(size = 1, alpha = 0.7) +
    scale_color_manual(values = c(
      Up = "red",
      Down = "blue",
      NS = "grey80"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype="dashed") +
    
    annotate("text",
             x = max(df$logFC, na.rm = TRUE),
             y = ymax_cut * 0.5,
             label = paste0("Up: ", up_n),
             hjust = 1, color = "red", size = 4) +
    
    annotate("text",
             x = min(df$logFC, na.rm = TRUE),
             y = ymax_cut * 0.5,
             label = paste0("Down: ", down_n),
             hjust = 0, color = "blue", size = 4) +
    
    coord_cartesian(ylim = c(0, ymax_cut)) +
    
    theme_classic() +
    labs(title = title,
         x = "logFC",
         y = "-log10(adj P)")
}

p1 <- plot_volcano(onlyFos_6hr_annotated, "onlyFos 6hr")
p2 <- plot_volcano(Mafb_6hr_annotated, "Mafb 6hr")
p3 <- plot_volcano(Mafb_Fos_6hr_annotated, "Mafb+Fos 6hr")
p4 <- plot_volcano(E269R_6hr_annotated, "E269R 6hr")
p5 <- plot_volcano(E269R_Fos_annotated_6hr, "E269R+Fos 6hr")

p6 <- plot_volcano(onlyFos_10hr_annotated, "onlyFos 10hr")
p7 <- plot_volcano(Mafb_10hr_annotated, "Mafb 10hr")
p8 <- plot_volcano(Mafb_Fos_10hr_annotated, "Mafb+Fos 10hr")
p9 <- plot_volcano(E269R_10hr_annotated, "E269R 10hr")
p10 <- plot_volcano(E269R_Fos_annotated_10hr, "E269R+Fos 10hr")

final_plot <- (p1 | p2 | p3 | p4 | p5) /
  (p6 | p7 | p8 | p9 | p10)

final_plot

ggsave("volcano_all_conditions_microarray.pdf",
       final_plot,
       width = 16,
       height = 6)

#Heatmap for all genes upregulated by Mafb, how do they change in E269R plus fos----

#For all significant genes heatmap 
Mafb_10hr_annotated_significant_genes<-Mafb_10hr_annotated_significant$SYMBOL
Mafb_6hr_annotated_significant_genes<-Mafb_6hr_annotated_significant$SYMBOL
Mafb_annotated_significant_genes<-union(Mafb_10hr_annotated_significant_genes,Mafb_6hr_annotated_significant_genes)

# --- 1. Combine expression data with annotation ---
# Extract the expression matrix
microarray_matrix <- exprs(expression_set_microarray_annotated)
# Convert to a data.frame if you prefer
Microarry_annotated_matrix <- as.data.frame(microarray_matrix)
feature_info <- fData(expression_set_microarray_annotated)
expr_df <- Microarry_annotated_matrix
expr_df$GeneSymbol <- feature_info$SYMBOL[match(rownames(expr_df), rownames(feature_info))]

# --- 2. Remove rows with missing gene names (optional) ---
expr_df <- expr_df[!is.na(expr_df$GeneSymbol) & expr_df$GeneSymbol != "", ]

# --- 3. Average expression per gene ---
# use dplyr for easy grouping
library(dplyr)
gene_expression_df <- expr_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
# --- 4. Make gene names rownames ---
gene_expression_df <- as.data.frame(gene_expression_df)
rownames(gene_expression_df) <- gene_expression_df$GeneSymbol
gene_expression_df$GeneSymbol <- NULL

gene_expression_Mafb_significant <- gene_expression_df[rownames(gene_expression_df) %in% Mafb_annotated_significant_genes, ]
gene_expression_Mafb_significant <- gene_expression_Mafb_significant[ ,c("wt_0h_rep_1","wt_0h_rep_2","wt_0h_rep_3",
                                                                         "wt_6h_rep_1","wt_6h_rep_2","wt_6h_rep_3",
                                                                         "wt_10h_rep_1","wt_10h_rep_2","wt_10h_rep_3",
                                                                         "E269R_plus_Fos_0h_rep_1","E269R_plus_Fos_0h_rep_2","E269R_plus_Fos_0h_rep_3",
                                                                         "E269R_plus_Fos_6h_rep_1","E269R_plus_Fos_6h_rep_2","E269R_plus_Fos_6h_rep_3",
                                                                         "E269R_plus_Fos_10h_rep_1","E269R_plus_Fos_10h_rep_2","E269R_plus_Fos_10h_rep_3")]



#Making heatmap
# --- Example: heatmap with clustering into groups ---
scaled_mat_gene_expression_Mafb_significant <- t(scale(t(gene_expression_Mafb_significant)))   # scale each gene (row) across samples

# Combine with scaled matrix
scaled_df <- as.data.frame(scaled_mat_gene_expression_Mafb_significant)
scaled_df$Cluster <- gene_clusters$Cluster

# Compute mean per cluster per sample
cluster_mean_scaled <- aggregate(. ~ Cluster, data = scaled_df, FUN = mean)
rownames(cluster_mean_scaled) <- paste0("Cluster_", cluster_mean_scaled$Cluster)
cluster_mean_scaled$Cluster <- NULL

# View result
cluster_mean_scaled

#Annotation row
conditions <- rep(c("NT", "6hr_Mafb", "10hr_Mafb","NT", "6hr_E269R_plus_Fos", "10hr_E269R_plus_Fos"), each = 3)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(gene_expression_Mafb_significant)
col_ha <- HeatmapAnnotation(df = annotation_col)


#annotation_row <- data.frame(Cluster = factor(gene_clusters$Cluster))
#rownames(annotation_row) <- rownames(gene_expression_Mafb_significant)
#row_ha <- rowAnnotation(df = annotation_row)

all_genes <- rownames(scaled_mat_gene_expression_Mafb_significant)

PBS_LPS_status <- ifelse(all_genes %in% rownames(LPS_wt_vs_PBS_wt.df_INCREASE), "Up",
                         ifelse(all_genes %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE), "Down", "NS"))

LPS_LPShet_status <- ifelse(all_genes %in% rownames(LPS_e269r_vs_LPS_wt.df_INCREASE), "Up",
                            ifelse(all_genes %in% rownames(LPS_e269r_vs_LPS_wt.df_DECREASE), "Down", "NS"))

PBS_PBShet_status <- ifelse(all_genes %in% rownames(PBS_e269r_vs_PBS_wt.df_INCREASE), "Up",
                            ifelse(all_genes %in% rownames(PBS_e269r_vs_PBS_wt.df_DECREASE), "Down", "NS"))

annotation_row <- data.frame(
  Cluster = factor(gene_clusters$Cluster),
  PBS_vs_LPS = PBS_LPS_status,
  LPS_vs_LPShet = LPS_LPShet_status,
  PBS_vs_PBShet = PBS_PBShet_status
)

rownames(annotation_row) <- rownames(scaled_mat_gene_expression_Mafb_significant)

ann_colors <- list(
  
  Cluster = c("1" = "#1B9E77",
              "2" = "#D95F02"),
  
  PBS_vs_LPS = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  ),
  
  LPS_vs_LPShet = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  ),
  
  PBS_vs_PBShet = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  )
)

row_ha <- rowAnnotation(
  df = annotation_row,
  col = ann_colors
)
#Labels
set.seed(123)
highlight_genes<-LPS_wt_vs_PBS_wt.df_DECREASE_cluster1

highlight_genes<-highlight_genes[highlight_genes %in% rownames(scaled_df)]
library(BuenColors)

cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)

 ht <- Heatmap(
  scaled_mat_gene_expression_Mafb_significant,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = 2,
  show_row_names = FALSE,
  show_column_names = TRUE,
  col = col_fun,
  left_annotation = row_ha,
  top_annotation = col_ha
)

# --- Add labels outside with connecting lines ---
lg <- rowAnnotation(
  labels = anno_mark(
    at = match(highlight_genes, rownames(scaled_mat_gene_expression_Mafb_significant)),
    labels = highlight_genes,
    labels_gp = gpar(fontsize = 12),
    link_gp = gpar(lwd = 1.2)
  )
)

# Draw everything
draw(ht + lg, merge_legend = TRUE)

#Make the boxplots 
df_long <- scaled_df %>%
  rownames_to_column("gene") %>%       # move rownames to a real column
  pivot_longer(
    cols = starts_with("wt_"),         # select the wt_ columns
    names_to = "sample",               
    values_to = "expression"
  )





#Heatmap for genes goign down after LPS (CLUSTER1), how are they affected in E269R plus Fos ----

# --- 1. Combine expression data with annotation ---
# Extract the expression matrix
microarray_matrix <- exprs(expression_set_microarray_annotated)
# Convert to a data.frame if you prefer
Microarry_annotated_matrix <- as.data.frame(microarray_matrix)
feature_info <- fData(expression_set_microarray_annotated)
expr_df <- Microarry_annotated_matrix
expr_df$GeneSymbol <- feature_info$SYMBOL[match(rownames(expr_df), rownames(feature_info))]

# --- 2. Remove rows with missing gene names (optional) ---
expr_df <- expr_df[!is.na(expr_df$GeneSymbol) & expr_df$GeneSymbol != "", ]

# --- 3. Average expression per gene ---
# use dplyr for easy grouping
library(dplyr)
gene_expression_df <- expr_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
# --- 4. Make gene names rownames ---
gene_expression_df <- as.data.frame(gene_expression_df)
rownames(gene_expression_df) <- gene_expression_df$GeneSymbol
gene_expression_df$GeneSymbol <- NULL

gene_expression_Mafb_significant <- gene_expression_df[rownames(gene_expression_df) %in% LPS_wt_vs_PBS_wt.df_DECREASE_cluster1, ]
gene_expression_Mafb_significant <- gene_expression_Mafb_significant[ ,c("wt_0h_rep_1","wt_0h_rep_2","wt_0h_rep_3",
                                                                         "wt_6h_rep_1","wt_6h_rep_2","wt_6h_rep_3",
                                                                         "wt_10h_rep_1","wt_10h_rep_2","wt_10h_rep_3",
                                                                         "E269R_plus_Fos_0h_rep_1","E269R_plus_Fos_0h_rep_2","E269R_plus_Fos_0h_rep_3",
                                                                         "E269R_plus_Fos_6h_rep_1","E269R_plus_Fos_6h_rep_2","E269R_plus_Fos_6h_rep_3",
                                                                         "E269R_plus_Fos_10h_rep_1","E269R_plus_Fos_10h_rep_2","E269R_plus_Fos_10h_rep_3")]



#Making heatmap
# --- Example: heatmap with clustering into groups ---
scaled_mat_gene_expression_Mafb_significant <- t(scale(t(gene_expression_Mafb_significant)))   # scale each gene (row) across samples

#Annotation row
conditions <- rep(c("NT", "6hr_Mafb", "10hr_Mafb","NT", "6hr_E269R_plus_Fos", "10hr_E269R_plus_Fos"), each = 3)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(gene_expression_Mafb_significant)
col_ha <- HeatmapAnnotation(df = annotation_col)


#annotation_row <- data.frame(Cluster = factor(gene_clusters$Cluster))
#rownames(annotation_row) <- rownames(gene_expression_Mafb_significant)
#row_ha <- rowAnnotation(df = annotation_row)

all_genes <- rownames(scaled_mat_gene_expression_Mafb_significant)

PBS_LPS_status <- ifelse(all_genes %in% rownames(LPS_wt_vs_PBS_wt.df_INCREASE), "Up",
                         ifelse(all_genes %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE), "Down", "NS"))

LPS_LPShet_status <- ifelse(all_genes %in% rownames(LPS_e269r_vs_LPS_wt.df_INCREASE), "Up",
                            ifelse(all_genes %in% rownames(LPS_e269r_vs_LPS_wt.df_DECREASE), "Down", "NS"))

PBS_PBShet_status <- ifelse(all_genes %in% rownames(PBS_e269r_vs_PBS_wt.df_INCREASE), "Up",
                            ifelse(all_genes %in% rownames(PBS_e269r_vs_PBS_wt.df_DECREASE), "Down", "NS"))

annotation_row <- data.frame(
  PBS_vs_LPS = PBS_LPS_status,
  LPS_vs_LPShet = LPS_LPShet_status,
  PBS_vs_PBShet = PBS_PBShet_status
)

rownames(annotation_row) <- rownames(scaled_mat_gene_expression_Mafb_significant)

ann_colors <- list(
  PBS_vs_LPS = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  ),
  
  LPS_vs_LPShet = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  ),
  
  PBS_vs_PBShet = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  )
)

row_ha <- rowAnnotation(
  df = annotation_row,
  col = ann_colors
)
#Labels
set.seed(123)
highlight_genes<-LPS_wt_vs_PBS_wt.df_DECREASE_cluster1

highlight_genes<-highlight_genes[highlight_genes %in% rownames(scaled_df)]
library(BuenColors)

cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)

ht <- Heatmap(
  scaled_mat_gene_expression_Mafb_significant,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = col_fun,
  left_annotation = row_ha,
  top_annotation = col_ha
)


# Draw everything
draw(ht)
scaled_df <- as.data.frame(scaled_mat_gene_expression_Mafb_significant)
df_long <- scaled_df %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = matches("wt_|E269R"),
    names_to = "sample",
    values_to = "expression"
  )

df_long <- df_long %>%
  mutate(
    condition = ifelse(grepl("^wt_", sample), "WT", "E269R"),
    
    timepoint = case_when(
      grepl("_0h_", sample) ~ "0h",
      grepl("_6h", sample) ~ "6h",
      grepl("_10h_", sample) ~ "10h"
    ),
    
    replicate = case_when(
      grepl("rep_1", sample) ~ "rep1",
      grepl("rep_2", sample) ~ "rep2",
      grepl("rep_3", sample) ~ "rep3"
    )
  )
df_long$timepoint <- factor(df_long$timepoint, levels = c("0h","6h","10h"))
summary_tab <- df_long %>%
  group_by(condition, timepoint) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    sd_expr = sd(expression, na.rm = TRUE),
    .groups = "drop"
  )

library(BuenColors)

pal <- jdb_palette("solar_flare", n = 100, type = "continuous")

ggplot(summary_tab,
       aes(x = timepoint,
           y = mean_expr,
           group = condition,
           color = mean_expr,
           fill = mean_expr)) +
  
  geom_line(size = 1.4) +
  geom_point(size = 3) +
  
  geom_ribbon(
    aes(ymin = mean_expr - sd_expr,
        ymax = mean_expr + sd_expr,
        group = condition),
    alpha = 0.25,
    color = NA
  ) +
  
  scale_color_gradientn(colours = pal, name = "Mean expression") +
  scale_fill_gradientn(colours = pal, name = "Mean expression") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 12)
  ) +
  
  labs(
    x = "Timepoint",
    y = "Mean scaled expression",
    title = "Mean expression trajectory of heatmap genes\nWT vs E269R"
  )

#boxplot
library(tidyverse)
library(BuenColors)

scaled_df <- as.data.frame(scaled_mat_gene_expression_Mafb_significant)

df_long <- scaled_df %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = matches("wt_|E269R"),
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    condition = ifelse(grepl("^wt_", sample), "WT", "E269R"),
    timepoint = case_when(
      grepl("_0h_", sample) ~ "0h",
      grepl("_6h_", sample) ~ "6h",
      grepl("_10h_", sample) ~ "10h"
    ),
    replicate = case_when(
      grepl("rep_1", sample) ~ "rep1",
      grepl("rep_2", sample) ~ "rep2",
      grepl("rep_3", sample) ~ "rep3"
    )
  )

df_long$timepoint <- factor(df_long$timepoint, levels = c("0h","6h","10h"))


ggplot(df_long, aes(x = timepoint, y = expression, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               position = position_dodge(width = 0.8), fill = "white") +
  
  scale_fill_manual(values = c("WT" = "#4DBBD5", "E269R" = "#E64B35")) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) +
  
  labs(
    x = "Timepoint",
    y = "Expression (scaled)",
    fill = "Condition",
    title = "Distribution of heatmap gene expression\nWT vs E269R"
  )






#Heatmap for genes goign down after LPS in Het mouse, how are they affected in E269R plus Fos in our system----

# --- 1. Combine expression data with annotation ---
# Extract the expression matrix
microarray_matrix <- exprs(expression_set_microarray_annotated)
# Convert to a data.frame if you prefer
Microarry_annotated_matrix <- as.data.frame(microarray_matrix)
feature_info <- fData(expression_set_microarray_annotated)
expr_df <- Microarry_annotated_matrix
expr_df$GeneSymbol <- feature_info$SYMBOL[match(rownames(expr_df), rownames(feature_info))]

# --- 2. Remove rows with missing gene names (optional) ---
expr_df <- expr_df[!is.na(expr_df$GeneSymbol) & expr_df$GeneSymbol != "", ]

# --- 3. Average expression per gene ---
# use dplyr for easy grouping
library(dplyr)
gene_expression_df <- expr_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
# --- 4. Make gene names rownames ---
gene_expression_df <- as.data.frame(gene_expression_df)
rownames(gene_expression_df) <- gene_expression_df$GeneSymbol
gene_expression_df$GeneSymbol <- NULL

gene_expression_Mafb_significant <- gene_expression_df[rownames(gene_expression_df) %in% LPS_wt_vs_PBS_wt.df_DECREASE_cluster1, ]
gene_expression_Mafb_significant <- gene_expression_Mafb_significant[ ,c("wt_0h_rep_1","wt_0h_rep_2","wt_0h_rep_3",
                                                                         "wt_6h_rep_1","wt_6h_rep_2","wt_6h_rep_3",
                                                                         "wt_10h_rep_1","wt_10h_rep_2","wt_10h_rep_3",
                                                                         "E269R_plus_Fos_0h_rep_1","E269R_plus_Fos_0h_rep_2","E269R_plus_Fos_0h_rep_3",
                                                                         "E269R_plus_Fos_6h_rep_1","E269R_plus_Fos_6h_rep_2","E269R_plus_Fos_6h_rep_3",
                                                                         "E269R_plus_Fos_10h_rep_1","E269R_plus_Fos_10h_rep_2","E269R_plus_Fos_10h_rep_3")]



#Making heatmap
# --- Example: heatmap with clustering into groups ---
scaled_mat_gene_expression_Mafb_significant <- t(scale(t(gene_expression_Mafb_significant)))   # scale each gene (row) across samples

#Annotation row
conditions <- rep(c("NT", "6hr_Mafb", "10hr_Mafb","NT", "6hr_E269R_plus_Fos", "10hr_E269R_plus_Fos"), each = 3)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(gene_expression_Mafb_significant)
col_ha <- HeatmapAnnotation(df = annotation_col)


#annotation_row <- data.frame(Cluster = factor(gene_clusters$Cluster))
#rownames(annotation_row) <- rownames(gene_expression_Mafb_significant)
#row_ha <- rowAnnotation(df = annotation_row)

all_genes <- rownames(scaled_mat_gene_expression_Mafb_significant)

PBS_LPS_status <- ifelse(all_genes %in% rownames(LPS_wt_vs_PBS_wt.df_INCREASE), "Up",
                         ifelse(all_genes %in% rownames(LPS_wt_vs_PBS_wt.df_DECREASE), "Down", "NS"))

LPS_LPShet_status <- ifelse(all_genes %in% rownames(LPS_e269r_vs_LPS_wt.df_INCREASE), "Up",
                            ifelse(all_genes %in% rownames(LPS_e269r_vs_LPS_wt.df_DECREASE), "Down", "NS"))

PBS_PBShet_status <- ifelse(all_genes %in% rownames(PBS_e269r_vs_PBS_wt.df_INCREASE), "Up",
                            ifelse(all_genes %in% rownames(PBS_e269r_vs_PBS_wt.df_DECREASE), "Down", "NS"))

annotation_row <- data.frame(
  PBS_vs_LPS = PBS_LPS_status,
  LPS_vs_LPShet = LPS_LPShet_status,
  PBS_vs_PBShet = PBS_PBShet_status
)

rownames(annotation_row) <- rownames(scaled_mat_gene_expression_Mafb_significant)

ann_colors <- list(
  PBS_vs_LPS = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  ),
  
  LPS_vs_LPShet = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  ),
  
  PBS_vs_PBShet = c(
    Up = "red",
    Down = "blue",
    NS = "grey90"
  )
)

row_ha <- rowAnnotation(
  df = annotation_row,
  col = ann_colors
)
#Labels
set.seed(123)
highlight_genes<-LPS_wt_vs_PBS_wt.df_DECREASE_cluster1

highlight_genes<-highlight_genes[highlight_genes %in% rownames(scaled_df)]
library(BuenColors)

cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)

ht <- Heatmap(
  scaled_mat_gene_expression_Mafb_significant,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = col_fun,
  left_annotation = row_ha,
  top_annotation = col_ha
)


# Draw everything
draw(ht)
scaled_df <- as.data.frame(scaled_mat_gene_expression_Mafb_significant)
df_long <- scaled_df %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = matches("wt_|E269R"),
    names_to = "sample",
    values_to = "expression"
  )

df_long <- df_long %>%
  mutate(
    condition = ifelse(grepl("^wt_", sample), "WT", "E269R"),
    
    timepoint = case_when(
      grepl("_0h_", sample) ~ "0h",
      grepl("_6h", sample) ~ "6h",
      grepl("_10h_", sample) ~ "10h"
    ),
    
    replicate = case_when(
      grepl("rep_1", sample) ~ "rep1",
      grepl("rep_2", sample) ~ "rep2",
      grepl("rep_3", sample) ~ "rep3"
    )
  )
df_long$timepoint <- factor(df_long$timepoint, levels = c("0h","6h","10h"))
summary_tab <- df_long %>%
  group_by(condition, timepoint) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    sd_expr = sd(expression, na.rm = TRUE),
    .groups = "drop"
  )

library(BuenColors)

pal <- jdb_palette("solar_flare", n = 100, type = "continuous")

ggplot(summary_tab,
       aes(x = timepoint,
           y = mean_expr,
           group = condition,
           color = mean_expr,
           fill = mean_expr)) +
  
  geom_line(size = 1.4) +
  geom_point(size = 3) +
  
  geom_ribbon(
    aes(ymin = mean_expr - sd_expr,
        ymax = mean_expr + sd_expr,
        group = condition),
    alpha = 0.25,
    color = NA
  ) +
  
  scale_color_gradientn(colours = pal, name = "Mean expression") +
  scale_fill_gradientn(colours = pal, name = "Mean expression") +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 12)
  ) +
  
  labs(
    x = "Timepoint",
    y = "Mean scaled expression",
    title = "Mean expression trajectory of heatmap genes\nWT vs E269R"
  )

#boxplot
library(tidyverse)
library(BuenColors)

scaled_df <- as.data.frame(scaled_mat_gene_expression_Mafb_significant)

df_long <- scaled_df %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = matches("wt_|E269R"),
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    condition = ifelse(grepl("^wt_", sample), "WT", "E269R"),
    timepoint = case_when(
      grepl("_0h_", sample) ~ "0h",
      grepl("_6h_", sample) ~ "6h",
      grepl("_10h_", sample) ~ "10h"
    ),
    replicate = case_when(
      grepl("rep_1", sample) ~ "rep1",
      grepl("rep_2", sample) ~ "rep2",
      grepl("rep_3", sample) ~ "rep3"
    )
  )

df_long$timepoint <- factor(df_long$timepoint, levels = c("0h","6h","10h"))


ggplot(df_long, aes(x = timepoint, y = expression, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               position = position_dodge(width = 0.8), fill = "white") +
  
  scale_fill_manual(values = c("WT" = "#4DBBD5", "E269R" = "#E64B35")) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) +
  
  labs(
    x = "Timepoint",
    y = "Expression (scaled)",
    fill = "Condition",
    title = "Distribution of heatmap gene expression\nWT vs E269R"
  )




#find overlapping genes----
library(dplyr)

get_DE <- function(df){
  df %>%
    filter(adj.P.Val < 0.05) %>%
    select(SYMBOL, logFC)
}

Mafb_all <- bind_rows(
  get_DE(Mafb_6hr_annotated) %>% mutate(time="6h"),
  get_DE(Mafb_10hr_annotated) %>% mutate(time="10h")
) %>%
  distinct(gene, .keep_all = TRUE)

E269R_Fos_all <- bind_rows(
  get_DE(E269R_6hr_annotated) %>% mutate(time="6h"),
  get_DE(E269R_10hr_annotated) %>% mutate(time="10h")
) %>%
  distinct(gene, .keep_all = TRUE)

overlap_microarray_E269R_fos_andMafb <- inner_join(
  Mafb_all %>% select(gene, Mafb_logFC = logFC),
  E269R_Fos_all %>% select(gene, E269R_Fos_logFC = logFC),
  by = "gene"
)


#Now seeing how Homodimer genes behave in hetrodimer condition----

get_DE <- function(df){
  df %>%
    tibble::rownames_to_column("gene") %>%
    select(SYMBOL, logFC, adj.P.Val)
}

Mafb6  <- get_DE(Mafb_6hr_annotated)
Mafb10 <- get_DE(Mafb_10hr_annotated)

E269R6  <- get_DE(E269R_Fos_annotated_6hr)
E269R10 <- get_DE(E269R_Fos_annotated_10hr)

library(dplyr)

Mafb_all <- bind_rows(
  Mafb6 %>% mutate(time="6h"),
  Mafb10 %>% mutate(time="10h")
) %>%
  filter(adj.P.Val < 0.05) %>%
  distinct(SYMBOL, .keep_all = TRUE)


E269R_all <- bind_rows(
  E269R6 %>% mutate(time="6h"),
  E269R10 %>% mutate(time="10h")
)

comparison_df <- comparison_df %>%
  mutate(
    E269R_direction = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Up",
      adj.P.Val < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "Not_significant"
    ),
    
    Mafb_direction = case_when(
      Mafb_logFC > 0 ~ "Up",
      Mafb_logFC < 0 ~ "Down",
      TRUE ~ "Neutral"
    ),
    
    status = case_when(
      adj.P.Val < 0.05 & Mafb_logFC * logFC > 0 ~ "Same",
      adj.P.Val < 0.05 & Mafb_logFC * logFC < 0 ~ "Reversed",
      TRUE ~ "Not_significant"
    )
  )
status_priority <- c("Reversed","Same","Not_significant")

comparison_df <- comparison_df %>%
  mutate(status = factor(status, levels=status_priority)) %>%
  group_by(SYMBOL) %>%
  slice_min(status) %>%
  ungroup()

comparison_df <- comparison_df %>%
  mutate(
    quadrant = case_when(
      Mafb_logFC > 0 & logFC > 0 ~ "Up_Up",
      Mafb_logFC < 0 & logFC < 0 ~ "Down_Down",
      Mafb_logFC > 0 & logFC < 0 ~ "Up_Down",
      Mafb_logFC < 0 & logFC > 0 ~ "Down_Up"
    )
  )

comparison_df <- as_tibble(comparison_df)

quad_counts <- comparison_df %>%
  dplyr::count(quadrant)
library(ggplot2)

ggplot(comparison_df, aes(Mafb_logFC, logFC)) +
  
  geom_point(aes(color = quadrant), alpha = 0.7, size = 1.5) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  
  annotate("text",
           x = max(comparison_df$Mafb_logFC)*0.5,
           y = max(comparison_df$logFC)*0.5,
           label = paste0("Up/Up: ", quad_counts$n[quad_counts$quadrant=="Up_Up"])) +
  
  annotate("text",
           x = min(comparison_df$Mafb_logFC)*0.5,
           y = min(comparison_df$logFC)*0.5,
           label = paste0("Down/Down: ", quad_counts$n[quad_counts$quadrant=="Down_Down"])) +
  
  annotate("text",
           x = max(comparison_df$Mafb_logFC)*0.5,
           y = min(comparison_df$logFC)*0.5,
           label = paste0("Up/Down: ", quad_counts$n[quad_counts$quadrant=="Up_Down"])) +
  
  annotate("text",
           x = min(comparison_df$Mafb_logFC)*0.5,
           y = max(comparison_df$logFC)*0.5,
           label = paste0("Down/Up: ", quad_counts$n[quad_counts$quadrant=="Down_Up"])) +
  
  theme_classic() +
  labs(
    x = "Mafb log2FC",
    y = "E269R + Fos logFC",
    title = "Transcriptional Direction Comparison"
  )


#Significant reversed genes

comparison_df <- as.data.frame(comparison_df)

bar_df <- comparison_df %>%
  dplyr::filter(Mafb_direction %in% c("Up","Down")) %>%
  dplyr::count(Mafb_direction, status)


library(ggplot2)

ggplot(bar_df, aes(x = Mafb_direction, y = n, fill = status)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  theme_classic() +
  labs(
    x = "MAFB regulation",
    y = "Number of genes",
    fill = "E269R effect"
  )








quad_genes <- split(comparison_df$SYMBOL, comparison_df$quadrant)
library(clusterProfiler)
library(org.Mm.eg.db)

compare_GO <- compareCluster(
  geneCluster = quad_genes,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)
dotplot(compare_GO, showCategory = 3)



