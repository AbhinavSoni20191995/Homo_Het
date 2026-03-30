#Figure 1 (library)----
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

#Fig 1a----
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
                                                                         "wt_10h_rep_1","wt_10h_rep_2","wt_10h_rep_3")]



#Making heatmap
# --- Example: heatmap with clustering into groups ---


scaled_mat_gene_expression_Mafb_significant <- t(scale(t(gene_expression_Mafb_significant)))   # scale each gene (row) across samples

# Compute distance matrix (e.g., correlation-based)
dist_mat <- as.dist(1 - cor(t(scaled_mat_gene_expression_Mafb_significant)))  # correlation distance

# Hierarchical clustering
hc <- hclust(dist_mat, method = "ward.D2")  # Ward's method

# Cut tree into 4 clusters
clusters <- cutree(hc, k = 2)

# Assign cluster info
gene_clusters <- data.frame(Cluster = factor(clusters))
rownames(gene_clusters) <- rownames(scaled_mat_gene_expression_Mafb_significant)

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
conditions <- rep(c("NT", "6hr_Mafb", "10hr_Mafb"), each = 3)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(gene_expression_Mafb_significant)
col_ha <- HeatmapAnnotation(df = annotation_col)


annotation_row <- data.frame(Cluster = factor(gene_clusters$Cluster))
rownames(annotation_row) <- rownames(gene_expression_Mafb_significant)
row_ha <- rowAnnotation(df = annotation_row)


#Labels
set.seed(123)
highlight_genes<-c("Tim4", "C1qa", "C1qb", "C1qc", "Folr2", "Cebpb","Vsig4", "Pf4", "Cd5l", "Cbr2", "Tgfbi", "Rasgrp3", "Tmem37", "Tmem117", "Fcna", "Il4", "Il13", "Retnlg","Il6st","Lbp","Il27","Tnfsf9","Il15","Nfkb1","Irak2","Irf8","Ptgs2","Il6","Cxcl9","Tnfsf15", "Cxcl1","Arg2","Il1a", "Cxcl10", "Arg1", "Mmp8", "Cd209g","Il4","Mmp25","Mmp9","Retnlg","Akt1","Cebp","Chd1","Cirh1a","Cited2","Dppa3","Eed","Klf2","Klf4","Myc","Nfya","Nfyb","RhoJ","Stat3","Suz12","Ube2f")

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
  row_split = 2,                      # cut rows into 2 clusters
  show_row_names = FALSE,             # hide all labels
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

#Make summary plots----
#compute mean expression per cluster & sample
mean_tab <- df_long %>%
  group_by(Cluster, sample) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

# Step 2: join the mean back
df_long2 <- df_long %>%
  left_join(mean_tab, by = c("Cluster", "sample"))

# Step 3: define the palette
pal <- jdb_palette("solar_flare", n = 100, type = "continuous")

# Step 4: ggplot
ggplot(df_long2, aes(x = sample, y = expression, fill = mean_expr)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Cluster, scales = "free_y") +
  scale_fill_gradientn(
    colours = pal,
    name = "Mean expression"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "Condition (sample)",
    y = "Expression",
    title = "Expression by Condition per Cluster\ncoloured by mean expression"
  )

 #Make other kind of plots----
 # Compute mean and SD per cluster × sample
summary_tab <- df_long %>%
  group_by(Cluster, sample) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    sd_expr = sd(expression, na.rm = TRUE),
    .groups = "drop"
  )

# Define solar_flare palette
pal <- jdb_palette("solar_flare", n = length(unique(summary_tab$Cluster)), type = "continuous")

sample_order <- c(
  "wt_0h_rep_1", "wt_0h_rep_2", "wt_0h_rep_3",
  "wt_6h_rep_1", "wt_6h_rep_2", "wt_6h_rep_3",
  "wt_10h_rep_1", "wt_10h_rep_2", "wt_10h_rep_3"
)

# Set factor levels
summary_tab$sample <- factor(summary_tab$sample, levels = sample_order)

library(ggplot2)
library(BuenColors)

ggplot(summary_tab, aes(x = sample, y = mean_expr, group = Cluster, color = mean_expr, fill = mean_expr)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr), alpha = 0.3, color = NA) +
  facet_wrap(~ Cluster, scales = "free_y") +
  scale_color_gradientn(colours = pal, name = "Mean expression") +
  scale_fill_gradientn(colours = pal, name = "Mean expression") +
  theme_minimal(base_size = 14) +   # clean theme
  theme(
    panel.grid = element_blank(),   # remove grey grid lines
    panel.border = element_blank(), # remove panel border
    axis.line = element_line(color = "black"), # optional axes
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    x = "Condition (sample)",
    y = "Expression",
    title = "Mean expression per cluster with SD\ncoloured by mean expression"
  )


# Define custom colors for clusters
cluster_colors <- c(
  "1" = "#0C2C84",
  "2" = "#FFFFCC"
)

# Plot
ggplot(summary_tab, aes(x = sample, y = mean_expr, group = Cluster, color = mean_expr, fill = Cluster)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr), alpha = 0.3, color = NA) +
  scale_color_gradientn(colours = pal, name = "Mean expression") +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    x = "Condition (sample)",
    y = "Expression",
    title = "Mean expression per cluster with SD"
  )



#Gene ontology for Cluster1 and cluster 2 and all genes----
Cluster1_EntrezID <- mapIds(org.Mm.eg.db, Cluster1_genes, "ENTREZID", "SYMBOL")
Cluster2_EntrezID <- mapIds(org.Mm.eg.db, Cluster2_genes, "ENTREZID", "SYMBOL")
Mafb_siginificant_ENTREZID<-mapIds(org.Mm.eg.db,Maf_significant_genes, "ENTREZID", "SYMBOL")
GO_cluster1 <- enrichGO(gene = Cluster1_EntrezID,
                                                     keyType = "ENTREZID",
                                                     OrgDb         = org.Mm.eg.db,
                                                     ont           = "BP",
                                                     pAdjustMethod = "BH",
                                                     pvalueCutoff  = 0.05,
                                                     qvalueCutoff  = 0.05,
                                                     readable      = TRUE)

GO_cluster2 <- enrichGO(gene = Cluster2_EntrezID,
                        keyType = "ENTREZID",
                        OrgDb         = org.Mm.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)


GO_allMafb_sginificant_genes <- enrichGO(gene = Mafb_siginificant_ENTREZID,
                        keyType = "ENTREZID",
                        OrgDb         = org.Mm.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)


GO_cluster1.df<-as.data.frame(GO_cluster1)
write.csv(GO_cluster1.df, file="GO_cluster1_Figure1.csv")
GO_cluster2.df<-as.data.frame(GO_cluster2)
write.csv(GO_cluster2.df, file="GO_cluster2_Figure1.csv")
GO_allMafb_sginificant_genes.df<-as.data.frame(GO_allMafb_sginificant_genes)
write.csv(GO_allMafb_sginificant_genes.df, file="GO_allMafb_sginificant_genes.csv")


edox_cluster2 <- setReadable(GO_cluster2, 'org.Mm.eg.db', 'ENTREZID')
heatplot(edox_cluster2, foldChange=geneList, showCategory=5)

dotplot(GO_cluster2, showCategory=10)
dotplot(GO_cluster1, showCategory=10)


#Isolate terms that contain certain keyword (to do)

keywords <- c("coagulation","homeostasis","wound healing","cell cycle","apoptotic cell clearance","receptor-mediated endocytosis")

filtered_go_cluster1 <- GO_cluster1.df %>%
  filter(str_detect(tolower(Description), 
                    paste(tolower(keywords), collapse = "|"))) %>%
  arrange(p.adjust) %>%      # sort by significance
  slice_head(n = 10)     


 #Run emaplot----
GO_cluster1<- pairwise_termsim(GO_cluster1)
emapplot(GO_cluster1, label_format = 0) +
  scale_color_gradientn(colors = jdb_palette("brewer_yes")) + theme_void()



 #run edox----
edox <- setReadable(GO_cluster1, 'org.Mm.eg.db', 'ENTREZID')
edox_sub = edox[edox$ID %in% filtered_go_cluster1, asis=T]

#getting log fold change via differential value at 10hr
Mafb_10hr_annotated_summary <- Mafb_10hr_annotated  %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
geneList <- Mafb_10hr_annotated_summary$logFC
names(geneList) <-  Mafb_10hr_annotated_summary$SYMBOL

cnetplot(
  edox_sub,
  color_category = 'firebrick',
  foldChange=geneList,
  node_label="all"
) + scale_color_gradientn(
  colors = jdb_palette("solar_flare"),
  limits = c(-1.5, 1.5),       # force min / max of scale
  oob = scales::squish         # values outside will be squished to limits
)


GO_tagged <- GO_allMafb_sginificant_genes.df %>%
  mutate(keyword_group = map_chr(
    Description,
    ~ keywords[str_detect(tolower(.x), tolower(keywords))] %>% 
      first() %>% 
      coalesce(NA_character_)
  ))

GO_tagged <- GO_tagged %>% filter(!is.na(keyword_group))
filtered_go <- GO_tagged %>%
  group_by(keyword_group) %>%
  arrange(p.adjust) %>%     # ascending → most significant first
  slice_head(n = 5) %>%
  ungroup()




edox_all <- setReadable(GO_allMafb_sginificant_genes, 'org.Mm.eg.db', 'ENTREZID')
edox_sub_2 = edox_all[edox_all$ID %in% filtered_go$ID, asis=T]




cnetplot(
  edox_sub,
  node_label = "all",
  showCategory = 12,
  color_category = 'firebrick',
) 





GO_cluster1<- pairwise_termsim(GO_cluster1)
GO_all_mafb<- pairwise_termsim(GO_allMafb_sginificant_genes)

p<-emapplot(GO_cluster1, showCategory= 200)
library(BuenColors)
p <- emapplot(GO_all_mafb, showCategory=200, label_format = 0) +
  scale_color_gradientn(colors = jdb_palette("brewer_yes")) + theme_void()

#Take out the genes from each cluster -> perform GSEA and then plots (To be done)---- 

Cluster1_genes <- rownames(scaled_df[scaled_df$Cluster == '1', ])
Cluster2_genes <- rownames(scaled_df[scaled_df$Cluster == '2', ])

#using GSEA in cluster profiler
# are also retrieved as tibbles
msigdbr_species()
ms_gsea <- msigdbr(species = "Mus musculus") #gets all collections/signatures with human gene IDs
#take a look at the categories and subcategories of signatures available to you
ms_gsea %>% 
  dplyr::distinct(gs_cat, gs_subcat) %>% 
  dplyr::arrange(gs_cat, gs_subcat)
# choose a specific msigdb collection/subcollection
# since msigdbr returns a tibble, we'll use dplyr to do a bit of wrangling
ms_gsea_H <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                     category = "H") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 

ms_gsea_C5 <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                      category = "C5") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 

ms_gsea_C2 <- msigdbr(species = "Mus musculus", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 


# Pull out just the columns corresponding to gene symbols and LogFC
t_DKO_Mafb_DOX_10_GSEA <- as.data.frame(cbind(Mafb_10hr_annotated$SYMBOL, Mafb_10hr_annotated$logFC)) 
t_DKO_Mafb_DOX_10_GSEA <- t_DKO_Mafb_DOX_10_GSEA[!duplicated(t_DKO_Mafb_DOX_10_GSEA$V1), ]
t_DKO_Mafb_DOX_10_GSEA <- t_DKO_Mafb_DOX_10_GSEA[!duplicated(t_DKO_Mafb_DOX_10_GSEA$V2), ]
t_DKO_Mafb_DOX_10_GSEA <- na.omit(t_DKO_Mafb_DOX_10_GSEA)

t_DKO_Mafb_DOX_10_GSEA_rnk <- as.numeric(t_DKO_Mafb_DOX_10_GSEA$V2)
names(t_DKO_Mafb_DOX_10_GSEA_rnk)<-as.character(t_DKO_Mafb_DOX_10_GSEA$V1)
t_DKO_Mafb_DOX_10_GSEA_rnk <- sort(t_DKO_Mafb_DOX_10_GSEA_rnk, decreasing = TRUE)


t_DKO_Mafb_DOX_6_GSEA <- as.data.frame(cbind(Mafb_6hr_annotated$SYMBOL, Mafb_6hr_annotated$logFC)) 
t_DKO_Mafb_DOX_6_GSEA <- t_DKO_Mafb_DOX_6_GSEA[!duplicated(t_DKO_Mafb_DOX_6_GSEA$V1), ]
t_DKO_Mafb_DOX_6_GSEA <- t_DKO_Mafb_DOX_10_GSEA[!duplicated(t_DKO_Mafb_DOX_6_GSEA$V2), ]
t_DKO_Mafb_DOX_6_GSEA <- na.omit(t_DKO_Mafb_DOX_6_GSEA)

t_DKO_Mafb_DOX_6_GSEA_rnk <- as.numeric(t_DKO_Mafb_DOX_6_GSEA$V2)
names(t_DKO_Mafb_DOX_6_GSEA_rnk)<-as.character(t_DKO_Mafb_DOX_6_GSEA$V1)
t_DKO_Mafb_DOX_6_GSEA_rnk <- sort(t_DKO_Mafb_DOX_6_GSEA_rnk, decreasing = TRUE)
#Run GSEA
library(tidyverse)
termtogene_homohet <- data.frame(gs_name = c(rep("Activated", length(Cluster1_genes)), rep("Repressed", length(Cluster2_genes))),
                                 gs_symbol = c(Cluster1_genes, Cluster2_genes))

termtogene_homohet<-na.omit(termtogene_homohet)

GSEA.results.CH.clusterprofiler<- GSEA(t_DKO_Mafb_DOX_10_GSEA_rnk, TERM2GENE=ms_gsea_H, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.CH.clusterprofiler_t<-as.tibble(GSEA.results.CH.clusterprofiler)

GSEA.results.C5.clusterprofiler<- GSEA(t_DKO_Mafb_DOX_10_GSEA_rnk, TERM2GENE=ms_gsea_C5, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.C5.clusterprofiler_t<-as.tibble(GSEA.results.C5.clusterprofiler)

dotplot(GSEA.results.CH.clusterprofiler,showCategory=5)


#GSEA for Kang et al related immunosuppresive genes, here we use only DE associated genes----

Kang_DE_genes <- as.data.framec(
  "Sepp1","Cd28","Klhl13","Mro","Tmem37","Ror2","Gpr34","Znf704","Crhbp","Thbs1",
  "Armc9","Deptor","Fblim1","Pdk4","Lgmn","Marco","Cnrip1","Cables1","Mtss1",
  "Ephb3","Scamp5","Ptgfrn","Pid1","Kiaa1147","Kiaa1671",
  "Itga9","Rnase4","Bmp2","Ltbp2","Kif26b","Pltp","Psca","Pmp22","Itsn1","Celsr1",
  "Rhob","Eepd1","Adam11","Pex11g","Lpar5","Hrh1","Slmo1","Ankh","Pdgfc",
  "Pik3ip1","Cryl1","Fuca1","Mmd","Asrgl1","Zc3h12d",
  "Klf2","Ldlrad4","Mertk","Nisch","Foxp1","Tesk2","Ap2a2","Dcbld1","Mgat4a",
  "Maml3","Maf","Tmem173","Ppm1l","Ssbp2","Tanc2","Sgpl1","Mgat5","B3gnt2",
  "Crem","Cpq","Fam198b","Apmap","Ntan1","Trps1","Sms",
  "Raph1","Tnfsf8","Zbtb4","Stard13","Xylt1","Thada","C20orf194","Socs6",
  "Zfp36l1","Def6","Mafb","Map3k1","Ly86","Lpar6","Swap70","Ap1b1","Msr1",
  "Adam9","Zadh2","Alox5","C22orf29","Kiaa0195","Dab2","Adcy9"
)


termtogene_kang_DE <- data.frame(
  gs_name   = rep("INTERFERRON_Disassembled_genes", length(Kang_DE_genes)),
  gs_symbol = Kang_DE_genes
)

termtogene_kang_DE<-na.omit(termtogene_kang_DE)

GSEA.results.Kang.DE.clusterprofiler_10hr<- GSEA(t_DKO_Mafb_DOX_10_GSEA_rnk, TERM2GENE=termtogene_kang_DE, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.Kang.DE.clusterprofiler_10hr<-as.tibble(GSEA.results.Kang.DE.clusterprofiler)
gseaplot2(GSEA.results.Kang.DE.clusterprofiler_10hr, geneSetID = 1,color= "#99000D")

GSEA.results.Kang.DE.clusterprofiler_6hr<- GSEA(t_DKO_Mafb_DOX_6_GSEA_rnk, TERM2GENE=termtogene_kang_DE, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.Kang.DE.clusterprofiler_6hr<-as.tibble(GSEA.results.Kang.DE.clusterprofiler_6hr)
gseaplot2(GSEA.results.Kang.DE.clusterprofiler_6hr, geneSetID = 1,color= "#FC9272")

leading_KANG_DE_10hr<-strsplit(GSEA.results.Kang.DE.clusterprofiler_10hr@result[["core_enrichment"]],"/")
leading_KANG_DE_10hr<-leading_KANG_DE_10hr[[1]]

leading_KANG_DE_6hr<-strsplit(GSEA.results.Kang.DE.clusterprofiler_6hr@result[["core_enrichment"]],"/")
leading_KANG_DE_6hr<-leading_KANG_DE_6hr[[1]]

leading_KANG_DE <- intersect(leading_KANG_DE_6hr,leading_KANG_DE_10hr)

#take out values from microarray dataset
library(dplyr)
gene_expression_df <- expr_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
# --- 4. Make gene names rownames ---
gene_expression_df <- as.data.frame(gene_expression_df)
rownames(gene_expression_df) <- gene_expression_df$GeneSymbol
gene_expression_df$GeneSymbol <- NULL

gene_expression_Mafb_Kang_de_leading_edge <- gene_expression_df[rownames(gene_expression_df) %in% leading_KANG_DE, ]
gene_expression_Mafb_Kang_de_leading_edge  <- gene_expression_Mafb_Kang_de_leading_edge [ ,c("wt_0h_rep_1","wt_0h_rep_2","wt_0h_rep_3",
                                                                         "wt_6h_rep_1","wt_6h_rep_2","wt_6h_rep_3",
                                                                         "wt_10h_rep_1","wt_10h_rep_2","wt_10h_rep_3")]
#intersection with cluster 1 genes
Kang_DE_c1 <- intersect(Cluster1_genes,Kang_DE_genes)

#Make a corresponding heatplot
gene_expression_Mafb_Kang_de_leading_edge<-t(scale(t(gene_expression_Mafb_Kang_de_leading_edge)))

#Annotation row
conditions <- rep(c("NT", "6hr_Mafb", "10hr_Mafb"), each = 3)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(gene_expression_Mafb_Kang_de_leading_edge)
col_ha_kang <- HeatmapAnnotation(df = annotation_col)


cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)


ht <- Heatmap(
  gene_expression_Mafb_Kang_de_leading_edge,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,             # hide all labels
  show_column_names = FALSE,
  col = col_fun,
  top_annotation = col_ha_kang
)


#Perform same for Self renewal gene set and TAM or immunosuppresive signature , need to choose according to biology (metabolic fitness? biosynthetic activity?)----

Self_renewal_genes<-c("Akt1","Cebpz","Chd1","Cited2","Dppa3","Eed","Klf2","Klf4","Myc","Nfya","Nfyb","Suz12","Ube2f")

termtogene_self_renewal <- data.frame(
  gs_name   = rep("Self_renewal", length(Self_renewal_genes)),
  gs_symbol = Self_renewal_genes
)

termtogene_self_renewal<-na.omit(termtogene_self_renewal)

GSEA.results.self_renewal.clusterprofiler_10hr<- GSEA(t_DKO_Mafb_DOX_10_GSEA_rnk, TERM2GENE=termtogene_self_renewal, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.self_renewal.clusterprofiler_10hr<-as.tibble(GSEA.results.self_renewal.clusterprofiler_10hr)
gseaplot2(GSEA.results.self_renewal.clusterprofiler_10hr, geneSetID = 1,color= "#99000D")

GSEA.results.self_renewal.clusterprofiler_6hr<- GSEA(t_DKO_Mafb_DOX_6_GSEA_rnk, TERM2GENE=termtogene_self_renewal, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.self_renewal.clusterprofiler_6hr<-as.tibble(GSEA.results.Kang.DE.clusterprofiler_6hr)
gseaplot2(GSEA.results.self_renewal.clusterprofiler_6hr, geneSetID = 1,color= "#FC9272")

leading_self_renewal_10hr<-strsplit(GSEA.results.self_renewal.clusterprofiler_10hr@result[["core_enrichment"]],"/")
leading_self_renewal_10hr<-leading_self_renewal_10hr[[1]]

leading_KANG_DE <- intersect(leading_KANG_DE_6hr,leading_KANG_DE_10hr)

#take out values from microarray dataset
library(dplyr)
gene_expression_df <- expr_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
# --- 4. Make gene names rownames ---
gene_expression_df <- as.data.frame(gene_expression_df)
rownames(gene_expression_df) <- gene_expression_df$GeneSymbol
gene_expression_df$GeneSymbol <- NULL

gene_expression_self_renewal_de_leading_edge <- gene_expression_df[rownames(gene_expression_df) %in% leading_self_renewal_10hr, ]
gene_expression_self_renewal_de_leading_edge  <- gene_expression_self_renewal_de_leading_edge[ ,c("wt_0h_rep_1","wt_0h_rep_2","wt_0h_rep_3",
                                                                                             "wt_6h_rep_1","wt_6h_rep_2","wt_6h_rep_3",
                                                                                             "wt_10h_rep_1","wt_10h_rep_2","wt_10h_rep_3")]
#intersection with cluster 1 genes
Self_renewal_c1 <- intersect(Cluster2_genes,Self_renewal_genes)

#Make a corresponding heatplot
gene_expression_self_renewal_de_leading_edge<-t(scale(t(gene_expression_self_renewal_de_leading_edge)))

#Annotation row
conditions <- rep(c("NT", "6hr_Mafb", "10hr_Mafb"), each = 3)
annotation_col <- data.frame(Condition = conditions)
rownames(annotation_col) <- colnames(gene_expression_self_renewal_de_leading_edge)
col_ha_kang <- HeatmapAnnotation(df = annotation_col)


cols <- jdb_palette("solar_flare")  # or jdb_palettes("solar_flare")
col_fun <- colorRamp2(
  seq(-2, 2, length.out = length(cols)), 
  cols
)


ht <- Heatmap(
  gene_expression_self_renewal_de_leading_edge,
  name = "Expr",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,             # hide all labels
  show_column_names = FALSE,
  col = col_fun,
  top_annotation = col_ha_kang
)

draw(ht)

#Genes Mafb dependent in LPMs----
library(clusterProfiler)
library(dplyr)
library(enrichplot)
library(ComplexHeatmap)
library(circlize)

#----------------------------
# 1. Define gene sets
#----------------------------
#LPM_Mafb_genes <- mmc4_Mafb_2026_Immunity$...1
#KC_Mafb_genes  <- mmc4_Mafb_2026_Immunity$...1
#IM_Mafb_genes  <- mmc4_Mafb_2026_Immunity$...1
#MG_Mafb_genes  <- mmc4_Mafb_2026_Immunity$...1

gene_sets <- list(
  LPM_Mafb = na.omit(LPM_Mafb_genes),
  KC_Mafb  = na.omit(KC_Mafb_genes),
  IM_Mafb  = na.omit(IM_Mafb_genes),
  MG_Mafb  = na.omit(MG_Mafb_genes)
)

# convert to TERM2GENE
termtogene <- stack(gene_sets)
colnames(termtogene) <- c("gs_symbol","gs_name")
termtogene <- termtogene[,c("gs_name","gs_symbol")]

#----------------------------
# 2. Run GSEA (6hr and 10hr)
#----------------------------
gsea_6hr <- GSEA(t_DKO_Mafb_DOX_6_GSEA_rnk,
                 TERM2GENE = termtogene,
                 verbose = FALSE,
                 eps = 0,
                 pvalueCutoff = 1,
                 by="fgsea",
                 minGSSize = 5)

gsea_10hr <- GSEA(t_DKO_Mafb_DOX_10_GSEA_rnk,
                  TERM2GENE = termtogene,
                  verbose = FALSE,
                  eps = 0,
                  pvalueCutoff = 1,
                  by="fgsea",

#----------------------------
# 3. Plot all GSEA curves
#----------------------------

p1 <- gseaplot2(gsea_6hr,
                geneSetID = c("LPM_Mafb","KC_Mafb","IM_Mafb","MG_Mafb"),
                color=c("#99000D","#FC9272","#6BAED6","#08519C"),
                title="6hr Mafb")

p2 <- gseaplot2(gsea_10hr,
                geneSetID = c("LPM_Mafb","KC_Mafb","IM_Mafb","MG_Mafb"),
                color=c("#99000D","#FC9272","#6BAED6","#08519C"),
                title="10hr Mafb")

p1
p2
#----------------------------
# 4. Extract leading edge genes
#----------------------------
# extract leading edge genes
res6 <- gsea_6hr@result
leading_edges <- lapply(res6$core_enrichment, function(x) strsplit(x,"/")[[1]])
names(leading_edges) <- res6$ID

# union of all leading edge genes
leading_union <- Reduce(union, leading_edges)

#----------------------------
# 5. Prepare expression matrix
#----------------------------
gene_expression_df <- expr_df %>%
  group_by(GeneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

gene_expression_df <- as.data.frame(gene_expression_df)
rownames(gene_expression_df) <- gene_expression_df$GeneSymbol
gene_expression_df$GeneSymbol <- NULL

expr_leading <- gene_expression_df[rownames(gene_expression_df) %in% leading_union,]

expr_leading <- expr_leading[,c("wt_0h_rep_1","wt_0h_rep_2","wt_0h_rep_3",
                                "wt_6h_rep_1","wt_6h_rep_2","wt_6h_rep_3",
                                "wt_10h_rep_1","wt_10h_rep_2","wt_10h_rep_3")]

expr_leading <- t(scale(t(expr_leading)))

#----------------------------
# 6. Heatmap
#----------------------------
conditions <- rep(c("NT","6hr_Mafb","10hr_Mafb"), each=3)
annotation_col <- data.frame(Condition=conditions)
rownames(annotation_col) <- colnames(expr_leading)

col_ha <- HeatmapAnnotation(df=annotation_col)

cols <- jdb_palette("solar_flare")
col_fun <- colorRamp2(seq(-2,2,length.out=length(cols)),cols)

ht <- Heatmap(expr_leading,
              name="Expr",
              cluster_rows=TRUE,
              cluster_columns=FALSE,
              show_row_names=TRUE,
              show_column_names=FALSE,
              col=col_fun,
              top_annotation=col_ha)

draw(ht)

#Take out peaks for Mafb for cluster 1 genes----

Annotate_FLAG_WTMAFB_peaks_gr_df<-as.data.frame(Annotate_FLAG_WTMAFB_peaks_gr)

Annotate_FLAG_WTMAFB_peaks_gr_df_cluster1<-Annotate_FLAG_WTMAFB_peaks_gr_df %>% filter(SYMBOL %in% Cluster1_genes & V9>=10)
Annotate_FLAG_WTMAFB_peaks_gr_df_cluster2<-Annotate_FLAG_WTMAFB_peaks_gr_df %>% filter(SYMBOL %in% Cluster2_genes)

Makegranges<-function(df)
{
  makeGRangesFromDataFrame(df,
                           keep.extra.columns=FALSE,
                           ignore.strand=TRUE,
                           seqnames.field=c("seqnames", "seqname",
                                            "chromosome", "chrom",
                                            "chr", "chromosome_name",
                                            "seqid"),
                           start.field="start",
                           end.field=c("end", "stop"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE)
}

Annotate_FLAG_WTMAFB_peaks_cluster1_gr<-Makegranges(Annotate_FLAG_WTMAFB_peaks_gr_df_cluster1)
Annotate_FLAG_WTMAFB_peaks_cluster2_gr<-Makegranges(Annotate_FLAG_WTMAFB_peaks_gr_df_cluster2)

export(Annotate_FLAG_WTMAFB_peaks_cluster1_gr,"Annotate_FLAG_WTMAFB_peaks_cluster1_gr.bed")
export(Annotate_FLAG_WTMAFB_peaks_cluster2_gr,"Annotate_FLAG_WTMAFB_peaks_cluster2_gr.bed")

#Deeptools based plots in unix----

#Compute matrix for both basal and ttx condition, 500bp
computeMatrix reference-point \
--referencePoint center \
-b 2000 \
-a 2000 \
-R Annotate_FLAG_WTMAFB_peaks_cluster1_gr.bed \
-S H3K27ac_DKO.bigWig \
DKO_Mafb_OE_H3K27ac.bigWig \
H3K27ac_WTBMmacs.bigWig \
AlvMac_H3K27ac.bigWig \
-o Annotate_FLAG_WTMAFB_peaks_cluster1_gr_computematrix.gz

computeMatrix reference-point \
--referencePoint center \
-b 2000 \
-a 2000 \
-R Annotate_FLAG_WTMAFB_peaks_cluster2_gr.bed \
-S H3K27ac_DKO.bigWig \
DKO_Mafb_OE_H3K27ac.bigWig \
H3K27ac_WTBMmacs.bigWig \
AlvMac_H3K27ac.bigWig \
-o Annotate_FLAG_WTMAFB_peaks_cluster2_gr_computematrix.gz

#P300

computeMatrix reference-point \
--referencePoint center \
-b 2000 \
-a 2000 \
-R Annotate_FLAG_WTMAFB_peaks_cluster1_gr.bed \
-S P300_DKO.bigWig \
DKO_Mafb_OE_P300.bigWig \
P300_WTBMmacs.bigWig \
AlvMac_P300.bigWig \
-o Annotate_FLAG_WTMAFB_peaks_cluster1_gr_computematrix_P300.gz

computeMatrix reference-point \
--referencePoint center \
-b 2000 \
-a 2000 \
-R Annotate_FLAG_WTMAFB_peaks_cluster2_gr.bed \
-S P300_DKO.bigWig \
DKO_Mafb_OE_P300.bigWig \
P300_WTBMmacs.bigWig \
AlvMac_P300.bigWig \
-o Annotate_FLAG_WTMAFB_peaks_cluster2_gr_computematrix_P300.gz




#See how the plotprofile and profile heatmap
plotProfile -m Annotate_FLAG_WTMAFB_peaks_cluster1_gr_computematrix.gz   \
-out Annotate_FLAG_WTMAFB_peaks_cluster1_gr_computematrix.pdf \
--perGroup

plotProfile -m Annotate_FLAG_WTMAFB_peaks_cluster2_gr_computematrix.gz   \
-out Annotate_FLAG_WTMAFB_peaks_cluster2_gr_computematrix.pdf \
--perGroup

plotHeatmap -m Annotate_FLAG_WTMAFB_peaks_cluster1_gr_computematrix.gz  \
--colorMap Blues \
--hclust 4  \
-out Annotate_FLAG_WTMAFB_peaks_cluster1_gr_computematrix_heatmap_4cluster.pdf 


#P300
plotProfile -m Annotate_FLAG_WTMAFB_peaks_cluster1_gr_computematrix_P300.gz   \
-out Annotate_FLAG_WTMAFB_peaks_cluster1_gr_computematrix_P300.pdf \
--perGroup

#Make a comprehensive visualization 


