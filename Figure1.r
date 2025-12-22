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

library(pheatmap)
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

#Run GSEA
library(tidyverse)
termtogene_homohet <- data.frame(gs_name = c(rep("Activated", length(Cluster1_genes)), rep("Repressed", length(Cluster2_genes))),
                                 gs_symbol = c(Cluster1_genes, Cluster2_genes))

termtogene_homohet<-na.omit(termtogene_homohet)

GSEA.results.CH.clusterprofiler<- GSEA(t_DKO_Mafb_DOX_10_GSEA_rnk, TERM2GENE=ms_gsea_H, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.CH.clusterprofiler_t<-as.tibble(GSEA.results.CH.clusterprofiler)

GSEA.results.C5.clusterprofiler<- GSEA(t_DKO_Mafb_DOX_10_GSEA_rnk, TERM2GENE=ms_gsea_C5, verbose=FALSE, eps=0, pvalueCutoff = 1,by="fgsea",minGSSize = 5)
GSEA.results.C5.clusterprofiler_t<-as.tibble(GSEA.results.C5.clusterprofiler)



edox_CH <- setReadable(GSEA.results.CH.clusterprofiler, 'org.Mm.eg.db', 'ENTREZID')
p<-cnetplot(edox_CH, node_label="category",showCategory = 543 )

?cnetplot


gseaplot2(GSEA.results.C8.CLUSTER7.clusterprofiler, 
          geneSetID = c(2:3), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          color= c("#ECCBAE", "#046C9A"))#can also turn off this title




#Gene ontology for Cluster1 and cluster 2 and all genes----
library(org.Mm.eg.db)
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

GO_allMafb_sginificant_genes <- enrichGO(gene = Mafb_siginificant_ENTREZID,
                        keyType = "ENTREZID",
                        OrgDb         = org.Mm.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)


GO_cluster1.df<-as.data.frame(GO_cluster1)
GO_allMafb_sginificant_genes.df<-as.data.frame(GO_allMafb_sginificant_genes)


#Isolate terms that contain certain keyword
library(dplyr)
library(stringr)
keywords <- c("coagulation","homeostasis","wound healing","innate immune response","T cell","cell cycle")

filtered_go <- GO_allMafb_sginificant_genes.df %>%
  filter(str_detect(tolower(Description), 
                    paste(tolower(keywords), collapse = "|"))) %>%
  arrange(p.adjust) %>%      # sort by significance
  slice_head(n = 25)         


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


#getting log fold change via differential value at 10hr
Mafb_10hr_annotated_summary <- Mafb_10hr_annotated  %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
geneList <- Mafb_10hr_annotated_summary$logFC
names(geneList) <-  Mafb_10hr_annotated_summary$SYMBOL

#run edox
edox <- setReadable(GO_cluster1, 'org.Mm.eg.db', 'ENTREZID')
edox_sub = edox[edox$ID %in% go_terms, asis=T]

edox_all <- setReadable(GO_allMafb_sginificant_genes, 'org.Mm.eg.db', 'ENTREZID')
edox_sub_2 = edox_all[edox_all$ID %in% filtered_go$ID, asis=T]

cnetplot(
  edox_sub,
  node_label = "all",
  showCategory = 12,
  color_category = 'firebrick',
) 

cnetplot(
  edox_sub_2,
  showCategory = 20,
  color_category = 'firebrick',
  foldChange=geneList,
  node_label="all"
) + scale_color_gradientn(
  colors = jdb_palette("solar_flare"),
  limits = c(-1.5, 1.5),       # force min / max of scale
  oob = scales::squish         # values outside will be squished to limits
)



GO_cluster1<- pairwise_termsim(GO_cluster1)
GO_all_mafb<- pairwise_termsim(GO_allMafb_sginificant_genes)

p<-emapplot(GO_cluster1, showCategory= 200)
library(BuenColors)
p <- emapplot(GO_all_mafb, showCategory=200, label_format = 0) +
  scale_color_gradientn(colors = jdb_palette("brewer_yes")) + theme_void()

#GSEA for Kang et al related immunosuppresive genes, here we use only DE associated genes

Kang_DE_genes <- c(
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




