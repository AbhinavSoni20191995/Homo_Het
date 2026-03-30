#CHIP-Seq data analysis to resolve Homo vs Het data

 #Convert the CHIP files into bed files----
#Upload files for Group B and Group Y CHIP-Seq peaks
#Convert them to the bed files

GroupB_peaks<-as.data.frame(groupB_Flag.E269R_fos_vs_Flag.wtMafB)
GroupY_peaks<-as.data.frame(groupY_flag_wtMafB_vs_flag_E269R_Fos)

GroupB_peaks_short <- GroupB_peaks %>%
  select(chromosome_name, position, enrichment, neg_log10_qvalue, rank, gene_name)

GroupY_peaks_short <- GroupY_peaks %>%
  select(chromosome_name, position, enrichment, neg_log10_qvalue, rank, gene_name)

names(GroupB_peaks_short)[which(names(GroupB_peaks_short) == "chromosome_name")] <- "chr"
names(GroupB_peaks_short)[which(names(GroupB_peaks_short) == "position")] <- "start"

names(GroupY_peaks_short)[which(names(GroupY_peaks_short) == "chromosome_name")] <- "chr"
names(GroupY_peaks_short)[which(names(GroupY_peaks_short) == "position")] <- "start"

GroupB_peaks_short$end<-GroupB_peaks_short$start
GroupY_peaks_short$end<-GroupY_peaks_short$start

GroupB_peaks_short$Interval_rownumber <- paste("Interval", 1:nrow(GroupB_peaks_short), sep="_")
GroupY_peaks_short$Interval_rownumber <- paste("Interval", 1:nrow(GroupY_peaks_short), sep="_")

GroupB_peaks_short_gr<-Makegranges(GroupB_peaks_short)
GroupY_peaks_short_gr<-Makegranges(GroupY_peaks_short)

library(rtracklayer)
library(dplyr)
chain_file_mm9_mm10 <- "mm9ToMm10.over.chain" 
chain_file_mm9_mm10<- import.chain(chain_file_mm9_mm10)
GroupB_mm10 <- as.data.frame(liftOver(GroupB_peaks_short_gr, chain_file_mm9_mm10))
GroupB_mm10 <- GroupB_mm10[, !colnames(GroupB_mm10) %in% c("group", "group_name")] 
GroupB_mm10 <- GroupB_mm10 %>% 
  rename(Chr = chr, Start = start, Strand = strand, End =end)

GroupY_mm10 <- as.data.frame(liftOver(GroupY_peaks_short_gr, chain_file_mm9_mm10))
GroupY_mm10 <- GroupY_mm10[, !colnames(GroupY_mm10) %in% c("group", "group_name")] 
GroupY_mm10 <- GroupY_mm10 %>% 
  rename(Chr = chr, Start = start, Strand = strand, End =end)

GroupB_mm10<-makeGRangesFromDataFrame(GroupB_mm10,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqnames.field=c("Chr"),
                                      start.field="Start",
                                      end.field=c("End"),
                                      starts.in.df.are.0based=FALSE)
GroupY_mm10<-makeGRangesFromDataFrame(GroupY_mm10,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqnames.field=c("Chr"),
                                      start.field="Start",
                                      end.field=c("End"),
                                      starts.in.df.are.0based=FALSE)

export(GroupY_mm10, "GroupY_mm10.bed")
export(GroupB_mm10,"GroupB_mm10.bed")


#See which peaks are overlapping
GroupY_mm10_resize <- resize(GroupY_mm10, width=150, fix="start")
GroupB_mm10_resize <- resize(GroupB_mm10, width=150, fix="start")

 #these above files were loaded and then everything was conducted below(archived)----

GroupB_mm10<-makeGRangesFromDataFrame(GroupB_mm10_resize,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqnames.field=c("V1"),
                                      start.field="V2",
                                      end.field=c("V3"),
                                      starts.in.df.are.0based=FALSE)
GroupY_mm10<-makeGRangesFromDataFrame(GroupY_mm10_resize,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqnames.field=c("V1"),
                                      start.field="V2",
                                      end.field=c("V3"),
                                      starts.in.df.are.0based=FALSE)

Annotate_GroupY_mm10<-annotatePeak(GroupY_mm10, TxDb = TxDb,annoDb = "org.Mm.eg.db")
Annotate_GroupB_mm10<-annotatePeak(GroupB_mm10, TxDb = TxDb,annoDb = "org.Mm.eg.db")


Annotate_GroupY_mm10_df<-as.data.frame(Annotate_GroupY_mm10)
Annotate_GroupB_mm10_df<-as.data.frame(Annotate_GroupB_mm10)


GroupY_genes<-unique(Annotate_GroupY_mm10_df$SYMBOL)
GroupB_genes<-unique(Annotate_GroupB_mm10_df$SYMBOL)


#I think the above called peaks are not good enough, so I will use peaks from my analysis


#lets make consensus peak set and then gene sets for specific genes, checks their enrichment using deeptools


Combined_peaks<-makeGRangesFromDataFrame(combined,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=TRUE,
                                      seqnames.field=c("V1"),
                                      start.field="V2",
                                      end.field=c("V3"),
                                      starts.in.df.are.0based=FALSE)

Combined_peaks<-annotatePeak(Combined_peaks, TxDb = TxDb,annoDb = "org.Mm.eg.db")
Combined_peaks_df<-as.data.frame(Combined_peaks)




 #Take out peaks related to Cluster1 Mafb genes in PBS to show Mafb in the state of Mafb-Het is destabilised or unable to bind(archived)----

Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1 <-Combined_peaks_df %>% filter(SYMBOL %in% PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1)
df<-Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1
df <- df |>
  dplyr::rename(
    chr = seqnames,
  )

gr_Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1 <- makeGRangesFromDataFrame(
  df,
  keep.extra.columns = TRUE,
  ignore.strand = FALSE,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end"
)

export.bed(Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1,"Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1.bed")

#done in shell using deeptools

computeMatrix scale-regions \
-S FLAG_WTMAFB.bigWig FLAG_E269R_1.bigWig \
-R Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1.bed \
-b 0 -a 0 \
--skipZeros \
-o Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1.gz


plotProfile \
-m Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1.gz \
--perGroup \
--legendLocation upper-right \
--plotTitle "Enrichment over regions" \
-out enrichment_profile_Combined_peaks_PBS_e269r_vs_PBS_wt.df_DECREASE_cluster1.pdf

 #Take out peaks related to Cluster1 Mafb genes in PBS to show Mafb in the state of Mafb-Het is destabilised or unable to bind(archived)----

Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2 <-Combined_peaks_df %>% filter(SYMBOL %in% PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2)
df<-Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2
df <- df |>
  dplyr::rename(
    chr = seqnames,
  )

gr_Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2 <- makeGRangesFromDataFrame(
  df,
  keep.extra.columns = TRUE,
  ignore.strand = FALSE,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end"
)

export.bed(Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2,"Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2.bed")

#done in shell using deeptools

computeMatrix scale-regions \
-S FLAG_WTMAFB.bigWig FLAG_E269R_1.bigWig \
-R Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2.bed \
-b 0 -a 0 \
--skipZeros \
-o Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2.gz


plotProfile \
-m Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2.gz \
--perGroup \
--legendLocation upper-right \
--plotTitle "Enrichment over regions" \
-out enrichment_profile_Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2.pdf



 #Take out peaks related to Cluster1 Mafb genes to show general change in Mafb----

Combined_peaks_Clusters <-Combined_peaks_df %>% filter(SYMBOL %in% Mafb_annotated_significant_genes)
df<-Combined_peaks_Clusters 
df <- df |>
  dplyr::rename(
    chr = seqnames,
  )

gr_Combined_peaks_Clusters <- makeGRangesFromDataFrame(
  df,
  keep.extra.columns = TRUE,
  ignore.strand = FALSE,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end"
)

export.bed(Combined_peaks_Clusters,"Combined_peaks_Clusters.bed")

#done in shell using deeptools

computeMatrix scale-regions \
-S FLAG_WTMAFB.bigWig FLAG_E269R_1.bigWig FLAG_WTMAFB_PLUS_FOS.bigWig FLAG_E269R_PLUS_FOS.bigWig  \
-R Combined_peaks_Clusters.bed \
-b 0 -a 0 \
--skipZeros \
-o Combined_peaks_Clusters.bed.gz


plotProfile \
-m Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2.gz \
--perGroup \
--legendLocation upper-right \
--plotTitle "Enrichment over regions" \
-out enrichment_profile_Combined_peaks_PBS_e269r_vs_PBS_wt.df_INCREASE_cluster2.pdf


plotHeatmap -m Combined_peaks_Clusters.bed.gz \
-out Combined_peaks_Clusters.pdf \
--colorMap RdBu \
--kmeans 4  





#Lets build seprate homo het peaks, from peaks I derived through nf-core----
#conceptually we keep things same as it was done to derive groupY_flag_wtMafB_vs_flag_E269R+Fos_p30_e3_100nt (group y homo peaks) or groupB1_flag_E269R+Fos_vs_flag_wtMafB_p30_e3_100nt  (group b hetero)
# see "/Volumes/crtd_sieweke/DATA/Former_Lab_Members/DD/Hao_2022_06/HomoHet_Carole Berruyer_Hao Huang/HOMO_HETERODKO/PRESENTATIONS_HOMOHET/PRESENTATIONS_FORMERS/Peak_analysis_of MafB_subsets_14fev2012.pptx"

FLAG_WTMAFB_peaks_gr<-makeGRangesFromDataFrame(FLAG_WTMAFB_peaks,
                                         keep.extra.columns=TRUE,
                                         ignore.strand=TRUE,
                                         seqnames.field=c("V1"),
                                         start.field="V2",
                                         end.field=c("V3"),
                                         starts.in.df.are.0based=FALSE)
FLAG_E269R_1_peaks_gr<-makeGRangesFromDataFrame(FLAG_E269R_1_peaks,
                                               keep.extra.columns=TRUE,
                                               ignore.strand=TRUE,
                                               seqnames.field=c("V1"),
                                               start.field="V2",
                                               end.field=c("V3"),
                                               starts.in.df.are.0based=FALSE)
FLAG_E269R_PLUS_FOS_peaks_gr<-makeGRangesFromDataFrame(FLAG_E269R_PLUS_FOS_peaks,
                                               keep.extra.columns=TRUE,
                                               ignore.strand=TRUE,
                                               seqnames.field=c("V1"),
                                               start.field="V2",
                                               end.field=c("V3"),
                                               starts.in.df.are.0based=FALSE)
FLAG_WTMAFB_PLUS_FOS_peaks<-makeGRangesFromDataFrame(FLAG_WTMAFB_PLUS_FOS_peaks,
                                               keep.extra.columns=TRUE,
                                               ignore.strand=TRUE,
                                               seqnames.field=c("V1"),
                                               start.field="V2",
                                               end.field=c("V3"),
                                               starts.in.df.are.0based=FALSE)

Annotate_FLAG_WTMAFB_peaks_gr<-annotatePeak(FLAG_WTMAFB_peaks_gr, TxDb = TxDb,annoDb = "org.Mm.eg.db")
Annotate_FLAG_E269R_1_peaks_gr<-annotatePeak(FLAG_E269R_1_peaks_gr, TxDb = TxDb,annoDb = "org.Mm.eg.db")
Annotate_FLAG_WTMAFB_PLUS_FOS_peaks<-annotatePeak(FLAG_WTMAFB_PLUS_FOS_peaks, TxDb = TxDb,annoDb = "org.Mm.eg.db")
Annotate_FLAG_E269R_PLUS_FOS_peaks_gr<-annotatePeak(FLAG_E269R_PLUS_FOS_peaks_gr, TxDb = TxDb,annoDb = "org.Mm.eg.db")

library(GenomicRanges)
library(ChIPpeakAnno)
library(VennDiagram)

peaklist <- list(
  Set1 = FLAG_WTMAFB_peaks_gr,
  Set2 = FLAG_E269R_1_peaks_gr,
  Set3 = FLAG_E269R_PLUS_FOS_peaks_gr,
  Set4 = FLAG_WTMAFB_PLUS_FOS_peaks
)

ol <- findOverlapsOfPeaks(peaklist)
ol$venn_cnt
makeVennDiagram(ol)

#Take out peaks in Set1 and not Set3

library(stringr)

merged <- ol$mergedPeaks
unique <- ol$uniquePeaks

peak_comb <- sapply(merged$peakNames, paste, collapse = "///")
length(peak_comb) == length(merged)

set1_not_set3 <- merged[
  grepl("Set1", peak_comb) &
    !grepl("Set3", peak_comb)
]


set3_not_set1 <- merged[
  grepl("Set3", peak_comb) &
    !grepl("Set1", peak_comb)
]

library(ChIPpeakAnno)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

annot_s1 <- annotatePeakInBatch(
  set1_not_set3,
  AnnotationData = genes(txdb),
  output = "nearestLocation",
  bindingRegion = c(-5000, 3000)
)

annot_s1 <- addGeneIDs(
  annot_s1,
  orgAnn = "org.Mm.eg.db",
  IDs2Add = "symbol",
  feature_id_type = "entrez_id"
)

?annotatePeakInBatch

annot_s3 <- annotatePeakInBatch(
  set3_not_set1,
  AnnotationData = genes(txdb),
  output = "nearestLocation",
  bindingRegion = c(-5000, 3000)
)

annot_s3 <- addGeneIDs(
  annot_s3,
  orgAnn = "org.Mm.eg.db",
  IDs2Add = "symbol",
  feature_id_type = "entrez_id"
)

table(annot_s1$insideFeature)
table(annot_s3$insideFeature)

prop.table(table(annot_s1$insideFeature))
prop.table(table(annot_s3$insideFeature))

df1 <- prop.table(table(annot_s1$insideFeature))
df3 <- prop.table(table(annot_s3$insideFeature))

barplot(rbind(df1, df3),
        beside = TRUE,
        legend = c("Set1_not_Set3", "Set3_not_Set1"),
        las = 2)

summary(abs(annot_s1$distancetoFeature))
summary(abs(annot_s3$distancetoFeature))

mean(abs(annot_s1$distancetoFeature) <= 3000)
mean(abs(annot_s3$distancetoFeature) <= 3000)


promoters_gr <- promoters(genes(txdb), upstream=3000, downstream=3000)

length(subsetByOverlaps(set1_not_set3, promoters_gr))
length(subsetByOverlaps(set3_not_set1, promoters_gr))

pie1(table(annot_s1$insideFeature))

?toGRanges
annotation_data <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
binOverFeature(annot_s1, 
               featureSite = "FeatureStart",
               nbins = 20,
               annotationData = annotation_data,
               xlab = "peak distance from TSS (bp)", 
               ylab = "peak count", 
               main = "Distribution of aggregated peak numbers around TSS")

genomicElementDistribution(annot_s1, 
                           TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene)


genomicElementDistribution(annot_s3, 
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

annot_s1_df<-as.data.frame(annot_s1)
annot_s3_df<-as.data.frame(annot_s3)

export(annot_s1,"Homo_peaks_chipanno.bed")
export(annot_s3,"Hetro_peaks_chipanno.bed")

#Deeptools based plots in unix----

#Compute matrix for both basal and ttx condition, 500bp
computeMatrix reference-point --referencePoint center \
-R  Homo_peaks_chipanno.bed \
-S  FLAG_WTMAFB.bigWig \
FLAG_E269R_1.bigWig \
FLAG_WTMAFB_PLUS_FOS.bigWig \
FLAG_E269R_PLUS_FOS.bigWig \
-a 1000 -b 1000 \
-o omo_peaks_chipanno_computematrix.gz 

computeMatrix reference-point --referencePoint center \
-R  Hetro_peaks_chipanno.bed \
-S  FLAG_WTMAFB.bigWig \
FLAG_E269R_1.bigWig \
FLAG_WTMAFB_PLUS_FOS.bigWig \
FLAG_E269R_PLUS_FOS.bigWig \
-a 1000 -b 1000 \
-o Hetro_peaks_chipanno_computematrix.gz 

#See how the plotprofile and profile heatmap
plotProfile -m Homo_peaks_chipanno_computematrix.gz   \
-out Homo_peaks_chipanno_computematrix_profileplot.pdf \
--perGroup

plotHeatmap -m omo_peaks_chipanno_computematrix.gz  \
--colorMap RdBu_r \
-out Homo_peaks_chipanno_computematrix.pdf 

plotHeatmap -m Hetro_peaks_chipanno_computematrix.gz   \
--colorMap RdBu_r \
-out Hetro_peaks_chipanno_computematrix_heatmap.pdf 

#the peaks look good, so homo het evaluation is nice

#motif analysis
#we want to capture homo het like peaks

findMotifsGenome.pl Homo_peaks_chipanno.bed mm10 homer_s1/ -size 200 -mask 
findMotifsGenome.pl Hetro_peaks_chipanno.bed mm10 homer_s3/ -size 200 -mask

#also run with homer2
findMotifsGenome.pl Homo_peaks_chipanno.bed mm10 homer_s1/ -size 200 -mask -homer2
findMotifsGenome.pl Hetro_peaks_chipanno.bed mm10 homer_s3/ -size 100 -mask -homer2

#Make a comprehensive visualization 


#Comparative Gene ontology for homodimers and hetrodimer genes----

# Extract genes from annotation
genes_homo_s1 <- annot_s1_df %>%
  pull(symbol) %>%
  unique()

genes_hetro_s1 <- annot_s3_df %>%
  pull(symbol) %>%
  unique()

library(clusterProfiler)
library(org.Mm.eg.db)

genes_homo_s1_entrez <- bitr(
  genes_homo_s1,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

genes_hetro_s1_entrez <- bitr(
  genes_hetro_s1,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

gene_list <- list(
  Homodimer = genes_homo_s1_entrez$ENTREZID,
  Heterodimer = genes_hetro_s1_entrez$ENTREZID
)

ego_compare_homo_het <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP"
)

dotplot(ego_compare_homo_het) 

#Take out peaks which are common with Cluster 1 

library(dplyr)

annot_s1_df <- annot_s1_df %>%
  mutate(
    Mafb_cluster = case_when(
      symbol %in% Cluster1_genes ~ "Cluster1",
      symbol %in% Cluster2_genes ~ "Cluster2",
      TRUE ~ NA_character_
    )
  )

annot_s3_df <- annot_s3_df %>%
  mutate(
    Mafb_cluster = case_when(
      symbol %in% Cluster1_genes ~ "Cluster1",
      symbol %in% Cluster2_genes ~ "Cluster2",
      TRUE ~ NA_character_
    )
  )

#Number of homodimer associated with cluster 1 genes

n_homo_c1_peaks <- annot_s1_df %>%
  filter(Mafb_cluster == "Cluster1")

genes_homo_c1 <- n_homo_c1_peaks %>%
  filter(Mafb_cluster == "Cluster1") %>%
  pull(symbol) %>%
  unique()

n_homo_c2_peaks <- annot_s1_df %>%
  filter(Mafb_cluster == "Cluster2")

genes_homo_c2 <- n_homo_c2_peaks %>%
  filter(Mafb_cluster == "Cluster2") %>%
  pull(symbol) %>%
  unique()



#Motif plot----

#Heatmap for motif erichment according to size of coverage and motif 

knownResults$X..of.Target.Sequences.with.Motif<-gsub("%","",knownResults$X..of.Target.Sequences.with.Motif)
knownResults$X..of.Target.Sequences.with.Motif<-as.numeric(knownResults$X..of.Target.Sequences.with.Motif)
knownResults$Motif.Name<- sub("/.*","",knownResults$Motif.Name)

library(ggplot2)
library(ggrepel)

# Updated plot
plot <- ggplot(knownResults, aes(x = -1 * Log.P.value, y = X..of.Target.Sequences.with.Motif)) +
  # Points with black border and colored fill
  geom_point(
    aes(size = -1 * Log.P.value, color = X..of.Target.Sequences.with.Motif),
    shape = 21, stroke = 1  # Adds a black border while keeping the fill
  ) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +  # Custom color scale
  theme_minimal() +                                             # Clean theme
  labs(
    title = "Motif Coverage and Significance",
    x = "-log10(p-value)",
    y = "% Sequence Covered",
    color = "% Sequence Covered",
    size = "-log10(p-value)"
  ) +
  # Add horizontal and vertical threshold lines
  geom_vline(xintercept = 5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black") +
  # Add text labels for motifs above thresholds with connecting lines
  geom_label_repel(
    data = subset(knownResults, X..of.Target.Sequences.with.Motif > 10 & -1 * Log.P.value > 20),
    aes(label = Motif.Name),
    fontface = "bold",                # Make text bold
    size = 4,                         # Increase text size
    box.padding = 1.5,                # Padding around text box
    point.padding = 1.5,              # Padding around points
    segment.color = "black"           # Black connecting lines
  ) +
  # Theme with larger axis labels and visible ticks
  theme(
    panel.background = element_rect(fill = "white", color = "black"),  # White background with black border
    panel.grid.minor = element_blank(),                                # No minor grid lines
    plot.title = element_text(hjust = 0.5, face = "bold"),             # Centered bold title
    axis.text = element_text(color = "black", size = 12),              # Larger axis numbers
    axis.title = element_text(face = "bold", size = 14),               # Larger axis labels
    axis.ticks = element_line(color = "black", size = 0.5),            # Add ticks to axes
    legend.key = element_rect(fill = "white", color = NA)              # White legend background
  )

# Display the plot
ggplot(knownResults, aes(x = -1 * Log.P.value, y = X..of.Target.Sequences.with.Motif)) +
  # Points with black border and colored fill
  geom_point(
    aes(size = -1 * Log.P.value, color = X..of.Target.Sequences.with.Motif),
    shape = 21, stroke = 1  # Adds a black border while keeping the fill
  ) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +  # Custom color scale
  theme_minimal() +                                             # Clean theme
  labs(
    title = "Motif Coverage and Significance",
    x = "-log10(p-value)",
    y = "% Sequence Covered",
    color = "% Sequence Covered",
    size = "-log10(p-value)"
  ) +
  # Add horizontal and vertical threshold lines
  geom_vline(xintercept = 5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black") +
  
  # Theme with larger axis labels and visible ticks
  theme(
    panel.background = element_rect(fill = "white", color = "black"),  # White background with black border
    panel.grid.minor = element_blank(),                                # No minor grid lines
    plot.title = element_text(hjust = 0.5, face = "bold"),             # Centered bold title
    axis.text = element_text(color = "black", size = 12),              # Larger axis numbers
    axis.title = element_text(face = "bold", size = 14),               # Larger axis labels
    axis.ticks = element_line(color = "black", size = 0.5),            # Add ticks to axes
    legend.key = element_rect(fill = "white", color = NA)              # White legend background
  )


#Check signficant genes Homodimer target genes bound 
