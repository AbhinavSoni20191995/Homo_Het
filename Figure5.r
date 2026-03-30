library(ProjecTILs)
library(BuenColors)
ref <- load.reference.map(). #this used for the moment

#Isolate Tcell from integrated dataset----
library(Seurat)
library(dplyr)
library(BuenColors)
Tcell_integrated <- subset(Integrated_sc, subset = SingleR.labels %in% c("T cells","Tgd","NKT"))

refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd")
DimPlot(ref, label = T, cols = refCols)

Tcell_integrated@meta.data <- Tcell_integrated@meta.data %>%
  mutate(
    samplid_combined = case_when(
      grepl("^LPS_24H", sampleid) ~ "LPS_24H",
      grepl("^LPS_4H",  sampleid) ~ "LPS_4H",
      grepl("^LPS_1H",  sampleid) ~ "LPS_1H",
      grepl("^NT_0H",   sampleid) ~ "NT_0H",
      TRUE ~ NA_character_
    )
  )

#run query reference 
Tcell_integrated <- Run.ProjecTILs(Tcell_integrated, ref = ref)
Tcell_integrated_0 <- subset(Tcell_integrated, subset = samplid_combined %in% c("NT_0H"))
Tcell_integrated_1 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_1H"))
Tcell_integrated_4 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_4H"))
Tcell_integrated_24 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_24H"))

#save query reference plot
p1<-plot.projection(ref, Tcell_integrated_0, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p2<-plot.projection(ref, Tcell_integrated_1, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p3<-plot.projection(ref, Tcell_integrated_4, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p4<-plot.projection(ref, Tcell_integrated_24, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p1+p2+p3+p4

#fraction plot
p1<-plot.statepred.composition(ref, Tcell_integrated_0, metric = "Percent",cols = jdb_palette("lawhoops"))
p2<-plot.statepred.composition(ref, Tcell_integrated_1, metric = "Percent",cols = jdb_palette("lawhoops"))
p3<-plot.statepred.composition(ref, Tcell_integrated_4, metric = "Percent",cols = jdb_palette("lawhoops"))
p4<-plot.statepred.composition(ref, Tcell_integrated_24, metric = "Percent",cols = jdb_palette("lawhoops"))

p1+p2+p3+p4


#Another way to plot it

library(dplyr)
library(ggplot2)

meta <- Tcell_integrated@meta.data

# count + fraction
df <- meta %>%
  group_by(samplid_combined, functional.cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(samplid_combined) %>%
  mutate(frac = n / sum(n))

df$samplid_combined <- factor(
  df$samplid_combined,
  levels = c("NT_0H","LPS_1H","LPS_4H","LPS_24H")
)

ggplot(df, aes(x = samplid_combined, y = frac, 
               group = functional.cluster, 
               color = functional.cluster)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  ylab("Fraction of cells") +
  xlab("Condition") +
  scale_color_manual(values = jdb_palette("lawhoops")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  facet_wrap(~functional.cluster)

#Compare signatures
library(SignatuR)
data(SignatuR)

programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)

library(UCell)
ref <- AddModuleScore_UCell(ref, features = programs, assay = "RNA", name = NULL)
Tcell_integrated_0 <- AddModuleScore_UCell(Tcell_integrated_0, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_1 <- AddModuleScore_UCell(Tcell_integrated_1, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_4 <- AddModuleScore_UCell(Tcell_integrated_4, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_24 <- AddModuleScore_UCell(Tcell_integrated_24, features = programs, assay = "RNA",name = NULL)

plot.states.radar(ref, query = Tcell_integrated_0, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_1, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_4, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_0, meta4radar = names(programs))


#visualize on the original umap 
Tcell_integrated <- ProjecTILs.classifier(query = Tcell_integrated, ref = ref)
Tcell_integrated$samplid_combined <- factor(
  Tcell_integrated$samplid_combined,
  levels = c("NT_0H", "LPS_1H", "LPS_4H", "LPS_24H")
)
#do it again with cca reduction
DimPlot(Tcell_integrated, group.by = "functional.cluster", split.by="samplid_combined",cols = jdb_palette("lawhoops"))

#save RDS obejct for future
Tcell_integrated_PTIL<-Tcell_integrated
saveRDS(Tcell_integrated_PTIL,file="Tcell_integrated_PTIL.rds")

#See genes which are different to connect to mechanism and cell chat part

discriminantGenes_0vs1 <- find.discriminant.genes(ref = ref, query = Tcell_integrated_1, query.control = Tcell_integrated_0,state = "Th1")
discriminantGenes_0vs4 <- find.discriminant.genes(ref = ref, query = Tcell_integrated_4, query.control = Tcell_integrated_0,state = "Th1")
discriminantGenes_0vs24 <- find.discriminant.genes(ref = ref, query = Tcell_integrated_24, query.control = Tcell_integrated_0,state = "Th1")

head(discriminantGenes_0vs1, n = 10)
head(discriminantGenes_0vs4, n = 10)
head(discriminantGenes_0vs24, n = 10)


#Take this information and use to subset CD4 cells and then analyze further---- 
#this is done for CD4 or CD8 specific reference annotation
Tcell_integrated_TIL_ms1<-readRDS("Tcell_integrated_PTIL.rds")
Tcell_integrated_CD4<-subset(Tcell_integrated_TIL_ms1, subset = Tcell_integrated_TIL_ms1$functional.cluster %in% c("Th1","CD4_NaiveLike","Treg","Tfh"))
CD4_cell<-colnames(Tcell_integrated_CD4)
Tcell_integrated<-subset(Integrated_sc, cells=CD4_cell)

list_of_ref.maps <- get.reference.maps(collection = "mouse")
ref <- list_of_ref.maps$mouse$Virus_CD4T #repeated using this

library(Seurat)
library(dplyr)
library(BuenColors)

refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd")
DimPlot(ref, label = T, cols = refCols)

Tcell_integrated@meta.data <- Tcell_integrated@meta.data %>%
  mutate(
    samplid_combined = case_when(
      grepl("^LPS_24H", sampleid) ~ "LPS_24H",
      grepl("^LPS_4H",  sampleid) ~ "LPS_4H",
      grepl("^LPS_1H",  sampleid) ~ "LPS_1H",
      grepl("^NT_0H",   sampleid) ~ "NT_0H",
      TRUE ~ NA_character_
    )
  )
#run query reference 
Tcell_integrated <- Run.ProjecTILs(Tcell_integrated, ref = ref)
Tcell_integrated_0 <- subset(Tcell_integrated, subset = samplid_combined %in% c("NT_0H"))
Tcell_integrated_1 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_1H"))
Tcell_integrated_4 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_4H"))
Tcell_integrated_24 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_24H"))

#save query reference plot
p1<-plot.projection(ref, Tcell_integrated_0, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p2<-plot.projection(ref, Tcell_integrated_1, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p3<-plot.projection(ref, Tcell_integrated_4, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p4<-plot.projection(ref, Tcell_integrated_24, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p1+p2+p3+p4

#fraction plot
p1<-plot.statepred.composition(ref, Tcell_integrated_0, metric = "Percent",cols = jdb_palette("lawhoops"))
p2<-plot.statepred.composition(ref, Tcell_integrated_1, metric = "Percent",cols = jdb_palette("lawhoops"))
p3<-plot.statepred.composition(ref, Tcell_integrated_4, metric = "Percent",cols = jdb_palette("lawhoops"))
p4<-plot.statepred.composition(ref, Tcell_integrated_24, metric = "Percent",cols = jdb_palette("lawhoops"))

p1+p2+p3+p4


#Another way to plot it

library(dplyr)
library(ggplot2)

meta <- Tcell_integrated@meta.data

# count + fraction
df <- meta %>%
  group_by(samplid_combined, functional.cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(samplid_combined) %>%
  mutate(frac = n / sum(n))

df$samplid_combined <- factor(
  df$samplid_combined,
  levels = c("NT_0H","LPS_1H","LPS_4H","LPS_24H")
)

ggplot(df, aes(x = samplid_combined, y = frac, 
               group = functional.cluster, 
               color = functional.cluster)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  ylab("Fraction of cells") +
  xlab("Condition") +
  scale_color_manual(values = jdb_palette("lawhoops")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  facet_wrap(~functional.cluster)

#Compare signatures
library(SignatuR)
data(SignatuR)

programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)

library(UCell)
ref <- AddModuleScore_UCell(ref, features = programs, assay = "RNA", name = NULL)
Tcell_integrated_0 <- AddModuleScore_UCell(Tcell_integrated_0, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_1 <- AddModuleScore_UCell(Tcell_integrated_1, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_4 <- AddModuleScore_UCell(Tcell_integrated_4, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_24 <- AddModuleScore_UCell(Tcell_integrated_24, features = programs, assay = "RNA",name = NULL)

plot.states.radar(ref, query = Tcell_integrated_0, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_1, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_4, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_0, meta4radar = names(programs))


#visualize on the original umap 
Tcell_integrated <- ProjecTILs.classifier(query = Tcell_integrated, ref = ref,skip.normalize=TRUE)
Tcell_integrated$samplid_combined <- factor(
  Tcell_integrated$samplid_combined,
  levels = c("NT_0H", "LPS_1H", "LPS_4H", "LPS_24H")
)
#do it again with cca reduction
DimPlot(Tcell_integrated, group.by = "functional.cluster", split.by="samplid_combined",cols = jdb_palette("lawhoops"))

#save RDS obejct for future
Tcell_integrated_PTIL<-Tcell_integrated
saveRDS(Tcell_integrated_PTIL,file="Tcell_integrated_PTIL_cd4_reference.rds")

#CD8 reference----
Tcell_integrated_TIL_ms1<-readRDS("Tcell_integrated_PTIL.rds")
Tcell_integrated_CD8<-subset(Tcell_integrated_TIL_ms1, subset = Tcell_integrated_TIL_ms1$functional.cluster %in% c("CD8_Tex","CD8_Tpex","CD8_EffectorMemory","CD8_EarlyActiv","CD8_NaiveLike"))
CD8_cell<-colnames(Tcell_integrated_CD8)
Tcell_integrated<-subset(Integrated_sc, cells=CD8_cell)

list_of_ref.maps <- get.reference.maps(collection = "mouse")
ref <- list_of_ref.maps$mouse$Virus_CD8T #repeated using this

library(Seurat)
library(dplyr)
library(BuenColors)

refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd")
DimPlot(ref, label = T, cols = refCols)

Tcell_integrated@meta.data <- Tcell_integrated@meta.data %>%
  mutate(
    samplid_combined = case_when(
      grepl("^LPS_24H", sampleid) ~ "LPS_24H",
      grepl("^LPS_4H",  sampleid) ~ "LPS_4H",
      grepl("^LPS_1H",  sampleid) ~ "LPS_1H",
      grepl("^NT_0H",   sampleid) ~ "NT_0H",
      TRUE ~ NA_character_
    )
  )
#run query reference 
Tcell_integrated <- Run.ProjecTILs(Tcell_integrated, ref = ref)
Tcell_integrated_0 <- subset(Tcell_integrated, subset = samplid_combined %in% c("NT_0H"))
Tcell_integrated_1 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_1H"))
Tcell_integrated_4 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_4H"))
Tcell_integrated_24 <- subset(Tcell_integrated, subset = samplid_combined %in% c("LPS_24H"))

#save query reference plot
p1<-plot.projection(ref, Tcell_integrated_0, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p2<-plot.projection(ref, Tcell_integrated_1, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p3<-plot.projection(ref, Tcell_integrated_4, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p4<-plot.projection(ref, Tcell_integrated_24, linesize = 0.5, pointsize = 0.5,cols = jdb_palette("lawhoops"))
p1+p2+p3+p4

#fraction plot
p1<-plot.statepred.composition(ref, Tcell_integrated_0, metric = "Percent",cols = jdb_palette("lawhoops"))
p2<-plot.statepred.composition(ref, Tcell_integrated_1, metric = "Percent",cols = jdb_palette("lawhoops"))
p3<-plot.statepred.composition(ref, Tcell_integrated_4, metric = "Percent",cols = jdb_palette("lawhoops"))
p4<-plot.statepred.composition(ref, Tcell_integrated_24, metric = "Percent",cols = jdb_palette("lawhoops"))

p1+p2+p3+p4


#Another way to plot it

library(dplyr)
library(ggplot2)

meta <- Tcell_integrated@meta.data

# count + fraction
df <- meta %>%
  group_by(samplid_combined, functional.cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(samplid_combined) %>%
  mutate(frac = n / sum(n))

df$samplid_combined <- factor(
  df$samplid_combined,
  levels = c("NT_0H","LPS_1H","LPS_4H","LPS_24H")
)

ggplot(df, aes(x = samplid_combined, y = frac, 
               group = functional.cluster, 
               color = functional.cluster)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  ylab("Fraction of cells") +
  xlab("Condition") +
  scale_color_manual(values = jdb_palette("lawhoops")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  facet_wrap(~functional.cluster)

#Compare signatures
library(SignatuR)
data(SignatuR)

programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)

library(UCell)
ref <- AddModuleScore_UCell(ref, features = programs, assay = "RNA", name = NULL)
Tcell_integrated_0 <- AddModuleScore_UCell(Tcell_integrated_0, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_1 <- AddModuleScore_UCell(Tcell_integrated_1, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_4 <- AddModuleScore_UCell(Tcell_integrated_4, features = programs, assay = "RNA",name = NULL)
Tcell_integrated_24 <- AddModuleScore_UCell(Tcell_integrated_24, features = programs, assay = "RNA",name = NULL)

plot.states.radar(ref, query = Tcell_integrated_0, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_1, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_4, meta4radar = names(programs))
plot.states.radar(ref, query = Tcell_integrated_0, meta4radar = names(programs))


#visualize on the original umap 
Tcell_integrated <- ProjecTILs.classifier(query = Tcell_integrated, ref = ref,skip.normalize=TRUE)
Tcell_integrated$samplid_combined <- factor(
  Tcell_integrated$samplid_combined,
  levels = c("NT_0H", "LPS_1H", "LPS_4H", "LPS_24H")
)
#do it again with cca reduction
DimPlot(Tcell_integrated, group.by = "functional.cluster", split.by="samplid_combined",cols = jdb_palette("lawhoops"))

#save RDS obejct for future
Tcell_integrated_PTIL<-Tcell_integrated
saveRDS(Tcell_integrated_PTIL,file="Tcell_integrated_PTIL_cd4_reference.rds")


