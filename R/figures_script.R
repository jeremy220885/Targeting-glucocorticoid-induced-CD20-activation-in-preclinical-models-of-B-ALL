# =============================================================================
# Figure generation script
# Targeting glucocorticoid-induced CD20 activation in preclinical models of B-ALL
# Input: All_merge1_integrated.RData (Seurat object, post QC/clustering/integration)
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(ggridges)
library(EnhancedVolcano)
library(UCell)
library(ggpubr)
library(dplyr)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
load("All_merge1_integrated.RData")
All_merge1.integrated[["RNA"]] <- JoinLayers(All_merge1.integrated[["RNA"]])

# =============================================================================
# FIGURE — UMAP / DimPlot (Fig. 1)
# =============================================================================

# Global UMAP, colored by cluster
DimPlot(All_merge1.integrated, label = TRUE)

# Remove cluster 12 (mouse cells contamination)
All_merge1.integrated <- subset(All_merge1.integrated, idents = c(0:11))
DimPlot(All_merge1.integrated, label = TRUE)

# UMAP split by the 4 sample types (Diag_NT, Diag_T, Rechute_NT, Rechute_T)
DimPlot(All_merge1.integrated, label = TRUE, split.by = "SampleType")

# =============================================================================
# FIGURE — Volcano plots (EnhancedVolcano) 
# Differential expression between conditions
# =============================================================================

Idents(All_merge1.integrated) <- All_merge1.integrated$SampleType

# Diag_T vs Diag_NT
Diag_TvsDiag_NT <- FindMarkers(All_merge1.integrated,
                                ident.1 = "Diag_T",
                                ident.2 = "Diag_NT")
write.csv2(Diag_TvsDiag_NT, file = "Diag_TvsDiag_NT.csv", row.names = TRUE)

EnhancedVolcano(Diag_TvsDiag_NT,
                lab = rownames(Diag_TvsDiag_NT),
                x = "avg_log2FC", y = "p_val_adj",
                pCutoff = 0.05, FCcutoff = 2,
                pointSize = 0.1, labSize = 2.5,
                title = "Diag_T vs Diag_NT")

# Rechute_T vs Rechute_NT
rechute_Tvsrechute_NT <- FindMarkers(All_merge1.integrated,
                                      ident.1 = "Rechute_T",
                                      ident.2 = "Rechute_NT")
write.csv2(rechute_Tvsrechute_NT, file = "rechute_Tvsrechute_NT.csv", row.names = TRUE)

EnhancedVolcano(rechute_Tvsrechute_NT,
                lab = rownames(rechute_Tvsrechute_NT),
                x = "avg_log2FC", y = "p_val_adj",
                pCutoff = 0.05, FCcutoff = 1,
                pointSize = 1, labSize = 2.5,
                title = "Rechute_T vs Rechute_NT")

# =============================================================================
# FIGURE — VlnPlot / FeaturePlot (Fig. 2)
# =============================================================================

# Violin plot — CD20 (MS4A1) grouped by sample type
VlnPlot(All_merge1.integrated, features = "MS4A1",
        group.by = "SampleType", pt.size = FALSE)

# Violin plot — CD34 grouped by sample type
VlnPlot(All_merge1.integrated, features = "CD34",
        group.by = "SampleType", pt.size = FALSE)

# FeaturePlot CD20 and CD34 split by sample type
Idents(All_merge1.integrated) <- "seurat_clusters"
FeaturePlot(object = All_merge1.integrated, label = TRUE,
            features = c("MS4A1"), split.by = "SampleType") # CD20
FeaturePlot(object = All_merge1.integrated, label = TRUE,
            features = c("CD34"), split.by = "SampleType")  # CD34

# =============================================================================
# FIGURE — Cell cycle scoring (Fig. S3A)
# =============================================================================

cc.genes_tinyatlas <- list(
  s.genes_tinyatlas = c("UBR7", "RFC2", "RAD51", "MCM2", "TIPIN", "MCM6", "UNG", "POLD3", "WDR76", "CLSPN", "CDC45", "CDC6", "MSH2", "MCM5", "POLA1", "MCM4", "RAD51AP1", "GMNN", "RPA2", "CASP8AP2", "HELLS", "E2F8", "GINS2", "PCNA", "NASP", "BRIP1", "DSCC1", "DTL", "CDCA7", "CENPU", "ATAD2", "CHAF1B", "USP1", "SLBP", "RRM1", "FEN1", "RRM2", "EXO1", "CCNE2", "TYMS", "BLM", "PRIM1", "UHRF1"),
  g2m.genes_tinyatlas = c("NCAPD2", "ANLN", "TACC3", "HMMR", "GTSE1", "NDC80", "AURKA", "TPX2", "BIRC5", "G2E3", "CBX5", "RANGAP1", "CTCF", "CDCA3", "TTK", "SMC4", "ECT2", "CENPA", "CDC20", "NEK2", "CENPF", "TMPO", "HJURP", "CKS2", "DLGAP5", "PIMREG", "TOP2A", "PSRC1", "CDCA8", "CKAP2", "NUSAP1", "KIF23", "KIF11", "KIF20B", "CENPE", "GAS2L3", "KIF2C", "NUF2", "ANP32E", "LBR", "MKI67", "CCNB2", "CDC25C", "HMGB2", "CKAP2L", "BUB1", "CDK1", "CKS1B", "UBE2C", "CKAP5", "AURKB", "CDCA2", "TUBB4B", "JPT1"))

s.genes_tinyatlas <- cc.genes_tinyatlas$s.genes_tinyatlas
g2m.genes_tinyatlas <- cc.genes_tinyatlas$g2m.genes_tinyatlas

All_merge1.integrated <- RunPCA(All_merge1.integrated,
                                 features = VariableFeatures(All_merge1.integrated),
                                 ndims.print = 6:10, nfeatures.print = 10)
DimHeatmap(All_merge1.integrated, dims = c(8, 10))

All_merge1.integrated <- CellCycleScoring(All_merge1.integrated,
                                           s.features   = s.genes_tinyatlas,
                                           g2m.features = g2m.genes_tinyatlas,
                                           set.ident    = TRUE)

# View cell cycle scores and phase assignments
head(All_merge1.integrated[[]])

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells
# separate entirely by phase
All_merge1.integrated <- RunPCA(All_merge1.integrated,
                                 features = c(s.genes_tinyatlas, g2m.genes_tinyatlas))
DimPlot(All_merge1.integrated)
DimPlot(All_merge1.integrated, split.by = "SampleType")
All_merge1.integrated$Phase <- as.character(All_merge1.integrated$Phase)
All_merge1.integrated$Phase[All_merge1.integrated$Phase == "G1"] <- "G0/G1"
DimPlot(All_merge1.integrated, group.by = "Phase")

# Fig. S3A — FeaturePlot of MKI67 expression, split by condition
FeaturePlot(All_merge1.integrated, features = "MKI67", split.by = "SampleType")

# Restore cluster identities for downstream figures
Idents(All_merge1.integrated) <- All_merge1.integrated$seurat_clusters

# =============================================================================
# FIGURE — Gene signature scoring (UCell) 
# =============================================================================

# -----------------------------------------------------------------------------
# Gene signature lists
# -----------------------------------------------------------------------------

markers_quies <- list()

# MRD signature — up-regulated genes (Ebinger et al., Cancer Cell 2016)
# DEG between primary diagnosis and primary MRD (111 genes)
markers_quies$UP_DEG_between_primary_diagnosis_and_primary_MRD_111 <- c(
  "RIN3", "CHRM3-AS2", "P2RX5", "SLA", "C16orf45", "SORT1", "ARHGAP18",
  "CPM", "GPR183", "IL10RA", "DFNA5", "TMEM100", "CD163", "TYROBP",
  "LGALS3", "SIGLEC15", "HBA1", "IL1R2", "GSN", "CTC-378H22.2", "LILRA1",
  "IGHA1", "DEPTOR", "SGK1", "LINC00961", "ADAM8", "MYRIP", "MIR22HG",
  "VOPP1", "FKBP5", "TNFRSF18", "NT5E", "IL10", "S100A6", "RAB20",
  "F2RL3", "UPP1", "ENSG00000117289", "JDP2", "MPEG1", "PFKP", "DUSP5",
  "TTN", "SIRPA", "SIK1", "IGLC3", "ENSG00000183748", "SLCO4A1", "RNASE6",
  "HMOX1", "CSGALNACT1", "CDC42EP3", "F3", "HBB", "CTTNBP2NL", "GBP2",
  "IFNGR1", "HBA2", "CEBPB", "SSH1", "S100A11", "FGD2", "ITPKB", "CD1C",
  "SEL1L3", "MXD4", "SCN1B", "MS4A1", "SMAP2", "CCDC107", "ARRDC3",
  "GLIPR2", "IQSEC1", "HERPUD1", "RP11-458D21.5", "APBB1", "SPRY1",
  "UNC50", "ETV5", "RHBDF2", "OGFRL1", "SELM", "SNX8", "LSP1", "NOTCH2",
  "TFEB", "ENDOD1", "ITGB2", "CARD19", "CALHM2", "PCTP", "PSTPIP1",
  "COMMD5", "SMARCA2", "FHOD1", "NCOA7", "NEDD9", "SERINC1", "ZNF331",
  "TACC1", "HCST", "CYTH4", "CD53", "NEAT1", "MALAT1", "HVCN1", "CTSB",
  "ITM2C", "MVP", "NPTN", "MEF2A"
)

# MRD signature — down-regulated genes (Ebinger et al., Cancer Cell 2016)
# DEG between primary diagnosis and primary MRD (95 genes)
markers_quies$down_DEG_between_primary_diagnosis_and_primary_MRD_95 <- c(
  "MZT2A", "EIF3F", "RP11-673C5.1", "DNMT1", "HMGB1", "SMARCA4", "SSRP1",
  "NUCB2", "HNRNPA1", "LDHB", "LIG1", "HMGA1", "H2AFY", "IDH2", "LGALS9",
  "CENPM", "RRM1", "MCM3", "ERG", "TUBA1A", "TUBA1B", "HIST1H1C", "PSMA6P1",
  "NASP", "TMEM263", "KIAA0101", "CHAF1A", "BIRC5", "FEN1", "SPTBN1",
  "HMGN2", "DUT", "TK1", "MCM7", "CTPS1", "MYL6B", "FHIT", "MAD2L1",
  "CDCA4", "CDK2", "MYB", "TRMT5", "C5orf56", "SNHG1", "TKT", "MND1",
  "CPXM1", "MCM5", "SOX4", "MDK", "GAMT", "IGLL1", "MME", "FAM129C",
  "MSH6", "CYTL1", "CHAF1B", "SCMH1", "CTD-2006C1.2", "C12orf75", "SMIM24",
  "CDCA5", "KIF11", "HPS4", "CDT1", "DNTT", "DBN1", "E2F1", "ADA",
  "PKMYT1", "TRH", "C4orf46", "AP005530.2", "KLF4", "MCM4", "ZWINT",
  "ENSG00000168274", "STMN1", "PCNA", "MCM6", "APBB2", "TUBBP1", "POLE",
  "CENPH", "KCNK12", "CDC45", "MCM2", "CBX2", "ENSG00000034063",
  "PLEKHG4B", "DHFR", "MACROD2", "NREP", "TYMS", "CDCA7"
)

# Cell-quiescence signature — up-regulated genes (Min et al., PLoS Biol 2019, 70 genes)
markers_quies$Quiescence_up_Mingwei_plos_2019_70 <- c(
  "CFLAR", "CD9", "CALCOCO1", "ROS1", "COL17A1", "FXYD3", "GPATCH2L",
  "YPEL3", "CST3", "NIPAL2", "MAN2B1", "RPS19", "SERINC1", "PERP",
  "CLIP4", "PCYOX1", "TMEM59", "RGS2", "IRF6", "PPL", "YPEL5", "AHNAK",
  "SAT1", "INADL", "DSC2", "DSG3", "CD63", "KIAA1109", "CDH13", "SYTL1",
  "GSN", "CCL28", "MR1", "CD109", "NMNAT2", "CYB5R1", "AZGP1", "NCSTN",
  "ZFYVE1", "SCN3B", "TMC4", "LENG8", "KRT15", "DMXL1", "NADSYN1",
  "CHST2", "GOLGA8A", "EPS8L2", "ALS2CL", "PTTG1IP", "TACSTD2", "ANO9",
  "AHNAK2", "MIR22HG", "PSAP", "S100A6", "NA", "NA", "GOLGA8A", "NA",
  "NA", "NA", "NA", "NEAT1", "MALAT1", "CHMP1B", "NOTCH2NL", "TXNIP",
  "MTRNR2L2", "NBPF14"
)

# Cell-quiescence signature — down-regulated genes (Min et al., PLoS Biol 2019, 128 genes)
markers_quies$Quiescence_down_Mingwei_plos_2019_128 <- c(
  "NCAPD2", "PTBP1", "MPHOSPH9", "NUCKS1", "TCOF1", "SMC1A", "MCM2",
  "SART3", "SNRPA", "KIF22", "HSP90AA1", "WBP11", "CAD", "SF3B2",
  "KHSRP", "WDR76", "NUP188", "HSP90AB1", "HNRNPM", "SMARCB1", "MCM5",
  "PNN", "RBBP7", "NPRL3", "USP10", "MCM4", "SGTA", "MRPL4", "PSMD3",
  "KPNB1", "CBX1", "LRRC59", "TMEM97", "WHSC1", "PRPF19", "PTGES3",
  "CPSF6", "SRSF3", "MCM3", "TCERG1", "SMC4", "EIF4G1", "ZNF142",
  "MSH6", "MRPL37", "SFPQ", "STMN1", "ARID1A", "PROSER1", "DDX39A",
  "EXOSC9", "USP22", "DEK", "DUT", "ILF3", "DNMT1", "PCNA", "NASP",
  "NA", "SRRM1", "CCNB1", "GNL2", "RNF138", "SRSF1", "TRA2B", "SMPD4",
  "ANP32B", "HMGA1", "MDC1", "HADH", "TP53", "ARHGDIA", "PRCC", "HDGF",
  "SF3B4", "UBAP2L", "ILF2", "PARP1", "LBR", "RQCD1", "SKP2", "MMS22L",
  "PPRC1", "SSRP1", "CCT5", "DLAT", "HNRNPU", "LARP1", "SCAF4", "RRP1B",
  "RRP1", "CHCHD4", "GMPS", "RFC4", "SLBP", "CDC25A", "PSIP1", "HNRNPK",
  "SKA3", "DIS3L", "USP39", "GPS1", "PA2G4", "HCFC1", "SLC19A1", "ETV4",
  "RAD23A", "DCTPP1", "RCC1", "EWSR1", "ALYREF", "PTMA", "HMGB1",
  "POM121", "MCMBP", "TEAD4", "TFDP1", "CHAMP1", "TOP1", "PRRC2A",
  "NA", "RBM14", "NA", "NA", "NA", "POM121C", "NA", "UHRF1"
)

# BCR signaling pathway — KEGG database (75 genes)
markers_quies$Kegg_Bcell_receptor_pathway <- c(
  "AKT1", "AKT2", "AKT3", "BCL10", "BLNK", "BTK", "CARD11", "CD19",
  "CD22", "CD72", "CD79A", "CD79B", "CD81", "CHP1", "CHP2", "CHUK",
  "CR2", "DAPP1", "FCGR2B", "FOS", "GRB2", "GSK3B", "HRAS", "IFITM1",
  "IKBKB", "IKBKG", "INPP5D", "JUN", "KRAS", "LILRB3", "LYN", "MALT1",
  "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "NFAT5", "NFATC1", "NFATC2",
  "NFATC3", "NFATC4", "NFKB1", "NFKBIA", "NFKBIB", "NFKBIE", "NRAS",
  "PIK3AP1", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2",
  "PIK3R3", "PIK3R5", "PLCG2", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1",
  "PPP3R2", "PRKCB", "PTPN6", "RAC1", "RAC2", "RAC3", "RAF1", "RASGRP3",
  "RELA", "SOS1", "SOS2", "SYK", "VAV1", "VAV2", "VAV3"
)

# BCR signaling pathway — Gene Ontology database (78 genes)
# GO:0050853 — B cell receptor signaling pathway
# Gene list provided by Anne Quillet
markers_quies$GO_BCR0050853 <- c(
  "ABL1", "BANK1", "BAX", "BCAR1", "BCL2", "BLK", "BLNK", "BMX", "BTK",
  "CD19", "CD22", "CD300A", "CD38", "CD79A", "CD79B", "CD81", "CMTM3",
  "CTLA4", "FCGR2B", "FCMR", "FCRL3", "FOXP1", "GCSAM", "GCSAML", "GPS2",
  "GRB2", "IGHA1", "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3",
  "IGHG4", "IGHM", "IGKC", "IGLC1", "IGLC3", "IGLC6", "IGLC7", "ITK",
  "KLHL6", "LAT2", "LCK", "LIME1", "LPXN", "LYN", "MAPK1", "MEF2C",
  "MIR18A", "MIR19A", "MIR34A", "MNDA", "MS4A1", "NCKAP1L", "NFAM1",
  "NFATC2", "NFKB1", "NFKBIA", "NFKBIZ", "PIK3CD", "PLCG2", "PLCL2",
  "PLEKHA1", "PRKCB", "PRKCH", "PTPN22", "PTPN6", "PTPRC", "RFTN1",
  "SH2B2", "SLC39A10", "SOS1", "STAP1", "SYK", "TEC", "TMIGD2", "VAV3"
)

# -----------------------------------------------------------------------------
# UCell scoring
# -----------------------------------------------------------------------------

All_merge1.integrated <- AddModuleScore_UCell(All_merge1.integrated,
                                               features = markers_quies,
                                               name = "_Ucell")
signature.names <- paste0(names(markers_quies), "_Ucell")

# -----------------------------------------------------------------------------
# Fig. S3C / S3D — RidgePlot by cluster
# -----------------------------------------------------------------------------

listDiagNT <- grep("Diag_NT", colnames(All_merge1.integrated), value = TRUE)
All_merge1_DNT.integrated <- subset(All_merge1.integrated, cells = listDiagNT)

for (feature_name in signature.names) {
  p <- RidgePlot(All_merge1_DNT.integrated,
                 group.by = "seurat_clusters",
                 features = feature_name) +
    theme(plot.title = element_text(size = 10),
          axis.text  = element_text(size = 8)) +
    ggtitle(paste("Diag_NT —", feature_name))
  print(p)
}

# -----------------------------------------------------------------------------
# Fig. S3E — Fold enrichment of cluster 10 after CHEMO (Diagnosis vs Relapse)
# -----------------------------------------------------------------------------

meta <- All_merge1.integrated@meta.data
meta$is_cluster10 <- meta$seurat_clusters == 10

prop_cluster10 <- meta %>%
  group_by(SampleType) %>%
  summarise(prop = sum(is_cluster10) / n(), .groups = "drop")

fc_diag <- prop_cluster10$prop[prop_cluster10$SampleType == "Diag_T"] /
           prop_cluster10$prop[prop_cluster10$SampleType == "Diag_NT"]
fc_rech <- prop_cluster10$prop[prop_cluster10$SampleType == "Rechute_T"] /
           prop_cluster10$prop[prop_cluster10$SampleType == "Rechute_NT"]

df_S3E <- data.frame(
  Timepoint = factor(c("Diagnosis", "Relapse"), levels = c("Diagnosis", "Relapse")),
  log2FC    = log2(c(fc_diag, fc_rech))
)

p_S3E <- ggplot(df_S3E, aes(x = Timepoint, y = log2FC, fill = Timepoint)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Diagnosis" = "#B07BB5", "Relapse" = "#B07BB5")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(title = "Cluster 10",
       y = "Fold enrichment after CHEMO (log2FC)",
       x = "")
print(p_S3E)

# -----------------------------------------------------------------------------
# Fig. S3F — RidgePlot of MRD up, Quiescence up, BCR GO and BCR KEGG
# split by MS4A1 status (MS4A1+ vs MS4A1-) on Diag_NT
# -----------------------------------------------------------------------------

All_merge1_DNT.integrated$MS4A1_status <- ifelse(
  GetAssayData(All_merge1_DNT.integrated, assay = "RNA", layer = "data")["MS4A1", ] > 1,
  "MS4A1+",
  "MS4A1-"
)

sig_S3F <- c(
  "UP_DEG_between_primary_diagnosis_and_primary_MRD_111_Ucell",
  "Quiescence_up_Mingwei_plos_2019_70_Ucell",
  "GO_BCR0050853_Ucell"
  
)

for (feature_name in sig_S3F) {
  p <- RidgePlot(All_merge1_DNT.integrated,
                 group.by = "MS4A1_status",
                 features = feature_name) +
    theme(plot.title = element_text(size = 10),
          axis.text  = element_text(size = 8)) +
    ggtitle(paste("MS4A1 status —", feature_name))
  print(p)
}

# =============================================================================
# FIGURE — SPIB vs CD20 (MS4A1) correlation 
#   - jitter added to spread points along discrete expression values
# =============================================================================

df <- FetchData(All_merge1.integrated, vars = c("SPIB", "MS4A1", "SampleType"))
df <- na.omit(df)

# Keep only cells positive for both SPIB and MS4A1
df <- subset(df, SPIB > 0 & MS4A1 > 0)

df$SampleType <- factor(df$SampleType,
                        levels = c("Diag_NT", "Diag_T", "Rechute_NT", "Rechute_T"),
                        labels = c("Diagnosis Vehicle", "Diagnosis CHEMO",
                                   "Relapse Vehicle",   "Relapse CHEMO"))

df$SPIB_log  <- log1p(df$SPIB)
df$MS4A1_log <- log1p(df$MS4A1)

p_corr <- ggplot(df, aes(x = SPIB_log, y = MS4A1_log)) +
  geom_jitter(width = 0.02, height = 0.02, alpha = 0.25, size = 0.4, color = "grey30") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick3", linewidth = 0.8) +
  stat_cor(method = "pearson",
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           label.x.npc = "left", label.y.npc = "top", size = 4) +
  facet_wrap(~ SampleType, ncol = 2) +
  theme_bw(base_size = 15) +
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    strip.text       = element_text(size = 16, face = "bold"),
    panel.border     = element_rect(color = "black", linewidth = 0.7),
    axis.title       = element_text(size = 15),
    axis.text        = element_text(size = 12)
  ) +
  labs(x = "SPIB expression (log1p, SPIB+ cells only)",
       y = "MS4A1 expression (log1p, MS4A1+ cells only)",
       title = "Correlation between SPIB and CD20 expression")

print(p_corr)
ggsave(filename = "SPIB_vs_CD20_correlation.pdf", plot = p_corr, width = 8, height = 8)

# Pearson correlation per group
for (grp in levels(df$SampleType)) {
  r <- with(subset(df, SampleType == grp), cor(SPIB, MS4A1, method = "pearson"))
  cat(grp, ": r =", r, "\n")
}

#Figure S7B
FeaturePlot(object = All_merge1.integrated, label = TRUE,
            features = c("SPIB"), split.by = "SampleType")  # SPIB

# =============================================================================
# FIGURE 5D — Normalized expression of 10 transcription factors
# =============================================================================

library(tidyr)

tf_genes <- c("MS4A1", "SPIB", "SOX7", "AHRR", "DPF3",
              "MAF", "HIC1", "TCF7", "SCML4", "ZNF831")

Idents(All_merge1.integrated) <- All_merge1.integrated$SampleType

avg_tf <- AverageExpression(All_merge1.integrated,
                             features = tf_genes,
                             group.by = "SampleType",
                             assays = "RNA",
                             slot = "data")$RNA

# Rename columns (Seurat may replace underscores with dashes)
colnames(avg_tf) <- gsub("-", "_", colnames(avg_tf))
avg_tf <- avg_tf[, c("Diag_NT", "Diag_T", "Rechute_NT", "Rechute_T")]

# Reshape to long format
df_tf <- as.data.frame(avg_tf)
df_tf$gene <- rownames(df_tf)

df_tf_long <- pivot_longer(df_tf,
                            cols = c("Diag_NT", "Diag_T", "Rechute_NT", "Rechute_T"),
                            names_to = "Condition",
                            values_to = "Expression")

df_tf_long$Condition <- factor(df_tf_long$Condition,
                                levels = c("Diag_NT", "Diag_T", "Rechute_NT", "Rechute_T"),
                                labels = c("Diagnosis Vehicle", "Diagnosis CHEMO",
                                           "Relapse Vehicle", "Relapse CHEMO"))

df_tf_long$gene <- factor(df_tf_long$gene, levels = tf_genes)

# Shared Y axis limit across all genes
y_max <- max(df_tf_long$Expression) * 1.05

condition_colors <- c("Diagnosis Vehicle" = "#E07B7B",
                       "Diagnosis CHEMO"   = "#6BB06B",
                       "Relapse Vehicle"   = "#6BAED6",
                       "Relapse CHEMO"     = "#9E7FC7")

# One plot per gene, same Y scale
for (g in tf_genes) {
  df_gene <- subset(df_tf_long, gene == g)

  p <- ggplot(df_gene, aes(x = Condition, y = Expression, fill = Condition)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = condition_colors) +
    scale_y_continuous(limits = c(0, y_max)) +
    theme_bw(base_size = 12) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title   = element_text(face = "italic", size = 12),
          legend.title = element_blank()) +
    labs(title = g, y = "Normalized expression")

  print(p)
}

############################################################
# Analysis:Fig. S7G
# Expression heatmap of MRD-upregulated genes associated
# with published SPIB occupancy
############################################################


############################################################
# Gene list
############################################################

spib_mrd_genes <- c(
  "ECE1",
  "ITGAM",
  "AHRR",
  "KREMEN2",
  "CR2",
  "DPF3",
  "RAPGEF5",
  "TCF7",
  "MS4A1",
  "SLC22A23",
  "CR1L",
  "MGAT3",
  "AK5",
  "GRASP",
  "PGF",
  "SSH3",
  "GNAZ",
  "PTGDR",
  "PLEKHG1",
  "NCAM2",
  "MYO6",
  "VAT1L",
  "ACTA2",
  "SLC2A6",
  "RNF157",
  "SUCLG2",
  "GBP1",
  "SPOCK2",
  "APOBEC3H",
  "LIPE",
  "EFNB2",
  "ALDH1A1",
  "TTC9",
  "SYTL2",
  "HTR3A",
  "SHROOM3",
  "PRX",
  "MFAP5",
  "FCRL1",
  "CD300A",
  "SCML4",
  "EPHA4",
  "NOS3",
  "ANPEP",
  "IL18RAP",
  "MAP3K9",
  "CD70",
  "FSTL4"
)

############################################################
# Basic checks
############################################################

length(spib_mrd_genes)
any(duplicated(spib_mrd_genes))
spib_mrd_genes[duplicated(spib_mrd_genes)]

genes_present <- spib_mrd_genes[
  spib_mrd_genes %in% rownames(All_merge1.integrated)
]

genes_missing <- spib_mrd_genes[
  !spib_mrd_genes %in% rownames(All_merge1.integrated)
]

length(genes_present)
genes_missing

write.csv(
  data.frame(
    gene = spib_mrd_genes,
    present_in_object = spib_mrd_genes %in% rownames(All_merge1.integrated)
  ),
  file = "01_gene_presence_check.csv",
  row.names = FALSE
)

############################################################
# Check condition names
############################################################

table(All_merge1.integrated$SampleType)

############################################################
# Average expression per condition
############################################################

Idents(All_merge1.integrated) <- All_merge1.integrated$SampleType

avg_expr <- AverageExpression(
  All_merge1.integrated,
  features = genes_present,
  group.by = "SampleType",
  assays = "RNA",
  slot = "data"
)$RNA

# Seurat may replace underscores by dashes in column names
colnames(avg_expr)

# Reorder columns
avg_expr <- avg_expr[, c("Diag-NT", "Diag-T", "Rechute-NT", "Rechute-T")]

# Rename columns for readability
colnames(avg_expr) <- c("Diag_NT", "Diag_T", "Rechute_NT", "Rechute_T")

############################################################
#  Heatmap with row scaling
############################################################
install.packages("pheatmap")
library(pheatmap)

pheatmap(
  avg_expr,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 7,
  fontsize_col = 11,
  angle_col = 45,
  main = "SPIB-bound MRD-upregulated genes",
  width = 6,
  height = 10
)

#############################################################
#vlnplot SPIB and MS4A1 Fif. S7C
#############################################################

library(ggplot2)

# condition
condition_to_plot <- "Diag_T"

# Subset
obj_cond <- subset(
  All_merge1.integrated,
  subset = SampleType == condition_to_plot
)

# Récupérer MS4A1, SPIB et les clusters
df <- FetchData(
  obj_cond,
  vars = c("MS4A1", "SPIB", "seurat_clusters")
)

df_long <- rbind(
  data.frame(
    seurat_clusters = df$seurat_clusters,
    gene = "MS4A1",
    expression = df$MS4A1
  ),
  data.frame(
    seurat_clusters = df$seurat_clusters,
    gene = "SPIB",
    expression = df$SPIB
  )
)

df_long$seurat_clusters <- factor(
  df_long$seurat_clusters,
  levels = sort(unique(as.numeric(as.character(df_long$seurat_clusters))))
)

df_long$gene <- factor(df_long$gene, levels = c("MS4A1", "SPIB"))

# Plot
p <- ggplot(
  df_long,
  aes(x = seurat_clusters, y = expression, fill = gene)
) +
  geom_violin(
    position = position_dodge(width = 0.8),
    scale = "width",
    trim = TRUE,
    linewidth = 0.2
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey85", color = "black")
  ) +
  labs(
    title = condition_to_plot,
    y = "Expression level"
  )
p