# =========================================================================
# Title: Integrated GO and KEGG Pathway Analysis for ST-CS Selected Proteins
# Description: 
#   This script performs GO and KEGG enrichment analyses on all proteins
#   selected by ST-CS in intrahepatic cholangiocarcinoma (CPTAC PDC000356) 
#   and glioblastoma (CPTAC PDC000446). 
#   Results are visualized in combined plots.
# =========================================================================

# ----------------------------
# Step 1: Load Libraries
# ----------------------------
library(clusterProfiler)    # For enrichment analysis
library(org.Hs.eg.db)       # Human gene annotation database
library(enrichplot)         # Visualization of enrichment results
library(ggplot2)            # Plot customization
library(cowplot)            # Combining multiple plots

# ----------------------------
# Step 2: Input Protein Lists (All Proteins, N≥1)
# ----------------------------

# Intrahepatic Cholangiocarcinoma (PDC000356): All selected proteins
intrahepatic_cholangiocarcinoma <- c(
  "AADAT", "ACSM2A", "ACSM2B", "ECI2", "GAS2", "GLS2", "MMUT", "NDRG2",
  "RMDN2", "SCP2", "SULT1A1", "SULT1A2", "SULT2A1", "TKFC", "TTC36",
  "ACSM3", "ANGPTL3", "CLEC4G", "CYB5A", "DDT", "FABP1", "LDHD", "PALMD",
  "STBD1", "ACMSD", "ACSM1", "EXOC3L4", "S100A11", "SELENBP1", "SULT1A3",
  "ACAA1", "ACBD4", "AK4", "ALDH4A1", "CES1", "ETFRF1", "FXYD1", "GLYCTK",
  "MMAB", "PCK2", "PIPOX", "REEP6", "RIDA", "S100A6", "TFR2", "ACSL1",
  "ALDH7A1", "AS3MT", "ASGR1", "ASGR2", "BDH1", "CAPG", "CLEC4M", "CYP2C9",
  "CYP2E1", "DHRS11", "ESPN", "HAO1", "LAMB1", "LYVE1", "OIT3", "PAH",
  "PGRMC1", "PHYH", "PLIN5", "PPP1R1A", "PVALB", "SPINT1", "VSNL1"
)

# Glioblastoma (PDC000446): All selected proteins
glioblastoma <- c(
  "ADAP1", "ATP5F1E", "BIN1", "GOT2", "TINAGL1", "ATP5IF1", "DLG2", "FTL",
  "NDUFA2", "ARRB1", "AUH", "EFHD1", "MAP6D1", "NAPEPLD", "PDE1A", "ADCY5",
  "AK5", "ASAH1", "ATP6V1G2", "CPNE6", "FBXO2", "HDAC11", "IQSEC2", "MAD2L1",
  "MBLAC2", "MYO1D", "PCSK2", "PLCL1", "PLEKHA1", "PPM1H", "RBP7", "SEPT4",
  "SIGLEC1", "STXBP6", "SUCLA2", "TPRG1L", "ACO2", "ALDH2", "ATP1B1", "ATP6V1D",
  "ATP6V1E1", "ATP8A1", "CAMKK1", "CARNS1", "CD55", "CD93", "CEP55", "ENPP2",
  "FAAH", "FTH1", "GJA1", "GLS", "HK1", "HSPA12A", "IDH3A", "IGSF8", "IQSEC1",
  "JAM3", "KCNAB2", "LYNX1", "LYRM1", "LYRM9", "ME3", "MKI67", "NCAPD2", "NECAB1",
  "NQO1", "NSF", "OGDHL", "PCNA", "PCSK1N", "PEX5L", "PIP4K2A", "PODXL", "PRKAR2B",
  "PRODH", "RASGRF2", "RGS14", "SH3GLB2", "SMC2", "TGIF1", "TMEM167A", "TPPP"
)

# ----------------------------
# Step 3: Convert Gene Symbols to Entrez ID (with error handling)
# ----------------------------

convert_to_entrez <- function(gene_symbols) {
  result <- tryCatch({
    bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  }, error = function(e) {
    message("Error in gene symbol conversion: ", e$message)
    return(NULL)
  })
  return(result)
}

# Convert protein lists
chol_entrez <- convert_to_entrez(intrahepatic_cholangiocarcinoma)
glio_entrez <- convert_to_entrez(glioblastoma)

# ----------------------------
# Step 4: Perform Enrichment Analysis (with error handling)
# ----------------------------

perform_enrichment <- function(entrez_ids, analysis_type = "GO", title_suffix) {
  if (is.null(entrez_ids)) {
    message("Skipping analysis: Empty Entrez IDs.")
    return(NULL)
  }
  
  if (analysis_type == "GO") {
    result <- enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,  # Relaxed threshold
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  } else if (analysis_type == "KEGG") {
    result <- enrichKEGG(
      gene = entrez_ids$ENTREZID,
      organism = "hsa",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    )
  }
  
  if (nrow(result) == 0) {
    message("No significant terms found for ", title_suffix)
    return(NULL)
  }
  return(result)
}

# ----------------------------
# Step 5: Run Analyses for Both Datasets
# ----------------------------

# Intrahepatic Cholangiocarcinoma
go_chol <- perform_enrichment(chol_entrez, "GO", "Intrahepatic Cholangiocarcinoma")
kegg_chol <- perform_enrichment(chol_entrez, "KEGG", "Intrahepatic Cholangiocarcinoma")

# Glioblastoma
go_glio <- perform_enrichment(glio_entrez, "GO", "Glioblastoma")
kegg_glio <- perform_enrichment(glio_entrez, "KEGG", "Glioblastoma")

# ----------------------------
# Step 6: Visualize Combined Results
# ----------------------------

# Custom plotting function
plot_enrichment <- function(enrich_result, title, plot_type = "dot") {
  if (is.null(enrich_result)) {
    return(NULL)
  }
  
  if (plot_type == "dot") {
    p <- dotplot(enrich_result, showCategory = 10, font.size = 10) + 
      ggtitle(title) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  } else if (plot_type == "bar") {
    p <- barplot(enrich_result, showCategory = 10, font.size = 10) + 
      ggtitle(title) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }
  return(p)
}

# Generate plots
plot_list <- list(
  plot_enrichment(go_chol, "GO: Intrahepatic Cholangiocarcinoma", "dot"),
  plot_enrichment(kegg_chol, "KEGG: Intrahepatic Cholangiocarcinoma", "bar"),
  plot_enrichment(go_glio, "GO: Glioblastoma", "dot"),
  plot_enrichment(kegg_glio, "KEGG: Glioblastoma", "bar")
)

# Remove NULL elements (failed analyses)
plot_list <- plot_list[!sapply(plot_list, is.null)]

# Combine plots using cowplot
combined_plot <- plot_grid(plotlist = plot_list, ncol = 2, labels = "AUTO")

# Save combined plot
ggsave("Outputs/Figures/Combined_Enrichment.tiff", combined_plot, 
       width = 14, height = 10, dpi = 300, compression = "lzw")

# ----------------------------
# Step 7: Save Results to CSV
# ----------------------------

write.csv(go_chol@result, "Outputs/Results/GO_Cholangiocarcinoma.csv")
write.csv(go_glio@result, "Outputs/Results/GO_Glioblastoma.csv")
write.csv(kegg_chol@result, "Outputs/Results/KEGG_Cholangiocarcinoma.csv")
write.csv(kegg_glio@result, "Outputs/Results/KEGG_Glioblastoma.csv")