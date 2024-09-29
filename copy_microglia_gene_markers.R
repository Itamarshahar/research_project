get_microglia_markers_copy <- function(full = FALSE) {
  if (full) {
    atlas <- list(
      # Microglia grouped subsets with marker genes and genes of interest
      microglia = list(
        General = list(genes = c("C3")),

        Profilerative = list(
          states = list(Mic.1 = c("ARHGAP11B", "TOP2A", "HMMR"))
        ),

        Surveilling = list(
          genes = c("P2RY12", "CX3CR1", "FRMD4A", "SYNDIG1"),
          states = list(Mic.2 = c(), Mic.3 = c(), Mic.4 = c(), Mic.5 = c())
        ),

        Reacting = list(
          genes = c("TMEM163", "HAMP", "SPP1"),
          states = list(Mic.6 = c("NAMPT", "RGS1", "ACSL1"), Mic.7 = c(), Mic.8 = c())
        ),

        Redox = list(
          genes = unique(c("AIF1", "APOC1", "RPL19", "RPLP2", "RPS11", "RPS6", "FTL", "FTH1", "APOE", "C1QA", "C1QB", "C1QC", "RPS19")),
          states = list(Mic.9 = c("P2RY12", "SORL1", "DOCK8"), Mic.10 = c("TMEM163", "SPP1", "HAMP"))
        ),

        Disease = list(
          genes = c("TPRG1", "PTPRG", "CPM", "DIRC3", "GPNMB", "APOE", "MYO1E", "TREM2"),
          states = list(Mic.12 = c("PTPRG", "TPRG1", "OLR1", "ADAMTS17", "BIN1", "ADAM10"),
                        Mic.13 = c("MS4A6A", "CD163", "MS4A4A", "MRC1", "SELENOP", "MCTP1", "HLA-DRA"))
        ),

        InterferonResponse = list(
          genes = c("IFI44L", "IFI6", "IFIT3", "PTPRG"),
          states = list(Mic.14 = c("IFI44L", "IFI6", "IFIT3", "PTPRG"))
        ),

        Other = list(
          states = list(Mic.11 = c("A2M"), Mic.15 = c("NLRP3", "IL1B", "IL18"), Mic.16 = c())
        ),

        Monocytes = list(
          genes = c("VCAN", "S100A4", "CD63")
        ),

        Macrophages = list(
          genes = c("MRC1", "CD63", "CD74", "CDK1")
        ),

        Homeostatic = list(
          genes = unique(c("FRMD4A", "P2RY12", "SYNDIG1", "CX3CR1", "SH3RF3", "TMEM163", "SOX5", "JUN", "SIGLECH", "SELPLG", "BTG2", "FOS", "TMEM119", "EGR1", "KLF2"))
        ),

        Inflammatory = list(
          genes = unique(c("GPNMB", "CD83", "IL1B", "NFKB1", "CCL3", "CCL4", "SLC2A3", "SRGN", "GPR183", "HSPA1A", "APOE", "IGF1", "APOC4", "CXCL2", "LYZ"))
        ),

        Stress = list(
          genes = unique(c("FOS", "JUNB", "HSPA1A", "HSPA1B", "HSP90AA1", "SLC2A3", "SRGN", "GPR183", "HIST1H2AC", "DDIT4"))
        ),

        Proliferating = list(
          genes = unique(c("ARHGAP11B", "TOP2A", "HMMR", "BRIP1", "MKI67", "FANCI"))
        ),

        DAM1_down = list(
          genes = c("P2RY12", "P2RY13", "SERINC3", "CX3CR1", "TGFBR1", "TMEM119", "TXNIP", "GLUL")
        ),

        DAM1 = list(
          genes = c("APOE", "B2M", "CSTB", "TYROBP", "TIMP2", "HLA-A", "HLA-B", "HLA-C", "FTH1", "LYZ", "CTSB", "CTSD")
        ),

        DAM2 = list(
          genes = c("ANK1", "SPP1", "AXL", "CSF1", "CST7", "CD9", "CADM1", "CLEC7A", "CCL6", "ITGAX", "CD63", "CD68", "CTSA", "LPL", "GUSB", "SERPINE2", "CTSZ", "CD52", "CTSL", "HIF1A")
        ),

        AD1.c7 = list(
          genes = c("DSCAM", "IPCEF1", "SOCS6", "ADAMTS17")
        ),

        AD1.c9 = list(
          genes = c("ITGAX", "MYO1E", "GLDN", "DTNA", "SPATS2L", "TPRG1", "PTPRG", "LPL")
        ),

        AD1.c10 = list(
          genes = c("P2RY12", "CX3CR1", "STARD13", "CPM", "EYA2", "GPNMB", "APOE")
        ),

        AD2 = list(
          genes = c("GRID2", "ADGRB3", "DPP10", "DDX17", "NAV2", "FOXP2")
        ),

        Dividing_Microglia = list(
          genes = c("STMN1", "TUBA1B", "PCLAF", "H2AZ1", "HMGB2", "BIRC5", "UBE2C", "TUBB5", "RAN", "CKS1B")
        ),

        Migrating_Microglia = list(
          genes = c("LGALS3", "NPY", "VIM", "PLIN2", "ADAM8", "LGALS1", "FABP5", "SPP1", "GPNMB", "CCL9")
        ),

        Interferon_Myeloid = list(
          genes = c("ISG15", "IFITM3", "IFIT3", "IRF7", "RSAD2", "IFI204", "RTP4", "BST2", "CCL5", "CXCL10")
        )
      )
    )
  }
  else {
    atlas <- list(
      # Microglia grouped subsets with marker genes and genes of interest
      microglia = list(
        General = list(genes = c("C3")),

        Profilerative = list(
          # states = list(Mic.1 = c("ARHGAP11B", "TOP2A", "HMMR"))
        ),

        Surveilling = list(
          genes = c("P2RY12", "CX3CR1", "FRMD4A", "SYNDIG1"),
          # states = list(Mic.2 = c(), Mic.3 = c(), Mic.4 = c(), Mic.5 = c())
        ),

        Reacting = list(
          genes = c("TMEM163", "HAMP", "SPP1"),
          # states = list(Mic.6 = c("NAMPT", "RGS1", "ACSL1"), Mic.7 = c(), Mic.8 = c())
        ),

        Redox = list(
          genes = unique(c("AIF1", "APOC1", "RPL19", "RPLP2", "RPS11", "RPS6", "FTL", "FTH1", "APOE", "C1QA", "C1QB", "C1QC", "RPS19")),
          # states = list(Mic.9 = c("P2RY12", "SORL1", "DOCK8"), Mic.10 = c("TMEM163", "SPP1", "HAMP"))
        ),

        Disease = list(
          genes = c("TPRG1", "PTPRG", "CPM", "DIRC3", "GPNMB", "APOE", "MYO1E", "TREM2"),
#           states = list(Mic.12 = c("PTPRG", "TPRG1", "OLR1", "ADAMTS17", "BIN1", "ADAM10"),
                        Mic.13 = c("MS4A6A", "CD163", "MS4A4A", "MRC1", "SELENOP", "MCTP1", "HLA-DRA"))
        ),

        InterferonResponse = list(
          genes = c("IFI44L", "IFI6", "IFIT3", "PTPRG"),
#           states = list(Mic.14 = c("IFI44L", "IFI6", "IFIT3", "PTPRG"))
        ),

        Other = list(
          # states = list(Mic.11 = c("A2M"), Mic.15 = c("NLRP3", "IL1B", "IL18"), Mic.16 = c())
        ),

        Monocytes = list(
          genes = c("VCAN", "S100A4", "CD63")
        ),

        Macrophages = list(
          genes = c("MRC1", "CD63", "CD74", "CDK1")
        ),

        Homeostatic = list(
          genes = unique(c("FRMD4A", "P2RY12", "SYNDIG1", "CX3CR1", "SH3RF3", "TMEM163", "SOX5", "JUN", "SIGLECH", "SELPLG", "BTG2", "FOS", "TMEM119", "EGR1", "KLF2"))
        ),

        Inflammatory = list(
          genes = unique(c("GPNMB", "CD83", "IL1B", "NFKB1", "CCL3", "CCL4", "SLC2A3", "SRGN", "GPR183", "HSPA1A", "APOE", "IGF1", "APOC4", "CXCL2", "LYZ"))
        ),

        Stress = list(
          genes = unique(c("FOS", "JUNB", "HSPA1A", "HSPA1B", "HSP90AA1", "SLC2A3", "SRGN", "GPR183", "HIST1H2AC", "DDIT4"))
        ),

        Proliferating = list(
          genes = unique(c("ARHGAP11B", "TOP2A", "HMMR", "BRIP1", "MKI67", "FANCI"))
        ),

        DAM1_down = list(
          genes = c("P2RY12", "P2RY13", "SERINC3", "CX3CR1", "TGFBR1", "TMEM119", "TXNIP", "GLUL")
        ),

        DAM1 = list(
          genes = c("APOE", "B2M", "CSTB", "TYROBP", "TIMP2", "HLA-A", "HLA-B", "HLA-C", "FTH1", "LYZ", "CTSB", "CTSD")
        ),

        DAM2 = list(
          genes = c("ANK1", "SPP1", "AXL", "CSF1", "CST7", "CD9", "CADM1", "CLEC7A", "CCL6", "ITGAX", "CD63", "CD68", "CTSA", "LPL", "GUSB", "SERPINE2", "CTSZ", "CD52", "CTSL", "HIF1A")
        ),

        AD1.c7 = list(
          genes = c("DSCAM", "IPCEF1", "SOCS6", "ADAMTS17")
        ),

        AD1.c9 = list(
          genes = c("ITGAX", "MYO1E", "GLDN", "DTNA", "SPATS2L", "TPRG1", "PTPRG", "LPL")
        ),

        AD1.c10 = list(
          genes = c("P2RY12", "CX3CR1", "STARD13", "CPM", "EYA2", "GPNMB", "APOE")
        ),

        AD2 = list(
          genes = c("GRID2", "ADGRB3", "DPP10", "DDX17", "NAV2", "FOXP2")
        ),

        Dividing_Microglia = list(
          genes = c("STMN1", "TUBA1B", "PCLAF", "H2AZ1", "HMGB2", "BIRC5", "UBE2C", "TUBB5", "RAN", "CKS1B")
        ),

        Migrating_Microglia = list(
          genes = c("LGALS3", "NPY", "VIM", "PLIN2", "ADAM8", "LGALS1", "FABP5", "SPP1", "GPNMB", "CCL9")
        ),

        Interferon_Myeloid = list(
          genes = c("ISG15", "IFITM3", "IFIT3", "IRF7", "RSAD2", "IFI204", "RTP4", "BST2", "CCL5", "CXCL10")
        )
      )
  }
  return(atlas)
}
