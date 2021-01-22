

#' @title get_MitoGeneSet
#'
#' @description returns a set of genes. 
#' @details https://www.genedx.com/wp-content/uploads/crm_docs/Mito-Gene-List.pdf
#' @param Name of of either c("mito1"); mito1 = www.genedx.com list.
#' @return a character vector of genes
#' @export
get_MitoGeneSet <- function(Name = NULL){
  
  if(is.null(Name)) stop("Name is NULL")
  
  if(Name == "mito1") {
    out = c("AARS", "AARS2", "ABCB11", "ABCB4", "ABCB7", "ABCD4", "ACAD9", "ACADM", 
            "ACADVL", "ACO2", "ACSF3", "ADCK3", "ADCK4", "AFG3L2", "AGK", "AGL", "AIFM1", 
            "ALAS2", "ALDOA", "ALDOB", "ALG1", "ALG11", "ALG13", "ALG2", "ALG3", "ALG6", "ALG9", 
            "AMACR", "APOPT1", "APTX", "ARG1", "ASL", "ASS1", "ATP5A1", "ATP5E", "ATP7B", "ATP8B1",
            "ATPAF2", "AUH", "B4GALT1", "BCKDHA", "BCKDHB", "BCS1L", "BOLA3", "C10orf2", "C12orf65", 
            "C19orf12", "CA5A*", "CARS2", "CHKB", "CISD2", "CLPB", "COA5", "COA6", "COASY", "COG4", 
            "COG5", "COG6", "COG7", "COG8", "COQ2", "COQ4", "COQ6", "COQ9", "COX10", "COX14", "COX15", 
            "COX20", "COX4I2", "COX6A1", "COX6B1", "COX7B", "CPS1", "CPT1A", "CPT2", "CYC1", "DARS", 
            "DARS2", "DBT", "DDHD1", "DDHD2", "DDOST", "DGUOK", "DLAT", "DLD", "DMGDH", "DNA2", "DNAJC19", 
            "DNM1L", "DNM2", "DOLK", "DPAGT1", "DPM1", "DPM3", "EARS2", "ECHS1", "ELAC2", "ENO3", 
            "ETFA", "ETFB", "ETFDH", "ETHE1", "FAH", "FARS2", "FASTKD2", "FBP1", "FBXL4", "FDX1L", 
            "FH", "FLAD1", "FOXRED1", "G6PC", "GAA", "GAMT", "GARS", "GATM", "GBE1", "GCDH", "GFER",
            "GFM1", "GFM2", "GLRX5", "GMPPA", "GSS", "GTPBP3", "GYG1", "GYG2", "GYS1", "GYS2", "HADHA", 
            "HADHB", "HARS2", "HCFC1", "HIBCH", "HLCS", "HMGCL", "HMGCS2", "HSD17B10", "HSPD1", "IARS2", 
            "IBA57", "ISCA2", "ISCU", "IVD", "LAMP2", "LARS", "LARS2", "LDHA", "LIAS", "LIPT1", "LMBRD1", 
            "LRPPRC", "LYRM4", "LYRM7", "MARS", "MARS2", "MCCC1", "MCCC2", "MCEE", "MFF", "MFN2", "MGAT2",
            "MGME1", "MICU1", "MLYCD", "MMAA", "MMAB", "MMACHC", "MMADHC", "MOGS", "MPC1", "MPDU1", "MPI", 
            "MPV17", "MRPL12", "MRPL3", "MRPL44", "MRPS16", "MRPS22", "MRPS7", "MTFMT", "MTO1", "MTPAP", 
            "MTR", "MTRR", "MUT", "NADK2", "NAGS", "NARS2", "NDUFA1", "NDUFA10", "NDUFA11", "NDUFA12", 
            "NDUFA2", "NDUFA4", "NDUFA9", "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", 
            "NDUFAF6 (C8ORF38)", "NDUFAF7", "NDUFB3", "NDUFB9", "NDUFS1", "NDUFS2", "NDUFS3",
            "NDUFS4", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NFS1", "NFU1", 
            "NGLY1", "NR2F1", "NUBPL", "OPA1", "OPA3", "OTC", "PARS2", "PC", "PCCA",
            "PCCB", "PDHA1", "PDHB", "PDHX", "PDP1", "PDSS1", "PDSS2", "PET100", 
            "PFKM", "PGAM2", "PGM1", "PHKA1", "PHKA2", "PHKB", "PHKG2", "PMM2", "PNPT1", 
            "POLG", "POLG2", "PRKAG2", "PRPS1", "PTRH2", "PUS1", "PYGM", "QARS", "RANBP2", 
            "RARS", "RARS2", "REEP1", "RFT1", "RMND1", "RRM2B", "SARS2", "SCO1", "SCO2", "SDHA*", 
            "SDHAF1", "SERAC1", "SFXN4", "SLC19A2", "SLC19A3", "SLC22A5", "SLC25A1", "SLC25A13", 
            "SLC25A15", "SLC25A19", "SLC25A20", "SLC25A22", "SLC25A3", "SLC25A38", "SLC25A4", 
            "SLC2A2", "SLC35A1", "SLC35A2", "SLC35C1", "SLC37A4", "SLC6A8*", "SLC7A7", "SPAST", 
            "SPG7", "SPTLC1", "SRD5A3", "SSR4", "STT3A", "STT3B", "STXBP1", "SUCLA2", "SUCLG1", 
            "SURF1", "TACO1", "TARS2", "TAZ", "TIMM8A", "TK2", "TMEM126A", "TMEM165", "TMEM70", 
            "TPK1", "TRIT1", "TRMU", "TRNT1", "TSFM", "TTC19", "TUFM", "TYMP", "UQCC2", "UQCC3", 
            "UQCRB", "UQCRC2", "UQCRQ", "VARS2", "WDR45", "WFS1", "YARS2")
  }
  
  return(out)
  
}

