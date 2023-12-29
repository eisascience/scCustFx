
DevMode = T


if(!DevMode){
  #start with getting alingments from STAR etc, you need sample fastq_1,_2 a ref .fa and a gtf
  
  # featureCounts can be run seperately but this is the R version of doing the same thing
  
  #in the bam folders you can run this script example to compare NCBI vs CAT genes 
  
  # Set the current working directory to the desired path
  setwd(".")
  
  # List all files in the current directory
  files <- list.files(pattern = "Aligned.sortedByCoord.out.bam")
  
  # Set the paths to GTF files
  GTF_NCBI <- "/home/groups/ConradLab/common_resources/genomes/Macaca_mulatta/Macaca_mulatta.Mmul_10.107.gtf"
  GTF_CAT <- "/home/groups/ConradLab/common_resources/genomes/Macaca_mulatta/UNMC_rhesus_annotation_v7.6.8.gtf"
  
  # Specify output directory
  OUTPUT_DIR <- "."
  
  # Load the Rsubread package
  library(Rsubread)
  
  # Function to perform feature counting
  perform_feature_counts <- function(bam_file, gtf_file, prefix, saveOut = FALSE) {
    counts <- featureCounts(
      bam_file,
      annot.ext = gtf_file,
      isGTFAnnotationFile = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = "exon_id", #exon_id gene_id transcript_id
      requireBothEndsMapped = TRUE,
      useMetaFeatures = TRUE,
      nthreads = 8,
      strandSpecific = 0, #0 , 1 2
      isLongRead = FALSE,
      isPairedEnd = TRUE,
      countReadPairs = TRUE,
      verbose = TRUE,
      minMQS = 0
    )
    
    # Save the counts to a file if saveOut is TRUE
    if (saveOut) {
      write.table(counts, file = paste0(OUTPUT_DIR, "/", prefix, "_featureCounts.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
      cat("Counts saved to", paste0(prefix, "_featureCounts.txt"), "\n")
    }
    
    # Print a message indicating completion
    cat("featureCounts for", prefix, "completed.\n")
    
    # Return the counts
    return(counts)
  }
  
  # Filter files based on whether they contain "ncbi" or "cat" in their names
  ncbi_files <- grep("ncbi", files, value = TRUE)
  cat_files <- grep("cat", files, value = TRUE)
  
  # counts <- perform_feature_counts(ncbi_files[1], GTF_NCBI, gsub("_Aligned.sortedByCoord.out.bam", "", ncbi_files[1]), saveOut = F)
  # names(counts)
  
  
  ncbiLS = list()
  catLS = list()
  # Perform feature counting for ncbi files
  for (ncbi_file in ncbi_files) {
    ncbiLS[[ncbi_file]] = perform_feature_counts(ncbi_file, GTF_NCBI, gsub("_Aligned.sortedByCoord.out.bam", "", ncbi_file), saveOut = TRUE)
  }
  
  # Perform feature counting for cat files
  for (cat_file in cat_files) {
    catLS[[cat_file]] = perform_feature_counts(cat_file, GTF_CAT, gsub("_Aligned.sortedByCoord.out.bam", "", cat_file), saveOut = TRUE)
  }
  
  
  
  
  
  
  # aggregate_data_script.r
  
  # Load the required libraries
  library(dplyr)
  library(data.table)
  
  # Get the current path
  current_path <- getwd()
  
  # Initialize an empty list to store dataframes
  df_list <- list()
  
  # Find all files with the name pattern "counts_paired.txt.summary" in the current path
  files <- list.files(path = current_path, pattern = "counts_paired.txt.summary", full.names = TRUE)
  
  # Iterate through each file and read its content into a dataframe
  for (file in files) {
    # Read the file into a dataframe, assuming it's a tab-separated file
    df <- read.table(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    # Transpose the dataframe
    df <- as.data.frame(t(df))
    
    # Extract the relevant part of the dataframe
    df_part <- df[2, -1]
    
    # Add a column for the file name (last part of the path)
    df_part$File <- gsub("_counts_paired.txt.summary", "", basename(file))
    
    # Convert numerical columns to numeric
    df_part[, -ncol(df_part)] <- sapply(df_part[, -ncol(df_part)], as.numeric)
    
    # Append the dataframe to the list
    df_list[[length(df_list) + 1]] <- df_part
  }
  
  # Combine all dataframes into a single dataframe using rbindlist
  df_aggregated <- rbindlist(df_list, fill = TRUE)
  
  # Set column names
  colnames(df_aggregated) <- c("Assigned", "Unassigned_Unmapped", "Unassigned_Read_Type", 
                               "Unassigned_Singleton", "Unassigned_MappingQuality", 
                               "Unassigned_Chimera", "Unassigned_FragmentLength", 
                               "Unassigned_Duplicate", "Unassigned_MultiMapping", 
                               "Unassigned_Secondary", "Unassigned_NonSplit", 
                               "Unassigned_NoFeatures", "Unassigned_Overlapping_Length", 
                               "Unassigned_Ambiguity", "File")
  
  # Display the aggregated dataframe
  print(df_aggregated)
  
  # Save the aggregated dataframe as an RDS file in the current path
  saveRDS(df_aggregated, file.path(current_path, "aggregated_data.rds"))
}


