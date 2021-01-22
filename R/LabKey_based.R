
#' @title download_SeuratHashingCalls
#'
#' @description Download the hashing calls for a seurat object based on LoupeId
#' @param loupeId LoupID e.g. 169823
#' @param outPath A path to save the hash files
#' @param lkBaseUrl setting of labkey
#' @param lkDefaultFolder setting of labkey
#' @return saves hashing calls in outPath as .3.hashing.calls.txt
#' @export
download_SeuratHashingCalls <- function(loupeId, outPath="./", 
                                       lkBaseUrl = "https://prime-seq.ohsu.edu",
                                       lkDefaultFolder = "/Labs/Bimber") {
  
  # loupeId = 169823
  rows <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="readset",
    colFilter=makeFilter(c("rowid", "EQUAL", loupeId)),
    containerFilter=NULL,
    colNameOpt="rname"
  )
  
  if (nrow(rows) != 1) {
    stop('Loupe output not found')
  }
  
  readset <- unique(rows$readset)
  
  rows <- labkey.selectRows(
    baseUrl=lkBaseUrl,
    folderPath=lkDefaultFolder,
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="rowid,",
    colFilter=makeFilter(c("readset", "EQUAL", readset), c("category", "EQUAL", "Seurat Data")),
    containerFilter=NULL,
    colNameOpt="rname"
  )
  
  if (nrow(rows) == 0) {
    stop('Seurat output not found')
  }
  
  #Take latest one
  seuratId <- max(rows$rowid)
  
  hashingFile <- paste0(outPath, '/hashing', '.', readset, '.', seuratId, '.calls.txt')
  DownloadOutputFile(outputFileId = seuratId, outFile = hashingFile, overwrite = T, pathTranslator = function(x){
    return(gsub(x, pattern = '.seurat.rds', replacement = '.3.hashing.calls.txt'))
  })
}
