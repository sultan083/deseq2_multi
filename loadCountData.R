#' Load count files
#'
#' Load one raw count file
#'
#' @param target target \code{data.frame} of the project returned by \code{loadTargetFile()}
#' @param counts.file is a raw count file
#' @param header a logical value indicating whether the file contains the names of the variables as its first line
#' @param skip number of lines of the data file to skip before beginning to read data
#' @param featuresToRemove vector of feature Ids (or character string common to feature Ids) to remove from the counts
#' @return The \code{matrix} of raw counts with row names corresponding to the feature Ids and column names to the sample names as provided in the first column of the target.
#' @details If \code{featuresToRemove} is equal to \code{"rRNA"}, all the features containing the character string "rRNA" will be removed from the counts.
#' @author Marie-Agnes Dillies and Hugo Varet

loadCountData <- function(target, counts.file, header=TRUE, skip=0, featuresToRemove){

  labels <- as.character(target[,1])

  rawCounts <- read.table(counts.file, sep="\t", quote="\"", header=header, skip=skip)
  
  colnames(rawCounts) <- c("Id", labels)
  
  if (any(duplicated(rawCounts$Id))) stop("Duplicated feature names in ", counts.file)
  
  for (i in labels) {
  cat(i,": ",length(rawCounts[,i])," rows and ",sum(rawCounts[,i]==0)," null count(s)\n",sep="")
  }
  
  rawCounts[is.na(rawCounts)] <- 0
  counts <- as.matrix(rawCounts[,-1])
  rownames(counts) <- rawCounts[,1]
  counts <- counts[order(rownames(counts)),]
  
  cat("\nFeatures removed:\n")
  
  for (f in setdiff(featuresToRemove,"")){
    match <- grep(f, rownames(counts))
    if (length(match)>0){
	  cat(rownames(counts)[match],sep="\n")
	  counts <- counts[-match,]
	}
  }

  cat("\nTop of the counts matrix:\n")
  print(head(counts))
  cat("\nBottom of the counts matrix:\n")
  print(tail(counts))
  return(counts)
}
