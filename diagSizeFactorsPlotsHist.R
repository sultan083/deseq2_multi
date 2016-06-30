#' Assess the estimations of the size factors
#'
#' Plots to assess the estimations of the size factors
#'
#' @param dds a \code{DESeqDataSet} object
#' @param outfile TRUE to export the figure in a png file
#' @param plots vector of plots to generate
#' @return Two files in the figures directory: diagSizeFactorsHist.png containing one histogram per sample and diagSizeFactorsTC.png for a plot of the size factors vs the total number of reads
#' @author Marie-Agnes Dillies and Hugo Varet

diagSizeFactorsPlotsHist <- function(dds, output.file="diagSizeFactorsHist.png"){
  # histograms
    nrow <- ceiling(sqrt(ncol(counts(dds))))
    ncol <- ceiling(ncol(counts(dds))/nrow)
    png(filename=output.file, width=1400*max(ncol,nrow), height=1400*min(ncol,nrow), res=300)
    par(mfrow=sort(c(nrow,ncol)))
    geomeans <- exp(rowMeans(log(counts(dds))))
    samples <- colnames(counts(dds))
    counts.trans <- log2(counts(dds)/geomeans)
    xmin <- min(counts.trans[is.finite(counts.trans)],na.rm=TRUE)
    xmax <- max(counts.trans[is.finite(counts.trans)],na.rm=TRUE)
    for (j in 1:ncol(dds)){
      hist(log2(counts(dds)[,j]/geomeans), nclass=100, xlab=expression(log[2] ~ (counts/geometric~mean)), las=1, xlim=c(xmin,xmax),
           main=paste0("Size factors diagnostic - Sample ",samples[j]),col="skyblue")
      abline(v = log2(sizeFactors(dds)[j]), col="red", lwd=1.5)
    }
    dev.off()

  }


