#!/usr/bin/Rscript

################################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Upendra Devisetty
### June 16th, 2016
################################################################################

# Load libraries
library(DESeq2)
library(genefilter)
library(devtools)
library(SARTools)
library(getopt)
library(knitr)
library(RColorBrewer)
library(genefilter)


args<-commandArgs(TRUE)

options<-matrix(c('project',  'pn', 1,  "character",    	# project name
                  'author', 'au', 1,  "character",      	# author of the statistical analysis/report  
                  'Dir',  'r',  2,  "character",        	# path to the directory containing raw counts files
		              'rawCounts', 'rc', 2, "character",		  # path to the combined raw count file
                  'OutDir',  'w',  1,  "character",     	# path to the output file 
                  'target', 't',  1,  "character",      	# path to the design/target file
                  'features', 'fe', 2,  "character",    	# names of the features to be removed (specific HTSeq-count information and rRNA for example)
                  'varInt', 'v',  2,  "character",      	# factor of interest
                  'condRef',  'c',  2,  "character",    	# reference biological condition
                  'batch',  'b',  2,  "character",      	# blocking factor: NULL (default) or "batch" for example
		              'locfunc', 'l', 2, "character", 		    # "median" (default) or "shorth" to estimate the size factors
		              'fitType', 'f',	2, "character", 	      # mean-variance relationship: "parametric" (default) or "local"
		              'cooksCutoff', 'cc',	2, "logical", 		# outliers detection threshold (TRUE to let DESeq2 choosing it or FALSE to disable the outliers detection)
		              'independentFiltering', 'if', 2, "logical", 	# TRUE or FALSE to perform the independent filtering or not
		              'typeTrans', 'tt', 2, "character",		# transformation method for PCA/clustering with DESeq2: "VST" or "rlog"
                  'alpha',  'a',  2,  "double",         	# threshold of statistical significance
                  'pAdjust',  'p',  2,  "character",    	# p-value adjustment method: "BH" (default) or "BY"
                  'colors', 'co', 2,  "character",      	# vector of colors of each biological condition on the plots
                  'help',   'h',    0,      "logical"),
                            ncol=4,byrow=TRUE)

ret.opts<-getopt(options,args)

if ( !is.null(ret.opts$help) ) {
  cat(getopt(options, usage=TRUE));
  q(status=1);
}

projectName <- ret.opts$project
author  <-  ret.opts$author
targetFile <- ret.opts$target
rawDir <- ret.opts$Dir
rawCounts <- ret.opts$rawCounts
OutDir <- ret.opts$OutDir
featuresToRemove <- ret.opts$features
varInt  <- ret.opts$varInt
condRef <- ret.opts$condRef
batch <- ret.opts$batch
fitType <- ret.opts$filtType
cooksCutoff <- ret.opts$cooksCutoff
independentFiltering <- ret.opts$independentFiltering
locfunc <- ret.opts$locfunc
alpha <- ret.opts$alpha
pAdjustMethod <- ret.opts$pAdjust
typeTrans <- ret.opts$typeTrans
col <- ret.opts$colors


# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
if (!is.null(ret.opts$Dir)) {
counts <- loadCountData(target=target, rawDir=rawDir, header=TRUE, skip=0, featuresToRemove=featuresToRemove)}
print(class(counts))


if (!is.null(ret.opts$rawCounts)) {
source("/Users/upendra_35/Documents/git.repos/deseq2_multi/loadCountData.R")
counts <- loadCountData(target=target, rawCounts=rawCounts, header=TRUE, skip=0, featuresToRemove=featuresToRemove)}


# description plots
source("/Users/upendra_35/Documents/git.repos/deseq2_multi/descriptionPlots.r")
majSequences <- descriptionPlots(counts=counts, n=3, group=target[,varInt], output.file=output.file, col=col)



