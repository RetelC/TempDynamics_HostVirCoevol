#############################
## Cas Retel
## cas.retel@eawag.ch
## 2017.04.05
#############################
## Takes a .sync file as input and outputs a file of similar
## dimensions, substituting allelic counts lower than min_depth
## and higher than max_depth 
## at every position with 0:0:0:0:0:0. 
## usage example if in.sync has four samples and no header: 
## "Rscript ~/syncFilterDepth.R in.sync FALSE 6 100,80,40,30 in.flt.sync"
#############################

require(magrittr)
source("/cluster/project/gdc2/special/shared/p368/scripts/cas/sumAF.R")

args <- commandArgs(TRUE)
vargs <- strsplit(args, ",")

#fname_in <- "D00-99_pbcv1.bwaal.sync"; header_present <- FALSE; min_depth=10; max_depth <- 0; fname_out <- "D00-99_pbcv1.bwaal.flt10.sync"

## input filename
fname_in <- vargs[[1]]
## does first row of input filename consist of column names?
header_present <- as.logical(vargs[[2]])
## minimum and maximum read depth; specify 0 not to use
## must be a numeric of length 1, or of the same length as the number
## of samples
min_depth <- vargs[[3]]; min_depth <- as.numeric(min_depth)
max_depth <- vargs[[4]]; max_depth <- as.numeric(max_depth)
## output filename
fname_out <- vargs[[5]]

## read in file
in_sync <- read.table(
  fname_in, sep="\t", header=header_present
)
print(paste("###### FILE ", fname_in, " READ IN, ", ncol(in_sync)-3, 
            " SAMPLES PRESENT ######", sep=""))

## check if min and max depths are specified correctly
if(!(length(min_depth) %in% c(1, ncol(in_sync)-3)) | 
   !(length(max_depth) %in% c(1, ncol(in_sync)-3))){
  stop("Input depth vectors do not all have the right length")
}
## elongate to vector if one value is provided
if(length(min_depth)==1){
  min_depth <- rep(min_depth, ncol(in_sync)-3)
}
if(length(max_depth)==1){
  max_depth <- rep(max_depth, ncol(in_sync)-3)
}

###################
## actual code; per sample in input .sync
for(col in 4:ncol(in_sync)){
  ## report progress
  print(
    paste(
      "###### FILTERING SAMPLE NUMBER ", col-3, 
      " FOR ALLELE FREQUENCIES BELOW ", min_depth[col-3], 
      " AND ABOVE ", max_depth[col-3], " #####", sep=""
    )
  )

  ## make sure the entries are no longer factors
  in_sync[, col] <- as.character(in_sync[, col])
  ## calculate sequencing depth per position
  tmp_af <- sumAF(in_sync[, col])

  ## determine which positions don't meet requirements
  ## (these if statements could be improved)
  ind_filtermax <- ind_filtermin <- logical(nrow(in_sync))
  if(max_depth[col-3] > 0){
    ind_filtermax <- (!is.na(tmp_af) & tmp_af > max_depth[col-3])
  }
  if(min_depth[col-3] > 0){
    ind_filtermin <- (!is.na(tmp_af) & tmp_af < min_depth[col-3])
  }
  ## replace those positions
  in_sync[(ind_filtermin | ind_filtermax), col] <- "0:0:0:0:0:0"
  rm(ind_filtermin, ind_filtermax, tmp_af)
}
## write to new file
write.table(
  in_sync, file=fname_out, quote=F, sep="\t", 
  row.names=F, col.names=F
)
rm(in_sync)