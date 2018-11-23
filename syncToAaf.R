#############################
## Cas Retel
## cas.retel@eawag.ch
## 2017.04.05
#############################
## Takes a .sync file as input and outputs a file of similar
## dimensions, substituting allelic counts at every position
## for the proportion equal to reference (=column 3 in .sync). 
## usage example: 
## "Rscript ~/syncToAaf.R D00-99_nc64a.bwaal.flt6.2of5.sync \
## FALSE D00-99_nc64a.bwaal.2of5.aaf"
#############################

require(magrittr)
source("/cluster/project/gdc2/special/shared/p368/scripts/cas/calcRefAF.R")

args <- commandArgs(TRUE)

## input filename
fname_in <- args[1]
## does first row of input filename consist of column names?
header_present <- as.logical(args[2])
## output filename (often identical to input with tag ".aaf" instead of ".sync")
fname_out <- args[3]

## read in file
in_sync <- read.table(
  fname_in, sep="\t", header=header_present
)

print(paste("###### FILE ", fname_in, " READ IN, ", ncol(in_sync)-3, 
            " SAMPLES PRESENT ######", sep=""))

out_sync <- in_sync[, 1:3]

for(col in 4:ncol(in_sync)){
  print(
    paste("###### CALCULATING ANCESTRAL ALLELE FREQUENCY OF SAMPLE NUMBER ", 
          col-3, " ######", sep=""))
  out_sync <- cbind(
    out_sync, calcRefAF(in_sync[, col], ref=in_sync[, 3])[3, ]
  )
}

print(paste0("###### syncToAaf.R DONE, WRITING TO ", fname_out))
write.table(
  out_sync, file=fname_out, quote=F, sep="\t", 
  row.names=F, col.names=F
)
rm(in_sync, out_sync)