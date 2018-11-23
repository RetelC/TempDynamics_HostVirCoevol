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
source("/cluster/project/gdc2/special/shared/p368/scripts/cas/fs_syncCalculations.R")

args <- commandArgs(TRUE)

## input filename
fname_in <- args[1]
## output filename (often identical to input with tag ".daf" instead of ".sync")
fname_out <- args[2]

## read in file
in_sync <- read.table(
  fname_in, sep="\t", header=FALSE
)
colnames(in_sync) <- c("chrom", "pos", "ref", paste0("af", 1:(ncol(in_sync)-3)))

print(paste("###### FILE ", fname_in, " READ IN, ", ncol(in_sync)-3, 
            " SAMPLES PRESENT ######", sep=""))

## calculate derived allele frequencies
out_dafs <- syncToDafs(in_sync)

## write to new file
write.table(
  cbind(in_sync[, 1:3], out_dafs), file=fname_out, quote=F, sep="\t", 
  row.names=F, col.names=F
)
rm(in_sync, out_dafs)

