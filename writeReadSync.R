
writeSync <- function(sync_in, fname, colnames_v=FALSE){
  ## calls write.table() with sensible default values. 
  ## check input
  if(!is.data.frame(sync_in)){
    stop("CR writeSync(): sync_in is not a data frame")
  }
  if(!is.character(fname)){
    stop("CR writeSync(): fname is not a character")
  }
  ## write file
  write.table(
    sync_in, quote=FALSE, row.names=FALSE, col.names=colnames_v, sep="\t", 
    file=fname
  )
}

readSync <- function(fname, tps_in=numeric(), header_v=FALSE){
  ## calls read.table() with sensible default values. Assigns column 
  ## names corresponding to the .sync file format (see Popoolation manual). 
  library(magrittr)
  ## check input
  if(!file.exists(fname)){
    stop("CR readSync(): fname does not exist")
  }
  if(!is.numeric(tps_in)){
    stop("CR readSync(): tps is not a numeric")
  }
  ## create time points if they're not provided
  if(length(tps_in)==0){
    print("CR readSync(): No time points provided!")
    ## obtain number of time points, i.e. number of fields in 
    ## first row of input file. 
    nt <- readLines(fname)[1] %>% 
      (function(ent) strsplit(ent, split="\t")[[1]]) %>% 
      (function(ent2) length(ent2) - 3)
    tps_in <- 1:nt
  }
  
  ## read file; ac stands for allele calls
  read.table(
    fname, header=header_v, stringsAsFactors=FALSE, sep="\t", 
    col.names=c("chrom", "pos", "ref", paste0("ac.", tps_in))
  )
}
