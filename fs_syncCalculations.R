#######################################
## sumAF(): from a 1-length character in .sync format (output by popoolation2), 
## calculate the sum allele frequency (i.e. coverage)
#######################################
## calcMajorDerivAF(): 
## from a vector of characters in .sync format (output by popoolation2), 
## calculate and output the allele frequency of the most abundant
## allele that is not equal to the reference
#######################################
## calcRefAF(): 
## from a vector in .sync format (output by popoolation2), 
## calculate the reference allele frequency
#######################################
## syncToDafs(): 
## from a .sync data.frame read into R using read.table(), calculate 
## derived or ancestral allele frequency. 
#######################################
## recodeDafs(): 
## takes an input matrix of (numeric) allele frequency values, and
## if the first element of the i'th row af_i,1 is higher than 0.5, changes that
## entire row's values to (1 - af_i ). 
#######################################
## syncPosToLogical()
## takes a correctly read in .sync file, a character vector scaffs
## and a numeric vector posits, and outputs a logical vector of length
## nrow(sync), with TRUEs if the sync entry matches a scaff-posit 
## combination. 
#######################################
## mergePositionsDaf()
## from a matrix (or data.frame) with chromosome names in column 1, 
## positions in column 2 and derived allele frequencies in the other columns, 
## merge the rows specified by {chroms} and {positions}. 
#######################################

sumAF <- function(x){
  ## from a 1-length character in .sync format (output by popoolation2), 
  ## formatted A-count:T-count:C-count:G-count:N-count:deletion-count
  ## calculate the sum allele frequency (i.e. coverage)
  
  require(magrittr)
  x <- as.character(x)
  ## split into six elements
  x <- strsplit(x, ":")[[1]] %>% as.numeric()
  if(length(x) != 6){
    stop("calcRefAF(): input doesn't contain six elements separated by colons")
  }
  return(sum(x))
}
sumAF <- Vectorize(sumAF, vectorize.args="x", USE=F)


# asd <- c("0:18:12:0:0:0", "0:2:5:0:0:0", "0:11:0:0:0:0", "0:0:12:0:0:0")
# asdf <- c("0:18:12:0:0:0", "0:2:5:0:0:0", NA, "0:11:0:0:0:0", "0:0:12:0:0:0")
# asdfg <- c("0:18:12:0:0:0", "0:5:0:7:0:0", NA, "0:0:0:1:0:0", "0:0:0:12:0:0")
# 
# calcMajorDerivAF(asd, ref="T")
# calcMajorDerivAF(asdf, ref="T")
# calcMajorDerivAF(asdfg, ref="T")

calcMajorDerivAF <- function(x, ref="N"){
  ## !!! DONT USE WITH APPLY !!!
  
  ## from a vector of characters in .sync format (output by popoolation2), 
  ## formatted A-count:T-count:C-count:G-count:N-count:deletion-count
  ## calculate and output the allele frequency of the most abundant
  ## allele that is not equal to the reference
  ## Most abundant = highest mean frequency across all samples. 
  ## outputs a matrix of as much columns as input vector has elements and 
  ## two rows: the first containing major derived allele frequencies and
  ## the second containing the major derived allele. 
  ## 23.10: added some checks of right input
  require(magrittr)
  if(!(ref %in% c("A", "T", "G", "C", "N"))){
    stop("Invalid reference value")
  }
  
  ## transform to character vector
  x <- sapply(x, as.character)
  
  ## generate output vector; fill in na elements
  xout <- rep(0, length(x))
  xout[is.na(x)] <- NA
  
  x <- x[!is.na(x)]
  
  ## split into six elements
  xl <- lapply(x, function(xel) as.numeric(strsplit(xel, ":")[[1]]))
  
  ## find which element represents reference
  bases_vect <- c("A", "T", "C", "G", "N", "d")
  ind_ref <- which(ref == bases_vect)
  
  ## calculate allele frequencies of all six variants
  x_afs <- sapply(xl, function(xel) xel / sum(xel))
  
  ## find which non-reference element has highest average frequency
  average_af <- rowMeans(x_afs, na.rm=T)
  average_af[ind_ref] <- 0
  ind_majorderiv <- which.max(average_af)
  
  ## fill in xout
  xout[!is.na(xout)] <- x_afs[ind_majorderiv, ]
  
  ## if sum of counts is 0, I want NA output, not NaN: 
  xout[is.nan(xout)] <- NA
  
  ## returns a data frame with columns major derived af and reference base
  ## this seems the most straightforward way (though not the most efficient)
  return(data.frame(
    maj.der.af=xout, 
    maj.der.base=rep(bases_vect[ind_majorderiv], length(xout)), 
    stringsAsFactors=FALSE
  ))
  ## !!! DONT USE WITH APPLY !!!
}
# calcRefAF <- Vectorize(calcRefAF, vectorize.args=c("x", "ref"))
## edit 2018.10.23: Just figured out why vectorizing gives problems: 
## a function that outputs a vector after vectorizing outputs a matrix of
## one column. This can't be done for the data frame of calcMajorDerivAF
## because one element is a character and the other is numeric. If I were
## to change it to a list of lists, this would work. 


calcRefAF <- function(x, ref="N"){
  ## from a .sync format (output by popoolation2), 
  ## formatted A-count:T-count:C-count:G-count:N-count:deletion-count
  ## calculate the reference allele frequency
  ## x = character in .sync-format ( e.g. "0:4:94:1:0:0" )
  ## ref = reference base at this genomic position
  
  ## outputs a numeric of length 3, containing
  ## - number of reference alleles
  ## - number of alternative alleles
  ## - reference allele frequency
  
  require(magrittr)
  x <- as.character(x)
  ## split into six elements; check input
  x <- strsplit(x, ":")[[1]] %>% as.numeric()
  if(length(x) != 6){
    stop("calcRefAF(): input doesn't contain six elements separated by colons")
  }
  if(!(ref %in% c("A", "T", "C", "G", "ind", "N"))){
    stop("calcRefAF(): invalid reference base")
  }
  
  ## if all entries are 0, raf can't be calculated
  if(all(x == 0)) return(t(c(0, 0, NA)))
  
  ## find which element represents reference
  ref_ind <- (ref=="A") + 2*(ref=="T") + 3*(ref=="C") + 4*(ref=="G") + 
    5*(ref=="N") + 6*(ref=="d")
  
  ## count reference, nonreference and reference allele frequency
  n_ref <- x[ref_ind]
  n_nonref <- sum(x[-ref_ind])
  return(t(c(n_ref, n_nonref, n_ref/sum(c(n_ref, n_nonref)))))
}
calcRefAF <- Vectorize(calcRefAF, vectorize.args=c("x", "ref"), USE=FALSE)

syncToDafs <- function(sync_in, ancestral=FALSE){
  ## from a .sync data.frame read into R using read.table(), calculate 
  ## derived or ancestral allele frequency. Returns a matrix of dimensions
  ## nrow(sync_in) * (ncol(sync_in) - 3)
  ## Basically a wrapper function around calcMajorDerivAF() and calcRefAF()
  ## ancestral = logical determining to output derived or ancestral frequency
  
  if(!all(colnames(sync_in)[1:3] == c("chrom", "pos", "ref"))){
    stop("Input sync file doesn't have the right column names")
  }
  if(ncol(sync_in) < 4){
    stop("Input sync file does not have enough columns")
  }
  if(nrow(sync_in) < 1){
    stop("Input sync file contains no loci")
  }
  ## create output vector
  afs_out <- numeric()
  ## per row, calculate derived allele frequency
  for(i in 1:nrow(sync_in)){
    if(!ancestral){
      afs_out <- rbind(
        afs_out, calcMajorDerivAF(
          sync_in[i, -(1:3)], ref=sync_in[i, 3]
        )[, 1])
    }else{
      afs_out <- rbind(
        afs_out, calcRefAF(
          sync_in[i, -(1:3)], ref=sync_in[i, 3]
        )[3, ]
      )
    }
  }
  rownames(afs_out) <- rownames(sync_in)
  colnames(afs_out) <- colnames(sync_in)[-(1:3)]
  return(afs_out)
}

recodeDafs <- function(dafs_in){
  ## takes an input matrix of (numeric) allele frequency values, and
  ## if the first element of the i'th row af_i,1 is higher than 0.5, changes that
  ## entire row's values to (1 - af_i ). Outputs a numeric matrix of same
  ## dimensions as daf_in
  if(!is.numeric(dafs_in)){
    stop("recodeDafs: Input matrix is not numeric")
  }
  if(any((dafs_in >= 1) | (dafs_in <= 0))){
    warning("recodeDafs: frequencies are not all on the interval ( 0, 1 )")
  }
  
  ## check if input is a matrix or just one row; if the latter, recode dafs
  ## for the one vector
  if(is.null(dim(dafs_in))){
    if(dafs_in[!is.na(dafs_in)][1] <= 0.5){
      return(dafs_in)
    }else{
      return(1 - dafs_in)
    }
  }
  
  ## if a matrix, change rows if first observed value is > 0.5
  dafs_in[apply(dafs_in, 1, (function(x) x[!is.na(x)][1] > 0.5)), ] <- 
    1 - dafs_in[apply(dafs_in, 1, (function(x) x[!is.na(x)][1] > 0.5)), ]
  return(dafs_in)
}

syncPosToLogical <- function(sync, scaffs, posits){
  ## takes a correctly read in .sync file, a character vector scaffs
  ## and a numeric vector posits, and outputs a logical vector of length
  ## nrow(sync), with TRUEs if the sync entry matches a scaff-posit 
  ## combination. 
  if(!(is.character(scaffs) & is.numeric(posits))){
    stop("syncPosToLogical(): scaffs and/or posits are of wrong class")
  }
  if(length(scaffs) != length(posits)){
    stop("syncPosToLogical(): scaffs and posits don't have the same length")
  }
  if(!(is.character(sync[, 1]) & is.numeric(sync[, 2]))){
    stop("syncPosToLogical(): input sync has wrong column classes")
  }
  if(!is.data.frame(sync)){
    stop("CR: input sync is not a data frame")
  }
  ## create output vector; initially contains only FALSES
  out <- logical(nrow(sync))
  ## walk over every scaff-posit combo
  for(i in 1:length(scaffs)){
    out[(scaffs[i] == sync[, 1]) & (posits[i] == sync[, 2])] <- TRUE
  }
  return(out)
}

mergePositionsDaf <- function(input, chroms, posits){
  ## from a matrix (or data.frame) with chromosome names in column 1, 
  ## positions in column 2 and derived allele frequencies in the other columns, 
  ## merge the rows specified by {chroms} and {positions}. 
  ## used to merge SNP calls that are very close to each other into one
  ## polymorphism
  ## chroms = list with every element a character of length 1, 
  ##   containing chromosome names of "clusters"
  ## positions = list with every element a numeric vector of length > 1, 
  ##   containing corresponding positions to merge into one row
  if(!(is.list(chroms) & is.list(posits))){
    stop("Either {chroms} or {posits} is not a list")
  }
  if(length(chroms) != length(posits)){
    stop("{chroms} and {posits} don't have equal length")
  }
  if(!all(unlist(chroms) %in% input[, 1])){
    stop("Not all {chroms} entries are present in first column of {input}")
  }
  for(i in 1:length(chroms)){
    if(sum(( input[, 1] == chroms[[i]] )  & ( input[, 2] %in% posits[[i]] )) != 
       length(posits[[i]])){
      stop(paste(
        "Chromosome ", chroms[[i]], 
        " (entry ", i, ") has a position entry that does not exist", sep=""
      ))
    }
  }
  
  ## create output object
  out <- input
  
  ## for every element of chroms and posits, replace the top row values 
  ## by the average daf values
  idx_remove <- idx_replace <- numeric()
  mean_afs <- matrix(nrow=0, ncol=ncol(input) - 2)
  
  for(i in 1:length(chroms)){
    ## elements to replace: first per chrom-posit combo
    idx_replace <- c(
      idx_replace, 
      which(( input[, 1] == chroms[[i]] ) & ( input[, 2] == posits[[i]][1] ))
    )
    ## elements to remove: everything except first per chrom-posit combo
    idx_remove <- c(
      idx_remove, 
      which(( input[, 1] == chroms[[i]] ) & ( input[, 2] %in% posits[[i]][-1] ))
    )
    ## calculate per-col average: 
    mean_afs <- rbind(
      mean_afs, 
      apply(
        input[
          ( input[, 1] == chroms[[i]] ) & ( input[, 2] %in% posits[[i]] ), 
          -c(1, 2)
          ], 2, (function(x) mean(x, na.rm=TRUE))
      )
    )
  }
  ## if whole column consists of NA's, mean(na.rm=TRUE) becomes NaN
  mean_afs[is.nan(mean_afs)] <- NA
  
  ## finally, remove and replace
  out[idx_replace, -c(1,2)] <- mean_afs
  out <- out[-idx_remove, ]
  return(out)
}
