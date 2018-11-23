#######################################
## calcAverageAF(): 
## from a matrix with numeric values, return colMeans(input)
## (from a vector with numeric values, return input)
#######################################
## checkFlankingRegs(): 
## if input position {pos} is within {flank.reg} base pairs of any 
## position given in {dels}, return TRUE. 
#######################################
## varAnc(): returns T if first element is > cutoff
#######################################
## varFirstObs(): returns T if first non-NA element is > cutoff
#######################################
## hasNeighbouringPos(): 
## Calculates if any genomic positions specified in input {chrom} and {pos}
## are within {flank.reg} bp of each other. 
#######################################
calcAverageAF <- function(afmat){
  ## from a vector with numeric values, return input
  ## from a matrix with numeric values, return colMeans(input)
  if(is.null(dim(afmat))){
    out <- afmat
  }else{ 
    out <- colMeans(afmat, na.rm=TRUE)
  }
  return(out)
}

checkFlankingRegs <- function(pos, dels, flank.reg=10){
  ## if input position {pos} is within {flank.reg} base pairs of any 
  ## position given in {dels}, return TRUE. 
  
  ## pos = one row of a read in .sync file, with element one scaffold name
  ##   and element two position
  ## dels = data frame of a read in .sync file, with column one scaffold names
  ##   and column two positions
  ## flank.reg = number of bp adjacent to a del to output TRUE for
  dels <- data.frame(dels[, 1:2])
  dels[, 2] <- as.numeric(as.character(dels[, 2]))
  
  ## becomes TRUE if scaffold names match and difference in positions <= flank.reg: 
  return(any(
    (dels[, 1]==as.character(pos[1])) & (abs(as.numeric(pos[2])-dels[, 2]) <= flank.reg)
  ))
}


varAnc <- function(daf_row, cutoff=0.1, also_derived=TRUE){
  ## returns T if first element is > cutoff
  ## return F otherwise (first element is <= cutoff or NA)
  ## also_derived: if set to TRUE, function returns TRUE if 
  ##    min(daf_row[1], 1 - daf_row[1]) > cutoff
  if(!also_derived) {
    return(!(is.na(daf_row[1]) | daf_row[1] < cutoff))
  }else{
    return(!(is.na(daf_row[1]) | (min(daf_row[1], 1 - daf_row[1]) < cutoff)))
  }
}

varFirstObs <- function(daf_row, cutoff=0.1, also_derived=TRUE){
  ## returns T if first non-NA element is > cutoff
  ## return F if first non-NA element is <= cutoff
  ## also_derived: if set to TRUE, function returns TRUE if 
  ##    min(daf_row[!is.na(daf_row)[1]], 1 - daf_row[!is.na(daf_row)[1]]) > 
  ##      cutoff
  
  if(mean(is.na(daf_row)) == 1){
    stop("No non-NA values observed")
  }
  
  ## find first observed element 
  daf_fo <- daf_row[which(!is.na(daf_row))[1]]
  if(!also_derived){
    return(daf_fo > cutoff)
  }else{
    return(min(daf_fo, 1 - daf_fo) > cutoff)
  }
}

hasNeighbouringPos <- function(chrom, pos, flank.reg=1000){
  ## Calculates if any genomic positions specified in input {chrom} and {pos}
  ## are within {flank.reg} bp of each other. Returns a logical vector of
  ## same length as {chrom} and {pos}, with TRUE if element has a neighbour
  ## and FALSE if it doesn't. 
  ## ! assumes {chrom} and {pos} are ordered !
  
  ## chrom = character vector containing chromosome names
  ## pos = character vector of same length containing positions
  ## flank.reg = maximum number of base pairs between positions that are 
  ##   still considered neighbours
  if(!all(c(is.character(chrom), is.numeric(pos)))){
    stop("{chrom} or {length} is of the wrong class")
  }
  if(length(chrom) < 2){
    stop("you can't do this with only one value")
  }
  if(length(chrom) != length(pos)){
    stop("{chrom} and {length} have different lengths")
  }
  if(!is.numeric(flank.reg) | (flank.reg < 1)){
    stop("invalid {flank.reg} specification")
  }
  
  
  ## determine number of entries
  p <- length(chrom)
  ## for elements 1 until {p-1}, check if they neighbour the next element
  tmp1 <- (chrom[1:(p - 1)] == chrom[2:p]) & 
    ((pos[2:p] - pos[1:(p - 1)]) <= flank.reg )
  ## if the i'th element of tmp1 is true, set the {i+1}'th element the 
  ## output vector to true
  out <- c(tmp1, FALSE) | c(FALSE, tmp1)
  return(out)
}


# chrt <- paste("chrom_", rep(c(1, 4, 7), c(4, 2, 6)))
# post <- c(4, 18, 26, 90, 92, 120, 3, 4, 5, 6, 16, 27)
# hasNeighbouringPos(chrt, post, flank.reg=10)
