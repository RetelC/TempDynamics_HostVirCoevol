round_10e3 <- function(inp, n_digits=1){
  ## rounds large numbers to the nearest smaller third power of 10; 
  ## to show it's behaviour, run:
  ## round_10e3(c(0.0003742, 0.004, 12, 73820188830), n_digits=1)
  
  
  ## store input
  inp <- as.numeric(unlist(inp))
  
  ## vector to store "exxxx"-appendix
  outp_power <- character(length(inp))
  ## check which power is highest
  poss_power <- (-33:33)*3
  max_power <- sapply(inp, function(x) max(poss_power[(10^poss_power) < x]))
  
  ## paste a "0" in front of max_power when it's smaller than 10
  outp_power <- paste("e", max_power, sep="")
  outp_power[max_power==0] <- ""
  outp_power[(max_power>0) & (max_power<10)] <- 
    paste("e0", max_power[(max_power>0) & (max_power<10)], sep="")
  outp_power[(max_power<0) & (max_power>-10)] <- 
    paste(
      "e-0", substr(
        max_power[(max_power<0) & (max_power>-10)], start=2, stop=3
      ), sep=""
    )
  # outp_power <- as.character(outp_power)
  
  ## calculate base
  outp_base <- round(inp / 10^max_power, digits=n_digits)
  
  return(paste(outp_base, outp_power, sep=""))
}
