require(magrittr); require(ggplot2)
# setwd("~/Documents/HVInt/Chlorella/tmp2016/ppl_dir/")

# ## simulate .syncaf file
# nloc <- 1000
# set.seed(2303)
# sim_syncaf <- data.frame(
#   chrom=rep("PBCV1_scaffold_1", nloc), 
#   pos=8351:(8350+nloc), 
#   ref=sample(c("A", "C", "C", "T", "G", "G"), replace=T, size=nloc), 
#   D00=1-rexp(nloc, rate=100), 
#   D05=1-rexp(nloc, rate=100), 
#   D10=1-rexp(nloc, rate=100), 
#   D15=1-rexp(nloc, rate=100), 
#   D20=1-rexp(nloc, rate=100)
# )
# for(col in 4:ncol(sim_syncaf)){
#   change <- runif(nloc)<.005
#   sim_syncaf[change, col] <- runif(sum(change))
# }
# rm(col, change)
# write.table(sim_syncaf, file="sim.syncaf", quote=F, sep="\t", 
#             row.names=F, col.names=F)

findVariantRows <- function(fname_in, af_change=0.1, af_timepoints=numeric(), 
                            plot_afs=TRUE, plot_afs_indiv=FALSE, 
                            fname_out=character(), fname_plot=character()){
  ## from a tab-delimited file containing CHROM, POS and REF as first
  ## three fields, then ancestral allele frequencies of pooled
  ## sequencing samples, return the rows that where the absolute 
  ## difference between the minimum and maximum allele frequency 
  ## is >= af_change
  ## 
  ## - af_change: minimum observed difference between timepoints
  ## - plot_afs: create a plot? (may not be informative if the number of 
  ##   loci is large)
  ## - timepoints_afs: vector of time points (assumed)
  
  require(magrittr); require(ggplot2)
  ## load in object
  obj_aaf <- read.table(fname_in, sep="\t", header=F)
  
  ## to prevent warnings later, remove positions that containg only na values
  ind_allna <- apply(
    obj_aaf[, 4:ncol(obj_aaf)], 1, (function(x) all(is.na(x)))
  )
  cat(paste(
    "##### REMOVING ", sum(ind_allna), " POSITIONS BECAUSE THEY CONTAIN ONLY NA VALUES #####\n", sep=""
  ))
  obj_aaf <- obj_aaf[!ind_allna, ]
  
  ## if not specified, generate output filename
  if(length(fname_out)==0){
    fname_out <- gsub(
      ".aaf", x=fname_in, 
      replacement=paste0(".variant", af_change * 100, ".aaf")
    )
  }
  ## if not specified, generate plot filename
  if(plot_afs & length(fname_plot)==0){
    fname_plot <- paste( 
      strsplit(fname_in, split="[.]")[[1]][
        -length(strsplit(fname_in, split="[.]")[[1]])
      ], 
      ".variant", af_change * 100, ".pdf", sep=""
    )
  }
  ## if not specified, generate time points
  if(plot_afs & length(af_timepoints)==0){
    af_timepoints=0:(ncol(obj_aaf)-4)
  }
  
  ## find rows with af differences of at least af_change
  cat(paste(
    "##### LOOKING FOR ROWS WITH ALLELE FREQUENCY DIFFERENCES OF AT LEAST ", 
    af_change * 100, " PERCENT #####\n", sep=""
  ))
  rows_variant <- apply(
    obj_aaf[, 4:ncol(obj_aaf)], 1, 
    (function(x) abs(max(x[!is.na(x)]) - min(x[!is.na(x)])) >= af_change)
  )
  ## write these rows to new file
  cat(paste("##### WRITING THEM TO FILE ( TOTAL LINES =", sum(rows_variant), 
              ") #####\n"))
  write.table(
    obj_aaf[rows_variant, ], file=fname_out, quote=F, sep="\t", 
    row.names=F, col.names=F
  )
  
  ## plot allele frequency trajectories in one plot
  if(plot_afs){
    cat(paste("##### PLOTTING ANCESTRAL ALLELE FREQUENCY TRAJECTORIES #####\n"))
    pdf(fname_plot, width=9, height=4)
    par(mar=c(4,4,1,1))
    plot(c(0, max(af_timepoints)), c(0, 1), type='n', ylab="", xlab="")
    mtext(side=1, text="Time", 2.2)
    mtext(side=2, text="Ancestral allele frequency", 2.4)
    for(i in 1:sum(rows_variant)){
      lines(af_timepoints, obj_aaf[which(rows_variant)[i], 4:ncol(obj_aaf)], 
            col=rainbow(sum(rows_variant))[i])
    }
    dev.off()
  }
  
  ## plot allele frequency trajectories in individual plots
  if(plot_afs_indiv){
    cat(paste("##### AGAIN, NOW IN SEPARATE PLOTS #####\n"))
    for(i in 1:sum(rows_variant)){
      ## create filename
      fname_plot_i <- paste( 
        strsplit(fname_in, split="[.]")[[1]][
          -length(strsplit(fname_in, split="[.]")[[1]])
          ], 
        ".variant_traj", i, ".pdf", sep=""
      )
      ## create title ("main")
      main_plot_i <- paste(
        obj_aaf[which(rows_variant)[i], 1], "; position ", 
        obj_aaf[which(rows_variant)[i], 2], sep=""
      )
      ## plot
      pdf(fname_plot_i, width=9, height=4)
      par(mar=c(4,4,3,1))
      plot(c(0, max(af_timepoints)), c(0, 1), type='n', ylab="", xlab="")
      mtext(side=1, text="Time", 2.2)
      mtext(side=2, text="Ancestral allele frequency", 2.4)
      mtext(side=3, text=main_plot_i, 1, cex=1.4)
      lines(af_timepoints, obj_aaf[which(rows_variant)[i], 4:ncol(obj_aaf)], 
            col=rainbow(sum(rows_variant))[i])
      dev.off()
    }
  }
  
}

## process input arguments
args <- commandArgs(TRUE)
Afname_in <- args[1]
Aaf_change <- as.numeric(args[2])
Aaf_timepoints <- as.numeric(args[3])
Aplot_afs <- as.logical(args[4])
Aplot_afs_indiv <- as.logical(args[5])

findVariantRows(
  fname_in=Afname_in, af_change=Aaf_change, af_timepoints=Aaf_timepoints, 
  plot_afs=Aplot_afs, plot_afs_indiv=Aplot_afs_indiv
)

###################
## example use: 
## cd to location containing complete .aaf file
# Rscript ~/cluster/project/gdc2/special/shared/p368/scripts/cas/findVariantAlleles.R sim.aaf 0.1 0,8,12,14,21,27,29,35,51,83,99 TRUE FALSE


