#!/bin/R
######################################
## Filtering of ancestral allele frequencies: PBCV1
## Coevolutionary treatments II.2, III.2 and IV.2
## name: Cas Retel
## date: 2017.10.13
## e-mail: cas.retel@eawag.ch
######################################
## Reads were mapped with bwa, 
## filtered per position for minimum coverage of 10X and
## maximum coverage of { mean + 3*sd }, 
## and a few quality criteria
## This script starts with a .sync file containing information
## on all positions that change derived allele frequency by at least 5%. 
## This initial filtering criterium doesn't remove anything that is not
## removed at a later stage anyway, and makes the dataset handleable in
## terms of size. 


rm(list=ls())
pkgs <- list("ggplot2", "magrittr", "RColorBrewer")
sapply(pkgs, require, character.only=T)
source('~/Documents/Functions/gg_multiplot.R')
source('~/Documents/Functions/round_10e3.R')
source('~/Documents/HVInt/scripts/fs_syncCalculations.R')
source('~/Documents/HVInt/scripts/fs_syncFiltering.R')
source('~/Documents/HVInt/scripts/plotAfs.R')
source('~/Documents/HVInt/scripts/writeReadSync.R')

for(tag_treat in c("2dot2_", "3dot2_", "4dot2_")){
  ## set unique experimental tag that defines treatment
# tag_treat <- "2dot2_" ## to run per replicate

exp2017_dir <- "/Users/reteladmin/Documents/HVInt/Exp2017/Exp2017_30/"  
tag_in <- paste(tag_treat, "D12-99_pbcv1.bwaal.variant5", sep="")
tag_out <- paste(tag_treat, "D12-99_pbcv1", sep="")

## extract virus population sizes (observed and smoothed)
virhos_psdf <- read.csv(
  paste(exp2017_dir, "demog/", tag_treat, "hosvir0_logcounts_spline.csv", sep="")
)
## create a vector that has same scale for host and virus
virhos_psdf$log.ps.smooth.samescale <- virhos_psdf$log.ps.smooth - 3*(virhos_psdf$species=="PBCV1")
## create a scaled population size vector
virhos_psdf$ps.smooth.scaled <- c(
  with(subset(virhos_psdf, species=="PBCV1"), 
       10^(log.ps.smooth - max(log.ps.smooth))), 
  with(subset(virhos_psdf, species=="NC64A"), 
       10^(log.ps.smooth - max(log.ps.smooth)))
)

## set time points with genomic data
if(tag_treat == "2dot2_"){
  tps <- c(16, 21, 27, 29, 41, 51, 64, 69, 83, 99)
}else if(tag_treat == "3dot2_"){
  tps <- c(15, 21, 27, 29, 41, 51, 64, 69, 83, 99)
}else if(tag_treat == "4dot2_"){
  tps <- c(15, 21, 27, 29, 35, 51, 64, 70, 83, 99)
}
## read in sync file
sync_var5 <- readSync(
  paste(exp2017_dir, "genom/ppl_dir/", tag_in, ".sync", sep=""), tps
)
## sync files now have trailing whitespaces if position number 
## is shorter: 
sync_var5[, "pos"] <- as.integer(
  gsub("^\\s+|\\s+$", "", sync_var5[, "pos"])
)


#############################
## this is raw sequencing calls; set calls below minimum and maximum depth
## to NA (I already did this before subsetting to _var5, but decided
## download the raw frequency calls)
if(tag_treat == "2dot2_"){
  mindepth <- c(10, 64, 125, 79, 10, 96, 114, 94, 120, 88)
  maxdepth <- c(2040, 1939, 1877, 1921, 2022, 1905, 1885, 1906, 1882, 1914)
}else if(tag_treat == "3dot2_"){
  mindepth <- c(114, 46, 90, 17, 10, 24, 31, 10, 60, 60)
  maxdepth <- c(1885, 1950, 1908, 1982, 1994, 1976, 1968, 2018, 1940, 1938)
}else if(tag_treat == "4dot2_"){
  mindepth <- c(10, 10, 10, 10, 10, 10, 62, 114, 10, 10)
  maxdepth <- c(2002, 2091, 2174, 2242, 2068, 2228, 1939, 1886, 2160, 2170)
}

## calculate coverage, check which calls are outside confidence intervals
cov_var5 <- apply(sync_var5[, 4:ncol(sync_var5)], 2, sumAF)
apply(cov_var5, 2, mean)
## plot a histogram
## apply(cov_var5, 2, (function(x) hist(as.numeric(x), breaks=50)))

flt_var5 <- t(apply(
  cov_var5, 1, (function(x) (x <= mindepth) | (x > maxdepth))
))
sum(flt_var5, na.rm=TRUE); mean(flt_var5, na.rm=TRUE)
## 22% resp. 31 and 13% of calls are removed.. 
## average coverage per column is much lower than 1000X => I checked, and 
## did not mess up the downsampling. The variable positions are 
## often in more difficult to map regions. 

sync_var5[cbind(
  matrix(rep(FALSE, 3*nrow(sync_var5)), ncol=3), flt_var5
)] <- "0:0:0:0:0:0"


#############################
## remove time points when average coverage is too low, i.e.
## if average depth is under 10X
idx_tps_lowcov <- 
  apply(cov_var5, 2, (function(x) mean(x, na.rm=TRUE))) < 10
sum(!idx_tps_lowcov) ## should not remove any timepoints
tps_hicov <- tps[!idx_tps_lowcov]
sync_var5_tmp <- sync_var5[, c(T, T, T, !idx_tps_lowcov)]

## calculate derived allele frequencies
dafs_var5 <- syncToDafs(sync_var5_tmp)
## I have a few positions that are only called at the low-coverage
## time points. Need to remove these full-NA rows from my dataset
sync_var5_hicov <- sync_var5_tmp[
  apply(dafs_var5, 1, (function(x) !all(is.na(x)))), 
  ]
dafs_var5 <- dafs_var5[apply(dafs_var5, 1, (function(x) !all(is.na(x)))), ]
rm(sync_var5_tmp)


## dafs now include alleles fixed for non-reference at D12. For now, 
## the frequency that's lowest at first non-NA time point is taken as derived 
dafs_var5[apply(dafs_var5, 1, (function(x) x[!is.na(x)][1] > 0.5)), ] <- 
  1- dafs_var5[apply(dafs_var5, 1, (function(x) x[!is.na(x)][1] > 0.5)), ]


###################
###################
## calculate derived allele frequencies that change at least 25%
idx_var5 <- apply(
  dafs_var5, 1,
  (function(x) (max(x, na.rm=T) - min(x, na.rm=T)) >= .05)
)
sum(idx_var5, na.rm=T)

dafs_var25 <- dafs_var5[idx_var5, ]
sync_var25 <- sync_var5_hicov[idx_var5, ]
rownames(dafs_var25) <- rownames(sync_var25) <- 1:sum(idx_var5)
## 275 resp. 167 and 431



###################
## find del calls and adjacent positions
idx_var25_dels <- apply(
  sync_var25, 1, 
  (function(x) calcMajorDerivAF(x[4:length(x)], ref=x[3])[1, 2])
) == "d"
idx_var25_delregs <- apply(
  sync_var25, 1, (function(x) checkFlankingRegs(
    x, dels=sync_var25[idx_var25_dels, ], flank.reg=10
  ))
)
sum(!idx_var25_delregs); mean(!idx_var25_delregs)
## 163 left, resp. 133 and 308
# View(sync_var25[idx_var25_dels, ])

sync_var25_nodels <- sync_var25[!idx_var25_delregs, ]
dafs_var25_nodels <- dafs_var25[!idx_var25_delregs, ]



###################
## find positions that are already variable at day 12
idx_var25_varanc <- apply(
  dafs_var25_nodels, 1, 
  (function(x) varFirstObs(x, cutoff=0.01, also_derived=TRUE))
)
sum(!idx_var25_varanc); mean(!idx_var25_varanc)
## 34 left, resp. 39 and 81
sync_var25_denovo <- sync_var25_nodels[!idx_var25_varanc, ]
dafs_var25_denovo <- dafs_var25_nodels[!idx_var25_varanc, ]

###################
## INTERMEZZO: WRITE TO FILE TO ASSESS BETWEEN-REP REPRODUCIBILITY
# write.table(
#   sync_var25_nodels[idx_var25_varanc, ], quote=FALSE, row.names=FALSE,
#   col.names=FALSE, sep="\t",
#   file=paste0(exp2017_dir, "genom/ppl_dir/", tag_out, "_varanc.sync")
# )
## INTERMEZZO: WRITE TO FILE TO ASSESS BETWEEN-REP REPRODUCIBILITY
###################

###################
## check if derived allele frequency reaches the detection limit at 
## least twice. (only works like this in conjunction with  _varanc)
idx_var25_1x <- apply(
  dafs_var25_denovo, 1, (function(x) sum(x > .05, na.rm=TRUE) <= 1)
)
sum(!idx_var25_1x); mean(!idx_var25_1x)
## 12 resp. 11 and 16
sync_var25_2x <- sync_var25_denovo[!idx_var25_1x, ]
dafs_var25_2x <- dafs_var25_denovo[!idx_var25_1x, ]


###################
## _4nas: remove positions if more than 1 NA's
nas_perpos <- table(apply(dafs_var25_2x, 1, (function(x) sum(is.na(x)))))
plot(as.numeric(names(nas_perpos)), nas_perpos, type='p', pch=16, cex=2, 
     xlab="Number of NA values", ylab="Locus count", 
     main=paste0("NA values per locus (p = ", nrow(dafs_var25), ")"))
text(
  x=as.numeric(names(nas_perpos)), y=.5*max(nas_perpos), 
  labels=paste0("p = \n", nas_perpos)
)
nas_max <- 1 ## ? 
abline(v=.5 + nas_max, lty=2, col=2) ## remove what's right of the line

## 
idx_var25_namax <- apply(
  dafs_var25_2x, 1, (function(x) sum(is.na(x)) > nas_max)
)
mean(!idx_var25_namax); sum(!idx_var25_namax)
## Leaves 10 resp. 7 and 10
sync_var25_nonas <- sync_var25_2x[!idx_var25_namax, ]
dafs_var25_nonas <- dafs_var25_2x[!idx_var25_namax, ]

###################
## check temporal autocorrelation
# idx_var25_noacf <- apply(
#   dafs_var25_nonas, 1,
#   (function(x) abs(acf(x, na.action=na.pass, lag.max=1, plot=F)$acf[2]) < 0.2)
# )
# sum(!idx_var25_noacf, na.rm=T); mean(!idx_var25_noacf, na.rm=T)
# ## left, resp. and



###################
## plot datasets with two instead of three filtering criteria: 
# pdf(paste0(exp2017_dir, "genom/figs/fig_", tag_out, "_var5_filteringsteps.pdf"),
#     width=9, height=12)
# par(mfrow=c(4, 1))
# plotAfs(
#   dafs_var25_nodels, tps=tps_hicov, ylim_v=c(0, 1),
#   col_v=rainbow(nrow(dafs_var25))[as.numeric(rownames(dafs_var25_nodels))],
#   title_v=paste0(tag_out, "; _var5 and deletions removed")
# )
# plotAfs(
#   dafs_var25_denovo, tps=tps_hicov, ylim_v=c(0, 1),
#   col_v=rainbow(nrow(dafs_var25))[as.numeric(rownames(dafs_var25_denovo))],
#   title_v=paste0(tag_out, "; _var5, deletions removed and de novo")
# )
# plotAfs(
#   dafs_var25_2x, tps=tps_hicov, ylim_v=c(0, 1),
#   col_v=rainbow(nrow(dafs_var25))[as.numeric(rownames(dafs_var25_2x))],
#   title_v=paste0(tag_out, "; _var5, deletions removed, de novo and 2times")
# )
# plotAfs(
#   dafs_var25_nonas, tps=tps_hicov, ylim_v=c(0, 1),
#   col_v=rainbow(nrow(dafs_var25))[as.numeric(rownames(dafs_var25_nonas))],
#   title_v=paste0(tag_out, "; _var5, deletions removed, de novo, 2times and <2 NA's")
# )
# dev.off(); system(paste0("open ", exp2017_dir, "genom/figs/fig_", tag_out, "_var5_filteringsteps.pdf"))

## if there's only one position left, manually turn numeric into matrix
if(nrow(sync_var25_nonas) <= 1){
  sync_var25_nonas <- matrix(sync_var25_nonas, nrow=1)
  dafs_var25_nonas <- matrix(dafs_var25_nonas, nrow=1)
}

###################
## Check if these positions fall within 1000bp from each other ? 
## Not the case
###################


#############################
## Optional: plot with ggplot2, presentation quality
## generate plot-dataframe
dafs_plotdf <- data.frame(
  tp=rep(tps, ifelse(is.matrix(dafs_var25_nonas),
                     nrow(dafs_var25_nonas), 1)),
  daf=as.numeric(t(dafs_var25_nonas)),
  position=factor(rep(
    sync_var25_nonas$pos,
    each=ifelse(is.matrix(dafs_var25_nonas),
                ncol(dafs_var25_nonas),
                length(dafs_var25_nonas))
  ))
)

## create a colour palette factor that colours consequently across
## the three replicates (only 5 SNPS)
allele_pos_total <- factor(c(35117, 257553, 257902, 258323, 258479))
col_pal <- brewer.pal(n=length(allele_pos_total), name="Set1")[
  allele_pos_total %in% sync_var25_nonas$pos
  ]
set.seed(2303)
colyr <- colorRampPalette(colors=brewer.pal(n=9, "YlOrRd"), space="rgb")(nrow(dafs_var25_nonas))[sample(1:nrow(dafs_var25_nonas))]
# plot(1:nrow(dafs_var25_nonas), rep(0, nrow(dafs_var25_nonas)), pch=16, cex=1, col=colyr)


## generate the population-size plot
gplot_ps <- ggplot() + geom_line(
  data=virhos_psdf, aes(x=day, y=log.ps.smooth.samescale, colour=species),
  size=1.2
) + scale_colour_manual(values=c("#009E73", "#E69F00")) +
  scale_x_continuous(breaks=10*(0:12)) +
  scale_y_continuous(breaks=1*(2:10)) +
  coord_cartesian(x=c(0, 100), y=c(5.5, 9)) +
  theme_bw() + theme(
    panel.grid=element_blank(),
    # panel.border=element_rect(colour="#E69F00"),
    # legend.title=element_text(colour="#E69F00"),
    axis.text.y=element_text(size=20),#, colour="#E69F00"),
    axis.ticks.y=element_line(size=1),#, colour="#E69F00"),
    axis.text.x=element_text(size=20),#, colour="#E69F00"),
    axis.ticks.x=element_line(size=1)#, colour="#E69F00"),
  ) + labs(x="Time (days)", y="10log( Population size )") + 
  geom_vline(
    data=data.frame(tp=tps), aes(xintercept=tp), colour="red",
    alpha=.6, linetype="dashed", size=.4
  )

## generate the derived allele frequency plot
gplot_dafs <- ggplot() + geom_line(
  data=dafs_plotdf[!is.na(dafs_plotdf$daf), ],
  aes(x=tp, y=daf, colour=position), size=1, alpha=.95
) +
  geom_point(
    data=dafs_plotdf, aes(x=tp, y=daf, colour=position), size=3
  ) + scale_colour_manual(values=colyr) +
  scale_x_continuous(breaks=10*(0:12)) +
  scale_y_continuous(breaks=.2*(0:5)) +
  coord_cartesian(x=c(0, 100), y=c(0, 1)) +
  theme_bw() + theme(
    panel.grid=element_blank(),
    # panel.border=element_rect(colour="#E69F00"),
    # legend.title=element_text(colour="#E69F00"),
    axis.text.y=element_text(size=20),#, colour="#E69F00"),
    axis.ticks.y=element_line(size=1),#, colour="#E69F00"),
    axis.text.x=element_text(size=20),#, colour="#E69F00"),
    axis.ticks.x=element_line(size=1)#, colour="#E69F00"),
  ) + labs(x="Time (days)", y="Derived allele frequency") +
  geom_vline(
    data=data.frame(tp=tps), aes(xintercept=tp), colour="red",
    alpha=.6, linetype="dashed", size=.4
  )

pdf(paste(exp2017_dir, "genom/figs/fig_", tag_out, "_var5.pdf", sep=""),
    width=9, height=6)
# pdf(paste0("~/Documents/Study/AVW9/fig_", tag_out, "_var5.pdf"),
#     width=9, height=6)
gg_multiplot(
  gplot_ps + theme(plot.margin=unit(c(.2, .2, .2, .8), "cm")),
  gplot_dafs + theme(plot.margin=unit(c(.2, .2, .2, .2), "cm"))
  # layout=matrix(rep(1:2, c(2, 2)), ncol=1)
)
dev.off(); system(paste("open ", exp2017_dir, "genom/figs/fig_", tag_out, "_var5.pdf", sep=""))

}

###################
## write to file; needed for the next step
## add a column ac.12
write.table(
  sync_var25_nonas, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",
  file=paste(exp2017_dir, "genom/ppl_dir/", tag_out, "_var5_filtered.sync", sep="")
)
write.table(
  cbind(sync_var25_nonas[, 1:3], dafs_var25_nonas),
  quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t",
  file=paste(exp2017_dir, "genom/daf_dir/", tag_out, "_var5_filtered.daf", sep="")
)

}
