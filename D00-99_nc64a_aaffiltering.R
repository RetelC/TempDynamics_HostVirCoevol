#!/bin/R
######################################
## Filtering of ancestral allele frequencies: NC64A
## Coevolutionary treatments II.2, III.2 and IV.2
## name: Cas Retel
## date: 2017.10.13
## e-mail: cas.retel@eawag.ch
######################################

## Reads were mapped with bwa, 
## filtered per position for minimum coverage of 10X and
## maximum coverage of { mean + 3*sd }, 
## and a few quality criteria
## Files with ancestral allele frequencies that go below
## 95% at at least one time-point were written to the 
## file read in below. 
#######################################
## Here, I filter ancestral allele frequencies for variation at D00 
## instead of at D12. To achieve this, the first thing to do is download
## information from all .variant25-positions for 4dot2_D00: 
#!/bin/bash
# awk 'FNR==NR{a[$1,$2]++;next}a[$1,$2]{print}' \
#   ppl_dir/2dot2_D00-99_nc64a.bwaal.variant25.sync \
#   ppl_dir/4dot2_D00-99_nc64a.bwaal.sync | cut -f1-4 > \
#   ppl_dir/D00_nc64a.bwaal.variant25.2dot2.sync
# awk 'FNR==NR{a[$1,$2]++;next}a[$1,$2]{print}' \
#   ppl_dir/3dot2_D00-99_nc64a.bwaal.variant25.sync \
#   ppl_dir/4dot2_D00-99_nc64a.bwaal.sync | cut -f1-4 > \
#   ppl_dir/D00_nc64a.bwaal.variant25.3dot2.sync
#######################################

#!/bin/R

rm(list=ls())
pkgs <- list("ggplot2", "magrittr", "RColorBrewer")
sapply(pkgs, require, character.only=T)
## files below can also be found at github.com/RetelC/TempDynamics_HostVirCoevolut
source('~/Documents/Functions/gg_multiplot.R')
source('~/Documents/Functions/round_10e3.R')
source('~/Documents/HVInt/scripts/fs_syncCalculations.R')
source('~/Documents/HVInt/scripts/fs_syncFiltering.R')
source('~/Documents/HVInt/scripts/plotAfs.R')
source('~/Documents/HVInt/scripts/writeReadSync.R')

for(tag_treat in c("2dot2_", "3dot2_", "4dot2_")){
## set unique experimental tag that defines treatment
# tag_treat <- "4dot2_"

exp2017_dir <- "/Users/reteladmin/Documents/HVInt/Exp2017/Exp2017_30/"
tag_in <- paste(tag_treat, "D00-99_nc64a.bwaal.variant25", sep="")
tag_out <- paste(tag_treat, "D00-99_nc64a", sep="")

## extract virus population sizes (observed and smoothed)
virhos_psdf <- read.csv(
  paste(exp2017_dir, "demog/", tag_treat, "hosvir0_logcounts_spline.csv", sep="")
)

## set time points with genomic data
if(tag_treat == "2dot2_"){
  tps0 <- c(7, 12, 14, 21, 27, 29, 41, 51, 64, 69, 83, 99)
}else if(tag_treat == "3dot2_"){
  tps0 <- c(8, 12, 15, 21, 27, 29, 41, 51, 64, 69, 83, 99)
}else if(tag_treat == "4dot2_"){
  tps0 <- c(0, 8, 12, 14, 21, 27, 29, 35, 51, 83, 99)
}
## read in sync file; I call these _var5, but actually they're
## {change 25% in aaf }
sync0_var5 <- readSync(
  paste(exp2017_dir, "genom/ppl_dir/", tag_in, ".sync", sep=""), tps0
)
## sync files now have trailing whitespaces if position number 
## is shorter: 
sync0_var5[, "pos"] <- as.integer(
  gsub("^\\s+|\\s+$", "", sync0_var5[, "pos"])
)

## for 2dot2 and 3dot2, add D00: 
if(tag_treat %in% c("2dot2_", "3dot2_")){
  sync_D00 <- readSync(paste0(
    exp2017_dir, "genom/ppl_dir/D00_nc64a.bwaal.variant25.", 
    substr(tag_treat, 1, 5), ".sync"), 0
  )
  if(identical(sync_D00$chrom, sync0_var5$chrom) & 
     identical(sync_D00$pos, sync0_var5$pos)){
    sync_var5 <- cbind(sync_D00, sync0_var5[, -(1:3)])
    tps <- c(0, tps0); # rm(sync_var5, tps0)
  }
}else{
  sync_var5 <- sync0_var5
  tps <- tps0; # rm(sync_var5, tps0)
}

#############################
## this is raw sequencing calls; set calls below minimum and maximum depth
## to NA (I already did this before subsetting to _var25, but decided to
## download the raw frequency calls). Max coverage (including D00 ,07 and 08): 
## 2dot2: 109, 181, 88, 41, 269, 55, 41, 109, 201, 130, 229, 360
## 3dot2: 176, 287, 70, 30, 728, 761, 48, 142, 392, 91, 146, 117
## 4dot2: 298, 444, 822, 332, 33, 709, 747, 33, 209, 361, 427

if(tag_treat == "2dot2_"){
  mindepth <- rep(10, length(tps))
  maxdepth <- c(298, 109, 181, 88, 41, 269, 55, 41, 109, 201, 130, 229, 360)
}else if(tag_treat == "3dot2_"){
  mindepth <- rep(10, length(tps))
  maxdepth <- c(298, 176, 287, 70, 30, 728, 761, 48, 142, 392, 91, 146, 117)
}else if(tag_treat == "4dot2_"){
  mindepth <- rep(10, length(tps))
  maxdepth <- c(298, 444, 822, 332, 33, 709, 747, 33, 209, 361, 427)
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
## D12: 45%, resp. 45% and 34% (!) of calls are removed
## D00: , resp 41% and 32 of calls are removed

## change sync entries outside coverage CI to "0:0:0:0:0:0"
sync_var5[cbind(
  matrix(rep(FALSE, 3*nrow(sync_var5)), ncol=3), flt_var5
)] <- "0:0:0:0:0:0"


#############################
## remove time points when average coverage is too low, i.e.
## if average depth is under 10X
# apply(cov_var25, 2, (function(x) mean(x, na.rm=TRUE)))
# apply(dafs_var25, 2, (function(x) mean(is.na(x))))
idx_tps_lowcov <- 
  apply(cov_var5, 2, (function(x) mean(x, na.rm=TRUE))) < 10
tps_hicov <- tps[!idx_tps_lowcov]
sum(!idx_tps_lowcov); sum(idx_tps_lowcov)
sync_var5_tmp <- sync_var5[, c(T, T, T, !idx_tps_lowcov)]
## D12: 8, resp. 8 and 7 time points left
## D00: 10 resp. 10 and 9 time points left

## calculate derived allele frequencies, recode to make ancestral 
## allele frequency 0
dafs_var5 <- syncToDafs(sync_var5_tmp)
# apply(aafs_var25_3nas, 2, (function(x) mean(is.na(x)))) < 0.25
## I have a few positions that are only called at the low-coverage
## time points. Need to remove these full-NA rows from my dataset
sync_var5_hicov <- sync_var5_tmp[
  apply(dafs_var5, 1, (function(x) !all(is.na(x)))), 
]
dafs_var5 <- dafs_var5[apply(dafs_var5, 1, (function(x) !all(is.na(x)))), ]
rm(sync_var5_tmp)

dafs_var5 <- recodeDafs(dafs_var5)

###################
###################
## calculate derived allele frequencies that change at least 25%
idx_var25 <- apply(
  dafs_var5, 1, 
  (function(x) (max(x, na.rm=T) - min(x, na.rm=T)) >= .25)
)
sum(idx_var25, na.rm=T)
dafs_var25 <- dafs_var5[idx_var25, ]
sync_var25 <- sync_var5_hicov[idx_var25, ]
## D12: 9296 resp. 8643 and 9342 loci change derived frequency by 25%
## D00: 10804 resp. 9623 and 13038 loci change derived frequency by 25%

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
## D12: 7944 resp. 7607 and 7958 after removing deletions
## D00: 9232 resp. 8403 and 11126 after removing deletions
sync_var25_nodels <- sync_var25[!idx_var25_delregs, ]
dafs_var25_nodels <- dafs_var25[!idx_var25_delregs, ]


###################
## find positions that are already variable at day 0
idx_var25_varanc <- apply(
  dafs_var25_nodels, 1, 
  (function(x) varFirstObs(x, cutoff=0.01, also_derived=TRUE))
)
sum(!idx_var25_varanc); mean(!idx_var25_varanc)
## D12: 1394 resp. 1409 and 1580 after removing variation in ancestor
## D00: 2326 resp. 2177 and 2262 after removing variation in ancestor
sync_var25_denovo <- sync_var25_nodels[!idx_var25_varanc, ]
dafs_var25_denovo <- dafs_var25_nodels[!idx_var25_varanc, ]

###################
## INTERMEZZO: WRITE TO FILE TO ASSESS BETWEEN-REP REPRODUCIBILITY
write.table(
  sync_var25_nodels[idx_var25_varanc, ], quote=FALSE, row.names=FALSE,
  col.names=FALSE, sep="\t",
  file=paste0(exp2017_dir, "genom/ppl_dir/", tag_out, "_varanc.sync")
)
## INTERMEZZO: WRITE TO FILE TO ASSESS BETWEEN-REP REPRODUCIBILITY
###################


###################
## check if derived allele frequency reaches the 
## detection limit at least twice. (only makes sense in conjunction with 
## _varanc)
idx_var25_1x <- apply(
  dafs_var25_denovo, 1, (function(x) sum(x > .25, na.rm=TRUE) <= 1)
)
sum(!idx_var25_1x); mean(!idx_var25_1x)
## D12: 89 resp. 129 and 281 reach detection limit twice
## D00: 418 resp. 470 and 646 reach detection limit twice
sync_var25_2x <- sync_var25_denovo[!idx_var25_1x, ]
dafs_var25_2x <- dafs_var25_denovo[!idx_var25_1x, ]


###################
## _4nas: remove positions if more than 4 NA's
nas_perpos <- table(apply(dafs_var25_2x, 1, (function(x) sum(is.na(x)))))
plot(as.numeric(names(nas_perpos)), nas_perpos, type='p', pch=16, cex=2, 
     xlab="Number of NA values", ylab="Locus count", 
     main=paste0("NA values per locus (p = ", nrow(dafs_var25), ")"))
text(
  x=as.numeric(names(nas_perpos)), y=.5*max(nas_perpos), 
  labels=paste0("p = \n", nas_perpos)
)
nas_max <- 4 ## or length(tps_hicov) - 3
abline(v=.5 + nas_max, lty=2, col=2) ## remove what's right of the line
## With two extra time points available, allowing for 4 NA values makes 
## sense for all three replicates

idx_3nas <- apply(
  dafs_var25_2x, 1, (function(x) sum(is.na(x)) > nas_max)
)
mean(!idx_3nas); sum(!idx_3nas)
## D12: Leaves 66 resp. 119 and 279 positions
## D00: Leaves 384 resp. 454 and 636 positions
dafs_var25_3nas <- dafs_var25_2x[!idx_3nas, ]
sync_var25_3nas <- sync_var25_2x[!idx_3nas, ]
# aafs_var25_3nas <- aafs_var25[!idx_3nas, ]


###################
## check temporal autocorrelation
# idx_var25_noacf <- apply(
#   dafs_var25_3nas, 1, 
#   (function(x) abs(acf(x, na.action=na.pass, lag.max=1, plot=F)$acf[2]) < 0.2)
# )
# sum(!idx_var25_noacf, na.rm=T); mean(!idx_var25_noacf, na.rm=T)
##  left, resp. 81 and 180


###################
## plot dataset after every filtering step: 
# pdf(paste0(exp2017_dir, "genom/figs/fig_", tag_out, "_var25_filteringsteps.pdf"),
#     width=9, height=12)
# par(mfrow=c(4, 1))
# plotAfs(
#   dafs_var25_nodels, tps=tps_hicov, xlim_v=c(0, 100), ylim_v=c(0, 1),
#   col_v=rainbow(nrow(dafs_var25))[as.numeric(rownames(dafs_var25_nodels))],
#   title_v=paste0(tag_out, "; _var25; deletions removed")
# )
# plotAfs(
#   dafs_var25_denovo, tps=tps_hicov, xlim_v=c(0, 100), ylim_v=c(0, 1),
#   col_v=rainbow(nrow(dafs_var25))[as.numeric(rownames(dafs_var25_denovo))],
#   title_v=paste0(tag_out, "; _var25; deletions removed and de novo")
# )
# plotAfs(
#   dafs_var25_2x, tps=tps_hicov, xlim_v=c(0, 100), ylim_v=c(0, 1),
#   col_v=rainbow(nrow(dafs_var25))[as.numeric(rownames(dafs_var25_2x))],
#   title_v=paste0(tag_out, "; _var25; deletions removed, de novo and 2times")
# )
# plotAfs(
#   dafs_var25_3nas, tps=tps_hicov, xlim_v=c(0, 100), ylim_v=c(0, 1),
#   col_v=rainbow(nrow(dafs_var25))[as.numeric(rownames(dafs_var25_3nas))],
#   title_v=paste0(tag_out, "; _var25; deletions removed, de novo, 2times and <= 3 NA's")
# )
# 
# dev.off(); system(paste0("open ", exp2017_dir, "genom/figs/fig_", tag_out, "_var25_filteringsteps.pdf"))

## if there's only one position left, manually turn numeric into matrix
if(nrow(sync_var25_3nas) <= 1){
  sync_var25_3nas <- matrix(sync_var25_3nas, nrow=1)
  dafs_var25_3nas <- matrix(dafs_var25_3nas, nrow=1)
}

###################
## Check if these positions fall within 1000bp from each other
nbs_var25_3nas <- with(
  sync_var25_3nas,
  hasNeighbouringPos(chrom=chrom, pos=pos, flank.reg=1000)
)
nbs_var25_3nas
## 2dot2: yes
## 3dot2: yes
## 4dot2: yes

## to check if these are sometimes one polymorphism, we need to plot them
## per chromosome
# set.seed(1990)
# pdf(paste0(exp2017_dir, "genom/tmp/fig_", tag_out, "_dafsperchrom.pdf"), 
#     width=7, height=5)
# for(i in 1:length(unique(sync_var25_3nas$chrom))){
#   idx_chr <- (sync_var25_3nas$chrom ==unique(sync_var25_3nas$chrom)[i]) & 
#     nbs_var25_3nas
#   if(!any(idx_chr)){ next }
#   tmp_col <- sample(rainbow(sum(idx_chr)))
#   plotAfs(
#     dafs_var25_3nas[idx_chr, ], tps=tps_hicov, col_v=tmp_col, 
#     title_v=paste0(
#       tag_out, ": SNPs on ", unique(sync_var25_3nas$chrom)[i], 
#       " and within 1000bp"
#     )
#   )
#   legend(
#     "topleft", legend=sync_var25_3nas[idx_chr, 2], fill=tmp_col, bty='n', 
#     ncol=max(c(1, round(sum(idx_chr) / 4, 0)))
#   )
# }
# dev.off(); system(paste0("open ", exp2017_dir, "genom/tmp/fig_", tag_out, "_dafsperchrom.pdf"))
# write.table(
#   x=cbind(sync_var25_3nas[, 1:3], dafs_var25_3nas),
#   file=paste("~/Downloads/", tag_out, "_checknb.dafs", sep=""), quote=F, 
#   row.names=F
# ); system(paste("open -a TextWrangler ~/Downloads/", tag_out, "_checknb.dafs", sep=""))
# }

if(tag_treat=="2dot2_"){
  mergechroms <- as.list(
    paste0(
      "NC64A_scaffold_",
      c(1, 1, 1, 1, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 11, 11, 12, 14, 16, 17, 18, 19, 19, 22, 23, 24 ,26 , 27, 32, 32), "_nuc"
    )
  )
  mergeposits <- list(
    c(240463, 240472, 240716), c(1808755, 1808756), 
    c(2479027, 2479028, 2479033, 2479034), c(3071221, 3071236), 
    c(255097, 255098, 255102), c(1010340, 1010346), 
    c(2195059, 2195068, 2195074, 2195083, 2195086, 2195089, 2195090, 2195091), 
    c(2195203, 2195209), c(1064328, 1064333, 1064334), 
    c(1676873, 1676879), c(1954491, 1954494), 
    c(247548, 247554), c(517236, 517239), 
    c(1414824, 1414827), c(387209, 387221, 387227, 387228, 387229, 387290), 
    c(1496798, 1496799), c(1036240, 1036243), 
    c(1140105, 1140115, 1140120), c(344562, 344564), 
    c(517996, 518343, 518346), c(820342, 820348, 820372, 820373), 
    c(55742, 55746), c(683146, 683149), 
    c(214708, 214711), c(212877, 212879, 212898, 212903), 
    c(218396, 218399, 218400), c(477054, 477060, 477075, 477078), 
    c(312336, 312338), c(262852, 262853, 262871), 
    c(280834, 280840, 280873, 280876, 280878)
  )
}else if(tag_treat=="3dot2_"){
  mergechroms <- as.list(
    paste0(
      "NC64A_scaffold_",
      c(1, 1, 1, 1, 1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 9, 9, 10, 10, 10, 10, 11, 12, 13, 14, 14, 16, 17, 17, 19, 19, 22, 23, 23, 23, 25, 26, 26, 31, 32, 32, 32, 32, 37, 38, 140, 440), "_nuc"
    )
  )
  mergeposits <- list(
    c(240463, 240472, 240716), c(452706, 452710, 452713, 452720), 
    c(612205, 612206), c(1808755, 1808756), 
    c(2479027, 2479033, 2479034), c(3071221, 3071236), 
    c(882430, 882439), c(255097, 255098, 255102), 
    c(616592, 616598), c(1010340, 1010346), 
    c(1029832, 1029836, 1029840), 
    c(2195059, 2195068, 2195074, 2195086, 2195089, 2195090, 2195091, 2195203, 2195209), 
    c(658904, 659082, 659091), c(1676873, 1676879), 
    c(247548, 247554), c(517230, 517239), 
    c(1414824, 1414827), c(1165229, 1165255, 1165270), 
    c(1371433, 1371454), c(1686946, 1686948), 
    c(349328, 349338), c(1089892, 1090286), 
    c(2837, 2858), c(3701, 3710, 3713), 
    c(4355, 4358), c(63934, 63936), 
    c(1496798, 1496799), c(1036240, 1036243), 
    c(1255208, 1255221), c(452988, 453000, 453003, 453024), 
    c(1140105, 1140115, 1140120), c(344562, 344564), 
    c(518343, 518346), c(1053226, 1053232), 
    c(55742, 55746), c(680125, 680128, 680131), 
    c(214708, 214711), c(212879, 212898, 212903), 
    c(266797, 266801, 266803, 266812, 266836), c(789048, 789049, 789052), 
    c(142702, 142717, 142718, 142725, 143095, 143098, 143107), 
    c(477054, 477060, 477075, 477078), c(477912, 477913), 
    c(298796, 298799), c(262852, 262853, 262871), 
    c(280834, 280840, 280873, 280876, 280878), c(281544, 281545), 
    c(359281, 359284, 359290, 359293), c(120794, 120796, 120799), 
    c(156289, 156292), c(1285, 1289, 1290, 1293, 1294, 1297, 1298, 1308), 
    c(961, 963, 969, 971)
    
    
  )
}else if(tag_treat=="4dot2_"){
  mergechroms <- as.list(
    paste0(
      "NC64A_scaffold_",
      c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 12, 12, 13, 14, 15, 15, 17, 18, 18, 19, 19, 20, 22, 23, 23, 23, 24, 25, 26, 31, 31, 32, 32, 32, 37, 38, 140, 746), "_nuc"
    )
  )
  mergeposits <- list(
    c(240463, 240472, 240716), c(277150, 277154), 
    c(452706, 452710, 452713, 452720), c(612205, 612206), 
    c(1808755, 1808756), c(2156841, 2156856), 
    c(2159062, 2159071), c(2159128, 2159152, 2159155), 
    c(2479027, 2479033, 2479034), c(3071221, 3071236), 
    c(882157, 882400, 882409, 882430, 882439), c(255097, 255098, 255102), 
    c(616592, 616598), c(2267240, 2267243, 2267252), 
    c(291371, 291732), 
    c(699012, 699015, 699274, 699278, 699282, 699286, 699949, 699957, 699964, 699970, 699979, 699985), 
    c(1010340, 1010346), c(1029836, 1029840), 
    c(1151004, 1151005), c(2195059, 2195068, 2195074, 2195083, 2195086, 2195089, 2195090, 2195091), 
    c(2195203, 2195209), c(357209, 357215, 357224, 357227), 
    c(533183, 533185), c(658904, 659082, 659091, 659118, 659536, 659545), 
    c(1676873, 1676879), c(1987277, 1987280), 
    c(300065, 300066), c(1414824, 1414827), 
    c(863444, 863448, 863454, 863456), c(869541, 869549, 869557), 
    c(1165255, 1165270), c(1371433, 1371454), 
    c(149652, 149664), c(1390980, 1390981, 1390986), 
    c(1686946, 1686948), c(349328, 349338), 
    c(1090286, 1090301), c(1613948, 1613950), 
    c(2837, 2858), c(3701, 3710, 3713), 
    c(4355, 4358), c(63934, 63936), 
    c(404682, 404691, 404698, 404699, 404701), 
    c(651541, 651543, 651547, 651562, 651574, 651583, 651674, 651676), 
    c(1496798, 1496799), c(216832, 216835, 216839), 
    c(1036240, 1036243), c(1254845, 1254848, 1255208, 1255221), 
    c(1140105, 1140115, 1140120), c(94241, 94242), 
    c(592888, 592909), c(518343, 518346), 
    c(245134, 245185, 245186), c(933138, 933142), 
    c(55742, 55746), c(680125, 680128), 
    c(225616, 225617), c(214708, 214711), 
    c(212877, 212879, 212898, 212903), c(272896, 272898), 
    c(789351, 789352, 789355), c(212064, 212070), 
    c(142702, 142717, 142718, 142725, 143083, 143095, 143098, 143107), 
    c(477054, 477060, 477075, 477078), c(343824, 344050), 
    c(418398, 418403), c(187628, 187629, 187652, 187654), 
    c(280834, 280840, 280873, 280876, 280878, 281544, 281545, 281547), 
    c(359281, 359284, 359290, 359293), c(120794, 120796, 120799), 
    c(156289, 156292), c(1285, 1287, 1289, 1290, 1291, 1293, 1294, 1297, 1298, 1308), 
    c(574, 581)
  )
}

## write these to file to check reproducibility between replicates; 
## these mutations are only considered reproducible if they're completely
## identical SNPs.
# if(file.exists(paste0(
#   exp2017_dir, "genom/snpol_dir/", tag_out, "_mergeSNPs.pos"
# ))){
#   system(paste0("rm ", exp2017_dir, "genom/snpol_dir/", tag_out, "_mergeSNPs.pos"))
# }
# ## to make *_multipos.R run, the object below has to have the same 
# ## number of columns: 
# lpos <- 14
# # lpos <- sapply(mergeposits, length) %>% max
# 
# for(i in 1:length(mergechroms)){
#   write.table(
#     cbind(
#       mergechroms[[i]],
#       t(c(mergeposits[[i]], rep(0, lpos - length(mergeposits[[i]]))))
#     ), append=TRUE, quote=FALSE,
#     row.names=FALSE, col.names=FALSE, sep="\t",
#     file=paste0(exp2017_dir, "genom/snpol_dir/", tag_out, "_mergeSNPs.pos")
#   )
# }
#system(paste0("open ~/Downloads/", tag_out, "_mergeSNPs.pos"))


## merge positions
dafs_var25_3nas_mr <- mergePositionsDaf(
  input=cbind(sync_var25_3nas[, 1:2],
              dafs_var25_3nas),
  chroms=mergechroms, posits=mergeposits
)
## find corresponding .sync rows
sync_var25_3nas_mr <- data.frame()
for(i in 1:nrow(dafs_var25_3nas_mr)){
  sync_var25_3nas_mr <- rbind(
    sync_var25_3nas_mr, subset(
      sync_var25_3nas, chrom==dafs_var25_3nas_mr[i, 1] &
        pos==dafs_var25_3nas_mr[i, 2]
    )
  )
}

## remove chrom and pos from dafs_*
dafs_var25_3nas_mr <- as.matrix(dafs_var25_3nas_mr[, -(1:2)])
dim(dafs_var25_3nas); dim(dafs_var25_3nas_mr)


#############################
## plot with ggplot2, presentation quality
colyr <- colorRampPalette(colors=c("yellow", "dark red"), space="rgb")
set.seed(2303)
colbg <- colorRampPalette(colors=brewer.pal(n=9, "Greens"), space="rgb")(nrow(dafs_var25_3nas))[sample(1:nrow(dafs_var25_3nas))]
# plot(1:102, rep(0, 102), pch=16, cex=1, col=colbg)

dafs_plotdf <- data.frame(
  tp=rep(tps_hicov, nrow(dafs_var25_3nas)),
  daf=as.numeric(t(dafs_var25_3nas)),
  allele=factor(rep(
    1:nrow(dafs_var25_3nas), each=ncol(dafs_var25_3nas)
  ))
)

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
    axis.ticks.x=element_line(size=1),#, colour="#E69F00"),
  ) +
  # geom_vline(
  #   data=data.frame(tp=tps_hicov), aes(xintercept=tp), colour="red",
  #   alpha=.6, linetype="dashed", size=.4
  # ) +
  labs(x="", y="")
  # labs(x="Time (days)", y="10log( Population size )")

gplot_dafs <- ggplot() + geom_line(
  data=dafs_plotdf[!is.na(dafs_plotdf$daf), ],
  aes(x=tp, y=daf, colour=allele), size=.5, alpha=.95
) + scale_colour_manual(values=colbg) +
  geom_point(
    data=dafs_plotdf, aes(x=tp, y=daf, colour=allele), size=3
  ) +
  scale_x_continuous(breaks=10*(0:12)) +
  scale_y_continuous(breaks=.2*(0:5)) +
  coord_cartesian(x=c(0, 100), y=c(0, 1)) +
  guides(colour=FALSE) + theme_bw() + theme(
    panel.grid=element_blank(),
    # panel.border=element_rect(colour="#E69F00"),
    # legend.title=element_text(colour="#E69F00"),
    axis.text.y=element_text(size=20),#, colour="#E69F00"),
    axis.ticks.y=element_line(size=1),#, colour="#E69F00"),
    axis.text.x=element_text(size=20),#, colour="#E69F00"),
    axis.ticks.x=element_line(size=1)#, colour="#E69F00"),
  ) +
  geom_vline(
    data=data.frame(tp=tps_hicov), aes(xintercept=tp), colour="red",
    alpha=.6, linetype="dashed", size=.4
  ) +
  labs(x="Time (days)", y="Derived allele frequency")
  labs(x="", y="")

pdf(paste(exp2017_dir, "genom/figs/fig_", tag_out, "_var25.pdf", sep=""),
    width=9, height=6)
# pdf(paste0("~/Documents/Study/AVW9/fig_", tag_out, "_var25.pdf"), 
#     width=9, height=6)
gg_multiplot(
  gplot_ps + theme(plot.margin=unit(c(.2, .2, .2, .8), "cm")),
  gplot_dafs + theme(plot.margin=unit(c(.2, 2.7, .2, .2), "cm"))
)
dev.off(); system(paste("open ", exp2017_dir, "genom/figs/fig_", tag_out, "_var25.pdf", sep=""))

}


#############################
## write to file; needed for the next step, to plot d(af) versus d(popsize)
## add a column ac.12
write.table(
  sync_var25_3nas, quote=FALSE, row.names=FALSE,
  col.names=FALSE, sep="\t",
  file=paste(exp2017_dir, "genom/ppl_dir/", tag_out, "_var25_filtered.sync", sep="")
)
write.table(
  cbind(sync_var25_3nas[, 1:3], dafs_var25_3nas),
  quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t",
  file=paste(exp2017_dir, "genom/daf_dir/", tag_out, "_var25_filtered.daf", sep="")
)
write.table(
  cbind(sync_var25_3nas_mr[, 1:3], dafs_var25_3nas_mr),
  quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t",
  file=paste0(exp2017_dir, "genom/daf_dir/", tag_out, "_var25_filtered_mr.daf")
)

}
