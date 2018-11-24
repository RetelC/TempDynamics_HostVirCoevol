#######################################
## Cas Retel
## 2018.04.13
## cas.retel@eawag.ch
#######################################
## After filtering for artefacts, read in derived allele frequencies
## of all three replicates, check if there's overlap between replicates. 
## Generate an 'allpos'-object containing a row for every position 
## filtered in at least one replicate, with information on in which
## replicate it is found, if it is single-position, it's annotation etc. 

## All annotation information comes from 
## http://genome.jgi.doe.gov/ChlNC64A_1/ChlNC64A_1.home.html

## Comes after nc64a_aaffiltering.R
rm(list=ls())
pkgs <- list("ggplot2", "magrittr", "RColorBrewer")
sapply(pkgs, require, character.only=T)
## files below can all be found at github.com/RetelC/TempDynamics_HostVirCoevol
source('~/Documents/Functions/gg_multiplot.R')
source('~/Documents/Functions/round_10e3.R')
source('~/Documents/HVInt/scripts/fs_syncCalculations.R')
source('~/Documents/HVInt/scripts/fs_syncFiltering.R')
source('~/Documents/HVInt/scripts/plotAfs.R')
source('~/Documents/HVInt/scripts/writeReadSync.R')

exp2017_dir <- "/Users/reteladmin/Documents/HVInt/Exp2017/Exp2017_30/"

## extract virus population sizes (observed and smoothed)
virhos_psdf <- read.csv(
  paste0(exp2017_dir, "demog/2dot2_hosvir0_logcounts_spline.csv")
)

## set timepoints with high enough coverage
tps_22 <- c(0, 7, 12, 14, 27, 51, 64, 69, 83, 99)
tps_32 <- c(0, 8, 12, 27, 29, 51, 64, 69, 83, 99)
tps_42 <- c(0, 8, 12, 14, 27, 29, 51, 83, 99)

## download positions that were merged into one polymorphism; 
## see 234dot2_D00-99_nc64a_aaffiltering.R
mergepos_22 <- read.table(paste0(
  exp2017_dir, "genom/snpol_dir/2dot2_D00-99_nc64a_mergeSNPs.pos"
), sep="\t", stringsAsFactors=FALSE)
mergepos_32 <- read.table(paste0(
  exp2017_dir, "genom/snpol_dir/3dot2_D00-99_nc64a_mergeSNPs.pos"
), sep="\t", stringsAsFactors=FALSE)
mergepos_42 <- read.table(paste0(
  exp2017_dir, "genom/snpol_dir/4dot2_D00-99_nc64a_mergeSNPs.pos"
), sep="\t", stringsAsFactors=FALSE)

mergepos_22$chrompos <- with(mergepos_22, paste(V1, V2, sep="::"))
mergepos_32$chrompos <- with(mergepos_32, paste(V1, V2, sep="::"))
mergepos_42$chrompos <- with(mergepos_42, paste(V1, V2, sep="::"))

## read in derived allele frequencies; not going via .sync because 
## multpos-dafs are averaged
dafs_22 <- readSync(
  paste0(exp2017_dir, "genom/daf_dir/2dot2_D00-99_nc64a_var25_filtered_mr.daf"), 
  tps=tps_22, header_v=TRUE
)
dafs_22 <- cbind(chrompos=with(dafs_22, paste(chrom, pos, sep="::")), dafs_22)
dafs_32 <- readSync(
  paste0(exp2017_dir, "genom/daf_dir/3dot2_D00-99_nc64a_var25_filtered_mr.daf"), 
  tps=tps_32, header_v=TRUE
)
dafs_32 <- cbind(chrompos=with(dafs_32, paste(chrom, pos, sep="::")), dafs_32)
dafs_42 <- readSync(
  paste0(exp2017_dir, "genom/daf_dir/4dot2_D00-99_nc64a_var25_filtered_mr.daf"), 
  tps=tps_42, header_v=TRUE
)
dafs_42 <- cbind(chrompos=with(dafs_42, paste(chrom, pos, sep="::")), dafs_42)

## create a data frame with one (unique) row per position found
## create variables chromnum for sorting and chrompos for matching
allpos <- rbind.data.frame(dafs_22[, 1:3], dafs_32[, 1:3], dafs_42[, 1:3])
allpos$pos <- as.numeric(allpos$pos)
allpos$chromnum <- gsub(
  "(NC64A_scaffold_)(\\d+)(_[npm][uli][cat])", replacement="\\2", 
  allpos$chrom, perl=TRUE
) %>% as.numeric
allpos <- allpos[!duplicated(allpos$chrompos), ]
allpos <- allpos[order(allpos$chromnum, allpos$pos), ]

## create vectors called found.[234]2 tracking in which replicate
## an allele is found
allpos$found.42 <- allpos$found.32 <- allpos$found.22 <- rep(NA, nrow(allpos))
for(i in 1:nrow(allpos)){
  allpos$found.22[i] <- c(FALSE, TRUE)[
    1+(any((allpos$chrom[i] == dafs_22$chrom) & (allpos$pos[i] == dafs_22$pos)))
    ]
  allpos$found.32[i] <- c(FALSE, TRUE)[
    1+(any((allpos$chrom[i] == dafs_32$chrom) & (allpos$pos[i] == dafs_32$pos)))
    ]
  allpos$found.42[i] <- c(FALSE, TRUE)[
    1+(any((allpos$chrom[i] == dafs_42$chrom) & (allpos$pos[i] == dafs_42$pos)))
    ]
}

## create vectors called multpos.[234]2 that become TRUE when the 
## position reflects a multi-position polymorphism. This information 
## comes from mergepos
allpos$multpos.22 <- (allpos$chrompos %in% mergepos_22$chrompos)
allpos$multpos.32 <- (allpos$chrompos %in% mergepos_32$chrompos)
allpos$multpos.42 <- (allpos$chrompos %in% mergepos_42$chrompos)

## Now let's see if multiple mutations are identical. 
identmp_2242 <- 
  matrix(rep(F, nrow(mergepos_22)* nrow(mergepos_42)), nrow=nrow(mergepos_42))
identmp_2232 <- 
  matrix(rep(F, nrow(mergepos_22)* nrow(mergepos_32)), nrow=nrow(mergepos_32))
identmp_3242 <- 
  matrix(rep(F, nrow(mergepos_32)* nrow(mergepos_42)), nrow=nrow(mergepos_42))

for(i in 1:nrow(mergepos_22)){
  for(j in 1:nrow(mergepos_32)){
    identmp_2232[j, i] <- all(mergepos_22[i, ] == mergepos_32[j, ])
  }
  for(j in 1:nrow(mergepos_42)){
    identmp_2242[j, i] <- all(mergepos_22[i, ] == mergepos_42[j, ])
  }
}
for(i in 1:nrow(mergepos_32)){
  for(j in 1:nrow(mergepos_42)){
    identmp_3242[j, i] <- all(mergepos_32[i, ] == mergepos_42[j, ])
  }
}

## create three vectors ol (for overlapping) that reflect whether
## a position is shared between replicates; for single-pos, this is 
## identical to found.[234]2; for multi-pos, this becomes TRUE if all
## merged positions exactly overlap. 
allpos$ol.34 <- allpos$ol.24 <- allpos$ol.23 <- F
allpos$ol.23[with(
  allpos, (found.22 & !multpos.22) & (found.32 & !multpos.32)
)] <- TRUE
allpos$ol.23[with(
  allpos, which(multpos.22)[apply(identmp_2232, 2, any)]
)] <- TRUE
allpos$ol.24[with(
  allpos, (found.22 & !multpos.22) & (found.42 & !multpos.42)
)] <- TRUE
allpos$ol.24[with(
  allpos, which(multpos.22)[apply(identmp_2242, 2, any)]
)] <- TRUE
allpos$ol.34[with(
  allpos, (found.32 & !multpos.32) & (found.42 & !multpos.42)
)] <- TRUE
allpos$ol.34[with(
  allpos, which(multpos.32)[apply(identmp_3242, 2, any)]
)] <- TRUE


## number of multi-position polymorphisms
sapply(list(mergepos_22, mergepos_32, mergepos_42), nrow)
with(allpos, c(table(multpos.22), table(multpos.32), table(multpos.42)))
## number of overlapping positions; see also 
## snpol_dir/234dot2_overlappingpositions.key
sapply(list(dafs_22, dafs_32, dafs_42), nrow)
with(allpos, c(sum(found.22), sum(found.32), sum(found.42)))
with(allpos, table(found.22, found.42, found.32))
with(allpos, table(found.22, found.32))
with(allpos, table(found.22, found.42))
with(allpos, table(found.32, found.42))

## create a vector with the number of replicates in which a locus is variable
allpos$reps.found <- with(allpos, found.22 + found.32 + found.42)

## select dafs found in multiple replicates
dafs_22_1rep <- dafs_22[
  dafs_22$chrompos %in% with(allpos, chrompos[found.22 & !found.32 & !found.42]), 
  ]
dafs_22_2rep <- dafs_22[
  dafs_22$chrompos %in% with(allpos, chrompos[found.22 & xor(found.32, found.42)]), 
  ]
dafs_22_2or3rep <- dafs_22[
  dafs_22$chrompos %in% with(allpos, chrompos[found.22 & (found.32 | found.42)]), 
  ]
dafs_22_3rep <- dafs_22[
  dafs_22$chrompos %in% with(allpos, chrompos[found.22 & found.32 & found.42]), 
  ]
## next replicate
dafs_32_1rep <- dafs_32[
  dafs_32$chrompos %in% with(allpos, chrompos[found.32 & !found.22 & !found.42]), 
  ]
dafs_32_2rep <- dafs_32[
  dafs_32$chrompos %in% with(allpos, chrompos[found.32 & xor(found.22, found.42)]), 
  ]
dafs_32_2or3rep <- dafs_32[
  dafs_32$chrompos %in% with(allpos, chrompos[found.32 & (found.22 | found.42)]), 
  ]
dafs_32_3rep <- dafs_32[
  dafs_32$chrompos %in% with(allpos, chrompos[found.22 & found.32 & found.42]), 
  ]
## next replicate
dafs_42_1rep <- dafs_42[
  dafs_42$chrompos %in% with(allpos, chrompos[found.42 & !found.22 & !found.32]), 
  ]
dafs_42_2rep <- dafs_42[
  dafs_42$chrompos %in% with(allpos, chrompos[found.42 & xor(found.22, found.32)]), 
  ]
dafs_42_2or3rep <- dafs_42[
  dafs_42$chrompos %in% with(allpos, chrompos[found.42 & (found.22 | found.32)]), 
  ]
dafs_42_3rep <- dafs_42[
  dafs_42$chrompos %in% with(allpos, chrompos[found.22 & found.32 & found.42]), 
  ]

#######################################
## Add annotations to allpos-object: predicted phenotypic effect was 
## assessed using snpeff, which outputs vcf files read in below: 
vcf_22 <- read.table(
  paste0(exp2017_dir, "genom/snpeff_dir/2dot2_D00-99_nc64a_var25_se_mp.vcf"), 
  header=T, sep="\t", stringsAsFactors=FALSE
)
vcf_32 <- read.table(
  paste0(exp2017_dir, "genom/snpeff_dir/3dot2_D00-99_nc64a_var25_se_mp.vcf"), 
  header=T, sep="\t", stringsAsFactors=FALSE
)
vcf_42 <- read.table(
  paste0(exp2017_dir, "genom/snpeff_dir/4dot2_D00-99_nc64a_var25_se_mp.vcf"), 
  header=T, sep="\t", stringsAsFactors=FALSE
)
vcf_22 <- cbind(
  chrompos=with(vcf_22, paste(CHROM, POS, sep="::")), vcf_22
)
vcf_32 <- cbind(
  chrompos=with(vcf_32, paste(CHROM, POS, sep="::")), vcf_32
)
vcf_42 <- cbind(
  chrompos=with(vcf_42, paste(CHROM, POS, sep="::")), vcf_42
)

## match vcf annotations with allpos-rows
head(allpos)
allpos$muts.42 <- allpos$muts.32 <- allpos$muts.22 <- 
  allpos$mutc.42 <- allpos$mutc.32 <- allpos$mutc.22 <- rep(NA, nrow(allpos))
allpos$mutc.22[allpos$chrompos %in% vcf_22$chrompos] <- vcf_22$MUT.CLASS
allpos$mutc.32[allpos$chrompos %in% vcf_32$chrompos] <- vcf_32$MUT.CLASS
allpos$mutc.42[allpos$chrompos %in% vcf_42$chrompos] <- vcf_42$MUT.CLASS

allpos$muts.22[allpos$chrompos %in% vcf_22$chrompos] <- vcf_22$MUT.SEVERITY
allpos$muts.32[allpos$chrompos %in% vcf_32$chrompos] <- vcf_32$MUT.SEVERITY
allpos$muts.42[allpos$chrompos %in% vcf_42$chrompos] <- vcf_42$MUT.SEVERITY

#############################
## add protein ID 
head(allpos)
allpos$proteinId <- allpos$gff.entry <- rep(NA, nrow(allpos))

## gff file downloaded from http://genome.jgi.doe.gov/ChlNC64A_1/ChlNC64A_1.home.html
nc64a_gff <- read.table(
  "~/Documents/HVInt/Chlorella/ChlNC64A_1/genes_JF2016_GENESONLY.gff", 
  sep="\t", header=F, stringsAsFactors=F
)

## for every variable position, check if it's inside a gene. 
## If yes, store the whole ninth column of the .gff, and extract
## protein ID from this. 
## ! note that there are two mutations that match (the same) 
## two genes; gene ID's 50233 and 6284 both have METABOLISM KOG-function
for(i in 1:nrow(allpos)){
  idx_insidegene <- which( (allpos[i, "chrom"] == nc64a_gff$V1) & 
                             (allpos[i, "pos"] >= nc64a_gff$V4) & 
                             (allpos[i, "pos"] <= nc64a_gff$V5) )
  
  if(length(idx_insidegene) > 0){
    if(length(idx_insidegene) > 1){
      print(paste0("WARNING: position number ", i, " matches multiple genes"))
      allpos$gff.entry[i] <- nc64a_gff[idx_insidegene[1], 9]
    }else{
      allpos$gff.entry[i] <- nc64a_gff[idx_insidegene, 9]
    }
    
    ## extract protein ID number from this entry
    allpos$proteinId[i] <- 
      strsplit(allpos$gff.entry[i], split=";")[[1]] %>%
      (function(x1) grep("ID=", x=x1, value=TRUE)) %>% 
      (function(x2) strsplit(x2, split="_")[[1]]) %>% tail(n=1) %>% as.numeric
  }
}


allpos$proteinId
# writeLines(
#   as.character(allpos$proteinId), 
#   paste0(exp2017_dir, "genom/annot_dir/234dot2_D00-99_nc64a_protIDs.txt")
# )



#######################################
## PLOT

par(mfrow=c(3, 1))

# all(as.character(dafs_22_3rep$chrompos) == as.character(dafs_42_3rep$chrompos))
set.seed(2303)
col_3rep <- rainbow(nrow(dafs_22_3rep))[sample(nrow(dafs_22_3rep))]
col_2or3rep <- rainbow(sum(allpos$reps.found >= 2))[
  sample(sum(allpos$reps.found >= 2))
  ]
names(col_2or3rep) <- with(allpos, chrompos[allpos$reps.found >= 2])
col_22_1rep <- colorRampPalette(colors=brewer.pal(n=9, "Greens"), space="rgb")(
  nrow(dafs_22_1rep))[sample(1:nrow(dafs_22_1rep))]
col_32_1rep <- colorRampPalette(colors=brewer.pal(n=9, "Greens"), space="rgb")(
  nrow(dafs_32_1rep))[sample(1:nrow(dafs_32_1rep))]
col_42_1rep <- colorRampPalette(colors=brewer.pal(n=9, "Greens"), space="rgb")(
  nrow(dafs_42_1rep))[sample(1:nrow(dafs_42_1rep))]



pdf(paste0(exp2017_dir, "genom/figs/fig_234dot2_nc64a_D00-99_overlapping.pdf"), 
    width=5, height=12)
par(mfrow=c(3, 1))
plotAfs(
  as.matrix(dafs_22[, -(1:4)]), tps=tps_22, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="II.2: Variants found in two out of three replicates"
)
plotAfs(
  as.matrix(dafs_22_2or3rep[, -(1:4)]), tps=tps_22, alpha_v=.8, 
  col_v=col_2or3rep[as.character(dafs_22_2or3rep$chrompos)], 
  lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_32[, -(1:4)]), tps=tps_32, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="III.2: Variants found in two out of three replicates"
)
plotAfs(
  as.matrix(dafs_32_2or3rep[, -(1:4)]), tps=tps_32, alpha_v=.8, 
  col_v=col_2or3rep[as.character(dafs_32_2or3rep$chrompos)], 
  lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_42[, -(1:4)]), tps=tps_42, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="IV.2: Variants found in two out of three replicates"
)
plotAfs(
  as.matrix(dafs_42_2or3rep[, -(1:4)]), tps=tps_42, alpha_v=.8, 
  col_v=col_2or3rep[as.character(dafs_42_2or3rep$chrompos)], 
  lines_lwd=2, dots_cex=1, add=TRUE
)
dev.off(); system(paste0("open ", exp2017_dir, "genom/figs/fig_234dot2_nc64a_D00-99_overlapping.pdf"))

pdf(paste0(exp2017_dir, "genom/figs/fig_234dot2_nc64a_D00-99_threereps.pdf"), 
    width=5, height=12)
par(mfrow=c(3, 1))
plotAfs(
  as.matrix(dafs_22[, -(1:4)]), tps=tps_22, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="II.2: Variants found in all three replicates"
)
plotAfs(
  as.matrix(dafs_22_3rep[, -(1:4)]), tps=tps_22, alpha_v=.8, 
  col_v=col_3rep, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_32[, -(1:4)]), tps=tps_32, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="III.2: Variants found in all three replicates"
)
plotAfs(
  as.matrix(dafs_32_3rep[, -(1:4)]), tps=tps_32, alpha_v=.8, 
  col_v=col_3rep, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_42[, -(1:4)]), tps=tps_42, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="IV.2: Variants found in all three replicates"
)
plotAfs(
  as.matrix(dafs_42_3rep[, -(1:4)]), tps=tps_42, alpha_v=.8, 
  col_v=col_3rep, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
dev.off(); system(paste0("open ", exp2017_dir, "genom/figs/fig_234dot2_nc64a_D00-99_threereps.pdf"))

## colour only the unique variants
pdf(paste0(exp2017_dir, "genom/figs/fig_234dot2_nc64a_D00-99_unique.pdf"), 
    width=5, height=12)
par(mfrow=c(3, 1))
plotAfs(
  as.matrix(dafs_22[, -(1:4)]), tps=tps_22, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="II.2: Variants found in only this replicate"
)
plotAfs(
  as.matrix(dafs_22_1rep[, -(1:4)]), tps=tps_22, alpha_v=.8, 
  col_v=col_22_1rep, lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_32[, -(1:4)]), tps=tps_32, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="III.2: Variants found in only this replicate"
)
plotAfs(
  as.matrix(dafs_32_1rep[, -(1:4)]), tps=tps_32, alpha_v=.8, 
  col_v=col_32_1rep, lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_42[, -(1:4)]), tps=tps_42, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="IV.2: Variants found in only this replicate"
)
plotAfs(
  as.matrix(dafs_42_1rep[, -(1:4)]), tps=tps_42, alpha_v=.8, 
  col_v=col_42_1rep, lines_lwd=2, dots_cex=1, add=TRUE
)
dev.off(); system(paste0("open ", exp2017_dir, "genom/figs/fig_234dot2_nc64a_D00-99_unique.pdf"))

## write dafs to file; 
writeSync(dafs_22_2or3rep[, -1], fname=paste0(exp2017_dir, "genom/daf_dir/2dot2_D00-99_nc64a_overlapping.daf"), colnames_v=TRUE)
writeSync(dafs_32_2or3rep[, -1], fname=paste0(exp2017_dir, "genom/daf_dir/3dot2_D00-99_nc64a_overlapping.daf"), colnames_v=TRUE)
writeSync(dafs_42_2or3rep[, -1], fname=paste0(exp2017_dir, "genom/daf_dir/4dot2_D00-99_nc64a_overlapping.daf"), colnames_v=TRUE)

writeSync(dafs_22_3rep[, -1], fname=paste0(exp2017_dir, "genom/daf_dir/2dot2_D00-99_nc64a_threereps.daf"), colnames_v=TRUE)
writeSync(dafs_32_2or3rep[, -1], fname=paste0(exp2017_dir, "genom/daf_dir/3dot2_D00-99_nc64a_threereps.daf"), colnames_v=TRUE)
writeSync(dafs_42_2or3rep[, -1], fname=paste0(exp2017_dir, "genom/daf_dir/4dot2_D00-99_nc64a_threereps.daf"), colnames_v=TRUE)

## proteinId can be matched with KOG, KEGG and GO information. 
