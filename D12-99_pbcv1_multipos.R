#######################################
## Cas Retel
## 2018.07.31
## cas.retel@eawag.ch
#######################################
## After filtering for artefacts, read in derived allele frequencies
## of all three replicates, check if there's overlap between replicates. 
## Generate an allpos object containing a row for every position 
## filtered in at least one replicate, with information on in which
## replicate it is found, if it is single-position, it's annotation etc. 

rm(list=ls())
pkgs <- list("ggplot2", "magrittr", "RColorBrewer")
sapply(pkgs, require, character.only=T)
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
tps_v22 <- c(16, 21, 27, 29, 41, 51, 64, 69, 83, 99)
tps_v32 <- c(15, 21, 27, 29, 41, 51, 64, 69, 83, 99)
tps_v42 <- c(15, 21, 27, 29, 35, 51, 64, 70, 83, 99)


######################################
## VIRUS GENOMICS: read in .sync files of filtered SNP datasets
sync_v22 <- readSync(
  paste0(exp2017_dir, 'genom/ppl_dir/2dot2_D12-99_pbcv1_var5_filtered.sync'), tps_v22
)
sync_v32 <- readSync(
  paste0(exp2017_dir, 'genom/ppl_dir/3dot2_D12-99_pbcv1_var5_filtered.sync'), tps_v32
)
sync_v42 <- readSync(
  paste0(exp2017_dir, 'genom/ppl_dir/4dot2_D12-99_pbcv1_var5_filtered.sync'), tps_v42
)
## calculate derived allele frequencies
dafs_v22 <- cbind(sync_v22[, 1:3], recodeDafs(syncToDafs(sync_v22)))
dafs_v32 <- cbind(sync_v32[, 1:3], recodeDafs(syncToDafs(sync_v32)))
dafs_v42 <- cbind(sync_v42[, 1:3], recodeDafs(syncToDafs(sync_v42)))

## create chrompos column to match across replicates
options(stringsAsFactors = FALSE)
sync_v22 <- cbind(chrompos=with(sync_v22, paste(chrom, pos, sep="::")), sync_v22)
sync_v32 <- cbind(chrompos=with(sync_v32, paste(chrom, pos, sep="::")), sync_v32)
sync_v42 <- cbind(chrompos=with(sync_v42, paste(chrom, pos, sep="::")), sync_v42)
dafs_v22 <- cbind(chrompos=with(dafs_v22, paste(chrom, pos, sep="::")), dafs_v22)
dafs_v32 <- cbind(chrompos=with(dafs_v32, paste(chrom, pos, sep="::")), dafs_v32)
dafs_v42 <- cbind(chrompos=with(dafs_v42, paste(chrom, pos, sep="::")), dafs_v42)


## create allpos object with one row per derived allele frequency found
allpos <- rbind.data.frame(sync_v22[, 1:3], sync_v32[, 1:3], sync_v42[, 1:3])
allpos$pos <- as.numeric(allpos$pos)
## chromnum is obsolete since PBCV1 only has one scaffold. 
allpos$chromnum <- gsub(
  "(PBCV1_scaffold_)(\\d)", replacement="\\2", 
  allpos$chrom, perl=TRUE
) %>% as.numeric
## remove duplicate rows
allpos <- allpos[!duplicated(allpos$chrompos), ]
allpos <- allpos[order(allpos$chromnum, allpos$pos), ]
nrow(allpos) ## 21 SNPs total, checks out

## create vectors called found.[234]2 tracking in which replicate
## an allele is found
allpos$found.42 <- allpos$found.32 <- allpos$found.22 <- rep(NA, nrow(allpos))
for(i in 1:nrow(allpos)){
  allpos$found.22[i] <- c(FALSE, TRUE)[
    1+(any((allpos$chrom[i] == dafs_v22$chrom) & (allpos$pos[i] == dafs_v22$pos)))
    ]
  allpos$found.32[i] <- c(FALSE, TRUE)[
    1+(any((allpos$chrom[i] == dafs_v32$chrom) & (allpos$pos[i] == dafs_v32$pos)))
    ]
  allpos$found.42[i] <- c(FALSE, TRUE)[
    1+(any((allpos$chrom[i] == dafs_v42$chrom) & (allpos$pos[i] == dafs_v42$pos)))
    ]
}
head(allpos)

## create three vectors ol (for overlapping) that reflect whether
## a position is shared between replicates
allpos$ol.34 <- allpos$ol.24 <- allpos$ol.23 <- F
allpos$ol.23[with(
  allpos, (found.22 & found.32)
)] <- TRUE
allpos$ol.24[with(
  allpos, (found.22 & found.42)
)] <- TRUE
allpos$ol.34[with(
  allpos, (found.32 & found.42)
)] <- TRUE

## number of overlapping positions; see also 
## snpol_dir/234dot2_overlappingpositions.key
sapply(list(dafs_v22, dafs_v32, dafs_v42), nrow)
with(allpos, c(sum(found.22), sum(found.32), sum(found.42)))
with(allpos, table(found.22, found.42, found.32))
with(allpos, table(found.22, found.32))
with(allpos, table(found.22, found.42))
with(allpos, table(found.32, found.42))

## create a vector with the number of replicates in which a locus is variable
allpos$reps.found <- with(allpos, found.22 + found.32 + found.42)

## select dafs found in multiple replicates
dafs_v22_2reps <- dafs_v22[
  dafs_v22$chrompos %in% with(allpos, chrompos[reps.found >= 2]), 
]
dafs_v32_2reps <- dafs_v32[
  dafs_v32$chrompos %in% with(allpos, chrompos[reps.found >= 2]), 
]
dafs_v42_2reps <- dafs_v42[
  dafs_v42$chrompos %in% with(allpos, chrompos[reps.found >= 2]), 
]

dafs_v22_3reps <- dafs_v22[
  dafs_v22$chrompos %in% with(allpos, chrompos[reps.found >= 3]), 
  ]
dafs_v32_3reps <- dafs_v32[
  dafs_v32$chrompos %in% with(allpos, chrompos[reps.found >= 3]), 
  ]
dafs_v42_3reps <- dafs_v42[
  dafs_v42$chrompos %in% with(allpos, chrompos[reps.found >= 3]), 
  ]

## see 234dot2_D00-99_writeVCFs_plotAnns.R; 
## predicted phenotypic effects were calculated with snpeff, and are matched
## to allele frequencies below. 
vcf_v22 <- read.table(
  paste0(exp2017_dir, "genom/snpeff_dir/2dot2_D12-99_pbcv1_var5_se_mp.vcf"), 
  header=T, sep="\t", stringsAsFactors=FALSE
)
vcf_v32 <- read.table(
  paste0(exp2017_dir, "genom/snpeff_dir/3dot2_D12-99_pbcv1_var5_se_mp.vcf"), 
  header=T, sep="\t", stringsAsFactors=FALSE
)
vcf_v42 <- read.table(
  paste0(exp2017_dir, "genom/snpeff_dir/4dot2_D12-99_pbcv1_var5_se_mp.vcf"), 
  header=T, sep="\t", stringsAsFactors=FALSE
)
vcf_v22 <- cbind(
  chrompos=with(vcf_v22, paste(CHROM, POS, sep="::")), vcf_v22
)
vcf_v32 <- cbind(
  chrompos=with(vcf_v32, paste(CHROM, POS, sep="::")), vcf_v32
)
vcf_v42 <- cbind(
  chrompos=with(vcf_v42, paste(CHROM, POS, sep="::")), vcf_v42
)

## match vcf annotations with allpos-rows
head(allpos)
allpos$muts.42 <- allpos$muts.32 <- allpos$muts.22 <- 
  allpos$mutc.42 <- allpos$mutc.32 <- allpos$mutc.22 <- rep(NA, nrow(allpos))
allpos$mutc.22[allpos$chrompos %in% vcf_v22$chrompos] <- vcf_v22$MUT.CLASS
allpos$mutc.32[allpos$chrompos %in% vcf_v32$chrompos] <- vcf_v32$MUT.CLASS
allpos$mutc.42[allpos$chrompos %in% vcf_v42$chrompos] <- vcf_v42$MUT.CLASS

allpos$muts.22[allpos$chrompos %in% vcf_v22$chrompos] <- vcf_v22$MUT.SEVERITY
allpos$muts.32[allpos$chrompos %in% vcf_v32$chrompos] <- vcf_v32$MUT.SEVERITY
allpos$muts.42[allpos$chrompos %in% vcf_v42$chrompos] <- vcf_v42$MUT.SEVERITY

#############################
## add ID, dbxref geneID, gene name
head(allpos)
allpos$gene.name <- allpos$dbxref.geneID <- allpos$gene.ID <- 
  allpos$gff.entry <- rep(NA, nrow(allpos))

pbcv1_gff <- read.table(
  "~/Documents/HVInt/Chlorella/ChlNC64A_1/PBCV1_genes.gff", 
  sep="\t", header=F, stringsAsFactors=F
)
pbcv1_gff <- pbcv1_gff[pbcv1_gff$V3 == "gene", ]

## for every variable position, check if it's inside a gene. 
## If yes, store the whole ninth column of the .gff, and extract
## the other columns from it
## ! positions 18839 and 24315 match multiple gene entries..: 
# pbcv1_gff[which( (allpos[2, "pos"] >= pbcv1_gff$V4) & (allpos[2, "pos"] <= pbcv1_gff$V5) ), ]
# pbcv1_gff[which( (allpos[3, "pos"] >= pbcv1_gff$V4) & (allpos[3, "pos"] <= pbcv1_gff$V5) ), ]
## always a forward and reverse strand gene. 
## the second position actually is inside two overlapping reverse-strand genes. 

for(i in 1:nrow(allpos)){
  idx_insidegene <- which( (allpos[i, "pos"] >= pbcv1_gff$V4) & 
                             (allpos[i, "pos"] <= pbcv1_gff$V5) )
  
  if(length(idx_insidegene) > 0){
    if(length(idx_insidegene) > 1){
      print(paste0("WARNING: position number ", i, " matches multiple genes"))
      allpos$gff.entry[i] <- pbcv1_gff[idx_insidegene[1], 9]
    }else{
      allpos$gff.entry[i] <- pbcv1_gff[idx_insidegene, 9]
    }
    ## extract gene ID number from this entry
    allpos$gene.ID[i] <- strsplit(allpos$gff.entry[i], split=";")[[1]] %>%
      (function(x1) grep("ID=", x=x1, value=TRUE)) %>% 
      (function(x2) strsplit(x2, split="gene")[[1]]) %>% tail(n=1) %>% as.numeric
    ## extract gene name from this entry
    allpos$gene.name[i] <- strsplit(allpos$gff.entry[i], split=";")[[1]] %>% 
      (function(x1) grep("Name=", x=x1, value=TRUE)) %>%
      (function(x2) strsplit(x2, split="=")[[1]]) %>% tail(n=1)
    ## extract dbxref gene ID from this entry
    allpos$dbxref.geneID[i] <- strsplit(allpos$gff.entry[i], split=";")[[1]] %>%
      (function(x1) grep("Dbxref", x=x1, value=TRUE)) %>% 
      (function(x2) strsplit(x2, split="GeneID:")[[1]]) %>% tail(n=1) %>% as.numeric
  }
}

allpos$gene.ID
allpos$gene.name
table(allpos$gene.name)
## 9 out of 21 mutations are in the A540L gene. 
# writeSync(allpos, fname=paste0(
#   exp2017_dir, "genom/snpol_dir/234dot2_pbcv1_D00-99_allvariablepositions.txt"
# ), colnames_v=TRUE)


#######################################
## PLOT

set.seed(105); col_2reps <- sample(rainbow(sum(allpos$reps.found >= 2)))
col_2reps_22 <- rep(NA, nrow(dafs_v22_2reps))
col_2reps_32 <- rep(NA, nrow(dafs_v32_2reps))
col_2reps_42 <- rep(NA, nrow(dafs_v42_2reps))
for(i in 1:nrow(dafs_v22_2reps)){
  col_2reps_22[i] <- col_2reps[
    (dafs_v22_2reps$chrompos[i] == allpos$chrompos[allpos$reps.found >= 2])
  ]
}
for(i in 1:nrow(dafs_v32_2reps)){
  col_2reps_32[i] <- col_2reps[
    (dafs_v32_2reps$chrompos[i] == allpos$chrompos[allpos$reps.found >= 2])
    ]
}
for(i in 1:nrow(dafs_v42_2reps)){
  col_2reps_42[i] <- col_2reps[
    (dafs_v42_2reps$chrompos[i] == allpos$chrompos[allpos$reps.found >= 2])
    ]
}

pdf(paste0(exp2017_dir, "genom/figs/fig_234dot2_pbcv1_overlapping.pdf"), width=5, height=12)
par(mfrow=c(3, 1))
plotAfs(
  as.matrix(dafs_v22[, -(1:4)]), tps=tps_v22, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="II.2: SNPs found in at least two replicates"
)
plotAfs(
  as.matrix(dafs_v22_2reps[, -(1:4)]), tps=tps_v22, col_v=col_2reps_22, alpha_v=.8, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_v32[, -(1:4)]), tps=tps_v32, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="III.2: SNPs found in at least two replicates"
)
plotAfs(
  as.matrix(dafs_v32_2reps[, -(1:4)]), tps=tps_v32, col_v=col_2reps_32, alpha_v=.8, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_v42[, -(1:4)]), tps=tps_v42, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="IV.2: SNPs found in at least two replicates"
)
plotAfs(
  as.matrix(dafs_v42_2reps[, -(1:4)]), tps=tps_v42, col_v=col_2reps_42, alpha_v=.8, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
dev.off(); system(paste0("open ", exp2017_dir, "genom/figs/fig_234dot2_pbcv1_overlapping.pdf"))

pdf(paste0(exp2017_dir, "genom/figs/fig_234dot2_pbcv1_threereps.pdf"), width=5, height=12)
par(mfrow=c(3, 1))
plotAfs(
  as.matrix(dafs_v22[, -(1:4)]), tps=tps_v22, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="II.2: SNPs found in at least two replicates"
)
plotAfs(
  as.matrix(dafs_v22_3reps[, -(1:4)]), tps=tps_v22, col_v=2, alpha_v=.8, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_v32[, -(1:4)]), tps=tps_v32, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="III.2: SNPs found in at least two replicates"
)
plotAfs(
  as.matrix(dafs_v32_3reps[, -(1:4)]), tps=tps_v32, col_v=2, alpha_v=.8, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
plotAfs(
  as.matrix(dafs_v42[, -(1:4)]), tps=tps_v42, col_v="dark grey", alpha_v=.6, 
  lines_lwd=2/3, dots_cex=2/3, xlim_v=c(0, 100), ylim_v=c(0, 1), 
  title_v="IV.2: SNPs found in at least two replicates"
)
plotAfs(
  as.matrix(dafs_v42_3reps[, -(1:4)]), tps=tps_v42, col_v=2, alpha_v=.8, 
  lines_lwd=2, dots_cex=1, add=TRUE
)
dev.off(); system(paste0("open ", exp2017_dir, "genom/figs/fig_234dot2_pbcv1_threereps.pdf"))

