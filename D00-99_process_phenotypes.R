#!/bin/R
######################################
## Processing of host-virus phenotype data
## Coevolutionary treatments II.2, III.2 and IV.2
## name: Cas Retel
## date: 2018.04.03
## e-mail: cas.retel@eawag.ch
######################################


## read in packages
rm(list=ls())
pkgs <- list("ggplot2", "magrittr", "RColorBrewer")
sapply(pkgs, require, character.only=T)
source('~/Documents/Functions/gg_multiplot.R')
source('~/Documents/Functions/round_10e3.R')

col_h <- "#009E73"
col_h2 <- rgb(t(col2rgb("#009E73")), alpha=204, maxColorValue=255)
col_v <- "#E69F00"
col_v2 <- rgb(t(col2rgb("#E69F00")), alpha=204, maxColorValue=255)

exp2017_dir <- "/Users/reteladmin/Documents/HVInt/Exp2017/Exp2017_30/"
db_dir <- "/Users/reteladmin/Dropbox/Talks_Algae_Virus/data/"
# file_in <- paste0(
#   "all_time_shift_data_", 
#   c("II_2", "III_2", "IV_2")[tag_treat==c("2dot2_", "3dot2_", "4dot2_")], 
#   ".csv"
# )
# tag_out <- paste(tag_treat, "D00-99", sep="")

## read in data
df_raw_22 <- read.csv(
  "/Users/reteladmin/Dropbox/Talks_Algae_Virus/data/all_time_shift_data_II_2.csv", sep=";", stringsAsFactors=FALSE, header=TRUE
)
df_raw_32 <- read.csv(
  "/Users/reteladmin/Dropbox/Talks_Algae_Virus/data/all_time_shift_data_III_2.csv", sep=";", stringsAsFactors=FALSE, header=TRUE
)
df_raw_42 <- read.csv(
  "/Users/reteladmin/Dropbox/Talks_Algae_Virus/data/all_time_shift_data_IV_2.csv", sep=";", stringsAsFactors=FALSE, header=TRUE
)

#######################################
## INTERMEZZO: check a few things, if Lutz made any copying errors on 
## the labelling etc. 
# with(df_raw_42, table(host.day, clone))
# with(df_raw_42, table(virus.day, host.day))
# 
# df_raw_22[df_raw_22$virus.day=="20.4", ] ## empty cell; changed in files
# View(with(df_raw_32, df_raw_32[virus.day=="100" & host.day=="75", ])) 
#   ## row 4308 is labelled neg neg, but that should actually still be 100 and 75; 
#   ## changed in files
# View(with(df_raw_32, df_raw_32[virus.day=="76" & host.day=="29", ]))
# 
# View(with(df_raw_42, df_raw_42[clone=="7" & virus.day=="neg", ]))
## there are always only three negative growth rates instead of four
## INTERMEZZO: check a few things
#######################################

#######################################
## DATA PREPARATION
## neg has a space in the factor level
df_raw_22$virus.ID[df_raw_22$virus.ID=="neg "] <- "neg"
## virus day 76 should be matched to host day 75, and I'm recoding virus
## day 1 as day 12 to get the right contemporary comparison
df_raw_22$virus.day[df_raw_22$virus.day=="76"] <- "75"
df_raw_22$virus.day[df_raw_22$virus.day=="1"] <- "12"

## virus day should be a factor, with day 100 coming last
df_raw_22$virus.day.f <- with(
  df_raw_22, factor(
    virus.day, levels=c("neg", unique(virus.day)[unique(virus.day) != "neg"])
  )
)
df_raw_22$host.day.f <- as.factor(df_raw_22$host.day)
## add a column with host individual
df_raw_22$host.indiv <- with(df_raw_22, paste(host.day, clone, sep="."))

## obtain rows containing resistance information per clone
df_inf_22 <- subset(df_raw_22, mean.growth != "")

## neg has a space in the factor level
df_raw_32$virus.ID[df_raw_32$virus.ID=="neg "] <- "neg"
## virus day 76 should be matched to host day 75, and I'm recoding virus
## day 1 as day 12 to get the right contemporary comparison
df_raw_32$virus.day[df_raw_32$virus.day=="76"] <- "75"
df_raw_32$virus.day[df_raw_32$virus.day=="1"] <- "12"
## virus day should be a factor, with day 100 coming last
df_raw_32$virus.day.f <- with(
  df_raw_32, factor(
    virus.day, levels=c("neg", unique(virus.day)[unique(virus.day) != "neg"])
  )
)
df_raw_32$host.day.f <- as.factor(df_raw_32$host.day)
## add a column with host individual
df_raw_32$host.indiv <- with(df_raw_32, paste(host.day, clone, sep="."))

## obtain rows containing resistance information per clone
df_inf_32 <- subset(df_raw_32, mean.growth != "")

## neg has a space in the factor level
df_raw_42$virus.ID[df_raw_42$virus.ID=="neg "] <- "neg"
## virus day 76 should be matched to host day 75, and I'm recoding virus
## day 1 as day 12 to get the right contemporary comparison
df_raw_42$virus.day[df_raw_42$virus.day=="76"] <- "75"
df_raw_42$virus.day[df_raw_42$virus.day=="1"] <- "12"

## virus day should be a factor, with day 100 coming last
df_raw_42$virus.day.f <- with(
  df_raw_42, factor(
    virus.day, levels=c("neg", unique(virus.day)[unique(virus.day) != "neg"])
  )
)
df_raw_42$host.day.f <- as.factor(df_raw_42$host.day)
## add a column with host individual
df_raw_42$host.indiv <- with(df_raw_42, paste(host.day, clone, sep="."))

## obtain rows containing resistance information per clone
df_inf_42 <- subset(df_raw_42, mean.growth != "")
## DATA PREPARATION
#######################################





#######################################
## FITNESS
## per-clone fitness; growth rate in presence of contemporary virus
## (or in absence at day 00)
df_contemp_22 <- rbind(
  df_inf_22[with(df_inf_22, virus.ID=="neg" & host.day=="1"), ], 
  with(df_inf_22, df_inf_22[(host.day==virus.day), ])
)
## I'm a bit scared of R switching factor levels on me
df_meanfit_22 <- with(df_contemp_22, data.frame(
  day = as.numeric(levels(factor(host.day))), 
  mean.fitness = tapply(mean.growth, INDEX=host.day, FUN=mean), 
  sd.fitness = tapply(mean.growth, INDEX=host.day, FUN=sd)
))

df_contemp_32 <- rbind(
  df_inf_32[with(df_inf_32, virus.ID=="neg" & host.day=="1"), ], 
  with(df_inf_32, df_inf_32[(host.day==virus.day), ])
)
## I'm a bit scared of R switching factor levels on me
df_meanfit_32 <- with(df_contemp_32, data.frame(
  day = as.numeric(levels(factor(host.day))), 
  mean.fitness = tapply(mean.growth, INDEX=host.day, FUN=mean), 
  sd.fitness = tapply(mean.growth, INDEX=host.day, FUN=sd)
))

df_contemp_42 <- rbind(
  df_inf_42[with(df_inf_42, virus.ID=="neg" & host.day=="1"), ], 
  with(df_inf_42, df_inf_42[(host.day==virus.day), ])
)
## I'm a bit scared of R switching factor levels on me
df_meanfit_42 <- with(df_contemp_42, data.frame(
  day = as.numeric(levels(factor(host.day))), 
  mean.fitness = tapply(mean.growth, INDEX=host.day, FUN=mean), 
  sd.fitness = tapply(mean.growth, INDEX=host.day, FUN=sd)
))



# plot fitness, i.e. host growth in presence of contemporary virus
# pdf("~/Dropbox/Talks_Algae_Virus/First paper/figures/extfig_fitness.pdf", width=12, height=4)
par(mfrow=c(1, 3))
par(mai=c(0.42, 0.62, 0.12, 0.12))
plot(c(0, 100), c(0, .4), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=20*(0:5), cex.axis=1.66)# ; axis(side=1, at=20*(0:5), cex.axis=1.66) 
axis(side=2, at=.1*(0:4), las=1, cex.axis=1.66)
with(df_contemp_22, points(host.day, mean.growth, col=col_h, pch=16))
with(df_meanfit_22, lines(day, mean.fitness, col=col_h, lty=2, lwd=2))

plot(c(0, 100), c(0, .4), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=20*(0:5), cex.axis=1.66)# ; axis(side=1, at=20*(0:5), cex.axis=1.66) 
axis(side=2, at=.1*(0:4), las=1, cex.axis=1.66)
with(df_contemp_32, points(host.day, mean.growth, col=col_h, pch=16))
with(df_meanfit_32, lines(day, mean.fitness, col=col_h, lty=2, lwd=2))

plot(c(0, 100), c(0, .4), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=20*(0:5), cex.axis=1.66)# ; axis(side=1, at=20*(0:5), cex.axis=1.66) 
axis(side=2, at=.1*(0:4), las=1, cex.axis=1.66)
with(df_contemp_42, points(host.day, mean.growth, col=col_h, pch=16))
with(df_meanfit_42, lines(day, mean.fitness, col=col_h, lty=2, lwd=2))
# dev.off(); system(paste0("open ~/Dropbox/Talks_Algae_Virus/First\\ paper/figures/extfig_fitness.pdf"))

# pdf("~/Dropbox/Talks_Algae_Virus/First paper/figures/extfig_fitness_sd.pdf", width=12, height=4)
par(mfrow=c(1, 3))
par(mai=c(0.42, 0.62, 0.12, 0.12))
plot(c(0, 100), c(0, .12), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=20*(0:5), cex.axis=1.66)# ; axis(side=1, at=20*(0:5), cex.axis=1.66) 
axis(side=2, at=.02*(0:6), las=1, cex.axis=1.66)
with(df_meanfit_22, lines(day, sd.fitness, col=col_h, lty=6, lwd=2))

plot(c(0, 100), c(0, .12), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=20*(0:5), cex.axis=1.66)# ; axis(side=1, at=20*(0:5), cex.axis=1.66) 
axis(side=2, at=.02*(0:6), las=1, cex.axis=1.66)
with(df_meanfit_32, lines(day, sd.fitness, col=col_h, lty=6, lwd=2))

plot(c(0, 100), c(0, .12), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=20*(0:5), cex.axis=1.66)# ; axis(side=1, at=20*(0:5), cex.axis=1.66) 
axis(side=2, at=.02*(0:6), las=1, cex.axis=1.66)
with(df_meanfit_42, lines(day, sd.fitness, col=col_h, lty=6, lwd=2))

# dev.off(); system(paste0("open ~/Dropbox/Talks_Algae_Virus/First\\ paper/figures/extfig_fitness_sd.pdf"))

## include two more covariables to see if there is a difference between
## the "sweep" timepoints and the rest
df_contemp_22$sweep.occurred <- 
  with(df_contemp_22, ifelse(host.day %in% c(27, 57), TRUE, FALSE))
df_contemp_32$sweep.occurred <- 
  with(df_contemp_32, ifelse(host.day %in% c(27, 57), TRUE, FALSE))
df_contemp_42$sweep.occurred <- 
  with(df_contemp_42, ifelse(host.day %in% c(27, 57), TRUE, FALSE))

plot(growth_rate ~ sweep.occurred, type='p', col=2, data=df_contemp_22)
points(growth_rate ~ sweep.occurred, col=3, data=df_contemp_32)
points(growth_rate ~ sweep.occurred, col=4, data=df_contemp_42)

summary(lm(growth_rate ~ sweep.occurred, data=df_contemp_22))
summary(lm(growth_rate ~ sweep.occurred, data=df_contemp_32))
summary(lm(growth_rate ~ sweep.occurred, data=df_contemp_42))
## not significant. If I look only at day 27, the growth rates are 
## significantly different, but the differences are not in the same direction. 
## FITNESS
#######################################

#######################################
## TRADEOFF BETWEEN RESISTANCE AND GROWTH
## calculate resistance range per individual
df_inf_22$host.indiv.f <- with(df_inf_22, factor(
  host.indiv, levels=unique(host.indiv)
))

df_inf_22$indiv.res <- with(df_inf_22, rep(tapply(
  resistant, INDEX=host.indiv.f, (function(x) sum(x, na.rm=T))
), table(host.indiv.f)))
## then extract rows of growth in absence of virus
df_troff_22 <- df_inf_22[df_inf_22$virus.ID == "neg", ]

df_inf_32$host.indiv.f <- with(df_inf_32, factor(
  host.indiv, levels=unique(host.indiv)
))

df_inf_32$indiv.res <- with(df_inf_32, rep(tapply(
  resistant, INDEX=host.indiv.f, (function(x) sum(x, na.rm=T))
), table(host.indiv.f)))
## then extract rows of growth in absence of virus
df_troff_32 <- df_inf_32[df_inf_32$virus.ID == "neg", ]

df_inf_42$host.indiv.f <- with(df_inf_42, factor(
  host.indiv, levels=unique(host.indiv)
))

df_inf_42$indiv.res <- with(df_inf_42, rep(tapply(
  resistant, INDEX=host.indiv.f, (function(x) sum(x, na.rm=T))
), table(host.indiv.f)))
## then extract rows of growth in absence of virus
df_troff_42 <- df_inf_42[df_inf_42$virus.ID == "neg", ]

with(df_troff_22, summary(lm(mean.growth ~ indiv.res)))
with(df_troff_32, summary(lm(mean.growth ~ indiv.res)))
with(df_troff_42, summary(lm(mean.growth ~ indiv.res)))


# pdf("~/Dropbox/Talks_Algae_Virus/First_paper/figures/extfig_tradeoff.pdf", 
#     width=8, height=4)
# pdf("~/Documents/HVInt/Exp2017/Exp2017_16/pheno/figs/extfig_tradeoff.pdf", width=8, height=4)
{
par(mfrow=c(1, 1), mai=c(0.42, 0.62, 0.12, 0.12))
  plot(c(0, 11), c(0, .4), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=2*(0:5), cex.axis=1.66)# ; axis(side=1, at=20*(0:5), cex.axis=1.66) 
axis(side=2, at=.1*(0:6), las=1, cex.axis=1.66)
with(df_troff_22, points(jitter(indiv.res, 0.5), mean.growth, col=col_h, pch=22))
with(df_troff_32, points(jitter(indiv.res, 0.5), mean.growth, col=col_h, pch=23))
with(df_troff_42, points(jitter(indiv.res, 0.5), mean.growth, col=col_h, pch=24))
with(df_troff_22, abline(lm(mean.growth ~ indiv.res), col=col_h, lwd=1.5, lty=2))
with(df_troff_32, abline(lm(mean.growth ~ indiv.res), col=col_h, lwd=1.5, lty=3))
with(df_troff_42, abline(lm(mean.growth ~ indiv.res), col=col_h, lwd=1.5, lty=4))
legend(x=9.4, y=0.42, legend=c("I", "II", "III"), col=col_h, xjust=0.5, 
       bty='n', pch=c(22, 23, 24), ncol=1)
legend(x=10.5, y=0.42, legend=c("I", "II", "III"), col=col_h, bty='o', 
       xjust=0.5, lty=c(2, 3, 4), ncol=1, box.col="white")
lines(c(8.75, 8.75), c(.31, .45)); lines(c(8.75, 11.5), c(.31, .31))
}
# dev.off(); system("open ~/Dropbox/Talks_Algae_Virus/First_paper/figures/extfig_tradeoff.pdf")
## 
#######################################



#####################################
## HEATMAP PER TIME POINT
col_v2h <- colorRampPalette(colors=c("#E69F00", "#009E73"), space="rgb")(12)
tps_v22 <- unique(df_inf_22$virus.day)[!unique(df_inf_22$virus.day)=="neg"] %>% as.integer
tps_h22 <- unique(df_inf_22$host.day)
## create the matrix; proud of this code
hm_tp_22 <- sapply(
  tps_v22, function(vtp) sapply(
    tps_h22, function(htp) mean(df_inf_22$resistant[
      (df_inf_22$virus.day == vtp) & (df_inf_22$host.day == htp)
      ], na.rm=TRUE)
  )
)
## make them recognizable
colnames(hm_tp_22) <- tps_v22
rownames(hm_tp_22) <- tps_h22

tps_v32 <- unique(df_inf_32$virus.day)[!unique(df_inf_32$virus.day)=="neg"] %>% as.integer
tps_h32 <- unique(df_inf_32$host.day)
## create the matrix; proud of this code
hm_tp_32 <- sapply(
  tps_v32, function(vtp) sapply(
    tps_h32, function(htp) mean(df_inf_32$resistant[
      (df_inf_32$virus.day == vtp) & (df_inf_32$host.day == htp)
      ], na.rm=TRUE)
  )
)
## make them recognizable
colnames(hm_tp_32) <- tps_v32
rownames(hm_tp_32) <- tps_h32

tps_v42 <- unique(df_inf_42$virus.day)[!unique(df_inf_42$virus.day)=="neg"] %>% as.integer
tps_h42 <- unique(df_inf_42$host.day)
## create the matrix; proud of this code
hm_tp_42 <- sapply(
  tps_v42, function(vtp) sapply(
    tps_h42, function(htp) mean(df_inf_42$resistant[
      (df_inf_42$virus.day == vtp) & (df_inf_42$host.day == htp)
      ], na.rm=TRUE)
  )
)
## make them recognizable
colnames(hm_tp_42) <- tps_v42
rownames(hm_tp_42) <- tps_h42

# pdf(paste0(exp2017_dir, "pheno/figs/fig_234dot2_hm_tps.pdf"), 
#     width=7, height=7)
par(mfrow=c(3, 1), mar=c(4, 4, 3, 2))
## FIRST plot: grey background
image(
  1:length(tps_h22), 1:length(tps_v22), axes=F, xlab="Host timepoint", 
  z=matrix(rep(0, length(tps_h22) * length(tps_v22)), nrow=length(tps_h22)), 
  ylab="Virus timepoint", col='dark grey', frame.plot=T, 
  main=paste0("II.2: Infection matrix (proportion of resistant clones)")
)
## overlay colour
image(
  1:length(tps_h22), 1:length(tps_v22), z=hm_tp_22, axes=F, add=TRUE, 
  col=col_v2h, zlim=c(0, 1)
)
abline(v=(-.5 + 1:length(tps_h22))); abline(h=(-.5 + 1:length(tps_v22)))
axis(side=1, at=1:length(tps_h22), cex.axis=1, font.axis=2, las=3, labels=tps_h22)
axis(side=2, at=1:length(tps_v22), cex.axis=1, font.axis=2, las=1, labels=tps_v22)
text(expand.grid(1:length(tps_h22), 1:length(tps_v22)),
     labels=round(hm_tp_22, 1))

image(
  1:length(tps_h32), 1:length(tps_v32), axes=F, xlab="Host timepoint", 
  z=matrix(rep(0, length(tps_h32) * length(tps_v32)), nrow=length(tps_h32)), 
  ylab="Virus timepoint", col='dark grey', frame.plot=T, 
  main=paste0("III.2: Infection matrix (proportion of resistant clones)")
)
## overlay colour
image(
  1:length(tps_h32), 1:length(tps_v32), z=hm_tp_32, axes=F, add=TRUE, 
  col=col_v2h, zlim=c(0, 1)
)
abline(v=(-.5 + 1:length(tps_h32))); abline(h=(-.5 + 1:length(tps_v32)))
axis(side=1, at=1:length(tps_h32), cex.axis=1, font.axis=2, las=3, labels=tps_h32)
axis(side=2, at=1:length(tps_v32), cex.axis=1, font.axis=2, las=1, labels=tps_v32)
text(expand.grid(1:length(tps_h32), 1:length(tps_v32)),
     labels=round(hm_tp_32, 1))


image(
  1:length(tps_h42), 1:length(tps_v42), axes=F, xlab="Host timepoint", 
  z=matrix(rep(0, length(tps_h42) * length(tps_v42)), nrow=length(tps_h42)), 
  ylab="Virus timepoint", col='dark grey', frame.plot=T, 
  main=paste0("IV.2: Infection matrix (proportion of resistant clones)")
)
## overlay colour
image(
  1:length(tps_h42), 1:length(tps_v42), z=hm_tp_42, axes=F, add=TRUE, 
  col=col_v2h, zlim=c(0, 1)
)
abline(v=(-.5 + 1:length(tps_h42))); abline(h=(-.5 + 1:length(tps_v42)))
axis(side=1, at=1:length(tps_h42), cex.axis=1, font.axis=2, las=3, labels=tps_h42)
axis(side=2, at=1:length(tps_v42), cex.axis=1, font.axis=2, las=1, labels=tps_v42)
text(expand.grid(1:length(tps_h42), 1:length(tps_v42)),
     labels=round(hm_tp_42, 1))
# dev.off(); system(paste0("open ", exp2017_dir, "pheno/figs/fig_234dot2_hm_tps.pdf"))
## HEATMAP PER TIME POINT
#####################################



#####################################
## HEATMAP PER INDIVIDUAL
tps_v22 <- unique(df_inf_22$virus.day)[!unique(df_inf_22$virus.day)=="neg"] %>% as.integer
cls_h22 <- unique(df_inf_22$host.indiv)
## clones per day: 
ncls_22 <- strsplit(cls_h22, split="[.]") %>% (function(lis) sapply(lis, "[[", 1)) %>% table

hm_cl_22 <- sapply(
  tps_v22, function(vtp) sapply(
    cls_h22, function(hcl) sum(df_inf_22$resistant[
      (df_inf_22$virus.day == vtp) & (df_inf_22$host.indiv == hcl)
      ])
  )
)
## make them recognizable
colnames(hm_cl_22) <- tps_v22
rownames(hm_cl_22) <- cls_h22

tps_v32 <- unique(df_inf_32$virus.day)[!unique(df_inf_32$virus.day)=="neg"] %>% as.integer
cls_h32 <- unique(df_inf_32$host.indiv)
## clones per day: 
ncls_32 <- strsplit(cls_h32, split="[.]") %>% (function(lis) sapply(lis, "[[", 1)) %>% table

hm_cl_32 <- sapply(
  tps_v32, function(vtp) sapply(
    cls_h32, function(hcl) sum(df_inf_32$resistant[
      (df_inf_32$virus.day == vtp) & (df_inf_32$host.indiv == hcl)
      ])
  )
)
## make them recognizable
colnames(hm_cl_32) <- tps_v32
rownames(hm_cl_32) <- cls_h32

tps_v42 <- unique(df_inf_42$virus.day)[!unique(df_inf_42$virus.day)=="neg"] %>% as.integer
cls_h42 <- unique(df_inf_42$host.indiv)
## clones per day: 
ncls_42 <- strsplit(cls_h42, split="[.]") %>% (function(lis) sapply(lis, "[[", 1)) %>% table

hm_cl_42 <- sapply(
  tps_v42, function(vtp) sapply(
    cls_h42, function(hcl) sum(df_inf_42$resistant[
      (df_inf_42$virus.day == vtp) & (df_inf_42$host.indiv == hcl)
      ])
  )
)
## make them recognizable
colnames(hm_cl_42) <- tps_v42
rownames(hm_cl_42) <- cls_h42


# pdf(paste0(exp2017_dir, "pheno/figs/fig_234dot2_hm_inds.pdf"), 
#     width=7, height=7)
# par(mfrow=c(3, 1), mar=c(4, 4, 3, 2))
pdf("~/Dropbox/Talks_Algae_Virus/figures/234dot2_pheno_hm.pdf", width=5, height=4)
par(mai=c(0.42, 0.82, 0.12, 0.12))
image(
  1:length(cls_h22), 1:length(tps_v22), axes=F, xlab="Host individual", 
  z=matrix(rep(0, length(cls_h22) * length(tps_v22)), nrow=length(cls_h22)), 
  ylab="Virus timepoint", col='dark grey', frame.plot=T
  # main=paste0("II.2: Infection matrix")
)
## overlay colour
image(
  1:length(cls_h22), 1:length(tps_v22), z=hm_cl_22, axes=F, add=TRUE, 
  col=col_v2h, zlim=c(0, 1)
)
abline(v=-.5+(0:sum(ncls_22)), lwd=.3)
abline(v=.5 + cumsum(ncls_22), lwd=1.5)
abline(h=(-.5 + 1:length(tps_v22)), lwd=1.5)
axis(side=1, at=-5+10*(1:length(tps_h22)), cex.axis=1, font.axis=2, 
     las=3, labels=tps_h22)
axis(side=2, at=1:length(tps_v22), cex.axis=1, font.axis=2, las=1, labels=tps_v22)

image(
  1:length(cls_h32), 1:length(tps_v32), axes=F, xlab="Host individual", 
  z=matrix(rep(0, length(cls_h32) * length(tps_v32)), nrow=length(cls_h32)), 
  ylab="Virus timepoint", col='dark grey', frame.plot=T
  # main=paste0("III.2: Infection matrix")
)
## overlay colour
image(
  1:length(cls_h32), 1:length(tps_v32), z=hm_cl_32, axes=F, add=TRUE, 
  col=col_v2h, zlim=c(0, 1)
)
abline(v=-.5+(0:sum(ncls_32)), lwd=.3)
abline(v=.5 + cumsum(ncls_32), lwd=1.5)
abline(h=(-.5 + 1:length(tps_v32)), lwd=1.5)
axis(side=1, at=-5+10*(1:length(tps_h32)), cex.axis=1, font.axis=2, 
     las=3, labels=tps_h32)
axis(side=2, at=1:length(tps_v32), cex.axis=1, font.axis=2, las=1, labels=tps_v32)

image(
  1:length(cls_h42), 1:length(tps_v42), axes=F, xlab="Host individual", 
  z=matrix(rep(0, length(cls_h42) * length(tps_v42)), nrow=length(cls_h42)), 
  ylab="Virus timepoint", col='dark grey', frame.plot=T
  # main=paste0("IV.2: Infection matrix")
)
## overlay colour
image(
  1:length(cls_h42), 1:length(tps_v42), z=hm_cl_42, axes=F, add=TRUE, 
  col=col_v2h, zlim=c(0, 1)
)
abline(v=-.5+(0:sum(ncls_42)), lwd=.3)
abline(v=.5 + cumsum(ncls_42), lwd=1.5)
abline(h=(-.5 + 1:length(tps_v42)), lwd=1.5)
axis(side=1, at=-5+10*(1:length(tps_h42)), cex.axis=1, font.axis=2, 
     las=3, labels=tps_h42)
axis(side=2, at=1:length(tps_v42), cex.axis=1, font.axis=2, las=1, labels=tps_v42)
dev.off(); system("open ~/Dropbox/Talks_Algae_Virus/figures/234dot2_pheno_hm.pdf")
# dev.off(); system(paste0("open ", exp2017_dir, "pheno/figs/fig_234dot2_hm_inds.pdf"))
## HEATMAP PER INDIVIDUAL
#####################################

#####################################
## PHENOTYPIC RANGE (HOST RESISTANCE & VIRUS INFECTIVITY)
## create data frame with one row per host individual
df_rrind_22 <- with(df_inf_22, data.frame(
  host.indiv = levels(as.factor(host.indiv)), 
  host.day.n = strsplit(levels(as.factor(host.indiv)), split="[.]") %>%
    (function(x) sapply(x, "[[", 1)) %>% as.numeric, 
  host.day.f = strsplit(levels(as.factor(host.indiv)), split="[.]") %>%
    (function(x) sapply(x, "[[", 1)) %>% as.factor, 
  res.range.indiv = tapply(resistant, host.indiv, (function(x) mean(x, na.rm=TRUE))), 
  is.generalist = tapply(
    resistant, INDEX=host.indiv, (function(res) (mean(res, na.rm=TRUE) == 1))
  )
))
## create data frame with one row per host time point
df_rr_22 <- with(df_rrind_22, data.frame(
  host.day.n = unique(host.day.n), 
  host.day.f = levels(host.day.f), 
  res.range.tp = tapply(res.range.indiv, INDEX=host.day.f, (function(x) mean(x, na.rm=TRUE))), 
  freq.generalist = colSums(matrix(as.numeric(is.generalist), ncol=length(unique(host.day.f))))
))
## create data frame with one row per host time point; maximum observed 
## resistance range; removing one outlier at day 15. 
df_rmax_22 <- with(df_rrind_22[df_rrind_22$host.indiv != "15.2", ], data.frame(
  host.day.n = unique(host.day.n), 
  host.day.f = levels(host.day.f), 
  res.max.tp = tapply(res.range.indiv, INDEX=host.day.f, (function(x) max(x, na.rm=TRUE)))
))
df_rmax_22 <- df_rmax_22[order(df_rmax_22$host.day.n), ]
## order the data frames by time point; careful to only do this after
## every tapply has been executed
df_rrind_22 <- df_rrind_22[order(df_rrind_22$host.day.n), ]
df_rr_22 <- df_rr_22[order(df_rr_22$host.day.n), ]

## create data frame with one row per virus time point
df_ir_22 <- with(df_inf_22, data.frame(
  virus.day.f = levels(virus.day.f), 
  inf.range = tapply(resistant, virus.day.f, (function(x) mean(1-x, na.rm=TRUE)))
))
df_ir_22 <- df_ir_22[df_ir_22$virus.day.f != "neg", ]
df_ir_22$virus.day.n <- as.numeric(as.character(df_ir_22$virus.day.f))
## Make a data frame that combines irange and rrange
df_irr_22 <- data.frame(
  day=c(df_rr_22$host.day.n, df_ir_22$virus.day.n) %>% unique %>% sort
)
df_irr_22$resist.range <- df_rr_22$res.range.tp[match(df_irr_22$day, df_rr_22$host.day.n)]
df_irr_22$infect.range <- df_ir_22$inf.range[match(df_irr_22$day, df_ir_22$virus.day.n)]
## now create a new match - mismatch variable. 
df_irr_22$match.range <- with(df_irr_22, infect.range > resist.range)
## now extrapolate the lines between the time points
df_irrext_22 <- data.frame(
  day = 1:99, 
  resist.range = with(df_irr_22, approx(x=day, y=resist.range, xout=1:99))$y, 
  infect.range = with(df_irr_22, approx(x=day, y=infect.range, xout=1:99))$y
)
## and for these, calculate if resistance is higher than infectivity
df_irrext_22$match <- with(df_irrext_22, resist.range <= infect.range)

## create data frame with one row per host individual
df_rrind_32 <- with(df_inf_32, data.frame(
  host.indiv = levels(as.factor(host.indiv)), 
  host.day.n = strsplit(levels(as.factor(host.indiv)), split="[.]") %>%
    (function(x) sapply(x, "[[", 1)) %>% as.numeric, 
  host.day.f = strsplit(levels(as.factor(host.indiv)), split="[.]") %>%
    (function(x) sapply(x, "[[", 1)) %>% as.factor, 
  res.range.indiv = tapply(resistant, host.indiv, (function(x) mean(x, na.rm=TRUE))), 
  is.generalist = tapply(
    resistant, INDEX=host.indiv, (function(res) (mean(res, na.rm=TRUE) == 1))
  )
))
## create data frame with one row per host time point
df_rr_32 <- with(df_rrind_32, data.frame(
  host.day.n = unique(host.day.n), 
  host.day.f = levels(host.day.f), 
  res.range.tp = tapply(res.range.indiv, INDEX=host.day.f, (function(x) mean(x, na.rm=TRUE))), 
  freq.generalist = colSums(matrix(as.numeric(is.generalist), ncol=length(unique(host.day.f))))
))
## create data frame with one row per host time point; maximum observed 
## resistance range; removing one outlier at day 15. 
df_rmax_32 <- with(df_rrind_32, data.frame(
  host.day.n = unique(host.day.n), 
  host.day.f = levels(host.day.f), 
  res.max.tp = tapply(res.range.indiv, INDEX=host.day.f, (function(x) max(x, na.rm=TRUE)))
))
df_rmax_32 <- df_rmax_32[order(df_rmax_32$host.day.n), ]
## create data frame with one row per virus time point
df_ir_32 <- with(df_inf_32, data.frame(
  virus.day.f = levels(virus.day.f), 
  inf.range = tapply(resistant, virus.day.f, (function(x) mean(1-x, na.rm=TRUE)))
))
df_ir_32 <- df_ir_32[df_ir_32$virus.day.f != "neg", ]
df_ir_32$virus.day.n <- as.numeric(as.character(df_ir_32$virus.day.f))
## order the data frames by time point
df_rrind_32 <- df_rrind_32[order(df_rrind_32$host.day.n), ]
df_rr_32 <- df_rr_32[order(df_rr_32$host.day.n), ]
## Make a data frame that combines irange and rrange
df_irr_32 <- data.frame(
  day=c(df_rr_32$host.day.n, df_ir_32$virus.day.n) %>% unique %>% sort
)
df_irr_32$resist.range <- df_rr_32$res.range.tp[match(df_irr_32$day, df_rr_32$host.day.n)]
df_irr_32$infect.range <- df_ir_32$inf.range[match(df_irr_32$day, df_ir_32$virus.day.n)]
## now create a new match - mismatch variable. 
df_irr_32$match.range <- with(df_irr_32, infect.range > resist.range)
## now extrapolate the lines between the time points
df_irrext_32 <- data.frame(
  day = 1:99, 
  resist.range = with(df_irr_32, approx(x=day, y=resist.range, xout=1:99))$y, 
  infect.range = with(df_irr_32, approx(x=day, y=infect.range, xout=1:99))$y
)
## and for these, calculate if resistance is higher than infectivity
df_irrext_32$match <- with(df_irrext_32, resist.range <= infect.range)


## create data frame with one row per host individual
df_rrind_42 <- with(df_inf_42, data.frame(
  host.indiv = levels(as.factor(host.indiv)), 
  host.day.n = strsplit(levels(as.factor(host.indiv)), split="[.]") %>%
    (function(x) sapply(x, "[[", 1)) %>% as.numeric, 
  host.day.f = strsplit(levels(as.factor(host.indiv)), split="[.]") %>%
    (function(x) sapply(x, "[[", 1)) %>% as.factor, 
  res.range.indiv = tapply(resistant, host.indiv, (function(x) mean(x, na.rm=TRUE))), 
  is.generalist = tapply(
    resistant, INDEX=host.indiv, (function(res) (mean(res, na.rm=TRUE) == 1))
  )
))
## create data frame with one row per host time point
df_rr_42 <- with(df_rrind_42, data.frame(
  host.day.n = unique(host.day.n), 
  host.day.f = levels(host.day.f), 
  res.range.tp = tapply(res.range.indiv, INDEX=host.day.f, (function(x) mean(x, na.rm=TRUE))), 
  freq.generalist = colSums(matrix(as.numeric(is.generalist), ncol=length(unique(host.day.f))))
))
## create data frame with one row per host time point; maximum observed 
## resistance range; removing one outlier at day 15. 
df_rmax_42 <- with(
  df_rrind_42[!(df_rrind_42$host.indiv %in% c("12.4", "15.10")), ], data.frame(
    host.day.n = unique(host.day.n), 
    host.day.f = levels(host.day.f), 
    res.max.tp = tapply(res.range.indiv, INDEX=host.day.f, (function(x) max(x, na.rm=TRUE)))
  )
)
df_rmax_42 <- df_rmax_42[order(df_rmax_42$host.day.n), ]
## create data frame with one row per virus time point
df_ir_42 <- with(df_inf_42, data.frame(
  virus.day.f = levels(virus.day.f), 
  inf.range = tapply(resistant, virus.day.f, (function(x) mean(1-x, na.rm=TRUE)))
))
df_ir_42 <- df_ir_42[df_ir_42$virus.day.f != "neg", ]
df_ir_42$virus.day.n <- as.numeric(as.character(df_ir_42$virus.day.f))
## order the data frames by time point
df_rrind_42 <- df_rrind_42[order(df_rrind_42$host.day.n), ]
df_rr_42 <- df_rr_42[order(df_rr_42$host.day.n), ]
## Make a data frame that combines irange and rrange
df_irr_42 <- data.frame(
  day=c(df_rr_42$host.day.n, df_ir_42$virus.day.n) %>% unique %>% sort
)
df_irr_42$resist.range <- df_rr_42$res.range.tp[match(df_irr_42$day, df_rr_42$host.day.n)]
df_irr_42$infect.range <- df_ir_42$inf.range[match(df_irr_42$day, df_ir_42$virus.day.n)]
## now create a new match - mismatch variable. 
df_irr_42$match.range <- with(df_irr_42, infect.range > resist.range)
## now extrapolate the lines between the time points
df_irrext_42 <- data.frame(
  day = 1:99, 
  resist.range = with(df_irr_42, approx(x=day, y=resist.range, xout=1:99))$y, 
  infect.range = with(df_irr_42, approx(x=day, y=infect.range, xout=1:99))$y
)
## and for these, calculate if resistance is higher than infectivity
df_irrext_42$match <- with(df_irrext_42, resist.range <= infect.range)


## plot

# pdf(paste0(exp2017_dir, "pheno/figs/fig_234dot2_phenot_range_tp.pdf"), width=9, height=3.4)
par(mfrow=c(1, 3), mai=c(.12, .12, .12, .12))
plot(c(0, 100), c(0, 1), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=2, at=.2*(0:5), labels=F)
axis(side=1, at=10*(0:10), labels=F)
with(df_ir_22, points(virus.day.n, inf.range, pch=16, cex=2, col=col_v))
with(df_ir_22, lines(virus.day.n, inf.range, lty=2, lwd=2, col=col_v))
with(df_rr_22, points(host.day.n, res.range.tp, pch=16, cex=2, col=col_h))
with(df_rr_22, lines(host.day.n, res.range.tp, lty=2, lwd=2, col=col_h))

plot(c(0, 100), c(0, 1), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=2, at=.2*(0:5), labels=F)
axis(side=1, at=10*(0:10), labels=F)
with(df_ir_32, points(virus.day.n, inf.range, pch=16, cex=2, col=col_v))
with(df_ir_32, lines(virus.day.n, inf.range, lty=2, lwd=2, col=col_v))
with(df_rr_32, points(host.day.n, res.range.tp, pch=16, cex=2, col=col_h))
with(df_rr_32, lines(host.day.n, res.range.tp, lty=2, lwd=2, col=col_h))

plot(c(0, 100), c(0, 1), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=2, at=.2*(0:5), labels=F)
axis(side=1, at=10*(0:10), labels=F)
with(df_ir_42, points(virus.day.n, inf.range, pch=16, cex=2, col=col_v))
with(df_ir_42, lines(virus.day.n, inf.range, lty=2, lwd=2, col=col_v))
with(df_rr_42, points(host.day.n, res.range.tp, pch=16, cex=2, col=col_h))
with(df_rr_42, lines(host.day.n, res.range.tp, lty=2, lwd=2, col=col_h))

# dev.off(); system(paste0("open ", exp2017_dir, "pheno/figs/fig_234dot2_phenot_range_tp.pdf"))


# pdf(paste0(exp2017_dir, "pheno/figs/fig_234dot2_phenot_range.pdf"))
# par(mfrow=c(3, 1), mar=c(4, 4, 3, 2))
# with(df_rmax_22, plot(
#   host.day.n, res.max.tp, type='s', lwd=2, col=col_h, xlim=c(0, 100), 
#   ylim=c(0, 1), xlab="Time (days)", ylab="Phenotype range", main="II.2"
# ))
# with(df_rrind_22, points(
#   jitter(host.day.n, 5), res.range.indiv, pch=16, cex=2, col=col_h
# ))
# par(new=TRUE)
# with(df_ir_22, plot(
#   virus.day.n, inf.range, type='s', lwd=2, col=col_v, xlim=c(0, 100), 
#   ylim=c(0, 1), xlab="Time (days)", ylab="Phenotype range", main="II.2"
# ))
# with(df_ir_22, points(
#   virus.day.n, inf.range, pch=16, cex=2.5, col=col_v
# ))
# 
# with(df_rmax_32, plot(
#   host.day.n, res.max.tp, type='s', lwd=2, col=col_h, xlim=c(0, 100), 
#   ylim=c(0, 1), xlab="Time (days)", ylab="Phenotype range", main="III.2"
# ))
# with(df_rrind_32, points(
#   jitter(host.day.n, 5), res.range.indiv, pch=16, cex=2, col=col_h
# ))
# par(new=TRUE)
# with(df_ir_32, plot(
#   virus.day.n, inf.range, type='s', lwd=2, col=col_v, xlim=c(0, 100), 
#   ylim=c(0, 1), xlab="Time (days)", ylab="Phenotype range", main="III.2"
# ))
# with(df_ir_32, points(
#   virus.day.n, inf.range, pch=16, cex=2.5, col=col_v
# ))
# 
# with(df_rmax_42, plot(
#   host.day.n, res.max.tp, type='s', lwd=2, col=col_h, xlim=c(0, 100), 
#   ylim=c(0, 1), xlab="Time (days)", ylab="Phenotype range", main="IV.2"
# ))
# with(df_rrind_42, points(
#   jitter(host.day.n, 5), res.range.indiv, pch=16, cex=2, col=col_h
# ))
# par(new=TRUE)
# with(df_ir_42, plot(
#   virus.day.n, inf.range, type='s', lwd=2, col=col_v, xlim=c(0, 100), 
#   ylim=c(0, 1), xlab="Time (days)", ylab="Phenotype range", main="IV.2"
# ))
# with(df_ir_42, points(
#   virus.day.n, inf.range, pch=16, cex=2.5, col=col_v
# ))
# dev.off(); system(paste0("open ", exp2017_dir, "pheno/figs/fig_234dot2_phenot_range.pdf"))
col_mm <- c(col_h, col_v)

# pdf(paste0(exp2017_dir, "pheno/figs/fig_234dot2_matchmismatch_extrapol.pdf"), 
#     width=9, height=.6)
par(mfrow=c(1, 3), mai=c(0.12, 0.12, 0.12, 0.12))
plot(c(0, 100), c(-.05, .05), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=10*(0:10), labels=F)
with(df_irrext_22, points(day, rep(0, length(day)), col=col_mm[1+as.numeric(match)], pch=15, cex=2))
with(subset(df_irrext_22, match & !is.na(match)), 
     points(day, rep(0, length(day)), col=col_mm[2], pch=15, cex=2))

plot(c(0, 100), c(-.05, .05), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=10*(0:10), labels=F)
with(df_irrext_32, points(day, rep(0, length(day)), col=col_mm[1+as.numeric(match)], pch=15, cex=2))
with(subset(df_irrext_32, match & !is.na(match)), 
     points(day, rep(0, length(day)), col=col_mm[2], pch=15, cex=2))

plot(c(0, 100), c(-.05, .05), xaxt='n', yaxt='n', type='n', ylab="", xlab="")
axis(side=1, at=10*(0:10), labels=F)
with(df_irrext_42, points(day, rep(0, length(day)), col=col_mm[1+as.numeric(match)], pch=15, cex=2))
with(subset(df_irrext_42, match & !is.na(match)), 
     points(day, rep(0, length(day)), col=col_mm[2], pch=15, cex=2))

# dev.off(); system(paste0("open ", exp2017_dir, "pheno/figs/fig_234dot2_matchmismatch_extrapol.pdf"))
## mismatch becomes match at one day, i.e. at day 45

## PHENOTYPIC RANGE (HOST RESISTANCE & VIRUS INFECTIVITY)
#####################################


#####################################
## ARD VS FSD PHASE
## Bar plots of resistance versus past:contemporary:future virus, 
## separated in time points before and after generalist. 
head(df_inf_22)
with(df_inf_22, table(virus.day, host.day))
## virus.day is a factor because of "neg"; extract actual resistance assays, 
## transform virus.day to numeric
df_phase_22 <- subset(df_inf_22, virus.day!="neg")[
  , c("virus.day", "host.day", "resistant", "virus.day.f", "host.day.f", "host.indiv.f")
]
df_phase_22$virus.day <- as.numeric(df_phase_22$virus.day)
## add phase
df_phase_22$phase <- with(df_phase_22, factor(
  1 + (host.day > 57), labels=c("ARD", "FSD")
))

## now add another factor reflecting comparison type
df_phase_22$comparison <- with(df_phase_22, factor(
  1 + (virus.day == host.day) + 2 * (virus.day > host.day), 
  labels = c("PAST", "CONTEMP", "FUTURE")
))
## and summarise in a table
df_respf_22 <- data.frame(
  phase = factor(rep(1:2, each=3), labels=c("ARD", "FSD")), 
  comparison = factor(rep(1:3, 2), labels=c("PAST", "CONTEMP", "FUTURE")), 
  resistance = c(with(subset(df_phase_22, phase == "ARD"), tapply(
    resistant, INDEX=comparison, (function(x) mean(x, na.rm=TRUE))
    )), 
    with(subset(df_phase_22, phase == "FSD"), tapply(
    resistant, INDEX=comparison, (function(x) mean(x, na.rm=TRUE))
    ))), 
  sterr = c(with(subset(df_phase_22, phase == "ARD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(
      mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE) / sum(!is.na(x))
    ))
  )), with(subset(df_phase_22, phase == "FSD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(
      mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE) / sum(!is.na(x))
    ))
  ))), 
  stdv = c(with(subset(df_phase_22, phase == "ARD"), tapply(
    resistant, INDEX=comparison, 
     (function(x) sqrt(mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE)))
  )), with(subset(df_phase_22, phase == "FSD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE)))
  )))
)


df_phase_32 <- subset(df_inf_32, virus.day!="neg")[
  , c("virus.day", "host.day", "resistant", "virus.day.f", "host.day.f", "host.indiv.f")
  ]
df_phase_32$virus.day <- as.numeric(df_phase_32$virus.day)
## add phase
df_phase_32$phase <- with(df_phase_32, factor(
  1 + (host.day > 57), labels=c("ARD", "FSD")
))

## now add another factor reflecting comparison type
df_phase_32$comparison <- with(df_phase_32, factor(
  1 + (virus.day == host.day) + 2 * (virus.day > host.day), 
  labels = c("PAST", "CONTEMP", "FUTURE")
))
## and summarise in a table
df_respf_32 <- data.frame(
  phase = factor(rep(1:2, each=3), labels=c("ARD", "FSD")), 
  comparison = factor(rep(1:3, 2), labels=c("PAST", "CONTEMP", "FUTURE")), 
  resistance = c(with(subset(df_phase_32, phase == "ARD"), tapply(
    resistant, INDEX=comparison, (function(x) mean(x, na.rm=TRUE))
  )), 
  with(subset(df_phase_32, phase == "FSD"), tapply(
    resistant, INDEX=comparison, (function(x) mean(x, na.rm=TRUE))
  ))), 
  sterr = c(with(subset(df_phase_32, phase == "ARD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(
      mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE) / sum(!is.na(x))
    ))
  )), with(subset(df_phase_32, phase == "FSD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(
      mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE) / sum(!is.na(x))
    ))
  ))), 
  stdv = c(with(subset(df_phase_32, phase == "ARD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE)))
  )), with(subset(df_phase_32, phase == "FSD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE)))
  )))
)

df_phase_42 <- subset(df_inf_42, virus.day!="neg")[
  , c("virus.day", "host.day", "resistant", "virus.day.f", "host.day.f", "host.indiv.f")
  ]
df_phase_42$virus.day <- as.numeric(df_phase_42$virus.day)
## add phase
df_phase_42$phase <- with(df_phase_42, factor(
  1 + (host.day > 57), labels=c("ARD", "FSD")
))

## now add another factor reflecting comparison type
df_phase_42$comparison <- with(df_phase_42, factor(
  1 + (virus.day == host.day) + 2 * (virus.day > host.day), 
  labels = c("PAST", "CONTEMP", "FUTURE")
))

## and summarise in a table
df_respf_42 <- data.frame(
  phase = factor(rep(1:2, each=3), labels=c("ARD", "FSD")), 
  comparison = factor(rep(1:3, 2), labels=c("PAST", "CONTEMP", "FUTURE")), 
  resistance = c(with(subset(df_phase_42, phase == "ARD"), tapply(
    resistant, INDEX=comparison, (function(x) mean(x, na.rm=TRUE))
  )), 
  with(subset(df_phase_42, phase == "FSD"), tapply(
    resistant, INDEX=comparison, (function(x) mean(x, na.rm=TRUE))
  ))), 
  sterr = c(with(subset(df_phase_42, phase == "ARD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(
      mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE) / sum(!is.na(x))
    ))
  )), with(subset(df_phase_42, phase == "FSD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(
      mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE) / sum(!is.na(x))
    ))
  ))), 
  stdv = c(with(subset(df_phase_42, phase == "ARD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE)))
  )), with(subset(df_phase_42, phase == "FSD"), tapply(
    resistant, INDEX=comparison, 
    (function(x) sqrt(mean(x == 0, na.rm=TRUE) * mean(x == 1, na.rm=TRUE)))
  )))
)

# ggplot(data=df_respf_22, aes(x=phase, y=resistance, fill=comparison)) + 
#   geom_bar(stat='identity', position=position_dodge()) + ylim(-2, 3) + 
#   geom_errorbar(
#     aes(ymin=resistance - stdv, ymax=resistance + stdv), width=.1, 
#     position=position_dodge(.9)
#   )
# brewer.pal(13, name="Greys")
col_gr <- c("#F7F7F7", "#969696", "#252525")

fig_phases_22 <- ggplot(
  data=df_respf_22, aes(x=phase, y=resistance, fill=comparison)
) + geom_bar(stat='identity', position=position_dodge(), 
             col=I('black'), size=I(.8)) + 
  geom_errorbar(
    aes(ymin=resistance - sterr, ymax=resistance + sterr), width=.2, size=.6, 
    position=position_dodge(.9)
  ) + 
  scale_colour_manual(values=col_gr) + 
  scale_fill_manual(name = "Comparison", values=col_gr) + 
  guides(fill=FALSE, colour=FALSE) + 
  scale_y_continuous(breaks=.2*(0:5), limits=c(0, 1.1)) + 
  labs(x="", y="") + theme_bw() + theme(
    plot.margin = unit(c(.2 ,.1 ,.2 ,0 ), "cm"), 
    panel.border = element_rect(fill=NA, colour = "black", size=1.2), 
    panel.grid=element_blank(), 
    # axis.text.y=element_text(size=24), axis.ticks.y=element_line(size=1), 
    axis.text.y=element_blank(), axis.ticks.y=element_line(size=1), 
    axis.text.x=element_blank(), axis.ticks.x=element_blank()
  )

fig_phases_32 <- ggplot(
  data=df_respf_32, aes(x=phase, y=resistance, fill=comparison)
) + geom_bar(stat='identity', position=position_dodge(), 
             col=I('black'), size=I(.8)) + 
  geom_errorbar(
    aes(ymin=resistance - sterr, ymax=resistance + sterr), width=.2, size=.6, 
    position=position_dodge(.9)
  ) + 
  scale_colour_manual(values=col_gr) + 
  scale_fill_manual(name = "Comparison", values=col_gr) + 
  guides(fill=FALSE, colour=FALSE) + 
  scale_y_continuous(breaks=.2*(0:5), limits=c(0, 1.1)) + 
  labs(x="", y="") + theme_bw() + theme(
    plot.margin = unit(c(.2 ,.1 ,.2 ,0 ), "cm"), 
    panel.border = element_rect(fill=NA, colour = "black", size=1.2), 
    panel.grid=element_blank(), 
    # axis.text.y=element_text(size=24), axis.ticks.y=element_line(size=1), 
    axis.text.y=element_blank(), axis.ticks.y=element_line(size=1), 
    axis.text.x=element_blank(), axis.ticks.x=element_blank()
  )

fig_phases_42 <- ggplot(
  data=df_respf_42, aes(x=phase, y=resistance, fill=comparison)
) + geom_bar(stat='identity', position=position_dodge(), 
             col=I('black'), size=I(.8)) + 
  geom_errorbar(
    aes(ymin=resistance - sterr, ymax=resistance + sterr), width=.2, size=.6, 
    position=position_dodge(.9)
  ) + 
  scale_colour_manual(values=col_gr) + 
  scale_fill_manual(name = "Comparison", values=col_gr) + 
  guides(fill=FALSE, colour=FALSE) + 
  scale_y_continuous(breaks=.2*(0:5), limits=c(0, 1.1)) + 
  labs(x="", y="") + theme_bw() + theme(
    plot.margin = unit(c(.2 ,.1 ,.2 ,0 ), "cm"), 
    panel.border = element_rect(fill=NA, colour = "black", size=1.2), 
    panel.grid=element_blank(), 
    # axis.text.y=element_text(size=24), axis.ticks.y=element_line(size=1), 
    axis.text.y=element_blank(), axis.ticks.y=element_line(size=1), 
    axis.text.x=element_blank(), axis.ticks.x=element_blank()
  )

# pdf(paste0(exp2017_dir, "pheno/figs/fig_234dot2_phases.pdf"), width=12, height=4)
gg_multiplot(
  fig_phases_22, fig_phases_32, fig_phases_42, cols=3
)
# dev.off(); system(paste0("open ", exp2017_dir, "pheno/figs/fig_234dot2_phases.pdf"))

with(within(
  subset(df_phase_22, phase=="ARD"), 
  comparison <- relevel(comparison, ref="CONTEMP")
), summary(glm(resistant ~ comparison, family='binomial')))
with(within(
  subset(df_phase_32, phase=="ARD"), 
  comparison <- relevel(comparison, ref="CONTEMP")
), summary(glm(resistant ~ comparison, family='binomial')))
with(within(
  subset(df_phase_42, phase=="ARD"), 
  comparison <- relevel(comparison, ref="CONTEMP")
), summary(glm(resistant ~ comparison, family='binomial')))

with(within(
  subset(df_phase_22, phase=="FSD"), 
  comparison <- relevel(comparison, ref="CONTEMP")
), summary(glm(resistant ~ comparison, family='binomial')))
with(within(
  subset(df_phase_32, phase=="FSD"), 
  comparison <- relevel(comparison, ref="CONTEMP")
), summary(glm(resistant ~ comparison, family='binomial')))
with(within(
  subset(df_phase_42, phase=="FSD"), 
  comparison <- relevel(comparison, ref="CONTEMP")
), summary(glm(resistant ~ comparison, family='binomial')))

## ARD VS FSD PHASE
#####################################

## write csv's for Carlos..
# write.csv(
#   hm_cl_22, quote=F, 
#   file=paste0(exp2017_dir, "pheno/repl2_infectionmat_hostindivs.csv")
# )
# write.csv(
#   hm_cl_32, quote=F, 
#   file=paste0(exp2017_dir, "pheno/repl3_infectionmat_hostindivs")
# )
# write.csv(
#   hm_cl_42, quote=F, 
#   file=paste0(exp2017_dir, "pheno/repl4_infectionmat_hostindivs.csv")
# )
# 
# write.csv(
#   df_inf_22[df_inf_22$virus.day != "neg", c("virus.day", "host.day", "host.indiv", "resistant")], 
#   quote=F, row.names=F, 
#   file=paste0(exp2017_dir, "pheno/repl2_infection_dataframe.csv")
# )
# write.csv(
#   df_inf_32[df_inf_32$virus.day != "neg", c("virus.day", "host.day", "host.indiv", "resistant")], 
#   quote=F, row.names=F, 
#   file=paste0(exp2017_dir, "pheno/repl3_infection_dataframe.csv")
# )
# write.csv(
#   df_inf_42[df_inf_42$virus.day != "neg", c("virus.day", "host.day", "host.indiv", "resistant")], 
#   quote=F, row.names=F, 
#   file=paste0(exp2017_dir, "pheno/repl4_infection_dataframe.csv")
# )

