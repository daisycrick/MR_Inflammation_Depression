rm(list = ls())
setwd("")

library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)

#load in data
dat <- read.table("test.dat.txt", header=T, sep=" ")

#subset data needed
MDD_BMI <- dat[, c('SNP', 'effect_allele.exposure', 'other_allele.exposure', 'beta.exposure', 'se.exposure', 'beta.outcome', 'se.outcome')]

head(MDD_BMI)

MDD_BMI <- rename(MDD_BMI, c(" " = "SNP", "allele_0" ="other_allele.exposure", "allele_1" = "effect_allele.exposure", "trait1_b" = "beta.exposure", 
"trait1_se" = "se.exposure", "trait2_b" = "beta.outcome",  "trait2_se" = "se.outcome"))


write.table(MDD_BMI,"/MDD_BMI.txt",sep=" ",row.names=FALSE, quotes=F)

MDD_BMI <- read.table("MDD_BMI.txt", header=T, sep=" ")
duplicates <- read.table("duplicate_SNPs.txt", header=F, sep=" ")
nrow(MDD_BMI)

MDD_BMI_2=MDD_BMI[MDD_BMI$SNP %in% duplicates$V1 ,]
MDD_BMI_2

MDD_BMI_dups=MDD_BMI_2[!duplicated(MDD_BMI_2[c('SNP')]), ]
MDD_BMI_dups

MDD_BMI_2=MDD_BMI[!MDD_BMI$SNP %in% duplicates$V1 ,]
nrow(MDD_BMI_2)
MDD_BMI_3=rbind(MDD_BMI_2, MDD_BMI_dups)
nrow(MDD_BMI_3)

write.table(dat,"MDD_BMI.txt",sep=" ",row.names=FALSE, quote=FALSE)





