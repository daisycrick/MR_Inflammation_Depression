library(TwoSampleMR)
library(tibble)
library(vroom)
library(data.table)
library(tidyverse)
#vignette("MVMR")

setwd("~/PhenoSPD")
file="trait_depression_tidy_outcome.txt.gz"
full_depression <- as_tibble(fread(file))

file="trait_BMI_tidy_outcome.txt.gz"
full_BMI <- as_tibble(fread(file))

setwd("/user/home/dc15053/PhenoSPD")

full_depression$outcome <- c("MDD")
full_depression$id.outcome <- c("MDD")


colnames(full_BMI)[7]<- "beta.exposure"
colnames(full_BMI)[8]<- "se.exposure"
colnames(full_BMI)[4]<- "effect_allele.exposure"
colnames(full_BMI)[5]<- "other_allele.exposure"
colnames(full_BMI)[6]<- "eaf.exposure"
colnames(full_BMI)[9]<- "pval.exposure"
full_BMI$exposure <- c("BMI")
full_BMI$id.exposure <- c("BMI")


glimpse(full_depression)
glimpse(full_BMI)

 dat <- harmonise_data(
 exposure_dat = full_BMI,
     outcome_dat = full_depression
)

write.table(dat,"/",sep=" ",row.names=FALSE)

