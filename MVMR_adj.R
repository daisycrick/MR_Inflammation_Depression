rm(list = ls())
setwd("C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/year two/SMFQ and GlycA/two sample") # Set the file path to your local directory

library(metafor)
library(plyr) 

library(meta) 
library(rmeta) 
library(TwoSampleMR)
library(MRInstruments)
library(MRPracticals)
library(MVMR)
library(TwoSampleMR)
library(MVMR)
library(tibble)
library(vroom)
library(data.table)
library(tidyverse)
library(dplyr)
#vignette("MVMR")


exposure_data<-extract_instruments("met-d-GlycA")
#SD
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "ieu-b-102")
H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
H_data$samplesize.exposure = 115078
H_data$samplesize.outcome = 500199
H_data<- power_prune(H_data, method = 1, dist.outcome = "binary")
mr_results_GlycA_dep<-mr(H_data, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_weighted_mode"))
mr_results_GlycA_dep
write.csv(mr_results_GlycA_dep, "mr_results_GlycA_dep.csv")


exposures<-mv_extract_exposures(c("met-d-GlycA","ukb-a-248") , clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = exposures$SNP, outcomes = "ieu-b-102")
H_data<-mv_harmonise_data(
  exposure_dat = exposures, 
  outcome_dat = outcome_data
)

F.data <- format_mvmr(BXGs = H_data$exposure_beta,
                      BYG = H_data$outcome_beta,
                      seBXGs = H_data$exposure_se,
                      seBYG = H_data$outcome_se,
                      RSID = rownames(H_data$exposure_beta))
if ( colnames(H_data$exposure_beta)[1] != "Exposure1") {
  colnames(F.data)=c("SNP", "betaYG", "sebetaYG", "betaX2", "betaX1", "sebetaX2", "sebetaX1")
  F.data <- F.data[, c(1:3,5,4,7,6)]
}
head(F.data)







## Step 3: Estimate causal effects

tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=2) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}

tidy_mvmr_output <- function(mvmr_res) {
  #  tidy up MVMR returned output
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate,
           se="Std. Error",
           pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

res <- ivw_mvmr(r_input = F.data)

res_2 <-mv_multiple(H_data)
result_2smr_GlycA_dep <- res_2$result %>%
  split_outcome() %>%
  separate(outcome, "outcome", sep="[(]") %>% 
  mutate(outcome=stringr::str_trim(outcome))%>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome) %>% 
  tidy_pvals()
result_2smr_GlycA_dep

write.csv(result_2smr_GlycA_dep, "result_2smr_GlycA_dep.csv")


mvmrcovmatrix <- cbind(c(1, -0.133620263291244),c(-0.133620263291244, 1))
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,6:7])
#if your exposures are from difference cohorts gencov =0
sres <- strength_mvmr(r_input = F.data, gencov = Xcovmat)



## Step 4: Test for horizontal pleiotropy using conventional Q-statistic estimation
#if your exposures are from difference cohorts gencov =0
pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)


write.csv(result_2smr_GlycA_dep, "result_2smr_GlycA_dep.csv")

dat<-read.table("dat.txt")

exposure_data<-extract_instruments("met-d-GlycA")
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = "ukb-a-248")
outcome_data_subset<- outcome_data[c("SNP", "outcome", "id.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome")]
outcome_data_subset <- outcome_data_subset %>% rename (exposure = outcome, id.exposure = id.outcome, effect_allele.exposure = effect_allele.outcome, other_allele.exposure = other_allele.outcome, beta.exposure = beta.outcome, eaf.exposure = eaf.outcome, se.exposure = se.outcome, pval.exposure = pval.outcome)
exposure_data_subset<- exposure_data[c("SNP", "exposure", "id.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")]
combined_subset<-rbind(outcome_data_subset, exposure_data_subset)


outcomes <- extract_outcome_data(
  snps = combined_subset$SNP, outcomes = "ieu-b-102")


H_data<-mv_harmonise_data(
  exposure_dat = combined_subset, 
  outcome_dat = outcomes
)
F.data <- format_mvmr(BXGs = H_data$exposure_beta,
                      BYG = H_data$outcome_beta,
                      seBXGs = H_data$exposure_se,
                      seBYG = H_data$outcome_se,
                      RSID = rownames(H_data$exposure_beta))
if ( colnames(H_data$exposure_beta)[1] != "Exposure1") {
  colnames(F.data)=c("SNP", "betaYG", "sebetaYG", "betaX2", "betaX1", "sebetaX2", "sebetaX1")
  F.data <- F.data[, c(1:3,5,4,7,6)]
}
head(F.data)

## Step 3: Estimate causal effects

tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=2) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}

tidy_mvmr_output <- function(mvmr_res) {
  #  tidy up MVMR returned output
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate,
           se="Std. Error",
           pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

res <- ivw_mvmr(r_input = F.data)

res <- ivw_mvmr(r_input = F.data)
res
res_2 <-mv_multiple(H_data)
result_2smr_GlycA_dep_adj <- res_2$result %>%
  split_outcome() %>%
  separate(outcome, "outcome", sep="[(]") %>% 
  mutate(outcome=stringr::str_trim(outcome))%>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome) %>% 
  tidy_pvals()
result_2smr_GlycA_dep_adj
write.csv(result_2smr_GlycA_dep_adj,"result_2smr_adj_GlycA_dep.csv")


mvmrcovmatrix <- cbind(c(1, -0.133620263291244),c(-0.133620263291244, 1))
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,6:7])
#if your exposures are from difference cohorts gencov =0
sres <- strength_mvmr(r_input = F.data, gencov = Xcovmat)



## Step 4: Test for horizontal pleiotropy using conventional Q-statistic estimation
#if your exposures are from difference cohorts gencov =0
pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat)



################################################



exposure_data<-extract_instruments("ieu-b-102")
#SD
outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP, outcomes = "met-d-GlycA")

H_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
H_data  <- power_prune(dat, method = 1, dist.outcome = "continuous")
mr_results_dep_G<-mr(H_data, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_weighted_mode"))
mr_results_dep_G
write.csv(mr_results_dep_G, "mr_results_dep_GlycA.csv")


exposures<-mv_extract_exposures(c("ieu-b-102","ukb-a-248") ,clump_r2 = 0.001, clump_kb = 10000, pval_threshold = 5e-08)
outcome_data <- extract_outcome_data(
  snps = exposures$SNP, outcomes = "met-d-GlycA")
H_data<-mv_harmonise_data(
  exposure_dat = exposures, 
  outcome_dat = outcome_data
)
F.data <- format_mvmr(BXGs = H_data$exposure_beta,
                      BYG = H_data$outcome_beta,
                      seBXGs = H_data$exposure_se,
                      seBYG = H_data$outcome_se,
                      RSID = rownames(H_data$exposure_beta))
if ( colnames(H_data$exposure_beta)[1] != "Exposure1") {
  colnames(F.data)=c("SNP", "betaYG", "sebetaYG", "betaX2", "betaX1", "sebetaX2", "sebetaX1")
  F.data <- F.data[, c(1:3,5,4,7,6)]
}
head(F.data)



## Step 3: Estimate causal effects

tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=2) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}

tidy_mvmr_output <- function(mvmr_res) {
  #  tidy up MVMR returned output
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate,
           se="Std. Error",
           pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

res <- ivw_mvmr(r_input = F.data)

res <- ivw_mvmr(r_input = F.data)
res
res_2 <-mv_multiple(H_data)
result_2smr_dep_G <- res_2$result %>%
  split_outcome() %>%
  separate(outcome, "outcome", sep="[(]") %>% 
  mutate(outcome=stringr::str_trim(outcome))%>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome) %>% 
  tidy_pvals()
result_2smr_dep_G
write.csv(result_2smr_dep_G, "result_2smr_dep_GlycA.csv")

mvmrcovmatrix <- cbind(c(1, 0.0548744039452972),c(0.0548744039452972, 1))
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,6:7])

#if your exposures are from difference cohorts gencov =0
sres <- strength_mvmr(r_input = F.data, gencov=Xcovmat)

## Step 4: Test for horizontal pleiotropy using conventional Q-statistic estimation
#if your exposures are from difference cohorts gencov =0
pres <- pleiotropy_mvmr(r_input = F.data, gencov=Xcovmat)




dat<-read.table("dat_dep.txt")

exposure_data<-extract_instruments("ieu-b-102")
outcome_data <- extract_outcome_data(snps = exposure_data$SNP, outcomes = "ukb-a-248")
outcome_data_subset<- outcome_data[c("SNP", "outcome", "id.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome")]
outcome_data_subset <- outcome_data_subset %>% rename (exposure = outcome, id.exposure = id.outcome, effect_allele.exposure = effect_allele.outcome, other_allele.exposure = other_allele.outcome, beta.exposure = beta.outcome, eaf.exposure = eaf.outcome, se.exposure = se.outcome, pval.exposure = pval.outcome)
exposure_data_subset<- exposure_data[c("SNP", "exposure", "id.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")]
combined_subset <- rbind(outcome_data_subset, exposure_data_subset)

outcome_data <- extract_outcome_data(
  snps = combined_subset$SNP, outcomes = "met-d-GlycA")


H_data<-mv_harmonise_data(
  exposure_dat = combined_subset, 
  outcome_dat = outcome_data
)

F.data <- format_mvmr(BXGs = H_data$exposure_beta,
                      BYG = H_data$outcome_beta,
                      seBXGs = H_data$exposure_se,
                      seBYG = H_data$outcome_se,
                      RSID = rownames(H_data$exposure_beta))
if ( colnames(H_data$exposure_beta)[1] != "Exposure1") {
  colnames(F.data)=c("SNP", "betaYG", "sebetaYG", "betaX2", "betaX1", "sebetaX2", "sebetaX1")
  F.data <- F.data[, c(1:3,5,4,7,6)]
}
head(F.data)

## Step 3: Estimate causal effects

tidy_pvals<-function(df){
  # round up output values and keep p-vals in scientific notation
  df %>% 
    mutate(pval= as.character(pval)) %>% 
    mutate_if(is.numeric, round, digits=2) %>% 
    mutate(pval=as.numeric(pval),
           pval=scales::scientific(pval, digits = 2),
           pval=as.numeric(pval))
}

tidy_mvmr_output <- function(mvmr_res) {
  #  tidy up MVMR returned output
  mvmr_res %>%
    as.data.frame() %>% 
    rownames_to_column("exposure") %>% 
    rename(b=Estimate,
           se="Std. Error",
           pval="Pr(>|t|)") %>% 
    select(-c(`t value`)) %>% 
    TwoSampleMR::generate_odds_ratios()
}

res <- ivw_mvmr(r_input = F.data)

res <- ivw_mvmr(r_input = F.data)
res
res_2 <-mv_multiple(H_data)
result_2smr_dep_G_adj <- res_2$result %>%
  split_outcome() %>%
  separate(outcome, "outcome", sep="[(]") %>% 
  mutate(outcome=stringr::str_trim(outcome))%>% 
  generate_odds_ratios() %>% 
  select(-id.exposure, -id.outcome) %>% 
  tidy_pvals()
result_2smr_dep_G_adj
write.csv(result_2smr_dep_G_adj, "result_2smr_adj_dep_GlycA.csv")

mvmrcovmatrix <- cbind(c(1, 0.0548744039452972),c(0.0548744039452972, 1))
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,6:7])

#if your exposures are from difference cohorts gencov =0
sres <- strength_mvmr(r_input = F.data, gencov=Xcovmat)

## Step 4: Test for horizontal pleiotropy using conventional Q-statistic estimation
#if your exposures are from difference cohorts gencov =0
pres <- pleiotropy_mvmr(r_input = F.data, gencov=Xcovmat)

