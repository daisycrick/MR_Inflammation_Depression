rm(list = ls())
setwd("C:/Users/dc15053/OneDrive - University of Bristol/Documents/
      PhD/year two/SMFQ and GlycA/two sample") # Set the file path to your local directory

install.packages("remotes")
library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(MRInstruments)
library(purrr)
remotes::install_github("WSpiller/RadialMR")
library(RadialMR)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(dplyr) 
ao <- available_outcomes()

#get instruments from CB GWAS on UKBB
GlycA <- extract_instruments("met-d-GlycA") #61;

GlycA <- clump_data(
  GlycA,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)



#get outcome from PGC consortium 
depression <- extract_outcome_data(snps=GlycA$SNP, outcomes="ieu-b-102",proxies = T)
#get outcome from SSGAC consortium
dep_symp <- extract_outcome_data(snps=GlycA$SNP, outcomes="ieu-a-1000",proxies = T)
#get outcome from SSGAC consortium

#Make sure all SNP effects in the same direction? 
depression[,c("SNP","effect_allele.outcome","beta.outcome","se.outcome","pval.outcome")] 
summary(depression$beta.outcome) 


# HARMONISE THE DATA
dat <- harmonise_data(GlycA, depression, action = 2)



#Create mr-keep and drop snps
mr_keep <-dat[, c('mr_keep','SNP')]
dat_2 <- merge(dat, mr_keep,by="SNP")
dat_3 <- subset(dat_2, mr_keep == "TRUE")

dat_3$mr_keep <- dat_3$mr_keep.x

dat_3 <- dat_3[-c(dat_3$mr_keep.y, )]


# Check data
head(dat)
dim(dat)





dat<-power_prune(dat,method=1,dist.outcome="binary")


# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_GlycA_depression <- mr(dat_3, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode", "mr_wald_ratio"))
mr_GlycA_depression
GlycA_depression<-cbind.data.frame(mr_GlycA_depression$outcome,mr_GlycA_depression$nsnp,mr_GlycA_depression$method,mr_GlycA_depression$b,mr_GlycA_depression$se,mr_GlycA_depression$pval)
mr_GlycA_depression <- generate_odds_ratios(mr_GlycA_depression)
#Export results
write.csv(mr_GlycA_depression,"./GLYCA_depression.csv")
# Estimate odds ratio and 95% confidence interval
exp(mr_GlycA_depression$b[1])
exp(mr_GlycA_depression$b[1]-1.96*mr_GlycA_depression$se[1])
exp(mr_GlycA_depression$b[1]+1.96*mr_GlycA_depression$se[1])

exp(mr_GlycA_depression$b[2])
exp(mr_GlycA_depression$b[2]-1.96*mr_GlycA_depression$se[2])
exp(mr_GlycA_depression$b[2]+1.96*mr_GlycA_depression$se[2])

exp(mr_GlycA_depression$b[3])
exp(mr_GlycA_depression$b[3]-1.96*mr_GlycA_depression$se[3])
exp(mr_GlycA_depression$b[3]+1.96*mr_GlycA_depression$se[3])

exp(mr_GlycA_depression$b[4])
exp(mr_GlycA_depression$b[4]-1.96*mr_GlycA_depression$se[4])
exp(mr_GlycA_depression$b[4]+1.96*mr_GlycA_depression$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat) #Q should be large and p-values are relatively small (<0.05)

mr_pleiotropy_test(dat)

# Single SNP analysis - performs multiple analyses for each exposure-outcome combination using a different SNP each time
single_GlycA <- mr_singlesnp(dat)
single_GlycA #(p-value should be large and the intercept is close to 0

# Get R^2 and Fstat

dat$samplesize.exposure = 115078
dat$samplesize.outcome = 500199
dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
r2=directionality_test(dat)
exposure_r2=r2[1,"snp_r2.exposure"]

Samplesize=dat$samplesize.exposure[1]
nsnp=as.numeric(mr_GlycA_depression$nsnp[1])
f_stat= ((Samplesize-(nsnp-1))/nsnp)*(exposure_r2/(1-exposure_r2))
f_stat=cbind(f_stat, exposure_r2)


#MR_PRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
'force=TRUE'
library(MRPRESSO)

run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)


dat$MAF <- 1-(dat$eaf.exposure)
dat$r2 <- 2*dat$MAF*(1-dat$MAF)*(dat$beta.exposure^2)

dat$F.stat <- (dat$r2 * (dat$samplesize.exposure-2))/(1-dat$r2)

median(dat$F.stat)
min(dat$F.stat)
max(dat$F.stat)

###########################################
#VISUALIZE THE CAUSAL EFFECT OF GlycA ON depression#
###########################################
# Generate a scatter plot comparing the different methods
png("./GlycA_dep_scatter.png")
mr_scatter_plot(mr_GlycA_depression, dat)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./GlycA_dep_forest.png")
single_GlycA_plot <- mr_forest_plot(single_GlycA)
single_GlycA_plot
dev.off()

single_GlycA %>%
  mutate(SNP = reorder(SNP, b)) %>%
  ggplot(aes(y=b, x=SNP, ymin=b-(1.96*se), ymax=b+(1.96*se))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +
  geom_errorbarh(xmin="", xmax="", height=.1) +
  coord_flip() +
  ylab("Effect estimate of individual GlycA SNPs on MDD") +
  xlab("SNP")+
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.7) +
  ggtitle("Effect of Individual GlycA SNPs on MDD") +
  theme_classic() 



# Generate a funnel plot to check asymmetry
png("./GlycA_dep_funnel.png")
mr_funnel_plot(single_GlycA)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png("./GlycA_dep_loo.png")
mr_leaveoneout_plot(res_loo)
dev.off()

res_loo %>%
  mutate(SNP = reorder(SNP, b)) %>%
  ggplot(aes(y=b, x=SNP, ymin=b-(1.96*se), ymax=b+(1.96*se))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +
  geom_errorbarh(xmin="", xmax="", height=.1) +
  coord_flip() +
  ylab("Effect estimate of GlycA SNPs on MDD after exlusing individual SNP") +
  xlab("SNP")+
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.7) +
  ggtitle("Effect estimate with each SNP excluded in turn using GlycA SNPs on MDD") +
  theme_classic() 

#Radial MR
#???
jointdataset <- jointdataset[ , c ('SNP', 'exposure.x' , 'outcome', 'pos.exposure', 'chr.exposure', 
                                   'effect_allele.exposure', 'other_allele.exposure' , 'eaf.exposure' , 'b', 'se', 'p')]

write.csv(jointdataset, ".combined_data_GlycA_MDD.csv", row.names=FALSE)






res_single2 <- mr_singlesnp(dat)
res_single2

# Single SNP analysis - performs multiple analyses for each exposure-outcome combination using a different SNP each time


#HARMONISE THE DATA
dat2 <- harmonise_data(GlycA, dep_symp, action = 2)


mr_keep <-dat2[, c('mr_keep','SNP')]
dat2_2 <- merge(dat2, mr_keep,by="SNP")
dat2_3 <- subset(dat2_2, mr_keep == "TRUE")

dat2_3 <- dat2_3 %>% rename (mr_keep = mr_keep.x)
dat2_3 <- dat2_3[-c(mr_keep.y)]
dat2 <- dat2_3







# Check data
head(dat2)
dim(dat2)

# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_GLYcA_dep_symp <- mr(dat2_3, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_GLYcA_dep_symp
GLYcA_dep_symp<-cbind.data.frame(mr_GLYcA_dep_symp$outcome,mr_GLYcA_dep_symp$nsnp,mr_GLYcA_dep_symp$method,mr_GLYcA_dep_symp$b,mr_GLYcA_dep_symp$se,mr_GLYcA_dep_symp$pval)

#Export results
write.csv(GLYcA_dep_symp,"./GlyA_depressive_symptoms_results_.csv")
# Estimate odds ratio and 95% confidence interval
exp(mr_GLYcA_dep_symp$b[1])
exp(mr_GLYcA_dep_symp$b[1]-1.96*mr_GLYcA_dep_symp$se[1])
exp(mr_GLYcA_dep_symp$b[1]+1.96*mr_GLYcA_dep_symp$se[1])

exp(mr_GLYcA_dep_symp$b[2])
exp(mr_GLYcA_dep_symp$b[2]-1.96*mr_GLYcA_dep_symp$se[2])
exp(mr_GLYcA_dep_symp$b[2]+1.96*mr_GLYcA_dep_symp$se[2])

exp(mr_GLYcA_dep_symp$b[3])
exp(mr_GLYcA_dep_symp$b[3]-1.96*mr_GLYcA_dep_symp$se[3])
exp(mr_GLYcA_dep_symp$b[3]+1.96*mr_GLYcA_dep_symp$se[3])

exp(mr_GLYcA_dep_symp$b[4])
exp(mr_GLYcA_dep_symp$b[4]-1.96*mr_GLYcA_dep_symp$se[4])
exp(mr_GLYcA_dep_symp$b[4]+1.96*mr_GLYcA_dep_symp$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat2)

mr_pleiotropy_test(dat2)
res_single2 <- mr_singlesnp(dat2)
res_single2

# Single SNP analysis - performs multiple analyses for each exposure-outcome combination using a different SNP each time
single_GlycA_dep_symp <- mr_singlesnp(dat2)
single_GlycA


# Get F stat and R^2
dat2$samplesize.exposure = 115078
dat2$r.exposure <- get_r_from_pn(dat2$pval.exposure, dat2$samplesize.exposure)
r2_symp=directionality_test(dat2)
exposure_r2=r2_symp[1,"snp_r2.exposure"]

Samplesize_symp=dat2$samplesize.exposure[1]
nsnp_symp=as.numeric(mr_GLYcA_dep_symp$nsnp[1])
f_stat_symp= ((Samplesize_symp-(nsnp_symp-1))/nsnp_symp)*(exposure_r2/(1-exposure_r2))
f_stat_symp=cbind(f_stat_symp, exposure_r2)

#MR_PRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
'force=TRUE'
library(MRPRESSO)

run_mr_presso(dat2, NbDistribution = 1000, SignifThreshold = 0.05)

dat$MAF <- 1-(dat$eaf.exposure)
dat$r2 <- 2*dat$MAF*(1-dat$MAF)*(dat$beta.exposure^2)

dat$F.stat <- (dat$r2 * (dat$samplesize.exposure-2))/(1-dat$r2)



###########################################
#VISUALIZE THE CAUSAL EFFECT OF GlycA ON depressive#
###########################################
# Generate a scatter plot comparing the different methods
png("./GlycA_depressive_symp_scatter.png")
mr_scatter_plot(mr_GLYcA_dep_symp, dat2)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./GlycA_depressive~_sym_forest.png")
res_single2_plot<-mr_forest_plot(res_single2)
res_single2_plot
dev.off()



res_single2 %>%
  mutate(SNP = reorder(SNP, b)) %>%
  ggplot(aes(y=b, x=SNP, ymin=b-(1.96*se), ymax=b+(1.96*se))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +
  geom_errorbarh(xmin="", xmax="", height=.1) +
  coord_flip() +
  ylab("Effect estimate of individual GlycA SNP on depressive symptoms") +
  xlab("SNP")+
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.7) +
  ggtitle("Effect of Individual GlycA SNPs on depressive symptoms") +
  theme_classic() 



# Generate a funnel plot to check asymmetry
png("./GlycA_depressive_sym_funnel.png")
mr_funnel_plot(res_single2)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo2 <- mr_leaveoneout(dat2)
png("./GlycA_depressive_sym_loo.png")
GlycA_depress_symp_plot <- mr_leaveoneout_plot(res_loo2)
dev.off()

res_loo2 %>%
  mutate(SNP = reorder(SNP, b)) %>%
  ggplot(aes(y=b, x=SNP, ymin=b-(1.96*se), ymax=b+(1.96*se))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +
  geom_errorbarh(xmin="", xmax="", height=.1) +
  coord_flip() +
  ylab("Effect estimate of GlycA SNPs on depressive symptoms after exlusing individual SNP") +
  xlab("SNP")+
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.7) +
  ggtitle("Effect estimate with each SNP excluded in turn using GlycA SNPs on depressive symptoms") +
  theme_classic() 


jointdataset2 <- merge(res_loo2, GlycA, by = 'SNP')
jointdataset2 <- jointdataset2[ , c ('SNP', 'exposure.x' , 'outcome', 'pos.exposure', 'chr.exposure', 
                                   'effect_allele.exposure', 'other_allele.exposure' , 'eaf.exposure' , 'b', 'se', 'p')]

write.csv(jointdataset2, ".combined_data_GlycA_dep_symp.csv", row.names=FALSE)
#steiger
steiger<- steiger_filtering(dat)
steiger <- steiger %>% select(-(remove:pos))
steiger <- steiger %>% select(-(originalname.outcome:id.outcome))
steiger <- steiger %>% select(-(chr.exposure))
steiger <- steiger %>% select(-(id.exposure))
steiger <- steiger %>% select(-(exposure:mr_keep))
steiger <- steiger %>% select(-(MAF))
steiger <- steiger %>% select(-(units.outcome))
steiger <- steiger %>% select(-(units.exposure))
steiger <- steiger %>% select(-(effective_n.exposure))
steiger <- steiger %>% select(-(effective_n.outcome))
write.csv(steiger,"./steiger.csv")
steiger_dep<- steiger_filtering(dat2)
steiger_dep <- steiger_dep %>% select(-(remove:pos))
steiger_dep <- steiger_dep %>% select(-(originalname.outcome:id.outcome))
steiger_dep <- steiger_dep %>% select(-(chr.exposure))
steiger_dep <- steiger_dep %>% select(-(id.exposure))
steiger_dep <- steiger_dep %>% select(-(exposure:mr_keep))
steiger_dep <- steiger_dep %>% select(-(MAF))
steiger_dep <- steiger_dep %>% select(-(units.outcome))
steiger_dep <- steiger_dep %>% select(-(units.exposure))
steiger_dep <- steiger_dep %>% select(-(effective_n.exposure))
steiger_dep <- steiger_dep %>% select(-(effective_n.outcome))
write.csv(steiger_dep,"./steiger_GlycA_dep.csv")

