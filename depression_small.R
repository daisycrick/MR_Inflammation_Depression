###########################################
#Depression on GlycA#
###########################################
rm(list = ls())
setwd("C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/year two/SMFQ and GlycA/two sample/Depression_GlycA_Small") # Set the file path to your local directory

library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(MRInstruments)
library(purrr)

ao <- available_outcomes()
#get instruments from CB GWAS on UKBB
depression <-read_exposure_data("Depression_snps.txt") #47;

depression <- clump_data(
  depression,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

 #2;
GlycA <- extract_outcome_data(snps=depression$SNP, outcomes="met-d-GlycA",proxies = T)

#HARMONISE THE DATA
dat <- harmonise_data(depression, GlycA, action = 2)
# Check data
head(dat)
dim(dat)

mr_report(
  dat,
  output_path = "C:/Users/dc15053/OneDrive - University of Bristol/Documents/PhD/year two/SMFQ and GlycA/two sample/Depression_GlycA_Small",
  output_type = "html",
  author = "Analyst",
  study = "Two Sample MR",
  path = system.file("reports", package = "TwoSampleMR"),
  )


dat$beta.outcome<-(dat$beta.outcome*-1)

merged.data <- merge(dat, depression, by=c("SNP"))
merged.data <- subset(merged.data, select = c("SNP", "pos", "chr", "effect_allele.exposure.x", "other_allele.exposure.x", "eaf.outcome",  "outcome", "pval.exposure.x", "exposure.x", "beta.exposure.x", "se.exposure.x"))

write.xlsx(merged.data, ".combined_data_MDDGlycA_SA.xlsx" )

# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_depression_small_GlycA <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_depression_small_GlycA
Depression_small_GlycA<-cbind.data.frame(mr_depression_small_GlycA$outcome,mr_depression_small_GlycA$nsnp,mr_depression_small_GlycA$method,mr_depression_small_GlycA$b,mr_depression_small_GlycA$se,mr_depression_small_GlycA$pval)

#Export results
write.csv(Depression_small_GlycA,"./Depression_small_GlycA_results_final.csv")
# Estimate odds ratio and 95% confidence interval
(mr_depression_small_GlycA$b[1])
(mr_depression_small_GlycA$b[1]-1.96*mr_depression_small_GlycA$se[1])
(mr_depression_small_GlycA$b[1]+1.96*mr_depression_small_GlycA$se[1])

(mr_depression_small_GlycA$b[2])
(mr_depression_small_GlycA$b[2]-1.96*mr_depression_small_GlycA$se[2])
(mr_depression_small_GlycA$b[2]+1.96*mr_depression_small_GlycA$se[2])

(mr_depression_small_GlycA$b[3])
(mr_depression_small_GlycA$b[3]-1.96*mr_depression_small_GlycA$se[3])
(mr_depression_small_GlycA$b[3]+1.96*mr_depression_small_GlycA$se[3])

(mr_depression_small_GlycA$b[4])
(mr_depression_small_GlycA$b[4]-1.96*mr_depression_small_GlycA$se[4])
(mr_depression_small_GlycA$b[4]+1.96*mr_depression_small_GlycA$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat)

mr_pleiotropy_test(dat)

res_single_small <- mr_singlesnp(dat)
res_single_small
write.csv(res_single_small, "./singlesnp.csv")


# Get Fstat and R^2
dat$samplesize.exposure = 173005
dat$samplesize.outcome = 115078

dat$r.exposure <- get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_pn(dat$pval.outcome, dat$samplesize.outcome)
r2=directionality_test(dat)
exposure_r2=r2[1,"snp_r2.exposure"]

Samplesize=as.numeric(dat$samplesize.exposure[1])
nsnp=as.numeric(mr_depression_small_GlycA$nsnp[1])
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


##########################################
#VISUALIZE THE CAUSAL EFFECT OF GlycA ON depressive#
###########################################
# Generate a scatter plot comparing the different methods
png("./Depression_small_Glyca_scatter.png")
scatter_plot<-mr_scatter_plot(mr_depression_small_GlycA, dat)

dev.off()


pot <- ggplot(dat, aes())



res_single_small <- read.xlsx("singlesnp.xlsx", rowNames = F)




# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./Depression_small_GlycA_forest.png")
plot<-mr_forest_plot(res_single_small)
dev.off()
print(plot)
dev.off()



res_single_small %>%
  mutate(SNP = reorder(SNP, SNP2)) %>%
  ggplot(aes(y=b, x=SNP, ymin=b-(1.96*se), ymax=b+(1.96*se))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +
  geom_errorbarh(xmin="", xmax="", height=.1) +
  coord_flip() +
  ylab("Effect estimate of individual MDD SNP on GlycA") +
  xlab("SNP")+
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.7) +
  ggtitle("Effect of Individual MDD SNPs on GlycA for the Sensitivity Analyis") +
  theme_classic() 



# Generate a funnel plot to check asymmetry
png("./Depression_small_GlycA_funnel2.png")
mr_funnel_plot(res_single_small)
dev.off()




# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png("./Depression_small_GlycA_loo.png")
mr_leaveoneout_plot(res_loo)
dev.off()

res_loo %>%
  mutate(SNP = reorder(SNP, b)) %>%
  ggplot(aes(y=b, x=SNP, ymin=b-(1.96*se), ymax=b+(1.96*se))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +
  geom_errorbarh(xmin="", xmax="", height=.1) +
  coord_flip() +
  ylab("Effect estimate of MDD SNPs on GlycA after exlusing individual SNP") +
  xlab("SNP")+
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.7) +
  ggtitle("Effect estimate with each SNP excluded in turn using sensitivity MDD SNPs on GlycA") +
  theme_classic() 




#steiger
steiger<- steiger_filtering(dat)
dat_2 <-steiger[!(steiger$steiger_dir=="FALSE"),]
steiger <- steiger %>% select(-(remove:pos))
steiger <- steiger %>% select(-(originalname.outcome:data_source.outcome))
steiger <- steiger %>% select(-(outcome))
steiger <- steiger %>% select(-(id.exposure))
steiger <- steiger %>% select(-(exposure:mr_keep))
steiger <- steiger %>% select(-(eaf.exposure))
steiger <- steiger %>% select(-(units.outcome))
steiger <- steiger %>% select(-(units.exposure))
steiger <- steiger %>% select(-(effective_n.exposure))
steiger <- steiger %>% select(-(effective_n.outcome))
write.csv(steiger,"./steiger.csv")

#with SNPs removed that did not pass filtering

# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_depression_small_GlycA_2 <- mr(dat_2, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_depression_small_GlycA_2
Depression_small_GlycA_2<-cbind.data.frame(mr_depression_small_GlycA_2$outcome,mr_depression_small_GlycA_2$nsnp,mr_depression_small_GlycA_2$method,mr_depression_small_GlycA_2$b,mr_depression_small_GlycA_2$se,mr_depression_small_GlycA_2$pval)

#Export results
write.csv(Depression_small_GlycA_2,"./Depression_small_GlycA_results_final_2.csv")
# Estimate odds ratio and 95% confidence interval
(mr_depression_small_GlycA$b[1])
(mr_depression_small_GlycA$b[1]-1.96*mr_depression_small_GlycA$se[1])
(mr_depression_small_GlycA$b[1]+1.96*mr_depression_small_GlycA$se[1])

(mr_depression_small_GlycA$b[2])
(mr_depression_small_GlycA$b[2]-1.96*mr_depression_small_GlycA$se[2])
(mr_depression_small_GlycA$b[2]+1.96*mr_depression_small_GlycA$se[2])

(mr_depression_small_GlycA$b[3])
(mr_depression_small_GlycA$b[3]-1.96*mr_depression_small_GlycA$se[3])
(mr_depression_small_GlycA$b[3]+1.96*mr_depression_small_GlycA$se[3])

(mr_depression_small_GlycA$b[4])
(mr_depression_small_GlycA$b[4]-1.96*mr_depression_small_GlycA$se[4])
(mr_depression_small_GlycA$b[4]+1.96*mr_depression_small_GlycA$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat_2)

mr_pleiotropy_test(dat_2)

res_single_small <- mr_singlesnp(dat_2)
res_single_small
write.csv(res_single_small, "./singlesnp.csv")


# Get Fstat and R^2
dat_2$samplesize.exposure = 173005
dat_2$samplesize.outcome = 115078

dat_2$r.exposure <- get_r_from_pn(dat_2$pval.exposure, dat_2$samplesize.exposure)
dat_2$r.outcome <- get_r_from_pn(dat_2$pval.outcome, dat_2$samplesize.outcome)
r2=directionality_test(dat_2)
exposure_r2=r2[1,"snp_r2.exposure"]

Samplesize=as.numeric(dat$samplesize.exposure[1])
nsnp=as.numeric(mr_depression_small_GlycA$nsnp[1])
f_stat= ((Samplesize-(nsnp-1))/nsnp)*(exposure_r2/(1-exposure_r2))
f_stat=cbind(f_stat, exposure_r2)

