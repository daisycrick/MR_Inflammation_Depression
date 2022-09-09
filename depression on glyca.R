rm(list = ls())
setwd("C:/") # Set the file path to your local directory

library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
library(MRInstruments)
library(purrr)
library(ggplot2)
library(dplyr)



install.packages('writexl')
library(writexl)



ao <- available_outcomes()
#get instruments from CB GWAS on UKBB
depression <- extract_instruments("ieu-b-102") #50;

depression <- clump_data(
  depression,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)


GlycA <- extract_outcome_data(snps=depression$SNP, outcomes="met-d-GlycA",proxies = T)


#HARMONISE THE DATA
dat <- harmonise_data(depression, GlycA, action = 2)

#Create mr-keep and drop snps
mr_keep <-dat[, c('mr_keep','SNP')]
dat_2 <- merge(dat, mr_keep,by="SNP")
dat_3 <- subset(dat_2, mr_keep == "TRUE")

dat_3 <- dat_3 %>% rename (mr_keep = mr_keep.x)
dat_3 <- dat_3[-c(mr_keep.y)]

# Check data
head(dat_3)
dim(dat_3)
write.table(dat_3, 'dat_dep.txt')

# Let's use the MR-Base R package to estimate the effects using the IVW, MR-Egger, weighted median and weighted mode methods
# Have a look at the mr_method_list() function 
mr_depression_GlycA <- mr(dat_3, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
mr_depression_GlycA
Depression_GlycA<-cbind.data.frame(mr_depression_GlycA$outcome,mr_depression_GlycA$nsnp,mr_depression_GlycA$method,mr_depression_GlycA$b,mr_depression_GlycA$se,mr_depression_GlycA$pval)

#Export results
write.csv(Depression_GlycA,"./Depression_GlycA_results.csv")

# Estimate odds ratio and 95% confidence interval
(mr_depression_GlycA$b[1])
(mr_depression_GlycA$b[1]-1.96*mr_depression_GlycA$se[1])
(mr_depression_GlycA$b[1]+1.96*mr_depression_GlycA$se[1])

(mr_depression_GlycA$b[2])
(mr_depression_GlycA$b[2]-1.96*mr_depression_GlycA$se[2])
(mr_depression_GlycA$b[2]+1.96*mr_depression_GlycA$se[2])

(mr_depression_GlycA$b[3])
(mr_depression_GlycA$b[3]-1.96*mr_depression_GlycA$se[3])
(mr_depression_GlycA$b[3]+1.96*mr_depression_GlycA$se[3])

(mr_depression_GlycA$b[4])
(mr_depression_GlycA$b[4]-1.96*mr_depression_GlycA$se[4])
(mr_depression_GlycA$b[4]+1.96*mr_depression_GlycA$se[4])

# Run some sensitivity analyses
mr_heterogeneity(dat_3)

mr_pleiotropy_test(dat_3)

res_single <- mr_singlesnp(dat_3)
res_single


# Get Fstat and R^2
dat_3$samplesize.exposure = 500199 
dat_3$samplesize.outcome = 115078 

dat_3$r.exposure <- get_r_from_pn(dat_3$pval.exposure, dat_3$samplesize.exposure)
r2=directionality_test(dat_3)
exposure_r2=r2[1,"snp_r2.exposure"]

Samplesize=as.numeric(dat$samplesize.exposure[1])
nsnp=as.numeric(mr_depression_GlycA$nsnp[1])
f_stat= ((Samplesize-(nsnp-1))/nsnp)*(exposure_r2/(1-exposure_r2))
f_stat=cbind(f_stat, exposure_r2)

###individual F-stats
dat_3$MAF <- 1-(dat_3$eaf.exposure)
dat_3$r2 <- 2*dat_3$MAF*(1-dat_3$MAF)*(dat_3$beta.exposure^2)

dat_3$F.stat <- (dat_3$r2 * (dat_3$samplesize.exposure-2))/(1-dat_3$r2)





#MR_PRESSO
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("rondolab/MR-PRESSO")
'force=TRUE'
library(MRPRESSO)
run_mr_presso(dat, NbDistribution = 1000, SignifThreshold = 0.05)



##########################################
#VISUALIZE THE CAUSAL EFFECT OF GlycA ON depressive#
###########################################
# Generate a scatter plot comparing the different methods
png("./Depression_Glyca_scatter.png")
mr_scatter_plot(mr_depression_GlycA, dat)
dev.off()


# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./Depression_GlycA_forest.png")
depression_GlycA_plot<-mr_forest_plot(res_single)
dev.off()


res_single %>%
  mutate(SNP = reorder(SNP, b)) %>%
  ggplot(aes(y=b, x=SNP, ymin=b-(1.96*se), ymax=b+(1.96*se))) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +
  geom_errorbarh(xmin="", xmax="", height=.1) +
  coord_flip() +
  ylab("Effect estimate of individual MDD SNP on GlycA") +
  xlab("SNP")+
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.7) +
  ggtitle("Effect of Individual MDD SNPs on GlycA") +
  theme_classic() 





# Generate a funnel plot to check asymmetry
png("./Depression_GlycA_funnel2.png")
mr_funnel_plot(res_single)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png("./Depression_GlycA_loo.png")
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
  ggtitle("Effect estimate with each SNP excluded in turn using MDD SNPs on GlycA") +
  theme_classic() 



jointdataset <- merge(res_loo, depression, by = 'SNP')
jointdataset <- jointdataset[ , c ('SNP', 'exposure.x' , 'outcome', 'pos.exposure', 'chr.exposure', 
                                    'effect_allele.exposure', 'other_allele.exposure' , 'eaf.exposure' , 'b', 'se', 'p')]

write.csv(jointdataset, ".combined_data_MDD_GlycA.csv", row.names=FALSE)




#steiger
steiger<- steiger_filtering(dat_3)
steiger <- steiger %>% select(-(remove:id.outcome))
steiger <- steiger %>% select(-(originalname.outcome:pos.exposure))
steiger <- steiger %>% select(-(chr.exposure))
steiger <- steiger %>% select(-(id.exposure))
steiger <- steiger %>% select(-(exposure:mr_keep))
steiger <- steiger %>% select(-(MAF))
steiger <- steiger %>% select(-(units.outcome))
steiger <- steiger %>% select(-(units.exposure))
steiger <- steiger %>% select(-(effective_n.exposure))
steiger <- steiger %>% select(-(effective_n.outcome))
write.csv(steiger,"./steiger_Dep_GlycA.csv")

















