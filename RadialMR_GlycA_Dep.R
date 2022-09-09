install.packages(c("devtools", "knitr", "rmarkdown"))
library(devtools)
install_github(c("MRCIEU/TwoSampleMR","MRCIEU/MRInstruments"))
library(devtools)
install_github("WSpiller/MRPracticals",build_opts = c("--no-resave-data", "--no-manual"),build_vignettes = TRUE)
install_github("WSpiller/RadialMR")
library(RadialMR)
remotes::install_github("WSpiller/RadialMR")
radial_data <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
head(radial_data)
ivw.model <- ivw_radial(radial_data,0.05/nrow(radial_data),3,0.0001)
ivw.model$outliers
egger.model<-egger_radial(radial_data,0.05/nrow(radial_data),3)
egger.model$outliers

IVWplot1<-plot_radial(ivw.model,T,F,F) #radial plot with reference scale
IVWplot2<-plot_radial(ivw.model,F,T,F) #outliers
IVWplot3<-plotly_radial(ivw.model) #can see RSID for follow up

Eggerplot1<-plot_radial(egger.model,T,F,F)
Eggerplot2<-plot_radial(egger.model,F,T,F)
Eggerplot3<-plotly_radial(egger.model)

Comboplot1<-plot_radial(c(ivw.model,egger.model),T,F,F) #combined plot with IVW and egger
Comboplot2<-plot_radial(c(ivw.model,egger.model),F,T,F)

#use this to remove outliers identified in the IVW model.
out_rem<-radial_data[radial_data$SNP %in% ivw.model$outliers$SNP,]
radial_data2<-radial_data[-c(as.numeric(row.names(out_rem))),]
ivw.model2<-ivw_radial(radial_data2,0.05/nrow(radial_data2),3,0.0001)

out_rem<-radial_data[radial_data$SNP %in% egger.model$outliers$SNP,]
radial_data3<-radial_data[-c(as.numeric(row.names(out_rem))),]
egger.model2<-egger_radial(radial_data3,0.05/nrow(radial_data3),3,0.0001)


#It is worth emphasising, however, that such an approach requires justification 
#and can potentially lead to a loss of valuable information. Better practice would 
#be to look up outliers using a tool such as Phenoscanner, to see if there is a pattern 
#in the set of phenotypes with which they are associated. If we have instruments for a 
#suspected pleiotropic pathway, we can then fit a multivariable MR model.


radial_data_dep <- format_radial(dat2$beta.exposure, dat2$beta.outcome, dat2$se.exposure, dat2$se.outcome, dat2$SNP)
head(radial_data_dep)
ivw.model2 <- ivw_radial(radial_data_dep,0.05/nrow(radial_data_dep),3,0.0001)
ivw.model2$outliers
egger.model2<-egger_radial(radial_data_dep,0.05/nrow(radial_data_dep),3)
egger.model2$outliers

IVWplot1_2<-plot_radial(ivw.model2,T,F,F) #radial plot with reference scale
IVWplot2_2<-plot_radial(ivw.model2,F,T,F) #outliers
IVWplot3_2<-plotly_radial(ivw.model2) #can see RSID for follow up

Eggerplot1_2<-plot_radial(egger.model2,T,F,F)
Eggerplot2_2<-plot_radial(egger.model2,F,T,F)
IVWplot3_2<-plotly_radial(egger.model2)

Comboplot1_2<-plot_radial(c(ivw.model2,egger.model2),T,F,F) #combined plot with IVW and egger
Comboplot2_2<-plot_radial(c(ivw.model2,egger.model2),F,T,F)

#use this to remove outliers identified in the IVW model.
out_rem2<-radial_data2[radial_data_dep$SNP %in% ivw.model2$outliers$SNP,]

radial_data2<-radial_data2[-c(as.numeric(row.names(out_rem2))),]
ivw.model2<-ivw_radial2(radial_data2,0.05/nrow(radial_data2),3,0.0001)

#It is worth emphasising, however, that such an approach requires justification 
#and can potentially lead to a loss of valuable information. Better practice would 
#be to look up outliers using a tool such as Phenoscanner, to see if there is a pattern 
#in the set of phenotypes with which they are associated. If we have instruments for a 
#suspected pleiotropic pathway, we can then fit a multivariable MR model.