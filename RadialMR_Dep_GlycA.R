install.packages(c("devtools", "knitr", "rmarkdown"))
library(devtools)
install_github(c("MRCIEU/TwoSampleMR","MRCIEU/MRInstruments"))
library(devtools)
install_github("WSpiller/MRPracticals",build_opts = c("--no-resave-data", "--no-manual"),build_vignettes = TRUE)
install_github("WSpiller/RadialMR")
library(RadialMR)
remotes::install_github("WSpiller/RadialMR")


radial_data <- format_radial(dat2$beta.exposure, dat2$beta.outcome, dat2$se.exposure, dat2$se.outcome, dat2$SNP)
head(radial_data)
ivw.model <- ivw_radial(radial_data,0.05/nrow(radial_data),3,0.0001)
ivw.model$outliers
egger.model<-egger_radial(radial_data,0.05/nrow(radial_data),3)
egger.model$outliers

IVWplot1<-plot_radial(ivw.model,T,F,F) #radial plot with reference scale
IVWplot2<-plot_radial(ivw.model,F,T,F) #outliers
IVWplot2
IVWplot3<-plotly_radial(ivw.model) #can see RSID for follow up
IVWplot3

Eggerplot1<-plot_radial(egger.model,T,F,F)
Eggerplot2<-plot_radial(egger.model,F,T,F)
IVWplot3<-plotly_radial(egger.model)

Comboplot1<-plot_radial(c(ivw.model,egger.model),T,F,F) #combined plot with IVW and egger
Comboplot2<-plot_radial(c(ivw.model,egger.model),F,T,F)

#use this to remove outliers identified in the IVW model.
out_rem<-radial_data[radial_data$SNP %in% ivw.model$outliers$SNP,]
radial_data<-radial_data[-c(as.numeric(row.names(out_rem))),]
ivw.model2<-ivw_radial(radial_data,0.05/nrow(radial_data),3,0.0001)

out_rem<-radial_data[radial_data$SNP %in% egger.model$outliers$SNP,]
radial_data<-radial_data[-c(as.numeric(row.names(out_rem))),]
egger.model2<-egger_radial(radial_data,0.05/nrow(radial_data),3,0.0001)

