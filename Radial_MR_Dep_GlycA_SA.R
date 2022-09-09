radial_data <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
head(radial_data)
ivw.model <- ivw_radial(radial_data,0.05/nrow(radial_data),3,0.0001)
ivw.model$outliers
egger.model<-egger_radial(radial_data,0.05/nrow(radial_data),3)
egger.model$outliers
#ivw
0.07466473   -1.96*0.02823036     
0.07466473   +1.96*0.02823036

#mr egger
-0.02293311   -1.96*0.1085519     
-0.02293311   +1.96*0.1085519

IVWplot1<-plot_radial(ivw.model,T,F,F) #radial plot with reference scale
IVWplot2<-plot_radial(ivw.model,F,T,F) #outliers
IVWplot2
IVWplot3<-plotly_radial(ivw.model) #can see RSID for follow up
IVWplot3

Eggerplot1<-plot_radial(egger.model,T,F,F)
Eggerplot2<-plot_radial(egger.model,F,T,F)
Eggerplot3<-plotly_radial(egger.model)

Comboplot1<-plot_radial(c(ivw.model,egger.model),T,F,F) #combined plot with IVW and egger
Comboplot2<-plotly_radial(ivw.model)
Comboplot3<-plotly_radial(egger.model)

#use this to remove outliers identified in the IVW model.
out_rem<-radial_data[radial_data$SNP %in% egger.model$outliers$SNP,]
radial_data<-radial_data[-c(as.numeric(row.names(out_rem))),]
egger.model2<-egger_radial(radial_data,0.05/nrow(radial_data),3)