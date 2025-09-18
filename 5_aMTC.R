###File 4. calculat aMTC and run sensitivity analyses ####
#5th Sept 2025 - NSS analysis
#Author: Dr Danny L Buss
#1. load libraries ####
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(bayesplot)
  library(forcats)
  library(stringr)
  library(writexl)
  library(scico)
  library(reprex)
  library(nlme)
  library(lme4)
  library(visreg)
  library(MASS)
  library(stargazer)
  library(plyr)
})
#2. load manuscript functions ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../ThemePublish.R")
source("MANUSCRIPT_functions_b.R")

#3. Import data (if not carrying straihgt on from 4_Aoristic_simulations.R)
#NSS_sp_tmp_all <- readxl::read_xlsx(
#  "df_species_with_100yrtime_increments_Sept2025.xlsx"
#)

NSS_sp_tmp_all %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Total_NISP       = sum(NISP, na.rm = TRUE),
    Assemblage_count = dplyr::n_distinct(DB_Assemblage_ID),
    .groups = "drop"
  )

NSS_sp_tmp_all$Temp.mid<-NSS_sp_tmp_all$Trait_Temp_Fishbase_Mid
repeat_size = 1000
repeats <- paste0("rep", seq_len(repeat_size))
test = repeats
keep_cols <- c("GBIF_species", "NISP", "Temp.mid", repeats)
bin_size = 100

###Add data filters here ####
#1) Main Run - all assemblages must have at least 5 NISP
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >= 5) %>%                           
  ungroup()

tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]

sum(tmp$NISP)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 1000 (takes a few minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"NSS-DB; 100-yr time bins (All Data)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 500 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("NSS-DB; 100-yr time bins (mean + SD); 500 simulations") + xlim(50,1850)

##Sieved assemblages only ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >= 5,
                DB_sieved == TRUE) %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"NSS DB (Sieved only)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

#Marine only ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >= 5,
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Marine sp. only"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

## Freshwater, no aquculture only ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >= 5,
                !Known_Aquaculture_R_Hoffmann == TRUE,
                Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Freshwater only (no aquaculture)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + 
  ylim(c(6, 20))
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

## Aquaculture only ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >=5,
                Known_Aquaculture_R_Hoffmann == TRUE) %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Aquaculture only"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + 
  ylim(c(6, 20))
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

## Marine (Untraded) ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >=5,
                !Known_Major_Trade_J_Barrett == TRUE,
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Marine (not traded only)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + 
  ylim(c(6, 20))
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

## Freshwater, no aquaculture only ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >=5,
                Known_Major_Trade_J_Barrett == TRUE,
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Marine (Known Traded sp.)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + 
  ylim(c(6, 20))
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

## Marine (Untraded) ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >=5,
                !Known_Major_Trade_J_Barrett == TRUE,
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Marine (not traded only)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + 
  ylim(c(6, 20))
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

## Rural, marine, no aquaculture ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >=5,
                Rural_urban_or_neither == "rural",
                !Known_Aquaculture_R_Hoffmann == TRUE,
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Rural, (marine sp. only; no aquaculture)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + 
  ylim(c(6, 20))
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

## Urban, marine, no aquaculture ####
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total >=5,
                Rural_urban_or_neither == "urban",
                !Known_Aquaculture_R_Hoffmann == TRUE,
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
sum(NSS_sp_tmp_all_b$NISP)
tmp <- NSS_sp_tmp_all_b[ , intersect(keep_cols, names(NSS_sp_tmp_all_b))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Urban, (marine sp. only; no aquaculture)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + 
  ylim(c(6, 20))
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Marine only, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

#1. Britain & Ireland aMTC ####
#2. Subset columns required to run aMTC function and melt dataframe
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total > 2,
                !Known_Aquaculture_R_Hoffmann == TRUE,
                !Rural_urban_or_neither == "neither",
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
BI<-NSS_sp_tmp_all_b %>% dplyr::filter(Region == "Britain & Ireland")
sum(BI$NISP)
tmp <- BI[ , intersect(keep_cols, names(BI))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Marine species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Britain & Ireland, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

#Scandinavia
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total > 2,
                !Known_Aquaculture_R_Hoffmann == TRUE,
                !Rural_urban_or_neither == "neither",
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
SC<-NSS_sp_tmp_all_b %>% dplyr::filter(Region == "Scandinavia")
sum(SC$NISP)
tmp <- SC[ , intersect(keep_cols, names(SC))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Marine species; 100-yr time bins (Scandinavia), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Scandinavia, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 500 simulations") + 
  facet_wrap(~Run)

###Western Europe
NSS_sp_tmp_all_b <- NSS_sp_tmp_all %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(NISP_total > 2,
                !Known_Aquaculture_R_Hoffmann == TRUE,
                !Rural_urban_or_neither == "neither",
                !Trait_LifeHistory == "Potamodromous") %>%                           
  ungroup()
WE<-NSS_sp_tmp_all_b %>% dplyr::filter(Region == "Western Europe")
sum(WE$NISP)
tmp <- WE[ , intersect(keep_cols, names(WE))]
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
tmp_sp$value<-as.character(tmp_sp$value)
tmp_sp$variable<-as.character(tmp_sp$variable)
tmp_sp$Temp.mid<-as.character(tmp_sp$Temp.mid)
tmp_sp<-aggregate(NISP~GBIF_species + Temp.mid + variable + value, data=tmp_sp, FUN=sum)
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "Subsample", "time.mid.round", "NISP_count")
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
tmp_sp$time<-as.numeric(tmp_sp$time.mid.round)

#Round here at the same time period bin size that you chose for generate_random_time_Data above
lst<-unique(tmp_sp$Subsample)

#Create dataframe of aMTCs for each simulation using custom function (currently there are 500 (takes c.20 minutes))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC_old(spec_count_new[spec_count_new$time > 50 & spec_count_new$time < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Marine species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Scandinavia, mean + SD); 500 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("MTC; 100-yr time bins; 500 simulations") + 
  facet_wrap(~Run, scales="free_y")
final$sizebins<-"100"

#Save outputs for figures
write_xlsx(final, "MTC_outputs.xlsx")
