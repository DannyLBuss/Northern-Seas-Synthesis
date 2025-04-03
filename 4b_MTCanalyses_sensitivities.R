###MTC analysis - recovery method variation
##24th Jan 2025 - NSS analysis
#Author: Dr Danny L Buss
#1. load libraries ####
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

source("~/Documents/4-Oceans/3.NSS/2.MANUSCRIPT_Jan2025/Code_final/MANUSCRIPT_functions.R")

#2. setwd & load data ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5a_NSS_SpeciesData_Feb2025_withCEE.xlsx")
#df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5b_NSS_SpeciesData_Feb2025_noCEE.xlsx")

DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69")
df_species$time_bins<-factor(df_species$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_species$time_bins2<-factor(df_species$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))
df_species$Temp.mid[df_species$GBIF_species=="Esox lucius"]<- ((21 - 10) /2) + 10

#3. Remove rows without temperature trait data ####
NSS_sp<-df_species[!is.na(df_species$GBIF_species),]
NSS_sp<-NSS_sp[!NSS_sp$GBIF_species=="NA",]
NSS_sp<-NSS_sp[!NSS_sp$Trait_Temp_Max=="NA",]
NSS_sp<-NSS_sp[!NSS_sp$Trait_Temp_Min=="NA",]
#Only species without temperature data: Chondrostoma nasus

#4. Create time range variable for subsampling chronologies ####
NSS_sp$Time.Range<-as.numeric(NSS_sp$End_date_CE) - as.numeric(NSS_sp$Start_date_CE)
NSS_sp<-NSS_sp[!NSS_sp$Time.Range=="NA",]

#5. Check temperatures are numeric and remove species with zero counts or mixed life history strategies ####
NSS_sp$Trait_LifeHistory<-as.character(NSS_sp$Trait_LifeHistory)
NSS_sp<-NSS_sp[!NSS_sp$Trait_LifeHistory=="Unclassified; mixed strategies",]
NSS_sp<-NSS_sp[!NSS_sp$Trait_LifeHistory=="Amphidromous",]
NSS_sp<-NSS_sp[!is.na(NSS_sp$Trait_LifeHistory),]
NSS_sp<-NSS_sp[!NSS_sp$Trait_LifeHistory=="NA",]
NSS_sp<-NSS_sp %>%
  dplyr::mutate_at(c('Trait_Temp_Max', 'Trait_Temp_Min', 'Temp.mid', 'Time.Range'), as.numeric) 
NSS_sp$Temp.range<-NSS_sp$Trait_Temp_Max - NSS_sp$Trait_Temp_Min
NSS_sp<-NSS_sp[NSS_sp$Temp.range < 21,]
#NSS_sp<-NSS_sp[!NSS_sp$GBIF_species=="Esox lucius",]

#6. Species and NISP counts in remaining dataset ####
NSS_sp$count<-"1"
NSS_sp<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                 !NSS_sp$DB_Assemblage_ID == "NA", 
                                 !NSS_sp$NISP < 1.5)

NSS_sp_hand<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                 !NSS_sp$DB_Assemblage_ID == "NA", 
                                 !NSS_sp$NISP < 1.5,
                                 NSS_sp$Recovery=="Hand-collected only")

NSS_sp_sieved<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                 !NSS_sp$DB_Assemblage_ID == "NA", 
                                 !NSS_sp$NISP < 1.5,
                                 !NSS_sp$Recovery=="Hand-collected only",
                                 !NSS_sp$Recovery=="Unknown")

NSS_sp_nomarine<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                      !NSS_sp$DB_Assemblage_ID == "NA", 
                                      !NSS_sp$NISP < 1.5,
                                      NSS_sp$Trait_LifeHistory=="Potamodromous",
                                      NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")

NSS_sp_marine<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                          !NSS_sp$DB_Assemblage_ID == "NA", 
                                          !NSS_sp$NISP < 1.5,
                                          NSS_sp$Trait_LifeHistory=="Oceanodromous",
                                          NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")

NSS_sp_marinetrade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                        !NSS_sp$DB_Assemblage_ID == "NA", 
                                        !NSS_sp$NISP < 1.5,
                                        !NSS_sp$Trait_LifeHistory=="Potamodromous",
                                        NSS_sp$Known_Major_Trade_J_Barrett=="TRUE",
                                        NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")

NSS_sp_freshwaternotrade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                             !NSS_sp$DB_Assemblage_ID == "NA", 
                                             !NSS_sp$NISP < 1.5,
                                             NSS_sp$Trait_LifeHistory=="Potamodromous",
                                             NSS_sp$Known_Major_Trade_J_Barrett=="FALSE",
                                             NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")

NSS_sp_marinenotrade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                             !NSS_sp$DB_Assemblage_ID == "NA", 
                                             !NSS_sp$NISP < 1.5,
                                             !NSS_sp$Trait_LifeHistory=="Potamodromous",
                                             NSS_sp$Known_Major_Trade_J_Barrett=="FALSE",
                                             NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")

NSS_sp_freshwateraqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                             !NSS_sp$DB_Assemblage_ID == "NA", 
                                             !NSS_sp$NISP < 1.5,
                                             !NSS_sp$Trait_LifeHistory=="Oceanodromous",
                                             NSS_sp$Known_Major_Trade_J_Barrett=="FALSE",
                                             NSS_sp$Known_Aquaculture_R_Hoffmann=="TRUE")

NSS_sp_marinenoaqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                                !NSS_sp$DB_Assemblage_ID == "NA", 
                                                !NSS_sp$NISP < 1.5,
                                                !NSS_sp$Trait_LifeHistory=="Potamodromous",
                                                NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")


NSS_sp_nofreshwater<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                        !NSS_sp$DB_Assemblage_ID == "NA", 
                                        !NSS_sp$NISP < 1.5,
                                        !NSS_sp$Trait_LifeHistory=="Potamodromous")

NSS_sp_Denmark<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                          !NSS_sp$DB_Assemblage_ID == "NA", 
                                          !NSS_sp$NISP < 1.5,
                                          NSS_sp$Country=="Denmark")

NSS_sp_NorwaySweden<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                              !NSS_sp$DB_Assemblage_ID == "NA", 
                                              !NSS_sp$NISP < 1.5,
                                              NSS_sp$Region=="Scandinavia",
                                              !NSS_sp$Country=="Denmark")

NSS_sp_hand %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_sieved %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_nofreshwater %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_nomarine %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_NorwaySweden %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_Denmark %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_sieved %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_hand %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_nofreshwater %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_nomarine %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_Denmark %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_NorwaySweden %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_sieved %>%
  dplyr::group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_hand %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_sieved %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_nomarine %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_nofreshwater %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_Denmark %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_NorwaySweden %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_noaqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                        !NSS_sp$DB_Assemblage_ID == "NA", 
                                        !NSS_sp$NISP < 1.5,
                                        NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")

NSS_sp_notrade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                         !NSS_sp$DB_Assemblage_ID == "NA", 
                                         !NSS_sp$NISP < 1.5,
                                         NSS_sp$Known_Major_Trade_J_Barrett=="FALSE")

NSS_sp_notradeaqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                             !NSS_sp$DB_Assemblage_ID == "NA", 
                                             !NSS_sp$NISP < 1.5,
                                             NSS_sp$Known_Major_Trade_J_Barrett=="FALSE",
                                             NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")
NSS_sp_noaqua %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_notrade %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_notradeaqua %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp_notrade %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_nomarine %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_noaqua %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_notradeaqua %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp_notrade %>%
  dplyr::group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_noaqua %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_notrade %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_notradeaqua %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

####Recovery
#100 year intervals - Global and regional aMTC (Sieved sites only) ####
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_sieved
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "All species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "All species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland (n=612)",
                   "All species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe (n=22)",
                   "All species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia (n=162)",
                   "All species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe (n=22)")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)


png("../SupplementaryFiles/Figure_MTC_sieved.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Sieved only; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

#100 year intervals - Global and regional aMTC (Hand-collected sites only) ####
NSS_sp<-NSS_sp_hand
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"All species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "All species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "All species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland (n=612)",
                   "All species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe (n=22)",
                   "All species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia (n=162)",
                   "All species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe (n=22)")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)


png("../SupplementaryFiles/Figure_MTC_sieved.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Sieved only; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()


####Without traded and aquaculture species (removes LOTS of NISP)
#100 year intervals - Global and regional aMTC (no aquaculture) ####
###RUN1###
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_noaqua
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without aquaculture species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without aquaculture species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture specie; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture specie; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "All species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "No traded species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland",
                   "No traded species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe",
                   "No traded species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia",
                   "All species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)

png("../SupplementaryFiles/Figure_MTC_noaqua.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

writexl::write_xlsx(final,"../SupplementaryFiles/MTC_sensitivities/MTC_noaqua.xlsx")

#100 year intervals - Global and regional aMTC (no trade) ####
###RUN1###
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_notrade
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without traded species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without traded species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "Without traded species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "Without traded species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland",
                   "Without traded species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe",
                   "Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia",
                   "Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)

png("../SupplementaryFiles/Figure_MTC_notrade.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

png("../SupplementaryFiles/Figure_MTC_notradeGlobal.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[7]) +
  scale_fill_manual(values=DB[7]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

writexl::write_xlsx(final,"../SupplementaryFiles/MTC_sensitivities/MTC_notrade.xlsx")

#100 year intervals - Global and regional aMTC (no aquaculture or trade) ####
###RUN1###
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_notradeaqua
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture or traded species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without aquaculture species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without aquaculture or traded species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture specie; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture specie; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "Without aquaculture or traded species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "Without aquaculture species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland",
                   "Without aquaculture species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe",
                   "Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia",
                   "Without aquaculture species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)

png("../SupplementaryFiles/Figure_MTC_noaquaculturenotrade.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture or traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

png("../SupplementaryFiles/Figure_MTC_noaquaculturenotradeGlobal.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture or traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[7]) +
  scale_fill_manual(values=DB[7]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

writexl::write_xlsx(final,"../SupplementaryFiles/MTC_sensitivities/MTC_noaquanotrade.xlsx")


####Oceanodromous or potamodromous patterns
#100 year intervals - Global and regional aMTC (Oceanodromous only) ####
###RUN1###
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_nofreshwater
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without freshwater species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without freshwater species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without freshwater species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without freshwater species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without freshwater species; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without freshwater species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No freshwater species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without freshwater species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without freshwater species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No freshwater species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without freshwater species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without freshwater species; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "Without freshwater species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "Without freshwater species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland",
                   "Without freshwater species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe",
                   "Without freshwater species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia",
                   "Without freshwater species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)

png("../SupplementaryFiles/Figure_MTC_nofreshwater.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without freshwater species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

png("../SupplementaryFiles/Figure_MTC_nofreshwaterGlobal.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without freshwater species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[7]) +
  scale_fill_manual(values=DB[7]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

writexl::write_xlsx(final,"../SupplementaryFiles/MTC_sensitivities/MTC_nofreshwater.xlsx")

#100 year intervals - Global and regional aMTC (Potamodromous only) ####
###RUN1###
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_nomarine
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Freshwater species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without marine species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without marine species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Freshwater species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Freshwater species; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Freshwater species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Freshwater species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Freshwater species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "Freshwater species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "Freshwater species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland",
                   "Without traded species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe",
                   "Freshwater species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia",
                   "Freshwater species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)

png("../SupplementaryFiles/Figure_MTC_freshwaternoaqua.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!labs == "Global",
                               TimePeriod < 1800), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Freshwater species (No Aquaculture); 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

png("../SupplementaryFiles/Figure_MTC_FreshwaterNoAquaGlobal.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Freshwater species (No Aquaculture); 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[7]) +
  scale_fill_manual(values=DB[7]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

writexl::write_xlsx(final,"../SupplementaryFiles/MTC_sensitivities/MTC_freshwater.xlsx")


####Scandinavia split into upper and lower

####Scandinavia split into upper and lower
#100 year intervals - Global and regional aMTC (Traded species only) ####
###RUN1###
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_marinetrade
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Traded Marine sp.; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Global"
a$Type<-"Traded Marine sp."
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-bind_rows(final, a)


#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Traded Marine sp.; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Britain & Ireland"
a$Type<-"Traded Marine sp."
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture specie; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Traded Marine sp; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Western Europe"
a$Type<-"Traded Marine sp"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Traded Marine sp; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Scandinavia"
a$Type<-"Traded Marine sp"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Central & Eastern Europe"
a$Type<-"Marine sp. (No Aquaculture)"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Traded Marine sp; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC, color=labs, fill=labs)) + geom_line(aes(y=MTC, colour = labs,),  linewidth = 1.2, alpha=0.9) + 
#  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, colour = labs, fill = labs), size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~Type, scale="free_y")

ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  theme_bw() + 
  ggtitle("100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~Type, scale="free_y")


NSS_sp<-NSS_sp_freshwatertrade
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Traded Marine sp.; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Global"
a$Type<-"Traded Freshwater sp."
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-bind_rows(final, a)


#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Traded Marine sp.; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Britain & Ireland"
a$Type<-"Traded Marine sp."
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture specie; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Traded Marine sp; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Western Europe"
a$Type<-"Traded Marine sp"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Traded Marine sp; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Scandinavia"
a$Type<-"Traded Marine sp"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Central & Eastern Europe"
a$Type<-"Marine sp. (No Aquaculture)"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Traded Marine sp; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC, color=labs, fill=labs)) + geom_line(aes(y=MTC, colour = labs,),  linewidth = 1.2, alpha=0.9) + 
  #  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, colour = labs, fill = labs), size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~Type, scale="free_y")

ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  theme_bw() + 
  ggtitle("100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~Type, scale="free_y")


png("../SupplementaryFiles/Figure_MTC_noaquaculture_all.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  theme_bw() + 
  ggtitle("100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~Type, scale="free_y")
dev.off()

png("../SupplementaryFiles/Figure_MTC_noaquaculture_allGlobal.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(labs == "Global"), aes(TimePeriod, MTC)) + 
  geom_line(aes(y=MTC), colour= NA, linewidth = 1.2, alpha=0.9) +
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), colour = NA, size = 1.2, alpha=0.3) +
  theme_bw() + 
  ggtitle("100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[3]) +
  scale_fill_manual(values=DB[3]) +
  theme(legend.position = "none") +
  facet_wrap(~Type, scale="free_y", ncol=1)
dev.off()

writexl::write_xlsx(final,"../SupplementaryFiles/MTC_sensitivities/MTC_noaquanotrade.xlsx")

#100 year intervals - Global and regional aMTC (traded species only)
NSS_sp<-NSS_sp_freshwatertrade
this_label<-"NSS_sp_freshwatertrade"
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-paste(this_label, "100-yr time bins (Global Data (no CEE), mean + SD)")
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-
a$Type
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-bind_rows(final, a)
final$Type<-c(rep("Freshwater sp. (No Aquaculture)",77),rep("Marine sp. (No Aquaculture)",19))
final$labs[78:96]<-"Global"

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Britain & Ireland"
a$Type<-"Marine sp. (No Aquaculture)"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture specie; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Western Europe"
a$Type<-"Marine sp. (No Aquaculture)"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Scandinavia"
a$Type<-"Marine sp. (No Aquaculture)"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs~Type, scale="free_y")

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
a$labs<-"Central & Eastern Europe"
a$Type<-"Marine sp. (No Aquaculture)"
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture specie; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC, color=labs, fill=labs)) + geom_line(aes(y=MTC, colour = labs,),  linewidth = 1.2, alpha=0.9) + 
  #  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, colour = labs, fill = labs), size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~Type, scale="free_y")

ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  theme_bw() + 
  ggtitle("100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~Type, scale="free_y")

#100 year intervals - Global and regional aMTC (Denmark only) ####
###RUN1###
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_notrade
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without traded species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without traded species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "Without traded species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "Without traded species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland",
                   "Without traded species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe",
                   "Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia",
                   "Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)

png("../SupplementaryFiles/Figure_MTC_notrade.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

png("../SupplementaryFiles/Figure_MTC_notradeGlobal.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[7]) +
  scale_fill_manual(values=DB[7]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

writexl::write_xlsx(final,"../SupplementaryFiles/MTC_sensitivities/MTC_notrade.xlsx")

#100 year intervals - Global and regional aMTC (Scandinavia only) ####
###RUN1###
#1. 200 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
NSS_sp<-NSS_sp_notrade
repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
repeats = as.character(gsub(" ","",repeats))
test = repeats

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
  distinct()
NSS_sp_tmp<-tmp1
#NSS_sp<-NSS_sp[,c(1:48)]
for(i in 1:length(repeats)){
  NSS_sp_tmp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_tmp)){
    NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                         100, 1, floor(runif(1, min=1, max=50000)))
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
rm(tmp1, NSS_sp_tmp)

tmp<-NSS_sp_tmp_all %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Global Data (no CEE), mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins (All Data, mean + SD); 200 simulations")

final<-a
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without traded species; 100-yr time bins (All Data, mean + SD); 200 simulations")

ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + 
  theme_bw() + ggtitle("Without traded species; 100-yr time bins (All Data, mean + SD); 200 simulations") + xlim(50,1850)

#1. Britain & Ireland aMTC
#2. Subset columns required to run aMTC function and melt dataframe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Britian & Ireland") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Britain & Ireland, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins (Britain & Ireland, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

#western Europe
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Western Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Western Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins (Western Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)

###Scandinavia
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Scandinavia") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("No traded species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

###CEE
tmp<-NSS_sp_tmp_all %>%
  dplyr::filter(Region=="Central & Eastern Europe") %>%
  dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
#Round here at the same time period bin size that you chose for generate_random_time_Data above
#tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)
tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)

#3. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","MTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD)"
names(a)<-c("TimePeriod","MTC","SDMTC","Run")
ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "darkred", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw()
ggplot(a, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + geom_ribbon(aes(y=MTC, ymin = MTC-0.0001, ymax = MTC + 0.0001), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD); 200 simulations")

final<-bind_rows(final, a)
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Run)
final$sizebins<-"100"

final$labs<-recode(final$Run,
                   "Without traded species; 100-yr time bins (Global Data (no CEE), mean + SD)" = "Global",
                   "Without traded species; 100-yr time bins (Britain & Ireland, mean + SD)" = "Britain & Ireland",
                   "Without traded species; 100-yr time bins (Western Europe, mean + SD)" = "Western Europe",
                   "Without aquaculture species; 100-yr time bins (Scandinavia, mean + SD)" = "Scandinavia",
                   "Without traded species; 100-yr time bins (Cental & Eastern Europe, mean + SD)" = "Central & Eastern Europe")

ggplot(final %>% dplyr::filter(!Run == "All species; 100-yr time bins (Global Data (no CEE), mean + SD)"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB) +
  scale_fill_manual(values=DB) +
  facet_wrap(~labs)

png("../SupplementaryFiles/Figure_MTC_notrade.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(!labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

png("../SupplementaryFiles/Figure_MTC_notradeGlobal.png", width = 7.5, height = 5, units = "in", res = 600)
ggplot(final %>% dplyr::filter(labs == "Global"), aes(TimePeriod, MTC)) + geom_line(aes(y=MTC, colour = labs), linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = labs), size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("Without traded species; 100-yr time bins; 200 simulations") + 
  scale_color_manual(values=DB[7]) +
  scale_fill_manual(values=DB[7]) +
  theme(legend.position = "none") +
  facet_wrap(~labs, scales="free_y")
dev.off()

writexl::write_xlsx(final,"../SupplementaryFiles/MTC_sensitivities/MTC_notrade.xlsx")


####Plots ####
ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 25-yr time bins; 200 simulations") + 
  facet_grid(sizebins~Run, scales="free")

ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 25-yr time bins; 200 simulations") + 
  facet_wrap(~Run, scales="free") + xlim(50,1850)

P1<-ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "#00CCCC", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("All species; 25-yr time bins; 200 simulations") + 
  facet_wrap(~Run, scales="free") + xlim(50,1850)

###Supplementary Plots
P1<-ggplot(final, aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "grey75", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "grey89", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("NSS-DB (no CEE), 200 simulations") + 
  facet_wrap(~Run, scales="free") + xlim(50,1850)

P1 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)


###Figure 2 MTC panels
final$plot_names<-recode(final$Run,
                         `All species; 100-yr time bins (Global Data (no CEE), mean + SD)` = "NSS-DB (no CEE) 100-yr",
                         `All species; 100-yr time bins (Britain & Ireland, mean + SD)` = "Britain & Ireland 100-yr",
                         `All species; 100-yr time bins (Western Europe, mean + SD)` = "Western Europe 100-yr",
                         `All species; 100-yr time bins (Scandinavia, mean + SD)` = "Scandinavia 100-yr",
                         `All species; 25-yr time bins (Global Data (no CEE), mean + SD)` = "NSS-DB (no CEE) 25-yr",
                         `All species; 25-yr time bins (Britain & Ireland, mean + SD)` = "Britain & Ireland 25-yr",
                         `All species; 25-yr time bins (Western Europe, mean + SD)` = "Western Europe 25-yr",
                         `All species; 25-yr time bins (Scandinavia, mean + SD)` = "Scandinavia 25-yr",
)
final$level<-recode(final$Run,
                    `All species; 100-yr time bins (Global Data (no CEE), mean + SD)` = "Global",
                    `All species; 100-yr time bins (Britain & Ireland, mean + SD)` = "Regional",
                    `All species; 100-yr time bins (Western Europe, mean + SD)` = "Regional",
                    `All species; 100-yr time bins (Scandinavia, mean + SD)` = "Regional",
                    `All species; 25-yr time bins (Global Data (no CEE), mean + SD)` = "Global",
                    `All species; 25-yr time bins (Britain & Ireland, mean + SD)` = "Regional",
                    `All species; 25-yr time bins (Western Europe, mean + SD)` = "Regional",
                    `All species; 25-yr time bins (Scandinavia, mean + SD)` = "Regional")
write_xlsx(final, "../Data/MTC_output.xlsx")

P1<-ggplot(final[final$sizebins=="25" &
                   final$level=="Regional",], aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "grey45", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "grey79", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("NSS-DB, 200 simulations") + 
  facet_wrap(~plot_names, scales="free", ncol=1) + xlim(50,1850)

P1 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)

P2<-ggplot(final[final$sizebins=="100" &
                   final$level=="Regional",], aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "grey45", linewidth = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "grey79", size = 1.2, alpha=0.3) + theme_bw() + 
  ggtitle("NSS-DB, 200 simulations") + 
  facet_wrap(~plot_names, scales="free", ncol=1) + xlim(50,1850)

P2 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)

P3<-ggplot(final[final$sizebins=="100" &
                   final$level=="Global",], aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "grey45", linewidth = 1.2, alpha=0.9) +                          # Main line
  geom_line(aes(y = MTC - SDMTC), linetype = "dashed", color = "grey55", size=1.5) +                          # Main line
  geom_line(aes(y = MTC + SDMTC), linetype = "dashed", color = "grey55", size=1.5) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "grey79", size = 1.2, alpha=0.3) + 
  theme_bw() +
  ggtitle("NSS-DB, 200 simulations") + 
  facet_wrap(~plot_names, scales="free", ncol=1) + xlim(50,1850)


P3 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)

P4<-ggplot(final[final$sizebins=="25" &
                   final$level=="Global",], aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "grey45", linewidth = 1.2, alpha=0.9) +                          # Main line
  geom_line(aes(y = MTC - SDMTC), linetype = "dashed", color = "grey79", size=0.5) +                          # Main line
  geom_line(aes(y = MTC + SDMTC), linetype = "dashed", color = "grey79", size=0.5) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "grey79", size = 1.2, alpha=0.3) + 
  theme_publish() +
  ggtitle("MTC, 200 simulations") + 
  labs(tag="A)") +
  facet_wrap(~plot_names, scales="free", ncol=1) + xlim(50,1850)


P4 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)


P4_A<-P4 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)


