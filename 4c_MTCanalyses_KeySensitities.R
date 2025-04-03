##24th Feb 2025 - NSS analysis
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
#df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5a_NSS_SpeciesData_Feb2025_withCEE.xlsx")
df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5b_NSS_SpeciesData_Feb2025_noCEE.xlsx")

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

#6. Species and NISP counts in subsetted datasets ####
NSS_sp$count<-"1"
NSS_sp<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                 !NSS_sp$DB_Assemblage_ID == "NA", 
                                 !NSS_sp$NISP < 1.5)

NSS_sp_marinetrade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                             !NSS_sp$DB_Assemblage_ID == "NA", 
                                             !NSS_sp$NISP < 1.5,
                                             NSS_sp$Trait_LifeHistory=="Oceanodromous",
                                             NSS_sp$Known_Major_Trade_J_Barrett=="TRUE",
                                             NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")

NSS_sp_freshwatertrade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                             !NSS_sp$DB_Assemblage_ID == "NA", 
                                             !NSS_sp$NISP < 1.5,
                                             !NSS_sp$Trait_LifeHistory=="Oceanodromous",
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

NSS_sp_marineaqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                              !NSS_sp$DB_Assemblage_ID == "NA", 
                                              !NSS_sp$NISP < 1.5,
                                              NSS_sp$Trait_LifeHistory=="Oceanodromous",
                                              NSS_sp$Known_Major_Trade_J_Barrett=="FALSE",
                                              NSS_sp$Known_Aquaculture_R_Hoffmann=="TRUE")

NSS_sp_freshwaternoaqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                                !NSS_sp$DB_Assemblage_ID == "NA", 
                                                !NSS_sp$NISP < 1.5,
                                                !NSS_sp$Trait_LifeHistory=="Oceanodromous",
                                                !NSS_sp$Known_Aquaculture_R_Hoffmann=="TRUE")

NSS_sp_marinenoaqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                            !NSS_sp$DB_Assemblage_ID == "NA", 
                                            !NSS_sp$NISP < 1.5,
                                            NSS_sp$Trait_LifeHistory=="Oceanodromous",
                                            !NSS_sp$Known_Aquaculture_R_Hoffmann=="TRUE")


NSS_sp_marine<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                              !NSS_sp$DB_Assemblage_ID == "NA", 
                                              !NSS_sp$NISP < 1.5,
                                              NSS_sp$Trait_LifeHistory=="Oceanodromous")


NSS_sp_freshwater<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                              !NSS_sp$DB_Assemblage_ID == "NA", 
                                              !NSS_sp$NISP < 1.5,
                                              !NSS_sp$Trait_LifeHistory=="Oceanodromous")

NSS_sp_Denmark<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                         !NSS_sp$DB_Assemblage_ID == "NA", 
                                         !NSS_sp$NISP < 1.5,
                                         NSS_sp$Country=="Denmark")

NSS_sp_NorwaySweden<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                              !NSS_sp$DB_Assemblage_ID == "NA", 
                                              !NSS_sp$NISP < 1.5,
                                              NSS_sp$Region=="Scandinavia",
                                              !NSS_sp$Country=="Denmark")

NSS_sp_Norway<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                              !NSS_sp$DB_Assemblage_ID == "NA", 
                                              !NSS_sp$NISP < 1.5,
                                              NSS_sp$Country=="Norway")

NSS_sp_Sweden<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                              !NSS_sp$DB_Assemblage_ID == "NA", 
                                              !NSS_sp$NISP < 1.5,
                                              NSS_sp$Country=="Sweden")

NSS_sp_Netherlands<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                        !NSS_sp$DB_Assemblage_ID == "NA", 
                                        !NSS_sp$NISP < 1.5,
                                        NSS_sp$Country=="Netherlands")

NSS_sp_Belgium<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                             !NSS_sp$DB_Assemblage_ID == "NA", 
                                             !NSS_sp$NISP < 1.5,
                                             NSS_sp$Country=="Belgium")

NSS_sp_Scotland<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                         !NSS_sp$DB_Assemblage_ID == "NA", 
                                         !NSS_sp$NISP < 1.5,
                                         NSS_sp$Country=="Scotland")

NSS_sp_England<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                         !NSS_sp$DB_Assemblage_ID == "NA", 
                                         !NSS_sp$NISP < 1.5,
                                         NSS_sp$Country=="England")

NSS_sp_Germany<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                         !NSS_sp$DB_Assemblage_ID == "NA", 
                                         !NSS_sp$NISP < 1.5,
                                         NSS_sp$Country=="Germany")

NSS_sp_Ireland<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                         !NSS_sp$DB_Assemblage_ID == "NA", 
                                         !NSS_sp$NISP < 1.5,
                                         NSS_sp$Country=="Ireland")

NSS_sp_N_Ireland<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                         !NSS_sp$DB_Assemblage_ID == "NA", 
                                         !NSS_sp$NISP < 1.5,
                                         NSS_sp$Country=="Northern Ireland")

#Loop through datasets and calculate MTC ####
alldataframes<-list(NSS_sp_marine, NSS_sp_marineaqua, NSS_sp_marinenoaqua, NSS_sp_marinetrade, NSS_sp_marinenotrade,
                    NSS_sp_freshwater, NSS_sp_freshwateraqua, NSS_sp_freshwaternoaqua, NSS_sp_freshwatertrade, NSS_sp_freshwaternotrade)
names(alldataframes)<-c("NSS_sp_marine","NSS_sp_marineaqua","NSS_sp_marinenoaqua","NSS_sp_marinetrade","NSS_sp_marinenotrade",
                        "NSS_sp_freshwater","NSS_sp_freshwateraqua","NSS_sp_freshwaternoaqua","NSS_sp_freshwatertrade","NSS_sp_freshwaternotrade")

# Apply the function to each dataframe and combine results
summary_results <- lapply(names(alldataframes), function(name) summarise_df(alldataframes[[name]], name))

# Combine into a single dataframe
final_summary_table <- dplyr::bind_rows(summary_results)

# Print the summary table
print(final_summary_table)
writexl::write_xlsx(final_summary_table,"../SupplementaryFiles/MTC_sensitivities/MTC_summarystates_counts.xlsx")

alldataframes<-list(NSS_sp_marine, NSS_sp_marinenoaqua, NSS_sp_marinetrade, NSS_sp_marinenotrade,
                    NSS_sp_freshwater, NSS_sp_freshwateraqua, NSS_sp_freshwaternoaqua, NSS_sp_freshwatertrade, NSS_sp_freshwaternotrade)
names(alldataframes)<-c("NSS_sp_marine","NSS_sp_marinenoaqua","NSS_sp_marinetrade","NSS_sp_marinenotrade",
                        "NSS_sp_freshwater","NSS_sp_freshwateraqua","NSS_sp_freshwaternoaqua","NSS_sp_freshwatertrade","NSS_sp_freshwaternotrade")

###loop through dataframes to calculate MTC globally and per region
final <- data.frame(TimePeriod = character(), MTC = numeric(), SDMTC = numeric(), Run = character(), labs = character(), Type = character())
MTC_results <- lapply(names(alldataframes), function(name) MTC_df(alldataframes[[name]], name))
final_MTC_table <- dplyr::bind_rows(MTC_results)
print(final_MTC_table)
writexl::write_xlsx(final_MTC_table,"../SupplementaryFiles/MTC_sensitivities/MTC_marine_freshwater.xlsx")

alldataframes<-list(NSS_sp_Denmark, NSS_sp_NorwaySweden)
names(alldataframes)<-c("NSS_sp_Denmark","NSS_sp_NorwaySweden")
final <- data.frame(TimePeriod = character(), MTC = numeric(), SDMTC = numeric(), Run = character(), labs = character(), Type = character())
summary_results <- lapply(names(alldataframes), function(name) summarise_df(alldataframes[[name]], name))
final_summary_table_scandinavia <- dplyr::bind_rows(summary_results)
writexl::write_xlsx(final_summary_table_scandinavia,"../SupplementaryFiles/MTC_sensitivities/MTC_scandi_counts.xlsx")
MTC_results <- lapply(names(alldataframes), function(name) MTC_df_global(alldataframes[[name]], name))
final_MTC_table_scandi <- dplyr::bind_rows(MTC_results)
print(final_MTC_table_scandi)
writexl::write_xlsx(final_MTC_table_scandi,"../SupplementaryFiles/MTC_sensitivities/MTC_scandi.xlsx")

final_MTC_table$Year<-final_MTC_table$TimePeriod
final_MTC_table$TimePeriod<-as.numeric(final_MTC_table$TimePeriod)
final_MTC_table$water<-recode(final_MTC_table$Type,
                              "NSS_sp_marinetrade" = "Marine",
                              "NSS_sp_marinenotrade" = "Marine",
                              "NSS_sp_marinenoaqua" = "Marine",
                              "NSS_sp_marine" = "Marine",
                              "NSS_sp_freshwater" = "Freshwater",
                              "NSS_sp_freshwaternotrade" = "Freshwater",
                              "NSS_sp_freshwateraqua" = "Freshwater",
                              "NSS_sp_freshwaternoaqua" = "Freshwater",
                              "NSS_sp_freshwatertrade" = "Freshwater")
final_MTC_table$subset<-recode(final_MTC_table$Type,
                              "NSS_sp_marinetrade" = "Trade",
                              "NSS_sp_marinenotrade" = "No_Trade",
                              "NSS_sp_marinenoaqua" = "No_Aquaculture",
                              "NSS_sp_marine" = "All_Species",
                              "NSS_sp_freshwater" = "All_Species",
                              "NSS_sp_freshwaternotrade" = "No_Trade",
                              "NSS_sp_freshwateraqua" = "Aquaculture",
                              "NSS_sp_freshwaternoaqua" = "No_Aquaculture",
                              "NSS_sp_freshwatertrade" = "Trade")
ggplot(final_MTC_table, aes(TimePeriod, MTC)) + 
  geom_line(aes(y = MTC), colour = "#00CCCC", linewidth = 1.2, alpha = 0.9) + 
  geom_ribbon(aes(y = MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha = 0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_grid(water ~ subset + labs, scale = "free_y")

png(filename="../SupplementaryFiles/MTC_freshwater_sensitivity.png", 
    width=8.5, height=5.5, units="in",res=600)
ggplot(final_MTC_table %>% dplyr::filter(water=="Freshwater"), aes(TimePeriod, MTC, color = subset)) + 
  geom_line(aes(y = MTC, color = subset, fill = subset), linewidth = 1.2, alpha = 0.9) + 
  geom_ribbon(aes(y = MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = subset),color = NA, size = 1.2, alpha = 0.3) + 
  theme_bw() + ggtitle("Marine; 100-yr time bins; 200 simulations") + 
  facet_grid(labs~subset, scale = "free_y") + 
  scale_color_manual(values=DB[2:6]) + 
  scale_fill_manual(values=DB[2:6])
dev.off()
png(filename="../SupplementaryFiles/MTC_marine_sensitivity.png", 
    width=8.5, height=5.5, units="in",res=600)
ggplot(final_MTC_table %>% dplyr::filter(water=="Marine"), aes(TimePeriod, MTC, color = subset)) + 
  geom_line(aes(y = MTC, color = subset, fill = subset), linewidth = 1.2, alpha = 0.9) + 
  geom_ribbon(aes(y = MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = subset),color = NA, size = 1.2, alpha = 0.3) + 
  theme_bw() + ggtitle("Marine; 100-yr time bins; 200 simulations") + 
  facet_grid(labs~subset, scale = "free_y") + 
  scale_color_manual(values=DB[2:6]) + 
  scale_fill_manual(values=DB[2:6])
dev.off()

final_MTC_table_scandi$Year<-final_MTC_table_scandi$TimePeriod
final_MTC_table_scandi$TimePeriod<-as.numeric(final_MTC_table_scandi$TimePeriod)
ggplot(final_MTC_table_scandi, aes(TimePeriod, MTC)) + 
  geom_line(aes(y = MTC), colour = "#00CCCC", linewidth = 1.2, alpha = 0.9) + 
  geom_ribbon(aes(y = MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha = 0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(labs ~ Type, scale = "free_y")

###All countries ####
alldataframes<-list(NSS_sp_Denmark, NSS_sp_Norway, NSS_sp_Belgium, NSS_sp_England, NSS_sp_Germany, 
                    NSS_sp_Ireland, NSS_sp_Netherlands, NSS_sp_Scotland, NSS_sp_Sweden)
names(alldataframes)<-c("NSS_sp_Denmark","NSS_sp_Norway", "NSS_sp_Belgium", "NSS_sp_England", "NSS_sp_Germany",
                        "NSS_sp_Ireland", "NSS_sp_Netherlands", "NSS_sp_Scotland",
                        "NSS_sp_Sweden")
final <- data.frame(TimePeriod = character(), MTC = numeric(), SDMTC = numeric(), Run = character(), labs = character(), Type = character())
summary_results <- lapply(names(alldataframes), function(name) summarise_df(alldataframes[[name]], name))
final_summary_table_countries <- dplyr::bind_rows(summary_results)
writexl::write_xlsx(final_summary_table_scandinavia,"../SupplementaryFiles/MTC_sensitivities/MTC_country_counts.xlsx")

MTC_results <- lapply(names(alldataframes), function(name) MTC_df_global(alldataframes[[name]], name))
final_MTC_table_countries <- dplyr::bind_rows(MTC_results)
print(final_MTC_table_countries)
writexl::write_xlsx(final_MTC_table_scandi,"../SupplementaryFiles/MTC_sensitivities/MTC_countries.xlsx")

final_MTC_table_countries$Year<-final_MTC_table_countries$TimePeriod
final_MTC_table_countries$TimePeriod<-as.numeric(final_MTC_table_countries$TimePeriod)
ggplot(final_MTC_table_countries, aes(TimePeriod, MTC)) + 
  geom_line(aes(y = MTC), colour = "#00CCCC", linewidth = 1.2, alpha = 0.9) + 
  geom_ribbon(aes(y = MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "#00CCCC", size = 1.2, alpha = 0.3) + 
  theme_bw() + ggtitle("All species; 100-yr time bins; 200 simulations") + 
  facet_wrap(~ Type, scale = "free_y")

png(filename="../SupplementaryFiles/MTC_scandi.png", 
    width=10.5, height=5.5, units="in",res=600)
ggplot(final_MTC_table_countries, aes(TimePeriod, MTC, color = Type)) + 
  geom_line(aes(y = MTC, color = Type), linewidth = 1.2, alpha = 0.9) + 
  geom_ribbon(aes(y = MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = Type),color = NA, size = 1.2, alpha = 0.3) + 
  theme_bw() + ggtitle("Scandinavian countries; 100-yr time bins; 200 simulations") + 
  facet_wrap(~Type, scale = "free_y") + 
  scale_color_manual(values=DB[2:6]) + 
  scale_fill_manual(values=DB[2:6])
dev.off()

####Regardless of life history
NSS_sp_noaqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                             !NSS_sp$DB_Assemblage_ID == "NA", 
                                             !NSS_sp$NISP < 1.5,
                                             NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE")

NSS_sp_notrade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                                 !NSS_sp$DB_Assemblage_ID == "NA", 
                                                 !NSS_sp$NISP < 1.5,
                                                 NSS_sp$Known_Major_Trade_J_Barrett=="FALSE")

NSS_sp_noaquanotrade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                         !NSS_sp$DB_Assemblage_ID == "NA", 
                                         !NSS_sp$NISP < 1.5,
                                         NSS_sp$Known_Aquaculture_R_Hoffmann=="FALSE",
                                         NSS_sp$Known_Major_Trade_J_Barrett=="FALSE")

NSS_sp_trade<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                                   !NSS_sp$DB_Assemblage_ID == "NA", 
                                                   !NSS_sp$NISP < 1.5,
                                                   NSS_sp$Known_Major_Trade_J_Barrett=="TRUE")

NSS_sp_aqua<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                               !NSS_sp$DB_Assemblage_ID == "NA", 
                                               !NSS_sp$NISP < 1.5,
                                               NSS_sp$Known_Aquaculture_R_Hoffmann=="TRUE")

alldataframes<-list(NSS_sp_aqua, NSS_sp_trade, NSS_sp_noaqua, NSS_sp_notrade, NSS_sp_noaquanotrade)
names(alldataframes)<-c("NSS_sp_aqua","NSS_sp_trade", "NSS_sp_noaqua", "NSS_sp_notrade", "NSS_sp_noaquanotrade")
final <- data.frame(TimePeriod = character(), MTC = numeric(), SDMTC = numeric(), Run = character(), labs = character(), Type = character())
summary_results <- lapply(names(alldataframes), function(name) summarise_df(alldataframes[[name]], name))
final_summary_table_all <- dplyr::bind_rows(summary_results)
writexl::write_xlsx(final_summary_table_all,"../SupplementaryFiles/MTC_sensitivities/MTC_trade_aqua_counts.xlsx")
MTC_results <- lapply(names(alldataframes), function(name) MTC_df_global(alldataframes[[name]], name))
final_MTC_table_all <- dplyr::bind_rows(MTC_results)
print(final_MTC_table_countries)
writexl::write_xlsx(final_MTC_table_all,"../SupplementaryFiles/MTC_sensitivities/MTC_alldata_tradeaqua.xlsx")

final_MTC_table_all$Year<-final_MTC_table_all$TimePeriod
final_MTC_table_all$TimePeriod<-as.numeric(final_MTC_table_all$TimePeriod)
final_MTC_table_all$Dataset<-recode(final_MTC_table_all$Type,
                                    "NSS_sp_trade" = "Commercially traded sp.",
                                    "NSS_sp_notrade" = "Not commercially traded sp.",
                                    "NSS_sp_noaquanotrade" = "Not commercially traded or farmed sp.",
                                    "NSS_sp_noaqua" = "Not commercially farmed sp.",
                                    "NSS_sp_aqua" = "Commercially farmed sp.")
final_MTC_table_all$Dataset <- factor(final_MTC_table_all$Dataset, 
                                      levels = c(
                                        "Commercially traded sp.",
                                        "Commercially farmed sp.",
                                        "Not commercially traded sp.",
                                        "Not commercially farmed sp.",
                                        "Not commercially traded or farmed sp."
                                      ))
library(cowplot)
library(patchwork)
ggplot(final_MTC_table_all, aes(TimePeriod, MTC, color = Dataset)) + 
  geom_line(aes(y = MTC, color = Dataset), linewidth = 1.2, alpha = 0.9) + 
  geom_ribbon(aes(y = MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = Dataset),color = NA, size = 1.2, alpha = 0.3) + 
  theme_bw() + ggtitle("MTC of commercially farmed and traded species - 200 simulations") + 
  facet_wrap(~Dataset, ncol=2) + 
  scale_color_manual(values=DB[2:6]) + 
  scale_fill_manual(values=DB[2:6]) +
  theme(legend.position = "none")

p<-ggplot(final_MTC_table_all, aes(TimePeriod, MTC, color = Dataset)) + 
  geom_line(aes(y = MTC, color = Dataset), linewidth = 1.2, alpha = 0.9) + 
  geom_ribbon(aes(y = MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill = Dataset),color = NA, size = 1.2, alpha = 0.3) + 
  theme_bw() + ggtitle("Median Temperature of the Catch") + 
  facet_wrap(~Dataset, ncol=2) + 
  scale_color_manual(values=DB[2:6]) + 
  scale_fill_manual(values=DB[2:6]) + theme_publish()

p<-p %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1600, fill = "skyblue", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 1750, xmax = 1950, fill = "darkred", alpha = 0.15) 
png(filename="../SupplementaryFiles/MTC_farmed.png", 
    width=7.2, height=8.7, units="in",res=600)
p +  
  theme(legend.position = "none",strip.text = element_text(size=10)) +
  annotate("text", x = 1375, y = 7, label = "Little Ice Age", color = "grey5", size = 3, 
           family="Arial", fontface="italic") +
  annotate("text", x = 999, y = 7, label = "MCA", color = "grey5", size = 3,family="Arial", fontface="italic") +
#  annotate("text", x = 280, y = 7, label = "RWP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 420, y = 7, label = "DACP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 510, y = 15.97, label = "**", color = "grey1", size = 5.2,family="Arial", fontface="italic") + ylim(6.75,16.75)
dev.off()

