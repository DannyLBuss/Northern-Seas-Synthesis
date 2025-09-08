###Plot species-level data summaries using analysis datasets
####14th Aug 2025 - NSS analysis
#Author: Dr Danny L Buss
###4. Plot summary statistics of aggregated data, sieved vs. hand-collected, and dataset following outlier removals 
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
})
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("ThemePublish.R")
source("~/Documents/2.Academic_Work/3.NSS/2.MANUSCRIPT_Jan2025/Code_final/MANUSCRIPT_functions.R")

#2. load data ####
df_species<-read_excel("NSS_SpeciesData_Aug2025_noPE.xlsx")
df_speciesb<-read_excel("NSS_SpeciesData_Aug2025_withPE.xlsx")
DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69")
df_species$time_bins<-factor(df_species$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_species$time_bins2<-factor(df_species$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))

#3. Plot counts ####
df_species$DB_Assemblage_ID<-df_species$Lumped_NSS_IDs
df_species %>% group_by(Country) %>% 
  summarize(n_distinct(DB_Assemblage_ID))

df_species %>% group_by(DB_sieved) %>% 
  summarize(n_distinct(DB_Assemblage_ID))

df_species %>% group_by(Region) %>% 
  summarize(n_distinct(DB_Assemblage_ID))

ggplot(df_species, aes(time_bins, as.numeric(Time.mid), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=4) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none",
                                       axis.text.x = element_text(angle=45, hjust=1)) + 
  labs(y="", x="Time Period", title = "Number of assemblages in each time period") 

ggplot(df_species, aes(time_bins2, as.numeric(Time.mid), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=4) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none",
                                       axis.text.x = element_text(angle=45, hjust=1)) + 
  labs(y="", x="Time Period", title = "Number of assemblages in each time period") 

ggplot(df_species, aes(time_bins2, as.numeric(Time.mid), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=4) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none",
                                       axis.text.x = element_text(angle=45, hjust=1)) + 
  labs(y="", x="Time Period", title = "Number of assemblages in each time period") 

df_species$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df_species$DB_Assemblage_ID %>% unique %>% length() %>%
  sum(na.rm=T)

ggplot(df_species, aes(time_bins2, as.numeric(Time.mid), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=2) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none") + 
  labs(y="", x="Time Period", title = "Number of assemblages in each time period (Species-data, clean)") 

df_species$count<-"1"
ggplot(df_species, aes(time_bins2, as.numeric(count), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=2) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none") + labs(y="", x="Time Period", title = "Number of assemblages in each time period (Species-data, clean)") 

###Plot counts of assemblages overtime (mid-point)
df_species$time_bin3<-factor(cut(as.numeric(df_species$Time.mid), breaks=c(-800,-500,-300,1,101,201,301,401,501,601,701,801,901,1001,1101,1201,
                                                                           1301,1401,1501,1601,1701,1801,1901,2001,2200))) 
ggplot(df_species[as.numeric(df_species$Time.mid) > 1,], aes(time_bin3, as.numeric(count), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=2) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none") + labs(y="", x="Time Period", title = "Number of assemblages in each time period (Species-data, clean)") 

###Count assemblages & NISP
df_species$NISP_counts<-as.numeric(df_species$NISP)
df_species %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

df_species %>% 
  dplyr::select(NISP_counts, GBIF_species,DB_Assemblage_ID, count) %>% 
  filter(!GBIF_species == "NA") %>%
  group_by(count) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_species %>% 
  dplyr::select(GBIF_family,GBIF_species, GBIF_class, GBIF_order, GBIF_genus, count) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count) %>%
  dplyr::summarize(Families = length(unique(GBIF_family)),,
            Class = length(unique(GBIF_class)),
            Order = length(unique(GBIF_order)),
            Genera = length(unique(GBIF_genus)),
            Species = length(unique(GBIF_species)))

df_species %>% 
  dplyr::select(GBIF_family,GBIF_species, GBIF_class, GBIF_order, GBIF_genus, count, Region) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count, Region) %>%
  dplyr::summarize(Families = length(unique(GBIF_family)),,
            Class = length(unique(GBIF_class)),
            Order = length(unique(GBIF_order)),
            Genera = length(unique(GBIF_genus)),
            Species = length(unique(GBIF_species)))

#Prop sieved
df_species %>% 
  dplyr::select(NISP_counts, GBIF_family, DB_sieved, count, Region) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count, Region) %>%
  dplyr::group_by(DB_sieved, count, Region) %>%
  summarize(NISP = sum(NISP_counts),
            All = 690892,
            prop = NISP/All)

#Prop sieved by region
df_species %>%
  dplyr::select(NISP_counts, GBIF_family, DB_sieved, count, Region) %>%
  dplyr::filter(!is.na(GBIF_family), GBIF_family != "NA") %>%
  dplyr::group_by(DB_sieved, count, Region) %>%
  dplyr::summarize(NISP = sum(NISP_counts, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::group_by(Region) %>%
  dplyr::mutate(
    Region_total = sum(NISP),
    prop = NISP / Region_total
  ) %>%
  ungroup()

###Taxa-level counts:
df_species %>% 
  dplyr::select(time_bins2, NISP_counts, GBIF_species, Region) %>% 
  filter(!GBIF_species == "NA") %>%
  group_by(time_bins2,Region, GBIF_species) %>%
  summarize(NISP = sum(NISP_counts))

df_species %>% 
  dplyr::select(time_bins2, NISP_counts, GBIF_species,DB_Assemblage_ID, Region) %>% 
  filter(!GBIF_species == "NA") %>%
  group_by(time_bins2,Region) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_species %>% 
  dplyr::select(time_bins2, NISP_counts, GBIF_species,DB_Assemblage_ID, Region) %>% 
  filter(!GBIF_species == "NA") %>%
  group_by(time_bins2,Region) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_species %>% 
  dplyr::select(time_bins, NISP_counts, GBIF_species,DB_Assemblage_ID, Region) %>% 
  filter(!GBIF_species == "NA") %>%
  group_by(time_bins,Region) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_species %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

#Plot Figure 1 ####



