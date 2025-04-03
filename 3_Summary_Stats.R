###Plot species and family-level data summaries using analysis datasets
####24th Jan 2025 - NSS analysis
#Author: Dr Danny L Buss
###4. Plot summary statistics of aggregated data, sieved vs. hand-collected, and dataset following outlier removals 
#1. load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(reshape2)
library(bayesplot)
library(forcats)
library(stringr)
library(writexl)
library(scico)
source("~/Documents/4-Oceans/3.NSS/2.MANUSCRIPT_Jan2025/Code/ThemePublish.R")

#2. setwd & load data ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5a_NSS_SpeciesData_Feb2025_withCEE.xlsx")
df_speciesb<-read_excel("../SupplementaryFiles/SupplementaryTable5b_NSS_SpeciesData_Feb2025_noCEE.xlsx")
df_family<-read_excel("../SupplementaryFiles/SupplementaryTable6a_NSS_FamilyData_Feb2025_withCEE.xlsx")
df_familyb<-read_excel("../SupplementaryFiles/SupplementaryTable6b_NSS_FamilyData_Feb2025_noCEE.xlsx")
DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69")
df_family$time_bins<-factor(df_family$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_species$time_bins<-factor(df_species$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_family$time_bins2<-factor(df_family$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))
df_species$time_bins2<-factor(df_species$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))

#3. Plot counts ####
df_family$NISP %>% as.numeric() %>%
  sum(na.rm=T)
df_family$DB_Assemblage_ID %>% unique %>% length() %>%
  sum(na.rm=T)
df_family %>% group_by(Region) %>% 
  summarize(n_distinct(DB_Assemblage_ID))
df_family %>% group_by(Country) %>% 
  summarize(n_distinct(DB_Assemblage_ID))
df_species %>% group_by(Country) %>% 
  summarize(n_distinct(DB_Assemblage_ID))

df_species %>% group_by(DB_sieved) %>% 
  summarize(n_distinct(DB_Assemblage_ID))

ggplot(df_family, aes(time_bins2, as.numeric(Time.mid), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
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

df_family$count<-"1"
ggplot(df_family, aes(time_bins2, as.numeric(count), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=2) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none") + 
  labs(y="", x="Time Period", title = "Number of assemblages in each time period (Family-data, clean)") 

df_species$count<-"1"
ggplot(df_species, aes(time_bins2, as.numeric(count), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=2) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none") + labs(y="", x="Time Period", title = "Number of assemblages in each time period (Species-data, clean)") 

###Plot counts of assemblages overtime (mid-point)
df_family$time_bin3<-factor(cut(as.numeric(df_family$Time.mid), breaks=c(-800,-500,-300,1,101,201,301,401,501,601,701,801,901,1001,1101,1201,
                                                                         1301,1401,1501,1601,1701,1801,1901,2001,2200))) 
df_species$time_bin3<-factor(cut(as.numeric(df_species$Time.mid), breaks=c(-800,-500,-300,1,101,201,301,401,501,601,701,801,901,1001,1101,1201,
                                                                           1301,1401,1501,1601,1701,1801,1901,2001,2200))) 
ggplot(df_family[as.numeric(df_family$Time.mid) > 1,], aes(time_bin3, as.numeric(count), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=4) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle=45, hjust=1)) + 
  labs(y="", x="Time Period", title = "Number of assemblages in each time period (Family-data, clean)") 

ggplot(df_species[as.numeric(df_species$Time.mid) > 1,], aes(time_bin3, as.numeric(count), color=Region, fill=Region)) + geom_bar(stat="identity", position = "stack") + 
  facet_wrap(~Region, nrow=2) + theme_publish() + 
  scale_color_manual(values=DB) + 
  scale_fill_manual(values=DB) + theme(legend.position = "none") + labs(y="", x="Time Period", title = "Number of assemblages in each time period (Species-data, clean)") 

###Count assemblages & NISP, family-level data
df_family$NISP_counts<-as.numeric(df_family$NISP)
df_family %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

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

df_family %>% 
  dplyr::select(NISP_counts, GBIF_family,DB_Assemblage_ID, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(count) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_family %>% 
  select(NISP_counts, GBIF_family,DB_Assemblage_ID, Region, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(Region, count) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_species %>% 
  dplyr::select(NISP_counts, GBIF_family,DB_Assemblage_ID, Region, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(Region, count) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_species %>% 
  dplyr::select(GBIF_family,GBIF_species, GBIF_class, GBIF_order, GBIF_genus, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(count) %>%
  summarize(Families = length(unique(GBIF_family)),,
            Class = length(unique(GBIF_class)),
            Order = length(unique(GBIF_order)),
            Genera = length(unique(GBIF_genus)),
            Species = length(unique(GBIF_species)))

df_species %>% 
  dplyr::select(GBIF_family,GBIF_species, GBIF_class, GBIF_order, GBIF_genus, count, Region) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(count, Region) %>%
  summarize(Families = length(unique(GBIF_family)),,
            Class = length(unique(GBIF_class)),
            Order = length(unique(GBIF_order)),
            Genera = length(unique(GBIF_genus)),
            Species = length(unique(GBIF_species)))

df_family %>% 
  dplyr::select(GBIF_family,GBIF_class, GBIF_order, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(count) %>%
  summarize(Families = length(unique(GBIF_family)),,
            Class = length(unique(GBIF_class)),
            Order = length(unique(GBIF_order)))

df_family %>% 
  dplyr::select(GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(count, Region) %>%
  summarize(Families = length(unique(GBIF_family)),,
            Class = length(unique(GBIF_class)),
            Order = length(unique(GBIF_order)))

#Prop sieved
df_family %>% 
  dplyr::select(NISP_counts, GBIF_family, DB_sieved, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(DB_sieved, count) %>%
  summarize(NISP = sum(NISP_counts),
            All = 1097819,
            prop = NISP/All)

df_species %>% 
  dplyr::select(NISP_counts, GBIF_family, DB_sieved, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(DB_sieved, count) %>%
  summarize(NISP = sum(NISP_counts),
            All = 829411,
            prop = NISP/All)

#Prop sieved per Region
df_species %>% 
  dplyr::select(NISP_counts, GBIF_family, DB_sieved, Region, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(DB_sieved, Region, count) %>%
  summarize(NISP = sum(NISP_counts)) %>%
  group_by(Region, count) %>%
  summarize(All = sum(NISP),
            prop = NISP/All)

df_family %>% 
  dplyr::select(NISP_counts, GBIF_family, DB_sieved, Region, count) %>% 
  filter(!GBIF_family == "NA") %>%
  group_by(DB_sieved, Region, count) %>%
  summarize(NISP = sum(NISP_counts)) %>%
  group_by(Region, count) %>%
  summarize(All = sum(NISP),
            prop = NISP/All)


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
  select(time_bins, NISP_counts, GBIF_species,DB_Assemblage_ID) %>% 
  filter(!GBIF_species == "NA") %>%
  group_by(time_bins) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_species %>% 
  select(time_bins2, NISP_counts, GBIF_species,DB_Assemblage_ID) %>% 
  filter(!GBIF_species == "NA") %>%
  group_by(time_bins2) %>%
  summarize(NISP = sum(NISP_counts),
            Assemblages = length(unique(DB_Assemblage_ID)))

df_species %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

###aMTC using time bins
df_family %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

