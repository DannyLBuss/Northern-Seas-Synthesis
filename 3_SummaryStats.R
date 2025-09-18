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
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(tidyverse)
  library(ggpubr)
})
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../ThemePublish.R")
source("MANUSCRIPT_functions.R")

#2. load data ####
df_species<-read_xlsx("NSS_SpeciesData_Aug2025_noPE.xlsx", col_types = "text")
df_speciesb<-read_xlsx("NSS_SpeciesData_Aug2025_withPE.xlsx", col_types = "text")
DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69")
df_species$time_bins<-factor(df_species$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_species$time_bins2<-factor(df_species$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))

#3. Plot counts ####
#df_species$DB_Assemblage_ID<-df_species$Lumped_NSS_IDs
df_species %>% dplyr::group_by(Database) %>% 
  dplyr::summarize(dplyr::n_distinct(DB_Assemblage_ID))

df_species %>% dplyr::group_by(Country) %>% 
  dplyr::summarize(dplyr::n_distinct(DB_Assemblage_ID))

df_species %>% dplyr::group_by(DB_sieved) %>% 
  dplyr::summarize(n_distinct(DB_Assemblage_ID))

df_species %>% dplyr::group_by(Region) %>% 
  dplyr::summarize(n_distinct(DB_Assemblage_ID))

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

tot<-df_species %>% 
  dplyr::summarize(NISP = sum(NISP_counts))

#Prop sieved
df_species %>% 
  dplyr::select(NISP_counts, GBIF_family, DB_sieved, count, Region) %>% 
  dplyr::filter(!GBIF_family == "NA") %>%
  dplyr::group_by(count, Region) %>%
  dplyr::group_by(DB_sieved, count, Region) %>%
  summarize(NISP = sum(NISP_counts),
            All =  tot$NISP,
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
df<-read_excel("~/Documents/2.Academic_Work/3.NSS/2.MANUSCRIPT_Jan2025/Final_code_runthrough/NSS_Workspace/Sept2025/All_Data_with_Chronology_Aug2025.xlsx")
head(df)
str(df)
sum(as.numeric(df$NISP))
t<-as.numeric(df$NISP)
sum(t, na.rm=TRUE)

###Figure 1 ####
DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69"
)

##Plot map with three region colours
df$region<-dplyr::recode(as.character(df$Country),"England" = 'Britain & Ireland',
                         "Scotland" = 'Britain & Ireland',
                         "Norway" = 'Scandinavia',
                         "Sweden" = 'Scandinavia',
                         "Denmark" = 'Scandinavia',
                         "Germany" = 'Western Europe',
                         "Belgium" = 'Western Europe',
                         "Netherlands" = 'Western Europe',
                         "Ireland" = 'Britain & Ireland',
                         "Estonia" = 'Poland & Estonia',
                         "Poland" = 'Poland & Estonia',
                         "Northern Ireland" = 'Britain & Ireland')

europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
plot_eur<-df %>% dplyr::select(`Decimal Longitude`, `Decimal Latitude`,region) %>% dplyr::distinct()
plot_eur$`Decimal Longitude`<-as.numeric(plot_eur$`Decimal Longitude`)/1000
plot_eur$`Decimal Latitude`<-as.numeric(plot_eur$`Decimal Latitude`)/1000
names(plot_eur)<-c("longitude","latitude","Region")

seas <- data.frame(
  lon = c(3.0, 19.0, 5.0), 
  lat = c(56.0, 58.5, 65.0), 
  name = c("North Sea", "Baltic", "Norwegian Sea") # City names
)

Oceans <- data.frame(
  lon = c(-15.5), 
  lat = c(50.0), 
  name = c("NORTH ATLANTIC") 
)

ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  
  geom_point(data = plot_eur, 
             aes(x = longitude, y = latitude, color = Region),
             size = 1) +  
  scale_color_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                                "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) +
  theme_minimal() +
  labs(title = "A)",
       x = "Longitude", y = "Latitude", color = "Region") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "none") 
##A - Map by region
Panel_A<-ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_eur, 
             aes(x = longitude, y = latitude, color = Region),
             size = 1.7) +  # Points with categories
  scale_color_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                                "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) +
  theme_minimal() +
  labs(title = "",
       x = "Longitude", y = "Latitude", color = "Region") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "none") +
  geom_text(data = seas, aes(x = lon, y = lat, label = name), color = "grey15", size = 3, fontface="italic", family="serif") +
  geom_text(data = Oceans, aes(x = lon, y = lat, label = name), color = "grey30", size = 3, fontface="italic", family="serif")

##Panel B - proportion of sites at each region
df_tmp<-df %>% dplyr::select(Lumped_NSS_IDs, region) %>% distinct()
nrow(df_tmp)
df_tmp$all<-"N = 2168"
df_tmp %>%
  ggplot(aes(x = all, fill = region, color = region)) +
  geom_bar(stat="count", position = "fill") + theme_publish() +
  labs(title = "C)", x = "", y = "Proportion of assemblages") +
  scale_color_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                                "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) +
  scale_fill_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                               "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) 
df$Site<-as.factor(df$`Site name`)
df_tmp<-df %>% dplyr::select(Site, region) %>% distinct()
nrow(df_tmp)
df_tmp$all<-"N = 912"
df_tmp %>%
  ggplot(aes(x = all, fill = region)) +
  geom_bar(stat="count", position = "fill") + theme_publish() +
  labs(title = "", x = "", y = "Proportion of archaeological sites", fill="") +
  scale_fill_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                               "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) +
  theme(legend.position = "right")

Panel_B<-df_tmp %>%
  ggplot(aes(x = all, fill = region)) +
  geom_bar(stat="count", position = "fill") + theme_publish() +
  labs(title = "", x = "", y = "Proportion of archaeological sites", fill="") +
  scale_fill_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                               "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) +
  theme(legend.position = "none")

###Panel C - NISP at each region
df_tmp<-df %>% dplyr::select(NISP, region, Lumped_NSS_IDs)
df_tmp$NISP<-as.numeric(df_tmp$NISP)
df_tmp<-df_tmp[!is.na(df_tmp$NISP),]
Panel_C<-df_tmp[df_tmp$NISP > 50,] %>%
  ggplot(aes(x = region, y = log(NISP), fill=region, color=region)) +
  geom_jitter(alpha=0.25,aes(fill=region, color=region)) + 
  geom_boxplot(aes(fill=region), color="grey30", outliers=F) +
  theme_publish() +
  labs(title = "", x = "", y = "log(NISP)", fill="") +
  scale_fill_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                               "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) +
  scale_color_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                                "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) +
  theme(legend.position = "none",axis.text.x = element_text(angle=45, hjust=1))


##D - Chron uncertainity by region
df$Start<-as.numeric(df$`Start date CE`)
df$End<-as.numeric(df$`End date CE`)
df$time.range<-as.numeric(df$`End date CE`) - as.numeric(df$`Start date CE`)
df$Assemblage<-as.factor(df$`Assemblage or sub-assemblage`)
df$all<-"all"

df_tmp<-df %>% dplyr::select(Start, End, `Decimal Latitude`, `Decimal Longitude`, region, time.range) %>% distinct()
df_tmp<-df_tmp %>%
  filter(!is.na(time.range)) %>%
  filter(Start > 0) %>%
  filter(End > 0) %>%
  filter(time.range < 401) %>%
  filter(time.range > 10)
df_tmp$Assemblage<-paste(df_tmp$`Decimal Latitude`,df_tmp$`Decimal Longitude`)
str(df_tmp)
df_tmp<-df_tmp[order(df_tmp$Start),]
Panel_D<-ggplot(df_tmp, aes(color=region,fill=region)) +
  geom_segment(linewidth=6,alpha=0.8,
               aes(color=region,
                   x=Start,
                   xend=End,
                   y=fct_reorder(Assemblage,Start, min),
                   yend=fct_reorder(Assemblage, Start, min))) + 
  labs(x="Year", y ="Assemblage ID", title="") +
  scale_fill_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                               "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) +
  scale_color_manual(values = c("Britain & Ireland" = "#14517B", "Western Europe" = "#6A8A73", 
                                "Scandinavia" = "#FE994F", "Poland & Estonia" = "#8e7d69")) + theme_publish() +
  theme(axis.text.y = element_text(size=0.3),
        legend.position="none")
Panel_D<-Panel_D + facet_wrap(~region, nrow=1)
Panel_D

top_left <- ggarrange(
  Panel_A, Panel_B,
  ncol = 2, widths = c(3, 1), 
  labels = c("A", "B"), align = "h"
)

left_col <- ggarrange(
  top_left, Panel_D,
  ncol = 1, heights = c(1, 1),
  labels = c("", "D"), align = "v"
)

final_plot <- ggarrange(
  left_col, Panel_C,
  ncol = 2, widths = c(2, 1),
  labels = c("", "C")
)

layout_matrix = rbind(c(1, 1, 1, 2, NA, NA),
                      c(1, 1, 1, 2, 4, 4),
                      c(1, 1, 1, 2, 4, 4),
                      c(1, 1, 1, 2, 4, 4),
                      c(3, 3, 3, NA, 4, 4),
                      c(3, 3, 3, NA, 4, 4))
png(filename="Figure_1_Sept25.png", 
    width=11, height=5, units="in",res=600)
gridExtra::grid.arrange(Panel_A, Panel_B, Panel_D, Panel_C,
                        layout_matrix = layout_matrix)
dev.off()
