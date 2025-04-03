###General Additive Modelling for fisheries
#Plot figure 1 usng all data - K = 4
library(ggplot2)
library(readxl)
library(forcats)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyverse)
library(scico)
library(gganimate)
library(reshape2)
library(zoo)
library(writexl)
library(mgcv)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
df_meta<-read_xlsx("~/Documents/4-Oceans/3.NSS/2.MANUSCRIPT_Jan2025/Data/SupplementaryTable5b_NSS_SpeciesData_Feb2025_noCEE.xlsx", col_names = T)
df_Assemblage<-read_xlsx("~/Documents/4-Oceans/3.NSS/2.MANUSCRIPT_Jan2025/Data/RecommenderSystems/Species_RS_outputs/Xmat5_4_30looprun.xlsx", col_names=T)
names(df_Assemblage)<-c("DB_Assemblage_ID","V1","V2","V3","V4")
df_Species<-read_xlsx("~/Documents/4-Oceans/3.NSS/2.MANUSCRIPT_Jan2025/Data/RecommenderSystems/Species_RS_outputs/Ymat5_4_30looprun.xlsx", col_names=T)
names(df_Species)<-c("GBIF_species","V1","V2","V3","V4")
df_Species<-df_Species[-1,]
source("~/Documents/4-Oceans/3.NSS/2.MANUSCRIPT_Jan2025/Code_final/MANUSCRIPT_functions.R")

df_Species<-df_Species[complete.cases(df_Species),]
df_meta_spec<-df_meta %>% dplyr::select(GBIF_species,GBIF_family,GBIF_order,GBIF_class,GBIF_phylum,GBIF_kingdom, Trait_Habitat, 
                                        Trait_LifeHistory, Trait_Temp_Max, Trait_Temp_Min,
                                        Known_Major_Trade_J_Barrett, Known_Aquaculture_R_Hoffmann) %>% distinct()
GBIF_species_missing<-df_meta_spec$GBIF_species[!df_meta_spec$GBIF_species%in%df_Species$GBIF_species]
df_meta_spec$GBIF_species<-gsub(" ","_",df_meta_spec$GBIF_species)
df_meta_spec$GBIF_species<-gsub("\\/",".",df_meta_spec$GBIF_species)
df_Ymat<-df_Species %>% left_join(df_meta_spec)
df_meta_assem<-df_meta %>% dplyr::select(DB_Assemblage_ID, Site_name,Settlement_modern_name,Country,Region,Decimal_Latitude,Decimal_Longitude,
                                         time_bins2,Start_date_CE, End_date_CE, Rural_urban_or_neither, Time.mid, time_bins, Time_range,
                                         Earliest_relevant_urban_date_in_NSS,Historic_Urban_Centre,Local_site_type, Recovery, DB_sieved) %>% distinct()
df_meta_assem$DB_Assemblage_ID[!df_meta_assem$DB_Assemblage_ID%in%df_Assemblage$DB_Assemblage_ID]
df_Xmat<-df_Assemblage %>% left_join(df_meta_assem)
df_Species_out<-df_Ymat

###Figure 1 
DB<-c("#FE994F","#6A8A73","#14517B","#8e7d69",
      "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
      "#EDC948", "#AF7AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
      "#A0CBE8", "#FFBE7D", "#FF9D9A", "#8CD17D", "#B6992D",
      "#D4A6C8", "#FABFD2", "#79706E", "#D37295", "#C3C3C3"
)

##Plot map with three region colours
df_Assemblage_out<-df_Xmat
europe <- ne_countries(scale = "medium", continent = "Europe", returnclass = "sf")
df_Assemblage_out$Decimal_Longitude<-as.character(df_Assemblage_out$Decimal_Longitude)
df_Assemblage_out$Decimal_Latitude<-as.character(df_Assemblage_out$Decimal_Latitude)
str(df_Assemblage_out)
df_Assemblage_out<-df_Assemblage_out %>%
  dplyr::mutate(across(c(V1,V2,V3,V4),as.numeric))
df_Assemblage_out<-df_Assemblage_out %>%
  dplyr::filter(!is.na(V1))
plot_eur<-melt(df_Assemblage_out %>% dplyr::select(DB_Assemblage_ID, Decimal_Longitude, Decimal_Latitude, Region, time_bins2, Recovery, DB_sieved, V1, V2, V3, V4))
names(plot_eur)[8]<-"fishery"
plot_eur$Decimal_Longitude<-as.numeric(plot_eur$Decimal_Longitude)
plot_eur$Decimal_Latitude<-as.numeric(plot_eur$Decimal_Latitude)
plot_eur<-plot_eur[,c(2:9)]
names(plot_eur)<-c("longitude","latitude","region","time_period","recovery","sieved","fishery","value")
plot_eur$fishery<-as.character(plot_eur$fishery)
plot_eur$time_period<-factor(plot_eur$time_period, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))

seas <- data.frame(
  lon = c(3.0, 19.0, 5.0), # Longitude
  lat = c(56.0, 58.5, 65.0), # Latitude
  name = c("North Sea", "Baltic", "Norwegian Sea") # City names
)

Oceans <- data.frame(
  lon = c(-15.5), # Longitude
  lat = c(50.0), # Latitude
  name = c("NORTH ATLANTIC") # City names
)

#Increasing
plot_eur<-plot_eur[!plot_eur$time_period=="NA",]
ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_eur %>%
               filter(value > 0.01), 
             aes(x = longitude, y = latitude, color = value),
             size = 1) +  # Points with categories
  scale_color_gradient(low = "white", high = "darkred") +  
  theme_minimal() +
  labs(title = "Propensity of sites for specific fishery types",
       x = "Longitude", y = "Latitude", color = "fishery") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "bottom") + facet_grid(time_period~fishery)
#Decreasing
ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_eur %>%
               filter(value < -0.01), 
             aes(x = longitude, y = latitude, color = value),
             size = 1) +  # Points with categories
  scale_color_gradient(low = "navyblue", high = "white") +  
  theme_minimal() +
  labs(title = "Propensity of sites for specific fishery types",
       x = "Longitude", y = "Latitude", color = "fishery") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "bottom") + facet_grid(time_period~fishery)

###Animation of Fishery Topic 4
plot_anim<-melt(df_Assemblage_out %>% 
                  dplyr::select(DB_Assemblage_ID, 
                                Start_date_CE, End_date_CE, 
                                Decimal_Longitude, Decimal_Latitude, Rural_urban_or_neither, Historic_Urban_Centre,
                                Region, time_bins2, Recovery, DB_sieved, V1, V2, V3, V4))
names(plot_anim)[12]<-"fishery"
plot_anim$Decimal_Longitude<-as.numeric(plot_anim$Decimal_Longitude)
plot_anim$Decimal_Latitude<-as.numeric(plot_anim$Decimal_Latitude)
plot_anim$start<-as.numeric(plot_anim$Start_date_CE)
plot_anim$end<-as.numeric(plot_anim$End_date_CE)
names(plot_anim)<-c("id","start","end","longitude","latitude","Rural_Urban","Historic_Urban_Centre","region","time_period","recovery","DB_sieved","fishery","value","StartDate","EndDate")
plot_anim$fishery<-as.character(plot_anim$fishery)
plot_anim$time_period<-factor(plot_anim$time_period, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))
plot_anim<-plot_anim[complete.cases(plot_anim),]
all_years<-seq(min(plot_anim$StartDate),max(plot_anim$EndDate), by=1)
plot_anim<-plot_anim[complete.cases(plot_anim),]
plot_anim$split<-cut(plot_anim$value,
                     breaks = c(-Inf,0,Inf),
                     labels = c("unlikely","likely"),
                     right = FALSE)

ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_anim,
             aes(x = longitude, y = latitude, color = value),
             size = 0.5) +  # Points with categories
  scale_color_gradient(low = "blue", high = "red") +  
  theme_minimal() +
  labs(title = "Propensity of sites for specific fishery types",
       x = "Longitude", y = "Latitude", color = "fishery") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "bottom") + facet_grid(fishery~split)
plot_anim$grouping<-paste(plot_anim$split, plot_anim$time_period)
plot_anim$grouping<-factor(plot_anim$grouping, 
                           levels=c("likely <500","likely 500-700","likely 700-900",
                                    "likely 900-1100","likely 1100-1300","likely 1300-1500","likely 1500-1700","likely >1700",
                                    "unlikely <500","unlikely 500-700","unlikely 700-900",
                                    "unlikely 900-1100","unlikely 1100-1300","unlikely 1300-1500","unlikely 1500-1700","unlikely >1700"))

ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_anim,
             aes(x = longitude, y = latitude, color = value),
             size = 0.5) +  # Points with categories
  scale_color_gradient(low = "blue", high = "red") +  
  theme_minimal() +
  labs(title = "Propensity of sites for specific fishery types",
       x = "Longitude", y = "Latitude", color = "fishery") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "bottom") + facet_grid(fishery~split)

ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_anim,
             aes(x = longitude, y = latitude, color = value),
             size = 0.5) +  # Points with categories
  scale_color_gradient(low = "blue", high = "red") +  
  theme_minimal() +
  labs(title = "Propensity of sites for specific fishery types",
       x = "Longitude", y = "Latitude", color = "fishery") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "bottom") + facet_grid(fishery~recovery)

ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_anim,
             aes(x = longitude, y = latitude, color = value),
             size = 0.5) +  # Points with categories
  scale_color_gradient(low = "blue", high = "red") +  
  theme_minimal() +
  labs(title = "Propensity of sites for specific fishery types",
       x = "Longitude", y = "Latitude", color = "fishery") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "bottom") + facet_grid(fishery~DB_sieved)

Fig_1A.1<-ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_anim[plot_anim$split=="likely",],
             aes(x = longitude, y = latitude, color = value),
             size = 0.5) +  # Points with categories
  scale_color_gradient(low = "white", high = "red") +  
  theme_minimal() +
  labs(title = "Propensity of sites for specific fishery types",
       x = "Longitude", y = "Latitude", color = "fishery") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "bottom") + facet_grid(fishery~grouping)
Fig_1A.1

Fig_1A.2<-ggplot(data = europe) +
  geom_sf(fill = "grey90", color = "black") +  # Europe map
  geom_point(data = plot_anim[plot_anim$split=="unlikely",],
             aes(x = longitude, y = latitude, color = value),
             size = 0.5) +  # Points with categories
  scale_color_gradient(low = "blue", high = "white") +  
  theme_minimal() +
  labs(title = "Propensity of sites for specific fishery types",
       x = "Longitude", y = "Latitude", color = "fishery") +
  coord_sf(xlim = c(-23, 33), ylim = c(37, 74)) + theme(panel.grid.major = element_blank(), # Remove major gridlines
                                                        panel.grid.minor = element_blank(),
                                                        legend.position = "bottom") + facet_grid(fishery~grouping)
Fig_1A.2

plot_anim<-plot_anim %>%
  dplyr::group_by(fishery) %>%
  dplyr::mutate(scaled_value = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

Fig_1B.1<-ggplot(plot_anim,aes(time_period, scaled_value, fill=fishery, color=fishery)) + 
  geom_point() + 
  geom_boxplot() + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(region~fishery)

plot_anim$rur_reg<-paste(plot_anim$Rural_Urban, plot_anim$region)
ggplot(plot_anim[!plot_anim$Rural_Urban=="neither",],aes(time_period, scaled_value, fill=fishery, color=fishery, shape=Rural_Urban)) + 
  geom_point() + 
  geom_boxplot() + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(rur_reg~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB)

ggplot(plot_anim[!plot_anim$Historic_Urban_Centre=="NA",],aes(time_period, scaled_value, fill=fishery, color=fishery, shape=Historic_Urban_Centre)) + 
  geom_point() + 
  geom_boxplot() + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(Historic_Urban_Centre~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB)

plot_anim$his_reg<-paste(plot_anim$Historic_Urban_Centre, plot_anim$region)
ggplot(plot_anim[!plot_anim$Historic_Urban_Centre=="NA",],aes(time_period, scaled_value, fill=fishery, color=fishery, shape=his_reg)) + 
  geom_point() + 
  geom_boxplot() + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(his_reg~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB)

###Expand rowwise and plot by year
plot_anim_expanded<- plot_anim %>%
  dplyr::rowwise() %>%
  dplyr::mutate(year = list(seq(StartDate, EndDate, by=25))) %>%
  unnest(year) %>%
  ungroup()
names(plot_anim_expanded)

df_moving_avg<-plot_anim_expanded[,c(6:12,13,18,16,21)] %>%
  group_by(fishery) %>%
  arrange(year) %>%
  mutate(moving_avg = rollapply(scaled_value, width=3, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 25 == 0) %>%
  ungroup()

df_moving_avg<-df_moving_avg[complete.cases(df_moving_avg),]
df_moving_avg$his_reg<-paste(df_moving_avg$Historic_Urban_Centre, df_moving_avg$region)
ggplot(df_moving_avg[!df_moving_avg$Historic_Urban_Centre=="NA",],
       aes(year, moving_avg, fill=fishery, color=fishery)) + 
  geom_point(alpha=0.1, shape=19, size=0.95) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(his_reg~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) 


df_moving_avg2 <- plot_anim_expanded %>%
  dplyr::select(fishery, scaled_value, his_reg, year) %>% 
  group_by(his_reg,fishery) %>%
  arrange(year) %>%
  mutate(moving_avg = rollapply(scaled_value, width=100, FUN=mean, fill=NA, align="center")) %>%
  #  filter(year %% 100 == 0) %>%
  ungroup()

ggplot(df_moving_avg[!df_moving_avg$Historic_Urban_Centre=="NA",],
       aes(year, moving_avg, fill=fishery, color=fishery)) + 
  geom_point(alpha=0.1, shape=19, size=0.95) + 
  geom_point(alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(his_reg~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) + geom_smooth( method = "loess",
                                               data=df_moving_avg2, aes(x = year, y = moving_avg, color = fishery), size = 1, span=0.2
  )
plot_anim_expanded$year<-as.numeric(plot_anim_expanded$year)
df_moving_avg3 <- plot_anim_expanded %>%
  dplyr::select(fishery, scaled_value, region, year) %>% 
  group_by(region, fishery) %>%
  arrange(year) %>%
  mutate(moving_avg = rollapply(scaled_value, width=20, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg3$year<-as.character(df_moving_avg3$year)
df_moving_avg3<-df_moving_avg3 %>% 
  dplyr::group_by(fishery, region, year) %>%
  summarise(
    meanavg = mean(moving_avg, na.rm = TRUE),
    sdavg = sd(moving_avg, na.rm = TRUE)
  )
df_moving_avg3$year<-as.numeric(df_moving_avg3$year)

df_moving_avgtmp<-df_moving_avg[!df_moving_avg$Historic_Urban_Centre=="NA",]
str(df_moving_avgtmp)
df_moving_avg3<-as.data.frame(df_moving_avg3[,c(1:5)])
a<-ggplot() + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=moving_avg, fill=fishery, color=fishery), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=moving_avg), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(region~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) + geom_point(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, y = meanavg, color = fishery), size = 1)

a + geom_errorbar(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg, color = fishery), size = 1, span=0.2) 

ggplot(df_moving_avg[!df_moving_avg$Historic_Urban_Centre=="NA",],
       aes(year, moving_avg, fill=fishery, color=fishery)) + 
  geom_point(alpha=0.1, shape=19, size=0.95) + 
  geom_point(alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(his_reg~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) + geom_smooth( method = "loess",
                                               data=df_moving_avg2, aes(x = year, y = moving_avg, color = fishery), size = 1, span=0.2
  )

##Plot for Rural Urban ####
df_moving_avg3 <- plot_anim_expanded %>%
  dplyr::select(fishery, scaled_value, Rural_Urban, year) %>% 
  group_by(Rural_Urban, fishery) %>%
  arrange(year) %>%
  mutate(moving_avg = rollapply(scaled_value, width=20, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg3$year<-as.character(df_moving_avg3$year)
df_moving_avg3<-df_moving_avg3 %>% 
  dplyr::group_by(fishery, Rural_Urban, year) %>%
  summarise(
    meanavg = mean(moving_avg, na.rm = TRUE),
    sdavg = sd(moving_avg, na.rm = TRUE)
  )
df_moving_avg3$year<-as.numeric(df_moving_avg3$year)

df_moving_avgtmp<-df_moving_avg[!df_moving_avg$Rural_Urban=="NA",]
str(df_moving_avgtmp)
df_moving_avg3<-as.data.frame(df_moving_avg3[,c(1:5)])
b<-ggplot() + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=moving_avg, fill=fishery, color=fishery), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=moving_avg), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(Rural_Urban~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) + geom_point(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, y = meanavg, color = fishery), size = 1)

b + geom_errorbar(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg, color = fishery), size = 1, span=0.2) 

##Plot for historic urban centres ####
df_moving_avg3 <- plot_anim_expanded %>%
  dplyr::select(fishery, scaled_value, Rural_Urban, year) %>% 
  group_by(Rural_Urban, fishery) %>%
  arrange(year) %>%
  mutate(moving_avg = rollapply(scaled_value, width=20, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg3$year<-as.character(df_moving_avg3$year)
df_moving_avg3<-df_moving_avg3 %>% 
  dplyr::group_by(fishery, Rural_Urban, year) %>%
  summarise(
    meanavg = mean(moving_avg, na.rm = TRUE),
    sdavg = sd(moving_avg, na.rm = TRUE)
  )
df_moving_avg3$year<-as.numeric(df_moving_avg3$year)

df_moving_avgtmp<-df_moving_avg[!df_moving_avg$Rural_Urban=="NA",]
df_moving_avg3<-as.data.frame(df_moving_avg3[,c(1:5)])
b<-ggplot() + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=moving_avg, fill=fishery, color=fishery), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=moving_avg), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(Rural_Urban~fishery) +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) + geom_point(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, y = meanavg, color = fishery), size = 1)

b + geom_errorbar(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg, color = fishery), size = 1, span=0.2) 

##Plot for splits using weighted averages ####
plot_anim_expanded<- plot_anim %>%
  dplyr::rowwise() %>%
  dplyr::mutate(year = list(seq(StartDate, EndDate, by=25))) %>%
  unnest(year) %>%
  ungroup()
names(plot_anim_expanded)

plot_anim_expanded<-plot_anim_expanded %>%
  group_by(year) %>%
  mutate(weight = n()) %>%
  ungroup()

window_size <- 3

df_moving_avg<-plot_anim_expanded[,c(6:12,13,18,16,21,22)] %>%
  group_by(fishery, year) %>%
  arrange(year) %>%
  mutate(
    weighted_rolling_avg = rollapply(
      data = value,
      width = window_size,  # Window size for rolling average
      FUN = function(x) weighted.mean(x, w = rep(weight[1], length(x))),
      fill = NA,  # Handle missing values in rolling windows
      align = "center"  # Align rolling window to center
    )
  ) %>%
  ungroup()
df_moving_avg<-df_moving_avg[complete.cases(df_moving_avg),]
df_moving_avg$his_reg<-paste(df_moving_avg$Historic_Urban_Centre, df_moving_avg$region)

window_size <- 20
df_moving_avg3<-plot_anim_expanded %>%
  dplyr::select(fishery, value, year, weight) %>% 
  arrange(year) %>%
  mutate(
    weighted_rolling_avg = rollapply(
      data = value,
      width = window_size,  # Window size for rolling average
      FUN = function(x) weighted.mean(x, w = rep(weight[1], length(x))),
      fill = NA,  # Handle missing values in rolling windows
      align = "center"  # Align rolling window to center
    )
  ) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg3<-df_moving_avg3[complete.cases(df_moving_avg3),]

df_moving_avg3$year<-as.character(df_moving_avg3$year)
df_moving_avg3<-df_moving_avg3 %>% 
  dplyr::group_by(fishery, year) %>%
  summarise(
    meanavg = mean(weighted_rolling_avg, na.rm = TRUE),
    sdavg = sd(weighted_rolling_avg, na.rm = TRUE)
  )
df_moving_avg3$year<-as.numeric(df_moving_avg3$year)
df_moving_avgtmp<-df_moving_avg[!df_moving_avg$split=="NA",]
str(df_moving_avgtmp)
df_moving_avg3<-as.data.frame(df_moving_avg3[,c(1:4)])
b<-ggplot() + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=weighted_rolling_avg, fill=fishery, color=fishery), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=weighted_rolling_avg), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~fishery, scales="free") +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) + geom_point(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, y = meanavg, color = fishery), size = 1)

b + geom_errorbar(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg, color = fishery), size = 1, span=0.2) + theme_publish() +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5, color = "darkred") +
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),  # Major grid lines
        panel.grid.minor = element_line(color = "grey", size = 0.2),  # Minor grid lines
        panel.grid.major.x = element_blank(),  # Removing major grid lines on x-axis
        panel.grid.minor.x = element_blank(),  # Removing minor grid lines on x-axis
        panel.background = element_rect(fill = "white", color = "white"),  # White background
        plot.background = element_rect(fill = "white", color = "white"),  # White plot background
        axis.text = element_text(size = 12),  # Adjusting axis text size
        axis.title = element_text(size = 14),  # Adjusting axis titles size
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
        legend.position = "none") + labs(x = "Year (CE)", y = "Propensity for fishery") +
  scale_x_continuous(expand= c(0,0), limits=c(1,1850))

###Split by region
df_moving_avg3<-plot_anim_expanded %>%
  dplyr::select(fishery, value, year, weight, region) %>% 
  arrange(year) %>%
  mutate(
    weighted_rolling_avg = rollapply(
      data = value,
      width = window_size,  # Window size for rolling average
      FUN = function(x) weighted.mean(x, w = rep(weight[1], length(x))),
      fill = NA,  # Handle missing values in rolling windows
      align = "center"  # Align rolling window to center
    )
  ) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg3<-df_moving_avg3[complete.cases(df_moving_avg3),]

df_moving_avg3$year<-as.character(df_moving_avg3$year)
df_moving_avg3<-df_moving_avg3 %>% 
  dplyr::group_by(fishery, year, region) %>%
  summarise(
    meanavg = mean(weighted_rolling_avg, na.rm = TRUE),
    sdavg = sd(weighted_rolling_avg, na.rm = TRUE)
  )
df_moving_avg3$year<-as.numeric(df_moving_avg3$year)
df_moving_avgtmp<-df_moving_avg[!df_moving_avg$split=="NA",]
str(df_moving_avgtmp)
str(df_moving_avg3)
df_moving_avg3<-as.data.frame(df_moving_avg3[,c(1:5)])
b<-ggplot() + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=weighted_rolling_avg, fill=fishery, color=fishery), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avgtmp,
             aes(x=year, y=weighted_rolling_avg), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(region~fishery, scales="free") +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) + geom_point(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, y = meanavg, color = fishery), size = 1)
b + geom_errorbar(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg, color = fishery), size = 1, span=0.2) + theme_publish() +
  geom_hline(yintercept = 0, linetype = "solid", size = 1.5, color = "darkred") +
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),  # Major grid lines
        panel.grid.minor = element_line(color = "grey", size = 0.2),  # Minor grid lines
        panel.grid.major.x = element_blank(),  # Removing major grid lines on x-axis
        panel.grid.minor.x = element_blank(),  # Removing minor grid lines on x-axis
        panel.background = element_rect(fill = "white", color = "white"),  # White background
        plot.background = element_rect(fill = "white", color = "white"),  # White plot background
        axis.text = element_text(size = 12),  # Adjusting axis text size
        axis.title = element_text(size = 14),  # Adjusting axis titles size
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
        legend.position = "none") + labs(x = "Year (CE)", y = "Propensity for fishery") +
  scale_x_continuous(expand= c(0,0), limits=c(1,1850)) + facet_grid(region~fishery)

###GAM
spatial_gam<-gam(value ~ s(longitude, latitude, bs="tp"),
                 data = plot_anim_expanded %>% filter(
                   fishery == "V1"
                 ), method = "REML")

grid <- expand.grid(
  longitude = seq(min(plot_anim_expanded$longitude), max(plot_anim_expanded$longitude), length.out = 100),
  latitude = seq(min(plot_anim_expanded$latitude), max(plot_anim_expanded$latitude), length.out = 100)
)

grid$prediction <- predict(spatial_gam, newdata = grid, type = "response")

ggplot() +
  geom_tile(data = grid, aes(x = longitude, y = latitude, fill = prediction)) +
  geom_point(data = plot_anim_expanded, aes(x = longitude, y = latitude), color = "black", size = 1) +
  scale_fill_viridis_c(option = "magma", name = "Predicted Value") +
  labs(title = "Spatial GAM Predictions", x = "Longitude", y = "Latitude") +
  theme_minimal()

write.csv(plot_anim_expanded, "../Data/expanded_RS_GAMS.csv")
df_plot<-plot_anim_expanded[,c(4,5,13,9)]

