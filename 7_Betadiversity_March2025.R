##24th Feb 2025 - NSS analysis - ecological diversity changes
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
library(tidyverse)
library(adespatial)
library(ade4)
library(magrittr)
library(vegan)
library(boot)
library(gridExtra)
library(car)
library(purrr)
library(sf)
library(zoo)

source("~/Documents/4-Oceans/3.NSS/2.MANUSCRIPT_Jan2025/Code_final/MANUSCRIPT_functions.R")

#2. setwd & load data ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5a_NSS_SpeciesData_Feb2025_withCEE.xlsx")
df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5b_NSS_SpeciesData_Feb2025_noCEE.xlsx")
df_species<-read_excel("../SupplementaryFiles/SupplementaryTable6b_NSS_FamilyData_Feb2025_noCEE.xlsx")

DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69")
DB2<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69",
      "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
      "#EDC948", "#AF7AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
      "#A0CBE8", "#FFBE7D", "#FF9D9A", "#8CD17D", "#B6992D",
      "#D4A6C8", "#FABFD2", "#79706E", "#D37295", "#C3C3C3"
)
df_species$time_bins<-factor(df_species$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_species$time_bins2<-factor(df_species$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))
NSS_sp<-df_species[!is.na(df_species$GBIF_species),]
NSS_sp<-NSS_sp[!NSS_sp$GBIF_species=="NA",]

#Create time range variable for subsampling chronologies
NSS_sp$Time.Range<-as.numeric(NSS_sp$End_date_CE) - as.numeric(NSS_sp$Start_date_CE)
NSS_sp<-NSS_sp[!NSS_sp$Time.Range=="NA",]
NSS_sp<-NSS_sp %>%
  dplyr::mutate_at('Time.Range', as.numeric)

NSS_sp$count<-"1"
NSS_sp<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                 !NSS_sp$DB_Assemblage_ID == "NA")

length(unique(NSS_sp$DB_Assemblage_ID))
name<-"NSS_sp"
summarise_df(NSS_sp, name)

NSS_sp %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  group_by(Region) %>%
  distinct(GBIF_species) %>%
  tally()

NSS_sp<-NSS_sp %>%
  dplyr::mutate_at('Total_Fish_NISP', as.numeric) 

tmp<-NSS_sp %>%
  dplyr::select("DB_Assemblage_ID","NISP","GBIF_species") %>%
  dplyr::mutate_at('NISP', as.numeric) 

fish_sp_JC<-dcast(NSS_sp %>%
                    dplyr::select("DB_Assemblage_ID","GBIF_species"), DB_Assemblage_ID ~ GBIF_species)
fish_sp_BC<-dcast(tmp, DB_Assemblage_ID ~ GBIF_species, value.var="NISP", sum)
rm(tmp)
fish_sp_JC_sht<-fish_sp_JC[,c(-1)]
fish_sp_BC_sht<-fish_sp_BC[,c(-1)]

set.seed(147)
#species_richness
richness_result <- boot(fish_sp_JC_sht, statistic = species_richness, R = 2000)
#shannon_diversity
shannon_result <- boot(fish_sp_JC_sht, statistic = shannon_diversity, R = 2000)
shannon_index <- diversity(fish_sp_JC_sht, index = "shannon")
rich_index <- specnumber(fish_sp_JC_sht)

# Extract mean richness and confidence intervals
richness_mean <- colMeans(richness_result$t)
richness_ci <- apply(richness_result$t, 2, quantile, probs = c(0.025, 0.975))
shannon_mean <- colMeans(shannon_result$t) 
shannon_ci <- apply(shannon_result$t, 2, quantile, probs = c(0.025, 0.975))
richness_group<-"All_Data"
richness_test<-"Jaccard"
richness_df <- data.frame(
  Site = rownames(fish_sp_JC_sht),
  Mean_Richness = richness_mean,
  Lower_CI = richness_ci[1, ],
  Upper_CI = richness_ci[2, ],
  Group = rep(richness_group, length(richness_ci[2, ])) ,
  Test = rep(richness_test, length(richness_ci[2, ]))
)

total_richness_per_site <- rowSums(fish_sp_JC_sht > 0)  # Count presence of species per site
hist(total_richness_per_site, main = "Local Alpha Diversity per Assemblage", xlab = "Number of Sp.")
richness_df$TotalRich<-total_richness_per_site

richness_df$DB_Assemblage_ID<-fish_sp_JC$DB_Assemblage_ID
richness_df<-richness_df %>% left_join(NSS_sp %>% 
  dplyr::select(DB_Assemblage_ID, Region, Recovery, time_bins, time_bins2,
                Rural_urban_or_neither, Start_date_CE, End_date_CE) %>% distinct())
richness_df$shannon<-shannon_index
richness_df$rich<-rich_index

richness_df %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Mean_Var1 = mean(rich, na.rm = TRUE),
    SD_Var1 = sd(rich, na.rm = TRUE),
    Min_Var1 = min(rich, na.rm = TRUE),
    Max_Var1 = max(rich, na.rm = TRUE),
    Mean_Var2 = mean(shannon, na.rm = TRUE),
    SD_Var2 = sd(shannon, na.rm = TRUE),
    .groups = "drop"
)

richness_df %>%
  dplyr::summarise(
    Mean_Var1 = mean(rich, na.rm = TRUE),
    SD_Var1 = sd(rich, na.rm = TRUE),
    Min_Var1 = min(rich, na.rm = TRUE),
    Max_Var1 = max(rich, na.rm = TRUE),
    Mean_Var2 = mean(shannon, na.rm = TRUE),
    SD_Var2 = sd(shannon, na.rm = TRUE),
    .groups = "drop"
  )

richness_df %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Mean_Var1 = mean(rich, na.rm = TRUE),
    SD_Var1 = sd(rich, na.rm = TRUE),
    Min_Var1 = min(rich, na.rm = TRUE),
    Max_Var1 = max(rich, na.rm = TRUE),
    Mean_Var2 = mean(shannon, na.rm = TRUE),
    SD_Var2 = sd(shannon, na.rm = TRUE),
    .groups = "drop"
  )

richness_df %>%
  dplyr::summarise(
    Mean_Var1 = mean(rich, na.rm = TRUE),
    SD_Var1 = sd(rich, na.rm = TRUE),
    Min_Var1 = min(rich, na.rm = TRUE),
    Max_Var1 = max(rich, na.rm = TRUE),
    Mean_Var2 = mean(shannon, na.rm = TRUE),
    SD_Var2 = sd(shannon, na.rm = TRUE),
    .groups = "drop"
  )

NSS_sp %>%
  group_by(time_bins2) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  group_by(time_bins2) %>%
  distinct(GBIF_family) %>%
  tally()

richness_df %>%
  dplyr::group_by(time_bins2) %>%
  dplyr::summarise(
    Mean_Var1 = mean(rich, na.rm = TRUE),
    SD_Var1 = sd(rich, na.rm = TRUE),
    Min_Var1 = min(rich, na.rm = TRUE),
    Max_Var1 = max(rich, na.rm = TRUE),
    Mean_Var2 = mean(shannon, na.rm = TRUE),
    SD_Var2 = sd(shannon, na.rm = TRUE),
    .groups = "drop"
  )

##Family-level alpha diversity ####
df_species<-read_excel("../SupplementaryFiles/SupplementaryTable6b_NSS_FamilyData_Feb2025_noCEE.xlsx")
df_species$time_bins<-factor(df_species$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_species$time_bins2<-factor(df_species$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))
NSS_sp<-df_species[!is.na(df_species$GBIF_family),]
NSS_sp<-NSS_sp[!NSS_sp$GBIF_family=="NA",]
names(NSS_sp)

#Create time range variable for subsampling chronologies
NSS_sp$Time.Range<-as.numeric(NSS_sp$End_date_CE) - as.numeric(NSS_sp$Start_date_CE)
NSS_sp<-NSS_sp[!NSS_sp$Time.Range=="NA",]
NSS_sp<-NSS_sp %>%
  dplyr::mutate_at('Time.Range', as.numeric)

NSS_sp$count<-"1"
NSS_sp<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                 !NSS_sp$DB_Assemblage_ID == "NA")

length(unique(NSS_sp$DB_Assemblage_ID))
name<-"NSS_sp"
summarise_df(NSS_sp, name)

NSS_sp %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  distinct(GBIF_family) %>%
  tally()

NSS_sp %>%
  group_by(Region) %>%
  distinct(GBIF_family) %>%
  tally()

NSS_sp<-NSS_sp %>%
  dplyr::mutate_at('Total_Fish_NISP', as.numeric) 

tmp<-NSS_sp %>%
  dplyr::select("DB_Assemblage_ID","NISP","GBIF_family") %>%
  dplyr::mutate_at('NISP', as.numeric) 

fish_sp_JC<-dcast(NSS_sp %>%
                    dplyr::select("DB_Assemblage_ID","GBIF_family"), DB_Assemblage_ID ~ GBIF_family)
fish_sp_BC<-dcast(tmp, DB_Assemblage_ID ~ GBIF_family, value.var="NISP", sum)
rm(tmp)
fish_sp_JC_sht<-fish_sp_JC[,c(-1)]
fish_sp_BC_sht<-fish_sp_BC[,c(-1)]

set.seed(147)
#species_richness
richness_result <- boot(fish_sp_JC_sht, statistic = species_richness, R = 2000)
#shannon_diversity
shannon_result <- boot(fish_sp_JC_sht, statistic = shannon_diversity, R = 2000)
shannon_index <- diversity(fish_sp_JC_sht, index = "shannon")
rich_index <- specnumber(fish_sp_JC_sht)

# Extract mean richness and confidence intervals
richness_mean <- colMeans(richness_result$t)
richness_ci <- apply(richness_result$t, 2, quantile, probs = c(0.025, 0.975))
shannon_mean <- colMeans(shannon_result$t) 
shannon_ci <- apply(shannon_result$t, 2, quantile, probs = c(0.025, 0.975))
richness_group<-"Family_Data"
richness_test<-"Jaccard"
richness_df <- data.frame(
  Site = rownames(fish_sp_JC_sht),
  Mean_Richness = richness_mean,
  Lower_CI = richness_ci[1, ],
  Upper_CI = richness_ci[2, ],
  Mean_Shannon = shannon_mean,
  Lower_Shannon_CI = shannon_ci[1, ],
  Upper_Shannon_CI = shannon_ci[2, ],
  Group = rep(richness_group, length(richness_ci[2, ])) ,
  Test = rep(richness_test, length(richness_ci[2, ]))
)

total_richness_per_site <- rowSums(fish_sp_JC_sht > 0)  # Count presence of species per site
hist(total_richness_per_site, main = "Local Alpha Diversity per Assemblage", xlab = "Number of Sp.")
richness_df$TotalRich<-total_richness_per_site

richness_df$DB_Assemblage_ID<-fish_sp_JC$DB_Assemblage_ID
richness_df<-richness_df %>% left_join(NSS_sp %>% 
                                         dplyr::select(DB_Assemblage_ID, Region, Recovery, time_bins, time_bins2,
                                                       Rural_urban_or_neither, Start_date_CE, End_date_CE) %>% distinct())
richness_df$shannon<-shannon_index
richness_df$rich<-rich_index

richness_df %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Mean_Var1 = mean(rich, na.rm = TRUE),
    SD_Var1 = sd(rich, na.rm = TRUE),
    Min_Var1 = min(rich, na.rm = TRUE),
    Max_Var1 = max(rich, na.rm = TRUE),
    Mean_Var2 = mean(shannon, na.rm = TRUE),
    SD_Var2 = sd(shannon, na.rm = TRUE),
    .groups = "drop"
  )

richness_df %>%
  dplyr::summarise(
    Mean_Var1 = mean(rich, na.rm = TRUE),
    SD_Var1 = sd(rich, na.rm = TRUE),
    Min_Var1 = min(rich, na.rm = TRUE),
    Max_Var1 = max(rich, na.rm = TRUE),
    Mean_Var2 = mean(shannon, na.rm = TRUE),
    SD_Var2 = sd(shannon, na.rm = TRUE),
    .groups = "drop"
  )

NSS_sp %>%
  group_by(time_bins2) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  group_by(time_bins2) %>%
  distinct(GBIF_family) %>%
  tally()

richness_df %>%
  dplyr::group_by(time_bins2) %>%
  dplyr::summarise(
    Mean_Var1 = mean(rich, na.rm = TRUE),
    SD_Var1 = sd(rich, na.rm = TRUE),
    Min_Var1 = min(rich, na.rm = TRUE),
    Max_Var1 = max(rich, na.rm = TRUE),
    Mean_Var2 = mean(shannon, na.rm = TRUE),
    SD_Var2 = sd(shannon, na.rm = TRUE),
    .groups = "drop"
  )

# Rarefaction - all data ####
rarefaction_results <- specaccum(fish_sp_JC_sht, method = "rarefaction", permutations = 100)
rarefaction_data <- data.frame(
  Sites = rarefaction_results$sites,
  Richness = rarefaction_results$richness,
  Upper_CI = rarefaction_results$richness + rarefaction_results$sd,  # Upper bound of CI
  Lower_CI = rarefaction_results$richness - rarefaction_results$sd   # Lower bound of CI
)
ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = "purple", size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "purple", alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve with Confidence Intervals: Species Richness vs. Number of Sites",
       x = "Number of Sites Sampled",
       y = "Species Richness") +
  theme(axis.text = element_text(size = 12))

richness_df$sites<-seq(1:length(richness_df$Mean_Richness))
#Check whether curve has reached asymptote
model <- nls(Mean_Richness ~ a * sites / (b + sites), 
             data = richness_df, 
             start = list(a = max(richness_df$Mean_Richness), b = 1))
summary(model)
fitted_values <- predict(model, newdata = data.frame(sites = richness_df$Sites))
plot_data <- data.frame(sites = richness_df$sites, 
                        richness = richness_df$Mean_Richness, 
                        fitted_richness = fitted_values)

fixed_b_model <- nls(Mean_Richness ~ a * sites / (1 + sites), 
                     data = richness_df, 
                     start = list(a = max(richness_df$Mean_Richness)))

rarefaction_plot <- ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = DB[2], size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill =  DB[2], alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve",
       x = "Number of Sites Sampled",
       y = "Species Richness",
       tag = "A)") +
  theme(axis.text = element_text(size = 12)) + theme_publish()

model_plot <- ggplot(plot_data, aes(x = sites, y = richness)) +
  geom_point(color =  DB[3]) +
  geom_line(aes(y = fitted_richness), color = DB[2], linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Michaelis-Menten Model Fit", x = "Sites", y = "Species Richness", tag="B)") +
  theme(axis.text = element_text(size = 12), plot.title = element_text(size = 14)) + theme_publish()

png(filename="../SupplementaryFiles/Rarefaction_alldata.png", 
    width=9.5, height=4.5, units="in",res=600)
grid.arrange(rarefaction_plot, model_plot, ncol = 2)
dev.off()

####Calculate moving avergae and plot changes overtime - species ####
richness_df$Start_date_CE<-as.numeric(richness_df$Start_date_CE)
richness_df$End_date_CE<-as.numeric(richness_df$End_date_CE)
richness_expanded<- richness_df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(year = list(seq(Start_date_CE, End_date_CE, by=25))) %>%
  unnest(year) %>%
  ungroup()
names(richness_expanded)

df_moving_avg<-richness_expanded %>%
  arrange(year) %>%
  mutate(moving_avg_rich = rollapply(rich, width=3, FUN=mean, fill=NA, align="center"),
         moving_avg_shannon = rollapply(shannon, width=3, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 25 == 0) %>%
  ungroup()

df_moving_avg<-df_moving_avg[complete.cases(df_moving_avg),]
ggplot(df_moving_avg,
       aes(year, moving_avg_rich)) + 
  geom_point(alpha=0.1, shape=19, size=0.95) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) 

ggplot(df_moving_avg,
       aes(year, moving_avg_shannon)) + 
  geom_point(alpha=0.1, shape=19, size=0.95) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) 

df_moving_avg2 <- richness_expanded %>%
  dplyr::select(year, rich, shannon) %>% 
  arrange(year) %>%
  mutate(moving_avg = rollapply(rich, width=100, FUN=mean, fill=NA, align="center"),
         moving_avg_shannon = rollapply(shannon, width=100, FUN=mean, fill=NA, align="center")) %>%
  #  filter(year %% 100 == 0) %>%
  ungroup()

ggplot(df_moving_avg,
       aes(year, moving_avg_rich)) + 
  geom_point(alpha=0.1, shape=19, size=0.95) + 
  geom_point(alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values=DB) +
  scale_color_manual(values=DB) + geom_smooth( method = "loess",
                                               data=df_moving_avg2, aes(x = year, y = moving_avg), size = 1, span=0.2
  )

#Plot species richness overtime ####
df_moving_avg3 <- richness_expanded %>%
  dplyr::select(rich, year) %>% 
  arrange(year) %>%
  mutate(moving_avg = rollapply(rich, width=20, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg3$year<-as.character(df_moving_avg3$year)
df_moving_avg3<-df_moving_avg3 %>% 
  dplyr::group_by(year) %>%
  dplyr::summarise(
    meanavg = mean(moving_avg, na.rm = TRUE),
    sdavg = sd(moving_avg, na.rm = TRUE)
  )
df_moving_avg3$year<-as.numeric(df_moving_avg3$year)
a<-ggplot() + 
  geom_point(data = df_moving_avg,
             aes(x=year, y=moving_avg_rich),fill="#76c0c1", color="#76c0c1", alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avg,
             aes(x=year, y=moving_avg_rich), alpha=0.5, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Species Richness over time", x = "Time Period", y = "Richness") +
  geom_point(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
                                                  aes(x = year, y = meanavg),color="#76c0c1", size = 1)

png(filename="../SupplementaryFiles/Richness_trends_NSSDB.png", 
    width=6.5, height=4.5, units="in",res=600)
a + geom_errorbar(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg),color = "#76c0c1", size = 1, span=0.2) + theme_publish() + theme(legend.position="none")
dev.off()

#Plot shannon indices overtime
df_moving_avg3 <- richness_expanded %>%
  dplyr::select(shannon, year) %>% 
  arrange(year) %>%
  mutate(moving_avg = rollapply(shannon, width=20, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg3$year<-as.character(df_moving_avg3$year)
df_moving_avg3<-df_moving_avg3 %>% 
  dplyr::group_by(year) %>%
  dplyr::summarise(
    meanavg = mean(moving_avg, na.rm = TRUE),
    sdavg = sd(moving_avg, na.rm = TRUE)
  )
df_moving_avg3$year<-as.numeric(df_moving_avg3$year)
a<-ggplot() + 
  geom_point(data = df_moving_avg,
             aes(x=year, y=moving_avg_shannon, fill=DB[4], color=DB[4]), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avg,
             aes(x=year, y=moving_avg_shannon), alpha=0.5, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(title = "Shannon index over time", x = "Time Period", y = "Shannon index") +
  geom_point(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
             aes(x = year, y = meanavg, color=DB[4]), size = 1)

png(filename="../SupplementaryFiles/Shannon_trends_NSSDB.png", 
    width=6.5, height=4.5, units="in",res=600)
a + geom_errorbar(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg, color = DB[4]), size = 1, span=0.2) + theme_publish() + theme(legend.position="none")
dev.off()


#plot_anim_expanded$year<-as.numeric(plot_anim_expanded$year)
df_moving_avg3 <- richness_expanded %>%
  dplyr::select(rich, Region, year) %>% 
  group_by(Region) %>%
  arrange(year) %>%
  mutate(moving_avg = rollapply(rich, width=20, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg3$year<-as.character(df_moving_avg3$year)
df_moving_avg3<-df_moving_avg3 %>% 
  dplyr::group_by(Region, year) %>%
  dplyr::summarise(
    meanavg = mean(moving_avg, na.rm = TRUE),
    sdavg = sd(moving_avg, na.rm = TRUE)
  )
df_moving_avg3$year<-as.numeric(df_moving_avg3$year)
a<-ggplot() + 
  geom_point(data = df_moving_avg,
             aes(x=year, y=moving_avg_rich, fill=Region, color=Region), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avg,
             aes(x=year, y=moving_avg_rich), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~Region) +
  labs(title = "Regional Species Richness over time", x = "Time Period", y = "Richness") +
  scale_fill_manual(values=DB[2:6]) +
  scale_color_manual(values=DB[2:6]) + geom_point(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
                                             aes(x = year, y = meanavg, color = Region), size = 1)

png(filename="../SupplementaryFiles/Richness_trendsperregion.png", 
    width=9, height=4.5, units="in",res=600)
a + geom_errorbar(data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg3[!is.na(df_moving_avg3$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg, color = Region), size = 1, span=0.2) + theme_publish()
dev.off()

###Shannon ####
df_moving_avg4 <- richness_expanded %>%
  dplyr::select(Mean_Shannon, Region, year) %>% 
  group_by(Region) %>%
  arrange(year) %>%
  mutate(moving_avg = rollapply(Mean_Shannon, width=20, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 100 == 0) %>%
  ungroup()
df_moving_avg4$year<-as.character(df_moving_avg4$year)
df_moving_avg4<-df_moving_avg4 %>% 
  dplyr::group_by(Region, year) %>%
  dplyr::summarise(
    meanavg = mean(moving_avg, na.rm = TRUE),
    sdavg = sd(moving_avg, na.rm = TRUE)
  )
df_moving_avg4$year<-as.numeric(df_moving_avg4$year)
b<-ggplot() + 
  geom_point(data = df_moving_avg,
             aes(x=year, y=moving_avg_shannon, fill=Region, color=Region), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avg,
             aes(x=year, y=moving_avg_shannon), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~Region) +
  labs(title = "Regional Shannon Indices over time", x = "Time Period", y = "Shannon Index") +
  scale_fill_manual(values=DB[2:6]) +
  scale_color_manual(values=DB[2:6]) + geom_point(data=df_moving_avg4[!is.na(df_moving_avg4$sdavg),], 
                                                  aes(x = year, y = meanavg, color = Region), size = 1)

png(filename="../SupplementaryFiles/Shannon_trendsperregion.png", 
    width=9, height=4.5, units="in",res=600)
b + geom_errorbar(data=df_moving_avg4[!is.na(df_moving_avg4$sdavg),], aes(x = year, ymin = meanavg - sdavg, 
                                                                          ymax = meanavg + sdavg), width=1.5, color="grey15") + 
  geom_smooth( method = "loess",
               data=df_moving_avg4[!is.na(df_moving_avg4$sdavg),], 
               se = TRUE, span = 0.2,
               aes(x = year, y = meanavg, color = Region), size = 1, span=0.2) + theme_publish()
dev.off()

###Plot Global Shannon/Richness (need to melt richness_expanded first and then use this df throughout)
str(richness_expanded)
richness_expanded_melt<-reshape2::melt(richness_expanded[,c(11,2,5,12,13,1,15,16,19)], id.vars = c("DB_Assemblage_ID","Region",
                                                                                                   "Recovery","time_bins2",
                                                                                                   "Rural_urban_or_neither","year"))

richness_expanded_melt<-reshape2::melt(richness_df %>% dplyr::select(
  shannon, rich, Region, Recovery, time_bins2, DB_Assemblage_ID, Rural_urban_or_neither), id.vars = c("DB_Assemblage_ID","Region",
                                                                                                   "Recovery","time_bins2",
                                                                                                   "Rural_urban_or_neither")) %>% distinct()
richness_expanded_melt<-recode()

df_moving_avg<-richness_expanded %>%
  arrange(year) %>%
  mutate(moving_avg_rich = rollapply(Mean_Richness, width=3, FUN=mean, fill=NA, align="center"),
         moving_avg_shannon = rollapply(Mean_Shannon, width=3, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 25 == 0) %>%
  ungroup()

df_moving_avg<-df_moving_avg[complete.cases(df_moving_avg),]
df_moving_avg$year<-as.character(df_moving_avg$year)
df_moving_avg_global<-df_moving_avg %>%
  dplyr::select(Mean_Richness, time_bins2, year, Recovery) %>% reshape2::melt()
df_moving_avg_global$year<-as.numeric(df_moving_avg_global$year)
df_moving_avg_global2 <- df_moving_avg_global %>%
  dplyr::select(year, variable, value, Recovery) %>% 
  arrange(year) %>%
  mutate(moving_avg = rollapply(value, width=100, FUN=mean, fill=NA, align="center")) %>%
  ungroup()

df_moving_avg_global3 <- df_moving_avg_global %>%
  dplyr::select(year, variable, value, Recovery) %>% 
  group_by(variable) %>%
  arrange(year) %>%
  mutate(moving_avg = rollapply(value, width=20, FUN=mean, fill=NA, align="center")) %>%
  filter(year %% 100 == 0) %>%
  ungroup()

df_moving_avg_global3$year<-as.character(df_moving_avg_global3$year)

df_moving_avg_global3<-df_moving_avg_global3 %>% 
  dplyr::group_by(variable, year, Recovery) %>%
  dplyr::summarise(
    meanavg = mean(moving_avg, na.rm = TRUE),
    sdavg = sd(moving_avg, na.rm = TRUE)
  )
df_moving_avg_global3$year<-as.numeric(df_moving_avg_global3$year)

ggplot() + 
  geom_point(data = df_moving_avg_global,
             aes(x=year, y=value, fill=Recovery, color=Recovery), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avg_global,
             aes(x=year, y=value), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(Recovery~variable, scales="free") +
  labs(title = "Alpha Diversity over time", x = "Time Period") +
  scale_fill_manual(values=DB[2:7]) +
  scale_color_manual(values=DB[2:7]) + geom_point(data=df_moving_avg_global3[!is.na(df_moving_avg_global3$sdavg),], 
                                                  aes(x = year, y = meanavg, color = variable), size = 1)

str(richness_expanded_melt)
ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
             aes(x=time_bins2, y=value, fill=Recovery, color=Recovery), alpha=0.1, shape=19, size=0.95) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(Region~variable, scales="free") +
  labs(title = "Alpha Diversity over time", x = "Time Period") +
  scale_fill_manual(values=DB[2:7]) +
  scale_color_manual(values=DB[2:7])

ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
               aes(x=time_bins2, y=value, fill=Region, color=Region)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_grid(~variable, scales="free") +
  labs(title = "Alpha Diversity over time", x = "Time Period") +
  scale_fill_manual(values=DB[2:7]) +
  scale_color_manual(values=DB[2:7])

ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
               aes(x=time_bins2, y=value, fill=Region, color=Region)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~variable, scales="free") +
  labs(title = "Alpha Diversity over time", x = "Time Period") +
  scale_fill_manual(values=DB[2:7]) +
  scale_color_manual(values=DB[2:7])

ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
               aes(x=time_bins2, y=value, fill=Recovery, color=Recovery)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~variable, scales="free") +
  labs(title = "Alpha Diversity over time", x = "Time Period") +
  scale_fill_manual(values=DB[2:7]) +
  scale_color_manual(values=DB[2:7])

ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
               aes(x=Region, y=value, fill=Recovery, color=Recovery)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~variable, scales="free") +
  labs(title = "Alpha Diversity by recovery type", x = "Region") +
  scale_fill_manual(values=DB[2:7]) +
  scale_color_manual(values=DB[2:7]) + theme_publish()

DB3<-c("#3F007D","#FEC287FF","#807DBA","#DADAEB")
png(filename="../SupplementaryFiles/DiversitybyRecoveryType.png", 
    width=7.2, height=3.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
               aes(x=Region, y=value, fill=Recovery, color=Recovery)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~variable2, scales="free") +
  labs(title = "Alpha Diversity by recovery type", x = "Region") +
  scale_fill_manual(values=DB3) +
  scale_color_manual(values=DB3) + theme_publish()
dev.off()

KW1<-kruskal.test(value ~ Region, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="rich"))

KW2<-kruskal.test(value ~ time_bins2, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="rich"))

KW3<-kruskal.test(value ~ Region, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="shannon"))

KW4<-kruskal.test(value ~ time_bins2, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="shannon"))

KW5<-kruskal.test(value ~ Recovery, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="rich"))

KW6<-kruskal.test(value ~ Recovery, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="shannon"))

KW7<-kruskal.test(value ~ Rural_urban_or_neither, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="rich"))

KW8<-kruskal.test(value ~ Rural_urban_or_neither, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="shannon"))

richness_expanded_melt$groups<-paste(richness_expanded_melt$Region, richness_expanded_melt$time_bins2)

KW9<-kruskal.test(value ~ groups, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="rich"))

KW10<-kruskal.test(value ~ groups, data = richness_expanded_melt %>% 
                    dplyr::filter(variable=="shannon"))

library(forcats)
richness_expanded_melt$variable2<-recode_factor(richness_expanded_melt$variable,
                                                levels = c(
                                         "rich" = 'Species Richness',
                                         "shannon" = 'Shannon Index'))
richness_expanded_melt$variable2<-gsub("shannon","Shannon Index",richness_expanded_melt$variable2)
richness_expanded_melt$variable2<-gsub("rich", "Species Richness",richness_expanded_melt$variable2)
ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
               aes(x=time_bins2, y=value, fill=Recovery, color=Recovery)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~variable, scales="free") +
  labs(title = "Alpha Diversity by recovery type", x = "Region") +
  scale_fill_manual(values=DB2) +
  scale_color_manual(values=DB2) + theme_publish()

DB4<-c("#F8A31B","#AA3929","#8E9CA3")
png(filename="../SupplementaryFiles/DiversitybySiteType.png", 
    width=7.2, height=3.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
               aes(x=Region, y=value, fill=Rural_urban_or_neither, color=Rural_urban_or_neither)) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~variable2, scales="free") +
  labs(title = "Alpha Diversity by site type", x = "Region") +
  scale_fill_manual(values=DB4) +
  scale_color_manual(values=DB4) + theme_publish()
dev.off()

ggplot() + 
  geom_point(data = df_moving_avg_global %>%
               dplyr::filter(variable=="Mean_Richness"),
             aes(x=year, y=value, fill=Recovery, color=Recovery), alpha=0.1, shape=19, size=0.95) + 
  geom_point(data = df_moving_avg_global,
             aes(x=year, y=value), alpha=0.1, shape=1, size=0.95, color="grey50") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~variable, scales="free_y") +
  labs(title = "Alpha Diversity over time", x = "Time Period") +
  scale_fill_manual(values=DB[2:6]) +
  scale_color_manual(values=DB[2:6]) + geom_point(data=df_moving_avg_global3[!is.na(df_moving_avg_global3$sdavg),], 
                                                  aes(x = year, y = meanavg, color = variable), size = 1)

png(filename="../SupplementaryFiles/DiversitybyTimeBins.png", 
    width=5.5, height=7.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_expanded_melt,
               aes(x=time_bins2, y=value, fill=time_bins2, color=time_bins2)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~variable2, scales="free", ncol=1) +
  labs(title = "Alpha Diversity by site type", x = "Region") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off()

#LCBD ####
richness_df$LCBD<-"Beta Diversity - nestedness"
richness_df$LCBD2<-"Beta Diversity - replacement"
names(richness_df)
str(richness_df)

#LCBD - Jaccard replacement ####
png(filename="../SupplementaryFiles/LCBDbyTimeBins.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=time_bins2, y=local.repl_J, fill=time_bins2, color=time_bins2)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~LCBD2, scales="free", ncol=1) +
  labs(title = "", x = "Region", y = "Replacement (Jaccard)") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off()

png(filename="../SupplementaryFiles/LCBDbyRecovery.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=Region, y=local.rich_J, fill=Recovery, color=Recovery)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~LCBD, scales="free", ncol=1) +
  labs(title = "", x = "Region", y = "Richness (Jaccard)") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off

png(filename="../SupplementaryFiles/LCBDbyRuralUrban.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=Region, y=local.repl_J, fill=Rural_urban_or_neither, color=Rural_urban_or_neither)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~LCBD2, scales="free", ncol=1) +
  labs(title = "", x = "Region", y = "Replacement (Jaccard)") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off

#Plot beta diversity nestedness by groups - Jaccard ####
png(filename="../SupplementaryFiles/LCBDbyTimeBins.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=time_bins2, y=local.rich_J, fill=time_bins2, color=time_bins2)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~local.repl_BC, scales="free", ncol=1) +
  labs(title = "Alpha Diversity by site type", x = "Region") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off()

png(filename="../SupplementaryFiles/LCBDbyRecovery.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=Region, y=local.repl_BC, fill=Recovery, color=Recovery)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~local.repl_BC, scales="free", ncol=1) +
  labs(title = "Alpha Diversity by site type", x = "Region") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off

png(filename="../SupplementaryFiles/LCBDbyRuralUrban.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=Region, y=local.repl_BC, fill=Rural_urban_or_neither, color=Rural_urban_or_neither)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~local.repl_BC, scales="free", ncol=1) +
  labs(title = "Alpha Diversity by site type", x = "Region") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off

#Plot beta diversity replacement by groups - bray-curtis ####
png(filename="../SupplementaryFiles/LCBDbyTimeBins.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=time_bins2, y=local.repl_BC, fill=time_bins2, color=time_bins2)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~local.repl_BC, scales="free", ncol=1) +
  labs(title = "Alpha Diversity by site type", x = "Region") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off()

png(filename="../SupplementaryFiles/LCBDbyRecovery.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=Region, y=local.repl_BC, fill=Recovery, color=Recovery)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~local.repl_BC, scales="free", ncol=1) +
  labs(title = "Alpha Diversity by site type", x = "Region") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off

png(filename="../SupplementaryFiles/LCBDbyRuralUrban.png", 
    width=5.5, height=4.5, units="in",res=600)
ggplot() + 
  geom_boxplot(data = richness_df,
               aes(x=Region, y=local.repl_BC, fill=Rural_urban_or_neither, color=Rural_urban_or_neither)) + 
  geom_jitter() +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + facet_wrap(~local.repl_BC, scales="free", ncol=1) +
  labs(title = "Alpha Diversity by site type", x = "Region") +
  scale_fill_manual(values=DB2[2:10]) +
  scale_color_manual(values=DB2[2:10]) + theme_publish()
dev.off

###Rarefaction - per region ####
#1. Britain & Ireland
NSS_sp_BI<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                 !NSS_sp$DB_Assemblage_ID == "NA",
                                 NSS_sp$Region=="Britian & Ireland")

###All Data:Create species matrix (long format, each row = assemblage; each column = species), presence and absence for jaccard and Bray-Curtis
tmp<-NSS_sp_BI %>%
  dplyr::select("DB_Assemblage_ID","NISP","GBIF_species") %>%
  dplyr::mutate_at('NISP', as.numeric) 
fish_sp_JC<-dcast(NSS_sp_BI %>%
                    dplyr::select("DB_Assemblage_ID","GBIF_species"), DB_Assemblage_ID ~ GBIF_species)
fish_sp_BC<-dcast(tmp, DB_Assemblage_ID ~ GBIF_species, value.var="NISP", sum)
rm(tmp)
fish_sp_JC_sht<-fish_sp_JC[,c(-1)]
fish_sp_BC_sht<-fish_sp_BC[,c(-1)]

set.seed(12)
boot_result <- boot(fish_sp_JC_sht, statistic = species_richness, R = 1000)
# Extract mean richness and confidence intervals
richness_mean <- colMeans(boot_result$t)  # Mean richness across bootstraps
richness_ci <- apply(boot_result$t, 2, quantile, probs = c(0.025, 0.975))  
richness_df <- data.frame(
  Site = rownames(fish_sp_JC_sht),
  Mean_Richness = richness_mean,
  Lower_CI = richness_ci[1, ],
  Upper_CI = richness_ci[2, ]
)

total_richness_per_site <- rowSums(fish_sp_JC_sht > 0)  # Count presence of species per site

# Create a combined dataframe for plotting mean and total richness
richness_comparison_df <- data.frame(
  Richness_Type = rep(c("Mean Richness", "Total Richness"), each = length(richness_df$Mean_Richness)),
  Richness = c(richness_df$Mean_Richness, total_richness_per_site)
)

# Violin plot of both mean and total richness
ggplot(richness_comparison_df, aes(x = Richness_Type, y = Richness, fill = Richness_Type)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Violin Plot: Mean and Total Species Richness Across All Sites",
       x = "Richness Type",
       y = "Species Richness") +
  scale_fill_manual(values = c("blue", "green")) +
  theme(axis.text.x = element_text(size = 12))

# Rarefaction curve (sampling curve) for species richness
rarefaction_results <- specaccum(fish_sp_JC_sht, method = "rarefaction", permutations = 100)
rarefaction_data <- data.frame(
  Sites = rarefaction_results$sites,
  Richness = rarefaction_results$richness,
  Upper_CI = rarefaction_results$richness + rarefaction_results$sd,  # Upper bound of CI
  Lower_CI = rarefaction_results$richness - rarefaction_results$sd   # Lower bound of CI
)

ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = "purple", size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "purple", alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve with Confidence Intervals: Species Richness vs. Number of Sites",
       x = "Number of Sites Sampled",
       y = "Species Richness") +
  theme(axis.text = element_text(size = 12))

richness_df$sites<-seq(1:length(richness_df$Mean_Richness))
#Check whether curve has reached asymptote
richness_df$sites<-seq(1:1509)
model <- nls(rich ~ a * sites / (b + sites), 
             data = richness_df, 
             start = list(a = 38, b = 1))
summary(model)
fitted_values <- predict(model, newdata = data.frame(sites = richness_df$sites))
plot_data <- data.frame(sites = richness_df$sites, 
                        richness = richness_df$rich, 
                        fitted_richness = fitted_values)

rarefaction_plot <- ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = DB[4], size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill =  DB[4], alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve",
       x = "Number of Sites Sampled",
       y = "Species Richness",
       tag = "A)") +
  theme(axis.text = element_text(size = 12)) + theme_publish()

model_plot <- ggplot(plot_data, aes(x = sites, y = richness)) +
  geom_point(color =  DB[3], alpha=0.7) +
  geom_line(aes(y = fitted_richness), color = DB[4], linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Michaelis-Menten Model Fit", x = "Sites", y = "Species Richness", tag="B)") +
  theme(axis.text = element_text(size = 12), plot.title = element_text(size = 14)) + theme_publish()

#chao2_result
specpool(fish_sp_JC_sht)

png(filename="../SupplementaryFiles/Rarefaction_NSSDB.png", 
    width=9.5, height=4.5, units="in",res=600)
rarefaction_plot
dev.off()

#Jaccard Index - species level #
fish.bd.presence<-beta.div.comp(fish_sp_JC_sht,coef="J",quant=F)
fish.bd.presence$part

#Quantitative form of the sorensen index
fish.bd.abundance<-beta.div.comp(fish_sp_BC_sht,coef="S", quant=T)
fish.bd.abundance$part

#Local contribution of each site to BD - replacement
local.repl_BC<-LCBD.comp(fish.bd.abundance$repl, sqrt.D = T)
richness_df$local.repl_BC<-local.repl_BC$LCBD
hist(local.repl_BC$LCBD)
#Almost normal distribution of replacement across all sites, suggesting no particular sites are driving the replacement

#Local contribution of each site to BD - richness
local.rich_BC<-LCBD.comp(fish.bd.abundance$rich, sqrt.D = T)
richness_df$local.rich_BC<-local.rich_BC$local.rich_BC
hist(local.repl_BC$LCBD)

#Local contribution of each site to BD - replacement
local.repl_J<-LCBD.comp(fish.bd.presence$repl, sqrt.D = T)
richness_df$local.repl_J<-local.repl_J$LCBD
hist(local.repl_J$LCBD)
#Almost normal distribution of replacement across all sites, suggesting no particular sites are driving the replacement

#Local contribution of each site to BD - richness
local.rich_J<-LCBD.comp(fish.bd.presence$rich, sqrt.D = T)
richness_df$local.rich_J<-local.rich_J$LCBD
hist(local.repl_BC$LCBD)

#Local contribution of each site to BD - replacement
local.repl<-LCBD.comp(fish.bd.abundance$repl, sqrt.D = T)
local.repl
hist(local.repl$LCBD)
#Almost normal distribution of replacement across all sites, suggesting no particular sites are driving the replacement

#Local contribution of each site to BD - richness
local.rich<-LCBD.comp(fish.bd.abundance$rich, sqrt.D = T)
local.rich
hist(local.rich$LCBD)
#Skewed distribution of richness, suggesting some sites have higher than average richness than others

#2. Scandinavia 
NSS_sp_Sc<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                    !NSS_sp$DB_Assemblage_ID == "NA",
                                    NSS_sp$Region=="Scandinavia")

tmp<-NSS_sp_Sc %>%
  dplyr::select("DB_Assemblage_ID","NISP","GBIF_species") %>%
  dplyr::mutate_at('NISP', as.numeric) 
fish_sp_JC<-dcast(NSS_sp_Sc %>%
                    dplyr::select("DB_Assemblage_ID","GBIF_species"), DB_Assemblage_ID ~ GBIF_species)
fish_sp_BC<-dcast(tmp, DB_Assemblage_ID ~ GBIF_species, value.var="NISP", sum)
rm(tmp)
fish_sp_JC_sht<-fish_sp_JC[,c(-1)]
fish_sp_BC_sht<-fish_sp_BC[,c(-1)]

set.seed(12)
boot_result <- boot(fish_sp_JC_sht, statistic = species_richness, R = 1000)
# Extract mean richness and confidence intervals
richness_mean <- colMeans(boot_result$t)  # Mean richness across bootstraps
richness_ci <- apply(boot_result$t, 2, quantile, probs = c(0.025, 0.975))  
richness_df <- data.frame(
  Site = rownames(fish_sp_JC_sht),
  Mean_Richness = richness_mean,
  Lower_CI = richness_ci[1, ],
  Upper_CI = richness_ci[2, ]
)

total_richness_per_site <- rowSums(fish_sp_JC_sht > 0)  # Count presence of species per site

# Create a combined dataframe for plotting mean and total richness
richness_comparison_df <- data.frame(
  Richness_Type = rep(c("Mean Richness", "Total Richness"), each = length(richness_df$Mean_Richness)),
  Richness = c(richness_df$Mean_Richness, total_richness_per_site)
)

# Violin plot of both mean and total richness
ggplot(richness_comparison_df, aes(x = Richness_Type, y = Richness, fill = Richness_Type)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Violin Plot: Mean and Total Species Richness Across All Sites",
       x = "Richness Type",
       y = "Species Richness") +
  scale_fill_manual(values = c("blue", "green")) +
  theme(axis.text.x = element_text(size = 12))

# Rarefaction curve (sampling curve) for species richness
rarefaction_results <- specaccum(fish_sp_JC_sht, method = "rarefaction", permutations = 100)
rarefaction_data <- data.frame(
  Sites = rarefaction_results$sites,
  Richness = rarefaction_results$richness,
  Upper_CI = rarefaction_results$richness + rarefaction_results$sd,  # Upper bound of CI
  Lower_CI = rarefaction_results$richness - rarefaction_results$sd   # Lower bound of CI
)
ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = "purple", size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "purple", alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve with Confidence Intervals: Species Richness vs. Number of Sites",
       x = "Number of Sites Sampled",
       y = "Species Richness") +
  theme(axis.text = element_text(size = 12))

richness_df$sites<-seq(1:length(richness_df$Mean_Richness))
#Check whether curve has reached asymptote
model <- nls(Mean_Richness ~ a * sites / (b + sites), 
             data = richness_df, 
             start = list(a = max(richness_df$Mean_Richness), b = 1))
summary(model)
fitted_values <- predict(model, newdata = data.frame(sites = richness_df$Sites))
plot_data <- data.frame(sites = richness_df$sites, 
                        richness = richness_df$Mean_Richness, 
                        fitted_richness = fitted_values)

fixed_b_model <- nls(Mean_Richness ~ a * sites / (1 + sites), 
                     data = richness_df, 
                     start = list(a = max(richness_df$Mean_Richness)))

rarefaction_plot <- ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = DB[1], size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill =  DB[1], alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve",
       x = "Number of Sites Sampled",
       y = "Species Richness",
       tag = "A)") +
  theme(axis.text = element_text(size = 12)) + theme_publish()

model_plot <- ggplot(plot_data, aes(x = sites, y = richness)) +
  geom_point(color =  DB[3], alpha=0.7) +
  geom_line(aes(y = fitted_richness), color = DB[1], linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Michaelis-Menten Model Fit", x = "Sites", y = "Species Richness", tag="B)") +
  theme(axis.text = element_text(size = 12), plot.title = element_text(size = 14)) + theme_publish()

png(filename="../SupplementaryFiles/Rarefaction_Scandinavia.png", 
    width=9.5, height=4.5, units="in",res=600)
grid.arrange(rarefaction_plot, model_plot, ncol = 2)
dev.off()

#Jaccard Index #
fish.bd.presence<-beta.div.comp(fish_sp_JC_sht,coef="J",quant=F)
fish.bd.presence$part

#Quantitative form of the sorensen index
fish.bd.abundance<-beta.div.comp(fish_sp_BC_sht,coef="J",quant=T)
fish.bd.abundance$part

#Local contribution of each site to BD - replacement
local.repl<-LCBD.comp(fish.bd.abundance$repl, sqrt.D = T)
local.repl
hist(local.repl$LCBD)
#Almost normal distribution of replacement across all sites, suggesting no particular sites are driving the replacement

#Local contribution of each site to BD - richness
local.rich<-LCBD.comp(fish.bd.abundance$rich, sqrt.D = T)
local.rich
hist(local.rich$LCBD)

#3. Western Europe 
NSS_sp_We<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                    !NSS_sp$DB_Assemblage_ID == "NA",
                                    NSS_sp$Region=="Western Europe")

###All Data:Create species matrix (long format, each row = assemblage; each column = species), presence and absence for jaccard and Bray-Curtis ####
tmp<-NSS_sp_We %>%
  dplyr::select("DB_Assemblage_ID","NISP","GBIF_species") %>%
  dplyr::mutate_at('NISP', as.numeric) 
fish_sp_JC<-dcast(NSS_sp_We %>%
                    dplyr::select("DB_Assemblage_ID","GBIF_species"), DB_Assemblage_ID ~ GBIF_species)
fish_sp_BC<-dcast(tmp, DB_Assemblage_ID ~ GBIF_species, value.var="NISP", sum)
rm(tmp)
fish_sp_JC_sht<-fish_sp_JC[,c(-1)]
fish_sp_BC_sht<-fish_sp_BC[,c(-1)]

set.seed(12)
boot_result <- boot(fish_sp_JC_sht, statistic = species_richness, R = 1000)
# Extract mean richness and confidence intervals
richness_mean <- colMeans(boot_result$t)  # Mean richness across bootstraps
richness_ci <- apply(boot_result$t, 2, quantile, probs = c(0.025, 0.975))  
richness_df <- data.frame(
  Site = rownames(fish_sp_JC_sht),
  Mean_Richness = richness_mean,
  Lower_CI = richness_ci[1, ],
  Upper_CI = richness_ci[2, ]
)

total_richness_per_site <- rowSums(fish_sp_JC_sht > 0)  # Count presence of species per site

# Create a combined dataframe for plotting mean and total richness
richness_comparison_df <- data.frame(
  Richness_Type = rep(c("Mean Richness", "Total Richness"), each = length(richness_df$Mean_Richness)),
  Richness = c(richness_df$Mean_Richness, total_richness_per_site)
)

# Violin plot of both mean and total richness
ggplot(richness_comparison_df, aes(x = Richness_Type, y = Richness, fill = Richness_Type)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Violin Plot: Mean and Total Species Richness Across All Sites",
       x = "Richness Type",
       y = "Species Richness") +
  scale_fill_manual(values = c("blue", "green")) +
  theme(axis.text.x = element_text(size = 12))

# Rarefaction curve (sampling curve) for species richness
rarefaction_results <- specaccum(fish_sp_JC_sht, method = "rarefaction", permutations = 100)
rarefaction_data <- data.frame(
  Sites = rarefaction_results$sites,
  Richness = rarefaction_results$richness,
  Upper_CI = rarefaction_results$richness + rarefaction_results$sd,  # Upper bound of CI
  Lower_CI = rarefaction_results$richness - rarefaction_results$sd   # Lower bound of CI
)
ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = "purple", size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "purple", alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve with Confidence Intervals: Species Richness vs. Number of Sites",
       x = "Number of Sites Sampled",
       y = "Species Richness") +
  theme(axis.text = element_text(size = 12))

richness_df$sites<-seq(1:length(richness_df$Mean_Richness))
#Check whether curve has reached asymptote
model <- nls(Mean_Richness ~ a * sites / (b + sites), 
             data = richness_df, 
             start = list(a = max(richness_df$Mean_Richness), b = 1))
summary(model)
fitted_values <- predict(model, newdata = data.frame(sites = richness_df$Sites))
plot_data <- data.frame(sites = richness_df$sites, 
                        richness = richness_df$Mean_Richness, 
                        fitted_richness = fitted_values)

fixed_b_model <- nls(Mean_Richness ~ a * sites / (1 + sites), 
                     data = richness_df, 
                     start = list(a = max(richness_df$Mean_Richness)))

rarefaction_plot <- ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = DB[5], size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill =  DB[5], alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve",
       x = "Number of Sites Sampled",
       y = "Species Richness",
       tag = "A)") +
  theme(axis.text = element_text(size = 12)) + theme_publish()

model_plot <- ggplot(plot_data, aes(x = sites, y = richness)) +
  geom_point(color =  DB[3], alpha=0.7) +
  geom_line(aes(y = fitted_richness), color = DB[5], linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Michaelis-Menten Model Fit", x = "Sites", y = "Species Richness", tag="B)") +
  theme(axis.text = element_text(size = 12), plot.title = element_text(size = 14)) + theme_publish()

png(filename="../SupplementaryFiles/Rarefaction_WesternEurope.png", 
    width=9.5, height=4.5, units="in",res=600)
grid.arrange(rarefaction_plot, model_plot, ncol = 2)
dev.off()

###CEE ####
df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5a_NSS_SpeciesData_Feb2025_withCEE.xlsx")
df_species$time_bins<-factor(df_species$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_species$time_bins2<-factor(df_species$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))
NSS_sp<-df_species[!is.na(df_species$GBIF_species),]
NSS_sp<-NSS_sp[!NSS_sp$GBIF_species=="NA",]

#Create time range variable for subsampling chronologies
NSS_sp$Time.Range<-as.numeric(NSS_sp$End_date_CE) - as.numeric(NSS_sp$Start_date_CE)
NSS_sp<-NSS_sp[!NSS_sp$Time.Range=="NA",]
NSS_sp<-NSS_sp %>%
  dplyr::mutate_at('Time.Range', as.numeric)

NSS_sp$count<-"1"
NSS_sp<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                 !NSS_sp$DB_Assemblage_ID == "NA")

length(unique(NSS_sp$DB_Assemblage_ID))
name<-"NSS_sp"
summarise_df(NSS_sp, name)

NSS_sp %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp_CEE<-NSS_sp %>% dplyr::filter(!NSS_sp$Temp.mid == "NA", 
                                    !NSS_sp$DB_Assemblage_ID == "NA",
                                    NSS_sp$Region=="Central & Eastern Europe")

tmp<-NSS_sp_CEE %>%
  dplyr::select("DB_Assemblage_ID","NISP","GBIF_species") %>%
  dplyr::mutate_at('NISP', as.numeric) 
fish_sp_JC<-reshape2::dcast(tmp %>%
                    dplyr::select("DB_Assemblage_ID","GBIF_species"), DB_Assemblage_ID ~ GBIF_species, fill = 0)
fish_sp_JC[fish_sp_JC != 0] <- 1
fish_sp_BC<-reshape2::dcast(tmp, DB_Assemblage_ID ~ GBIF_species, value.var="NISP", sum)
rm(tmp)
fish_sp_JC_sht<-fish_sp_JC[,c(-1)]
fish_sp_JC_sht[] <- lapply(fish_sp_JC_sht, function(x) {
  x <- as.numeric(x)   # Convert to numeric
  x[is.na(x)] <- 0     # Replace NAs with 0
  return(x)
})
fish_sp_BC_sht<-fish_sp_BC[,c(-1)]

set.seed(12)
boot_result <- boot(fish_sp_JC_sht, statistic = species_richness, R = 1000)
# Extract mean richness and confidence intervals
richness_mean <- colMeans(boot_result$t)  # Mean richness across bootstraps
richness_ci <- apply(boot_result$t, 2, quantile, probs = c(0.025, 0.975))  
richness_df <- data.frame(
  Site = rownames(fish_sp_JC_sht),
  Mean_Richness = richness_mean,
  Lower_CI = richness_ci[1, ],
  Upper_CI = richness_ci[2, ]
)

total_richness_per_site <- rowSums(fish_sp_JC_sht > 0)  # Count presence of species per site

# Create a combined dataframe for plotting mean and total richness
richness_comparison_df <- data.frame(
  Richness_Type = rep(c("Mean Richness", "Total Richness"), each = length(richness_df$Mean_Richness)),
  Richness = c(richness_df$Mean_Richness, total_richness_per_site)
)

# Violin plot of both mean and total richness
ggplot(richness_comparison_df, aes(x = Richness_Type, y = Richness, fill = Richness_Type)) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.6) +
  theme_minimal() +
  labs(title = "Violin Plot: Mean and Total Species Richness Across All Sites",
       x = "Richness Type",
       y = "Species Richness") +
  scale_fill_manual(values = c("blue", "green")) +
  theme(axis.text.x = element_text(size = 12))

# Rarefaction curve (sampling curve) for species richness
rarefaction_results <- specaccum(fish_sp_JC_sht, method = "rarefaction", permutations = 100)
rarefaction_data <- data.frame(
  Sites = rarefaction_results$sites,
  Richness = rarefaction_results$richness,
  Upper_CI = rarefaction_results$richness + rarefaction_results$sd,  # Upper bound of CI
  Lower_CI = rarefaction_results$richness - rarefaction_results$sd   # Lower bound of CI
)
ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = "purple", size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "purple", alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve with Confidence Intervals: Species Richness vs. Number of Sites",
       x = "Number of Sites Sampled",
       y = "Species Richness") +
  theme(axis.text = element_text(size = 12))

richness_df$sites<-seq(1:length(richness_df$Mean_Richness))
#Check whether curve has reached asymptote
model <- nls(Mean_Richness ~ a * sites / (b + sites), 
             data = richness_df, 
             start = list(a = max(richness_df$Mean_Richness), b = 1))
summary(model)
fitted_values <- predict(model, newdata = data.frame(sites = richness_df$Sites))
plot_data <- data.frame(sites = richness_df$sites, 
                        richness = richness_df$Mean_Richness, 
                        fitted_richness = fitted_values)

fixed_b_model <- nls(Mean_Richness ~ a * sites / (1 + sites), 
                     data = richness_df, 
                     start = list(a = max(richness_df$Mean_Richness)))

rarefaction_plot <- ggplot(rarefaction_data, aes(x = Sites, y = Richness)) +
  geom_line(color = DB[5], size = 1) +  # Line for richness
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill =  DB[5], alpha = 0.2) +  # Error ribbon for CI
  theme_minimal() +
  labs(title = "Rarefaction Curve",
       x = "Number of Sites Sampled",
       y = "Species Richness",
       tag = "A)") +
  theme(axis.text = element_text(size = 12)) + theme_publish()

model_plot <- ggplot(plot_data, aes(x = sites, y = richness)) +
  geom_point(color =  DB[3], alpha=0.7) +
  geom_line(aes(y = fitted_richness), color = DB[5], linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Michaelis-Menten Model Fit", x = "Sites", y = "Species Richness", tag="B)") +
  theme(axis.text = element_text(size = 12), plot.title = element_text(size = 14)) + theme_publish()

png(filename="../SupplementaryFiles/Rarefaction_CEE.png", 
    width=9.5, height=4.5, units="in",res=600)
grid.arrange(rarefaction_plot, model_plot, ncol = 2)
dev.off()

#Run NMDS - Family level ####
df_family<-read_excel("../SupplementaryFiles/SupplementaryTable6b_NSS_FamilyData_Feb2025_noCEE.xlsx")
df_family$time_bins<-factor(df_family$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_family$time_bins2<-factor(df_family$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))
NSS_sp<-df_family[!is.na(df_family$GBIF_family),]
NSS_sp<-NSS_sp[!NSS_sp$GBIF_family=="NA",]

#Create time range variable for subsampling chronologies
NSS_sp$Time.Range<-as.numeric(NSS_sp$End_date_CE) - as.numeric(NSS_sp$Start_date_CE)
NSS_sp<-NSS_sp[!NSS_sp$Time.Range=="NA",]
NSS_sp<-NSS_sp %>%
  dplyr::mutate_at('Time.Range', as.numeric)

NSS_sp$count<-"1"
NSS_sp<-NSS_sp %>% dplyr::filter(!NSS_sp$DB_Assemblage_ID == "NA")
length(unique(NSS_sp$DB_Assemblage_ID))
name<-"NSS_sp"

NSS_sp %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  group_by(Region) %>%
  distinct(GBIF_family) %>%
  tally()

NSS_sp<-NSS_sp %>%
  dplyr::mutate_at('Total_Fish_NISP', as.numeric) 

tmp<-NSS_sp %>%
  dplyr::select("DB_Assemblage_ID","NISP","GBIF_family","time_bins2") %>%
  dplyr::mutate_at('NISP', as.numeric) 

###Global diversity at the family level (NDY) ####
NSS_FAM<-read_excel("1_Data/NSS_Family_Data.xlsx", col_names = TRUE, col_types = "text")
names(NSS_FAM)
tmp<-NSS_FAM[,c(1,12,26)] %>%
  dplyr::mutate_at('NISP_counts', as.numeric) 
fish_sp_JC<-dcast(NSS_FAM[,c(1,12)], Assemblage_ID ~ GBIF_family)
fish_sp_BC<-dcast(tmp, Assemblage_ID ~ GBIF_family, value.var="NISP_counts", sum)
rm(tmp)
fish_sp_JC_sht<-fish_sp_JC[,c(-1)]
fish_sp_BC_sht<-fish_sp_BC[,c(-1)]

#Jaccard Index - Global (all species presence and absence data)
fish.bd.presence<-beta.div.comp(fish_sp_JC_sht,coef="J",quant=F)
fish.bd.presence$part

#Quantitative form of the sorensen index
fish.bd.abundance<-beta.div.comp(fish_sp_BC_sht,coef="J",quant=T)
fish.bd.abundance$part

#Local contribution of each site to BD - replacement
local.repl<-LCBD.comp(fish.bd.abundance$repl, sqrt.D = T)
local.repl
hist(local.repl$LCBD)
#Normal distribution of replacement across all sites, suggesting no particular sites are driving the replacement

#Local contribution of each site to BD - richness
local.rich<-LCBD.comp(fish.bd.abundance$rich, sqrt.D = T)
local.rich
hist(local.rich$LCBD)
#Skewed distribution of richness, suggesting some sites have higher than average richness than others

#Species contributions to beta diversity (This function will tell you if a site contributes significantly to beta diversity)
scbd<-beta.div(fish_sp_BC_sht, method="hellinger")
scbd
hist(scbd$SCBD)

####Jaccard Index ####
fish_sp_BC<-dcast(NSS_sp %>%
                    dplyr::select("DB_Assemblage_ID","GBIF_family","time_bins2"), DB_Assemblage_ID + time_bins2 ~ GBIF_family)
fish_sp_BC[,c(3:63)]<- ifelse(fish_sp_BC[,c(3:63)] > 0, 1, 0)
meta<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Region, time_bins2, Rural_urban_or_neither, Historic_Urban_Centre) %>% distinct()
threshold <- 0.015
fish_sp_BC_all<-fish_sp_BC
rownames(fish_sp_BC_all)<-fish_sp_BC_all$DB_Assemblage_ID
fish_sp_BC_all<-fish_sp_BC_all[,3:length(colnames(fish_sp_BC_all))] %>%
  mutate(across(everything(.), as.numeric))
fish_sp_BC_all_pruned<-fish_sp_BC_all[, colSums(fish_sp_BC_all, na.rm = T) > (threshold * nrow(fish_sp_BC_all))]
fish_sp_BC_all_pruned_jaccard<- vegdist(fish_sp_BC_all_pruned, method = "jaccard")
nmds_result_all <- metaMDS(fish_sp_BC_all_pruned_jaccard)
nmds_result_all <- as.data.frame(nmds_result_all$points)
nmds_result_all$DB_Assemblage_ID<-row.names(nmds_result_all)
nmds_result_all_out_jaccard<-merge(nmds_result_all, meta)
nmds_result_all_out_jaccard$seq<-seq(1:1511)
adonis_result <- adonis2(fish_sp_BC_all_pruned_jaccard ~ Region, data = meta, method = "jaccard")
adonis_result_b <- adonis2(fish_sp_BC_all_pruned_jaccard ~ time_bins2, data = meta, method = "jaccard")
adonis_result_c <- adonis2(fish_sp_BC_all_pruned_jaccard ~ time_bins2 + Region, data = meta, method = "jaccard")
adonis_result_d <- adonis2(fish_sp_BC_all_pruned_jaccard ~ time_bins2*Region, data = meta, method = "jaccard")

ggplot(nmds_result_all_out_jaccard,
       aes(x = MDS1, y = MDS2, color = Region)) +
  geom_point(size = 1.5, alpha=0.8) +
  scale_color_manual(values = DB2[2:6]) + 
  theme_minimal() +
  stat_ellipse(type='t',size = 1, level = 0.9) + 
  labs(title = "NMDS per time bin") + facet_wrap(~time_bins2, scales="free")

###Bray-curtis time*region
fish_sp_BC<-dcast(NSS_sp %>%
                    dplyr::select("DB_Assemblage_ID","GBIF_family","time_bins2"), DB_Assemblage_ID + time_bins2 ~ GBIF_family)
meta<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Region, time_bins2, Rural_urban_or_neither, Historic_Urban_Centre) %>% distinct()
threshold <- 0.015
fish_sp_BC_all<-fish_sp_BC
rownames(fish_sp_BC_all)<-fish_sp_BC_all$DB_Assemblage_ID
fish_sp_BC_all<-fish_sp_BC_all[,3:length(colnames(fish_sp_BC_all))] %>%
  mutate(across(everything(.), as.numeric))
fish_sp_BC_all_pruned<-fish_sp_BC_all[, colSums(fish_sp_BC_all, na.rm = T) > (threshold * nrow(fish_sp_BC_all))]
nmds_result_all <- metaMDS(fish_sp_BC_all_pruned, distance = "bray", trymax = 100)
nmds_result_all <- as.data.frame(nmds_result_all$points)
nmds_result_all$DB_Assemblage_ID<-row.names(nmds_result_all)
nmds_result_all_out<-merge(nmds_result_all, meta)
nmds_result_all_out$seq<-seq(1:1511)

ggplot(nmds_result_all_out %>%
         filter(!seq == 3,
                !seq == 809,
                !seq == 97,
                !seq == 1031,
                !seq == 45),
       aes(x = MDS1, y = MDS2, color = Region)) +
  geom_point(size = 1.5, alpha=0.8) +
  scale_color_manual(values = DB2[2:6]) + 
  theme_minimal() +
  stat_ellipse(type='t',size = 1, level = 0.95) + 
  labs(title = "NMDS per time bin") + facet_wrap(~time_bins2, scales="free")

ggplot(nmds_result_all_out %>%
         filter(!seq == 3,
                !seq == 809,
                !seq == 97,
                !seq == 1031,
                !seq == 45),
       aes(x = MDS1, y = MDS2, color = Region)) +
  geom_point(size = 1.5, alpha=0.8) +
  scale_color_manual(values = DB2[2:6]) + 
  theme_minimal() +
  stat_ellipse(type='t',size = 1, level = 0.75) + 
  labs(title = "NMDS per time bin") + facet_grid(Rural_urban_or_neither~time_bins2, scales="free")

###Calculate ellipse overlap
library(car)
nmds_result_all_out$group<-as.character(paste(nmds_result_all_out$Region, nmds_result_all_out$time_bins2))
nmds_result_all_out$count<-1
tm<-nmds_result_all_out %>% group_by(group) %>% tally()
sim_results <- run_simulation(nmds_result_all_out, iterations = 100)
f_test_output <- perform_f_test(sim_results)
print(sim_results$summary)
print(f_test_output)


#Jaccard Index #
fish.bd.presence<-beta.div.comp(fish_sp_JC_sht,coef="J",quant=F)
fish.bd.presence$part

#Quantitative form of the sorensen index
fish.bd.abundance<-beta.div.comp(fish_sp_BC_sht,coef="J",quant=T)
fish.bd.abundance$part

#Local contribution of each site to BD - replacement
local.repl<-LCBD.comp(fish.bd.abundance$repl, sqrt.D = T)
local.repl
hist(local.repl$LCBD)
#Almost normal distribution of replacement across all sites, suggesting no particular sites are driving the replacement

#Local contribution of each site to BD - richness
local.rich<-LCBD.comp(fish.bd.abundance$rich, sqrt.D = T)
local.rich
hist(local.rich$LCBD)

#Jaccard Index #
fish.bd.presence<-beta.div.comp(fish_sp_JC_sht,coef="J",quant=F)
fish.bd.presence$part

#Quantitative form of the sorensen index
fish.bd.abundance<-beta.div.comp(fish_sp_BC_sht,coef="J",quant=T)
fish.bd.abundance$part

#Local contribution of each site to BD - replacement
local.repl<-LCBD.comp(fish.bd.abundance$repl, sqrt.D = T)
local.repl
hist(local.repl$LCBD)
#Normal distribution of replacement across all sites, suggesting no particular sites are driving the replacement

#Local contribution of each site to BD - richness
local.rich<-LCBD.comp(fish.bd.abundance$rich, sqrt.D = T)
local.rich
hist(local.rich$LCBD)
#Skewed distribution of richness, suggesting some sites have higher than average richness than others

#Species contributions to beta diversity (This function will tell you if a site contributes significantly to beta diversity)
scbd<-beta.div(fish_sp_BC_sht, method="hellinger")
scbd
hist(scbd$SCBD)


#Jaccard Index #
fish.bd.presence<-beta.div.comp(fish_sp_JC_sht,coef="J",quant=F)
fish.bd.presence$part

#Quantitative form of the sorensen index
fish.bd.abundance<-beta.div.comp(fish_sp_BC_sht,coef="J",quant=T)
fish.bd.abundance$part

#Local contribution of each site to BD - replacement
local.repl<-LCBD.comp(fish.bd.abundance$repl, sqrt.D = T)
local.repl
hist(local.repl$LCBD)
#Almost normal distribution of replacement across all sites, suggesting no particular sites are driving the replacement

#Local contribution of each site to BD - richness
local.rich<-LCBD.comp(fish.bd.abundance$rich, sqrt.D = T)
local.rich
hist(local.rich$LCBD)

