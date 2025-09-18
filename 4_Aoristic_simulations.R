###File 4. simulate aoristic time bins ####
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
source("MANUSCRIPT_functions.R")

#3. setup data
#n_cols <- readxl::read_xlsx("NSS_SpeciesData_Sept2025_noPE.xlsx", n_max = 0) %>% ncol()
#df_species <- readxl::read_xlsx(
#  "NSS_SpeciesData_Sept2025_noPE.xlsx",
#  col_types = rep("text", n_cols)   # force all columns to character)
df_species<-df_new_with_cities_b

df_species$NISP<-as.numeric(df_species$NISP)
df_species<-df_species[!is.na(df_species$GBIF_species),]
df_species<-df_species[!df_species$GBIF_species=="NA",]

df_species %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

# columns that should be numeric
num_cols <- c(
  "Decimal_Latitude","Decimal_Longitude",
  "start_date_CE_final_num","end_date_CE_final_num",
  "min_sieve_size_final_num","NISP","total_NISP", 
  "Trait_Temp_Max","Trait_Temp_Min",
  "Trait_Temp_Fishbase_Median","Trait_Temp_Fishbase_Range",
  "Trait_Temp_Fishbase_Env_Low","Trait_Temp_Fishbase_Env_High",
  "Trait_Temp_Fishbase_Mid"
)

cols_to_fix<-c(  "Decimal_Latitude","Decimal_Longitude",
                 "start_date_CE_final_num","end_date_CE_final_num",
                 "min_sieve_size_final_num","NISP","total_NISP", 
                 "Trait_Temp_Max","Trait_Temp_Min",
                 "Trait_Temp_Fishbase_Median","Trait_Temp_Fishbase_Range",
                 "Trait_Temp_Fishbase_Env_Low","Trait_Temp_Fishbase_Env_High",
                 "Trait_Temp_Fishbase_Mid")
df_species <- df_species %>%
  dplyr::mutate(across(all_of(cols_to_fix), ~ as.numeric(gsub(",", ".", .x))))
df_species$NISP<-df_species$NISP %>% as.numeric()

# 2. Compute Time.mid and Time.Range ####
df_species$Time.mid <- round((as.numeric(df_species$start_date_CE_final_num) + as.numeric(df_species$end_date_CE_final_num)) / 2 , -2)

df_species$Time.Range <- as.numeric(df_species$end_date_CE_final_num) - as.numeric(df_species$start_date_CE_final_num)

df_species %>%
  dplyr::group_by(Region) %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69")
df_species$time_bins<-factor(df_species$time_bins, levels=c("<600","600-900","900-1200","1200-1500",">1500"))
df_species$time_bins2<-factor(df_species$time_bins2, levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"))

#3. Remove rows without temperature trait data ####
NSS_sp<-df_species
NSS_sp<-NSS_sp[!NSS_sp$Trait_Temp_Fishbase_Env_Low=="NA",]
NSS_sp<-NSS_sp[!NSS_sp$Trait_Temp_Fishbase_Env_High=="NA",]

NSS_sp %>% 
  dplyr::summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(DB_Assemblage_ID)
  )
NSS_sp<-NSS_sp[!is.na(NSS_sp$Region),]
NSS_sp %>%
  dplyr::group_by(Region) %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))
writexl::write_xlsx(NSS_sp, "NSS_sp_tmp.xlsx")

NSS_sp <- NSS_sp[setdiff(names(NSS_sp), 
                         c("Minimum_sieve_mesh_mm", "Has_taxa", "minimum_sieve_size_num"))]

#5. Check temperatures are numeric and remove outlier species or those with zero counts or mixed life history strategies ####
NSS_sp$Trait_LifeHistory<-as.character(NSS_sp$Trait_LifeHistory)
NSS_sp<-NSS_sp[!NSS_sp$Trait_LifeHistory=="Unclassified; mixed strategies",]
NSS_sp<-NSS_sp[!NSS_sp$Trait_LifeHistory=="Amphidromous",]
NSS_sp<-NSS_sp[!is.na(NSS_sp$Trait_LifeHistory),]
NSS_sp<-NSS_sp[!NSS_sp$Trait_LifeHistory=="NA",]
names(NSS_sp)

hist(as.numeric(NSS_sp$Trait_Temp_Fishbase_Range), breaks=12)
hist(as.numeric(NSS_sp$Trait_Temp_Fishbase_Mid), breaks=12)
species_ranges <- NSS_sp %>%
  dplyr::group_by(GBIF_species) %>%
  dplyr::summarise(Trait_Temp_Fishbase_Range = max(Trait_Temp_Fishbase_Range, na.rm = TRUE), .groups = "drop")

###REMOVE TEMP RANGE MAX ####
thr <- quantile(species_ranges$Trait_Temp_Fishbase_Range, probs = 0.975, na.rm = TRUE)
species_to_drop <- species_ranges %>%
  dplyr::filter(Trait_Temp_Fishbase_Range >= thr) %>%
  dplyr::pull(GBIF_species)
NSS_sp <- NSS_sp %>%
  filter(!GBIF_species %in% species_to_drop)
hist(NSS_sp$Trait_Temp_Fishbase_Range, breaks = 12)

NSS_sp %>%
  dplyr::group_by(Region) %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

#6. Species and NISP counts in remaining dataset ####
NSS_sp$count<-"1"

NSS_sp %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::group_by(count, Region) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))
NSS_sp %>% 
  dplyr::select(GBIF_species,GBIF_family,GBIF_class, GBIF_order, Region, count) %>% 
  dplyr::group_by(count) %>%
  dplyr::summarize(Species = length(unique(GBIF_species)),
                   Families = length(unique(GBIF_family)),
                   Class = length(unique(GBIF_class)),
                   Order = length(unique(GBIF_order)))

NSS_sp %>% 
  dplyr::select(time_bins2, NISP, GBIF_species, Region) %>% 
  dplyr::group_by(Region) %>%
  dplyr::summarize(NISP = sum(NISP))

NSS_sp %>% 
  dplyr::select(time_bins2, NISP, GBIF_species) %>% 
  dplyr::summarize(NISP = sum(NISP))

NSS_sp %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>% 
  dplyr::select(time_bins2, NISP, GBIF_species) %>% 
  dplyr::summarize(NISP = sum(NISP))

#100 year intervals - Global and regional aMTC (no PE) ####
###RUN1###
#1. 500 subsamples - 100-yr time bins (All Data)
# 1) Create the blank columns
#Create blank columns in dataframe where simulated timestamps will be generated. Change the number below to edit number of simulations for aMTC, currently set at 200
repeat_size = 1000
repeats <- paste0("rep", seq_len(repeat_size))
test = repeats
bin_size = 100

names(NSS_sp)
##Subset unique assemblages 
tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, start_date_CE_final_num, end_date_CE_final_num) %>%
  distinct()

NSS_sp %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Total_NISP       = sum(NISP, na.rm = TRUE),
    Assemblage_count = dplyr::n_distinct(DB_Assemblage_ID),
    .groups = "drop"
  )

tmp1 <- tmp1 %>% 
  dplyr::mutate(
    start_date_CE_final_num = as.numeric(start_date_CE_final_num),
    end_date_CE_final_num   = as.numeric(end_date_CE_final_num)
  )

NSS_sp_tmp<-tmp1

for (i in seq_along(repeats)) {
  NSS_sp_tmp[, repeats[i]] <- 0L
}

# 2) Populate with one simulated year per row/column (non-vectorised)
bin_size <- 100
for (i in seq_along(repeats)) {
  # optional: set.seed(1000 + i)
  
  for (j in seq_len(nrow(NSS_sp_tmp))) {
    NSS_sp_tmp[j, repeats[i]] <- generate_random_time_data(
      NSS_sp_tmp$start_date_CE_final_num[j],
      NSS_sp_tmp$end_date_CE_final_num[j],
      time_bin_size = bin_size,
      reps = 1
    )
  }
}

###Rejoin NSS_sp and NSS_sp_tmp
NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 4:1003)], by="DB_Assemblage_ID", all.x=T)
rm(tmp1, NSS_sp_tmp)

#Save species data with 100-yr subsets
save(NSS_sp_tmp_all, file = "df_species_with_100yrtime_increments_Sept2025.RData")
write_xlsx(NSS_sp_tmp_all, "df_species_with_100yrtime_increments_Sept2025.xlsx")

NSS_sp_tmp_all %>%
  dplyr::group_by(Region) %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

###For abba
write_xlsx(NSS_sp_tmp_all[,c(1:50)], "df_species_with_100yrtime_increments_Sept2025_Abba.xlsx")
