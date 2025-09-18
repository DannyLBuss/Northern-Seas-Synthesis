####13th Aug 2025 - NSS analysis
#Author: Dr Danny L Buss
###Merge datasets
#1. load libraries and setwd ####
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(bayesplot)
  library(stringr)
  library(writexl)
  library(readr)
  library(tidyr)
})
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../ThemePublish.R")
source("MANUSCRIPT_functions.R")

#2. load data ####
#df<-read_xlsx("All_Data_with_Chronology_Sept2025.xlsx", col_names = TRUE, col_types = "text")
df<-df_final
names(df)<-gsub(" ","_",names(df))
df<-df %>% dplyr::select(-any_of(c("Site_code","Year_reported")))

#df$Region<-dplyr::recode(df$Country,
#                             "England" = 'Britain & Ireland',
#                             "Scotland" = 'Britain & Ireland',
#                             "Norway" = 'Scandinavia',
#                             "Sweden" = 'Scandinavia',
#                             "Denmark" = 'Scandinavia',
#                             "Germany" = 'Western Europe',
#                             "Belgium" = 'Western Europe',
#                             "Netherlands" = 'Western Europe',
#                             "Ireland" = 'Britain & Ireland',
#                             "Estonia" = 'Poland & Estonia',
#                             "Poland" = 'Poland & Estonia',
#                             "Northern Ireland" = 'Britain & Ireland')
#
cols_to_fix<-c('Start_date_CE', 'End_date_CE',
               'Decimal_Latitude','Decimal_Longitude', 
               'ID_for_linked_rows_differing_only_by_recovery', 
               'NISP','Trait_Temp_Max','Trait_Temp_Min', 'Trait_Temp_Fishbase_Median',
               'Trait_Temp_Fishbase_Range','Trait_Temp_Fishbase_Env_Low','Trait_Temp_Fishbase_Env_High',
               'Trait_Temp_Fishbase_Mid')
df <- df %>%
  dplyr::mutate(across(all_of(cols_to_fix), ~ as.numeric(gsub(",", ".", .x))))
df$NISP<-df$NISP %>% as.numeric()

df %>%
  dplyr::group_by(region) %>%
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(Lumped_NSS_IDs))

df %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(Lumped_NSS_IDs))

df %>%
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(Lumped_NSS_IDs))

#Convert georeferences to three decimal places
df$Decimal_Latitude<-round(as.numeric(df$Decimal_Latitude),3)/1000
df$Decimal_Longitude<-round(as.numeric(df$Decimal_Longitude),3)/1000

#3. count number of assemblages ####
df %>% 
  summarise(UniqueCount = n_distinct(Lumped_NSS_IDs))

# Count total NISP
df$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df$NISP<-df$NISP %>% as.numeric()

df_merge<-df

df_merge$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df_merge$Recovery<-df_merge$`Recovery_method_(hand-collected;_sieved;_both;_unknown)`

df_merge$DB_ID3 <- with(
  df_merge,
  paste(
    Settlement_modern_name, Site_name, Decimal_Latitude, Decimal_Longitude,
    Country, Site_type_using_local_categories, Rural_urban_or_neither, Recovery, DB_sieved, Database, 
    sep = ";"
  )
)

#Unidentified specimens
df_old_short<-df_merge %>% dplyr::select(Lumped_NSS_IDs, Number_of_unidentified_fish_specimens) %>% distinct()
names(df_old_short)[2]<-"unidentified_fish_specimens"
df_old_short$unidentified_fish_specimens<-as.character(df_old_short$unidentified_fish_specimens)
df_old_short<-aggregate(unidentified_fish_specimens ~ Lumped_NSS_IDs, data = df_old_short, paste, collapse = ",")
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#Totals
df_old_short<-df_merge %>% dplyr::select(Lumped_NSS_IDs, Total_Fish_NISP) %>% distinct()
names(df_old_short)[2]<-"Total_NISP"
df_old_short$Total_NISP<-as.character(df_old_short$Total_NISP)
df_old_short<-aggregate(Total_NISP ~ Lumped_NSS_IDs, data = df_old_short, paste, collapse = ",")
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#Local site types
df_old_short<-df_merge %>% dplyr::select(Lumped_NSS_IDs, Site_type_using_local_categories) %>% distinct()
names(df_old_short)[2]<-"Local_site_type"
df_old_short$Local_site_type<-as.character(df_old_short$Local_site_type)
df_old_short<-aggregate(Local_site_type ~ Lumped_NSS_IDs, data = df_old_short, paste, collapse = ",")
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#Sieve_size_min
df_old_short<-df_merge %>% dplyr::select(Lumped_NSS_IDs, Minimum_sieve_mesh_mm) %>% distinct()
names(df_old_short)[2]<-"Minimum_sieve_mesh_mm"
df_old_short$Minimum_sieve_mesh_mm<-as.character(df_old_short$Minimum_sieve_mesh_mm)
df_old_short <- aggregate(Minimum_sieve_mesh_mm ~ Lumped_NSS_IDs, 
                          data = df_old_short, 
                          FUN = function(x) {
                            c(pasted = paste(x, collapse = ","),
                              minval = min(as.numeric(x)))
                          })

# Split the combined result into separate columns
df_old_short$Minimum_sieve_mesh_mm_agg <- df_old_short$Minimum_sieve_mesh_mm[, "minval"]
df_old_short$Minimum_sieve_mesh_mm <- df_old_short$Minimum_sieve_mesh_mm[, "pasted"]
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#Rural_urban_or_neither
df_old_short<-df_merge %>% dplyr::select(Lumped_NSS_IDs, Rural_urban_or_neither) %>% distinct()
names(df_old_short)[2]<-"Rural_or_urban"
df_old_short$Rural_or_urban<-as.character(df_old_short$Rural_or_urban)
df_old_short<-aggregate(Rural_or_urban ~ Lumped_NSS_IDs, data = df_old_short, paste, collapse = ",")
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#5. Convert numeric variables to text for aggregations ####
df_merge <-df_merge %>% mutate_at(c('Start_date_CE', 'End_date_CE',
                        'Decimal_Latitude','Decimal_Longitude', 
                        'ID_for_linked_rows_differing_only_by_recovery', 
                        'NISP','Trait_Temp_Max','Trait_Temp_Min', 'Trait_Temp_Fishbase_Median',
                        'Trait_Temp_Fishbase_Range','Trait_Temp_Fishbase_Env_Low','Trait_Temp_Fishbase_Env_High',
                        'Trait_Temp_Fishbase_Mid'), as.character)
#Remove special characters from column names
names(df_merge)<-gsub(" ","_", names(df_merge))
names(df_merge)<-gsub("\\(","", names(df_merge))
names(df_merge)<-gsub("\\)","", names(df_merge))

#Remove unwanted columns for stacking
remove_cols<-c("DB_ID2","DB_ID4","Assemblage_ID","NSS_unique_ID","Input_by","ID_number_within_dataset",
               "File_name","Assemblage_or_sub-assemblage","County_province_or_state","Georef_source",
               "Context_types","ID_for_linked_rows_differing_only_by_recovery","Part_of_sieve_stack",
               "Recovery_method_hand-collected;_sieved;_both;_unknown","Sieve_sizes_and_recovery_details_with_original_mesh_units_if_not_mm",
               "Maximum_sieve_mesh_mm","Lumped_by_phase_or_split_by_context",
               "Unpublished_zooarchaeology_reference_J_Arch_Sci_format","Published_zooarchaeology_reference_J_Arch_Sci_format",
               "Archaeological_reference_for_chronology_etc_not_in_zooarch_refs_J_Arch_Sci_format","Analyst_name",
               "Dataset_contact_name","IP_status","Skeletal_element_data_available_y_or_n","Measurement_or_fish_size_data_available_y_or_n",
               "Data_quality_green_amber_or_red","Explanation_of_amber_or_red_data_quality_noting_what_variables_and_why",
               "General_comments","Taxa_with_only_presence_data","Phase_or_equivalent_optional","Context_or_equivalent_optional",
               "Sample_optional","Sediment_volumes_available_y_or_n_optional","Associated_mammal_NISP_available_y;_n;_not_matching_or_enter_value_if_easily_to_hand_optional",
               "Associated_bird_NISP_available_y;_n;_not_matching_or_enter_value_if_easily_to_hand_optional","Extra_column_for_dataset-specific_variable_optional"
               ,"Extra_column_for_dataset-specific_variable_optional","Definition_of_extra_column_variable_optional",
               "Number_of_unidentified_fish_specimens","Total_Fish_NISP","Rural_urban_or_neither","Site_type_using_local_categories", "database","Original_name","Taxa_cleaned","Minimum_sieve_mesh",
               "Data_quality","unidentified_fish_specimens","Total_NISP","Minimum_sieve_mesh_mm_agg")

#Minimum_sieve_mesh_mm	Has_taxa	Trait_Temp_Fishbase_Median	Trait_Temp_Fishbase_Range	Trait_Temp_Fishbase_Source	Trait_Temp_Fishbase_Env_Low	Trait_Temp_Fishbase_Env_High	Trait_Temp_Fishbase_Mid	Trait_Temp_Fishbase_Source_Choice	minimum_sieve_size_num
df_merge<-df_merge %>% dplyr::select(-any_of(remove_cols))

df_merge$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_merge$Lumped_NSS_IDs))

####NISP BY DATA QUALITY ####
df_merge %>%
  dplyr::mutate(
    quality = as.character(quality),
    NISP = as.numeric(NISP)
  ) %>%
  dplyr::group_by(quality) %>%
  dplyr::summarise(total_NISP = sum(NISP, na.rm = TRUE))

### QUALITY FILTER HERE ####
df_merge<-df_merge[df_merge$quality==1,]

#set up aggregate function inputs ####
collapse_mode <- "span_ok"
group_keys <- c(
  "Settlement_modern_name","Site_name","Decimal_Latitude","Decimal_Longitude",
  "Country","Local_site_type","Rural_or_urban","Recovery","DB_sieved","Database"
)
taxa_cols <- c(
  "GBIF_species","GBIF_genus","GBIF_family","GBIF_order",
  "GBIF_chrondrichthyes_superorder","GBIF_chondricthyes_infraclass","GBIF__chondricthyes_subclass",
  "GBIF_class","GBIF_phylum","GBIF_kingdom","GBIF_level","Trait_Habitat","Trait_LifeHistory","Trait_endangered","Trait_Temp_Max",
  "Trait_Temp_Min","Known_Major_Trade_J_Barrett","Known_Aquaculture_R_Hoffmann","Trait_Temp_Fishbase_Median","Trait_Temp_Fishbase_Range",
  "Trait_Temp_Fishbase_Source","Trait_Temp_Fishbase_Env_Low","Trait_Temp_Fishbase_Env_High","Trait_Temp_Fishbase_Mid","Trait_Temp_Fishbase_Source_Choice"
)

# keep only keys that exist (prevents 'columns do not exist' errors)
group_keys <- intersect(group_keys, names(df_merge))
taxa_cols  <- intersect(taxa_cols,  names(df_merge))

#aggregate chron overlaps
df_out <- df_merge %>%
  dplyr::mutate(
    .start_num = readr::parse_double(base::as.character(Start_date_CE)),
    .end_num   = readr::parse_double(base::as.character(End_date_CE)),
    minimum_sieve_size_num = readr::parse_number(Minimum_sieve_mesh_mm),
    NISP_num = readr::parse_double(base::as.character(NISP))
  ) %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(group_keys))) %>%
  dplyr::mutate(
    .g_min_start = { x <- .start_num[!base::is.na(.start_num)]; if (base::length(x)) base::min(x) else NA_real_ },
    .g_max_end   = { x <- .end_num[!base::is.na(.end_num)];     if (base::length(x)) base::max(x) else NA_real_ },
    .g_max_start = { x <- .start_num[!base::is.na(.start_num)]; if (base::length(x)) base::max(x) else NA_real_ },
    .g_min_end   = { x <- .end_num[!base::is.na(.end_num)];     if (base::length(x)) base::min(x) else NA_real_ },
    .total_span  = .g_max_end - .g_min_start,
    .span_ok     = !base::is.na(.total_span) & .total_span <= 150,
    .overlap_start = .g_max_start,
    .overlap_end   = .g_min_end,
    .has_overlap   = !base::is.na(.overlap_start) & !base::is.na(.overlap_end) &
      (.overlap_start <= .overlap_end),
    aggregate = .has_overlap & .span_ok,
    collapse_group = dplyr::case_when(
      collapse_mode == "overlap" ~ (.has_overlap & .span_ok),
      collapse_mode == "span_ok" ~ (.span_ok),
      collapse_mode == "always"  ~ TRUE,
      TRUE                       ~ FALSE
    ),
    .in_overlap = (.has_overlap) &
      !base::is.na(.start_num) & !base::is.na(.end_num) &
      (.start_num <= .overlap_end) & (.end_num >= .overlap_start)
  ) %>%
  dplyr::mutate(
    NSS_IDs_overlap = dplyr::if_else(
      .has_overlap,
      {
        ids <- base::unique(stats::na.omit(Lumped_NSS_IDs[.in_overlap]))
        if (base::length(ids)) base::paste(ids, collapse = ";") else NA_character_
      },
      NA_character_
    ),
    NSS_IDs_group_all = {
      ids <- base::unique(stats::na.omit(Lumped_NSS_IDs))
      if (base::length(ids)) base::paste(ids, collapse = ";") else NA_character_
    }
  ) %>%
  dplyr::mutate(
    start_date_CE_final_num  = dplyr::if_else(aggregate, .g_min_start, .start_num),
    end_date_CE_final_num    = dplyr::if_else(aggregate, .g_max_end,   .end_num),
    min_sieve_size_final_num = dplyr::if_else(collapse_group,
                                              { x <- minimum_sieve_size_num[!base::is.na(minimum_sieve_size_num)];
                                              if (base::length(x)) base::min(x) else NA_real_ },
                                              minimum_sieve_size_num),
    NSS_IDs_final = dplyr::if_else(.has_overlap, NSS_IDs_overlap, NSS_IDs_group_all)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    -dplyr::starts_with(".g_"),
    -.total_span, -.span_ok,
    -dplyr::starts_with(".overlap"),
    -.has_overlap, -.in_overlap,
    -NSS_IDs_group_all,
    -.start_num, -.end_num
  )

group_keys <- intersect(group_keys, names(df_merge))
taxa_cols  <- intersect(taxa_cols,  names(df_merge))

remove <- c(
  "DB_ID4","lowest_start_overlap","highest_end_overlap",
  "DB_ID4_Group","DB_ID3","start_date_CE","end_date_CE"
)
remove <- intersect(remove, names(df_out))

group_keys <- intersect(group_keys, names(df_out))
taxa_cols  <- intersect(taxa_cols,  names(df_out))

df_base <- df_out %>%
  dplyr::select(-dplyr::any_of(remove)) %>%
  dplyr::mutate(
    NISP_num = ifelse(is.na(NISP_num),
                      readr::parse_double(as.character(NISP)),
                      NISP_num)
  ) %>%
  dplyr::distinct()

gt_keys <- c(group_keys, taxa_cols)
gt_keys <- intersect(gt_keys, names(df_base))

total_by_gt <- df_base %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(gt_keys))) %>%
  dplyr::summarise(
    Total_NISP_num = if (all(is.na(NISP_num))) NA_real_ else sum(NISP_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Total_NISP = dplyr::if_else(is.na(Total_NISP_num), NA_character_, as.character(Total_NISP_num))) %>%
  dplyr::select(-Total_NISP_num)

id_cols <- unique(c(
  group_keys, taxa_cols,
  "start_date_CE_final_num","end_date_CE_final_num","min_sieve_size_final_num","NSS_IDs_final"
))
id_cols <- intersect(id_cols, names(df_base))

collapsed <- df_base %>%
  dplyr::filter(aggregate %in% TRUE & collapse_group %in% TRUE) %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
  dplyr::summarise(
    NISP_num = if (all(is.na(NISP_num))) NA_real_ else sum(NISP_num, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  dplyr::left_join(total_by_gt, by = gt_keys) %>%
  dplyr::mutate(
    aggregate      = TRUE,
    collapse_group = TRUE
  ) %>%
  dplyr::rename(
    NISP            = NISP_num,
    Lumped_NSS_IDs  = NSS_IDs_final
  )

not_collapsed <- df_base %>%
  dplyr::filter(!(aggregate %in% TRUE & collapse_group %in% TRUE)) %>%
  dplyr::mutate(Total_NISP = NA_character_) %>% mutate(NISP = as.numeric(NISP))

#Final dataframe following new chron overlapped aggregations
df_result <- dplyr::bind_rows(collapsed, not_collapsed)

not_collapsed$NISP %>% as.numeric() %>%
  sum(na.rm=T)

length(unique(not_collapsed$Lumped_NSS_IDs))

collapsed$NISP %>% as.numeric() %>%
  sum(na.rm=T)

length(unique(collapsed$Lumped_NSS_IDs))

#Recalculate total NISP
df_result <- dplyr::relocate(df_result, NISP, Total_NISP, .after = NISP)

df_result$NISP %>% as.numeric() %>%
  sum(na.rm=T)

length(unique(df_result$Lumped_NSS_IDs))

df <- df %>%
  dplyr::group_by(Lumped_NSS_IDs) %>%
  dplyr::filter(sum(NISP, na.rm = TRUE) >= 1) %>%
  ungroup()

id_cols <- unique(c(
  "Settlement_modern_name","Site_name","Decimal_Latitude","Decimal_Longitude",
  "Country","Local_site_type","Rural_or_urban","Recovery","DB_sieved",
  "aggregate","start_date_CE_final_num","end_date_CE_final_num","min_sieve_size_final_num",
  "GBIF_species","GBIF_genus","GBIF_family","GBIF_order",
  "GBIF_chrondrichthyes_superorder","GBIF_chondricthyes_infraclass","GBIF__chondricthyes_subclass",
  "GBIF_class","GBIF_phylum","GBIF_kingdom","GBIF_level","Trait_Habitat","Trait_LifeHistory","Trait_endangered","Trait_Temp_Max",
  "Trait_Temp_Min","Known_Major_Trade_J_Barrett","Known_Aquaculture_R_Hoffmann", "NSS_IDs_final","DB_Assemblage_ID"
))

# Remove assemblages that are chronological outliers ####
df_result<-df_result[df_result$start_date_CE_final_num |> as.numeric() >= -1,]
df_result<-df_result[df_result$start_date_CE_final_num |> as.numeric() <= 1850,]
df_result<-df_result[df_result$end_date_CE_final_num |> as.numeric() <= 1900,]
ggplot(df_result, aes(as.numeric(start_date_CE_final_num))) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "white") +
  scale_x_continuous(limits = c(1, 1901), breaks = seq(0, 2000, 500)) +
  labs(title = "Start dates",
       x = "Values",
       y = "Count") +
  theme_publish()
ggplot(df_result, aes(as.numeric(end_date_CE_final_num))) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "white") +
  scale_x_continuous(limits = c(1, 1901), breaks = seq(0, 2000, 500)) +
  labs(title = "End Dates",
       x = "Values",
       y = "Count") +
  theme_publish()

df_result$NISP %>% as.numeric() %>%
  sum(na.rm=T)

length(unique(df_result$Lumped_NSS_IDs))

#N=2083
#NISP=1802193

# Remove assemblages with low or v. high fish counts ####
#Recalculate total fish NISP
df_summed <- df_result %>%
  dplyr::group_by(
    Settlement_modern_name,
    Site_name,
    Decimal_Latitude,
    Decimal_Longitude,
    Country,
    Local_site_type,
    Rural_or_urban,
    Recovery,
    DB_sieved,
    aggregate,
    start_date_CE_final_num,
    end_date_CE_final_num,
    min_sieve_size_final_num,
    NSS_IDs_final
  ) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%
  ungroup()

# Remove assemblages with chronological range of more than 300
df_tmp<- df_summed %>%
  filter((end_date_CE_final_num |> as.numeric() - start_date_CE_final_num |> as.numeric()) <= 350) %>%
  filter(complete.cases(Decimal_Longitude, Decimal_Latitude))
df_tmp$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_tmp$Lumped_NSS_IDs))
#N=1536
#NISP=1378018

df_tmp$region<-dplyr::recode(df_tmp$Country,
                          "England" = 'Britain & Ireland',
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

df_tmp %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    Total_NISP       = sum(NISP, na.rm = TRUE),
    Assemblage_count = dplyr::n_distinct(Lumped_NSS_IDs),
    .groups = "drop"
  )

df_tmp %>%
  dplyr::group_by(region) %>%
  dplyr::filter(region=="Britain & Ireland",
                !GBIF_species=="NA",
                !is.na(GBIF_species)) %>%
  dplyr::summarise(
    Total_NISP       = sum(NISP, na.rm = TRUE),
    Assemblage_count = dplyr::n_distinct(Lumped_NSS_IDs),
    .groups = "drop"
  )

df_tmp %>%
  dplyr::group_by(region) %>%
  dplyr::filter(region=="Western Europe",
                !GBIF_species=="NA",
                !is.na(GBIF_species)) %>%
  dplyr::summarise(
    Total_NISP       = sum(NISP, na.rm = TRUE),
    Assemblage_count = dplyr::n_distinct(Lumped_NSS_IDs),
    .groups = "drop"
  )

df_tmp %>%
  dplyr::group_by(region) %>%
  dplyr::filter(region=="Scandinavia",
                !GBIF_species=="NA",
                !is.na(GBIF_species)) %>%
  dplyr::summarise(
    Total_NISP       = sum(NISP, na.rm = TRUE),
    Assemblage_count = dplyr::n_distinct(Lumped_NSS_IDs),
    .groups = "drop"
  )

df_tmp %>% 
  summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(Lumped_NSS_IDs)
  )

df_tmp$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_tmp$Lumped_NSS_IDs))

#N=1521
#NISP=1,378,018

###Add regions, temperature and time midpoints, and timebins
df_species<-df_tmp
df_species$Temp.mid<-as.character((( as.numeric(df_species$Trait_Temp_Max) - as.numeric(df_species$Trait_Temp_Min) ) / 2) + as.numeric(df_species$Trait_Temp_Min))
df_species$Time.mid<-as.character(((as.numeric(df_species$end_date_CE_final_num) - as.numeric(df_species$start_date_CE_final_num)) / 2) + as.numeric(df_species$start_date_CE_final_num))
#Add custom time ranges and regions
df_species$time_bins<-factor(cut(as.numeric(df_species$Time.mid), breaks=c(-800,601,901,1201,1501,2200), labels=c("<600","600-900","900-1200","1200-1500",">1500"),
                                 levels =c ("<600","600-900","900-1200","1200-1500",">1500")))
df_species$time_bins2<-factor(cut(as.numeric(df_species$Time.mid), breaks=c(-800,501,701,901,1101,1301,1501,1701,2200), 
                                  labels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700"),
                                  levels=c("<500","500-700","700-900","900-1100","1100-1300","1300-1500","1500-1700",">1700")))

names(df_species)
remove <- c(
  "aggregate","collapse_group","quality",
  "Assemblage_date_as_locally_defined","Start_date_CE","End_date_CE",
  "Unpublished_zooarchaeology_reference","Published_zooarchaeology_reference",
  "Archaeological_reference_for_chronology_etc_not_in_zooarch_refs","database","Taxa_cleaned",
  "Minimum_sieve_mesh","Data_quality","unidentified_fish_specimens","NISP_num","NSS_IDs_overlap","NSS_IDs_final",
  "NISP_total"
)
df_species<-df_species %>% dplyr::select(-any_of(remove))

# Add Cities Data ####
df_cities<-read_xlsx("../../../SupplementaryFiles/Supplementary_Danny's_Towns_Jan2025.xlsx")
names(df_cities)[4]<-"Country"
df_new_with_cities<-df_species%>% left_join(df_cities, by=c("Settlement_modern_name", "Country"), relationship="many-to-many")
df_new_with_cities$Time_range<-as.numeric(df_new_with_cities$end_date_CE_final_num) - as.numeric(df_new_with_cities$start_date_CE_final_num)
#Counts
df_new_with_cities$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_new_with_cities$Lumped_NSS_IDs))

#Rename columns
names(df_new_with_cities)<-gsub("Lumped_NSS_IDs","DB_Assemblage_ID",names(df_new_with_cities))
names(df_new_with_cities)<-gsub("region","Region",names(df_new_with_cities))
names(df_new_with_cities)<-gsub("Rural_or_urban","Rural_urban_or_neither",names(df_new_with_cities))
df_new_with_cities<-df_new_with_cities %>% dplyr::select(-`Original Spelling (Danny's table)`)
df_new_with_cities<-df_new_with_cities %>% dplyr::select(-Total_NISP)

df_new_with_cities <- df_new_with_cities %>%
  mutate(
    Time.mid = dplyr::if_else(
      !is.na(start_date_CE_final_num) & !is.na(end_date_CE_final_num),
      (start_date_CE_final_num + end_date_CE_final_num) / 2,
      as.numeric(NA)
    ),
    Time.range = dplyr::if_else(
      !is.na(start_date_CE_final_num) & !is.na(end_date_CE_final_num),
      end_date_CE_final_num - start_date_CE_final_num,
      as.numeric(NA)
    )
  )

df_new_with_cities <- df_new_with_cities %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(as.numeric(NISP), na.rm = TRUE)) %>%
  ungroup()

df_new_with_cities %>% 
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(DB_Assemblage_ID)
  )
#Export species-level dataset for sesnsitivty analysis - with PE ####
write_xlsx(df_new_with_cities, "NSS_SpeciesData_Aug2025_withPE.xlsx", col_names = TRUE, format_headers = TRUE)

df_new_with_cities_b<-df_new_with_cities[!df_new_with_cities$Region=="Poland & Estonia",]

df_new_with_cities_b %>% 
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(DB_Assemblage_ID)
  )

df_new_with_cities_b %>% 
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

df_new_with_cities_b <- df_new_with_cities_b %>%
  dplyr::mutate(NISP = as.numeric(NISP)) %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(total_NISP = sum(NISP, na.rm = TRUE)) %>%
  ungroup()

rm_cols<-c("Time.range","NISP_total")
df_new_with_cities_b<-df_new_with_cities_b %>% dplyr::select(-all_of(rm_cols))

#Remove small and large assemblages ####
df_new_with_cities_b %>%
  dplyr::group_by(Region) %>%
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

# 1) totals per assemblage (one row per DB_Assemblage_ID)

assemblage_tmp <- df_new_with_cities_b %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::summarise(total_NISP_assemblage = sum(NISP, na.rm = TRUE), .groups = "drop")

qs <- stats::quantile(assemblage_tmp$total_NISP_assemblage, c(0.01, 0.99), na.rm = TRUE)

ids <- assemblage_tmp %>%
  dplyr::filter(total_NISP_assemblage < qs[2]) %>%
  dplyr::pull(DB_Assemblage_ID)

df_new_with_cities_b <- df_new_with_cities_b %>%
  dplyr::filter(DB_Assemblage_ID %in% ids)

df_new_with_cities_b %>%
  dplyr::group_by(Region) %>%
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

rm(assemblage_tmp)

#Export species-level dataset for analysis - without PE ####
write_xlsx(df_new_with_cities_b, "NSS_SpeciesData_Sept2025_noPE.xlsx", col_names = TRUE, format_headers = TRUE)

#all good quality data
df_new_with_cities_b %>% 
  dplyr::summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(DB_Assemblage_ID)
  )

#NISP = 860,520
#N = 1292

#species level only
df_new_with_cities_b %>%
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))
#NISP = 625,072
#N = 1445

#species level only (Assemblages < 5 NISP removed)
df_supp<-df_new_with_cities_b %>%
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

df_supp<-df_new_with_cities_b %>%
  dplyr::group_by(DB_Assemblage_ID) %>%
  dplyr::mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%  # sum per ID
  dplyr::filter(!NISP_total < 5,
                !GBIF_species == "NA") %>%                           
  ungroup()

df_supp %>%
  dplyr::filter(!GBIF_species == "NA") %>%
  dplyr::summarize(NISP = sum(NISP),
                   Assemblage_count = n_distinct(DB_Assemblage_ID))

#NISP = 624,799
#N = 1309

