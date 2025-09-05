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
source("ThemePublish.R")
source("~/Documents/2.Academic_Work/3.NSS/2.MANUSCRIPT_Jan2025/Code_final/MANUSCRIPT_functions.R")

#2. load data ####
df<-read_xlsx("All_Data_with_Chronology_Aug2025.xlsx", col_names = TRUE, col_types = "text")
names(df)<-gsub(" ","_",names(df))
df<-df %>% dplyr::select(-any_of(c("Site_code","Year_reported")))
#Remove those assemblages with a NISP of less than 1
df<-df[!df$NISP < 0.99,]
#Convert georeferences to three decimal places
df$Decimal_Latitude<-round(as.numeric(df$Decimal_Latitude),3)/1000
df$Decimal_Longitude<-round(as.numeric(df$Decimal_Longitude),3)/1000

#3. count number of assemblages ####
df %>% 
  summarise(UniqueCount = n_distinct(NSS_unique_ID))

# Count total NISP
df$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df_merge<-df

df_merge$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df_old_short<-df_merge %>% dplyr::select(DB_ID3, Number_of_unidentified_fish_specimens) %>% distinct()
names(df_old_short)[2]<-"unidentified_fish_specimens"
df_old_short$unidentified_fish_specimens<-as.character(df_old_short$unidentified_fish_specimens)
df_old_short<-aggregate(unidentified_fish_specimens ~ DB_ID3, data = df_old_short, paste, collapse = ",")
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#Totals
df_old_short<-df_merge %>% dplyr::select(DB_ID3, Total_Fish_NISP) %>% distinct()
names(df_old_short)[2]<-"Total_NISP"
df_old_short$Total_NISP<-as.character(df_old_short$Total_NISP)
df_old_short<-aggregate(Total_NISP ~ DB_ID3, data = df_old_short, paste, collapse = ",")
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#Local site types
df_old_short<-df_merge %>% dplyr::select(DB_ID3, Site_type_using_local_categories) %>% distinct()
names(df_old_short)[2]<-"Local_site_type"
df_old_short$Local_site_type<-as.character(df_old_short$Local_site_type)
df_old_short<-aggregate(Local_site_type ~ DB_ID3, data = df_old_short, paste, collapse = ",")
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#Rural_urban_or_neither
df_old_short<-df_merge %>% dplyr::select(DB_ID3, Rural_urban_or_neither) %>% distinct()
names(df_old_short)[2]<-"Rural_or_urban"
df_old_short$Rural_or_urban<-as.character(df_old_short$Rural_or_urban)
df_old_short<-aggregate(Rural_or_urban ~ DB_ID3, data = df_old_short, paste, collapse = ",")
df_merge<- df_merge %>% merge(df_old_short, all.x=T)

#5. Convert numeric variables to text for aggregations ####
df_merge <-df_merge %>% mutate_at(c('Start_date_CE', 'End_date_CE',
                        'Decimal_Latitude','Decimal_Longitude', 
                        'ID_for_linked_rows_differing_only_by_recovery', 
                        'NISP','Trait_Temp_Max','Trait_Temp_Min'), as.character)
#Remove special characters from column names
names(df_merge)<-gsub(" ","_", names(df_merge))
names(df_merge)<-gsub("\\(","", names(df_merge))
names(df_merge)<-gsub("\\)","", names(df_merge))

#Remove unwanted columns for stacking
remove_cols<-c("DB_ID2","DB_ID4","Assemblage_ID","NSS_unique_ID","Input_by","ID_number_within_dataset",
               "File_name","Assemblage_or_sub-assemblage","County_province_or_state","Georef_source",
               "Context_types","ID_for_linked_rows_differing_only_by_recovery","Part_of_sieve_stack",
               "Recovery_method_hand-collected;_sieved;_both;_unknown","Sieve_sizes_and_recovery_details_with_original_mesh_units_if_not_mm",
               "Minimum_sieve_mesh_mm","Maximum_sieve_mesh_mm","Lumped_by_phase_or_split_by_context",
               "Unpublished_zooarchaeology_reference_J_Arch_Sci_format","Published_zooarchaeology_reference_J_Arch_Sci_format",
               "Archaeological_reference_for_chronology_etc_not_in_zooarch_refs_J_Arch_Sci_format","Analyst_name",
               "Dataset_contact_name","IP_status","Skeletal_element_data_available_y_or_n","Measurement_or_fish_size_data_available_y_or_n",
               "Data_quality_green_amber_or_red","Explanation_of_amber_or_red_data_quality_noting_what_variables_and_why",
               "General_comments","Taxa_with_only_presence_data","Phase_or_equivalent_optional","Context_or_equivalent_optional",
               "Sample_optional","Sediment_volumes_available_y_or_n_optional","Associated_mammal_NISP_available_y;_n;_not_matching_or_enter_value_if_easily_to_hand_optional",
               "Associated_bird_NISP_available_y;_n;_not_matching_or_enter_value_if_easily_to_hand_optional","Extra_column_for_dataset-specific_variable_optional"
               ,"Extra_column_for_dataset-specific_variable_optional","Definition_of_extra_column_variable_optional",
               "Number_of_unidentified_fish_specimens","Total_Fish_NISP","Rural_urban_or_neither","Site_type_using_local_categories")
df_merge<-df_merge %>% dplyr::select(-any_of(remove_cols))

df_merge_new$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df_merge$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df_merge$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df_out <- df_merge %>%
  mutate(
    start_date_CE = as.numeric(Start_date_CE),
    end_date_CE   = as.numeric(End_date_CE),
    minimum_sieve_size_num = readr::parse_number(Minimum_sieve_mesh)
  ) %>%
  group_by(
    Settlement_modern_name, Site_name, Decimal_Latitude, Decimal_Longitude,
    Country, Local_site_type, Rural_or_urban, Recovery, DB_sieved
  ) %>%
  mutate(
    .grp_min_start = if (all(is.na(start_date_CE))) NA_real_ else min(start_date_CE, na.rm = TRUE),
    .grp_max_end   = if (all(is.na(end_date_CE)))   NA_real_ else max(end_date_CE,   na.rm = TRUE),
    .total_span = .grp_max_end - .grp_min_start,
    .span_ok    = !is.na(.total_span) & .total_span < 200,
    .overlap_start = if (all(is.na(start_date_CE))) NA_real_ else max(start_date_CE, na.rm = TRUE),
    .overlap_end   = if (all(is.na(end_date_CE)))   NA_real_ else min(end_date_CE,   na.rm = TRUE),
    .has_overlap   = !is.na(.overlap_start) & !is.na(.overlap_end) & (.overlap_start <= .overlap_end),
    aggregate = .has_overlap & .span_ok,
    min_sieve_size_group = if (all(is.na(minimum_sieve_size_num))) NA_real_
    else min(minimum_sieve_size_num, na.rm = TRUE),
    lowest_start_overlap = ifelse(aggregate, .overlap_start, NA_real_),
    highest_end_overlap  = ifelse(aggregate, .overlap_end,   NA_real_),
    .in_overlap = aggregate &
      !is.na(start_date_CE) & !is.na(end_date_CE) &
      (start_date_CE <= .overlap_end) & (end_date_CE >= .overlap_start)
  ) %>%
  mutate(
    NSS_IDs_overlap = ifelse(
      aggregate,
      {
        ids <- unique(na.omit(Lumped_NSS_IDs[.in_overlap]))
        if (length(ids) == 0) NA_character_ else paste(ids, collapse = ";")
      },
      NA_character_
    )
  ) %>%
  mutate(
    start_date_CE_final_num = ifelse(aggregate, .grp_min_start, start_date_CE),
    end_date_CE_final_num   = ifelse(aggregate, .grp_max_end,   end_date_CE),
    min_sieve_size_final_num = ifelse(aggregate, min_sieve_size_group, minimum_sieve_size_num),
    NSS_IDs_final = ifelse(
      aggregate,
      NSS_IDs_overlap,
      as.character(Lumped_NSS_IDs)
    )
  ) %>%
  select(
    Start_date_CE, End_date_CE, Minimum_sieve_mesh,
    start_date_CE, end_date_CE, minimum_sieve_size_num,
    everything(),
    start_date_CE_final_num, end_date_CE_final_num, min_sieve_size_final_num, NSS_IDs_final,
    min_sieve_size_group, lowest_start_overlap, highest_end_overlap, aggregate
  ) %>%
  select(
    -starts_with(".grp_"),
    -.total_span, -.span_ok,
    -starts_with(".overlap"),
    -.has_overlap, -.in_overlap
  ) %>%
  ungroup()

remove<-c("DB_ID4","min_sieve_size_group","lowest_start_overlap","highest_end_overlap","NSS_IDs_overlap","DB_ID4_Group","DB_ID3","start_date_CE","end_date_CE","minimum_sieve_size_num")
df_out<-df_out %>% dplyr::select(-any_of(remove))
#Convert all columns except NISP to character
df_out <- df_out %>%
  mutate(across(-all_of("NISP"), as.character)) %>%
  # Use one of the following lines for NISP:
  mutate(NISP = as.numeric(NISP)) 
names(df_out)
# ---- Count TRUE/FALSE in aggregate ----
table(df_out$aggregate)

# Columns to form the identifier (deduped)
id_cols <- unique(c(
  "Settlement_modern_name","Site_name","Decimal_Latitude","Decimal_Longitude",
  "Country","Local_site_type","Rural_or_urban","Recovery","DB_sieved",
  "aggregate","start_date_CE_final_num","end_date_CE_final_num","min_sieve_size_final_num",
  "GBIF_species","GBIF_genus","GBIF_family","GBIF_order",
  "GBIF_chrondrichthyes_superorder","GBIF_chondricthyes_infraclass","GBIF__chondricthyes_subclass",
  "GBIF_class","GBIF_phylum","GBIF_kingdom","GBIF_level","NSS_IDs_final"
))

# 1) Build a stable identifier (without altering originals)
id_tmp <- df_out %>%
  transmute(across(all_of(id_cols), to_chr)) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), "(NA)", .))) %>%
  unite("identifier", all_of(id_cols), sep = ";", remove = TRUE)

df_with_id <- df_out %>%
  mutate(identifier = id_tmp$identifier)

# Column order for the final output: identifier + NISP + all original columns
final_cols <- c("identifier", "NISP", names(df_out))

# 2) Aggregate ONLY rows where aggregate == TRUE
df_agg <- df_with_id %>%
  filter(aggregate == TRUE) %>%
  group_by(identifier) %>%
  summarise(
    # Sum NISP (robust to any stray text like "12 pcs")
    NISP = sum(readr::parse_number(as.character(NISP)), na.rm = TRUE),
    # Paste ALL other columns as strings so nothing is dropped
    across(.cols = setdiff(names(df_with_id), c("identifier", "NISP")),
           .fns  = collapse_col,
           .names = "{.col}"),
    .groups = "drop"
  ) %>%
  # make sure columns are in the same order as requested
  select(all_of(final_cols))

# 3) Keep NON-aggregating rows as rows (no grouping), but cast non-NISP to character for consistency
df_nonagg <- df_with_id %>%
  filter(is.na(aggregate) | aggregate == FALSE) %>%
  mutate(
    NISP = readr::parse_number(as.character(NISP)),
    across(.cols = setdiff(names(.), c("identifier","NISP")), .fns = to_chr)
  ) %>%
  select(all_of(final_cols))

# 4) Combine
df_result <- bind_rows(df_agg, df_nonagg)

# df_result now has:
# 1. one row per aggregated identifier for aggregate==TRUE groups (NISP summed, others pasted),
# 2. original rows for aggregate==FALSE/NA,
# 3. all original columns preserved or lumped if differences in ro, plus 'identifier'.
df_result$NISP %>% as.numeric() %>%
  sum(na.rm=T)

df_result$NISP %>% as.numeric() %>%
  sum(na.rm=T)

length(unique(df_result$Lumped_NSS_IDs))

# 8. Assemblages that are chronological outliers ####
df_result<-df_result[df_result$start_date_CE_final_num |> as.numeric() >= 1,]
df_result<-df_result[df_result$start_date_CE_final_num |> as.numeric() <= 1850,]
df_result<-df_result[df_result$end_date_CE_final_num |> as.numeric() <= 1900,]
ggplot(df_result, aes(as.numeric(start_date_CE_final_num))) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "white") +
  scale_x_continuous(limits = c(1, 1901), breaks = seq(0, 2000, 500)) +
  labs(title = "Histogram with x-axis 1–2000",
       x = "Values",
       y = "Count") +
  theme_publish()
ggplot(df_result, aes(as.numeric(end_date_CE_final_num))) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "white") +
  scale_x_continuous(limits = c(1, 1901), breaks = seq(0, 2000, 500)) +
  labs(title = "Histogram with x-axis 1–2000",
       x = "Values",
       y = "Count") +
  theme_publish()

df_result$NISP %>% as.numeric() %>%
  sum(na.rm=T)

length(unique(df_result$Lumped_NSS_IDs))

#N=2109
#NISP=1805727

# Remove assemblages with low or v. high fish counts ####
#Recalculate total fish NISP
df_summed <- df_result %>%
  group_by(
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
  mutate(NISP_total = sum(NISP, na.rm = TRUE)) %>%
  ungroup()

#Remove small and large assemblages ####
df_tmp<-df_summed %>% filter(NISP_total < quantile(df_summed$NISP_total, 0.975, na.rm=T))
df_tmp<-df_tmp %>% filter(NISP_total > quantile(df_tmp$NISP_total, 0.025, na.rm=T))
df_tmp$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_tmp$Lumped_NSS_IDs))
#N=1751
#NISP=1308730

# Remove assemblages with chronological range of more than 400
df_tmp<- df_tmp %>%
  filter((end_date_CE_final_num |> as.numeric() - start_date_CE_final_num |> as.numeric()) <= 400) %>%
  filter(complete.cases(Decimal_Longitude, Decimal_Latitude))
df_tmp$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_tmp$Lumped_NSS_IDs))
#N=1370
#NISP=1010617

df_tmp$region<-recode(df_tmp$Country,
                          "England" = 'Britian & Ireland',
                          "Scotland" = 'Britian & Ireland',
                          "Norway" = 'Scandinavia',
                          "Sweden" = 'Scandinavia',
                          "Denmark" = 'Scandinavia',
                          "Germany" = 'Western Europe',
                          "Belgium" = 'Western Europe',
                          "Netherlands" = 'Western Europe',
                          "Ireland" = 'Britian & Ireland',
                          "Estonia" = 'Poland & Estonia',
                          "Poland" = 'Poland & Estonia',
                          "Northern Ireland" = 'Britian & Ireland')
df_tmp %>% 
  group_by(region) %>%
  summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(Lumped_NSS_IDs)
  )

df_tmp %>% 
  summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(Lumped_NSS_IDs)
  )

###Subset to species level data
df_tmp<-df_tmp %>% filter(
  !is.na(GBIF_species)
)

df_tmp<-df_tmp %>% filter(
  !GBIF_species == "NA"
)

#Remove these columns
rm_cols<-c("Original_name","Has_taxa","GBIF_level")
df_tmp<-df_tmp %>% dplyr::select(-all_of(rm_cols))

df_tmp$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_tmp$Lumped_NSS_IDs))
#N=1364
#NISP=691,514

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
df_species$region<-recode(df_species$Country,
                          "England" = 'Britian & Ireland',
                          "Scotland" = 'Britian & Ireland',
                          "Norway" = 'Scandinavia',
                          "Sweden" = 'Scandinavia',
                          "Denmark" = 'Scandinavia',
                          "Germany" = 'Western Europe',
                          "Belgium" = 'Western Europe',
                          "Netherlands" = 'Western Europe',
                          "Ireland" = 'Britian & Ireland',
                          "Estonia" = 'Poland & Estonia',
                          "Poland" = 'Poland & Estonia',
                          "Northern Ireland" = 'Britian & Ireland')
# Add Cities Data ####
df_cities<-read_xlsx("../../SupplementaryFiles/Supplementary_Danny's_Towns_Jan2025.xlsx")
names(df_cities)[4]<-"Country"
df_new_with_cities<-df_species%>% left_join(df_cities)
df_new_with_cities$Time_range<-as.numeric(df_new_with_cities$End_date_CE) - as.numeric(df_new_with_cities$Start_date_CE)
#Counts
df_new_with_cities$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_new_with_cities$Lumped_NSS_IDs))

#Rename columns
names(df_new_with_cities)<-gsub("DB_ID3","DB_Assemblage_ID",names(df_new_with_cities))
names(df_new_with_cities)<-gsub("region","Region",names(df_new_with_cities))
names(df_new_with_cities)<-gsub("Rural_or_urban","Rural_urban_or_neither",names(df_new_with_cities))
df_new_with_cities<-df_new_with_cities %>% dplyr::select(-`Original Spelling (Danny's table)`)

df_new_with_cities %>% 
  group_by(Region) %>%
  summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(Lumped_NSS_IDs)
  )
#Export species-level dataset for sesnsitivty analysis - with PE ####
write_xlsx(df_new_with_cities, "NSS_SpeciesData_Aug2025_withPE.xlsx", col_names = TRUE, format_headers = TRUE)

df_new_with_cities_b<-df_new_with_cities[!df_new_with_cities$Region=="Poland & Estonia",]
df_new_with_cities_b$NISP %>% as.numeric() %>%
  sum(na.rm=T)
length(unique(df_new_with_cities_b$Lumped_NSS_IDs))

df_new_with_cities_b %>% 
  group_by(Region) %>%
  summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(Lumped_NSS_IDs)
  )

df_new_with_cities_b %>% 
  summarise(
    Total_NISP = sum(NISP, na.rm=TRUE),
    Assemblage_count = n_distinct(Lumped_NSS_IDs)
  )

#Export species-level dataset for analysis - without PE ####
write_xlsx(df_new_with_cities_b, "NSS_SpeciesData_Aug2025_noPE.xlsx", col_names = TRUE, format_headers = TRUE)

