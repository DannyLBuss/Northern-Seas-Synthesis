# NSS Script 1: data merge and harmonisation
# Author: Dr Danny L Buss
# Date: 2025-02-05

# 0. Load Packages ####
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(writexl)
  library(textclean)
  library(here)
  library(reshape2)
})

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../ThemePublish.R")

#3. load data
# Belgium ####
Belgium<-read_xlsx("../NSS_DBs_VORs/Belgium_VoR_14Aug2025.xlsx", col_names = TRUE, col_types = "text")
Belgium$`NSS unique ID`<-as.character(Belgium$`NSS unique ID`)
Belgium<- Belgium %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Belgium<-Belgium %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Belgium)<-gsub("¬†"," ",names(Belgium))
names(Belgium)<-trimws(iconv(names(Belgium), "UTF-8", "ASCII", sub = " "))
Belgium <- Belgium %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Belgium$Database<-"Belgium_VoR_14Aug2025"
names(Belgium)
Belgium<-Belgium[,c(1:45,184,46:183)]
length(unique(names(Belgium[46:183])))
length(unique(Belgium$`Site name`))
nrow(Belgium)
length(unique(Belgium$`NSS unique ID`))
Belgium$`Total Fish NISP`<-stringr::str_replace(Belgium$`Total Fish NISP`,"\\.0","")
Belgium$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Belgium$`NSS unique ID`<-stringr::str_replace(Belgium$`NSS unique ID`,"\\.0","")

# Denmark ####
Denmark<-read_xlsx("../NSS_DBs_VORs/Denmark_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Denmark<- Denmark %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Denmark<-Denmark %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Denmark)<-gsub("¬†"," ",names(Denmark))
names(Denmark)<-trimws(iconv(names(Denmark), "UTF-8", "ASCII", sub = " "))
Denmark <- Denmark %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Denmark$Database<-"Denmark_VoR_14Aug2025"
names(Denmark)
Denmark<-Denmark[,c(1:45,132,46:131)]
length(unique(names(Denmark[46:131])))
length(unique(Denmark$`Site name`))
nrow(Denmark)
Denmark$`Total Fish NISP`<-stringr::str_replace(Denmark$`Total Fish NISP`,"\\.0","")
Denmark$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Denmark$`NSS unique ID`<-stringr::str_replace(Denmark$`NSS unique ID`,"\\.0","")

# England ####
England<-read_xlsx("../NSS_DBs_VORs/England_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text") 
England<- England %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
England<-England %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(England)<-gsub("¬†"," ",names(England))
names(England)<-trimws(iconv(names(England), "UTF-8", "ASCII", sub = " "))
England<-England %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
England<-England %>%
  dplyr::mutate(across(where(is.character), ~ gsub("  ", " ", .x)))
England$Database<-"England_VoR_14Aug2025"
England<-England[,c(1:45,680,46:679)]
length(unique(names(England[46:679])))
length(unique(England$`Site name`))
nrow(England)
England$`Total Fish NISP`<-stringr::str_replace(England$`Total Fish NISP`,"\\.0","")
England$`Total Fish NISP` %>% as.numeric() %>%
  sum()
England$`NSS unique ID`<-stringr::str_replace(England$`NSS unique ID`,"\\.0","")

# Estonia ####
Estonia<-read_xlsx("../NSS_DBs_VORs/Estonia_VoR_13Aug2025.xlsx", col_names = TRUE,  col_types = "text") 
Estonia<- Estonia %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Estonia<-Estonia %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Estonia)<-gsub("¬†"," ",names(Estonia))
names(Estonia)<-trimws(iconv(names(Estonia), "UTF-8", "ASCII", sub = " "))
Estonia<-Estonia %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Estonia$Database<-"Estonia_VoR_13Aug2025"
Estonia<-Estonia[,c(1:45,76,46:75)]
length(unique(names(Estonia[46:75])))
length(unique(Estonia$`Site name`))
nrow(Estonia)
Estonia$`Total Fish NISP`<-stringr::str_replace(Estonia$`Total Fish NISP`,"\\.0","")
Estonia$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Estonia$`NSS unique ID`<-stringr::str_replace(Estonia$`NSS unique ID`,"\\.0","")

# Germany x 2 ####
Germany_Hol<-read_xlsx("../NSS_DBs_VORs/Germany_Holocene_DB_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Germany_Hol<- Germany_Hol %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Germany_Hol<-Germany_Hol %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Germany_Hol)<-gsub("¬†"," ",names(Germany_Hol))
names(Germany_Hol)<-trimws(iconv(names(Germany_Hol), "UTF-8", "ASCII", sub = " "))
Germany_Hol<-Germany_Hol %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Germany_Hol$Database<-"Germany_Holocene_DB_VoR_14Aug2025"
Germany_Hol<-Germany_Hol[,c(1:45,118,46:117)]
length(unique(names(Germany_Hol[46:117])))
length(unique(Germany_Hol$`Site name`))
nrow(Germany_Hol)
Germany_Hol<-Germany_Hol[complete.cases(Germany_Hol$`NSS unique ID`),]
Germany_Hol$`Total Fish NISP`<-stringr::str_replace(Germany_Hol$`Total Fish NISP`,"\\.0","")
Germany_Hol$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Germany_Hol$`NSS unique ID`<-stringr::str_replace(Germany_Hol$`NSS unique ID`,"\\.0","")

Germany_Hein<-read_xlsx("../NSS_DBs_VORs/Germany_Küchelmann_Heinrich_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Germany_Hein<- Germany_Hein %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Germany_Hein<-Germany_Hein %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Germany_Hein)<-gsub("¬†"," ",names(Germany_Hein))
names(Germany_Hein)<-trimws(iconv(names(Germany_Hein), "UTF-8", "ASCII", sub = " "))
Germany_Hein<-Germany_Hein %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Germany_Hein$Database<-"Germany_Küchelmann_Heinrich_VoR_14Aug2025"
Germany_Hein<-Germany_Hein[,c(1:45,139,46:138)]
length(unique(names(Germany_Hein[46:138])))
length(unique(Germany_Hein$`Site name`))
nrow(Germany_Hein)
Germany_Hein$`Total Fish NISP`<-stringr::str_replace(Germany_Hein$`Total Fish NISP`,"\\.0","")
Germany_Hein$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Germany_Hein$`NSS unique ID`<-stringr::str_replace(Germany_Hein$`NSS unique ID`,"\\.0","")

# Netherlands x 2 ####
Netherlands_RomSup<-read_xlsx("../NSS_DBs_VORs/Netherlands_Rom_suppl_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Netherlands_RomSup<- Netherlands_RomSup %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Netherlands_RomSup<-Netherlands_RomSup %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Netherlands_RomSup)<-gsub("¬†"," ",names(Netherlands_RomSup))
names(Netherlands_RomSup)<-trimws(iconv(names(Netherlands_RomSup), "UTF-8", "ASCII", sub = " "))
Netherlands_RomSup<-Netherlands_RomSup %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Netherlands_RomSup$Database<-"Netherlands_Rom_suppl_VoR_14Aug2025"
Netherlands_RomSup<-Netherlands_RomSup[,c(1:45,159,46:158)]
length(unique(names(Netherlands_RomSup[46:158])))
length(unique(Netherlands_RomSup$`Site name`))
nrow(Netherlands_RomSup)
Netherlands_RomSup$`Total Fish NISP`<-stringr::str_replace(Netherlands_RomSup$`Total Fish NISP`,"\\.0","")
Netherlands_RomSup$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Netherlands_RomSup$`NSS unique ID`<-stringr::str_replace(Netherlands_RomSup$`NSS unique ID`,"\\.0","")

Netherlands_Rom<-read_xlsx("../NSS_DBs_VORs/Netherlands_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Netherlands_Rom<- Netherlands_Rom %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Netherlands_Rom<-Netherlands_Rom %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Netherlands_Rom)<-gsub("¬†"," ",names(Netherlands_Rom))
names(Netherlands_Rom)<-trimws(iconv(names(Netherlands_Rom), "UTF-8", "ASCII", sub = " "))
Netherlands_Rom<-Netherlands_Rom %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Netherlands_Rom$Database<-"Netherlands_VoR_14Aug2025"
Netherlands_Rom<-Netherlands_Rom[,c(1:45,191,46:190)]
length(unique(names(Netherlands_Rom[46:190])))
length(unique(Netherlands_Rom$`Site name`))
nrow(Netherlands_Rom)
Netherlands_Rom$`Total Fish NISP`<-stringr::str_replace(Netherlands_Rom$`Total Fish NISP`,"\\.0","")
Netherlands_Rom$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Netherlands_Rom$`NSS unique ID`<-stringr::str_replace(Netherlands_Rom$`NSS unique ID`,"\\.0","")

# Ireland ####
Ireland<-read_xlsx("../NSS_DBs_VORs/Ireland_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Ireland<- Ireland %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Ireland<-Ireland %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Ireland)<-gsub("¬†"," ",names(Ireland))
names(Ireland)<-trimws(iconv(names(Ireland), "UTF-8", "ASCII", sub = " "))
Ireland<-Ireland %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Ireland$Database<-"Ireland_VoR_14Aug2025"
Ireland<-Ireland[,c(1:45,86,46:85)]
length(unique(names(Ireland[46:85])))
length(unique(Ireland$`Site name`))
nrow(Ireland)
Ireland$`Total Fish NISP`<-stringr::str_replace(Ireland$`Total Fish NISP`,"\\.0","")
Ireland$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Ireland$`NSS unique ID`<-stringr::str_replace(Ireland$`NSS unique ID`,"\\.0","")

# Norway ####
Norway<-read_xlsx("../NSS_DBs_VORs/Norway_VoR_13Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Norway<- Norway %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Norway<-Norway %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Norway)<-gsub("¬†"," ",names(Norway))
names(Norway)<-trimws(iconv(names(Norway), "UTF-8", "ASCII", sub = " "))
Norway<-Norway %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Norway$Database<-"Norway_VoR_13Aug2025"
Norway<-Norway[,c(1:45,122,46:121)]
length(unique(names(Norway[46:121])))
length(unique(Norway$`Site name`))
nrow(Norway)
Norway$`Total Fish NISP`<-stringr::str_replace(Norway$`Total Fish NISP`,"\\.0","")
Norway$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Norway$`NSS unique ID`<-stringr::str_replace(Norway$`NSS unique ID`,"\\.0","")

# Poland ####
Poland<-read_xlsx("../NSS_DBs_VORs/Poland_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Poland<- Poland %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Poland<-Poland %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Poland)<-gsub("¬†"," ",names(Poland))
names(Poland)<-trimws(iconv(names(Poland), "UTF-8", "ASCII", sub = " "))
Poland<-Poland %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Poland<-Poland %>%
  dplyr::mutate(across(where(is.character), ~ gsub("  ", " ", .x)))
Poland$Database<-"Poland_VoR_14Aug2025"
Poland<-Poland[,c(1:45,96,46:95)]
length(unique(names(Poland[46:95])))
length(unique(Poland$`Site name`))
nrow(Poland)
Poland$`Total Fish NISP`<-stringr::str_replace(Poland$`Total Fish NISP`,"\\.0","")
Poland$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Poland$`NSS unique ID`<-stringr::str_replace(Poland$`NSS unique ID`,"\\.0","")

# Scotland x 2 ####
Scotland<-read_xlsx("../NSS_DBs_VORs/Scotland_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Scotland<- Scotland %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Scotland<-Scotland %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Scotland)<-gsub("¬†"," ",names(Scotland))
names(Scotland)<-trimws(iconv(names(Scotland), "UTF-8", "ASCII", sub = " "))
Scotland<-Scotland %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Scotland<-Scotland %>%
  dplyr::mutate(across(where(is.character), ~ gsub("  ", " ", .x)))
Scotland$Database<-"Scotland_VoR_14Aug2025"
Scotland<-Scotland[,c(1:45,165,46:164)]
length(unique(names(Scotland[46:164])))
length(unique(Scotland$`Site name`))
Scotland<-Scotland[complete.cases(Scotland$`NSS unique ID`),]
nrow(Scotland)
Scotland$`Total Fish NISP`<-stringr::str_replace(Scotland$`Total Fish NISP`,"\\.0","")
Scotland$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Scotland$`NSS unique ID`<-stringr::str_replace(Scotland$`NSS unique ID`,"\\.0","")

Scotland_b<-read_xlsx("../NSS_DBs_VORs/Scotland_suppl_VoR_14Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Scotland_b<- Scotland_b %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Scotland_b<-Scotland_b %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Scotland_b)<-gsub("¬†"," ",names(Scotland_b))
names(Scotland_b)<-trimws(iconv(names(Scotland_b), "UTF-8", "ASCII", sub = " "))
Scotland_b<-Scotland_b %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Scotland_b<-Scotland_b %>%
  dplyr::mutate(across(where(is.character), ~ gsub("  ", " ", .x)))
Scotland_b$Database<-"Scotland_suppl_VoR_14Aug2025"
Scotland_b<-Scotland_b[,c(1:45,69,46:68)]
length(unique(names(Scotland_b[46:68])))
length(unique(Scotland_b$`Site name`))
nrow(Scotland_b)
Scotland_b$`Total Fish NISP`<-stringr::str_replace(Scotland_b$`Total Fish NISP`,"\\.0","")
Scotland_b$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Scotland_b$`NSS unique ID`<-stringr::str_replace(Scotland_b$`NSS unique ID`,"\\.0","")

# Sweden ####
Sweden<-read_xlsx("../NSS_DBs_VORs/Sweden_VoR_13Aug2025.xlsx", col_names = TRUE,  col_types = "text")
Sweden<-Sweden %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
Sweden<-Sweden %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(Sweden)<-gsub("¬†"," ",names(Sweden))
names(Sweden)<-trimws(iconv(names(Sweden), "UTF-8", "ASCII", sub = " "))
Sweden<-Sweden %>%
  dplyr::mutate(across(where(is.character), ~ trimws(iconv(.x, "UTF-8", "ASCII", sub = " "))))
Sweden$Database<-"Sweden_VoR_13Aug2025"
Sweden<-Sweden[,c(1:45,124,46:123)]
length(unique(names(Sweden[46:123])))
length(unique(Sweden$`Site name`))
nrow(Sweden)
Sweden$`Total Fish NISP`<-stringr::str_replace(Sweden$`Total Fish NISP`,"\\.0","")
Sweden$`Total Fish NISP` %>% as.numeric() %>%
  sum()
Sweden$`NSS unique ID`<-stringr::str_replace(Sweden$`NSS unique ID`,"\\.0","")

# Load Data for Taxonomic Harmonisation ####
META<-readxl::read_xlsx("S4_Taxon_Traits_Table.xlsx", col_names = TRUE)
str(META)
taxa_meta<-META
taxa_meta<-taxa_meta %>%
  dplyr::mutate(across(everything(), ~ trimws(.)))
taxa_meta<-taxa_meta %>%
  dplyr::mutate(across(everything(), ~ textclean::replace_non_ascii(.)))
names(taxa_meta)<-gsub("¬†"," ",names(taxa_meta))
taxa_meta<-taxa_meta %>%
  dplyr::mutate(across(everything(), ~ gsub("¬†"," ",.)))

#Merge datasets
taxa<-data.frame("Original_name" = c(names(Germany_Hein)[47:ncol(Germany_Hein)],
                                     names(Germany_Hol)[47:ncol(Germany_Hol)],
                                     names(Belgium)[47:ncol(Belgium)],
                                     names(Denmark)[47:ncol(Denmark)],
                                     names(England)[47:ncol(England)],
                                     names(Estonia)[47:ncol(Estonia)],
                                     names(Poland)[47:ncol(Poland)],
                                     names(Netherlands_Rom)[47:ncol(Netherlands_Rom)],
                                     names(Netherlands_RomSup)[47:ncol(Netherlands_RomSup)],
                                     names(Ireland)[47:ncol(Ireland)],
                                     names(Norway)[47:ncol(Norway)],
                                     names(Scotland)[47:ncol(Scotland)],
                                     names(Scotland_b)[47:ncol(Scotland_b)],
                                     names(Sweden)[47:ncol(Sweden)]),
                 "Original_DB" = c(rep("Germany_Küchelmann_Heinrich_VoR_14Aug2025", length(names(Germany_Hein)[47:ncol(Germany_Hein)])),
                                   rep("Germany_Holocene_DB_VoR_14Aug2025", length(names(Germany_Hol)[47:ncol(Germany_Hol)])),
                                   rep("Belgium_VoR_14Aug2025", length(names(Belgium)[47:ncol(Belgium)])),
                                   rep("Denmark_VoR_14Aug2025", length(names(Denmark)[47:ncol(Denmark)])),
                                   rep("England_VoR_14Aug2025", length(names(England)[47:ncol(England)])),
                                   rep("Estonia_VoR_13Aug2025", length(names(Estonia)[47:ncol(Estonia)])),
                                   rep("Poland_VoR_14Aug2025", length(names(Poland)[47:ncol(Poland)])),
                                   rep("Netherlands_VoR_14Aug2025", length(names(Netherlands_Rom)[47:ncol(Netherlands_Rom)])),
                                   rep("Netherlands_Rom_suppl_VoR_14Aug2025",length(names(Netherlands_RomSup)[47:ncol(Netherlands_RomSup)])),
                                   rep("Ireland_VoR_14Aug2025",length(names(Ireland)[47:ncol(Ireland)])),
                                   rep("Norway_VoR_13Aug2025", length(names(Norway)[47:ncol(Norway)])),
                                   rep("Scotland_VoR_14Aug2025",length(names(Scotland)[47:ncol(Scotland)])),
                                   rep("Scotland_suppl_VoR_14Aug2025",length(names(Scotland_b)[47:ncol(Scotland_b)])),
                                   rep("Sweden_VoR_13Aug2025", length(names(Sweden)[47:ncol(Sweden)]))))
taxa$Original_name<-trimws(taxa$Original_name)
taxa_meta$Original_name<-trimws(taxa_meta$Original_name)
taxa$Original_DB<-as.factor(taxa$Original_DB)
taxa$Original_name<-gsub("¬†"," ",taxa$Original_name)
taxa$Original_name<-gsub("\\¬†"," ",taxa$Original_name)
taxa$Original_name<-gsub("\\†"," ",taxa$Original_name)
taxa$Original_name<-gsub("\\¬"," ",taxa$Original_name)
taxa$Original_name<-trimws(iconv(taxa$Original_name, "UTF-8", "ASCII", sub = " "))
taxa$Original_name<-gsub("  "," ",taxa$Original_name)
taxa_meta$Original_name<-trimws(iconv(taxa_meta$Original_name, "UTF-8", "ASCII", sub = " "))
all_taxa_headers<-unique(as.character(taxa$Original_name))
meta<-unique(as.character(taxa_meta$Original_name))
all_taxa_headers[!all_taxa_headers%in%meta]
x<-meta[!meta%in%all_taxa_headers]

#All taxa across the databases are found within the taxa metadata column, and no additional taxa were present

#Convert to long format ####
str((Belgium)[47:ncol(Belgium)])
Belgium[47:ncol(Belgium)] <- sapply(Belgium[47:ncol(Belgium)],as.numeric)
Belgium$`NSS unique ID` <- as.character(Belgium$`NSS unique ID`)
Belgium_Meta<-Belgium[,c(1:46)]
Belgium_Taxa<-melt(Belgium, id.vars = 1:46)
names(Belgium_Taxa)[47]<-"Original_name"
names(Belgium_Taxa)[48]<-"NISP"
Belgium_Taxa$Original_name<-trimws(Belgium_Taxa$Original_name)
Belgium_Taxa<-Belgium_Taxa[,c(1,47,48)]
Belgium_Taxa<-Belgium_Taxa[!is.na(Belgium_Taxa$NISP),]
Belgium_Taxa$database<-"Belgium_VoR_14Aug2025"
Belgium_Meta<-Belgium_Meta %>%
  dplyr::mutate_all(as.character)
Belgium_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Belgium_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Denmark)[47:ncol(Denmark)])
Denmark[47:ncol(Denmark)] <- sapply(Denmark[47:ncol(Denmark)],as.numeric)
Denmark$`NSS unique ID` <- as.character(Denmark$`NSS unique ID`)
Denmark_Meta<-Denmark[,c(1:46)]
Denmark_Taxa<-melt(Denmark, id.vars = 1:46)
names(Denmark_Taxa)[47]<-"Original_name"
names(Denmark_Taxa)[48]<-"NISP"
Denmark_Taxa$Original_name<-trimws(Denmark_Taxa$Original_name)
Denmark_Taxa<-Denmark_Taxa[,c(1,47,48)]
Denmark_Taxa<-Denmark_Taxa[!is.na(Denmark_Taxa$NISP),]
Denmark_Taxa$database<-"Denmark_VoR_14Aug2025"
Denmark_Meta<-Denmark_Meta %>%
  dplyr::mutate_all(as.character)
names(Denmark_Meta)<-names(Belgium_Meta)
Denmark_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Denmark_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((England)[47:ncol(England)])
England[47:ncol(England)] <- sapply(England[47:ncol(England)],as.numeric)
England$`NSS unique ID` <- as.character(England$`NSS unique ID`)
England_Meta<-England[,c(1:46)]
England_Taxa<-melt(England, id.vars = 1:46)
names(England_Taxa)[47]<-"Original_name"
names(England_Taxa)[48]<-"NISP"
England_Taxa$Original_name<-trimws(England_Taxa$Original_name)
England_Taxa<-England_Taxa[,c(1,47,48)]
England_Taxa<-England_Taxa[!is.na(England_Taxa$NISP),]
England_Taxa$database<-"England_VoR_14Aug2025"
England_Meta<-England_Meta %>%
  dplyr::mutate_all(as.character)
names(England_Meta)<-names(Belgium_Meta)
England_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
England_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Estonia)[47:ncol(Estonia)])
Estonia[47:ncol(Estonia)] <- sapply(Estonia[47:ncol(Estonia)],as.numeric)
Estonia$`NSS unique ID` <- as.character(Estonia$`NSS unique ID`)
Estonia_Meta<-Estonia[,c(1:46)]
Estonia_Taxa<-melt(Estonia, id.vars = 1:46)
names(Estonia_Taxa)[47]<-"Original_name"
names(Estonia_Taxa)[48]<-"NISP"
Estonia_Taxa$Original_name<-trimws(Estonia_Taxa$Original_name)
Estonia_Taxa<-Estonia_Taxa[,c(1,47,48)]
Estonia_Taxa<-Estonia_Taxa[!is.na(Estonia_Taxa$NISP),]
Estonia_Taxa$database<-"Estonia_VoR_13Aug2025"
Estonia_Meta<-Estonia_Meta %>%
  dplyr::mutate_all(as.character)
names(Estonia_Meta)<-names(Belgium_Meta)
Estonia_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Estonia_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Germany_Hein)[47:ncol(Germany_Hein)])
Germany_Hein[47:ncol(Germany_Hein)] <- sapply(Germany_Hein[47:ncol(Germany_Hein)],as.numeric)
Germany_Hein$`NSS unique ID` <- as.character(Germany_Hein$`NSS unique ID`)
Germany_Hein_Meta<-Germany_Hein[,c(1:46)]
Germany_Hein_Taxa<-melt(Germany_Hein, id.vars = 1:46)
names(Germany_Hein_Taxa)[47]<-"Original_name"
names(Germany_Hein_Taxa)[48]<-"NISP"
Germany_Hein_Taxa$Original_name<-trimws(Germany_Hein_Taxa$Original_name)
Germany_Hein_Taxa<-Germany_Hein_Taxa[,c(1,47,48)]
Germany_Hein_Taxa<-Germany_Hein_Taxa[!is.na(Germany_Hein_Taxa$NISP),]
Germany_Hein_Taxa$database<-"Germany_Hein_VoR_14Aug2025"
Germany_Hein_Meta<-Germany_Hein_Meta %>%
  dplyr::mutate_all(as.character)
names(Germany_Hein_Meta)<-names(Belgium_Meta)
Germany_Hein_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Germany_Hein_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Germany_Hol)[47:ncol(Germany_Hol)])
Germany_Hol[47:ncol(Germany_Hol)] <- sapply(Germany_Hol[47:ncol(Germany_Hol)],as.numeric)
Germany_Hol$`NSS unique ID` <- as.character(Germany_Hol$`NSS unique ID`)
Germany_Hol_Meta<-Germany_Hol[,c(1:46)]
Germany_Hol_Taxa<-melt(Germany_Hol, id.vars = 1:46)
names(Germany_Hol_Taxa)[47]<-"Original_name"
names(Germany_Hol_Taxa)[48]<-"NISP"
Germany_Hol_Taxa$Original_name<-trimws(Germany_Hol_Taxa$Original_name)
Germany_Hol_Taxa<-Germany_Hol_Taxa[,c(1,47,48)]
Germany_Hol_Taxa<-Germany_Hol_Taxa[!is.na(Germany_Hol_Taxa$NISP),]
Germany_Hol_Taxa$database<-"Germany_Hol_VoR_14Aug2025"
Germany_Hol_Meta<-Germany_Hol_Meta %>%
  dplyr::mutate_all(as.character)
names(Germany_Hol_Meta)<-names(Belgium_Meta)
Germany_Hol_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Germany_Hol_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Netherlands_Rom)[47:ncol(Netherlands_Rom)])
Netherlands_Rom[47:ncol(Netherlands_Rom)] <- sapply(Netherlands_Rom[47:ncol(Netherlands_Rom)],as.numeric)
Netherlands_Rom$`NSS unique ID` <- as.character(Netherlands_Rom$`NSS unique ID`)
Netherlands_Rom_Meta<-Netherlands_Rom[,c(1:46)]
Netherlands_Rom_Taxa<-melt(Netherlands_Rom, id.vars = 1:46)
names(Netherlands_Rom_Taxa)[47]<-"Original_name"
names(Netherlands_Rom_Taxa)[48]<-"NISP"
Netherlands_Rom_Taxa$Original_name<-trimws(Netherlands_Rom_Taxa$Original_name)
Netherlands_Rom_Taxa<-Netherlands_Rom_Taxa[,c(1,47,48)]
Netherlands_Rom_Taxa<-Netherlands_Rom_Taxa[!is.na(Netherlands_Rom_Taxa$NISP),]
Netherlands_Rom_Taxa$database<-"Netherlands_Rom_VoR_14Aug2025"
Netherlands_Rom_Meta<-Netherlands_Rom_Meta %>%
  dplyr::mutate_all(as.character)
names(Netherlands_Rom_Meta)<-names(Belgium_Meta)
Netherlands_Rom_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Netherlands_Rom_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Netherlands_RomSup)[47:ncol(Netherlands_RomSup)])
Netherlands_RomSup[47:ncol(Netherlands_RomSup)] <- sapply(Netherlands_RomSup[47:ncol(Netherlands_RomSup)],as.numeric)
Netherlands_RomSup$`NSS unique ID` <- as.character(Netherlands_RomSup$`NSS unique ID`)
Netherlands_RomSup_Meta<-Netherlands_RomSup[,c(1:46)]
Netherlands_RomSup_Taxa<-melt(Netherlands_RomSup, id.vars = 1:46)
names(Netherlands_RomSup_Taxa)[47]<-"Original_name"
names(Netherlands_RomSup_Taxa)[48]<-"NISP"
Netherlands_RomSup_Taxa$Original_name<-trimws(Netherlands_RomSup_Taxa$Original_name)
Netherlands_RomSup_Taxa<-Netherlands_RomSup_Taxa[,c(1,47,48)]
Netherlands_RomSup_Taxa<-Netherlands_RomSup_Taxa[!is.na(Netherlands_RomSup_Taxa$NISP),]
Netherlands_RomSup_Taxa$database<-"Netherlands_RomSup_VoR_14Aug2025"
Netherlands_RomSup_Meta<-Netherlands_RomSup_Meta %>%
  dplyr::mutate_all(as.character)
names(Netherlands_RomSup_Meta)<-names(Belgium_Meta)
Netherlands_RomSup_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Netherlands_RomSup_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Norway)[47:ncol(Norway)])
Norway[47:ncol(Norway)] <- sapply(Norway[47:ncol(Norway)],as.numeric)
Norway$`NSS unique ID` <- as.character(Norway$`NSS unique ID`)
Norway_Meta<-Norway[,c(1:46)]
Norway_Taxa<-melt(Norway, id.vars = 1:46)
names(Norway_Taxa)[47]<-"Original_name"
names(Norway_Taxa)[48]<-"NISP"
Norway_Taxa$Original_name<-trimws(Norway_Taxa$Original_name)
Norway_Taxa<-Norway_Taxa[,c(1,47,48)]
Norway_Taxa<-Norway_Taxa[!is.na(Norway_Taxa$NISP),]
Norway_Taxa$database<-"Norway_VoR_13Aug2025"
Norway_Meta<-Norway_Meta %>%
  dplyr::mutate_all(as.character)
names(Norway_Meta)<-names(Belgium_Meta)
Norway_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Norway_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Ireland)[47:ncol(Ireland)])
Ireland[47:ncol(Ireland)] <- sapply(Ireland[47:ncol(Ireland)],as.numeric)
Ireland$`NSS unique ID` <- as.character(Ireland$`NSS unique ID`)
Ireland_Meta<-Ireland[,c(1:46)]
Ireland_Taxa<-melt(Ireland, id.vars = 1:46)
names(Ireland_Taxa)[47]<-"Original_name"
names(Ireland_Taxa)[48]<-"NISP"
Ireland_Taxa$Original_name<-trimws(Ireland_Taxa$Original_name)
Ireland_Taxa<-Ireland_Taxa[,c(1,47,48)]
Ireland_Taxa<-Ireland_Taxa[!is.na(Ireland_Taxa$NISP),]
Ireland_Taxa$database<-"Ireland_VoR_14Aug2025"
Ireland_Meta<-Ireland_Meta %>%
  dplyr::mutate_all(as.character)
names(Ireland_Meta)<-names(Belgium_Meta)
Ireland_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Ireland_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Poland)[47:ncol(Poland)])
Poland[47:ncol(Poland)] <- sapply(Poland[47:ncol(Poland)],as.numeric)
Poland$`NSS unique ID` <- as.character(Poland$`NSS unique ID`)
Poland_Meta<-Poland[,c(1:46)]
Poland_Taxa<-melt(Poland, id.vars = 1:46)
names(Poland_Taxa)[47]<-"Original_name"
names(Poland_Taxa)[48]<-"NISP"
Poland_Taxa$Original_name<-trimws(Poland_Taxa$Original_name)
Poland_Taxa<-Poland_Taxa[,c(1,47,48)]
Poland_Taxa<-Poland_Taxa[!is.na(Poland_Taxa$NISP),]
Poland_Taxa$database<-"Poland_VoR_14Aug2025"
Poland_Meta<-Poland_Meta %>%
  dplyr::mutate_all(as.character)
names(Poland_Meta)<-names(Belgium_Meta)
Poland_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Poland_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Scotland)[47:ncol(Scotland)])
Scotland[47:ncol(Scotland)] <- sapply(Scotland[47:ncol(Scotland)],as.numeric)
Scotland$`NSS unique ID` <- as.character(Scotland$`NSS unique ID`)
Scotland_Meta<-Scotland[,c(1:46)]
Scotland_Taxa<-melt(Scotland, id.vars = 1:46)
names(Scotland_Taxa)[47]<-"Original_name"
names(Scotland_Taxa)[48]<-"NISP"
Scotland_Taxa$Original_name<-trimws(Scotland_Taxa$Original_name)
Scotland_Taxa$Original_name<-gsub("  "," ",Scotland_Taxa$Original_name)
Scotland_Taxa<-Scotland_Taxa[,c(1,47,48)]
Scotland_Taxa<-Scotland_Taxa[!is.na(Scotland_Taxa$NISP),]
Scotland_Taxa$database<-"Scotland_VoR_14Aug2025"
Scotland_Meta<-Scotland_Meta %>%
  dplyr::mutate_all(as.character)
names(Scotland_Meta)<-names(Belgium_Meta)
Scotland_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Scotland_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Scotland_b)[47:ncol(Scotland_b)])
Scotland_b[47:ncol(Scotland_b)] <- sapply(Scotland_b[47:ncol(Scotland_b)],as.numeric)
Scotland_b$`NSS unique ID` <- as.character(Scotland_b$`NSS unique ID`)
Scotland_b_Meta<-Scotland_b[,c(1:46)]
Scotland_b_Taxa<-melt(Scotland_b, id.vars = 1:46)
names(Scotland_b_Taxa)[47]<-"Original_name"
names(Scotland_b_Taxa)[48]<-"NISP"
Scotland_b_Taxa$Original_name<-trimws(Scotland_b_Taxa$Original_name)
Scotland_b_Taxa$Original_name<-gsub("  "," ",Scotland_b_Taxa$Original_name)
Scotland_b_Taxa<-Scotland_b_Taxa[,c(1,47,48)]
Scotland_b_Taxa<-Scotland_b_Taxa[!is.na(Scotland_b_Taxa$NISP),]
Scotland_b_Taxa$database<-"Scotland_suppl_VoR_14Aug2025"
Scotland_b_Meta<-Scotland_b_Meta %>%
  dplyr::mutate_all(as.character)
names(Scotland_b_Meta)<-names(Belgium_Meta)
Scotland_b_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Scotland_b_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

str((Sweden)[47:ncol(Sweden)])
Sweden[47:ncol(Sweden)] <- sapply(Sweden[47:ncol(Sweden)],as.numeric)
Sweden$`NSS unique ID` <- as.character(Sweden$`NSS unique ID`)
Sweden_Meta<-Sweden[,c(1:46)]
Sweden_Taxa<-melt(Sweden, id.vars = 1:46)
names(Sweden_Taxa)[47]<-"Original_name"
names(Sweden_Taxa)[48]<-"NISP"
Sweden_Taxa$Original_name<-trimws(Sweden_Taxa$Original_name)
Sweden_Taxa<-Sweden_Taxa[,c(1,47,48)]
Sweden_Taxa<-Sweden_Taxa[!is.na(Sweden_Taxa$NISP),]
Sweden_Taxa$database<-"Sweden_VoR_13Aug2025"
Sweden_Meta<-Sweden_Meta %>%
  dplyr::mutate_all(as.character)
names(Sweden_Meta)<-names(Belgium_Meta)
Sweden_Meta %>%
  dplyr::summarize(Assemblage_count = n_distinct(`NSS unique ID`))
Sweden_Taxa %>%
  dplyr::summarise(total_nisp = sum(NISP))

df_Meta<-dplyr::bind_rows(list(Belgium_Meta, Denmark_Meta, 
                        England_Meta, 
                        Estonia_Meta, Germany_Hein_Meta,
                        Germany_Hol_Meta, Ireland_Meta, Netherlands_Rom_Meta, Netherlands_RomSup_Meta,
                        Norway_Meta, Poland_Meta, 
                        Scotland_Meta, Scotland_b_Meta, 
                        Sweden_Meta))
df_Meta$`Year reported`<-gsub("1995 \\(mentioned in cited article\\)",1995,df_Meta$`Year reported`)
df_Meta$`Number of unidentified fish specimens`<-gsub("\\?","",df_Meta$`Number of unidentified fish specimens`)
df_Meta$`NSS unique ID`
df_Meta <-df_Meta %>% dplyr::mutate_at(c('NSS unique ID','Start date CE', 'End date CE',
                                  'Decimal Latitude','Decimal Longitude', 'Year reported',
                                  'Number of unidentified fish specimens',
                                  'ID for linked rows differing only by recovery', 'Minimum sieve mesh mm',
                                  'Maximum sieve mesh mm', 'Total Fish NISP'), as.numeric)
df_Meta<-df_Meta[!is.na(df_Meta$`NSS unique ID`),]
#2. Merge taxa with taxa metadata and assign quality scores before assigning taxonomy
df_NISP<-dplyr::bind_rows(list(Belgium_Taxa, Denmark_Taxa, 
                        England_Taxa, 
                        Estonia_Taxa, Germany_Hein_Taxa,
                        Germany_Hol_Taxa, Ireland_Taxa, Netherlands_Rom_Taxa, Netherlands_RomSup_Taxa,
                        Norway_Taxa, Poland_Taxa, 
                        Scotland_Taxa, Scotland_b_Taxa, 
                        Sweden_Taxa))
sum(as.numeric(c(Belgium_Taxa$NISP,
                 Denmark_Taxa$NISP,
                 England_Taxa$NISP, 
                 Estonia_Taxa$NISP, 
                 Germany_Hein_Taxa$NISP,
                 Germany_Hol_Taxa$NISP, 
                 Ireland_Taxa$NISP, 
                 Netherlands_Rom_Taxa$NISP,
                 Netherlands_RomSup_Taxa$NISP,
                 Norway_Taxa$NISP, Poland_Taxa$NISP, 
                 Scotland_Taxa$NISP, Scotland_b_Taxa$NISP, 
                 Sweden_Taxa$NISP)), na.rm=T)
str(df_NISP)
df_NISP$`NSS unique ID`<-as.character(df_NISP$`NSS unique ID`)
x<-df_NISP %>% dplyr::group_by(`NSS unique ID`, NISP) %>% dplyr::summarise(NISP_all = sum(NISP))
x<-aggregate(NISP ~ `NSS unique ID`, data = df_NISP, FUN= sum)
y<-df_Meta[,c(1,46)]
y$`NSS unique ID`<-as.character(y$`NSS unique ID`)
tmp<-y[!y$`NSS unique ID` %in% x$`NSS unique ID`,]
tmp<-tmp$`NSS unique ID`
df_NISP$Has_taxa<-"TRUE"
df_NISP$`NSS unique ID`<-as.character(df_NISP$`NSS unique ID`)
df_NISP<-df_NISP %>% add_row(`NSS unique ID` = tmp, database = "Belgium_VoR_14Aug2025", Has_taxa = "FALSE")
df_NISP$`NSS unique ID`<-str_replace(df_NISP$`NSS unique ID`,"\\.0","")

###Create Unique ID Variable for Recommender systems
df_Meta$`Assemblage or sub-assemblage`<-gsub(" ", ",", df_Meta$`Assemblage or sub-assemblage`)
df_Meta$Assemblage_ID<-paste(df_Meta$Database, df_Meta$`Site name`, df_Meta$`Settlement modern name`, df_Meta$`Start date CE`, df_Meta$`End date CE`, df_Meta$`Decimal Latitude`, df_Meta$`Decimal Longitude`)
df_Meta$Assemblage_ID<-gsub(" ",";",df_Meta$Assemblage_ID)
df_Meta<-df_Meta %>%
  dplyr::mutate_all(as.character)
df_NISP<-df_NISP %>%
  dplyr::mutate_all(as.character)

sum(as.numeric(df_NISP$NISP), na.rm=T)

df_NISP %>%
  dplyr::group_by(database) %>%
  dplyr::summarise(
    total_NISP = sum(as.numeric(NISP), na.rm = TRUE)
  )

df_NISP %>%
  dplyr::group_by(database) %>%
  dplyr::summarise(
    total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
    Assemblage_count = n_distinct(`NSS unique ID`),
    Taxa = n_distinct(`Original_name`)
  )

df_NISP %>%
  dplyr::summarise(
    total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
    Assemblage_count = n_distinct(`NSS unique ID`),
    Taxa = n_distinct(`Original_name`)
  )

df_NISP %>%
  dplyr::summarise(
    total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
    Assemblage_count = n_distinct(`NSS unique ID`),
    Taxa = n_distinct(`Original_name`)
  )

df_NISP_new %>%
  dplyr::group_by(database) %>%
  dplyr::summarise(
    species_count = n_distinct(`GBIF_species`),
    genera_count = n_distinct(`GBIF_genus`),
    family_count = n_distinct(`GBIF_family`),
    order_count = n_distinct(`GBIF_order`),
    order_count = n_distinct(`GBIF_class`)
  )

df_NISP_new %>%
  dplyr::group_by(database) %>%
  dplyr::summarise(
    Assemblage_count = n_distinct(`NSS unique ID`),
    Taxa = n_distinct(`Taxa_cleaned`)
  )

saveRDS(c(df_NISP, df_Meta), file = "NSS_longformat_Aug2025.R")
write_xlsx(df_Meta, "NSS_metadata_Aug2025.xlsx", col_names = TRUE, format_headers = TRUE)
write_xlsx(df_NISP, "NSS_NISPData_Aug2025.xlsx", col_names = TRUE, format_headers = TRUE)
df_Meta <-df_Meta %>% dplyr::mutate_at(c('NSS unique ID'), as.character)
df_all<-df_Meta %>% left_join(df_NISP)
df_all<-df_all %>%
  dplyr::mutate_all(as.character)
sum(as.numeric(df_all$NISP), na.rm=T)
sum(as.numeric(df_NISP$NISP), na.rm=T)
sum(as.numeric(df_Meta$`Total Fish NISP`), na.rm=T)

###Join all datasets together with taxonomic information
taxa_new<-taxa_meta
df_NISP_new<-merge(df_all, taxa_new, all.x = T, by="Original_name")
sum(as.numeric(df_NISP_new$NISP), na.rm=T)
df_NISP_new<-df_NISP_new[!is.na(df_NISP_new$`NSS unique ID`),]
new<-df_NISP_new[is.na(df_NISP_new$Taxa_cleaned),]
df_NISP_new<-df_NISP_new[!is.na(df_NISP_new$Taxa_cleaned),]

###Summary stats of taxon (old)
#cols <- c("Taxa_cleaned", "GBIF_class", "GBIF_order", 
#          "GBIF_family", "GBIF_genus", "GBIF_species")
#
## Combine results
#results <- map_dfr(cols, ~get_counts(df_NISP_new, .x)) %>%
#  bind_rows(map_dfr(cols, ~get_overall(df_NISP_new, .x)))
#
#check<-data.frame(results)

#Summary stats of taxon and counts ####
totalcounts<-df_NISP_new %>% 
  dplyr::summarise(SiteCount = n_distinct(`Site name`),
            AssemblageCount = n_distinct(`Assemblage or sub-assemblage`),
            total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
            SpeciesCount = n_distinct(`GBIF_species`),
            GeneraCount = n_distinct(`GBIF_genus`),
            FamilyCount = n_distinct(`GBIF_family`),
            OrderCount = n_distinct(`GBIF_order`),
            ClassCount = n_distinct(`GBIF_class`))

databasecounts<-df_NISP_new %>% 
  dplyr::group_by(database) %>%
  dplyr::summarise(SiteCount = n_distinct(`Site name`),
            AssemblageCount = n_distinct(`Assemblage or sub-assemblage`),
            total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
            SpeciesCount = n_distinct(`GBIF_species`),
            GeneraCount = n_distinct(`GBIF_genus`),
            FamilyCount = n_distinct(`GBIF_family`),
            OrderCount = n_distinct(`GBIF_order`),
            ClassCount = n_distinct(`GBIF_class`))

write_xlsx(totalcounts, "totalcounts.xlsx")
write_xlsx(databasecounts, "databasecounts.xlsx")

###Sort out typos
df_NISP_new$`Assemblage date as locally defined`[df_NISP_new$`ID for linked rows differing only by recovery`=="2214"]<-"late 7th/early 8th to early-mid 8th"
df_NISP_new$`Assemblage date as locally defined`[df_NISP_new$`ID for linked rows differing only by recovery`=="2350"]<-"12th; EMED"
df_NISP_new$`Assemblage date as locally defined`[df_NISP_new$`Site name`=="Bostadh Beach"]<-"Late Iron Age/Norse"
df_NISP_new$`Year reported`[df_NISP_new$`NSS unique ID`=="60330"]<-"2016"
df_NISP_new$`Assemblage date as locally defined`[df_NISP_new$`NSS unique ID`==10009]<-"MED c.1100 to c.1350"

###Remove assemblages without start and end dates
df_all<-df_NISP_new[!is.na(df_NISP_new$`Start date CE`),]
df_all<-df_all[!is.na(df_all$`End date CE`),]
df_all<-df_all[!df_all$`Data quality (green amber or red)`=="red",]

###Sort out Assemblage metadata
df_old<-df_all
df_old$DB_ID2<-trimws(paste(df_old$`Site name`,";",df_old$`Settlement modern name`,";",df_old$`Assemblage date as locally defined`,";",df_old$`Start date CE`,";",df_old$`End date CE`,";",round(as.numeric(df_old$`Decimal Longitude`),3),";",round(as.numeric(df_old$`Decimal Latitude`),3)))
df_old$DB_ID2<-gsub(" ","",df_old$DB_ID2)
df_old_short<-df_old %>% dplyr::select(DB_ID2, `Recovery method (hand-collected; sieved; both; unknown)`) %>% distinct()
names(df_old_short)[2]<-"Recovery_method"
df_old_short<-dcast(df_old_short, DB_ID2 ~ Recovery_method, value.var="Recovery_method", fun.aggregate = length)
df_old_short$Recovery<-"NA"
names(df_old_short)[4]<-"hand"
df_old_short$Recovery <- df_old_short %>%
  dplyr::mutate(Recovery = case_when(
    both == 0 & complete == 0 & hand == 0 & sieved == 0 & unknown == 1 ~ "Unknown",
    both == 1 & complete %in% c(0,1) & hand %in% c(0,1) & sieved %in% c(0,1) & unknown %in% c(0,1) ~ "both",
    sieved == 1 & both == 0 & complete == 0 & hand == 0 & unknown == 0 ~ "Sieved only",
    hand == 1 & both == 0 & complete == 0 & sieved == 0 & unknown == 0 ~ "Hand-collected only",
    hand == 1 & both == 0 & complete == 0 & sieved == 1 & unknown == 0 ~ "both",
    hand == 0 & both == 0 & complete == 1 & sieved == 0 & unknown == 0 ~ "both",
    hand == 0 & both == 0 & complete == 0 & sieved == 1 & unknown == 1 ~ "both",
    hand == 1 & both == 0 & complete == 0 & sieved == 0 & unknown == 1 ~ "Unknown",
    hand == 1 & both == 0 & complete == 0 & sieved == 1 & unknown == 1 ~ "both",
    TRUE ~ "CHECK"
  ))
all_df<- df_old %>% merge(df_old_short, all.x=T)

#Add DB_Sieved column (True/False)
df_old_short<-df_old_short$Recovery[,c(1,8)] %>% distinct()
df_old_short<- df_old_short %>%
  dplyr::mutate(DB_sieved = case_when(
    Recovery %in% c("both","Sieved only") ~ TRUE,
    Recovery == "Hand-collected only" ~ FALSE,
    TRUE ~ NA
  ))
all_df<- df_old %>% merge(df_old_short, all.x=T)

#Minimum sieve size
all_df$DB_ID3<-trimws(paste(all_df$`Site name`,";",all_df$`Settlement modern name`,";",all_df$`Assemblage date as locally defined`,";",
                            all_df$Recovery,";",all_df$DB_sieved,";",all_df$`Start date CE`,";",all_df$`End date CE`,";",round(as.numeric(all_df$`Decimal Longitude`),3),";",round(as.numeric(all_df$`Decimal Latitude`),3)))
all_df$DB_ID3<-gsub(" ","",all_df$DB_ID3)
df_old_short<-all_df %>% dplyr::select(DB_ID3, `Minimum sieve mesh mm`, Recovery) %>% distinct()
names(df_old_short)[2]<-"Minimum_sieve_mesh"
df_old_short$Minimum_sieve_mesh<-as.numeric(df_old_short$Minimum_sieve_mesh)
df_old_short<-aggregate(Minimum_sieve_mesh ~ DB_ID3 + Recovery, df_old_short, min)
df_old_short$Minimum_sieve_mesh[df_old_short$Recovery == "Hand-collected only"]<-NA
all_df<- all_df %>% merge(df_old_short, all.x=T)

#List of NSS IDs
df_old_short<-all_df %>% dplyr::select(DB_ID3, `NSS unique ID`) %>% distinct()
names(df_old_short)[2]<-"Lumped_NSS_IDs"
df_old_short$Lumped_NSS_IDs<-as.character(df_old_short$Lumped_NSS_IDs)
df_old_short<-aggregate(Lumped_NSS_IDs ~ DB_ID3, data = df_old_short, paste, collapse = ",")
all_df<- all_df %>% merge(df_old_short, all.x=T)

#List of Green/Amber flags
df_old_short<-all_df %>% dplyr::select(DB_ID3, `Data quality (green amber or red)`) %>% distinct()
names(df_old_short)[2]<-"Data_quality"
df_old_short$Data_quality<-as.character(df_old_short$Data_quality)
df_old_short<-aggregate(Data_quality ~ DB_ID3, data = df_old_short, paste, collapse = ",")
all_df<- all_df %>% merge(df_old_short, all.x=T)

###Export data (Non-Aggregated) 
names(all_df)
df_final<-all_df[,c(1,3,57,5:54,2,4,56,59,60:87)] %>% distinct()
df_final<-df_final[!is.na(df_final$`NSS unique ID`),]
df_final$NISP<-as.numeric(df_final$NISP)
x<-df_final %>% 
  dplyr::group_by(Database) %>%
  dplyr::summarise(Sum_X = sum(NISP, na.rm = TRUE))
df_final %>% 
  dplyr::summarise(Sum_X = sum(NISP, na.rm = TRUE))

df_final$`Decimal Latitude`<-round(as.numeric(df_final$`Decimal Latitude`),3)
df_final$`Decimal Longitude`<-round(as.numeric(df_final$`Decimal Longitude`),3)

df_final %>% 
  dplyr::summarise(SiteCount = n_distinct(`Site name`),
                   AssemblageCount = n_distinct(`Assemblage or sub-assemblage`),
                   total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
                   SpeciesCount = n_distinct(`GBIF_species`),
                   GeneraCount = n_distinct(`GBIF_genus`),
                   FamilyCount = n_distinct(`GBIF_family`),
                   OrderCount = n_distinct(`GBIF_order`),
                   ClassCount = n_distinct(`GBIF_class`))

df_final %>% 
  dplyr::group_by(database) %>%
  dplyr::summarise(SiteCount = n_distinct(`Site name`),
                   AssemblageCount = n_distinct(`Assemblage or sub-assemblage`),
                   total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
                   SpeciesCount = n_distinct(`GBIF_species`),
                   GeneraCount = n_distinct(`GBIF_genus`),
                   FamilyCount = n_distinct(`GBIF_family`),
                   OrderCount = n_distinct(`GBIF_order`),
                   ClassCount = n_distinct(`GBIF_class`))

df_final$region<-dplyr::recode(df_final$Country,
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


df_final<-df_final %>% 
  dplyr::filter(!is.na(region))

df_final %>% 
  dplyr::filter(!is.na(GBIF_family),
         !GBIF_family=="NA") %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(SiteCount = n_distinct(`Site name`),
                   AssemblageCount = n_distinct(`Assemblage or sub-assemblage`),
                   total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
                   SpeciesCount = n_distinct(`GBIF_species`),
                   GeneraCount = n_distinct(`GBIF_genus`),
                   FamilyCount = n_distinct(`GBIF_family`),
                   OrderCount = n_distinct(`GBIF_order`),
                   ClassCount = n_distinct(`GBIF_class`))

df_final %>% 
  dplyr::summarise(SiteCount = n_distinct(`Site name`),
                   AssemblageCount = n_distinct(`Assemblage or sub-assemblage`),
                   total_NISP = sum(as.numeric(NISP), na.rm = TRUE),
                   SpeciesCount = n_distinct(`GBIF_species`),
                   GeneraCount = n_distinct(`GBIF_genus`),
                   FamilyCount = n_distinct(`GBIF_family`),
                   OrderCount = n_distinct(`GBIF_order`),
                   ClassCount = n_distinct(`GBIF_class`))

write_xlsx(df_NISP_new, "All_Data_with_Chronology_Sept2025.xlsx", col_names = TRUE, format_headers = TRUE)
saveRDS(df_NISP_new, file = "All_Data_with_Chronology_Sept2025.R")

