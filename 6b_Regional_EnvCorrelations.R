#Plot Figure 2 for NSS - Env Proxies 
#Author: D L Buss
#Load libraries/data ####
library(dplyr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(corrplot)
library(psych)
library(tidyverse)
library(reshape2)
DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69", "Gold3")
DB2<-c("#14517B","#FE994F","#6A8A73","#8e7d69","#F8F1E9", "Gold3")
DB3<-c("#6A8A73","#FE994F","#14517B","grey5","darkred","#006ba4","grey18","grey40","#D6604D","#2B547E",
       "gold2","#F39B7FFF",
       "#91D1C2FF","#8e7d69","#0072B2",
       "#6A8A73","#FE994F","#8491B4FF","#14517B","#FE994F","#6A8A73","#8e7d69","#F8F1E9", "Gold3")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("MANUSCRIPT_functions.R")

###NSS-DB Regional correlations with Env Var. ####
final<-read_xlsx("../Data/MTC_output.xlsx", col_names=T)
final2<-final %>%
  filter(sizebins == "100") %>% select(TimePeriod, MTC, plot_names)
final2$TimePeriod<-as.character(final2$TimePeriod)
final2<-dcast(final2, TimePeriod ~ plot_names, value.var = "MTC")
#final2$Year<-as.numeric(final2$Year)

#Load and reformat proxies
df_fell<-read_xlsx("../Data/Env_proxies/Felling_Dates/Felling dates_Ljungqvist et al 2022_2.xlsx", col_names=T)
df_fell_b<-read_xlsx("../Data/Env_proxies/Felling_Dates/Fellingcounts_Ljungqvist et al 2022_2_avg.xlsx", col_names=T)
df_fell$all<-1
overtime<-aggregate(all ~ Year, df_fell, FUN=sum)
overtime$all[overtime$all < 4]<-NA
overtime$title<-"Proxy of demographic change from Ljungqvist et al. (2022)"
names(overtime)[2]<-"Felling_Counts"
#overtime$TimePeriod<-as.character(overtime$Year)
overtime<-overtime[!is.na(overtime$Felling_Counts),]
overtime<-overtime %>%
  mutate(TimePeriod = (Year %/% 100) * 100) %>%
  group_by(TimePeriod) %>%
  summarize(Total_Count = sum(Felling_Counts), .groups = "drop") %>%
  arrange(TimePeriod)
names(overtime)[2]<-"Felling_Counts_Ljungqvist_2022"
final2<-merge(final2, overtime[,c(1:2)], all.x=T, by="TimePeriod")

df_SST<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Berner, K.S., N. Koç, F. Godtliebsen, and D. Divine. 2011..xlsx", sheet=2)
head(df_SST)
str(df_SST)
names(df_SST)[1]<-"TimePeriod"
names(df_SST)[3]<-"SST_25yr_av"
names(df_SST)[4]<-"SST_Berner_2011"
final2<-merge(final2, df_SST[,c(1,4)], all.x=T, by="TimePeriod")

df_CS<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Fohlmeister, J., Schröder-Ritzrau, A., Scholz...2012.xlsx",sheet=2)
head(df_CS)
str(df_CS)
names(df_CS)[1]<-"TimePeriod"
names(df_CS)[3]<-"CS_25yr_av"
names(df_CS)[4]<-"Precip_Fohlmeister_2012)"
final2<-merge(final2, df_CS[,c(1,4)], all.x=T, by="TimePeriod")
#oxygen_lab<-expression(paste(delta^18, "O (\u2030)"))

df_Forams<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Keigwin_1996.xlsx", sheet=2)
head(df_Forams)
str(df_Forams)
df_Forams$`25yr_avg`<-as.numeric(df_Forams$`25yr_avg`)
df_Forams$`100yr_avg`<-as.numeric(df_Forams$`100yr_avg`)
names(df_Forams)[1]<-"TimePeriod"
names(df_Forams)[2]<-"OI_25yr_av"
names(df_Forams)[4]<-"d18O_Keigwin_1996"
final2<-merge(final2, df_Forams[,c(1,4)], all.x=T, by="TimePeriod")

df_SST_Af<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_McGregor, H.V., M. Dima, H.W. Fischer, and S. Mulitza 2007.xlsx", sheet=2)
df_SST_Af$`25yr_avg`<-as.numeric(df_SST_Af$`25yr_avg`)
df_SST_Af$`100yr_avg`<-as.numeric(df_SST_Af$`100yr_avg`)
names(df_SST_Af)[1]<-"TimePeriod"
names(df_SST_Af)[4]<-"SST_25yr_av"
names(df_SST_Af)[5]<-"SST_McGregor_2011"
final2<-merge(final2, df_SST_Af[,c(1,5)], all.x=T, by="TimePeriod")

df_SST_SP<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Miettinen, A., D. Divine, N. Koç¸ F. Godtliebsen, and I.R. Hall. 2012..xlsx", sheet=2)
head(df_SST_SP)
str(df_SST_SP)
df_SST_SP$`25yr_avg`<-as.numeric(df_SST_SP$`25yr_avg`)
df_SST_SP$`100yr_avg`<-as.numeric(df_SST_SP$`100yr_avg`)
names(df_SST_SP)[1]<-"TimePeriod"
names(df_SST_SP)[3]<-"SST_25yr_av"
names(df_SST_SP)[4]<-"SST_Miettinen_2012"
final2<-merge(final2, df_SST_SP[,c(1,4)], all.x=T, by="TimePeriod")

df_SST_SP2<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Sicre, M.-A., I.R. Hall, J. Mignot, M. Khodri, U. Ezat.xlsx", sheet=2)
head(df_SST_SP2)
str(df_SST_SP2)
df_SST_SP2$`25yr_Rolling_Avg`<-as.numeric(df_SST_SP2$`25yr_Rolling_Avg`)
df_SST_SP2$`100yr_Rolling_Avg`<-as.numeric(df_SST_SP2$`100yr_Rolling_Avg`)
names(df_SST_SP2)[1]<-"TimePeriod"
names(df_SST_SP2)[3]<-"SST_25yr_av"
names(df_SST_SP2)[4]<-"SST_Sicre_2011"
final2<-merge(final2, df_SST_SP2[,c(1,4)], all.x=T, by="TimePeriod")

df_XRAY<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Tiljander, M., M. Saarnisto, A.E.K. Ojala, and T. Saarinen. 2003..xlsx", sheet=2)
head(df_XRAY)
str(df_XRAY)
df_XRAY$`25yr_Rolling_Avg`<-as.numeric(df_XRAY$`25yr_Rolling_Avg`)
df_XRAY$`100yr_Rolling_Avg`<-as.numeric(df_XRAY$`100yr_Rolling_Avg`)
names(df_XRAY)[1]<-"TimePeriod"
names(df_XRAY)[3]<-"XR_25yr_av"
names(df_XRAY)[4]<-"XR_Tiljander_2003"
final2<-merge(final2, df_XRAY[,c(1,4)], all.x=T, by="TimePeriod")

df_pages<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_Pages2K.xlsx", sheet=1)
names(df_pages)[1]<-"TimePeriod"
names(df_pages)[3]<-"Europe2K_25yr_av"
names(df_pages)[4]<-"Europe2K_100yr_av"
names(df_pages)[12]<-"Arctic2K_25yr_av"
names(df_pages)[13]<-"Arctic2K_100yr_av"
final2<-merge(final2, df_pages[,c(1,4)], all.x=T, by="TimePeriod")
final2<-merge(final2, df_pages[,c(1,13)], all.x=T, by="TimePeriod")

###Cross and Pearson's correlations with MTC (global) ####
names(final2)<-gsub(" 100-yr","", names(final2))
names(final2)<-gsub("\\(no CEE\\)","", names(final2))
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(3,7:15)) %>% as.matrix(.)
cor.matrix<-cor(cors)

corrplot::corrplot(cor.matrix, method = "number", type = "lower",
         tl.col = "black", tl.cex = 0.75,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))

png("../SupplementaryFiles/SuppFigureCorr_GlobalMTC.png", width = 6, height = 6, units = "in", res = 600)
corrplot::corrplot(cor.matrix, method = "number", type = "lower",
                   tl.col = "black", tl.cex = 0.75,
                   col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###Cross and Pearson's correlations with MTC (regional)
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(2,4,5,7:15)) %>% as.matrix(.)
cor.matrix<-cor(cors)
png("../SupplementaryFiles/SuppFigureCorr_RegionalMTC.png", width = 6, height = 6, units = "in", res = 600)
corrplot::corrplot(cor.matrix, method = "number", type = "lower",
                   tl.col = "black", tl.cex = 0.75,
                   col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()
###Cross and Pearson's correlations between MTCs
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(2:5)) %>% as.matrix(.)
cor.matrix<-cor(cors)
png("../SupplementaryFiles/SuppFigureCorr_MTConly.png", width = 6, height = 6, units = "in", res = 600)
corrplot::corrplot(cor.matrix, method = "number", type = "lower",
                   tl.col = "black", tl.cex = 0.75,
                   col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###Significant values only ####
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(3,7:15)) %>% as.matrix(.)
cor.matrix<-cor(cors)
p_values <- cor.mtest(cor.matrix)
png("../SupplementaryFiles/SuppFigureCorr_SigGlobalMTC.png", width = 6, height = 6, units = "in", res = 600)
corrplot(cor.matrix, method = "number",
         main = "",
         p.mat = p_values,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###Cross and Pearson's correlations with MTC (regional) ####
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(2:5,7:15)) %>% as.matrix(.)
cor.matrix<-cor(cors)
p_values <- cor.mtest(cor.matrix)
png("../SupplementaryFiles/SuppFigureCorr_SigRegionalMTC.png", width = 6, height = 6, units = "in", res = 600)
corrplot(cor.matrix, method = "number",
         main = "",
         p.mat = p_values,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###With Bonferonni corrections ####
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(3,7:15)) %>% as.matrix(.)
corr_results <- corr.test(cors, adjust = "bonferroni") 
cor_matrix <- corr_results$r
p_matrix <- corr_results$p  
png("../SupplementaryFiles/SuppFigureCorr_SigGlobalMTC.png", width = 7.5, height = 7.5, units = "in", res = 600)
corrplot(cor_matrix, method = "number",
         main = "",
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(2:5,7:15)) %>% as.matrix(.)
corr_results <- corr.test(cors, adjust = "bonferroni") 
cor_matrix <- corr_results$r
p_matrix <- corr_results$p  
png("../SupplementaryFiles/SuppFigureCorr_SigRegionalMTC.png", width = 7.5, height = 7.5, units = "in", res = 600)
corrplot(cor_matrix, method = "number",
         main = "",
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###Pearsons correlations split before and after the little ice age ####
final2$segment<-ifelse(as.numeric(final2$TimePeriod) < 1400, "First Millennium", "Second Millennium")
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100",
                        !segment == "Second Millennium") %>%
  select(c(3,7:15)) %>% as.matrix(.)
corr_results <- corr.test(cors, adjust = "bonferroni") 
cor_matrix <- corr_results$r
p_matrix <- corr_results$p  
png("../SupplementaryFiles/SuppFigureCorr_SigPreLittleIceAgeMTC.png", width = 7.5, height = 7.5, units = "in", res = 600)
corrplot(cor_matrix, method = "number",
         main = "",
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###Pearsons correlations split before and after the little ice age ####
final2$segment<-ifelse(as.numeric(final2$TimePeriod) >= 1400, "First Millennium", "Second Millennium")
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100",
                        !segment == "Second Millennium") %>%
  select(c(3,7:15)) %>% as.matrix(.)
corr_results <- corr.test(cors, adjust = "bonferroni") 
cor_matrix <- corr_results$r
p_matrix <- corr_results$p  
png("../SupplementaryFiles/SuppFigureCorr_SigPostLittleIceAgeMTC.png", width = 7.5, height = 7.5, units = "in", res = 600)
corrplot(cor_matrix, method = "number",
         main = "",
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###Regional Pearsons correlations split before and after the little ice age ####
final2$segment<-ifelse(as.numeric(final2$TimePeriod) < 1400, "First Millennium", "Second Millennium")
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100",
                        !segment == "Second Millennium") %>%
  select(c(2:5,7:15)) %>% as.matrix(.)
corr_results <- corr.test(cors, adjust = "bonferroni") 
cor_matrix <- corr_results$r
p_matrix <- corr_results$p  
png("../SupplementaryFiles/SuppFigureCorr_SigPreLittleIceAge_RegionalMTC.png", width = 7.5, height = 7.5, units = "in", res = 600)
corrplot(cor_matrix, method = "number",
         main = "",
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###Regional Pearsons correlations split before and after the little ice age ####
final2$segment<-ifelse(as.numeric(final2$TimePeriod) >= 1400, "First Millennium", "Second Millennium")
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100",
                        !segment == "Second Millennium") %>%
  select(c(2:5,7:15)) %>% as.matrix(.)
corr_results <- corr.test(cors, adjust = "bonferroni") 
cor_matrix <- corr_results$r
p_matrix <- corr_results$p  
png("../SupplementaryFiles/SuppFigureCorr_SigPostLittleIceAge_RegionalMTC.png", width = 7.5, height = 7.5, units = "in", res = 600)
corrplot(cor_matrix, method = "number",
         main = "",
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()

###Regional Pearsons correlations split before and after the little ice age ####
final2$segment<-ifelse(as.numeric(final2$TimePeriod) >= 1400, "First Millennium", "Second Millennium")
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(4,7:15)) %>% as.matrix(.)
corr_results <- corr.test(cors, adjust = "bonferroni") 
cor_matrix <- corr_results$r
p_matrix <- corr_results$p  
corrplot(cor_matrix, method = "number",
         main = "",
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))

###Cross correlations ####
# Get all variable pairs
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100") %>%
  select(c(2:5,7:15))
var_pairs <- combn(names(cors), 2, simplify = FALSE)
get_ccf(cors$`Britain & Ireland`,cors$`NSS-DB `)
ccf(cors$`Britain & Ireland`,cors$`NSS-DB `)

# Calculate CCF for all pairs
ccf_results <- purrr::map(var_pairs, ~ {
  ccf_data <- get_ccf(cors[[.x[1]]], cors[[.x[2]]])
  ccf_data$Var1 <- .x[1]
  ccf_data$Var2 <- .x[2]
  return(ccf_data)
})

# Combine into one data frame
ccf_df <- bind_rows(ccf_results)
conf_level <- 1.96 / sqrt(12)

ccf_df_lower <- ccf_df %>%
  filter(Var1 > Var2)

png("../SupplementaryFiles/CrossCorrelations_all.png", width = 7.5, height = 7.5, units = "in", res = 600)
ggplot(ccf_df_lower, aes(x = lag, y = ccf)) +
  geom_line() +
  geom_hline(yintercept = c(-conf_level, conf_level), linetype = "dashed", color = "blue") +
  facet_grid(Var1 ~ Var2) +  # Retains original shape, but only lower diagonal
  theme_minimal() +
  xlim(-4, 4) +
  ylim(-0.7, 0.7) +
  labs(title = "Cross-Correlations", x = "Lag", y = "CCF")
dev.off()
# Apply the function to all variable pairs and combine the results
ccf_results <- purrr::map_dfr(var_pairs, ~ {
  get_significant_ccf(cors[[.x[1]]], cors[[.x[2]]], .x, conf_level = conf_level)
})

ccf_results_pairs<-ccf_results[,c(1:2)] %>% distinct()
var_pairs <- combn(ccf_results_pairs, 2, simplify = FALSE)

plots_per_page <- 6
total_plots <- length(var_pairs[[1]]$Var1)
num_pages <- ceiling(total_plots / plots_per_page)

for (page in 1:num_pages) {
  # Define the filename for the PNG file, incorporating the page number
  png_filename <- paste0("CCF_plots_page_", page, ".png")
  
  # Open a PNG device for the current plot
  png(png_filename, width = 550, height = 650)
  
  # Set up the plotting layout for a 3x2 grid
  par(mfrow = c(3, 2), mar = c(4, 4, 4, 1))  # Adjust margins as needed
  
  # Determine the indices for the plots on this page
  start_index <- (page - 1) * plots_per_page + 1
  end_index <- min(page * plots_per_page, total_plots)
  
  # Generate and display CCF plots for the current page
  for (i in start_index:end_index) {
    # Extract variable names
    var1 <- var_pairs[[1]]$Var1[i]
    var2 <- var_pairs[[1]]$Var2[i]
    
    # Generate the CCF plot
    ccf_plot <- ccf(cors[[var1]], cors[[var2]], plot = TRUE, lag.max = 5, main = paste(var1, "vs", var2))
  }
  dev.off()
}

####Pre 1300s ####
###Cross correlations ####
# Get all variable pairs
final2$time<-as.numeric(final2$TimePeriod)
cors<-final2 %>% filter(!TimePeriod == "1900",
                        !TimePeriod == "100",
                        time > 1301) %>%
  select(c(2:5,7:15))
var_pairs <- combn(names(cors), 2, simplify = FALSE)

# Calculate CCF for all pairs
ccf_results <- purrr::map(var_pairs, ~ {
  ccf_data <- get_ccf(cors[[.x[1]]], cors[[.x[2]]])
  ccf_data$Var1 <- .x[1]
  ccf_data$Var2 <- .x[2]
  return(ccf_data)
})

# Combine into one data frame
ccf_df <- bind_rows(ccf_results)
conf_level <- 1.96 / sqrt(12)

ccf_df_lower <- ccf_df %>%
  filter(Var1 > Var2)

plots_per_page <- 9
total_plots <- length(var_pairs)
num_pages <- ceiling(total_plots / plots_per_page)

for (page in 1:num_pages) {
  # Define the filename for the PNG file, incorporating the page number
  png_filename <- paste0("CCF_plots_page_", page, ".png")
  
  # Open a PNG device for the current plot
  png(png_filename, width = 550, height = 650)
  
  # Set up the plotting layout for a 3x2 grid
  par(mfrow = c(3, 3), mar = c(4, 4, 4, 1))  # Adjust margins as needed
  
  # Determine the indices for the plots on this page
  start_index <- (page - 1) * plots_per_page + 1
  end_index <- min(page * plots_per_page, total_plots)
  
  # Generate and display CCF plots for the current page
  for (i in start_index:end_index) {
    # Extract variable names
    var1 <- var_pairs[[i]][1]
    var2 <- var_pairs[[i]][2]
    
    # Generate the CCF plot
    ccf_plot <- ccf(cors[[var1]], cors[[var2]], plot = TRUE, lag.max = 5, main = paste(var1, "vs", var2))
  }
  dev.off()
}

png("../SupplementaryFiles/CrossCorrelations_Pre1300.png", width = 7.5, height = 7.5, units = "in", res = 600)
ggplot(ccf_df_lower, aes(x = lag, y = ccf)) +
  geom_line() +
  geom_hline(yintercept = c(-conf_level, conf_level), linetype = "dashed", color = "blue") +
  facet_grid(Var1 ~ Var2) +  # Retains original shape, but only lower diagonal
  theme_minimal() +
  xlim(-4, 4) +
  ylim(-0.7, 0.7) +
  labs(title = "Cross-Correlations", x = "Lag", y = "CCF")
dev.off()
# Apply the function to all variable pairs and combine the results
ccf_results <- purrr::map_dfr(var_pairs, ~ {
  get_significant_ccf(cors[[.x[1]]], cors[[.x[2]]], .x, conf_level = conf_level)
})
ccf_df <- bind_rows(ccf_results)

plots_per_page <- 6
total_plots <- length(var_pairs[[1]]$Var1)
num_pages <- ceiling(total_plots / plots_per_page)

for (page in 1:num_pages) {
  # Define the filename for the PNG file, incorporating the page number
  png_filename <- paste0("CCF_plots_page_", page, ".png")
  
  # Open a PNG device for the current plot
  png(png_filename, width = 550, height = 650)
  
  # Set up the plotting layout for a 3x2 grid
  par(mfrow = c(3, 2), mar = c(4, 4, 4, 1))  # Adjust margins as needed
  
  # Determine the indices for the plots on this page
  start_index <- (page - 1) * plots_per_page + 1
  end_index <- min(page * plots_per_page, total_plots)
  
  # Generate and display CCF plots for the current page
  for (i in start_index:end_index) {
    # Extract variable names
    var1 <- var_pairs[[1]]$Var1[i]
    var2 <- var_pairs[[1]]$Var2[i]
    
    # Generate the CCF plot
    ccf_plot <- ccf(cors[[var1]], cors[[var2]], plot = TRUE, lag.max = 5, main = paste(var1, "vs", var2))
  }
  dev.off()
}
