##24th Feb 2025 - NSS analysis - ecological diversity changes
#Author: Dr Danny L Buss
remove.packages("TMB")
install.packages("glmmTMB", type="source")

#1. load libraries ####
library(readxl)
library(ggplot2)
library(dplyr)
library(reshape2)
library(writexl)
library(glmmTMB)
library(visreg)
library(broom)
library(MASS)
library(gridExtra)
library(gtools)

source("~/Documents/4-Oceans/3.NSS/2.MANUSCRIPT_Jan2025/Code_final/MANUSCRIPT_functions.R")

#2. setwd & load data ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5a_NSS_SpeciesData_Feb2025_withCEE.xlsx")
df_species<-read_excel("../SupplementaryFiles/SupplementaryTable5b_NSS_SpeciesData_Feb2025_noCEE.xlsx")

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
NSS_sp<-NSS_sp %>% dplyr::filter(!NSS_sp$DB_Assemblage_ID == "NA",
                                 !NSS_sp$Start_date_CE == "NA")

name<-"NSS_sp"
summarise_df(NSS_sp, name)

NSS_sp<-NSS_sp %>%
  dplyr::filter(!Rural_urban_or_neither == "neither",
                !Recovery == "Unknown")

NSS_sp<-NSS_sp %>%
  dplyr::filter(!Recovery == "Unknown")

NSS_sp %>%
  group_by(Region) %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  distinct(DB_Assemblage_ID) %>%
  tally()

NSS_sp %>%
  distinct(GBIF_species) %>%
  tally()

NSS_sp %>%
  dplyr::summarise(total_counts = sum(NISP))

props_grouped <- NSS_sp %>%
  dplyr::group_by(DB_Assemblage_ID, Recovery, Region, Rural_urban_or_neither, Start_date_CE, Time_range, Trait_LifeHistory) %>%
  dplyr::summarise(total_counts = sum(NISP)) %>%
  ungroup()

#3. Calculate proportion of counts within each site ####
props_grouped_b<-aggregate(total_counts ~ Trait_LifeHistory + DB_Assemblage_ID, data = props_grouped, FUN = sum)
props_grouped_c<-aggregate(total_counts ~ DB_Assemblage_ID, data = props_grouped, FUN = sum)
names(props_grouped_c)[2]<-"Totals"
props_grouped_d<-left_join(props_grouped_b,props_grouped_c)
props_grouped_d$proportion<-props_grouped_d$total_counts/props_grouped_d$Totals
hist(props_grouped_d$proportion[props_grouped_d$Trait_LifeHistory=="Oceanodromous"])
props_grouped<-left_join(props_grouped, props_grouped_d)

props_oceanodromous<-props_grouped %>%
  filter(Trait_LifeHistory == "Oceanodromous") %>% distinct()

dataset_non_oceanodromous <- props_grouped %>%
  distinct(DB_Assemblage_ID, Recovery, Region, Rural_urban_or_neither, Start_date_CE, Time_range, Trait_LifeHistory) %>%
  filter(!DB_Assemblage_ID %in% props_oceanodromous$DB_Assemblage_ID) %>%
  mutate(Trait_LifeHistory = "Oceanodromous", total_counts = 0)

props_grouped <- data.frame(bind_rows(props_oceanodromous, dataset_non_oceanodromous) %>% distinct())
str(props_grouped)

props_grouped<-props_grouped %>%
  dplyr::mutate_at('Start_date_CE', as.numeric,
                   'Time_range', as.numeric,
                   'Recovery', as.factor,
                   'Rural_urban_or_neither', as.factor,
                   'Region', as.factor) 

#4. beta regression (linear model) ####
#remove values with 0 and 1 (as they cannot be estimated using betareg)
tmp<-props_grouped %>% 
  dplyr::filter(proportion > 0.01)
tmp<-tmp %>% 
  dplyr::filter(proportion < 1)

br_M1 <- glmmTMB::glmmTMB(
  proportion ~ Start_date_CE + Rural_urban_or_neither + Region + Recovery + Time_range,
  family = beta_family(),
  data = tmp)
summary(br_M1)
AIC(br_M1)
visreg(br_M1)

drop1(br_M1, test = "Chi")

variables <- c("Start_date_CE", "Rural_urban_or_neither", "Region", "Recovery", "Time_range")
comb <- unlist(lapply(1:length(variables), function(i) combn(variables, i, simplify = FALSE)), recursive = FALSE)
comb <- unlist(lapply(1:length(variables), function(i) {
  combs <- combn(variables, i, simplify = FALSE)
  combs_with_interactions <- lapply(combs, function(comb_set) {
    if (length(comb_set) > 1) {
      interactions <- combn(comb_set, 2, FUN = function(x) paste(x[1], x[2], sep = ":"))
      return(c(comb_set, interactions))
    } else {
      return(comb_set)
    }
  })
  return(combs_with_interactions)
}), recursive = FALSE)

aic_values <- sapply(comb, function(x) {
  formula_str <- paste("proportion ~", paste(x, collapse = " + "))
  formula <- as.formula(formula_str)
  tmp_br <- glmmTMB(formula, family = beta_family(), data = tmp)
  return(AIC(tmp_br))
})

result_table <- data.frame(
  Combination = sapply(comb, function(x) paste(x, collapse = ", ")),
  AIC = aic_values
)

# Sort the results by AIC value
result_table <- result_table[order(result_table$AIC), ]
write_xlsx(result_table,"../SupplementaryFiles/GLM_Oceanodromous_AIC.xlsx")
print(result_table)

best_BR <- glmmTMB::glmmTMB(
  proportion ~ Start_date_CE*Region,
  family = beta_family(),
  data = tmp)
summary(best_BR)
AIC(best_BR)
visreg(best_BR, "Start_date_CE", by="Region",overlay=TRUE, partial=FALSE)
visreg(best_BR, "Start_date_CE", by="Region")
visreg(best_BR)

best_BR <- glmmTMB::glmmTMB(
  proportion ~ Start_date_CE*Region,
  family = beta_family(),
  data = tmp)
summary(best_BR)

BR_AIC_out<-bind_rows(result_table) %>%
  arrange(AIC)

summary_df <- tmp %>%
  dplyr::group_by(Region) %>%  # Replace group_variable with the name of your grouping variable (e.g., Region, Recovery)
  dplyr::summarise(
    mean_proportion = mean(proportion, na.rm = TRUE),
    sd_proportion = sd(proportion, na.rm = TRUE),
    min = min(proportion, na.rm = TRUE),
    max = max(proportion, na.rm = TRUE)) %>% ungroup

#Plot model coefficients
BR_data <- data.frame(Start_date_CE = seq(min(props_grouped$Start_date_CE), 
                                           max(props_grouped$Start_date_CE), length.out = 150),
                       Region = rep(unique(props_grouped$Region),50))
BR_data$predicted <- predict(best_BR, newdata = BR_data, type = "response",se.fit = TRUE)$fit
BR_data$se.fit <- predict(best_BR, newdata = BR_data, type = "response", se.fit = TRUE)$se.fit
BR_data$lower <- BR_data$predicted - 1.96 * BR_data$se.fit
BR_data$upper <- BR_data$predicted + 1.96 * BR_data$se.fit

ggplot(BR_data, aes(x = Start_date_CE, y = predicted, color = Region, ymin = lower, ymax = upper)) +
   geom_line() +
   geom_point() +
   geom_ribbon(aes(ymin = lower, ymax = upper, fill = Region), alpha = 0.2) + 
   labs(x = "Year", y = "Propotion of oceanodromous", color = "Region", fill = "Region",
        title = "Beta Regression") +
  ylim(0.1, 1.0) +
  xlim(1,1850) +
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) + theme_publish() + theme(legend.position = "none") 

BR1<-ggplot(BR_data, aes(x = Start_date_CE, y = predicted, color = Region, ymin = lower, ymax = upper)) +
  geom_line(linewidth=2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Region),color=NA, alpha = 0.2) + 
  labs(x = "Year", y = "Propotion of oceanodromous NISP", color = "Region", fill = "Region",
       title = "Beta Regression") +
  ylim(0.1, 1.0) +
  xlim(1,1850) +
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) + theme_publish() + theme(legend.position = "none") 

BR1<-BR1 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1600, fill = "skyblue", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 1700, xmax = 1850, fill = "darkred", alpha = 0.15)
png("../SupplementaryFiles/BR_Prop_Oceanodromous.png", width = 5.5, height = 4.5, units = "in", res = 600)
BR1 + theme(strip.text = element_text(size=10)) +
  annotate("text", x = 1375, y = 0.1, label = "Little Ice Age", color = "grey5", size = 3, 
           family="Arial", fontface="italic") +
  annotate("text", x = 999, y = 0.1, label = "MCA", color = "grey5", size = 3,family="Arial", fontface="italic") +
  #  annotate("text", x = 280, y = 7, label = "RWP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 420, y = 0.1, label = "DACP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 510, y = 0.99, label = "**", color = "grey1", size = 5.2,family="Arial", fontface="italic")
dev.off()

BetaReg<-BR1 + theme(strip.text = element_text(size=10)) +
  annotate("text", x = 1375, y = 0.05, label = "Little Ice Age", color = "grey5", size = 3, 
           family="Arial", fontface="italic") +
  annotate("text", x = 999, y = 0.05, label = "MCA", color = "grey5", size = 3,family="Arial", fontface="italic") +
  #  annotate("text", x = 280, y = 7, label = "RWP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 420, y = 0.05, label = "DACP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 510, y = 0.99, label = "**", color = "grey1", size = 5.2,family="Arial", fontface="italic")

# GAM (nonlinear model) ####
props_grouped$Recovery<-as.factor(props_grouped$Recovery)
props_grouped$Region<-as.factor(props_grouped$Region)
props_grouped$Rural_urban_or_neither<-as.factor(props_grouped$Rural_urban_or_neither)
model_results <- data.frame(
  model_formula = character(0),
  AIC_value = numeric(0)
)

gam_M1.1 <- mgcv::gam(
  proportion ~ s(Start_date_CE),
  data = props_grouped, method="REML")
AIC(gam_M1.1)
model_results <- rbind(model_results, data.frame(
  model_formula = "proportion ~ s(Start_date_CE)",
  AIC_value = AIC(gam_M1.1)
))

gam_M1.2 <- mgcv::gam(
  proportion ~ s(Start_date_CE) + Region,
  data = props_grouped, method="REML")
AIC(gam_M1.2)
model_results <- rbind(model_results, data.frame(
  model_formula = "proportion ~ s(Start_date_CE) + Region",
  AIC_value = AIC(gam_M1.2)
))

gam_M1.3 <- mgcv::gam(
  proportion ~ Region + s(Start_date_CE, by = Region),
  data = props_grouped, method="REML")
AIC(gam_M1.3)
model_results <- rbind(model_results, data.frame(
  model_formula = "proportion ~ Region + s(Start_date_CE, by = Region)",
  AIC_value = AIC(gam_M1.3)
))
summary(gam_M1.3)  
par(mfrow = c(2, 2))
plot(gam_M1.3)

gam_M1.4 <- mgcv::gam(
  proportion ~ Region + s(Start_date_CE, by = Region) + Rural_urban_or_neither + Recovery,
  data = props_grouped, method="REML")
AIC(gam_M1.4)
summary(gam_M1.4)
model_results <- rbind(model_results, data.frame(
  model_formula = "proportion ~ Region + s(Start_date_CE, by = Region) + Rural_urban_or_neither + Recovery",
  AIC_value = AIC(gam_M1.4)
))

gam_M1.5 <- mgcv::gam(
  proportion ~ Region + s(Start_date_CE, by = Region) + Rural_urban_or_neither + Recovery + s(Time_range),
  data = props_grouped, method="REML")
AIC(gam_M1.5)
summary(gam_M1.5)
model_results <- rbind(model_results, data.frame(
  model_formula = "proportion ~ Region + s(Start_date_CE, by = Region) + Rural_urban_or_neither + Recovery + Time_range",
  AIC_value = AIC(gam_M1.5)
))

gam_M1.6 <- mgcv::gam(
  proportion ~ Region + s(Start_date_CE, by = Region) + Rural_urban_or_neither + Recovery,
  data = props_grouped, method="REML")
AIC(gam_M1.6)
model_results <- rbind(model_results, data.frame(
  model_formula = "proportion ~ Region + s(Start_date_CE, by = Region) + Rural_urban_or_neither + Recovery",
  AIC_value = AIC(gam_M1.6)
))

variables <- c("s(Start_date_CE)", "Rural_urban_or_neither", "Region", "Recovery", "s(Time_range)")
comb <- unlist(lapply(1:length(variables), function(i) combn(variables, i, simplify = FALSE)), recursive = FALSE)

aic_values <- sapply(comb, function(x) {
  formula_str <- paste("proportion ~", paste(x, collapse = " + "))
  formula <- as.formula(formula_str)
  tmp_gam <- mgcv::gam(formula, data = props_grouped, method="REML")
  return(AIC(tmp_gam))
})

result_table_gams <- data.frame(
  Combination = sapply(comb, function(x) paste(x, collapse = ", ")),
  AIC = aic_values
)

result_table_gams$Combination<-gsub(", ","+ ",result_table_gams$Combination)
names(result_table_gams)<-c("model_formula","AIC_value")
GAMS_AIC_out<-bind_rows(result_table_gams,model_results) %>%
  arrange(AIC_value)
write_xlsx(GAMS_AIC_out,"../SupplementaryFiles/GAMS_AIC_out.xlsx")

gam_best <- mgcv::gam(
  proportion ~ Region + s(Start_date_CE, by = Region) + Rural_urban_or_neither + Recovery,
  data = props_grouped, select=TRUE
  )
AIC(gam_best)
summary(gam_best)

gam_best <- mgcv::gam(
  proportion ~ Region + s(Start_date_CE, by = Region) + Rural_urban_or_neither + Recovery,
  data = props_grouped, select=TRUE
)
AIC(gam_best)
summary(gam_best)

gam_data <- data.frame(Start_date_CE = seq(min(props_grouped$Start_date_CE), 
                                          max(props_grouped$Start_date_CE), length.out = 300),
                      Region = rep(unique(props_grouped$Region),100),
                      Recovery = rep("both",300),
                      Rural_urban_or_neither = rep("rural",300))

gam_data$predicted <- predict(gam_best, newdata = gam_data, type = "response",se.fit = TRUE)$fit
gam_data$se.fit <- predict(gam_best, newdata = gam_data, type = "response", se.fit = TRUE)$se.fit
gam_data$lower <- gam_data$predicted - 1.96 * gam_data$se.fit
gam_data$upper <- gam_data$predicted + 1.96 * gam_data$se.fit
gam_data$all<-"Generalized Additive Model"
summary(gam_best)
plot(gam_best)

ggplot(gam_data, aes(x = Start_date_CE, y = predicted, color = Region, ymin = lower, ymax = upper)) +
  geom_line(linewidth=2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Region),color=NA, alpha = 0.2) + 
  labs(x = "Year", y = "Propotion of oceanodromous", color = "Region", fill = "Region",
       title = "Generalized Additive Model") +
  ylim(0.01, 1.1) +
  xlim(1,1850) +
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) + theme_publish() + theme(legend.position = "none")

M2<-ggplot(gam_data, aes(x = Start_date_CE, y = predicted, color = Region, ymin = lower, ymax = upper)) +
  geom_line(linewidth=2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Region),color=NA, alpha = 0.2) + 
  labs(x = "Year", y = "Propotion of oceanodromous NISP", color = "Region", fill = "Region",
       title = "Generalized Additive Model") +
  ylim(0.01, 1.1) +
  xlim(1,1850) +
  scale_color_manual(values=DB[2:6]) +
  scale_fill_manual(values=DB[2:6]) + theme_publish() + theme(legend.position = "none")

###Add environmental periods
png("../SupplementaryFiles/GAM_Prop_Oceanodromous.png", width = 5.5, height = 4.5, units = "in", res = 600)
Gam1<-M2 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1600, fill = "skyblue", alpha = 0.15) %>%
  add_shaded_rectangle(xmin = 1700, xmax = 1850, fill = "darkred", alpha = 0.15)
Gam1 + theme(strip.text = element_text(size=10)) +
  annotate("text", x = 1375, y = 0.05, label = "Little Ice Age", color = "grey5", size = 3, 
           family="Arial", fontface="italic") +
  annotate("text", x = 999, y = 0.05, label = "MCA", color = "grey5", size = 3,family="Arial", fontface="italic") +
  #  annotate("text", x = 280, y = 7, label = "RWP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 420, y = 0.05, label = "DACP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 510, y = 0.99, label = "**", color = "grey1", size = 5.2,family="Arial", fontface="italic")
dev.off()

GAM<-Gam1 + theme(strip.text = element_text(size=10)) +
  annotate("text", x = 1375, y = 0.05, label = "Little Ice Age", color = "grey5", size = 3, 
           family="Arial", fontface="italic") +
  annotate("text", x = 999, y = 0.05, label = "MCA", color = "grey5", size = 3,family="Arial", fontface="italic") +
  #  annotate("text", x = 280, y = 7, label = "RWP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 420, y = 0.05, label = "DACP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 510, y = 0.99, label = "**", color = "grey1", size = 5.2,family="Arial", fontface="italic")

###Plot both together for main manuscript
png("../SupplementaryFiles/BetaReg_GAM_Figure.png", width = 9.5, height = 4.5, units = "in", res = 600)
grid.arrange(BetaReg,GAM, nrow=1)
dev.off()

#Statistically compare both models ####
log_likelihood_model<-logLik(best_BR) 
deviance_best_BR <- -2 * as.numeric(logLik(best_BR))
null_model <- update(best_BR, . ~ 1) 
null_log<-logLik(null_model)
null_deviance <- -2 * as.numeric(logLik(null_model))
r2_ml <- 1 - (log_likelihood_model / null_log)
print(r2_ml)

###For main figure
library(tidyr)
NSS_sp_long <- tmp %>%
  dplyr::select(proportion, Region, Recovery, Rural_urban_or_neither) %>%
  tidyr::gather(key = "group_type", value = "group", Region, Recovery, Rural_urban_or_neither)
NSS_sp_long$group_type<-gsub("Recovery","Recovery method",NSS_sp_long$group_type)
NSS_sp_long$group_type<-gsub("Rural_urban_or_neither","Site type",NSS_sp_long$group_type)

remove_outliers <- function(x) {
  # Calculate IQR
  IQR_x <- IQR(x)
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  
  # Identify outliers
  lower_bound <- Q1 - 1.5 * IQR_x
  upper_bound <- Q3 + 1.5 * IQR_x
  
  # Return values that are within bounds
  x[x >= lower_bound & x <= upper_bound]
}
NSS_sp_long <- NSS_sp_long %>%
  mutate(proportion = remove_outliers(proportion)) %>%
  filter(!is.na(proportion))

groups<-ggplot(NSS_sp_long, aes(x = group, y = proportion, color=group, fill=group)) +
  geom_jitter(alpha=0.35) +
  geom_boxplot(alpha=0.75) +
  facet_wrap(~ group_type, scales = "free") +
  scale_color_manual(values=DB2[2:12]) +
  scale_fill_manual(values=DB2[2:12]) + theme_publish() + theme(legend.position = "none") +
  labs(x = "Group", y = "Proportion", title = "")

png("../SupplementaryFiles/BetaReg_GAM_Figure.png", width = 9.5, height = 4.5, units = "in", res = 600)
grid.arrange(BetaReg,GAM, nrow=1)
dev.off()


png("../SupplementaryFiles/BetaReg_GAM_Figure.png", width = 8.5, height = 7, units = "in", res = 600)
layout_matrix <- matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE)
grid.arrange(groups, BetaReg,GAM, layout_matrix = layout_matrix)
dev.off()
