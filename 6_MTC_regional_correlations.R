#Cross correlations of regional MTC ####
#Author: D L Buss
library(magrittr)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(corrplot)
library(tidyverse)

#1. Felling Dates as proxy of human population size (1200 onwards)
DB<-c("#F8F1E9","#FE994F","#6A8A73","#14517B","#8e7d69", "Gold3")
DB2<-c("#14517B","#FE994F","#6A8A73","#8e7d69","#F8F1E9", "Gold3")
DB3<-c("#6A8A73","#FE994F","#14517B","grey5","darkred","#006ba4","grey18","grey40","#D6604D","#2B547E",
       "gold2","#F39B7FFF",
       "#91D1C2FF","#8e7d69","#0072B2",
       "#6A8A73","#FE994F","#8491B4FF","#14517B","#FE994F","#6A8A73","#8e7d69","#F8F1E9", "Gold3")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

###
final<-read_xlsx("../Data/MTC_output.xlsx", col_names=T)
#MAKE matrix
final<-dcast(final[final$sizebins=="100",], TimePeriod ~ plot_names , value.var = "MTC")
final<-final[complete.cases(final),]

# Get all variable pairs
var_pairs <- combn(names(final[,c(2:5)]), 2, simplify = FALSE)
names(final)[3]<-"NSS-DB 100-yr"
cor_matrix <- cor(final[,c(2:5)])
corrplot(cor_matrix, method = "number", type = "lower",
         tl.col = "black", tl.cex = 0.75,
         col = colorRampPalette(c("navyblue","grey80","darkred"))(200))

tmo<-final %>% dplyr::mutate(TimePeriod = as.factor(TimePeriod)) %>%
  melt() 
tmo$TimePeriod<-as.numeric(tmo$TimePeriod)
ggplot(tmo,aes(x = TimePeriod, y = value, color = variable)) +
  geom_line(linewidth=1.2, alpha=0.9) + 
  facet_wrap(~tmo$variable, ncol=1, scale="free") + theme_publish() + 
  scale_color_manual(values=DB)

# Sample data: Numeric matrix with multiple time series
# Get the number of columns
num_vars <- ncol(final[,c(2:5)])

# Loop through each pair of variables
par(mfrow = c(2, 3))  # Set plot layout
for (i in 1:(num_vars - 1)) {
  for (j in (i + 1):num_vars) {
    ccf(final[, i], final[, j], 
        main = paste("CCF:", colnames(final[,c(2:5)])[i], "vs", colnames(final[,c(2:5)])[j]),
        na.action = na.omit)
  }
}
par(mfrow = c(1, 1))

# Get all variable pairs
var_pairs <- combn(names(final[,c(2:5)]), 2, simplify = FALSE)
finalb<-final[,c(2:5)]
# Calculate CCF for all pairs

ccf_results <- map(var_pairs, ~ {
  ccf_data <- get_ccf(df[[.x[1]]], df[[.x[2]]])
  ccf_data$Var1 <- .x[1]
  ccf_data$Var2 <- .x[2]
  return(ccf_data)
})

# Combine into one data frame
ccf_df <- bind_rows(ccf_results)

# Plot the CCFs using ggplot2
ggplot(ccf_df, aes(x = lag, y = ccf)) +
  geom_line() +
  facet_grid(Var1 ~ Var2, scales = "free") +
  theme_minimal() +
  labs(title = "Cross-Correlation Functions", 
       x = "Lag", y = "CCF")



# Combine into one data frame
ccf_df <- bind_rows(ccf_results)
conf_level <- 1.96 / sqrt(15)

# Plot the CCFs using ggplot2
ggplot(ccf_df, aes(x = lag, y = ccf)) +
  geom_line() +
  geom_hline(yintercept = c(-conf_level, conf_level), 
             linetype = "dashed", color = "blue") +
  facet_grid(Var1 ~ Var2) +
  theme_minimal() +
  xlim(-6, 6) +
  labs(title = "Cross-Correlations of Env. Variables", 
       x = "Lag", y = "CCF")

##png("../SupplementaryFiles/SuppFigureCorr.png", width = 10.5, height = 10.5, units = "in", res = 600)
ggplot(ccf_df, aes(x = lag, y = ccf)) +
  geom_line() +
  geom_hline(yintercept = c(-conf_level, conf_level), 
             linetype = "dashed", color = "blue") +
  facet_grid(Var1 ~ Var2) +
  theme_minimal() +
  xlim(-4.5, 4.5) +
  labs(title = "Cross-Correlations of Env. Variables (100-yr)", 
       x = "Lag", y = "CCF") +
  theme(strip.text = element_text(size=7))
#dev.off()

cor_matrix <- cor(envs_all_new2[,c(2:9)])
#png("../SuppFigureCorr2.png", width = 7.5, height = 7.5, units = "in", res = 600)
corrplot(cor_matrix, method = "number", type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
#dev.off()



###with sig.
p_values <- cor.mtest(cor_matrix)
#png("../SupplementaryFiles/SuppFigureCorr2.png", width = 6.5, height = 6.5, units = "in", res = 600)
corrplot(cor_matrix, method = "number", 
         p.mat = p_values,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
#dev.off()

final2<-final %>%
  dplyr::filter(level == "Regional") %>%
  dplyr::filter(Region == "Britain & Ireland") %>%
  dplyr::filter(sizebins == "100") %>%
  dplyr::select(TimePeriod,MTC,Region) %>% distinct()
final2$Proxy<-"MTC"
names(final2)<-c("Year","Value","Region","Proxy")
final2$Region<-"*MTC (NSS-DB)"

envs<-bind_rows(df_SST[,c(1,4,5)])
envs$Region<-"Norwegian Sea SST (Berner et al. 2011)"
names(envs)<-c("Year","Value","Proxy","Region")
envs$Proxy<-"Subpolar SST"
envs_all<-bind_rows(envs,final2)

envs<-overtime
envs$Region<-"Felling counts (Ljungqvist et al. 2022)"
names(envs)<-c("Year","Value","Proxy","Region")
envs$Proxy<-"Pop. Size Proxy"
envs$Year<-envs$Year %>% as.numeric()
envs$Value2<-envs$Value
envs$Value<-NA
envs_all<-bind_rows(envs_all,envs)

envs<-bind_rows(df_SST_Af[,c(1,5)])
envs$Proxy<-"N.Africa SST"
envs$Region<-"N.Africa SST (McGregor et al. 2007)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

envs<-bind_rows(df_SST_SP[,c(1,4)])
envs$Proxy<-"Subpolar SST."
envs$Region<-"Subpolar SST (Miettinen et al. 2012)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

envs<-bind_rows(df_SST_SP2[,c(1,4)])
envs$Proxy<-"Polar SST"
envs$Region<-"Polar SST (Sicre et al. 2011)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

envs<-bind_rows(df_XRAY[,c(1,4)])
envs$Proxy<-"d18O"
envs$Region<-"N.Atl Precip. (Tiljander et al. 2003)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

envs<-bind_rows(df_Forams[,c(1,3)])
envs$Proxy<-"Ice.Vol."
envs$Region<-"N.Atl Ice Volume (Keigwin 1996)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

envs<-bind_rows(df_CS[,c(1,4)])
envs$Proxy<-"d18O."
envs$Region<-"N.Atl Precip. Rates (Fohlmeister et al. 2012)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

stats_sums<-envs_all %>%
  dplyr::group_by(Region, Proxy) %>%
  dplyr::summarise(
    min = min(Value, na.rm=T),
    max = max(Value, na.rm=T),
    mean = mean(Value, na.rm=T),
    sd = sd(Value, na.rm=T),
    range = diff(range(Value, na.rm=T)),
    count = sum(!is.na(Value))
  ) %>% ungroup()
stats_sums

stats_sumsb<-envs_all %>%
  dplyr::filter(level == "Regional") %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(
    min = min(Value, na.rm=T),
    max = max(Value, na.rm=T),
    mean = mean(Value, na.rm=T),
    sd = sd(Value, na.rm=T),
    range = diff(range(Value, na.rm=T)),
    count = sum(!is.na(Value))
  ) %>% ungroup()
stats_sums
