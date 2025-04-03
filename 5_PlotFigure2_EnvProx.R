#Plot Figure 2 for NSS - Env Proxies ####
#Author: D L Buss
#Load libraries/data ####
library("remotes")
remotes::install_github("nickmckay/geoChronR")
library(geoChronR)
library(magrittr)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(corrplot)
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

#1. Felling Dates as proxy of human population size (1200 onwards) ####
df_fell<-read_xlsx("../Data/Env_proxies/Felling_Dates/Felling dates_Ljungqvist et al 2022_2.xlsx", col_names=T)
df_fell_b<-read_xlsx("../Data/Env_proxies/Felling_Dates/Fellingcounts_Ljungqvist et al 2022_2_avg.xlsx", col_names=T)
df_fell$all<-1
str(df_fell)
overtime<-aggregate(all ~ Year, df_fell, FUN=sum)
overtime$all[overtime$all < 4]<-NA
overtime$title<-"Proxy of demographic change from Ljungqvist et al. (2022)"
ggplot(overtime, aes(Year, all)) + 
  geom_point(color="#FE994F", alpha=0.8) + 
  labs(y="Counts of trees felled", tag="B)") + 
  geom_smooth(color="grey30",level = 0.95) + 
  theme_publish() + facet_wrap(~title)
#write_xlsx(overtime, "~/Documents/4-Oceans/3.NSS/MANUSCRIPT_Oct2024/1_Data/EnvironmentalProxies/Fellingcounts_Ljungqvist et al 2022_2.xlsx")
B<-ggplot(overtime, aes(Year, all)) + 
  geom_point(color="#FE994F", alpha=0.8) + 
  labs(y="Counts of trees felled", tag="B)") + 
  geom_smooth(color="grey30",level = 0.95) + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850)
P4_B<-B %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)

#2. SST_Aug ####
df_SST<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Berner, K.S., N. Koç, F. Godtliebsen, and D. Divine. 2011..xlsx", sheet=2)
head(df_SST)
str(df_SST)
names(df_SST)[1]<-"Year"
names(df_SST)[3]<-"SST_25yr_av"
names(df_SST)[4]<-"SST_100yr_av"
df_SST$title<-"Proxy of Sea Surface Tempertaure from Berner et al. (2011)"
#100-yr rolling average
ggplot(df_SST, aes(Year, SST_100yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.8) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish() 

C<-ggplot(df_SST, aes(Year, SST_100yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.95, size=1.9) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish()  + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850)
P4_C<-C %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) + ylim(9.75,11.25)
P4_C

#3. Oxygen isotopes of Cave_Stalagamites ####
df_CS<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Fohlmeister, J., Schröder-Ritzrau, A., Scholz...2012.xlsx",sheet=2)
head(df_CS)
str(df_CS)
names(df_CS)[1]<-"Year"
names(df_CS)[3]<-"CS_25yr_av"
names(df_CS)[4]<-"CS_100yr_av"
df_CS$title<-"Proxy of precipitation rates in western Europe from Fohlmeister et al. (2012)"
oxygen_lab<-expression(paste(delta^18, "O (\u2030)"))
df_CS$CS_100yr_av<-1-df_CS$CS_100yr_av
df_CS$CS_25yr_av<-1-df_CS$CS_25yr_av

#25-yr rolling average
ggplot(df_CS, aes(Year, CS_25yr_av)) + 
  geom_point(color="#14517B", alpha=0.8) + 
  labs(y="1 - (Oxygen isotopes of cave stalagmites)") + 
  theme_publish() 
#100-yr rolling average
ggplot(df_CS, aes(Year, CS_100yr_av)) + 
  geom_point(color="#14517B", alpha=0.8) + 
  labs(y="1 - (Oxygen isotopes of cave stalagmites)") + 
  theme_publish() 
D<-ggplot(df_CS, aes(Year, CS_100yr_av)) + 
  geom_point(color="#14517B", alpha=0.95, size=1.9) + 
  labs(y=oxygen_lab) +
  theme_publish()  + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850) + 
  ylim(5.7,7.4)
P4_D<-D %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) 
P4_D

#4. Oxygen Isotopes of Sargasso Sea Forams ####
df_Forams<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Keigwin_1996.xlsx", sheet=2)
head(df_Forams)
str(df_Forams)
df_Forams$`25yr_avg`<-as.numeric(df_Forams$`25yr_avg`)
df_Forams$`100yr_avg`<-as.numeric(df_Forams$`100yr_avg`)
names(df_Forams)[1]<-"Year"
names(df_Forams)[2]<-"OI_25yr_av"
names(df_Forams)[3]<-"OI_100yr_av"
df_Forams$title<-"Proxy of Seasonal Ice Volumes in subpolar N.Atlantic from Keigwin (1996)"
#25-yr rolling average
ggplot(df_Forams, aes(Year, OI_25yr_av)) + 
  geom_point(color="#14517B", alpha=0.8) + 
  labs(y="1 - (Oxygen isotopes of sargasso sea foraminifera)") + 
  theme_publish() 
#100-yr rolling average
ggplot(df_Forams, aes(Year, OI_100yr_av)) + 
  geom_point(color="#14517B", alpha=0.8) + 
  labs(y="1 - (Oxygen isotopes of sargasso sea foraminifera)") + 
  theme_publish()

E<-ggplot(df_Forams, aes(Year, OI_100yr_av)) + 
  geom_point(color="#14517B", alpha=0.95, size=1.9) + 
  labs(y=oxygen_lab) +
  theme_publish()  + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850) + 
  ylim(-0.8,0.01)
P4_E<-E %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) 
P4_E

#5. Sea Surface Temp - Africa ####
df_SST_Af<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_McGregor, H.V., M. Dima, H.W. Fischer, and S. Mulitza 2007.xlsx", sheet=2)
head(df_SST_Af)
str(df_SST_Af)
df_SST_Af$`25yr_avg`<-as.numeric(df_SST_Af$`25yr_avg`)
df_SST_Af$`100yr_avg`<-as.numeric(df_SST_Af$`100yr_avg`)
names(df_SST_Af)[1]<-"Year"
names(df_SST_Af)[4]<-"SST_25yr_av"
names(df_SST_Af)[5]<-"SST_100yr_av"
df_SST_Af$title<-"Proxy of Sea Surface Tempertaure for N.Africa Atlantic from McGregor et al. (2011)"

#25-yr rolling average
ggplot(df_SST_Af, aes(Year, SST_25yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.8) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish() 
#100-yr rolling average
ggplot(df_SST_Af, aes(Year, SST_100yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.8) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish() 
PF<-ggplot(df_SST_Af, aes(Year, SST_100yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.95, size=1.9) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish()  + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850) + 
  ylim(19.0,20.4)
P4_F<-PF %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) 
P4_F

#6. Sea Surface Temp - Subpolar ####
df_SST_SP<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Miettinen, A., D. Divine, N. Koç¸ F. Godtliebsen, and I.R. Hall. 2012..xlsx", sheet=2)
head(df_SST_SP)
str(df_SST_SP)
df_SST_SP$`25yr_avg`<-as.numeric(df_SST_SP$`25yr_avg`)
df_SST_SP$`100yr_avg`<-as.numeric(df_SST_SP$`100yr_avg`)
names(df_SST_SP)[1]<-"Year"
names(df_SST_SP)[3]<-"SST_25yr_av"
names(df_SST_SP)[4]<-"SST_100yr_av"
df_SST_SP$title<-"Proxy of Sea Surface Tempertaure for subpolar N.Atlantic from Miettinen et al. (2012)"

#25-yr rolling average
ggplot(df_SST_SP, aes(Year, SST_25yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.8) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish() 
#100-yr rolling average
ggplot(df_SST_SP, aes(Year, SST_100yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.8) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish() 
PG<-ggplot(df_SST_SP, aes(Year, SST_100yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.95, size=1.9) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish()  + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850) + 
  ylim(12.2,13.8)
P4_G<-PG %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) 
P4_G

#7. SST - Subpolar 2 ####
df_SST_SP2<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Sicre, M.-A., I.R. Hall, J. Mignot, M. Khodri, U. Ezat.xlsx", sheet=2)
head(df_SST_SP2)
str(df_SST_SP2)
df_SST_SP2$`25yr_Rolling_Avg`<-as.numeric(df_SST_SP2$`25yr_Rolling_Avg`)
df_SST_SP2$`100yr_Rolling_Avg`<-as.numeric(df_SST_SP2$`100yr_Rolling_Avg`)
names(df_SST_SP2)[1]<-"Year"
names(df_SST_SP2)[3]<-"SST_25yr_av"
names(df_SST_SP2)[4]<-"SST_100yr_av"
df_SST_SP2$title<-"Proxy of Sea Surface Tempertaure for subpolar N.Atlantic from Sicre et al. (2011)"
#25-yr rolling average
ggplot(df_SST_SP2, aes(Year, SST_25yr_av)) + 
  geom_point(color="#8460A3", alpha=0.8) + 
  labs(y="Sea surface temperatures") + 
  theme_publish() 
#100-yr rolling average
ggplot(df_SST_SP2, aes(Year, SST_100yr_av)) + 
  geom_point(color="#8460A3", alpha=0.8) + 
  labs(y="Sea surface temperatures") + 
  theme_publish() 

PH<-ggplot(df_SST_SP2, aes(Year, SST_100yr_av)) + 
  geom_point(color="#6A8A73", alpha=0.95, size=1.9) + 
  labs(y=expression("SST ("*~degree*C*")")) +
  theme_publish()  + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850) + 
  ylim(6.2,9.8)
P4_H<-PH %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) 
P4_H

#8. XRay Fluorescence as an indicator of biogenic carbonate production, increased = warmer ####
df_XRAY<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_AVG_Tiljander, M., M. Saarnisto, A.E.K. Ojala, and T. Saarinen. 2003..xlsx", sheet=2)
head(df_XRAY)
str(df_XRAY)
df_XRAY$`25yr_Rolling_Avg`<-as.numeric(df_XRAY$`25yr_Rolling_Avg`)
df_XRAY$`100yr_Rolling_Avg`<-as.numeric(df_XRAY$`100yr_Rolling_Avg`)
names(df_XRAY)[1]<-"Year"
names(df_XRAY)[3]<-"XR_25yr_av"
names(df_XRAY)[4]<-"XR_100yr_av"
df_XRAY$title<-"Proxy for subpolar North Atlantic precipiation rates from Tiljander et al. (2003)"
#25-yr rolling average
ggplot(df_XRAY, aes(Year, XR_25yr_av)) + 
  geom_point(color="#8460A3", alpha=0.8) + 
  labs(y="X-Ray fluorescence (Proxy of oceanic productivity)") + 
  theme_publish() 
#100-yr rolling average
ggplot(df_XRAY, aes(Year, XR_100yr_av)) + 
  geom_point(color="#F5E9DC", alpha=0.8) + 
  labs(y="X-Ray fluorescence (Proxy of oceanic productivity)") + 
  theme_publish() 

PI<-ggplot(df_XRAY, aes(Year, XR_100yr_av)) + 
  geom_point(color="#8460A3", alpha=0.95, size=1.9) + 
  labs(y=expression("Precipitation rates")) +
  theme_publish()  + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850) + 
  ylim(55,105)
P4_I<-PI %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) 
P4_I

#9. Pages 2K ####
df_pages<-read_xlsx("../Data/Env_proxies/Env_Datasets/25_100_Pages2K.xlsx", sheet=1)
names(df_pages)[1]<-"Year"
names(df_pages)[3]<-"Europe2K_25yr_av"
names(df_pages)[4]<-"Europe2K_100yr_av"
names(df_pages)[12]<-"Arctic2K_25yr_av"
names(df_pages)[13]<-"Arctic2K_100yr_av"
df_pages$title<-"Proxy of Global Tempertaure Change from Pages 2K (2017) - Arctic"
#Arctic - 100 years
ggplot(df_pages, aes(Year, Arctic2K_100yr_av)) + 
  geom_point(color = "skyblue", alpha=0.8) + 
  labs(y=expression("Temp. ("*~degree*C*")")) + 
  theme_publish() 

#Europe - 100 years
ggplot(df_pages, aes(Year, Europe2K_100yr_av)) + 
  geom_point(color="#D6604D", alpha=0.8) + 
  labs(y=expression("Temp. ("*~degree*C*")")) +
  theme_publish() 

J<-ggplot(df_pages, aes(Year, Arctic2K_100yr_av)) + 
  geom_point(color = "skyblue", alpha=0.8) + 
  labs(y=expression("Temp. ("*~degree*C*")")) + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850)
P4_J<-J %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) + ylim(-1.5,0)
P4_J

K<-ggplot(df_pages, aes(Year, Europe2K_100yr_av)) + 
  geom_point(color = "#A2C8EC", alpha=0.8) + 
  labs(y=expression("Temp. ("*~degree*C*")")) + 
  theme_publish() + facet_wrap(~title) + xlim(50,1850)
P4_K<-K %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25) + ylim(-0.6,0.6)
P4_K

###Plot MTC with env correlates
final<-read_xlsx("../Data/MTC_output.xlsx", col_names=T)
P4<-ggplot(final[final$sizebins=="25" &
                   final$level=="Global",], aes(TimePeriod, MTC)) + geom_line(aes(y=MTC), colour = "grey45", linewidth = 1.2, alpha=0.9) +                          # Main line
  geom_line(aes(y = MTC - SDMTC), linetype = "dashed", color = "grey79", size=0.5) +                          # Main line
  geom_line(aes(y = MTC + SDMTC), linetype = "dashed", color = "grey79", size=0.5) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC), fill = "grey79", size = 1.2, alpha=0.3) + 
  theme_publish() +
  ggtitle("MTC, 200 simulations") + 
  labs(tag="A)") +
  facet_wrap(~plot_names, scales="free", ncol=1) + xlim(50,1850)


P4 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)


P4_A<-P4 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)

plots<-list(P4_A,P4_B,P4_C,P4_D,P4_E,P4_F,
         P4_G,P4_H,P4_I, P4_J,P4_K)
plots[[8]]
#for (i in seq_along(plots)) {
#  ggsave(filename = paste0("plot_", i, ".png"), plot = plots[[i]], width=6, height=4.5, units="in", dpi=600, device = "png")
#}

final$Region<-recode(final$Run,
                    `All species; 100-yr time bins (Global Data (no CEE), mean + SD)` = "Global",
                    `All species; 100-yr time bins (Britain & Ireland, mean + SD)` = "Britain & Ireland",
                    `All species; 100-yr time bins (Western Europe, mean + SD)` = "Western Europe",
                    `All species; 100-yr time bins (Scandinavia, mean + SD)` = "Scandinavia",
                    `All species; 25-yr time bins (Global Data (no CEE), mean + SD)` = "Global",
                    `All species; 25-yr time bins (Britain & Ireland, mean + SD)` = "Britain & Ireland",
                    `All species; 25-yr time bins (Western Europe, mean + SD)` = "Western Europe",
                    `All species; 25-yr time bins (Scandinavia, mean + SD)` = "Scandinavia")

###Plot only MTC (NSS DB without CEE)
A_v1<-ggplot(final[final$sizebins=="100" &
               final$level=="Global",], 
       aes(TimePeriod, MTC, color=Region)) + geom_line(aes(y=MTC, colour = Region, linetype=Region), linewidth = 1.8, alpha=0.9) +                          # Main line
  #  geom_line(aes(y = MTC - SDMTC, colour=Run), linetype = "dashed", size=0.5) +                          # Main line
  # geom_line(aes(y = MTC + SDMTC, colour=Run), linetype = "dashed", size=0.5) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill=Region, color=Region), color=NA, size = 1.2, alpha=0.1) + 
  theme_bw() +
  ggtitle("European temperature variability (inferred from MTC)") + 
  labs(tag="A)", y="Median Temp. Catch (MTC)", x="Year") + xlim(50,1850) + theme(legend.position = "none") + 
  scale_color_manual(values="Grey15") +
  scale_fill_manual(values="Grey15") + 
  scale_linetype_manual(values=c("solid","solid","solid")) +
  annotate("text", x = 1325, y = 7, label = "Little Ice Age", color = "grey5", size = 3, 
           family="Arial", fontface="italic") +
  annotate("text", x = 999, y = 7, label = "MCA", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 280, y = 7, label = "RWP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 420, y = 7, label = "DACP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 510, y = 14.97, label = "**", color = "grey1", size = 5.2,family="Arial", fontface="italic") + ylim(6.95,15)

A_v1 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)

Panel_A_MTC<-A_v1 %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)


F2A<-ggplot(final[final$sizebins=="100" &
                   final$level=="Regional",], 
       aes(TimePeriod, MTC, color=Region)) + geom_line(aes(y=MTC, colour = Region, linetype=Region), linewidth = 1.2, alpha=0.9) +                          # Main line
#  geom_line(aes(y = MTC - SDMTC, colour=Run), linetype = "dashed", size=0.5) +                          # Main line
# geom_line(aes(y = MTC + SDMTC, colour=Run), linetype = "dashed", size=0.5) + 
  geom_ribbon(aes(y=MTC, ymin = MTC - SDMTC, ymax = MTC + SDMTC, fill=Region, color=Region), color=NA, size = 1.2, alpha=0.1) + 
  theme_bw() +
  ggtitle("Regionional temperature variability") + 
  labs(tag="A)", y="Median Temp. Catch (MTC)", x="Year") + xlim(50,1850) + theme(legend.position = "bottom") + 
  scale_color_manual(values=DB2) +
  scale_fill_manual(values=DB2) + 
  scale_linetype_manual(values=c("solid","solid","solid")) +
  annotate("text", x = 1325, y = 7, label = "Little Ice Age", color = "grey5", size = 3, 
           family="Arial", fontface="italic") +
  annotate("text", x = 999, y = 7, label = "MCA", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 280, y = 7, label = "RWP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 420, y = 7, label = "DACP", color = "grey5", size = 3,family="Arial", fontface="italic") +
  annotate("text", x = 510, y = 14.97, label = "**", color = "grey1", size = 5.2,family="Arial", fontface="italic") + ylim(6.95,15)

F2A %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)

Panel_A<-F2A %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)

###Melt all environmental variables together 
final2<-final %>%
  dplyr::filter(level == "Global") %>%
  dplyr::filter(sizebins == "100") %>%
  dplyr::select(TimePeriod,MTC,Region) %>% distinct()
final2$Proxy<-"MTC"
names(final2)<-c("Year","Value","Region","Proxy")
final2$Region<-"*MTC (NSS-DB)"
head(final2)

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

envs<-bind_rows(df_pages[,c(1,4)])
envs$Proxy<-"Eur. Temp."
envs$Region<-"Global Temperature Model (Pages2K Con. 2017)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

envs<-bind_rows(df_pages[,c(1,13)])
envs$Proxy<-"Arctic Temp."
envs$Region<-"Global Temperature Model (Pages2K Con. 2017)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

##Plot
DB3<-c("#6A8A73","#FE994F","#14517B","grey5","darkred","#006ba4","grey18","grey40","#D6604D","#2B547E",
       "gold2","#F39B7FFF",
       "#91D1C2FF","#8e7d69","#0072B2",
       "#6A8A73","#FE994F","#8491B4FF","#14517B","#FE994F","#6A8A73","#8e7d69","#F8F1E9", "Gold3")
DB4<-c("#006BA4","#FF800E","#000000","#C85200","#5F9Ed1","#CFCFCF","#FFBc79","#00887d","#76c0c1","#014d64","#B09C85FF")

Fig2_plot<-ggplot(envs_all,
            aes(Year, Value, color=Region)) + 
  geom_line(aes(y=Value, colour = Region, linetype=Region), linewidth = 1.2, alpha=0.9) +                          # Main line
  theme_bw() +
  ggtitle("") + 
  labs(tag="B)", y="", colour="",linetype="") + xlim(50,1850) + theme(legend.position = "right") + 
  scale_color_manual(values=DB4) +
  scale_fill_manual(values=DB4) + 
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","solid","solid","solid","solid","solid","solid")) +
  scale_color_manual(values=DB4) +
  facet_grid(Proxy ~., scales="free_y") + 
  geom_point(data = envs_all, aes(Year, Value2),color="#FE994F", alpha=0) + 
  geom_smooth(data = envs_all, aes(Year, Value2), color="#FE994F",level = 0.95) 

Fig2_plot<-Fig2_plot +
  theme(
    panel.spacing = unit(0, "lines"),  # Removes space between panels
    plot.margin = margin(0, 0, 0, 0)  # Adjust plot margins if needed
  ) 

Fig2_plot<-Fig2_plot %>% add_shaded_rectangle(xmin = 350, xmax = 550, fill = "skyblue", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 490, xmax = 540, fill = "orange", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 900, xmax = 1100, fill = "darkred", alpha = 0.25) %>%
  add_shaded_rectangle(xmin = 1150, xmax = 1500, fill = "skyblue", alpha = 0.25)
Fig2_plot
Panel_B<-Fig2_plot + theme(strip.text = element_text(size=8))

ggarrange(Panel_A, Panel_B, ncol=2)
ggarrange(Panel_A_MTC, Panel_B, ncol=2)

png("../Figure2_v1.png", width = 11.5, height = 6.9, units = "in", res = 600)
ggarrange(Panel_A_MTC, Panel_B, ncol=2)
dev.off()

png("../Figure2_v2.png", width = 11.5, height = 6.9, units = "in", res = 600)
ggarrange(Panel_A, Panel_B, ncol=2)
dev.off()

###Panel A MTC only
png("../Figure2_v3.png", width = 8, height = 6.9, units = "in", res = 600)
Panel_A_MTC + labs(tag="")
dev.off()

png("../Figure2_v4.png", width = 8, height = 6.9, units = "in", res = 600)
Panel_A + labs(tag="")
dev.off()

png("../Figure2_v5.png", width = 6, height = 6.9, units = "in", res = 600)
Panel_B + labs(tag="")
dev.off()

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

df_fell_b$all2<-"Felling"
df_fell_b[!is.na(df_fell_b$`20yr SUM`),] %>%
  dplyr::group_by(all2) %>%
  dplyr::summarise(
    min = min(`100yr SUM`, na.rm=T),
    max = max(`100yr SUM`, na.rm=T),
    mean = mean(`100yr SUM`, na.rm=T),
    sd = sd(`100yr SUM`, na.rm=T),
    range = diff(range(`100yr SUM`, na.rm=T)),
    count = sum(!is.na(`100yr SUM`))
  ) %>% ungroup()

tmp<-df_fell_b[,c(1,2)] %>%
  dplyr::mutate(time_bin = cut(as.numeric(Year),
                               breaks = seq(1249,2000, by = 25),
                               right = FALSE,
                               include.lowest = TRUE))
tmp<-tmp %>%
  dplyr::group_by(time_bin) %>%
  dplyr::summarise(counts = sum(all, na.rm=T)) %>%
  ungroup()

tmp %>%
  dplyr::summarise(
    min = min(counts, na.rm=T),
    max = max(counts, na.rm=T),
    mean = mean(counts, na.rm=T),
    sd = sd(counts, na.rm=T),
    range = diff(range(counts, na.rm=T)),
    count = sum(!is.na(counts))
  ) %>% ungroup()


###Melt all environmental variables together (for cross-correlations)
final2<-final
final2$TimePeriod<-as.character(final2$TimePeriod)
final2<-melt(final2)
final2$Region<-final2$variable
names(final2)<-c("Year","Proxy","Value","Region")
final2$Year<-as.numeric(final2$Year)

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

envs<-bind_rows(df_pages[,c(1,4)])
envs$Proxy<-"Eur. Temp."
envs$Region<-"Global Temperature Model (Pages2K Con. 2017)"
names(envs)<-c("Year","Value","Proxy","Region")
envs_all<-bind_rows(envs_all,envs)

envs<-bind_rows(df_pages[,c(1,13)])
envs$Proxy<-"Arctic Temp."
envs$Region<-"Global Temperature Model (Pages2K Con. 2017)"
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

###Cross correlations between environmental variables (all years)
envs_all_corrs<-envs_all[,c(1,2,4)] %>% filter(
  !Region == "*MTC (NSS-DB)",
  !is.na(Value),
  Value > 99)
envs_all_new<-df_CS %>%
  dplyr::select(1,4) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(df_SST[,c(1,4)], by = "Year")
envs_all_new<-df_Forams %>%
  dplyr::select(1,3) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new, by = "Year")
names(envs_all_new)<-c("Year","Keigwin, L.D. 1996","Fohlmeister et al. 2012",
                       "Berner et al. 2011")
envs_all_new<-df_SST_Af %>%
  dplyr::select(1,5) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new, by = "Year")
names(envs_all_new)<-c("Year","McGregor et al. 2007","Keigwin, L.D. 1996","Fohlmeister et al. 2012",
                       "Berner et al. 2011")
envs_all_new<-df_SST_SP %>%
  dplyr::select(1,4) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new, by = "Year")
names(envs_all_new)<-c("Year","Miettinen et al. (2012)","McGregor et al. 2007","Keigwin, L.D. 1996","Fohlmeister et al. 2012",
                       "Berner et al. 2011")
envs_all_new<-df_SST_SP2 %>%
  dplyr::select(1,4) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new, by = "Year")
names(envs_all_new)<-c("Year","Sicre et al. (2011)","Miettinen et al. (2012)","McGregor et al. 2007","Keigwin, L.D. 1996","Fohlmeister et al. 2012",
                       "Berner et al. 2011")
envs_all_new<-df_XRAY %>%
  dplyr::select(1,4) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new, by = "Year")
names(envs_all_new)<-c("Year","Precip. Tiljander et al. (2003)",
                       "Sicre et al. (2011)",
                       "Miettinen et al. (2012)",
                       "McGregor et al. 2007",
                       "Keigwin, L.D. 1996",
                       "Fohlmeister et al. 2012",
                       "Berner et al. 2011")

envs_all_new<-df_pages %>%
  dplyr::select(Year, Europe2K_100yr_av) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new, by = "Year")
names(envs_all_new)<-c("Year","Precip. Tiljander et al. (2003)",
                       "Sicre et al. (2011)",
                       "Miettinen et al. (2012)",
                       "McGregor et al. 2007",
                       "Keigwin, L.D. 1996",
                       "Fohlmeister et al. 2012",
                       "Berner et al. 2011","Pages2K Europe 2017")

envs_all_new<-df_pages %>%
  dplyr::select(Year, Arctic2K_100yr_av) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new, by = "Year")
names(envs_all_new)<-c("Year","Precip. Tiljander et al. (2003)",
                       "Sicre et al. (2011)",
                       "Miettinen et al. (2012)",
                       "McGregor et al. 2007",
                       "Keigwin, L.D. 1996",
                       "Fohlmeister et al. 2012",
                       "Berner et al. 2011","Pages2K Europe 2017","Pages2K Arctic 2017")
envs_all_new$Year<-round(envs_all_new$Year,0)
envs_all_new<-envs_all_new %>% filter(
  complete.cases(.)) %>% distinct()

cor_matrix <- cor(envs_all_new[,c(2:10)])
corrplot(cor_matrix, method = "number", type = "lower",
         tl.col = "black", tl.cex = 0.75,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))

##Cross correlations with MTC (and felling dates)
envs_all_new2<-final2 %>%
  dplyr::select(1,2) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new, by = "Year")
#names(envs_all_new2)[2]<-"MTC NSS_DB (This Study)"
envs_all_new2<-envs_all_new2 %>% filter(
  complete.cases(.)) %>% distinct()
tmp<-envs_all_new2[,c(2:11)]
num_vars <- ncol(tmp)
par(mfrow = c(4, 5))  # Set plot layout
for (i in 1:(num_vars - 1)) {
  for (j in (i + 1):num_vars) {
    ccf(tmp[, i], tmp[, j], lag.max = 3,
        main = paste("CCF:", colnames(tmp)[i], "vs", colnames(tmp)[j]),
        na.action = na.omit)
  }
}
par(mfrow = c(1, 1))

###And felling (removed for now)
envs_all_new2<-df_fell_b %>%
  dplyr::select(1,5) %>% 
  dplyr::mutate(Year = as.numeric(Year) %>% round(0)) %>%
  left_join(envs_all_new2, by = "Year")
names(envs_all_new2)[2]<-"Pop Size. Proxy Ljungqvist et al. 2022"
envs_all_new2<-envs_all_new2 %>% filter(
  complete.cases(.)) %>% distinct()
tmp<-envs_all_new2[,c(2:11)]
names(tmp)[1]<-"NSS DB MTC"

# Get all variable pairs
var_pairs <- combn(names(tmp), 2, simplify = FALSE)
get_ccf(tmp$`Precip. Tiljander et al. (2003)`, tmp$`Sicre et al. (2011)`)
x<-ccf(tmp$`Precip. Tiljander et al. (2003)`, tmp$`Sicre et al. (2011)`)

# Calculate CCF for all pairs
ccf_results <- purrr::map(var_pairs, ~ {
  ccf_data <- get_ccf(tmp[[.x[1]]], tmp[[.x[2]]])
  ccf_data$Var1 <- .x[1]
  ccf_data$Var2 <- .x[2]
  return(ccf_data)
})

# Combine into one data frame
ccf_df <- bind_rows(ccf_results)
conf_level <- 1.96 / sqrt(16)

ccf_df_lower <- ccf_df %>%
  filter(Var1 > Var2)

ggplot(ccf_df_lower, aes(x = lag, y = ccf)) +
  geom_line() +
  geom_hline(yintercept = c(-conf_level, conf_level), linetype = "dashed", color = "blue") +
  facet_grid(Var1 ~ Var2) +  # Retains original shape, but only lower diagonal
  theme_minimal() +
  xlim(-5, 5) +
  ylim(-0.7, 0.7) +
  labs(title = "Cross-Correlations (Lower Triangle)", x = "Lag", y = "CCF")

ggplot(ccf_df, aes(x = lag, y = ccf)) +
  geom_line() +  # Plot CCF line
  geom_hline(yintercept = c(-conf_level, conf_level), 
             linetype = "dashed", color = "blue") +  # Confidence interval lines
  facet_wrap(~ interaction(Var1, Var2), scales = "free") +  # Separate plot for each pair
  theme_minimal() +
  labs(title = "Cross-Correlation Functions", x = "Lag", y = "CCF") +
  theme(strip.text = element_text(size = 10))  

ggplot(ccf_df, aes(x = lag, y = ccf)) +
  geom_line() +
  geom_hline(yintercept = c(-conf_level, conf_level), linetype = "dashed", color = "blue") +
  facet_wrap(~ interaction(Var1, Var2), scales = "free") +  # Avoid unnecessary duplicates
  theme_minimal() +
  xlim(-5, 5) +
  labs(title = "Cross-Correlations of Env. Variables", x = "Lag", y = "CCF")

cor_matrix <- cor(tmp)
par(mfrow=c(1,1))
corrplot(cor_matrix, method = "number", type = "lower",
         tl.col = "black", tl.cex = 0.75,
         col = colorRampPalette(c("navyblue","grey80","darkred"))(200))

p_values <- cor.mtest(cor_matrix)
png("../SupplementaryFiles/SuppFigureCorr_MTC.png", width = 6, height = 6, units = "in", res = 600)
corrplot(cor_matrix, method = "number",
         main = "",
         p.mat = p_values,
         sig.level = 0.05,
         insig = "blank",
         type = "lower",
         tl.col = "black", tl.cex = 0.85,
         col = colorRampPalette(c("navyblue","white","darkred"))(200))
dev.off()
###CCFs for first millennia


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
                     