###Random time assignment for NSS_sp DF 

###RUN1###
#1. 1000 subsamples - 100-yr time bins (All Data)
#Create blank columns in dataframe where simulated timestamps will be generated
repeats = trimws(paste(rep("rep",1000),as.character(seq(1:1000))))
repeats = as.character(gsub(" ","",repeats))
test = repeats
test = repeats
#repeats = repeats[1:10]
NSS_sp<-NSS_sp[,c(1:45)]
for(i in 1:length(repeats)){
  NSS_sp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp)){
    NSS_sp[,test[i]][j,] = generate_random_time_data_yr(NSS_sp$Start.date.CE[j], NSS_sp$End.date.CE[j],
                                                        100, 1, floor(runif(1, min=1, max=50000)))
  }
}

#3. Subset columns required to run aMTC function and melt dataframe
tmp_sp<-melt(NSS_sp[,c(14, 29, 46:245,42)], id.vars=c("GBIF_species","Trait_Temp_Mid","NISP_count"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
tmp_sp$time.mid.round<-round_any(tmp_sp$time.mid.round, 100)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)

#4. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","aMTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[!spec_count_new$time=="0" & !spec_count_new$time=="2000",])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"100-yr time bins (All Data, mean + SD)"
names(a)<-c("TimePeriod","aMTC","SDMTC")
#ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw() + ggtitle("100-yr time bins (All)")
ggplot(a, aes(TimePeriod, aMTC)) + geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins (mean + SD)")
final<-a

###RUN2###
#1000 subsamples - 200-yr time bins (All Data)
#1. Create blank columns in dataframe where simulated timestamps will be generated
repeats = trimws(paste(rep("rep",1000),as.character(seq(1:1000))))
repeats = as.character(gsub(" ","",repeats))
test = repeats
#repeats = repeats[1:10]
NSS_sp<-NSS_sp[,c(1:45)]
for(i in 1:length(repeats)){
  NSS_sp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp)){
    NSS_sp[,test[i]][j,] = generate_random_time_data_yr(NSS_sp$Start.date.CE[j], NSS_sp$End.date.CE[j],
                                                        200, 1, floor(runif(1, min=1, max=50000)))
  }
}

#3. Subset columns required to run aMTC function and melt dataframe
tmp_sp<-melt(NSS_sp[,c(14, 29, 46:245,42)], id.vars=c("GBIF_species","Trait_Temp_Mid","NISP_count"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
tmp_sp$time.mid.round<-round_any(tmp_sp$time.mid.round, 200)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)

#4. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","aMTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[!spec_count_new$time=="0" & !spec_count_new$time=="2000",])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"200-yr time bins (All Data, mean + SD)"
names(a)<-c("TimePeriod","aMTC","SDMTC")
#ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw() + ggtitle("200-yr time bins (All)")
ggplot(a, aes(TimePeriod, aMTC)) + geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("200-yr time bins (mean + SD)")
final<-bind_rows(final,a)

###RUN3###
#1000 subsamples - 300-yr time bins (All Data)
#1. Create blank columns in dataframe where simulated timestamps will be generated
repeats = trimws(paste(rep("rep",1000),as.character(seq(1:1000))))
repeats = as.character(gsub(" ","",repeats))
test = repeats
#repeats = repeats[1:10]
NSS_sp<-NSS_sp[,c(1:45)]
for(i in 1:length(repeats)){
  NSS_sp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp)){
    NSS_sp[,test[i]][j,] = generate_random_time_data_yr(NSS_sp$Start.date.CE[j], NSS_sp$End.date.CE[j],
                                                        300, 1, floor(runif(1, min=1, max=50000)))
  }
}

#3. Subset columns required to run aMTC function and melt dataframe
tmp_sp<-melt(NSS_sp[,c(14, 29, 46:245,42)], id.vars=c("GBIF_species","Trait_Temp_Mid","NISP_count"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
tmp_sp$time.mid.round<-round_any(tmp_sp$time.mid.round, 300)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)

#4. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","aMTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[!spec_count_new$time=="0" & !spec_count_new$time=="2000",])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"300-yr time bins (All Data, mean + SD)"
names(a)<-c("TimePeriod","aMTC","SDMTC")
#ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw() + ggtitle("300-yr time bins (All)")
ggplot(a, aes(TimePeriod, aMTC)) + geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("300-yr time bins (mean + SD)")
final<-bind_rows(final,a)

###RUN4###
#500 subsamples - 50-yr time bins (All Data)
###Run 4. Random time assignment for NSS_sp DF - 500 subsamples - 50-yr time bins (All data)
#1. Create blank columns in dataframe where simulated timestamps will be generated
repeats = trimws(paste(rep("rep",500),as.character(seq(1:500))))
repeats = as.character(gsub(" ","",repeats))
test = repeats
#repeats = repeats[1:10]
NSS_sp<-NSS_sp[,c(1:45)]
for(i in 1:length(repeats)){
  NSS_sp[,test[i]] = 0
}

#2. Populate columns with random subsampled time data 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp)){
    NSS_sp[,test[i]][j,] = generate_random_time_data_yr(NSS_sp$Start.date.CE[j], NSS_sp$End.date.CE[j],
                                                        50, 1, floor(runif(1, min=1, max=50000)))
  }
}

#3. Subset columns required to run aMTC function and melt dataframe
tmp_sp<-melt(NSS_sp[,c(14, 29, 46:245,42)], id.vars=c("GBIF_species","Trait_Temp_Mid","NISP_count"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
tmp_sp$time.mid.round<-round_any(tmp_sp$time.mid.round, 50)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)

#4. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","aMTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[!spec_count_new$time=="0" & !spec_count_new$time=="2000",])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
a<-merge(df_new,df_new2, by="TimePeriod")
a$run<-"50-yr time bins (All Data, mean + SD)"
names(a)<-c("TimePeriod","aMTC","SDMTC")
#ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw() + ggtitle("50-yr time bins (All)")
ggplot(a, aes(TimePeriod, aMTC)) + geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("50-yr time bins (mean + SD)")
final<-bind_rows(final,a)
names(final)[4]<-"run"
ggplot(final, aes(TimePeriod, aMTC)) + geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + facet_wrap(~run, scale="free_x")

###Plot all data together
ggplot(final, aes(TimePeriod, aMTC)) + 
  geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + 
  theme_bw() + facet_wrap(~run, ncol=2) + ggtitle("All Data")


###REGIONSPECIFIC DATA###

#1. British Isles Only (100-yr time bins) - 300 simulations
#Create blank columns in dataframe where simulated timestamps will be generated
repeats = trimws(paste(rep("rep",300),as.character(seq(1:300))))
repeats = as.character(gsub(" ","",repeats))
test = repeats
test = repeats
#repeats = repeats[1:10]
NSS_sp<-NSS_sp[,c(1:45)]
BI<-c("England","Scotland","Ireland")
NSS_sp_BI<-NSS_sp[NSS_sp$Country%in%BI,]
for(i in 1:length(repeats)){
  NSS_sp_BI[,test[i]] = 0
}

#2. Populate columns with random subsampled time data 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_BI)){
    NSS_sp_BI[,test[i]][j,] = generate_random_time_data_yr(NSS_sp_BI$Start.date.CE[j], NSS_sp_BI$End.date.CE[j],
                                                        100, 1, floor(runif(1, min=1, max=50000)))
  }
}

#3. Subset columns required to run aMTC function and melt dataframe
tmp_sp<-melt(NSS_sp_BI[,c(14, 29, 46:245,42)], id.vars=c("GBIF_species","Trait_Temp_Mid","NISP_count"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
tmp_sp$time.mid.round<-round_any(tmp_sp$time.mid.round, 100)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)

#4. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","aMTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[!spec_count_new$time=="0" & !spec_count_new$time=="2000",])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
b<-merge(df_new,df_new2, by="TimePeriod")
b$run<-"100-yr time bins (British Isles, mean + SD)"
names(b)<-c("TimePeriod","aMTC","SDMTC")
#ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw() + ggtitle("100-yr time bins (All)")
ggplot(b, aes(TimePeriod, aMTC)) + geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins (British Isles)")
final2<-b

#2. Scandinavia
#Create blank columns in dataframe where simulated timestamps will be generated
repeats = trimws(paste(rep("rep",300),as.character(seq(1:300))))
repeats = as.character(gsub(" ","",repeats))
test = repeats
test = repeats
#repeats = repeats[1:10]
NSS_sp<-NSS_sp[,c(1:45)]
SC<-c("Norway","Sweden","Denmark")
NSS_sp_SC<-NSS_sp[NSS_sp$Country%in%SC,]
for(i in 1:length(repeats)){
  NSS_sp_SC[,test[i]] = 0
}

#2. Populate columns with random subsampled time data 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_SC)){
    NSS_sp_SC[,test[i]][j,] = generate_random_time_data_yr(NSS_sp_SC$Start.date.CE[j], NSS_sp_SC$End.date.CE[j],
                                                           100, 1, floor(runif(1, min=1, max=50000)))
  }
}

#3. Subset columns required to run aMTC function and melt dataframe
tmp_sp<-melt(NSS_sp_SC[,c(14, 29, 46:245,42)], id.vars=c("GBIF_species","Trait_Temp_Mid","NISP_count"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
tmp_sp$time.mid.round<-round_any(tmp_sp$time.mid.round, 100)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)

#4. Create dataframe of aMTCs for each simulation (currently there are 300 but probably need to up to 500 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","aMTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[!spec_count_new$time=="0" & !spec_count_new$time=="2000",])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
b<-merge(df_new,df_new2, by="TimePeriod")
b$run<-"100-yr time bins (Scandinavia, mean + SD)"
names(b)<-c("TimePeriod","aMTC","SDMTC")
#ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw() + ggtitle("100-yr time bins (Scandinavia)")
ggplot(b, aes(TimePeriod, aMTC)) + geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins (Scandinavia)")
final2<-bind_rows(final2,b)

#3. Western Europe
#Create blank columns in dataframe where simulated timestamps will be generated
repeats = trimws(paste(rep("rep",300),as.character(seq(1:300))))
repeats = as.character(gsub(" ","",repeats))
test = repeats
test = repeats
#repeats = repeats[1:10]
NSS_sp<-NSS_sp[,c(1:45)]
MWE<-c("Germany","Netherlands","Belgium")
NSS_sp_MWE<-NSS_sp[NSS_sp$Country%in%MWE,]
for(i in 1:length(repeats)){
  NSS_sp_MWE[,test[i]] = 0
}

#2. Populate columns with random subsampled time data 
#(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
x<-1
for(i in 1:length(repeats)){
  for(j in 1:nrow(NSS_sp_MWE)){
    NSS_sp_MWE[,test[i]][j,] = generate_random_time_data_yr(NSS_sp_MWE$Start.date.CE[j], NSS_sp_MWE$End.date.CE[j],
                                                           100, 1, floor(runif(1, min=1, max=50000)))
  }
}

#3. Subset columns required to run aMTC function and melt dataframe
tmp_sp<-melt(NSS_sp_MWE[,c(14, 29, 46:245,42)], id.vars=c("GBIF_species","Trait_Temp_Mid","NISP_count"))
names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
tmp_sp$time.mid.round<-round_any(tmp_sp$time.mid.round, 100)
tmp_sp$time<-as.character(tmp_sp$time.mid.round)
lst<-unique(tmp_sp$Subsample)

#4. Create dataframe of aMTCs for each simulation (currently there are 200 but probably need to up to 2000 and rerun (takes a few hours))
df_out<-data.frame(rep(NA, length(test)),
                   rep(0.0, length(test)),
                   rep("Remove",length(test)))
names(df_out)<-c("TimePeriod","aMTC","Subsample")
for(i in 1:length(lst)){
  spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
  df_temps<-calculate_aMTC(spec_count_new[!spec_count_new$time=="0" & !spec_count_new$time=="2000",])
  df_temps$Subsample<-paste(i)
  df_out<-bind_rows(df_out,df_temps)
}

###Plot random time data
df_out<-df_out[!df_out$Subsample=="Remove",]
df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
b<-merge(df_new,df_new2, by="TimePeriod")
b$run<-"100-yr time bins (Scandinavia, mean + SD)"
names(b)<-c("TimePeriod","aMTC","SDMTC")
#ggplot(df_out, aes(TimePeriod, aMTC, group=Subsample, color=Subsample)) + geom_line(colour = "blue", size = 1.2, alpha=0.2) + theme_bw() + ggtitle("100-yr time bins (Western Europe)")
ggplot(b, aes(TimePeriod, aMTC)) + geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + theme_bw() + ggtitle("100-yr time bins (Western Europe)")
final2<-bind_rows(final2,b)
names(final2)[4]<-"run"

###Plot all regional data together
ggplot(final2, aes(TimePeriod, aMTC)) + 
  geom_line(aes(y=aMTC), colour = "darkred", size = 1.2, alpha=0.9) + 
  geom_ribbon(aes(y=aMTC, ymin = aMTC - SDMTC, ymax = aMTC + SDMTC), fill = "darkred", size = 1.2, alpha=0.3) + 
  theme_bw() + facet_wrap(~run, ncol=2) + ggtitle("Regional Data (300 simulations")
