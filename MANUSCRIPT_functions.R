###Function to generate random time data

generate_random_time_data<- function(start_year, end_year, time_bin_size, reps, seed){
  # Define the time bins (100-year intervals)
  time_bins <- seq(start_year, end_year, by = time_bin_size)
  # Generate a data frame with random time assignments
  set.seed(seed)  # Set a seed for reproducibility
  ifelse(length(time_bins)< 2,
         data <- data.frame(
           observation_id = seq(1:reps),
           year = rep(time_bins,reps)),
         data <- data.frame(
           observation_id = seq(1:reps),
           year = sample(time_bins, size = reps, replace = TRUE)))
  return(data$year)
}

###Function to calculate aMTC
calculate_aMTC<-function(spec_count_new){
  aMTC = rep(0, length(unique(spec_count_new$time)))
  tmp_times<-as.character(unique(spec_count_new$time.mid.round))
  for(i in 1:length(unique(spec_count_new$time))){
    df_tmp<-spec_count_new[spec_count_new$time%in%tmp_times[i],]
    df_tmp<-df_tmp[complete.cases(df_tmp$time.mid.round),]
    MTC = 0 
    MTC_tmp = 0
    for(j in 1:length(unique(df_tmp$GBIF_species))){
      spec<-df_tmp$GBIF_species
      MTC_tmp<-df_tmp[df_tmp$GBIF_species%in%spec[j],]$Trait_Temp_Mid*df_tmp[df_tmp$GBIF_species%in%spec[j],]$NISP_count
      MTC<-MTC + MTC_tmp
    }
    aMTC[i]<-MTC / sum(df_tmp$NISP_count)
  }
  df_out = data.frame(
    as.numeric(tmp_times),
    aMTC
  )
  names(df_out)<-c("TimePeriod","aMTC")
  return(df_out)
}