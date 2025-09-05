###Customised Functions for NSS-DB from Buss et al. 2025 ####
#Author: Dr Danny L Buss
#Date: 4th Sept 2025

###Function to generate random time data
round_to_interval <- function(x, interval) {
  round(x / interval) * interval
}

# Function for counting taxa per database
get_counts <- function(df, col) {
  df %>%
    group_by(Database) %>%
    summarise(UniqueCount = n_distinct(.data[[col]]), .groups = "drop") %>%
    mutate(Level = col, Scope = "By Database")
}

get_overall <- function(df, col) {
  df %>%
    summarise(UniqueCount = n_distinct(.data[[col]])) %>%
    mutate(Level = col, Database = "All", Scope = "Overall")
}

###Function to merge split assemblages
merge_assemblages_purrr <- function(df, group_cols, start_col, end_col,
                                    max_span = 200, inclusive_touch = TRUE) {
  start_sym <- ensym(start_col)
  end_sym   <- ensym(end_col)
  start_nm  <- as_string(start_sym)
  end_nm    <- as_string(end_sym)
  
  merge_group <- function(g) {
    g <- arrange(g, !!start_sym, !!end_sym)
    
    # Reduce: walk through each row, carrying a list of merged blocks
    blocks <- reduce(seq_len(nrow(g)), .init = list(), .f = function(acc, i) {
      row <- g[i, , drop = FALSE]
      if (length(acc) == 0) {
        return(list(row))
      }
      
      last <- acc[[length(acc)]]
      cur_start <- last[[start_nm]]
      cur_end   <- last[[end_nm]]
      ns <- row[[start_nm]]
      ne <- row[[end_nm]]
      
      overlaps <- if (inclusive_touch) ns <= cur_end else ns < cur_end
      if (overlaps) {
        new_start <- min(cur_start, ns, na.rm = TRUE)
        new_end   <- max(cur_end,   ne, na.rm = TRUE)
        if ((new_end - new_start) < max_span) {
          # replace the last block
          last[[start_nm]] <- new_start
          last[[end_nm]]   <- new_end
          acc[[length(acc)]] <- last
        } else {
          acc <- append(acc, list(row))
        }
      } else {
        acc <- append(acc, list(row))
      }
      acc
    })
    
    bind_rows(blocks)
  }
  
  df %>%
    arrange(across(all_of(c(group_cols, start_col, end_col)))) %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(~ merge_group(.x)) %>%
    ungroup()
}


generate_random_time_data<- function(start_year, end_year, time_bin_size, reps, seed){
  # Define the time bins (100-year intervals)
  start_year = round_to_interval(as.numeric(as.character(start_year)), time_bin_size)
  end_year = (round_to_interval(as.numeric(as.character(end_year)), time_bin_size) + 1)
  time_bins <- seq(start_year, end_year, by = time_bin_size)
  # Generate a data frame with random time assignments
  set.seed(seed)  # Set a seed for reproducibility
  year<-0
  ifelse(length(time_bins)< 2,
         data <- data.frame(
           observation_id = seq(1:reps),
           year = rep(time_bins,reps)),
         data <- data.frame(
           observation_id = seq(1:reps),
           year = round_to_interval(sample(time_bins, size = reps, replace = TRUE), time_bin_size)))
#  data$year = round_to_interval(year, time_bin_size)
  return(data$year)
}

generate_random_time_data(125, 175, 100, 100, 2)

####Example
#generate_random_time_data(2014, 2400, 50, 3, 2)

###Function to calculate aMTC
calculate_aMTC<-function(spec_count_new){
  aMTC = rep(0, length(unique(spec_count_new$time)))
  tmp_times<-as.character(unique(spec_count_new$time.mid.round))
  df_out = data.frame(
    TimePeriod = numeric(),
    aMTC = numeric()
  )
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
  df_out_tmp = data.frame(
    as.numeric(tmp_times),
    as.numeric(aMTC)
  )
  names(df_out_tmp)<-c("TimePeriod","aMTC")
  df_out<-bind_rows(df_out, df_out_tmp)
  rm(aMTC, MTC, MTC_tmp)
  return(df_out)
}


###Function to calculate alpha diversity for each time bin
calculate_AlphaDiversity<-function(spec_count_new){
  ADiv = rep(0, length(unique(spec_count_new$time)))
  tmp_times<-as.character(unique(spec_count_new$time.mid.round))
  for(i in 1:length(unique(spec_count_new$time))){
    df_tmp<-spec_count_new[spec_count_new$time%in%tmp_times[i],]
    Div = 0
    Div_tmp<-length(unique(df_tmp$GBIF_species))
    Div<-Div + Div_tmp
  }
  ADiv[i]<-Div
  df_out = data.frame(
    as.numeric(tmp_times),
    ADiv)
  names(df_out)<-c("TimePeriod","AlphaDiversity")
  return(df_out)
}

###Function to calculate beta diversity for each time bin whilst accounting for chronological uncertainty
calculate_beta_diversity<-function(spec_count_new){
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


theme_publish <- function(base_size = 12, base_family = "",
                          line_size = 0.25, ...) {
  half_line <- base_size / 2
  small_rel <- 0.8
  small_size <- small_rel * base_size
  
  theme_bw(base_size = base_size, base_family = base_family, ...) %+replace%
    theme(
      rect = element_rect(fill = "transparent", colour = NA, color = NA,
                          linewidth = 0, linetype = 0),
      text = element_text(family = base_family, face = "plain",
                          colour = "black", size = base_size, hjust = 0.5,
                          vjust = 0.5, angle = 0, lineheight = 0.9,
                          margin = ggplot2::margin(), debug = F),
      
      axis.text = element_text(size = small_size),
      axis.text.x = element_text(margin = ggplot2::margin(t = small_size/4),
                                 vjust = 1),
      axis.text.y = element_text(margin = ggplot2::margin(r = small_size/4), 
                                 hjust = 1),
      axis.title.x = element_text(margin = ggplot2::margin(t = small_size,
                                                           b = small_size)),
      axis.title.y = element_text(angle = 90,
                                  margin = ggplot2::margin(r = small_size,
                                                           l = small_size/4)),
      axis.ticks = element_line(colour = "black", linewidth = line_size),
      axis.ticks.length = unit(0.25, 'lines'),
      
      axis.line = element_line(colour = "black", linewidth = line_size),
      axis.line.x = element_line(colour = "black", linewidth = line_size), 
      axis.line.y = element_line(colour = "black", linewidth = line_size), 
      
      legend.spacing = unit(base_size/4, "pt"),
      legend.key = element_blank(),
      legend.key.size = unit(1 * base_size, "pt"),
      legend.key.width = unit(1.5 * base_size, 'pt'),
      legend.text = element_text(size = rel(small_rel)),
      legend.title = element_text(size = rel(small_rel), face = 'bold'),
      legend.position = 'bottom',
      legend.box = 'horizontal',
      
      panel.spacing = unit(1, "lines"),
      panel.background = element_blank(),
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      strip.text = element_text(size = base_size),
      strip.background = element_rect(fill = NA, colour = "black", linewidth = 0.125),
      strip.text.x = element_text(face = 'bold', hjust = 0,
                                  margin = ggplot2::margin(b = small_size/2,
                                                           t = small_size/4)),
      strip.text.y = element_text(angle = -90, face = 'bold',
                                  margin = ggplot2::margin(l = small_size/2,
                                                           r = small_size/4)),
      
      plot.margin = unit(c(5,5,0,0), "pt"),
      plot.background = element_blank(),
      plot.title = element_text(face = "bold", size = 1.2 * base_size, 
                                margin = ggplot2::margin(b = half_line),
                                hjust = 0)
    )
}

add_shaded_rectangle <- function(plot, xmin, xmax, ymin = -Inf, ymax = Inf, fill = "lightgrey", alpha = 0.3) {
  plot + 
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
             fill = fill, alpha = alpha)
}

cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

# Function to calculate CCF for two vectors
get_ccf <- function(x, y, lag.max = 10) {
  ccf_result <- ccf(x, y, plot = FALSE, lag.max = 4, na.action = na.omit)
  return(data.frame(lag = ccf_result$lag, ccf = ccf_result$acf))
}

get_ccf_plot <- function(var1, var2, var_names) {
  # Generate the CCF plot with a custom title
  ccf_result <- ccf(var1, var2, plot = TRUE, main = paste(var_names[1], "vs", var_names[2]))
}

get_significant_ccf <- function(var1, var2, var_names, conf_level, lag_range = -3:3) {
  ccf_result <- ccf(var1, var2, plot = FALSE, lag.max = max(abs(lag_range)))
  ccf_values <- ccf_result$acf[,,1]
  lags <- ccf_result$lag
  
  # Filter for significant correlations within the specified lag range
  significant_indices <- which(lags %in% lag_range & abs(ccf_values) > conf_level)  # Adjust threshold as needed
  if (length(significant_indices) > 0) {
    data.frame(
      Var1 = rep(var_names[1], length(significant_indices)),
      Var2 = rep(var_names[2], length(significant_indices)),
      Lag = lags[significant_indices],
      Correlation = ccf_values[significant_indices]
    )
  } else {
    NULL
  }
}

plot_significant_ccf <- function(ccf_data) {
  ggplot(ccf_data, aes(x = Lag, y = Correlation)) +
    geom_bar(stat = "identity") +
    labs(title = paste(ccf_data$Var1[1], "vs", ccf_data$Var2[1]),
         x = "Lag",
         y = "Correlation")
}

MTC_df <- function(df, name) {
  NSS_sp<-df
  repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
  repeats = as.character(gsub(" ","",repeats))
  test = repeats
  names(NSS_sp)
  ##Subset unique assemblages 
  tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
    distinct()
  NSS_sp_tmp<-tmp1
  #NSS_sp<-NSS_sp[,c(1:48)]
  for(i in 1:length(repeats)){
    NSS_sp_tmp[,test[i]] = 0
  }
  #2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
  #(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
  x<-1
  for(i in 1:length(repeats)){
    for(j in 1:nrow(NSS_sp_tmp)){
      NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                           100, 1, floor(runif(1, min=1, max=50000)))
    }
  }
  NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
  NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
  rm(tmp1, NSS_sp_tmp)
  tmp<-NSS_sp_tmp_all %>%
    dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
  tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
  names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
  tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
  tmp_sp$time<-as.character(tmp_sp$time.mid.round)
  lst<-unique(tmp_sp$Subsample)
  tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
  tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
  df_out<-data.frame(rep(NA, length(test)),
                     rep(0.0, length(test)),
                     rep("Remove",length(test)))
  names(df_out)<-c("TimePeriod","MTC","Subsample")
  for(i in 1:length(lst)){
    spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
    df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
    df_temps$Subsample<-paste(i)
    df_out<-bind_rows(df_out,df_temps)
  }
  df_out<-df_out[!df_out$Subsample=="Remove",]
  df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
  df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
  a<-merge(df_new,df_new2, by="TimePeriod")
  a$Run<-paste(as.character(name),"Global")
  names(a)<-c("TimePeriod","MTC","SDMTC","Run")
  a$labs<-"Global"
  a$Type <- as.character(name)
  a$TimePeriod<-as.character(a$TimePeriod)
  final<-bind_rows(final, a)
  #Britain & Ireland aMTC
  tmp<-NSS_sp_tmp_all %>%
    dplyr::filter(Region=="Britian & Ireland") %>%
    dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
  tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
  names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
  #Round here at the same time period bin size that you chose for generate_random_time_Data above
  #tmp_sp$time.mid.round<-plyr::round_any(tmp_sp$time.mid.round, 100)
  tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
  tmp_sp$time<-as.character(tmp_sp$time.mid.round)
  lst<-unique(tmp_sp$Subsample)
  tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
  tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
  df_out<-data.frame(rep(NA, length(test)),
                     rep(0.0, length(test)),
                     rep("Remove",length(test)))
  names(df_out)<-c("TimePeriod","MTC","Subsample")
  for(i in 1:length(lst)){
    spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
    df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
    df_temps$Subsample<-paste(i)
    df_out<-bind_rows(df_out,df_temps)
  }
  df_out<-df_out[!df_out$Subsample=="Remove",]
  df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
  df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
  a<-merge(df_new,df_new2, by="TimePeriod")
  a$Run<-paste(as.character(name),"Britain & Ireland")
  names(a)<-c("TimePeriod","MTC","SDMTC","Run")
  a$labs<-"Britain & Ireland"
  a$Type <- as.character(name)
  a$TimePeriod<-as.character(a$TimePeriod)
  final<-bind_rows(final, a)
  #western Europe
  tmp<-NSS_sp_tmp_all %>%
    dplyr::filter(Region=="Western Europe") %>%
    dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
  tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
  names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
  tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
  tmp_sp$time<-as.character(tmp_sp$time.mid.round)
  lst<-unique(tmp_sp$Subsample)
  tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
  tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
  df_out<-data.frame(rep(NA, length(test)),
                     rep(0.0, length(test)),
                     rep("Remove",length(test)))
  names(df_out)<-c("TimePeriod","MTC","Subsample")
  for(i in 1:length(lst)){
    spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
    df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
    df_temps$Subsample<-paste(i)
    df_out<-bind_rows(df_out,df_temps)
  }
  df_out<-df_out[!df_out$Subsample=="Remove",]
  df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
  df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
  a<-merge(df_new,df_new2, by="TimePeriod")
  a$Run<-paste(as.character(name),"Western Europe")
  names(a)<-c("TimePeriod","MTC","SDMTC","Run")
  a$labs<-"Western Europe"
  a$Type <- as.character(name)
  a$TimePeriod<-as.character(a$TimePeriod)
  final<-bind_rows(final, a)
  ###Scandinavia
  tmp<-NSS_sp_tmp_all %>%
    dplyr::filter(Region=="Scandinavia") %>%
    dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
  tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
  names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
  tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
  tmp_sp$time<-as.character(tmp_sp$time.mid.round)
  lst<-unique(tmp_sp$Subsample)
  tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
  tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
  df_out<-data.frame(rep(NA, length(test)),
                     rep(0.0, length(test)),
                     rep("Remove",length(test)))
  names(df_out)<-c("TimePeriod","MTC","Subsample")
  for(i in 1:length(lst)){
    spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
    df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
    df_temps$Subsample<-paste(i)
    df_out<-bind_rows(df_out,df_temps)
  }
  df_out<-df_out[!df_out$Subsample=="Remove",]
  df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
  df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
  a<-merge(df_new,df_new2, by="TimePeriod")
  a$Run<-paste(as.character(name),"Scandinavia")
  names(a)<-c("TimePeriod","MTC","SDMTC","Run")
  a$labs<-"Scandinavia"
  a$Type <- as.character(name)
  a$TimePeriod<-as.character(a$TimePeriod)
  final<-bind_rows(final, a)
  return(final)
}

MTC_df_global <- function(df, name) {
  NSS_sp<-df
  repeats = trimws(paste(rep("rep",200),as.character(seq(1:200))))
  repeats = as.character(gsub(" ","",repeats))
  test = repeats
  names(NSS_sp)
  ##Subset unique assemblages 
  tmp1<-NSS_sp %>% dplyr::select(DB_Assemblage_ID, Start_date_CE, End_date_CE, time_bins2, time_bins) %>%
    distinct()
  NSS_sp_tmp<-tmp1
  #NSS_sp<-NSS_sp[,c(1:48)]
  for(i in 1:length(repeats)){
    NSS_sp_tmp[,test[i]] = 0
  }
  #2. Populate columns with random subsampled time data. Choose size of time bins here. Below I have used 100. 
  #(note this slo prints the random seed chosen for each simulation so you can copy and reproduce the results if you want to). 
  x<-1
  for(i in 1:length(repeats)){
    for(j in 1:nrow(NSS_sp_tmp)){
      NSS_sp_tmp[,test[i]][j,] = generate_random_time_data(NSS_sp_tmp$Start_date_CE[j], NSS_sp_tmp$End_date_CE[j],
                                                           100, 1, floor(runif(1, min=1, max=50000)))
    }
  }
  NSS_sp_tmp<-NSS_sp_tmp %>% distinct()
  NSS_sp_tmp_all<-merge(NSS_sp, NSS_sp_tmp[,c(1, 6:205)], all.x=T)
  rm(tmp1, NSS_sp_tmp)
  tmp<-NSS_sp_tmp_all %>%
    dplyr::select(GBIF_species, NISP, Temp.mid, !!!repeats)
  tmp_sp<-reshape2::melt(tmp, id.vars=c("GBIF_species","Temp.mid","NISP"))
  names(tmp_sp)<-c("GBIF_species", "Trait_Temp_Mid", "NISP_count", "Subsample", "time.mid.round")
  tmp_sp$NISP_count<-as.numeric(tmp_sp$NISP)
  tmp_sp$time<-as.character(tmp_sp$time.mid.round)
  lst<-unique(tmp_sp$Subsample)
  tmp_sp$Subsample<-as.character(tmp_sp$Subsample)
  tmp_sp$Trait_Temp_Mid<-as.numeric(tmp_sp$Trait_Temp_Mid)
  df_out<-data.frame(rep(NA, length(test)),
                     rep(0.0, length(test)),
                     rep("Remove",length(test)))
  names(df_out)<-c("TimePeriod","MTC","Subsample")
  for(i in 1:length(lst)){
    spec_count_new<-aggregate(NISP_count ~ GBIF_species + time.mid.round + time + Trait_Temp_Mid, data=tmp_sp[tmp_sp$Subsample==lst[i],], FUN=sum)
    df_temps<-calculate_aMTC(spec_count_new[spec_count_new$time.mid.round > 50 & spec_count_new$time.mid.round < 1950,])
    df_temps$Subsample<-paste(i)
    df_out<-bind_rows(df_out,df_temps)
  }
  df_out<-df_out[!df_out$Subsample=="Remove",]
  df_new<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=mean)
  df_new2<-aggregate(aMTC ~ TimePeriod, data=df_out, FUN=sd)
  a<-merge(df_new,df_new2, by="TimePeriod")
  a$Run<-paste(as.character(name),"Global")
  names(a)<-c("TimePeriod","MTC","SDMTC","Run")
  a$labs<-"Global"
  a$Type <- as.character(name)
  a$TimePeriod<-as.character(a$TimePeriod)
  final<-bind_rows(final, a)
  return(final)
}

summarise_df <- function(df, name) {
  
  species_count <- df %>%
    filter(!is.na(GBIF_species)) %>%
    summarize(Species = n_distinct(GBIF_species),
              Families = n_distinct(GBIF_family),
              Class = n_distinct(GBIF_class),
              Order = n_distinct(GBIF_order)) %>%
    mutate(Dataset = name)
  
  assemblage_count <- df %>%
    distinct(DB_Assemblage_ID) %>%
    tally(name = "Assemblages") %>%
    mutate(Dataset = name)
  
  nisp_count <- df %>%
    filter(!is.na(NISP)) %>%
    summarize(NISP = sum(NISP, na.rm = TRUE)) %>%
    mutate(Dataset = name)
  
  # Merge all summary statistics
  summary_table <- full_join(species_count, assemblage_count, by = "Dataset") %>%
    full_join(nisp_count, by = "Dataset")
  
  return(summary_table)
}

species_richness <- function(data, indices) {
  sampled_data <- data[indices, ]  # Resample sites
  return(rowSums(sampled_data > 0))  # Species richness per site
}

shannon_diversity <- function(data, indices) {
  sampled_data <- data[indices, ]  # Resample sites
  proportions <- sampled_data / rowSums(sampled_data)
  diversity <- -rowSums(proportions * log(proportions), na.rm = TRUE) 
  return(diversity)  # Return Shannon diversity
}

get_ellipse <- function(data, level = 0.95, count = 8) {
  sampled_data <- data[sample(nrow(data), count), ]
  if (var(sampled_data$MDS1) == 0 | var(sampled_data$MDS2) == 0) {
    warning("Selected points have zero variance, returning empty dataframe")
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  df_ellipse <- tryCatch(
    {
      dataEllipse(
        x = sampled_data$MDS1, y = sampled_data$MDS2, levels = level,
        plot.points = FALSE, draw = FALSE, center.pch = NULL
      )
    },
    error = function(e) {
      warning("dataEllipse() failed, returning empty dataframe")
      return(matrix(numeric(0), ncol = 2))
    }
  )
  if (nrow(df_ellipse) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  df_ellipse <- rbind(df_ellipse, df_ellipse[1, ])
  return(data.frame(x = df_ellipse[, 1], y = df_ellipse[, 2]))
}

calculate_overlap <- function(ellipse1, ellipse2) {
  poly1 <- st_polygon(list(as.matrix(ellipse1)))
  poly2 <- st_polygon(list(as.matrix(ellipse2)))
  area1 <- st_area(st_sfc(poly1))
  area2 <- st_area(st_sfc(poly2))
  intersection <- st_area(st_intersection(st_sfc(poly1), st_sfc(poly2)))
  return(as.numeric(intersection / ((area1 + area2) / 2) * 100))  
}

run_pairwise_overlap <- function(data, iterations = 100) {
  groups <- as.character(unique(data$group))
  pairs <- as.data.frame(t(combn(groups, 2)))
  colnames(pairs) <- c("group1", "group2")
  results_list <- list()  
  for (i in 1:nrow(pairs)) {
    g1 <- pairs$group1[i]
    g2 <- pairs$group2[i]
    group1_data <- filter(data, group == g1)
    group2_data <- filter(data, group == g2)
    overlap_values <- c()
    for (j in 1:100) {
      ellipse1 <- get_ellipse(group1_data)
      ellipse2 <- get_ellipse(group2_data)
      overlap <- calculate_overlap(ellipse1, ellipse2)
      overlap_values <- c(overlap_values, overlap)
    }
    results_list[[i]] <- data.frame(
      group1 = g1, group2 = g2,
      mean_overlap = mean(overlap_values),
      sd_overlap = sd(overlap_values),
      overlaps = paste(overlap_values, collapse = ", ")
    )
  }
  return(do.call(rbind, results_list))
}

perform_f_test <- function(raw_data) {
    f_test_result <- aov(overlap ~ interaction(group1, group2), data = raw_data)
    summary(f_test_result)
}

to_chr <- function(x) {
  if (is.numeric(x)) {
    y <- format(x, trim = TRUE, scientific = FALSE, digits = 15)
    y[is.na(x)] <- NA_character_
    y
  } else {
    as.character(x)
  }
}

collapse_col <- function(x, sep = ";") {
  x <- to_chr(x)
  x <- x[!is.na(x) & x != ""]
  u <- unique(x)
  if (length(u) == 0) NA_character_ else paste(u, collapse = sep)
}