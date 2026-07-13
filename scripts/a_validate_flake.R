library(tidyverse)
library(ggpubr)

source('R/flake_functions.R')

flake_dir <- 'data/flake'

flake_results <- list.files(flake_dir, pattern = '*.rslt', recursive = T, full.names = T)
flake_nlms <- list.files(flake_dir, pattern = '*.nml', recursive = T, full.names = T)

# plotting function
all_flake <- flake_results |> 
  purrr::map(~plot_flake(.x,
                         var = 'Ts',
                         ylim = c(0,20),
                         xlim = NULL))
names(all_flake) <- str_remove(basename(flake_results), '.rslt')

ggarrange(plotlist = all_flake[1:(length(flake_results)/2)], 
          nrow = 3, ncol = 6)
ggarrange(plotlist = all_flake[(1 + (length(flake_results)/2)):length(flake_results)], 
          nrow = 3, ncol = 6)

# Read in data
val_data <- read_csv('data/Validation/temperature and oxygen.csv') |> 
  filter(variable == 'TEMP')

WIND <- c(47007, 47008) # these are the WBIDs for the NBAS and SBAS
WIND_WBID <- 29233

match_IDs <- readxl::read_xlsx('data/SITE ID_MULTIPLE DATA SOURCES_LD LAKES.xlsx') 

lake_IDs <- unique(str_extract(basename(flake_results), "(?<=).*?(?=_)"))

# Plots for validation -------------------------
for (i in 1:length(lake_IDs)) {
  
  lake_ID_use <- lake_IDs[i]  
  flake_IDs <- flake_results[str_detect(toupper(flake_results), lake_ID_use)]
  flake_nml <-  glmtools::read_nml(flake_nlms[str_detect(toupper(flake_nlms), lake_ID_use)])
  
  if (lake_ID_use == WIND_WBID) {
    lake_ID_use <- 47008 # change to the ID for the SBAS for validation data
  }
  
  lake_name_use <- match_IDs$LAKE_Lakes_Tour_Chem_TeOx.xlsx[which(match_IDs$`WBID_Lake District_UKCEH Portal data_raw.xlsx` == lake_ID_use)]
  
  # read the FLake results
  clear_df <- read_flake(flake_IDs[1]) |> 
    mutate(doy = yday(as_date(time)),
           year = year(as_date(time)))|> 
    # slice_tail(n = 365) |> 
    mutate(time = row_number())
  
  turbid_df <- read_flake(flake_IDs[2]) |> 
    mutate(doy = yday(as_date(time)),
           year = year(as_date(time)))|> 
    # slice_tail(n = 365) |> 
    mutate(time = row_number())
  
  mean_df <- as.data.frame(Map(function(x, y) {(x + y) / 2}, clear_df, turbid_df))
  
  # how to know the max - siple, take the max most common
  use_max_depth <- val_data |> 
    filter(site == lake_name_use) |> 
    group_by(depth) |> 
    summarise(n = n()) |> 
    filter(n > 20) |> 
    slice_max(depth) |> pull(depth)
  
  val_surface <- val_data |> 
    filter(site == lake_name_use ,
           depth < 1) |> 
    mutate(date = lubridate::parse_date_time2(date, 'dby', cutoff_2000 = 25),
           doy = yday(date),
           year = year(date))
  
  ggarrange(clear_df |> 
              ggplot(aes(x=doy, y = Ts, group = year)) + 
              geom_point(data = val_surface, aes(x=doy, y = value), alpha = 0.1) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = basename(flake_IDs[1])) ,
            turbid_df |> 
              ggplot(aes(x=doy, y = Ts, group = year))  + 
              geom_point(data = val_surface, aes(x=doy, y = value), alpha = 0.1) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = basename(flake_IDs[2])),
            # calculate a mean from the clear and turbid
            mean_df  |> 
              ggplot(aes(x=doy, y = Ts, group = year))  + 
              geom_point(data = val_surface, aes(x=doy, y = value), alpha = 0.1) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = 'mean'),
            
            nrow = 3) |> 
    ggsave(filename = file.path('./data/flake/plots', paste0(lake_IDs[i], '_val_surface_obs.png')),
           width = 7, height = 18, units = 'cm')
  
  val_bottom <- val_data |> 
    filter(site == lake_name_use ,
           depth >= use_max_depth) |> 
    reframe(.by = date, value = mean(value, na.rm = T)) |> 
    mutate(date = lubridate::parse_date_time2(date, 'dby', cutoff_2000 = 25),
           doy = yday(date),
           year = year(date))
  
  ggarrange(clear_df |> 
              ggplot(aes(x=doy, y = Tb, group = year)) + 
              geom_point(data = val_bottom, aes(x=doy, y = value), alpha = 0.1) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = basename(flake_IDs[1])) ,
            turbid_df |> 
              ggplot(aes(x=doy, y = Tb, group = year))  + 
              geom_point(data = val_bottom, aes(x=doy, y = value), alpha = 0.1) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = basename(flake_IDs[2])),
            # calculate a mean from the clear and turbid
            mean_df  |>
              ggplot(aes(x=doy, y = Tb, group = year))  +
              geom_point(data = val_bottom, aes(x=doy, y = value), alpha = 0.1) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = 'mean'),
            
            nrow = 3) |> 
    ggsave(filename = file.path(flake_dir, 'plots', paste0(lake_IDs[i], '_val_bottom_obs.png')),
           width = 7, height = 18, units = 'cm')

  
}

