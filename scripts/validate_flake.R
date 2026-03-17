library(tidyverse)
library(ggpubr)
source('R/flake_functions.R')

flake_results <- list.files('data/flake/', pattern = '*.rslt', recursive = T, full.names = T)

# plotting function
all_flake <- flake_results |> 
  purrr::map(~plot_flake(.x,
                         var = 'Ts',
                         ylim = c(0,20),
                         xlim = NULL))


ggarrange(plotlist = all_flake[1:(length(flake_results)/2)], 
          nrow = 3, ncol = 6)
ggarrange(plotlist = all_flake[(1 + (length(flake_results)/2)):length(flake_results)], 
          nrow = 3, ncol = 6)

# Read in data
val_data <- read_csv('data/Validation/temperature and oxygen.csv') |> 
  filter(variable == 'TEMP')


lake_IDs <- toupper(unique(str_extract(basename(flake_results), "(?<=).*?(?=_)")))

for (i in 1:length(lake_IDs)) {
  
  lake_ID_use <- lake_IDs[i]  
  flake_IDs <- flake_results[str_detect(toupper(flake_results), lake_ID_use)]
  if (lake_ID_use == 'WIND') {
    lake_ID_use <- 'SBAS'
  }

  # how to know the max - siple, take the max most common
  use_max_depth <- val_data |> 
    filter(site == lake_ID_use) |> 
    group_by(depth) |> 
    summarise(n = n()) |> 
    filter(n > 20) |> 
    slice_max(depth) |> pull(depth)
  
  val_surface <- val_data |> 
    filter(site == lake_ID_use ,
           depth < 1) |> 
    mutate(date = lubridate::parse_date_time2(date, 'dby', cutoff_2000 = 25),
           doy = yday(date),
           year = year(date))
  
  ggarrange(read_flake(flake_IDs[1]) |> 
              slice_head(n = 365) |> 
              ggplot(aes(x=time, y = Ts)) + 
              geom_point(data = val_surface, aes(x=doy, y = value), alpha = 0.1) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = basename(flake_IDs[1])) ,
            read_flake(flake_IDs[2]) |> 
              slice_head(n = 365) |> 
              ggplot(aes(x=time, y = Ts))  + 
              geom_point(data = val_surface, aes(x=doy, y = value), alpha = 0.1) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = basename(flake_IDs[2])),
            nrow = 2) |> 
    ggsave(filename = file.path('./data/flake/plots', paste0(lake_IDs[i], '_val_surface.png')),
           width = 7, height = 12, units = 'cm')
  
  val_bottom <- val_data |> 
    filter(site == lake_ID_use ,
           depth >= use_max_depth) |> 
    reframe(.by = date, value = mean(value, na.rm = T)) |> 
    mutate(date = lubridate::parse_date_time2(date, 'dby', cutoff_2000 = 25),
           doy = yday(date),
           year = year(date))
  
  ggarrange(read_flake(flake_IDs[1]) |> 
              slice_head(n = 365) |> 
              ggplot(aes(x=time, y = Tb)) + 
              geom_point(data = val_bottom, aes(x=doy, y = value), alpha = 0.2) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = basename(flake_IDs[1])) ,
            read_flake(flake_IDs[2]) |> 
              slice_head(n = 365) |> 
              ggplot(aes(x=time, y = Tb))  + 
              geom_point(data = val_bottom, aes(x=doy, y = value), alpha = 0.2) +
              geom_line() +
              coord_cartesian(ylim = c(0,30)) +
              theme_bw() +
              labs(title = basename(flake_IDs[2])),
            nrow = 2) |> 
    ggsave(filename = file.path('./data/flake/plots', paste0(lake_IDs[i], '_val_bottom.png')),
           width = 7, height = 12, units = 'cm')
  
  
  ggarrange(read_flake(flake_IDs[1]) |> 
              slice_tail(n = 365) |> 
              ggplot(aes(x=time, y = h_ML)) +
              geom_line() +
              theme_bw()+
              labs(title = basename(flake_IDs[1])),
            read_flake(flake_IDs[2]) |> 
              slice_tail(n = 365) |> 
              ggplot(aes(x=time, y = h_ML)) +
              geom_line() +
              theme_bw() +
              labs(title = basename(flake_IDs[2]))) |> 
    ggsave(filename = file.path('./data/flake/plots', paste0(lake_IDs[i], '_ML.png')),
           width = 12, height = 7, units = 'cm')
}
