#--------------------------------------#
## Project: Lake District pathways modelling
## Script purpose: model validation using "observed" nutrient loads
## Date: 2025-12-03
## Author: Freya Olsson
#--------------------------------------#
library(tidyverse)
library(ggpubr)

out_dir <- 'output/'
val_dir <- 'data/Validation'
start_date <- '2000-01-01'

# read in observations -----------------------
# make sure you have run 00_compile_obs.R
WIND <- c(47007, 47008) # these are the WBIDs for the NBAS and SBAS
WIND_WBID <- 29233

obs <- read_csv('data/Validation/output/LD_combined_database_core_variables.csv') |> 
  filter(date > ymd('1990-01-01')) |> 
  mutate(site_id = ifelse(site_id %in% WIND, WIND_WBID, site_id))


# Read in model output ----------------
pclake_results <- list.files(out_dir, pattern = 'turbid.csv', recursive = T, full.names = T)

subset <- T # run all lakes or not

# the lakes portal data has all the basic info we need
lakes_portal_df <- read_csv('data/lakes4PCLake.csv', show_col_types = F)

if (subset == T) {
  # LakeIDs to loop through
  lake_names_lookup <- distinct(obs, site_id, site_name) |> 
    filter(!is.na(site_id)) |> 
    mutate(site_id = ifelse(site_id %in% WIND, WIND_WBID, site_id),
           site_name = ifelse(site_id == WIND_WBID, "Windermere", site_name)) |> 
    distinct(site_id, site_name)
  # the NBAS and SBAS have seperate WBIDs but the one from the lakes portal has the combined one which is different
  
  pclake_results <- pclake_results[str_detect(string = pclake_results, paste(lake_names_lookup$site_id, collapse = "|"))]
  
} else {
  
  lake_names_lookup <- data.frame(site_id = lakes_portal_df$WBID, 
                                  site_name = lakes_portal_df$NAME)
  
  
}


## Read in the model output -----------------#
baseline_runs <- pclake_results |>  
  lapply(read_csv, show_col_types = F, id = 'filename') |> 
  bind_rows() |> 
  group_by(filename) |>
  slice_tail(n = 365*25) |>  # take the last 25 years
  mutate(time = row_number()) |> # renumber
  ungroup() |> 
  mutate(Date = as_date(time, origin = start_date), 
         site_id = parse_number(filename)) |> 
  left_join(lake_names_lookup,
            by = join_by(site_id))


# plots

(baseline_runs |> 
    select(Date, aSecchiT, site_name, site_id) |> 
    left_join(filter(obs, variable_code == 'SECC'),
              by = join_by(Date == date, site_id, site_name)) |> 
    filter(.by = c(site_name, site_id), any(!is.na(value))) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    ggplot(aes(x=doy, y = value)) +
    geom_point(size = 0.9, alpha = 0.3) +
    geom_line(aes(y = aSecchiT, group = year), alpha = 0.6, colour = 'goldenrod') + # PClake in mg/m3
    theme_bw() +
    facet_wrap(~site_id+site_name, scales = 'free_y')) #|> 
  # ggsave(filename = file.path(out_dir, 'plot', 'secchi_val_turbid.png'),
  #        height = 20, width = 20, units = 'cm')

(baseline_runs |> 
    select(Date, oChlaEpi, site_name, site_id) |> 
    left_join(filter(obs, variable_code == 'TOCA'),
              by = join_by(Date == date, site_id, site_name)) |> 
    filter(.by = c(site_name, site_id), any(!is.na(value))) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    ggplot(aes(x=doy, y = value)) + # obs are in ug/L
    geom_point(size = 0.9, alpha = 0.3) +
    geom_line(aes(y = oChlaEpi, group = year), alpha = 0.6, colour = 'seagreen') + # PClake in mg/m3
    theme_bw() +
    facet_wrap(~site_id+site_name, scales = 'free_y')) #|> 
  # ggsave(filename = file.path(out_dir, 'plot', 'chla_val_turbid.png'),
  #        height = 20, width = 20, units = 'cm')


(baseline_runs |> 
    select(Date, oPTotWEpi, site_name, site_id) |> 
    left_join(filter(obs, variable_code == 'TOTP'),
              by = join_by(Date == date, site_id, site_name)) |> 
    filter(.by = c(site_name, site_id), any(!is.na(value))) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    ggplot(aes(x=doy, y = value/1000)) + # obs are in ug/L
    geom_point(size = 0.9, alpha = 0.3) +
    geom_line(aes(y = oPTotWEpi, group = year), alpha = 0.6, colour = 'orchid4') + # PClake in g/m3
    theme_bw() +
    facet_wrap(~site_id+site_name, scales = 'free_y')) #|> 
  # ggsave(filename = file.path(out_dir, 'plot', 'TP_val_turbid.png'),
  #        height = 20, width = 20, units = 'cm')


(baseline_runs |> 
    select(Date, uQInEpi, site_id, site_name) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    
    ggplot(aes(x=doy, y = uQInEpi)) + # obs are in ug/L
    geom_line(aes(group = year), alpha = 0.6, colour = 'blue4') + # PClake in g/m3
    theme_bw() +
    facet_wrap(~site_id+site_name, scales = 'free_y')) 


# calculate the TN:TP ratio
#TN:TP (molar) = (TN / 14.007) ÷ (TP / 30.974) 
N_mol <- 14.007
P_mol <- 30.974

baseline_runs |> 
  select(Date, oNTotWEpi, oPTotWEpi, site_name, site_id) |> 
  mutate(N = oNTotWEpi / N_mol,
         P = oPTotWEpi / P_mol,
         `N:P` = N/P) |> 
  reframe(.by = c(site_name,site_id),
          mean_ratio = mean(`N:P`)) |> 
  arrange(mean_ratio)

baseline_runs |> 
  select(Date, oNTotWEpi, oPTotWEpi, site_name, site_id) |> 
  mutate(N = oNTotWEpi / N_mol,
         P = oPTotWEpi / P_mol,
         `N:P` = N/P) |> 
  mutate(lake = fct_reorder(site_name, desc(`N:P`), .fun='median')) |> 
  ggplot(aes(x=lake, y=`N:P`)) + 
  geom_boxplot() + 
  # coord_cartesian(ylim = c(0, 600)) + 
  theme_bw() +
  geom_hline(yintercept = 30, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90))

baseline_runs |> 
  select(Date, uNLoadEpi, uPLoadEpi, site_name, site_id) |> 
  mutate(N = uNLoadEpi / N_mol,
         P = uPLoadEpi / P_mol,
         `N:P` = N/P) |> 
  mutate(lake = fct_reorder(site_name, desc(`N:P`), .fun='median')) |> 
  ggplot(aes(x=lake, y=`N:P`)) + 
  geom_boxplot() + 
  theme_bw() +
  geom_hline(yintercept = 30, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90))
  


