#--------------------------------------#
## Project: Lake District nutrient critical limit estimation
## Script purpose: model validation using "observed" nutrient loads
## Date: 2025-12-03
## Author: Freya Olsson
#--------------------------------------#
library(tidyverse)
library(ggpubr)

out_dir <- 'output/'
val_dir <- 'data/Validation'
start_date <- '2000-01-01'


# Read in model output ----------------
pclake_results <- list.files(out_dir, pattern = 'turbid.csv', recursive = T, full.names = T)

subset <- T # run all lakes or not

# the lakes portal data has all the basic info we need
lakes_portal_df <- read_csv('data/lakes4PCLake.csv', show_col_types = F)

if (subset == T) {
  # LakeIDs to loop through
  lakeIDs <- readxl::read_xlsx('data/SITE ID_MULTIPLE DATA SOURCES_LD LAKES.xlsx') |> 
    filter(!is.na(site_chemistry.csv))   # remove Loughrigg because it's not right
  
  lake_names_lookup <- lakeIDs$`WBID_Lake District_UKCEH Portal data_raw.xlsx`
  names(lake_names_lookup) <- str_extract(lakeIDs$site_chemistry.csv, "....")
  
  WIND <- c(47007, 47008) # these are the WBIDs for the NBAS and SBAS
  WIND_WBID <- 29233
  
  lake_names_lookup[which(lake_names_lookup %in% WIND)] <- WIND_WBID
  names(lake_names_lookup)[which(lake_names_lookup == WIND_WBID)] <- "WIND"
  unique(lake_names_lookup)
  # the NBAS and SBAS have seperate WBIDs but the one from the lakes portal has the combined one which is different
  
  pclake_results <- pclake_results[str_detect(string = pclake_results, paste(lake_names_lookup, collapse = "|"))]
  
} else {
  
  lake_names_lookup <- lakes_portal_df$WBID
  names(lake_names_lookup) <- lakes_portal_df$NAME
  
  
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
         WBID = parse_number(filename)) |> 
  left_join(distinct(data.frame(WBID = lake_names_lookup,
                                NAME = names(lake_names_lookup))),
            by = join_by(WBID))

# read in observations -----------------------
# make sure you have run 00_compile_obs.R
obs <- read_csv('data/Validation/output/LD_combined_database_core_variables.csv') |> 
  filter(date > ymd('1990-01-01')) |> 
  mutate(site_id = ifelse(site_id %in% WIND, WIND_WBID, site_id))

# plots

(baseline_runs |> 
    select(Date, aSecchiT, NAME, WBID) |> 
    left_join(filter(obs, variable_code == 'SECC'),
              by = join_by(Date == date, WBID == site_id)) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    ggplot(aes(x=doy, y = value)) +
    geom_point(size = 0.9, alpha = 0.3) +
    geom_line(aes(y = aSecchiT, group = year), alpha = 0.6, colour = 'goldenrod') + # PClake in mg/m3
    theme_bw() +
    facet_wrap(~WBID+NAME, scales = 'free_y')) #|> 
  # ggsave(filename = file.path(out_dir, 'plot', 'secchi_val_turbid.png'),
  #        height = 20, width = 20, units = 'cm')

(baseline_runs |> 
    select(Date, oChlaEpi, NAME, WBID) |> 
    left_join(filter(obs, variable_code == 'TOCA'), 
              by = join_by(Date == date, WBID == site_id)) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    ggplot(aes(x=doy, y = value)) + # obs are in ug/L
    geom_point(size = 0.9, alpha = 0.3) +
    geom_line(aes(y = oChlaEpi, group = year), alpha = 0.6, colour = 'seagreen') + # PClake in mg/m3
    theme_bw() +
    facet_wrap(~WBID+NAME, scales = 'free_y')) #|> 
  # ggsave(filename = file.path(out_dir, 'plot', 'chla_val_turbid.png'),
  #        height = 20, width = 20, units = 'cm')


(baseline_runs |> 
    select(Date, oPTotWEpi, NAME, WBID) |> 
    left_join(filter(obs, variable_code == 'TOTP'), 
              by = join_by(Date == date, WBID == site_id)) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    
    ggplot(aes(x=doy, y = value/1000)) + # obs are in ug/L
    geom_point(size = 0.9, alpha = 0.3) +
    geom_line(aes(y = oPTotWEpi, group = year), alpha = 0.6, colour = 'orchid4') + # PClake in g/m3
    theme_bw() +
    facet_wrap(~WBID+NAME, scales = 'free_y')) #|> 
  # ggsave(filename = file.path(out_dir, 'plot', 'TP_val_turbid.png'),
  #        height = 20, width = 20, units = 'cm')


(baseline_runs |> 
    select(Date, uQInEpi, NAME, WBID) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    
    ggplot(aes(x=doy, y = uQInEpi)) + # obs are in ug/L
    geom_line(aes(group = year), alpha = 0.6, colour = 'blue4') + # PClake in g/m3
    theme_bw() +
    facet_wrap(~WBID+NAME, scales = 'free_y')) 


# calculate the TN:TP ratio
#TN:TP (molar) = (TN / 14.007) ÷ (TP / 30.974) 
N_mol <- 14.007
P_mol <- 30.974

baseline_runs |> 
  select(Date, oNTotWEpi, oPTotWEpi, NAME, WBID) |> 
  mutate(N = oNTotWEpi / N_mol,
         P = oPTotWEpi / P_mol,
         `N:P` = N/P) |> 
  reframe(.by = c(NAME,WBID),
          mean_ratio = mean(`N:P`)) |> 
  arrange(mean_ratio)

baseline_runs |> 
  select(Date, oNTotWEpi, oPTotWEpi, NAME, WBID) |> 
  mutate(N = oNTotWEpi / N_mol,
         P = oPTotWEpi / P_mol,
         `N:P` = N/P) |> 
  mutate(lake = fct_reorder(NAME, desc(`N:P`), .fun='median')) |> 
  ggplot(aes(x=lake, y=`N:P`)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0, 600)) + 
  theme_bw() +
  geom_hline(yintercept = 30, linetype = 'dashed')

baseline_runs |> 
  select(Date, uNLoadEpi, uPLoadEpi, NAME, WBID) |> 
  mutate(N = uNLoadEpi / N_mol,
         P = uPLoadEpi / P_mol,
         `N:P` = N/P) |> 
  mutate(lake = fct_reorder(NAME, desc(`N:P`), .fun='median')) |> 
  ggplot(aes(x=lake, y=`N:P`)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0, 600)) +
  theme_bw() +
  geom_hline(yintercept = 30, linetype = 'dashed')
  


