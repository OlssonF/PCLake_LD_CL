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
pclake_results <- list.files(out_dir, pattern = '*.csv', recursive = T, full.names = T)

subset <- T # run all lakes or not

# the lakes portal data has all the basic info we need
lakes_portal_df <- read_csv('data/lakes4PCLake.csv', show_col_types = F)

if (subset == T) {
  # LakeIDs to loop through
  lakeIDs <- readxl::read_xlsx('data/SITE ID_MULTIPLE DATA SOURCES_LD LAKES.xlsx') |> 
    filter(!is.na(site_chemistry.csv))   # remove Loughrigg because it's not right
  
  lake_names_lookup <- lakeIDs$`WBID_Lake District_UKCEH Portal data_raw.xlsx`
  names(lake_names_lookup) <- str_extract(lakeIDs$site_chemistry.csv, "....")
  
  lake_names_lookup[which(names(lake_names_lookup) == 'Wind')] <- 29233 
  # the NBAS and SBAS have seperate WBIDs but the one from the lakes portal has the combined one which is different
  
  pclake_results <- pclake_results[str_detect(string = pclake_results, paste(lake_names_lookup, collapse = "|"))]
  
} else {
  
  lake_names_lookup <- lakes_portal_df$WBID
  names(lake_names_lookup) <- lakes_portal_df$NAME
  
  
}

WIND <- c(47007, 47008) # these are the WBIDs for the NBAS and SBAS
WIND_WBID <- 29233

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
  left_join(data.frame(WBID = lake_names_lookup,
                       NAME = names(lake_names_lookup)),
            by = join_by(WBID))

# read in observations -----------------------

obs_secchi <- read_csv(file.path(val_dir, 'Secchi.csv'), show_col_types = F, name_repair = ) |> 
  rename_with(tolower) |> 
  mutate(date = ymd(date)) |> 
  mutate(site = ifelse(site %in% c('SBAS', 'NBAS'), 'WIND', site)) |> 
  left_join(data.frame(WBID = lake_names_lookup,
                       site = names(lake_names_lookup)),
            by = join_by(site)) |> 
  filter(between(date, as_date('2000-01-01'), as_date('2024-12-31')))

obs_tempDO <- read_csv(file.path(val_dir, 'temperature and oxygen.csv'), show_col_types = F) |> 
  rename_with(tolower) |> 
  mutate(date = dmy(date)) |> 
  mutate(site = ifelse(site %in% c('SBAS', 'NBAS'), 'WIND', site)) |> 
  left_join(data.frame(WBID = lake_names_lookup,
                       site = names(lake_names_lookup)),
            by = join_by(site)) |> 
  filter(between(date, as_date('2000-01-01'), as_date('2024-12-31')))

# also requires the sample.csv that determines where the sample was taken
# regular samples are integrated top 5/7 m at deepest location (7 in Windermere only), or DIP = shoreline

samples_chem <- read_csv(file.path(val_dir, 'samples.csv'), show_col_types = F) |> 
  mutate(site = ifelse(site %in% c('SBAS', 'NBAS'), 'WIND', site))

chem_vars <- c('TOCA') 

obs_chem <- read_csv(file.path(val_dir, 'chemistry.csv'), show_col_types = F) |>
  mutate(date = ymd(date)) |> 
  mutate(site = ifelse(site %in% c('SBAS', 'NBAS'), 'WIND', site)) |> 
  filter(between(date, as_date('2000-01-01'), as_date('2024-12-31'))) |> 
  left_join(data.frame(WBID = lake_names_lookup,
                       site = names(lake_names_lookup)),
            by = join_by(site)) |> 
  left_join(samples_chem,
            by = join_by(date, site), 
            relationship = 'many-to-many') |> 
  distinct() 


# plots

(baseline_runs |> 
    select(Date, aSecchiT, NAME, WBID) |> 
    left_join(obs_secchi, by = join_by(NAME == site, Date == date, WBID == WBID)) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    ggplot(aes(x=doy, y = diskvalu)) +
    geom_line(aes(y = aSecchiT, group = year), alpha = 0.3) +
    geom_point(size = 0.9, alpha = 0.3) +
    theme_bw() +
    facet_wrap(~WBID+NAME, scales = 'free_y', 
               nrow = 6, ncol = 3)) |> 
  ggsave(filename = file.path(out_dir, 'plot', 'secchi_val.png'), 
         height = 15, width = 20, units = 'cm')

(baseline_runs |> 
    select(Date, oChlaEpi, NAME, WBID) |> 
    left_join(filter(obs_chem, variable == 'TOCA'), 
              by = join_by(NAME == site, Date == date, WBID ==WBID)) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 
    ggplot(aes(x=doy, y = chemvalu)) + # obs are in ug/L
    geom_line(aes(y = oChlaEpi, group = year), alpha = 0.3) + # PClake in mg/m3
    geom_point(size = 0.9, alpha = 0.6) +
    theme_bw() +
    facet_wrap(~WBID+NAME, scales = 'free_y', 
               nrow = 6, ncol = 3)) |> 
  ggsave(filename = file.path(out_dir, 'plot', 'chla_val.png'), 
         height = 15, width = 20, units = 'cm')


(baseline_runs |> 
    select(Date, oPTotWEpi, NAME, WBID) |> 
    left_join(filter(obs_chem, variable == 'TOTP'), 
              by = join_by(NAME == site, Date == date, WBID ==WBID)) |> 
    mutate(doy = yday(Date),
           year = year(Date)) |> 

    ggplot(aes(x=doy, y = chemvalu/1000)) + # obs are in ug/L
    geom_line(aes(y = oPTotWEpi, group = year), alpha = 0.3) + # PClake in g/m3
    geom_point(size = 0.9, alpha = 0.6) +
    theme_bw() +
    facet_wrap(~WBID+NAME, scales = 'free_y', 
               nrow = 6, ncol = 3)) |> 
  ggsave(filename = file.path(out_dir, 'plot', 'TP_val.png'), 
         height = 20, width = 20, units = 'cm')


