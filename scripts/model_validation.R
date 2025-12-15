#--------------------------------------#
## Project: Lake District nutrient critical limit estimation
## Script purpose: model validation using "observed" nutrient loads
## Date: 2025-12-03
## Author: Freya Olsson
#--------------------------------------#
library(tidyverse)

start_date <- '1998-01-01'

# Read in baseline runs ----------------
baseline_runs <- list.files('output/', pattern = 'baseline', full.names = TRUE) |>  
  lapply(read_csv, show_col_types = F, id = 'NAME_short') |> 
  bind_rows() |> 
  mutate(NAME_short = str_to_upper(str_sub(NAME_short, 17, - 5)),
         Date = as_date(time, origin = start_date)) # figure out how to deal with leap years?
  

# read in observations -----------------
## Secchi depth -----------------------
obs_secchi <- read_csv('data/Validation/Secchi.csv') |> 
  mutate(Site = ifelse(Site %in% c('SBAS', 'NBAS'), 'WIND', Site))

ggplot(obs_secchi, aes(x=Date, y= Diskvalu)) + 
  geom_point() + 
  facet_wrap(~Site, scales = 'free_y') 

(baseline_runs |> 
  select(Date, aSecchiT, NAME_short) |> 
  left_join(obs_secchi, by = join_by(NAME_short == Site, Date == Date)) |> 
  ggplot(aes(x=Date, y = Diskvalu)) +
  geom_line(aes(y = aSecchiT, colour = 'model')) +
  geom_point() +
  facet_wrap(~NAME_short, scales = 'free_y')) |> 
  ggsave(filename = 'output/plots/val_secchi.png', width = 30, height = 20, units = 'cm')

## Chlorophyll -----------------------
# also requires teh sample.csv that determines where the sample was taken
# regular samples are integrated top 5/7 m at deepest location (7 in Windermere only), or DIP = shoreline

samples_chem <- read_csv('data/Validation/samples.csv', show_col_types = F) |> 
  mutate(site = ifelse(site %in% c('SBAS', 'NBAS'), 'WIND', site))

obs_chl <- read_csv('data/Validation/chemistry.csv', show_col_types = F) |>
  mutate(site = ifelse(site %in% c('SBAS', 'NBAS'), 'WIND', site)) |> 
  filter(variable == 'TOCA') |> # total chlorophyll
  left_join(samples_chem,
            by = join_by(date, site), 
            relationship = 'many-to-many') |> 
  distinct() 

ggplot(obs_chl, aes(x=date, y= chemvalu)) + 
  geom_point() + 
  facet_wrap(~site, scales = 'free_y') 

(baseline_runs |> 
  select(Date, oChlaEpi, NAME_short) |> 
  left_join(obs_chl, by = join_by(NAME_short == site, Date == date)) |> 
  ggplot(aes(x=Date, y = chemvalu)) +
  geom_line(aes(y = oChlaEpi, colour = 'model')) +
  geom_point(size = 0.9, alpha = 0.6) +
  facet_wrap(~NAME_short, scales = 'free_y')) |> 
  ggsave(filename = 'output/plots/val_chla.png', width = 30, height = 20, units = 'cm')

## Phosphorus -----------------------
obs_TP <- read_csv('data/Validation/chemistry.csv', show_col_types = F) |> 
  mutate(site = ifelse(site %in% c('SBAS', 'NBAS'), 'WIND', site)) |> 
  filter(variable == 'TOTP') |> # total phosphorus
  left_join(samples_chem,
            by = join_by(date, site), 
            relationship = 'many-to-many') |> 
  distinct() 


ggplot(obs_TP, aes(x=date, y= chemvalu)) + 
  geom_point() + 
  facet_wrap(~site, scales = 'free_y') 

(baseline_runs |> 
  select(Date, oPTotWEpi, NAME_short) |> 
  left_join(obs_TP, by = join_by(NAME_short == site, Date == date)) |> 
  ggplot(aes(x=Date, y = chemvalu)) +
  geom_line(aes(y = oPTotWEpi *1000, colour = 'model')) +
  geom_point(size = 0.9, alpha = 0.6) +
  facet_wrap(~NAME_short, scales = 'free_y')) |> 
  ggsave(filename = 'output/plots/val_TP.png', width = 30, height = 20, units = 'cm')

