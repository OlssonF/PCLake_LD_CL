#--------------------------------------#
## Project: Lake District PCLake
## Script purpose: Access the EA API to fetch the river discharge and utrient concentrations data to 
#                   generate scaling factors to convert annual SAGIS loads to monthly loads
## Date: 2026-06-05
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

library(httr)
library(jsonlite)
library(tidyverse)
library(purrr)

source('R/sagis_funs.R')

# to save the data locally or not?
archive <- TRUE
plot <-  FALSE
# within 50km of centre of Cumbria
lat  <- 54.51
long <- -3.16
dist <- 30   # 30 km radius

# --------------------------------------------#
# --------- 1. Extract SAGIS from GIS ------------
# --------------------------------------------#
nutrients <- c("Ammonia",
               "Nitrate",
               "Phosphate",
               "Total_Phosphorus")

# 1. Extract SAGIS data from the files 
sagis_nuts <- map(.x = nutrients, .f = extract_SAGIS, 
                  inlake_summary = mean, 
                  portion = 'MeanLdKGd') |> 
  list_rbind()

if (archive) {
  sagis_nuts |> write_csv("data/FW_ Request for information - Ref_ EIR2026_11092GC/lake_loads.csv")
}


if (plot) {
  # plot this by getting the polygons again
  ld_lakesportal <- readr::read_csv('data/Lake District_UKCEH Portal RT_data.csv') |> 
    select(WBID)
  polys <- st_read("data/uklakes/data/uklakes_v3_6_poly.gpkg") |> 
    filter(WBID %in% ld_lakesportal$WBID) # filter to just the lake district lakes
  
  polys_summary <- polys %>%
    left_join(sagis_nuts, by = join_by(WBID, NAME))
  
  
  ggpubr::ggarrange(
    polys_summary |> 
      filter(variable  == nutrients[1]) |> 
      ggplot() +
      geom_sf(aes(fill = value_kg_day), colour = NA) +
      scale_fill_viridis_c(option = "plasma", na.value = "grey") +
      theme_minimal() +
      labs(title = nutrients[1], 
           fill = "Mean loading per lake"),
    
    polys_summary |> 
      filter(variable  == nutrients[2]) |> 
      ggplot() +
      geom_sf(aes(fill = value_kg_day), colour = NA) +
      scale_fill_viridis_c(option = "viridis", na.value = "grey") +
      theme_minimal() +
      labs(title = nutrients[2],
           fill = "Mean loading per lake"),
    
    polys_summary |> 
      filter(variable  == nutrients[3]) |> 
      ggplot() +
      geom_sf(aes(fill = value_kg_day), colour = NA) +
      scale_fill_viridis_c(option = "mako", na.value = "grey") +
      theme_minimal() +
      labs(title = nutrients[3], 
           fill = "Mean loading per lake"),
    
    polys_summary |> 
      filter(variable  == nutrients[4]) |> 
      ggplot() +
      geom_sf(aes(fill = value_kg_day), colour = NA) +
      scale_fill_viridis_c(option = "rocket", na.value = "grey") +
      theme_minimal() +
      labs(title = nutrients[4], 
           fill = "Mean loading per lake per day")
  )
  
  
  polys_summary |> 
    ggplot() +
    geom_histogram(aes(x = value_kg_day),) +
    facet_wrap(~variable, scales = 'free') +
    theme_bw()
  
  
  polys_summary |> 
    ggplot() +
    geom_point(aes(x = POLY_AREA_HA, y = value_kg_day),) +
    facet_wrap(~variable, scales = 'free') +
    theme_bw()
  
}
# ----------------------------------------------- #
# ---- 2. Downscale to a monthly load -------------
# -----------------------------------------------#

# --------------------------------------------#
## ------------ 2a Discharge scaling -------------
# ----------------------------------------------- #

### Get the station info -------------------------
url_stations <- paste0("http://environment.data.gov.uk/hydrology/id/stations?_limit=1000",
                       "&dist=", dist, "&lat=", lat, "&long=", long)

# Fetch using httr
res <- GET(url_stations)

# Parse JSON
raw <- content(res, as = "text", encoding = "UTF-8")
json <- fromJSON(raw, flatten = TRUE)

stations <- json$items

# Filter stations that have FLOW measures
stations_with_flow <- stations %>%
  mutate(has_flow = map_lgl(measures, ~ any(.x$parameter == "flow"))) |> 
  filter(has_flow)

station_ids_flow <- stations_with_flow$stationGuid

##----------------------------------------------#

### Get the flow data -------------------------
# Apply to all station IDs
all_ld_flow <- map(station_ids_flow, get_flow_data) |> 
  list_rbind() |> 
  left_join(stations_with_flow, by = join_by(station_id == stationGuid))

if (plot) {
  all_ld_flow |> 
    mutate(doy = yday(datetime - months(9)), # set the first day as Oct 1
           year =  year(datetime - months(9)) + 1) |>
    ggplot(aes(x= as_date(doy + 274), # get the actual data
               y=value)) +
    # geom_line(aes(group = year)) +
    geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cc")) +
    facet_wrap(~label, scales = 'free_y') +
    scale_x_date(date_labels = "%d %b", name = 'Day of year') +
    theme_bw()
  
  all_ld_flow |> 
    mutate(doy = yday(datetime - months(9)), # set the first day as Oct 1
           year =  year(datetime - months(9)) + 1) |>
    ggplot(aes(x= as_date(doy + 274), # get the actual data
               y=value)) +
    geom_line(aes(group = year)) +
    facet_wrap(~label, scales = 'free_y') +
    scale_x_date(date_labels = "%d %b", name = 'Day of year') +
    theme_bw()
}


if (archive) {
  all_ld_flow |> 
    write_csv("data/EArivers_archive/flow_measurements.csv")
}

##----------------------------------------------#

### Calculate scaling factors -------------------
# ratio between average monthly discharge and average annual discharge
# calculate annual average discharge per river
annual_ld_flow <- all_ld_flow |>
  mutate(# set the first day as Oct 1 (of the next year)
    year = year(datetime)) |> # - months(9)) + 1) |> 
  reframe(.by = c(station_id, label, year),
          annual_average = mean(value, na.rm = T)) 

monthly_ld_flow <- all_ld_flow |>
  mutate(# set the first day as Oct 1
    year = year(datetime),# - months(9)) + 1,
    month = month(datetime)) |> # - months(9))) |> 
  reframe(.by = c(station_id, label, year, month),
          monthly_average = mean(value, na.rm = T)) 

discharge_ratios <- full_join(monthly_ld_flow, annual_ld_flow,
                              by = join_by(station_id, label, year)) |>
  filter(!is.na(year)) |> 
  mutate(ratio = monthly_average/annual_average) |> 
  reframe(.by = c(station_id, label, month),
          mean_ratio = mean(ratio))


if (plot) {
  discharge_ratios |> 
    # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
    ggplot(aes(x=month, y=mean_ratio)) +
    geom_col() +
    facet_wrap(~label)
  
  discharge_ratios |> 
    # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
    ggplot(aes(x=month, y=mean_ratio))  + 
    geom_point() + 
    geom_smooth(method = "gam", formula = y~s(x, bs = 'cc')) + 
    theme_bw() + 
    geom_hline(yintercept = 1, linetype = 'dashed') + 
    scale_x_continuous(labels = month.abb, breaks = 1:12)
  
}

##-------------------------------------------------#

# -------------------------------------------------# 
## ------------- 2b Nutrient scaling ----------------
# -------------------------------------------------# 

### Get sample location info -----------------------

# Can only do 250 at a time
skip <- 0 # start with first 250
station_ids <- NULL
run <- 0
rows <- 250

while (rows == 250) {
  message("skipping ", skip)
  # Fetch using httr
  res <- GET(paste0("https://environment.data.gov.uk/water-quality/sampling-point?",
                    "skip=", skip,
                    "&limit=250&latitude=54.51&longitude=-3.16&radius=30&samplingPointType=F6"))
  
  # Parse JSON
  raw <- content(res, as = "text", encoding = "UTF-8")
  json <- fromJSON(raw, flatten = TRUE)
  
  station_ids_wq <- bind_rows(station_ids, json$member)
  
  run <- run + 1
  rows <- nrow(json$member)
  message(rows, " rows in JSON")
  skip <- run * rows
} 

sampling_points <- station_ids_wq$notation

### Get sample measurements --------------------------
# want both determinands at all sites
point_det <- expand_grid(sampling_points, determinand = c('0348', '0117'))

all_samples <- map2(point_det$sampling_points,
                    point_det$determinand,
                    get_samples) |> 
  list_rbind() |> 
  left_join(station_ids_wq, by = join_by(notation))

if (plot) {
  all_samples |> 
    mutate(datetime = ymd_hms(datetime),
           date = ymd(format(datetime, "%Y-%m-%d"))) |>
    filter(variable == 'Nitrate-N') |> 
    ggplot(aes(x=date, y= value)) +
    geom_point() +
    facet_wrap(~notation, scales = 'free_y')
  
  all_samples |> 
    mutate(datetime = ymd_hms(datetime),
           date = ymd(format(datetime, "%Y-%m-%d"))) |>
    filter(variable == 'Phosphorus-P') |> 
    ggplot(aes(x=date, y= value)) +
    geom_point() +
    facet_wrap(~notation, scales = 'free_y')
}


if (archive) {
  all_samples |> 
    write_csv("data/EArivers_archive/nutrient_samples.csv")
}

##--------------------------------------------------------#

### Calculate scaling factors -----------------------------

# ratio between average monthly concentration and average annual concentration
# calculate annual average concentration per river
annual_samples <- all_samples |>
  mutate(# set the first day as Oct 1 (of the next year)
    year = year(datetime)) |> # - months(9)) + 1) |> 
  reframe(.by = all_of(c('notation', 'variable', 'altLabel', 'year')),
          annual_average = mean(value, na.rm = T)) 

monthly_samples <- all_samples |>
  mutate(# set the first day as Oct 1
    year = year(datetime),# - months(9)) + 1,
    month = month(datetime)) |> # - months(9))) |> 
  reframe(.by = all_of(c('notation', 'variable', 'altLabel', 'year', 'month')),
          monthly_average = mean(value, na.rm = T)) 

nut_ratios <- full_join(monthly_samples, annual_samples,
                        by = join_by(notation, altLabel, variable, year)) |>
  filter(!is.na(year)) |> 
  mutate(ratio = monthly_average/annual_average) |> 
  reframe(.by = c(notation, altLabel, variable, month),
          mean_ratio = mean(ratio))

if (plot) {
  
  nut_ratios |> 
    # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
    ggplot(aes(x=month, y=mean_ratio, fill = variable))  + 
    geom_bar(position="dodge", stat="identity") +
    facet_wrap(~notation)
  
  
  nut_ratios |> 
    # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
    ggplot(aes(x=month, y=mean_ratio, fill = variable))  + 
    geom_bar(position="dodge", stat="identity") +
    facet_wrap(~notation)
  
  nut_ratios |> 
    # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
    ggplot(aes(x=month, y=mean_ratio))  + 
    geom_point() + 
    geom_smooth(method = "gam", formula = y~s(x, bs = 'cc')) + 
    facet_wrap(~variable, nrow = 2) + 
    theme_bw() + 
    geom_hline(yintercept = 1, linetype = 'dashed')  + 
    scale_x_continuous(labels = month.abb[seq(1,12,2)], 
                       breaks = seq(1,12,2))
  
}


# -------------------------------------------------# 
## -------------- 2c Combine scalings ---------------
# -------------------------------------------------# 
if (plot) {
  nut_ratios |> 
    reframe(.by = c(month, variable),
            av_scaling = mean(mean_ratio, na.rm = T)) |> 
    ggplot(aes(x=month, y=av_scaling)) +
    geom_col() + 
    scale_x_continuous(labels = month.abb[seq(1,12,2)], 
                       breaks = seq(1,12,2)) +
    facet_wrap(~variable)
  
  discharge_ratios |> 
    reframe(.by = month,
            av_scaling = mean(mean_ratio, na.rm = T)) |> 
    ggplot(aes(x=month, y=av_scaling)) +
    geom_col() + 
    scale_x_continuous(labels = month.abb[seq(1,12,2)], 
                       breaks = seq(1,12,2))
}


# calculate the average concentration and discharge scalings
av_scaling_nuts <- nut_ratios |> 
  reframe(.by = c(month, variable),
          nut_scaling = mean(mean_ratio, na.rm = T))

av_scaling_discharge <- discharge_ratios |> 
  reframe(.by = month,
          discharge_scaling = mean(mean_ratio, na.rm = T))

# calculate the combined scaling by multiplying and "normalising" by proportional scaling
combined_scaling <- full_join(av_scaling_nuts, av_scaling_discharge, by = "month") |> 
  arrange(variable, month) |> 
  mutate(total_scaling = nut_scaling * discharge_scaling) |> 
  mutate(.by = variable,
         rescaled = scales::rescale(total_scaling), #rescaled 0-1
         prop_scale = total_scaling / sum(total_scaling))

if (plot) {
  combined_scaling |>  # the sum adds to 1
    ggplot(aes(x=month, y = prop_scale)) +
    geom_bar(stat = 'identity') +
    facet_wrap(~variable) + 
    scale_x_continuous(labels = month.abb[seq(1,12,2)], breaks = seq(1,12,2))
  
}

## ----------------------------------------------------#

## ----------------------------------------------------#
## ----------------- 2d Apply scalings ----------------
## ----------------------------------------------------#

# Apply scalings to  SAGIS annual loadings 
sagis_loads <- read.csv("data/FW_ Request for information - Ref_ EIR2026_11092GC/lake_loads.csv") |> 
  filter(variable %in% c('Nitrate', 'Total_Phosphorus')) |> 
  mutate(value_kg_year = value_kg_day * 365) # SAGIS data are in kg/day

sagis_monthly <- sagis_loads |> 
  mutate(variable = ifelse(variable == 'Nitrate', 'Nitrate-N', 
                           ifelse(variable == 'Total_Phosphorus', 'Phosphorus-P', NA))) |> 
  full_join(combined_scaling, by = join_by(variable),
            relationship = 'many-to-many') |> 
  mutate(value_kg_month = value_kg_year * prop_scale) 

if (plot) {
  sagis_monthly |> 
    ggplot(aes(x = month, y = value_kg_month, group = WBID)) +
    geom_point() +
    facet_wrap(~variable, scales = 'free')  + 
    scale_x_continuous(labels = month.abb[seq(1,12,2)], 
                       breaks = seq(1,12,2)) +
    theme_bw()
  
  
  sagis_monthly |> 
    filter(NAME %in% c("Ullswater", 'Windermere', 'Bowscale Tarn', 'Grasmere')) |> 
    ggplot(aes(x = NAME, y = value_kg_month, fill = as_factor(month))) + 
    geom_bar(position="stack", stat="identity") +
    facet_wrap(~variable, scales = 'free') +
    theme_bw() +
    geom_point(aes(y = value_kg_year))
  
}

if (archive) {
  sagis_monthly |> write_csv("data/FW_ Request for information - Ref_ EIR2026_11092GC/lake_loads_monthly.csv")
}

# convert from the month loads to a load per day, equal in the month
year_dates <- distinct(sagis_monthly, WBID, NAME, variable) |> 
  group_by(WBID, NAME, variable) |> 
  reframe(Date = as_date(seq.Date(as_date('2026-01-01'), as_date('2026-12-31'), 'day')))


sagis_daily <- sagis_monthly |> 
  mutate(Date = dmy(paste0('01-',month, '-2026')),
         n_days = days_in_month(Date),
         value_kg_day = value_kg_month/n_days) |> 
  select(Date, WBID, NAME, value_kg_day, variable) |> 
  full_join(year_dates, by = join_by(Date, WBID, NAME, variable)) |> 
  group_by(WBID, NAME, variable) |> 
  arrange(Date) |> 
  mutate(value_kg_day = imputeTS::na_locf(value_kg_day, na_remaining = 'rev')) |> 
  ungroup()|> 
  mutate(day = yday(Date))


# also try generating one with a smooth curve, interpolate to mid-month
sagis_daily_smoothed <- sagis_monthly |> 
  mutate(Date = dmy(paste0('01-',month, '-2026')),
         n_days = days_in_month(Date),
         value_kg_day = value_kg_month/n_days) |> 
  select(Date, WBID, NAME, value_kg_day, variable) |> 
  full_join(year_dates, by = join_by(Date, WBID, NAME, variable)) |> 
  group_by(WBID, NAME, variable) |> 
  arrange(Date) |> 
  mutate(value_kg_day = imputeTS::na_interpolation(value_kg_day, option = 'linear', na.rm = F)) |> 
  ungroup() |> 
  mutate(day = yday(Date))



if (plot) {
  sagis_daily |> 
    filter(NAME %in% c('Bowscale Tarn', 'Crummock Water', 'Bassenthwaite Lake')) |> 
    ggplot(aes(x=Date, y = value_kg_day)) +
    geom_line() +
    scale_x_datetime() +
    facet_grid(WBID+NAME~variable, scales = 'free')
}

sagis_daily |> 
  reframe(daily_av = mean(value_kg_day),
          .by = c('variable', 'NAME', 'WBID')) |> 
  mutate(variable = ifelse(variable == 'Nitrate-N', 'Nitrate',
                           ifelse(variable == 'Phosphorus-P', 
                                  'Total_phosphorus', NA))) |> 
  full_join(sagis_nuts) |> 
  mutate(test = ifelse((daily_av - value_kg_day)/value_kg_day < 0.01, T, F))|> 
  filter(test == F) 
# check that the downscaled daily average matches the daily average from SAGIS - this should be zero!

if (archive) {
  sagis_daily |> 
    select(day, variable, NAME, WBID, value_kg_day) |> 
    mutate(variable = ifelse(variable == 'Nitrate-N', 'nitrate',
                             ifelse(variable == 'Phosphorus-P', 
                                    'phosphorus', NA))) |> 
    write_csv("data/FW_ Request for information - Ref_ EIR2026_11092GC/lake_loads_daily.csv")
}
 