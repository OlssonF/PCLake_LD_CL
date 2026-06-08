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

# to save the data locally or not?
archive <- TRUE

# within 50km of centre of Cumbria
lat  <- 54.51
long <- -3.16
dist <- 30   # 30 km radius

# --------------------------------------------#
# ------------ Discharge scaling --------------
# --------------------------------------------#

url_stations <- paste0(
  "http://environment.data.gov.uk/hydrology/id/stations?_limit=1000",
  "&dist=", dist,
  "&lat=", lat,
  "&long=", long
)

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


#' get_flow_data
#'
#' @param station_id stationGuid from the stations endpoint
#'
#' @returns data frame
#' @export
#'
#' @examples
get_flow_data <- function(station_id) {
  
  message(station_id)
  # 1. Get station metadata
  station_url <- paste0(
    "https://environment.data.gov.uk/hydrology/id/measures/", 
    station_id,
    "-flow-m-86400-m3s-qualified/readings.json?mineq-date=2016-10-01&max-date=2025-10-01"
  )
  
  res <- GET(station_url)
  raw <- content(res, as = "text", encoding = "UTF-8")
  station_json <- fromJSON(raw, flatten = TRUE)
  
  # 2. Identify flow measure 
  measures <- station_json$items
  
  
  if (length(measures) != 0) {
    df <- measures %>%
      mutate(station_id = station_id,
             datetime = ymd_hms(dateTime)) %>%
      select(station_id, datetime, value, quality) 
    if (nrow(df) < 8*365) {
      return(NULL)
    } else {
      return(df)
    }
  } else {
    return(NULL)
  }
  
  
}

# Apply to all station IDs
all_ld_flow <- map(station_ids_flow, get_flow_data) |> 
  list_rbind() |> 
  left_join(stations_with_flow, by = join_by(station_id == stationGuid))

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

if (archive) {
  all_ld_flow |> 
    write_csv("data/EArivers_archive/flow_measurements.csv")
}


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


full_join(monthly_ld_flow, annual_ld_flow,
          by = join_by(station_id, label, year)) |>
  filter(!is.na(year)) |> 
  mutate(ratio = monthly_average/annual_average) |> 
  reframe(.by = c(station_id, label, month),
          mean_ratio = mean(ratio)) |> 
  # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
  ggplot(aes(x=month, y=mean_ratio)) +
  geom_col() +
  facet_wrap(~label)

discharge_ratios <- full_join(monthly_ld_flow, annual_ld_flow,
                              by = join_by(station_id, label, year)) |>
  filter(!is.na(year)) |> 
  mutate(ratio = monthly_average/annual_average) |> 
  reframe(.by = c(station_id, label, month),
          mean_ratio = mean(ratio))

discharge_ratios |> 
  # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
  ggplot(aes(x=month, y=mean_ratio))  + 
  geom_point() + 
  geom_smooth(method = "gam", formula = y~s(x, bs = 'cc')) + 
  theme_bw() + 
  geom_hline(yintercept = 1, linetype = 'dashed') + 
  scale_x_continuous(labels = month.abb, breaks = 1:12)

# -------------------------------------------------# 
# --------------- Nutrient scaling -----------------
# -------------------------------------------------# 


# Identify the sampling locations
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

# Step 2: function to grab a time series of sampling point and determinands
get_samples <- function(sp, determinand) {
  
  message(sp)
  skip <- 0 # start with first 250
  run <- 0
  rows <- 250
  
  while(rows == 250) { # use while loop to make sure we get sampling locations where there are more than 250 rows
    message("skipping ", skip)
    dat_full <- NULL
    
    url_samples <- paste0("http://environment.data.gov.uk/water-quality/sampling-point/", sp, "/observation?",
                          "skip=", skip,
                          "&limit=250",
                          "&dateFrom=", '2016-10-01',
                          "&dateTo=", '2025-10-01',
                          "&determinand=", determinand)
    res_samp <- GET(url_samples)
    samp <- fromJSON(content(res_samp, "text", encoding = "UTF-8"))
    
    if (length(samp$member) > 0) {
      # extract the results
      dat <- data.frame(datetime = samp$member$phenomenonTime,
                        value = samp$member$hasResult$numericValue) |> 
        mutate(notation = sp, 
               variable = samp$member$observedProperty$altLabel) |> 
        # where the sample is < take the upper bound
        mutate(BLQ = ifelse(str_detect(samp$member$hasSimpleResult, pattern = '<'), T, F),
               value = ifelse(is.na(value), samp$member$hasResult$upperBound, value))
      
      dat_full <- bind_rows(dat_full, dat)
      
      run <- run + 1
      rows <- nrow(dat)
      message(rows, " rows in JSON")
      skip <- run * rows
      
    } else {
      rows <- 0
    }
    
  }
  
  if(rows != 0){
    return(dat_full)
  }
  
}

# combine with the sampling point information
point_det <- expand_grid(sampling_points, determinand = c('0348', '0117'))

all_samples <- map2(point_det$sampling_points,
                    point_det$determinand,
                    get_samples) |> 
  list_rbind() |> 
  left_join(station_ids_wq, by = join_by(notation))


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

if (archive) {
  all_samples |> 
    write_csv("data/EArivers_archive/nutrient_samples.csv")
}


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


full_join(monthly_samples, annual_samples,
          by = join_by(notation, variable, altLabel, year)) |>
  filter(!is.na(year)) |> 
  mutate(ratio = monthly_average/annual_average) |> 
  reframe(.by = all_of(c('notation', 'variable', 'altLabel', 'month')),
          mean_ratio = mean(ratio)) |> 
  # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
  ggplot(aes(x=month, y=mean_ratio, fill = variable))  + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~notation)


full_join(monthly_samples, annual_samples,
          by = join_by(notation, altLabel, variable, year)) |>
  filter(!is.na(year)) |> 
  mutate(ratio = monthly_average/annual_average) |> 
  reframe(.by = all_of(c('notation', 'variable', 'altLabel', 'month')),
          mean_ratio = mean(ratio)) |> 
  # mutate(month = ((month + 8) %% 12) + 1) |> # convert back to the actual month
  ggplot(aes(x=month, y=mean_ratio, fill = variable))  + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~notation)

nut_ratios <- full_join(monthly_samples, annual_samples,
                        by = join_by(notation, altLabel, variable, year)) |>
  filter(!is.na(year)) |> 
  mutate(ratio = monthly_average/annual_average) |> 
  reframe(.by = c(notation, altLabel, variable, month),
          mean_ratio = mean(ratio))


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


# -------------------------------------------------# 
# --------------- Combine scalings -----------------
# -------------------------------------------------# 
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

# 1. calculate the average concentration and discharge scalings
av_scaling_nuts <- nut_ratios |> 
  reframe(.by = c(month, variable),
          nut_scaling = mean(mean_ratio, na.rm = T))

av_scaling_discharge <- discharge_ratios |> 
  reframe(.by = month,
          discharge_scaling = mean(mean_ratio, na.rm = T))

# 2. calculate the combined scaling by multiplying and "normalising" by proportional scaling
combined_scaling <- full_join(av_scaling_nuts, av_scaling_discharge, by = "month") |> 
  arrange(variable, month) |> 
  mutate(total_scaling = nut_scaling * discharge_scaling) |> 
  mutate(.by = variable,
         rescaled = scales::rescale(total_scaling), #rescaled 0-1
         prop_scale = total_scaling / sum(total_scaling))


combined_scaling |>  # the sum adds to 1
  ggplot(aes(x=month, y = prop_scale)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~variable) + 
  scale_x_continuous(labels = month.abb[seq(1,12,2)], breaks = seq(1,12,2))


# 3. apply scalings to SAGIS annual loadings
sagis_loads <- read_csv("data/FW_ Request for information - Ref_ EIR2026_11092GC/lake_loads.csv", show_col_types = F) |> 
  mutate(value = ifelse(is.na(inlake_summary), mean_value, inlake_summary)) |> 
  filter(variable %in% c('Nitrate', 'Total_Phosphorus'))

sagis_monthly <- sagis_loads |> 
  mutate(variable = ifelse(variable == 'Nitrate', 'Nitrate-N', 
                           ifelse(variable == 'Total_Phosphorus', 'Phosphorus-P', NA))) |> 
  full_join(combined_scaling, by = join_by(variable),
            relationship = 'many-to-many') |> 
  mutate(monthly_apportion = value * prop_scale) 

sagis_monthly |> 
  ggplot(aes(x = month, y = monthly_apportion, group = NAME)) +
  geom_point() +
  facet_wrap(~variable, scales = 'free')  + 
  scale_x_continuous(labels = month.abb[seq(1,12,2)], 
                     breaks = seq(1,12,2)) +
  theme_bw()


sagis_monthly |> 
  filter(NAME %in% c("Ullswater", 'Windermere', 'Bowscale Tarn', 'Grasmere')) |> 
  ggplot(aes(x = NAME, y = monthly_apportion, fill = as_factor(month))) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() +
  geom_point(aes(y = value))


sagis_monthly |> 
  filter(NAME %in% c("Ullswater", 'Windermere', 'Bowscale Tarn', 'Grasmere'))
