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

# Step 2: ffunction to grab a time series of sampling point and determinand
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
               variable = samp$member$observedProperty$altLabel)
      
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
all_samples <- map(sampling_points, get_samples, 
                   determinand = '0348') |> 
  list_rbind() |> 
  left_join(station_ids_wq, by = join_by(notation))


all_samples |> 
  mutate(datetime = ymd_hms(datetime),
         date = ymd(format(datetime, "%Y-%m-%d"))) |>
  ggplot(aes(x=date, y= value)) +
  geom_point() +
  facet_wrap(~notation, scales = 'free_y')

