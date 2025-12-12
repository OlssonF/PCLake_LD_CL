#--------------------------------------#
## Project: Lake District nutrient critical limit estimation
## Script purpose: Accessing the SAGIS data and matching to the lakes
## Date: 2025-12-09
## Author: Freya Olsson
#--------------------------------------#
library(tidyverse)
library(readxl)
# Read in data ---------------------------------

# Read in SAGIS data
P_dat <- read_csv('data/SAGIS/data/PO4_SAGIS.csv') |>  mutate(Nutrient = 'P') # these concentrations are in mg/L (or g/m3)
N_dat <- read_csv('data/SAGIS/data/NO3_SAGIS.csv') |>  mutate(Nutrient = 'N')

nut_dat <- bind_rows(P_dat, N_dat) 

lakes_portal_df <- read_csv('data/lakes4PCLake.csv')

lakeIDs <- read_xlsx('data/SITE ID_MULTIPLE DATA SOURCES_LD LAKES.xlsx') |> 
  filter(!is.na(LAKE_LakesTour2021_Zooplankton.csv),
         !str_detect(LAKE_LakesTour2021_Zooplankton.csv, 'Loughrigg'))   # remove Loughrigg because it's not right

lake_names_lookup <- str_extract(lakeIDs$LAKE_LakesTour2021_Zooplankton.csv, "....")
names(lake_names_lookup) <- lakeIDs$`WBID_Lake District_UKCEH Portal data_raw.xlsx`

names(lake_names_lookup)[which(lake_names_lookup == 'Wind')] <- 29233 
# the NBAS and SBAS have seperate WBIDs but the one from the lakes portal has the combined one which is different

## Match sites ----------------------------------------
# Match inflow (SAGIS) with Lakes (portal/tour data)

# Loop through the two nutrients and the inflows  
lake_inflows <- NULL

# look up the the inflows for each lake
for (i in 1:length(lake_names_lookup)) {
  
  lake_name <- lake_names_lookup[i]
  
  # get teh northing and easting
  northing_use <- lakes_portal_df |> 
    filter(str_detect(NAME,lake_name)) |> 
    select(WBNORTH) |> pull()
  
  easting_use <- lakes_portal_df |> 
    filter(str_detect(NAME, lake_name)) |> 
    select(WBEAST) |> pull()
  
  # First check if one is defined in the dataframe
  lake_ins <- nut_dat |>
    mutate(distanceSum = abs(X- easting_use) + abs(Y- northing_use)) |>
    filter(str_detect(FeatName, 'Lake In')) |> 
    filter(str_detect(FeatName, lake_name),
           distanceSum <= 10000) |> # max 10 km away? 
    mutate(WBID = as.numeric(str_split_i(gsub(pattern = "GB3..", x= FeatName, replacement = ""), ' - ', 2))) 
  
  # If not, look another way - look for the end of a reach but not in the lake
  if (nrow(lake_ins) == 0) {
    # message(lake_name)
    lake_ins <- nut_dat |>
      mutate(distanceSum = abs(X- easting_use) + abs(Y- northing_use)) |>
      # select(distanceSum, FeatName, ReachName, ReachNo, MeanConc_m, X, Y) |>
      filter(FeatName  == 'End of reach' & !str_detect(ReachName, 'Lake')) |>
      slice_min(distanceSum, n = 1)
      
  }
  
  if (nrow(lake_ins) == 0) {
    message(lake_name, ' missing')
  }
  
  lake_inflows <- bind_rows(lake_inflows, mutate(lake_ins, NAME_short = lake_name, WBID = as.numeric(names(lake_name))))
  
}


## Write output ------------------------------------
# write out the relevent columns for the baseline runs
lake_inflows |> 
  select(WBID, X, Y, NAME_short, MeanConc_m, LCLimMnCon, UCLimMnCon, Nutrient) |>
  write_csv('data/SAGIS/data/InflowNutrientConcs_SAGIS.csv')
