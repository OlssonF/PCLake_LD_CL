#' extract_SAGIS
#'
#' @param variable the name of nutrient table to query e.g. Total_Phosphorus
#' @param inlake_summary which summarising function to use (mean, median etc.)
#' @param portion which statistic to sample (meanLdKGd, or another apportionment)
#'
#' @returns
#' @export
#'
#' @examples
extract_SAGIS <- function(variable = 'Total_Phosphorus', 
                          inlake_summary = mean, 
                          portion = 'MeanLdKGd') {
  
  library(sf)
  library(tidyverse)
  library(readxl)
  # ------------------------------------------------------------#
  # 1. Load data
  # ------------------------------------------------------------#
  
  # Replace with your file paths
  ld_lakesportal <- readr::read_csv('data/Lake District_UKCEH Portal RT_data.csv', 
                                    show_col_types = F) |>
    select(WBID)
  
  polys <- st_read("data/uklakes/data/uklakes_v3_6_poly.gpkg", 
                   quiet = T) |> 
    filter(WBID %in% ld_lakesportal$WBID) |> # filter to just the lake district lakes
    select(WBID, NAME)
  
  catchments <- st_read("data/uklakes/data/uklakes_v3_7_catchments.gpkg", quiet = T) |> 
    filter(WBID %in% ld_lakesportal$WBID) |> # filter to just the lake district lakes
    select(WBID, NAME)
  
  SAGIS_files <- list.files("data/FW_ Request for information - Ref_ EIR2026_11092GC/", 
                            pattern = paste0("\\", variable, ".xlsx$"), full.names = TRUE)
  
  message("reading ", paste(basename(SAGIS_files), sep = ";"))
  
  pts <- map_df(SAGIS_files,
                read_excel) |> 
    select(all_of(c("X", "Y", "OBJECTID", portion))) |> 
    st_as_sf(coords = c("X", "Y"), 
             crs = 27700)
  
  # ------------------------------------------------------------#
  # 2. Spatial join: attach polygon attributes to each point
  # ------------------------------------------------------------#
  
  # This gives each point the polygon it falls inside
  joined <- st_join(pts, polys, left = F, join = st_within) 
  
  # ------------------------------------------------------------#
  # 3. Summary of intersecting points per polygon
  # ------------------------------------------------------------#
  
  summary_tbl <- joined %>%
    st_drop_geometry() %>%
    group_by(WBID, NAME) %>%
    summarise(#inlake_n   = sum(!is.na(get(portion))),
              inlake_summary = inlake_summary(get(portion)), .groups = 'drop') |> 
    pivot_longer(inlake_summary, names_to = 'source', values_to = 'value_kg_day') |> 
    mutate(distance_m = 0)
  
  # ------------------------------------------------------------#
  # 4. Identify polygons with NO intersecting points
  # ------------------------------------------------------------#
  
  no_points <- polys %>%
    anti_join(summary_tbl, by = "WBID") 
  # anti_join() return all rows from x without a match in y.
  
  
  # ------------------------------------------------------------#
  # 5. Filter SAGIS points to those inside catchments
  # ------------------------------------------------------------#
  
  pts_in_catch <- st_join(pts, catchments, join = st_within, left = FALSE) %>%
    select(OBJECTID, all_of(portion), WBID, NAME)
  
  # ------------------------------------------------------------#
  # 6. For each polygon with no points:
  #    Try nearest point inside catchment
  #    If none, fall back to nearest point from all SAGIS points
  # ------------------------------------------------------------#
  
  nearest_results <- map_df(no_points$WBID, function(wbid) {
    
    lake_poly <- polys %>% filter(WBID == wbid)
    
    # 6a. Points inside matching catchment
    pts_subset <- pts_in_catch %>% filter(WBID == wbid)
    
    if (nrow(pts_subset) > 0) {
      # nearest point inside catchment
      idx <- st_nearest_feature(lake_poly, pts_subset)
      nearest_pt <- pts_subset[idx, ]
      dist <- st_distance(st_geometry(lake_poly), st_geometry(nearest_pt))
      tibble(
        WBID = wbid,
        NAME = lake_poly$NAME,
       # nearest_id = nearest_pt$OBJECTID,
        value_kg_day = nearest_pt[[portion]],
        distance_m = as.numeric(dist),
        source = "nearest_catchment"
      )
    } else {
      # 6b. fallback: nearest point from ALL SAGIS points
      idx <- st_nearest_feature(lake_poly, pts)
      nearest_pt <- pts[idx, ]
      dist <- st_distance(st_geometry(lake_poly), st_geometry(nearest_pt))
      tibble(
        WBID = wbid,
        NAME = lake_poly$NAME,
        #nearest_id = nearest_pt$OBJECTID,
        value_kg_day = nearest_pt[[portion]],
        distance_m = as.numeric(dist),
        source = "nearest_all"
      )
    }
  })
  
  # ------------------------------------------------------------#
  # 7. Combine outputs
  # ------------------------------------------------------------#
  
  out <- bind_rows(summary_tbl, nearest_results) %>%
    mutate(variable = variable)
  
  return(out)
}



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


#' get_samples
#'
#' @param sp sample point (e.g. NW-SSNxxxx)
#' @param determinand code of the determinand (P = 0348, N = 0118)
#'
#' @returns dataframe
#' @export
#'
#' @examples
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
