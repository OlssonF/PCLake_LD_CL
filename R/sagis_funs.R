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
  # ------------------------------------------------------------
  # 1. Load data
  # ------------------------------------------------------------
  
  # Replace with your file paths
  ld_lakesportal <- readr::read_csv('data/Lake District_UKCEH Portal RT_data.csv', show_col_types = F) |> select(WBID)
  
  polys <- st_read("data/uklakes/data/uklakes_v3_6_poly.gpkg", 
                   quiet = T) |> 
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
  
  # ------------------------------------------------------------
  # 2. Spatial join: attach polygon attributes to each point
  # ------------------------------------------------------------
  
  # This gives each point the polygon it falls inside
  joined <- st_join(pts, polys, left = F, join = st_within) 
  
  # ------------------------------------------------------------
  # 3. Summary of intersecting points per polygon
  # ------------------------------------------------------------
  
  summary_tbl <- joined %>%
    st_drop_geometry() %>%
    group_by(WBID, NAME) %>%
    summarise(inlake_n   = sum(!is.na(get(portion))),
              inlake_summary = inlake_summary(get(portion)), .groups = 'drop')
  
  # ------------------------------------------------------------
  # 4. Identify polygons with NO intersecting points
  # ------------------------------------------------------------
  
  no_points <- polys %>%
    anti_join(summary_tbl, by = "WBID") 
  # anti_join() return all rows from x without a match in y.
  
  # ------------------------------------------------------------
  # 5. Find nearest point for polygons with no points
  # ------------------------------------------------------------
  
  # Compute nearest point index
  nearest_index <- st_nearest_feature(no_points, pts)
  
  # Extract the nearest points
  nearest_points <- pts[nearest_index, ]
  
  # Compute centroid-to-point distance
  distances <- distances <- st_distance(st_geometry(st_centroid(no_points)),
                                        st_geometry(nearest_points),
                                        by_element = T)
  
  # ------------------------------------------------------------
  # 6. Build final nearest-point table
  # ------------------------------------------------------------
  
  nearest_tbl <- no_points %>%
    st_drop_geometry() %>%
    mutate(#SAGIS_location    = nearest_points$OBJECTID,  
      mean_value = nearest_points$MeanLdKGd,
      distance_m = as.numeric(distances))
  
  
  # ------------------------------------------------------------
  # 7. Combine outputs from the two methods
  # ------------------------------------------------------------
  
  out <- full_join(nearest_tbl, summary_tbl,
                   by = join_by(WBID, NAME)) |> 
    mutate(variable = variable)
  
  return(out)
}



nutrients <- c("Ammonia",
               "Nitrate",
               "Phosphate",
               "Total_Phosphorus")

sagis_nuts <- map(nutrients, extract_SAGIS) |> 
  list_rbind()

sagis_nuts |> write_csv("data/FW_ Request for information - Ref_ EIR2026_11092GC/lake_loads.csv")

# plot this by getting the polygons again
ld_lakesportal <- readr::read_csv('data/Lake District_UKCEH Portal RT_data.csv') |> select(WBID)
polys <- st_read("data/uklakes/data/uklakes_v3_6_poly.gpkg") |> 
  filter(WBID %in% ld_lakesportal$WBID) # filter to just the lake district lakes

polys_summary <- polys %>%
  left_join(sagis_nuts, by = join_by(WBID, NAME)) |> 
  mutate(value = ifelse(is.na(inlake_summary), mean_value, inlake_summary))


ggpubr::ggarrange(
  polys_summary |> 
    filter(variable  == nutrients[1]) |> 
    ggplot() +
    geom_sf(aes(fill = value), colour = NA) +
    scale_fill_viridis_c(option = "plasma", na.value = "grey") +
    theme_minimal() +
    labs(title = nutrients[1], 
         fill = "Mean loading per lake"),
  
  polys_summary |> 
    filter(variable  == nutrients[2]) |> 
    ggplot() +
    geom_sf(aes(fill = value), colour = NA) +
    scale_fill_viridis_c(option = "viridis", na.value = "grey") +
    theme_minimal() +
    labs(title = nutrients[2],
         fill = "Mean loading per lake"),
  
  polys_summary |> 
    filter(variable  == nutrients[3]) |> 
    ggplot() +
    geom_sf(aes(fill = value), colour = NA) +
    scale_fill_viridis_c(option = "mako", na.value = "grey") +
    theme_minimal() +
    labs(title = nutrients[3], 
         fill = "Mean loading per lake"),
  
  polys_summary |> 
    filter(variable  == nutrients[4]) |> 
    ggplot() +
    geom_sf(aes(fill = value), colour = NA) +
    scale_fill_viridis_c(option = "rocket", na.value = "grey") +
    theme_minimal() +
    labs(title = nutrients[4], 
         fill = "Mean loading per lake")
)


polys_summary |> 
  ggplot() +
  geom_histogram(aes(x = value),) +
  facet_wrap(~variable, scales = 'free') +
  theme_bw()


polys_summary |> 
  ggplot() +
  geom_point(aes(x = POLY_AREA_HA, y = value),) +
  facet_wrap(~variable, scales = 'free') +
  theme_bw()
