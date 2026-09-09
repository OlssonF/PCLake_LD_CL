#--------------------------------------#
## Project: Lake District pathways modelling
## Script purpose: builds the combined Lake District water quality database
## Date: 2026-09-04; 2026-09-09
## Author: Sven Teurlincx (updates FO to match existing project directory)
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

# ============================================================================
# from 00_run_all.R
# Top-level script: builds the combined Lake District water quality database
# from all sources described in Validation/README.txt.
#
# INPUT (expected folder layout, under data_dir):
#   <data_dir>/
#     Lake District_UKCEH Portal data_raw.xlsx
#     SITE ID_MULTIPLE DATA SOURCES_LD LAKES_EBM.xlsx
#     Lake District_UKCEH Portal RT_data.csv
#     Validation/
#       metadata.csv, README.txt
#       EA_WQ/LDFiltered_20XX.csv            (one or more, any year range)
#       UKCEH_LakesTour/Lakes_Tour_Chem_TeOx.xlsx
#       UKCEH_Monitoring/chemistry.csv, samples.csv, Secchi.csv,
#                        temperature and oxygen.csv
#       UKCEH_UWMN/UWMN_Scoat_Burnmoor.xlsx
#       Hydroscape/hydroscape_chemistry_lakes_eidc.csv,
#                  hydroscape_chemistry_uplands_eidc.csv
#
# OUTPUT (written to <validation_dir>/output/):
#   LD_combined_database.csv   - the harmonised long-format time series
#   LD_lake_metadata.csv       - static per-lake attributes, keyed by site_id (WBID)
#   unmatched_sites_log.csv    - every site/sub-site code that could NOT be
#                                matched to a master WBID, with row counts and
#                                date ranges, for manual review
#   variable_mapping_lookup.csv - the variable/unit harmonisation table used,
#                                for review/editing without touching the code
#   duplicate_records_removed.csv - every row dropped as a cross-source
#                                duplicate (same site/date/variable/value
#                                reported by more than one raw database),
#                                with which source was kept vs removed
#   LD_combined_database_core_variables.csv - the same combined table,
#                                filtered to only the 10 variables documented
#                                in metadata.csv (TEMP, OXYG, SECC, ALKA,
#                                NH4N, TON, PO4P, TOTP, SIO2, TOCA)
#
# NOTE ON FILE NAMES: paths below use data_dir's actual local file names
# (with spaces, matching Sven's OneDrive folder). If any file gets renamed
# locally, update the corresponding path_* line below to match.
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(rlang)
})

# ---- paths (EDIT THESE if your folder names differ) ------------------------
project_dir <- here::here() 
data_dir       <- file.path(project_dir, "data") # project folder as the top level
validation_dir <- file.path(data_dir, "Validation")
output_dir     <- file.path(data_dir, "Validation/output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

path_sitemap <- file.path(data_dir, "SITE ID_MULTIPLE DATA SOURCES_LD LAKES_EBM.xlsx")
path_portal  <- file.path(data_dir, "Lake District_UKCEH Portal data_raw.xlsx")
path_rt      <- file.path(data_dir, "Lake District_UKCEH Portal RT_data.csv")

path_chemistry <- file.path(validation_dir, "UKCEH_Monitoring", "chemistry.csv")
path_secchi    <- file.path(validation_dir, "UKCEH_Monitoring", "Secchi.csv")
path_tempox    <- file.path(validation_dir, "UKCEH_Monitoring", "temperature and oxygen.csv")
path_lakestour <- file.path(validation_dir, "UKCEH_LakesTour", "Lakes_Tour_Chem_TeOx.xlsx")
path_uwmn      <- file.path(validation_dir, "UKCEH_UWMN", "UWMN_Scoat_Burnmoor.xlsx")
path_ea_folder <- file.path(validation_dir, "EA_WQ")
path_hydroscape_lakes   <- file.path(validation_dir, "Hydroscape", "hydroscape_chemistry_lakes_eidc.csv")
path_hydroscape_uplands <- file.path(validation_dir, "Hydroscape", "hydroscape_chemistry_uplands_eidc.csv")
path_metadata <- file.path(validation_dir, "metadata.csv")

# ---- load pipeline functions (all functions live in one file) --------------
source(file.path(project_dir, "R", "LD_database_functions.R"))

# ---- 1. master site lookup + lake metadata ----------------------------------
message("Building master site lookup...")
site_lookup <- build_site_lookup(path_sitemap)

message("Building lake metadata table...")
lake_metadata <- build_lake_metadata(path_portal, path_rt)

# ---- 2. variable/unit mapping + hydroscape sample-name mapping --------------
var_mapping <- build_variable_mapping()
hydroscape_sample_mapping <- build_hydroscape_sample_mapping()

# ---- 3. load every source ----------------------------------------------------
message("Loading UKCEH_Monitoring chemistry.csv...")
d_chem <- load_ukceh_chemistry(path_chemistry, site_lookup, var_mapping)

message("Loading UKCEH_Monitoring Secchi.csv...")
d_secchi <- load_ukceh_secchi(path_secchi, site_lookup, var_mapping)

message("Loading UKCEH_Monitoring temperature and oxygen.csv...")
d_tempox <- load_ukceh_tempox(path_tempox, site_lookup, var_mapping)

message("Loading UKCEH_LakesTour CHEM sheet...")
d_lt_chem <- load_lakestour_chem(path_lakestour, site_lookup, var_mapping, lakestour_chem_name_aliases)

message("Loading UKCEH_LakesTour TeOx sheet...")
d_lt_teox <- load_lakestour_teox(path_lakestour, site_lookup, var_mapping)

message("Loading UKCEH_UWMN...")
d_uwmn <- load_uwmn(path_uwmn, site_lookup, var_mapping)

message("Loading EA_WQ LDFiltered_YYYY.csv files...")
d_ea <- load_ea_wq_folder(path_ea_folder, site_lookup, var_mapping)

message("Loading Hydroscape lakes chemistry...")
d_hydro_lakes <- load_hydroscape(path_hydroscape_lakes, "lakes", site_lookup, lake_metadata, var_mapping, hydroscape_sample_mapping)

message("Loading Hydroscape uplands chemistry...")
d_hydro_uplands <- load_hydroscape(path_hydroscape_uplands, "uplands", site_lookup, lake_metadata, var_mapping, hydroscape_sample_mapping)

# ---- 4. combine + depth aggregation ------------------------------------------
message("Combining all sources...")
combined_raw <- bind_rows(d_chem, d_secchi, d_tempox, d_lt_chem, d_lt_teox, d_uwmn, d_ea,
                          d_hydro_lakes, d_hydro_uplands)

message("Aggregating depth-resolved profiles (surface / bottom / whole-lake mean)...")
combined_pre_dedup <- aggregate_depth_profiles(combined_raw)

message("Checking for cross-source duplicate observations...")
dedup_result <- deduplicate_records(combined_pre_dedup)
combined <- dedup_result$deduplicated %>%
  arrange(site_id, date, source, variable_code, depth_zone)
duplicate_records_removed <- dedup_result$removed_log

# ---- 5. unmatched-site log ---------------------------------------------------
unmatched_log <- collect_unmatched_log()

# ---- 5b. metadata.csv "core variables" subset --------------------------------
# metadata.csv only documents 10 of the ~140 variable codes in the combined
# database (the original UKCEH core suite: TEMP, OXYG, SECC, ALKA, NH4N, TON,
# PO4P, TOTP, SIO2, TOCA). This filters the final combined table down to just
# those, across every source - useful as a "known-quantity" subset separate
# from the many source-specific/unstandardised variables also present.
metadata_variables <- read.csv(path_metadata)$UKCEH_code
combined_core_variables <- combined %>% filter(variable_code %in% metadata_variables)

# ---- 6. write outputs ---------------------------------------------------------
message("Writing outputs to ", output_dir, " ...")
write_csv(combined, file.path(output_dir, "LD_combined_database.csv"))
write_csv(lake_metadata, file.path(output_dir, "LD_lake_metadata.csv"))
write_csv(unmatched_log, file.path(output_dir, "unmatched_sites_log.csv"))
write_csv(var_mapping, file.path(output_dir, "variable_mapping_lookup.csv"))
write_csv(duplicate_records_removed, file.path(output_dir, "duplicate_records_removed.csv"))
write_csv(combined_core_variables, file.path(output_dir, "LD_combined_database_core_variables.csv"))

# ---- 7. summary ---------------------------------------------------------------
message("\n=========== SUMMARY ===========")
message("Combined database rows: ", nrow(combined))
message("Date range: ", paste(range(combined$date, na.rm = TRUE), collapse = " to "))
message("Distinct lakes (matched, site_id not NA): ", n_distinct(na.omit(combined$site_id)))
message("Rows with NO matched site_id: ", sum(is.na(combined$site_id)),
        " (", round(100 * mean(is.na(combined$site_id)), 1), "%)")
message("Rows by source:")
print(combined %>% count(source, sort = TRUE))
message("Unmatched site/sub-site codes logged: ", nrow(unmatched_log),
        " (see unmatched_sites_log.csv)")
message("Cross-source duplicate observations removed: ", nrow(duplicate_records_removed),
        " (see duplicate_records_removed.csv)")
message("Core (metadata.csv) variable rows: ", nrow(combined_core_variables),
        " across ", n_distinct(combined_core_variables$variable_code), " variables",
        " (see LD_combined_database_core_variables.csv)")
message("Lake metadata rows: ", nrow(lake_metadata))
message("================================\n")
