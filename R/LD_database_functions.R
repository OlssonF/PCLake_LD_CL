#--------------------------------------#
## Project: Lake District pathways modelling
## Script purpose: unctions used to build the combined Lake District water quality
## Date: 2026-09-04
## Author: Sven Teurlincx
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

# ============================================================================
# LD_database_functions.R
# All functions used to build the combined Lake District water quality
# database. Called from 00_compile_obs.R - not meant to be run on its own.
#
# Sections:
#   1. Site lookup (master WBID table + Lakes Tour CHEM name aliases)
#   2. Lake metadata (static morphometry/catchment/retention table)
#   3. Variable/unit mapping lookup (cross-source harmonisation table)
#   4. Source loaders (one per file/sheet -> common long schema)
#   5. Depth aggregation (surface / bottom / whole-lake mean)
# ============================================================================


# ============================================================================
# 01_site_lookup.R
# Builds the master site lookup table (WBID = master site ID) from
# SITE_ID_MULTIPLE_DATA_SOURCES_LD_LAKES.xlsx, forward-filling lake-level
# identifiers down over the "extra sampling point" rows.
#
# Also defines a small manual alias table for the UKCEH Lakes Tour CHEM
# sheet, whose "Site" field uses slightly different name spellings/spacing
# than the master NAME column (e.g. "Brotherswater" vs "Brothers Water").
# ============================================================================

build_site_lookup <- function(sitemap_path) {

  raw <- read_excel(sitemap_path, sheet = "in")

  names(raw) <- c(
    "NAME", "WBID", "WBID_uklakesportal", "site_UWMN",
    "site_lakestour_chemteox", "lake_lakestour_zoo", "site_microbial",
    "EA_label", "EA_notation",
    "site_algae", "site_chemistry", "site_secchi", "site_tempox"
  )

  # NAME / WBID are only populated on the first row of each lake's block;
  # subsequent rows (extra EA sampling points under the same lake) repeat
  # the same lake but leave NAME/WBID blank. Forward-fill to attach them.
  lookup <- raw %>%
    tidyr::fill(NAME, WBID, WBID_uklakesportal, .direction = "down") %>%
    mutate(WBID = as.integer(WBID))

  lookup
}

# Manual alias table: UKCEH Lakes Tour CHEM sheet "Site" text -> master NAME
# (the LAKE_Lakes_Tour_Chem_TeOx.xlsx mapping column holds a short CODE that
# matches the TeOx sheet, not the full names used in the CHEM sheet, so the
# CHEM sheet has to be matched on NAME instead - these 7 lakes have spelling/
# spacing differences from the master NAME column and need manual aliasing)
lakestour_chem_name_aliases <- c(
  "Brotherswater"          = "Brothers Water",
  "Derwentwater"           = "Derwent Water",
  "Elterwater"             = "Elter Water or Elterwater",
  "Haweswater"             = "Haweswater Reservoir",
  "Wastwater"              = "Wast Water",
  "Windermere North Basin" = "Windermere (N Basin)",
  "Windermere South Basin" = "Windermere (S Basin)"
)


# ============================================================================
# 02_lake_metadata.R
# Builds the static lake-metadata lookup table (one row per WBID) from:
#   - Lake_District_UKCEH_Portal_data_raw.xlsx ("Combined" sheet): morphometry,
#     catchment, land cover, lake typology
#   - Lake_District_UKCEH_Portal_RT_data.csv: volume, discharge, retention time
#
# This is a separate table (LD_lake_metadata.csv), not merged into every row
# of the time-series database. Join on site_id (= WBID).
# ============================================================================

build_lake_metadata <- function(portal_path, rt_path) {

  portal <- read_excel(portal_path, sheet = "Combined") %>%
    mutate(WBID = as.integer(WBID)) %>%
    select(
      WBID, NAME, UKCNTRY, UKCOUNTY, WBLAT, WBLONG, WBALT,
      WBSAREA_km2 = WBSAREA, MNDP_m = MNDP, MXDP_m = MXDP, VOL_m3 = VOL,
      WBPERIM_KM, FETCH_KM, DIST2C_KM,
      UK_CORE_TYPE, UK_ALT_TYPE, UK_SIZE_TYPE, UK_DEPTH_TYPE,
      UK_GEOL_TYPE_, UK_HUMIC_TYPE,
      Alk_uEql, COL, DOC, COND,
      CTAREA_ha = CTAREA,
      HS_CATCHMENT_MEAN_ALTITUDE_M, HS_CATCHMENT_MEAN_SLOPE_DEG
    )

  # Robust column pickup (the DISCHARGE column name varies with m3/y vs m3.y.)
  rt_raw <- read.csv(rt_path)
  disch_col <- grep("^DISCHARGE", names(rt_raw), value = TRUE)[1]
  rt <- rt_raw %>%
    transmute(
      WBID = as.integer(WBID),
      VOL_m3_rt = VOLm3,
      DISCHARGE_m3y = .data[[disch_col]],
      RET_TIME_yrs = RET_TIMEyrs
    )

  lake_metadata <- portal %>%
    left_join(rt, by = "WBID") %>%
    distinct(WBID, .keep_all = TRUE) %>%
    arrange(WBID)

  lake_metadata
}


# ============================================================================
# 03_variable_mapping.R
# Defines the cross-source variable/unit harmonisation lookup table.
#
# This is deliberately an explicit, editable table rather than hidden logic:
# every raw variable name/label from every source is listed once, mapped to
# a standard variable_code + standard_unit + conversion_factor
# (value_standard = value_raw * conversion_factor).
#
# WHERE I WAS CONFIDENT the match is a genuine like-for-like measurement
# (same determinand, compatible unit), a real conversion factor is given.
#
# WHERE I WAS NOT CONFIDENT (different analytical basis, e.g. filtered vs
# unfiltered, or an alkalinity/ion reported in charge-equivalents that can't
# be safely converted to a mass unit without assumptions), the raw
# source-specific label is kept as its own variable_code with factor = 1 and
# a note explaining why it was NOT merged with an existing code. Please
# review rows with a non-empty "review" column.
#
# This table is also written out to output/variable_mapping_lookup.csv so it
# can be checked/edited without touching this script.
# ============================================================================

build_variable_mapping <- function() {

  tribble(
    ~source, ~raw_code, ~variable_code, ~description, ~raw_unit, ~standard_unit, ~conversion_factor, ~review,

    # ---- UKCEH_Monitoring: chemistry.csv (81 codes; only these documented
    #      in metadata.csv, unit assumed as given there - file has no unit column) ----
    "UKCEH_Monitoring_chemistry", "ALKA", "ALKA", "alkalinity",                        "unspecified", "unspecified", 1, "Unit per metadata.csv is ambiguous (labelled ugL but alkalinity is usually ueq/L or mg/L CaCO3) - verify before using quantitatively",
    "UKCEH_Monitoring_chemistry", "NH4N", "NH4N", "ammoniacal nitrogen",               "ugL",         "ug/L",        1, "",
    "UKCEH_Monitoring_chemistry", "PO4P", "PO4P", "soluble reactive phosphorus",       "ugL",         "ug/L",        1, "",
    "UKCEH_Monitoring_chemistry", "TOTP", "TOTP", "total phosphorus",                  "ugL",         "ug/L",        1, "",
    "UKCEH_Monitoring_chemistry", "SIO2", "SIO2", "dissolved reactive silicon",        "ugL",         "ug/L",        1, "",
    "UKCEH_Monitoring_chemistry", "TOCA", "TOCA", "total chlorophyll",                 "ugL",         "ug/L",        1, "",
    "UKCEH_Monitoring_chemistry", "NO3N", "NO3N", "nitrate nitrogen (not in metadata.csv, unit assumed ug/L to match sibling nutrient codes)", "ugL_assumed", "ug/L", 1, "Unit assumed, not in metadata.csv - verify",
    "UKCEH_Monitoring_chemistry", "NO2N", "NO2N", "nitrite nitrogen (not in metadata.csv, unit assumed ug/L)", "ugL_assumed", "ug/L", 1, "Unit assumed, not in metadata.csv - verify",
    "UKCEH_Monitoring_chemistry", "PH",   "PH",   "pH",                                "pH units",    "pH units",    1, "",
    "UKCEH_Monitoring_chemistry", "COND", "COND", "conductivity",                      "unspecified", "unspecified", 1, "",
    # all other chemistry.csv codes not listed here are passed through
    # unchanged (variable_code = raw code, unit = "unspecified") by the
    # loader - see 04_load_sources.R

    # ---- UKCEH_Monitoring: temperature_and_oxygen.csv ----
    "UKCEH_Monitoring_tempox", "TEMP", "TEMP", "temperature",       "degree_C",   "degree_C", 1, "",
    "UKCEH_Monitoring_tempox", "OXYG", "OXYG", "oxygen saturation", "saturation", "percent",  1, "",

    # ---- UKCEH_Monitoring: Secchi.csv ----
    "UKCEH_Monitoring_secchi", "SECC", "SECC", "secchi depth", "metre", "m", 1, "",

    # ---- UKCEH_LakesTour: CHEM sheet (wide -> long; column headers below) ----
    "UKCEH_LakesTour_CHEM", "Total P ug/L",       "TOTP", "total phosphorus",           "ug/L",  "ug/L", 1,    "",
    "UKCEH_LakesTour_CHEM", "PO4-P ug/L",         "PO4P", "soluble reactive phosphorus", "ug/L",  "ug/L", 1,    "",
    "UKCEH_LakesTour_CHEM", "NO3-N ug/L",         "NO3N", "nitrate nitrogen",            "ug/L",  "ug/L", 1,    "",
    "UKCEH_LakesTour_CHEM", "NH4-N ug/L",         "NH4N", "ammoniacal nitrogen",         "ug/L",  "ug/L", 1,    "",
    "UKCEH_LakesTour_CHEM", "SiO2 mg/L",          "SIO2", "dissolved reactive silicon",  "mg/L",  "ug/L", 1000, "",
    "UKCEH_LakesTour_CHEM", "pH",                  "PH",   "pH",                          "pH units", "pH units", 1, "",
    "UKCEH_LakesTour_CHEM", "Chla ugL-1",          "TOCA", "total chlorophyll",           "ug/L",  "ug/L", 1,    "",
    "UKCEH_LakesTour_CHEM", "Secchi depth m",      "SECC", "secchi depth",                "m",     "m",    1,    "",
    "UKCEH_LakesTour_CHEM", "Cond uS/cm",         "COND", "conductivity",                "uS/cm", "unspecified", 1, "Unit differs from other COND sources (uS/cm here vs unspecified elsewhere) - not converted",
    "UKCEH_LakesTour_CHEM", "NO3-N uE/L",         "NO3N_ueqL", "nitrate as charge-equivalents (ionic balance variable, not the nutrient conc.)", "uE/L", "ueq/L", 1, "Kept separate from NO3N (ug/L) - different basis, not the same as the nutrient concentration column",
    "UKCEH_LakesTour_CHEM", "Alky uE/L",          "ALKA_ueqL", "alkalinity (charge-equivalents)", "uE/L", "ueq/L", 1, "Kept separate from chemistry.csv's ALKA - unit basis of that one is unverified (see review note above)",
    "UKCEH_LakesTour_CHEM", "SO4 uE/L",           "SULP_ueqL", "sulphate (charge-equivalents)",   "uE/L", "ueq/L", 1, "",
    "UKCEH_LakesTour_CHEM", "Cl uE/L",            "CHLORIDE_ueqL", "chloride (charge-equivalents)", "uE/L", "ueq/L", 1, "",
    "UKCEH_LakesTour_CHEM", "Ca uE/L",            "CALC_ueqL", "calcium (charge-equivalents)",    "uE/L", "ueq/L", 1, "",
    "UKCEH_LakesTour_CHEM", "Mg uE/L",            "MAGN_ueqL", "magnesium (charge-equivalents)",  "uE/L", "ueq/L", 1, "",
    "UKCEH_LakesTour_CHEM", "Na uE/L",            "SODI_ueqL", "sodium (charge-equivalents)",     "uE/L", "ueq/L", 1, "",
    "UKCEH_LakesTour_CHEM", "K uE/L",             "POTA_ueqL", "potassium (charge-equivalents)",  "uE/L", "ueq/L", 1, "",

    # ---- UKCEH_LakesTour: TeOx sheet ----
    "UKCEH_LakesTour_TeOx", "TEMP",        "TEMP",       "temperature",       "degree_C", "degree_C", 1, "",
    "UKCEH_LakesTour_TeOx", "OXY(%)",      "OXYG",       "oxygen saturation", "percent",  "percent",  1, "",
    "UKCEH_LakesTour_TeOx", "OXY(mg/L)",   "OXYG_mgL",   "oxygen concentration (mass basis)",  "mg/L",   "mg/L",   1, "Not converted to %sat (requires temperature/pressure-dependent solubility) - kept as separate variable",
    "UKCEH_LakesTour_TeOx", "OXY(µmol/L)", "OXYG_umolL", "oxygen concentration (molar basis)", "umol/L", "umol/L", 1, "Not converted - kept as separate variable",

    # ---- UWMN_Scoat_Burnmoor.xlsx ----
    "UKCEH_UWMN", "Value", "TEMP", "temperature (assumed - no variable name given in source file; values/range consistent with lake temperature profiles)", "degree_C_assumed", "degree_C", 1, "Confirm with data provider - source file has no variable/units column",

    # ---- EA_WQ: LDFiltered_*.csv (determinand.label -> code) ----
    "EA_WQ", "Temp Water",    "TEMP", "temperature",                  "cel",  "degree_C", 1,    "",
    "EA_WQ", "O Diss %sat",   "OXYG", "oxygen saturation",            "%",    "percent",  1,    "",
    "EA_WQ", "Oxygen Diss",   "OXYG_mgL", "oxygen concentration (mass basis)", "mg/l", "mg/L", 1, "Not converted to %sat - kept as separate variable, matches UKCEH_LakesTour_TeOx OXYG_mgL",
    "EA_WQ", "SDiscm",        "SECC", "secchi depth",                 "cm",   "m",        0.01, "",
    "EA_WQ", "Chlorophylls",  "TOCA", "total chlorophyll",            "ug/l", "ug/L",     1,    "",
    "EA_WQ", "PhaeophytnAB",  "PHAEO", "phaeophytin",                 "ug/l", "ug/L",     1,    "EA-only variable, no equivalent in other sources",
    "EA_WQ", "Phosphorus-P",  "TOTP", "total phosphorus",             "mg/l", "ug/L",     1000, "",
    "EA_WQ", "Orthophospht",  "PO4P", "soluble reactive phosphorus (unfiltered)",           "mg/l", "ug/L", 1000, "",
    "EA_WQ", "OrthophsFilt",  "PO4P", "soluble reactive phosphorus (filtered)",             "mg/l", "ug/L", 1000, "Merged with unfiltered PO4P - Sven confirmed these are the same determinand for this dataset",
    "EA_WQ", "Ammonia(N)",    "NH4N", "ammoniacal nitrogen (unfiltered)",                    "mg/l", "ug/L", 1000, "",
    "EA_WQ", "NH3 filt N",    "NH4N", "ammoniacal nitrogen (filtered)",                     "mg/l", "ug/L", 1000, "Merged with unfiltered NH4N - Sven confirmed these are the same determinand for this dataset",
    "EA_WQ", "N Oxidised",    "TON",  "total oxidised nitrogen (unfiltered) - matches metadata.csv TON", "mg/l", "ug/L", 1000, "",
    "EA_WQ", "N Oxid Filt",   "TON",  "total oxidised nitrogen (filtered)",                 "mg/l", "ug/L", 1000, "Merged with unfiltered TON - Sven confirmed these are the same determinand for this dataset",
    "EA_WQ", "Nitrate-N",     "NO3N", "nitrate nitrogen (unfiltered)",                      "mg/l", "ug/L", 1000, "",
    "EA_WQ", "Nitrate Filt",  "NO3N", "nitrate nitrogen (filtered)",                        "mg/l", "ug/L", 1000, "Merged with unfiltered NO3N - Sven confirmed these are the same determinand for this dataset",
    "EA_WQ", "Nitrite-N",     "NO2N", "nitrite nitrogen (unfiltered)",                      "mg/l", "ug/L", 1000, "",
    "EA_WQ", "Nitrite Filt",  "NO2N", "nitrite nitrogen (filtered)",                        "mg/l", "ug/L", 1000, "Merged with unfiltered NO2N - Sven confirmed these are the same determinand for this dataset",
    "EA_WQ", "SiO2 Rv",       "SIO2", "dissolved reactive silicon",                         "mg/l", "ug/L", 1000, "",
    "EA_WQ", "pH",            "PH",   "pH",                                                 "phunits", "pH units", 1, "",
    "EA_WQ", "Alky pH 4.5",   "ALKA_mgL", "alkalinity (as CaCO3, mass basis)",              "mg/l", "mg/L", 1, "Kept separate from chemistry.csv's ALKA and CHEM sheet's ALKA_ueqL - different unit basis, not converted",
    "EA_WQ", "Cond @ 20C",    "COND", "conductivity",                                       "us/cm", "unspecified", 1, "Unit differs from other COND sources - not converted",
    "EA_WQ", "N-Kjeldahl",    "TOTN_KJEL", "total Kjeldahl nitrogen", "mg/l", "mg/L", 1, "EA-only variable",
    # All other EA_WQ determinands not listed here (metals, major ions,
    # hardness, colour, turbidity, TOC, coliforms, phenol, oil&grease,
    # solids) are passed through unchanged (variable_code = sanitised raw
    # label, unit = as reported) by the loader.

    # ---- Hydroscape: hydroscape_chemistry_lakes_eidc.csv /
    #      hydroscape_chemistry_uplands_eidc.csv (short 2017-2018 campaign
    #      dataset, wide format, one row per discrete grab sample) ----
    "Hydroscape", "Temperature", "TEMP", "temperature",             "(C)",             "degree_C", 1,    "",
    "Hydroscape", "pH",          "PH",   "pH",                      "pH units",        "pH units", 1,    "",
    "Hydroscape", "EC",          "COND_25C", "conductivity (normalised to 25C)", "uS/cm_at_25C", "unspecified", 1, "Reported at a standardised 25C reference temperature - not merged with other COND sources, which don't state a reference temperature",
    "Hydroscape", "NO2",         "NO2N", "nitrite nitrogen",        "mg_N/L",          "ug/L",     1000, "",
    "Hydroscape", "NO3",         "NO3N", "nitrate nitrogen",        "mg_N/L",          "ug/L",     1000, "",
    "Hydroscape", "NH4",         "NH4N", "ammoniacal nitrogen",     "mg_N/L",          "ug/L",     1000, "",
    "Hydroscape", "SRP",         "PO4P", "soluble reactive phosphorus", "mg_P/L",      "ug/L",     1000, "",
    "Hydroscape", "d15N_NO3",    "ISO_D15N_NO3", "nitrate isotope tracer (delta-15N)", "permil", "permil", 1, "New isotope variable, no equivalent in other sources - not converted/merged with anything",
    "Hydroscape", "d18O_NO3",    "ISO_D18O_NO3", "nitrate isotope tracer (delta-18O)", "permil", "permil", 1, "New isotope variable, no equivalent in other sources - not converted/merged with anything"
  )
}


# ============================================================================
# 04_load_sources.R
# One loader function per data source. Every loader returns rows in the
# SAME common long schema:
#
#   source        - which raw data source the row came from
#   site_id       - master WBID (integer) if matched, else NA
#   site_name     - master lake NAME if matched, else the raw source site code
#   source_site_code - the original site code/name as it appeared in the raw file
#   date          - Date (no time component)
#   datetime      - POSIXct if the source had a time component, else NA
#   variable_code - standardised variable code (see 03_variable_mapping.R)
#   value         - value, converted to standard_unit
#   unit          - standard_unit
#   depth_m       - sample depth in metres (NA if not depth-resolved)
#   depth_zone    - "surface" (<=1m), "bottom" (within 1m of the deepest
#                   sample taken that site-date-variable), "whole_lake_mean"
#                   for depth-resolved sources; NA for single-depth sources
#   censored      - TRUE/FALSE, EA_WQ "less-than" detection-limit results
#   raw_variable  - original variable code/label before mapping
#   raw_value     - original value before unit conversion
#   raw_unit      - original unit before conversion
#   notes         - any review flags carried over from the mapping table
# ============================================================================

empty_common_schema <- function() {
  tibble(
    source = character(), site_id = integer(), site_name = character(),
    source_site_code = character(), date = as.Date(character()),
    datetime = as.POSIXct(character()), variable_code = character(),
    value = double(), unit = character(), depth_m = double(),
    depth_zone = character(), censored = logical(), raw_variable = character(),
    raw_value = double(), raw_unit = character(), notes = character()
  )
}

# ---- helper: attach master site_id/site_name given a source-code column,
#      logging anything that doesn't match (with row counts, date range, and
#      an example lake name if we can guess one from context) ----
attach_site_id <- function(df, code_col, lookup, lookup_col, unmatched_log_name) {
  m <- lookup %>% select(WBID, NAME, !!lookup_col) %>% filter(!is.na(.data[[lookup_col]])) %>% distinct()
  out <- df %>% left_join(m, by = setNames(lookup_col, code_col))
  unmatched_rows <- out %>% filter(is.na(WBID))
  if (nrow(unmatched_rows) > 0) {
    log_unmatched(unmatched_log_name, unmatched_rows, code_col)
  }
  out %>% rename(site_id = WBID, site_name = NAME)
}

unmatched_registry <- new.env()

# code_col may be NULL when called from a context where the "code" is really
# a resolved lake name (e.g. LakesTour CHEM) rather than a raw site code
log_unmatched <- function(log_name, unmatched_rows, code_col = NULL) {
  code_values <- if (!is.null(code_col)) unmatched_rows[[code_col]] else unmatched_rows$site_name_raw
  date_col <- if ("date" %in% names(unmatched_rows)) unmatched_rows$date else as.Date(NA)
  summary_tbl <- tibble(unmatched_code = code_values, date = date_col) %>%
    group_by(unmatched_code) %>%
    summarise(
      n_rows = n(),
      first_date = suppressWarnings(min(date, na.rm = TRUE)),
      last_date = suppressWarnings(max(date, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(log_source = log_name, .before = 1)
  unmatched_registry[[log_name]] <- bind_rows(unmatched_registry[[log_name]], summary_tbl)
}

# call after running ALL loaders to get one combined, de-duplicated table of
# every unmatched site/sub-site code across every source, for manual review
collect_unmatched_log <- function() {
  logs <- mget(ls(unmatched_registry), envir = unmatched_registry)
  if (length(logs) == 0) return(tibble())
  bind_rows(logs) %>%
    group_by(log_source, unmatched_code) %>%
    summarise(
      n_rows = sum(n_rows),
      first_date = suppressWarnings(min(first_date, na.rm = TRUE)),
      last_date = suppressWarnings(max(last_date, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(log_source, desc(n_rows))
}

# ============================================================================
# Hydroscape site matching
# Hydroscape has no dedicated mapping column in the site file - "Sample_ID"
# is free text like "Bassenthwaite epilimnion" or "Angle 1 inflow". This is
# an explicit lookup (built by inspecting every unique Sample_ID in both
# hydroscape_chemistry_*_eidc.csv files) mapping each one to a base lake/tarn
# name (NA if it's a stream/beck/gill sample with no corresponding standing
# water body) and a depth/position descriptor.
#
# IMPORTANT: two names are genuinely ambiguous in the Lake District - there
# are two separate lakes both called "Angle Tarn" (WBID 29093 and 29179) and
# three separate lakes called "Blea Tarn" (WBID 29097, 29218, 29265). Name
# matching alone cannot tell these apart, so resolve_hydroscape_site() below
# disambiguates using the sample's own lat/lon against each candidate's
# coordinates in the lake metadata table, and flags every disambiguated row
# in "notes" so it can be spot-checked.
# ============================================================================
build_hydroscape_sample_mapping <- function() {
  tribble(
    ~sample_id, ~base_lake_name, ~depth_descriptor,
    # -- hydroscape_chemistry_lakes_eidc.csv --
    "Bassenthwaite epilimnion",  "Bassenthwaite Lake", "epilimnion",
    "Bassenthwaite hypolimnion", "Bassenthwaite Lake", "hypolimnion",
    "Bassenthwaite oxycline",    "Bassenthwaite Lake", "oxycline",
    "Black Beck 1", NA, "stream", "Black Beck 2", NA, "stream",
    "Cunsey Beck 1", NA, "stream", "Cunsey Beck 2", NA, "stream", "Cunsey Beck 3", NA, "stream",
    "Derwent epilimnion",  "Derwent Water", "epilimnion",
    "Derwent hypolimnion", "Derwent Water", "hypolimnion",
    "Derwent oxycline",    "Derwent Water", "oxycline",
    "Esthwaite M epilimnion",  "Esthwaite Water", "epilimnion",
    "Esthwaite M hypolimnion", "Esthwaite Water", "hypolimnion",
    "Esthwaite M oxycline",    "Esthwaite Water", "oxycline",
    "Esthwaite Middle",        "Esthwaite Water", NA,
    "Esthwaite N epilimnion",  "Esthwaite Water", "epilimnion",
    "Esthwaite N hypolimnion", "Esthwaite Water", "hypolimnion",
    "Esthwaite N oxycline",    "Esthwaite Water", "oxycline",
    "Esthwaite North",         "Esthwaite Water", NA,
    "Esthwaite S epilimnion",  "Esthwaite Water", "epilimnion",
    "Esthwaite S hypolimnion", "Esthwaite Water", "hypolimnion",
    "Esthwaite S oxycline",    "Esthwaite Water", "oxycline",
    "Esthwaite South",         "Esthwaite Water", NA,
    "Grasmere epilimnion",  "Grasmere", "epilimnion",
    "Grasmere hypolimnion", "Grasmere", "hypolimnion",
    "Grasmere oxycline",    "Grasmere", "oxycline",
    "Grasmere-Rydal 1", NA, "stream", "Grasmere-Rydal 2", NA, "stream",
    "Newlands Beck 1", NA, "stream", "Newlands Beck 2", NA, "stream",
    "River Derwent 1", NA, "stream", "River Derwent 2", NA, "stream",
    "River Derwent 3", NA, "stream", "River Derwent 4", NA, "stream",
    "River Derwent 5", NA, "stream", "River Derwent 6", NA, "stream",
    "Rothay 1", NA, "stream", "Rothay 2", NA, "stream",
    "Rothay 3", NA, "stream", "Rothay 4", NA, "stream",
    "Rydal epilimnion",  "Rydal Water", "epilimnion",
    "Rydal hypolimnion", "Rydal Water", "hypolimnion",
    "Rydal oxycline",    "Rydal Water", "oxycline",
    "Stonethwaite Beck", NA, "stream",
    "Watendlath Beck",   NA, "stream", # distinct from Watendlath Tarn below

    # -- hydroscape_chemistry_uplands_eidc.csv --
    # tarn networks: inflow / in-tarn / outflow all tagged to the same tarn
    # WBID (per instruction - matched like EA sub-points tagged to a lake)
    "Angle 1 inflow", "Angle Tarn", "inflow", "Angle inflow", "Angle Tarn", "inflow",
    "Angle 2 tarn", "Angle Tarn", "tarn", "Angle tarn", "Angle Tarn", "tarn", "Angle Tarn", "Angle Tarn", "tarn",
    "Angle 3 outflow", "Angle Tarn", "outflow", "Angle outflow", "Angle Tarn", "outflow",
    "Blea 1 inflow", "Blea Tarn", "inflow", "Blea inflow", "Blea Tarn", "inflow",
    "Blea 2 tarn", "Blea Tarn", "tarn", "Blea tarn", "Blea Tarn", "tarn", "Blea Tarn", "Blea Tarn", "tarn",
    "Blea 3 outflow", "Blea Tarn", "outflow", "Blea outflow", "Blea Tarn", "outflow",
    "Codale 1 inflow", "Codale Tarn", "inflow", "Codale inflow", "Codale Tarn", "inflow",
    "Codale 2 tarn", "Codale Tarn", "tarn", "Codale open water", "Codale Tarn", "tarn", "Codale Tarn", "Codale Tarn", "tarn",
    "Codale 3 outflow", "Codale Tarn", "outflow", "Codale outflow", "Codale Tarn", "outflow",
    "Easedale 1 inflow", "Easedale Tarn", "inflow", "Easedale inflow", "Easedale Tarn", "inflow",
    "Easedale 2 tarn", "Easedale Tarn", "tarn", "Easedale Tarn", "Easedale Tarn", "tarn",
    "Easedale 3 outflow", "Easedale Tarn", "outflow", "Easedale outflow", "Easedale Tarn", "outflow",
    "Far Easedale Gill 1", NA, "stream", "Far Easedale Gill 2", NA, "stream",
    "Far Easedale Gill 3", NA, "stream", "Far Easedale Gill 4", NA, "stream",
    "Langstrath Beck 1", NA, "stream", "Langstrath Beck 2", NA, "stream",
    "River Derwent 0", NA, "stream",
    "Ruddy Gill 1", NA, "stream", "Ruddy Gill 2", NA, "stream", "Ruddy Gill 3", NA, "stream",
    "Sourmilk Gill", NA, "stream", "Sourmilk Gill 1", NA, "stream", "Sourmilk Gill 2", NA, "stream",
    "Sprinkling (secondary) inflow", "Sprinkling Tarn", "inflow",
    "Sprinkling inflow",  "Sprinkling Tarn", "inflow",
    "Sprinkling Tarn",    "Sprinkling Tarn", "tarn",
    "Sprinkling outflow", "Sprinkling Tarn", "outflow",
    "Stake Beck 1", NA, "stream", "Stake Beck 2", NA, "stream",
    # "Styhead Gill" is NOT explicitly tied to the tarn's own inflow/outflow
    # naming (unlike "Styhead inflow"/"Styhead outflow"/"Styhead Outflow"),
    # so it's kept unmatched rather than assumed to be the same stream
    "Styhead Gill", NA, "stream",
    "Styhead Outflow", "Styhead Tarn", "outflow", "Styhead outflow", "Styhead Tarn", "outflow",
    "Styhead Tarn",    "Styhead Tarn", "tarn",
    "Styhead inflow",  "Styhead Tarn", "inflow",
    "Watendlath 1 tarn", "Watendlath Tarn", "tarn", "Watendlath tarn", "Watendlath Tarn", "tarn", "Watendlath Tarn", "Watendlath Tarn", "tarn",
    "Watendlath 2 outflow", "Watendlath Tarn", "outflow", "Watendlath outflow", "Watendlath Tarn", "outflow",
    "Watendlath inflow", "Watendlath Tarn", "inflow"
  )
}

# disambiguates a base_lake_name against the master site table using the
# sample's own coordinates when more than one WBID shares that name (e.g.
# "Angle Tarn" x2, "Blea Tarn" x3) - picks the nearest by simple planar
# distance in decimal degrees (fine at Lake District scale, candidates here
# are always >5km apart)
resolve_hydroscape_site <- function(base_names, lat, lon, site_lookup, lake_metadata) {
  candidates <- lake_metadata %>% select(WBID, NAME, WBLAT, WBLONG)

  purrr::pmap(list(base_names, lat, lon), function(nm, la, lo) {
    if (is.na(nm)) return(tibble(WBID = NA_integer_, site_name = NA_character_, n_candidates = 0L))
    matches <- candidates %>% filter(NAME == nm)
    if (nrow(matches) == 0) return(tibble(WBID = NA_integer_, site_name = nm, n_candidates = 0L))
    if (nrow(matches) == 1) return(tibble(WBID = matches$WBID, site_name = matches$NAME, n_candidates = 1L))
    # multiple candidates share this name - pick nearest by coordinates
    matches <- matches %>%
      mutate(dist = sqrt((WBLAT - la)^2 + (WBLONG - lo)^2)) %>%
      arrange(dist)
    tibble(WBID = matches$WBID[1], site_name = matches$NAME[1], n_candidates = nrow(matches))
  }) %>% bind_rows()
}

# ============================================================================
# Hydroscape loader - shared by both hydroscape_chemistry_lakes_eidc.csv and
# hydroscape_chemistry_uplands_eidc.csv (identical column layout). Row 1
# (index 0) of the raw file is a units row, not data - skipped on read.
# Longitude is given as positive "degrees West" ("W" column) - negated to
# standard signed decimal degrees to match the lake metadata table.
# ============================================================================
load_hydroscape <- function(path, source_label, site_lookup, lake_metadata, var_mapping, sample_mapping) {
  raw <- read_csv(path, locale = locale(encoding = "latin1"), skip = 2,
                   col_names = c("datetime_chr", "Sample_ID", "N", "W", "Temperature", "pH", "EC",
                                  "NO2", "NO3", "NH4", "SRP", "d15N_NO3", "d18O_NO3"),
                   col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    filter(!is.na(Sample_ID), Sample_ID != "") %>%
    mutate(across(c(N, W, Temperature, pH, EC, NO2, NO3, NH4, SRP, d15N_NO3, d18O_NO3), as.numeric))

  long <- raw %>%
    mutate(
      datetime = dmy_hm(datetime_chr, quiet = TRUE),
      date = as.Date(datetime),
      lon = -W
    ) %>%
    left_join(sample_mapping, by = c("Sample_ID" = "sample_id")) %>%
    pivot_longer(
      cols = c(Temperature, pH, EC, NO2, NO3, NH4, SRP, d15N_NO3, d18O_NO3),
      names_to = "raw_variable", values_to = "raw_value_chr"
    ) %>%
    mutate(raw_value = suppressWarnings(as.numeric(raw_value_chr))) %>%
    filter(!is.na(raw_value))

  site_resolved <- resolve_hydroscape_site(long$base_lake_name, long$N, long$lon, site_lookup, lake_metadata)
  long <- bind_cols(long, site_resolved)

  unmatched_rows <- long %>% filter(is.na(WBID))
  if (nrow(unmatched_rows) > 0) {
    log_unmatched(paste0("hydroscape_", source_label), unmatched_rows, "Sample_ID")
  }

  mapped <- long %>% apply_variable_mapping("raw_variable", "Hydroscape", var_mapping)

  mapped %>%
    transmute(
      source = paste0("Hydroscape_", source_label), site_id = WBID, site_name,
      source_site_code = Sample_ID, date, datetime, variable_code,
      value = raw_value * conversion_factor, unit = standard_unit,
      depth_m = NA_real_,
      depth_zone = case_when(
        depth_descriptor == "epilimnion" ~ "surface",
        depth_descriptor == "hypolimnion" ~ "bottom",
        depth_descriptor == "oxycline" ~ "oxycline",
        TRUE ~ NA_character_
      ),
      censored = FALSE, raw_variable, raw_value, raw_unit = NA_character_,
      notes = if_else(
        !is.na(n_candidates) & n_candidates > 1,
        paste0(notes, if_else(notes == "", "", "; "),
               "site name \"", base_lake_name, "\" matched ", n_candidates,
               " candidate lakes - resolved to WBID ", WBID, " by nearest sample coordinates, verify"),
        notes
      )
    )
}

# ---- helper: map raw variable to standard code/unit using the mapping table,
#      passing through unmapped variables unchanged ----
apply_variable_mapping <- function(df, raw_col, source_name, mapping) {
  map_sub <- mapping %>% filter(source == source_name) %>%
    select(raw_code, variable_code, standard_unit, conversion_factor, review)

  df %>%
    left_join(map_sub, by = setNames("raw_code", raw_col)) %>%
    mutate(
      variable_code = coalesce(variable_code, .data[[raw_col]]),
      conversion_factor = coalesce(conversion_factor, 1),
      standard_unit = coalesce(standard_unit, "unspecified"),
      notes = coalesce(review, "")
    )
}

# ============================================================================
# UKCEH_Monitoring: chemistry.csv
# ============================================================================
load_ukceh_chemistry <- function(path, site_lookup, var_mapping) {
  raw <- read.csv(path, colClasses = c(chemvalu = "numeric")) %>%
    mutate(variable = str_trim(variable))

  df <- raw %>%
    transmute(
      date = as.Date(date),
      source_site_code = site,
      raw_variable = variable,
      raw_value = chemvalu
    ) %>%
    apply_variable_mapping("raw_variable", "UKCEH_Monitoring_chemistry", var_mapping) %>%
    attach_site_id("source_site_code", site_lookup, "site_chemistry", "chemistry_sites") %>%
    transmute(
      source = "UKCEH_Monitoring_chemistry", site_id, site_name, source_site_code,
      date, datetime = as.POSIXct(NA), variable_code,
      value = raw_value * conversion_factor, unit = standard_unit,
      depth_m = NA_real_, depth_zone = NA_character_, censored = FALSE,
      raw_variable, raw_value, raw_unit = NA_character_, notes
    )
  df
}

# ============================================================================
# UKCEH_Monitoring: Secchi.csv
# ============================================================================
load_ukceh_secchi <- function(path, site_lookup, var_mapping) {
  raw <- read.csv(path)
  df <- raw %>%
    transmute(
      date = as.Date(Date),
      source_site_code = Site,
      raw_variable = "SECC",
      raw_value = Diskvalu
    ) %>%
    apply_variable_mapping("raw_variable", "UKCEH_Monitoring_secchi", var_mapping) %>%
    attach_site_id("source_site_code", site_lookup, "site_secchi", "secchi_sites") %>%
    transmute(
      source = "UKCEH_Monitoring_secchi", site_id, site_name, source_site_code,
      date, datetime = as.POSIXct(NA), variable_code,
      value = raw_value * conversion_factor, unit = standard_unit,
      depth_m = NA_real_, depth_zone = NA_character_, censored = FALSE,
      raw_variable, raw_value, raw_unit = NA_character_, notes
    )
  df
}

# ============================================================================
# UKCEH_Monitoring: temperature_and_oxygen.csv (depth-resolved)
# Dates are "DD-Mon-YY" with an AMBIGUOUS 2-digit year (spans 1947-2024) -
# century is inferred manually: yy <= 24 -> 20yy, else 19yy.
# ============================================================================
parse_ukceh_tempox_date <- function(date_str) {
  d <- dmy(date_str, quiet = TRUE)
  yy <- year(d) %% 100
  century <- if_else(yy <= 24, 2000, 1900)
  make_date(century + yy, month(d), day(d))
}

load_ukceh_tempox <- function(path, site_lookup, var_mapping) {
  raw <- read.csv(path, check.names = TRUE) # "depth " -> "depth."
  df <- raw %>%
    transmute(
      date = parse_ukceh_tempox_date(date),
      source_site_code = site,
      raw_variable = str_trim(variable),
      depth_m = depth,
      raw_value = value
    ) %>%
    apply_variable_mapping("raw_variable", "UKCEH_Monitoring_tempox", var_mapping) %>%
    attach_site_id("source_site_code", site_lookup, "site_tempox", "tempox_sites") %>%
    transmute(
      source = "UKCEH_Monitoring_tempox", site_id, site_name, source_site_code,
      date, datetime = as.POSIXct(NA), variable_code,
      value = raw_value * conversion_factor, unit = standard_unit,
      depth_m, depth_zone = NA_character_, censored = FALSE,
      raw_variable, raw_value, raw_unit = NA_character_, notes
    )
  df
}

# ============================================================================
# UKCEH_LakesTour: Lakes_Tour_Chem_TeOx.xlsx, CHEM sheet (wide -> long)
# Site names use full text names, some spelled differently to the master
# NAME column -> resolved via lakestour_chem_name_aliases (01_site_lookup.R)
# ============================================================================
load_lakestour_chem <- function(path, site_lookup, var_mapping, name_aliases) {
  raw <- read_excel(path, sheet = "CHEM", skip = 4)
  # column headers contain embedded \r\n line breaks (e.g. "Total P\r\nug/L") -
  # normalise to a single space so they match the mapping table's raw_code
  names(raw) <- str_squish(gsub("[\r\n]+", " ", names(raw)))
  names(raw)[1:5] <- c("site_name_raw", "date", "year", "month", "season")

  long <- raw %>%
    filter(!is.na(site_name_raw)) %>%
    mutate(date = as.Date(date)) %>%
    pivot_longer(
      cols = -c(site_name_raw, date, year, month, season),
      names_to = "raw_variable", values_to = "raw_value"
    ) %>%
    filter(!is.na(raw_value)) %>%
    mutate(
      site_name_resolved = if_else(
        site_name_raw %in% names(name_aliases),
        unname(name_aliases[site_name_raw]),
        site_name_raw
      )
    )

  mapped <- long %>%
    apply_variable_mapping("raw_variable", "UKCEH_LakesTour_CHEM", var_mapping)

  m <- site_lookup %>% select(WBID, NAME) %>% distinct()
  joined <- mapped %>% left_join(m, by = c("site_name_resolved" = "NAME"))
  unmatched_rows <- joined %>% filter(is.na(WBID))
  if (nrow(unmatched_rows) > 0) {
    log_unmatched("lakestour_chem_sites", unmatched_rows, "site_name_raw")
  }

  joined %>%
    transmute(
      source = "UKCEH_LakesTour_CHEM", site_id = WBID, site_name = site_name_resolved,
      source_site_code = site_name_raw, date, datetime = as.POSIXct(NA),
      variable_code, value = raw_value * conversion_factor, unit = standard_unit,
      depth_m = NA_real_, depth_zone = NA_character_, censored = FALSE,
      raw_variable, raw_value, raw_unit = NA_character_, notes
    )
}

# ============================================================================
# UKCEH_LakesTour: Lakes_Tour_Chem_TeOx.xlsx, TeOx sheet (depth-resolved)
# Site codes match the "site_lakestour_chemteox" mapping column directly.
# ============================================================================
# DATE column is character with MIXED formats: some "M/D/YYYY" strings, some
# Excel serial-date numbers stored as text (e.g. "30717") - handle both.
parse_lakestour_teox_date <- function(date_chr) {
  as_mdy <- suppressWarnings(mdy(date_chr, quiet = TRUE))
  serial <- suppressWarnings(as.numeric(date_chr))
  as_serial <- as.Date(serial, origin = "1899-12-30")
  if_else(!is.na(as_mdy), as_mdy, as_serial)
}

load_lakestour_teox <- function(path, site_lookup, var_mapping) {
  raw <- read_excel(path, sheet = "TeOx", col_types = "text")
  # last column name is "OXY(<mu>mol/L)" - matched by position to dodge
  # encoding issues with the mu character, then given a clean ASCII name
  names(raw)[12] <- "OXY_umolL"
  # source file uses the literal text "NA" for missing values in these
  # numeric columns (confirmed - no other non-numeric values present), which
  # trips R's "NAs introduced by coercion" warning even though the result is
  # correct; suppress it here rather than let it leak to the console.
  raw <- raw %>%
    mutate(across(c(DEPTH, TEMP, `OXY(mg/L)`, `OXY(%)`, OXY_umolL),
                   ~ suppressWarnings(as.numeric(na_if(.x, "NA")))))
  long <- raw %>%
    transmute(
      date = parse_lakestour_teox_date(DATE), source_site_code = LAKE, depth_m = DEPTH,
      TEMP = TEMP, `OXY(mg/L)` = `OXY(mg/L)`, `OXY(%)` = `OXY(%)`,
      `OXY(µmol/L)` = OXY_umolL
    ) %>%
    pivot_longer(cols = c(TEMP, `OXY(mg/L)`, `OXY(%)`, `OXY(µmol/L)`),
                 names_to = "raw_variable", values_to = "raw_value") %>%
    filter(!is.na(raw_value))

  mapped <- long %>% apply_variable_mapping("raw_variable", "UKCEH_LakesTour_TeOx", var_mapping)

  mapped %>%
    attach_site_id("source_site_code", site_lookup, "site_lakestour_chemteox", "lakestour_teox_sites") %>%
    transmute(
      source = "UKCEH_LakesTour_TeOx", site_id, site_name, source_site_code,
      date, datetime = as.POSIXct(NA), variable_code,
      value = raw_value * conversion_factor, unit = standard_unit,
      depth_m, depth_zone = NA_character_, censored = FALSE,
      raw_variable, raw_value, raw_unit = NA_character_, notes
    )
}

# ============================================================================
# UKCEH_UWMN: UWMN_Scoat_Burnmoor.xlsx (depth-resolved, Min/Mean/Max already
# summarised in the source file). Sven confirmed values = temperature.
# Site names ("Burnmoor Tarn","Scoat Tarn") match master NAME directly, but
# we route through site_UWMN mapping column for consistency/robustness.
# ============================================================================
load_uwmn <- function(path, site_lookup, var_mapping) {
  raw <- read_excel(path)
  df <- raw %>%
    transmute(
      date = as.Date(Date), source_site_code = Site, depth_m = Depth_cm / 100,
      metric = Metric, raw_value = Value, raw_variable = "Value"
    ) %>%
    apply_variable_mapping("raw_variable", "UKCEH_UWMN", var_mapping) %>%
    attach_site_id("source_site_code", site_lookup, "site_UWMN", "uwmn_sites") %>%
    transmute(
      source = paste0("UKCEH_UWMN_", metric), site_id, site_name, source_site_code,
      date, datetime = as.POSIXct(NA), variable_code,
      value = raw_value * conversion_factor, unit = standard_unit,
      depth_m, depth_zone = NA_character_, censored = FALSE,
      raw_variable = paste0("Value (", metric, ")"), raw_value,
      raw_unit = NA_character_, notes
    )
  df
}

# ============================================================================
# EA_WQ: LDFiltered_20XX.csv - loops over every file matching the pattern in
# the given folder. Only water-quality determinands are kept (weather codes,
# lab admin numbers, and grid references are dropped, per Sven's decision).
# Censored ("<" qualifier) results: numeric detection-limit value is kept,
# with censored = TRUE flag.
# ============================================================================
EA_NON_WQ_DETERMINANDS <- c(
  "OMR Smpl Num", "Weth-Visibty", "WethPresPrec", "WethPresTemp",
  "Weth7Dy-Temp", "Weth7Dy-Prec", "Weth7Dy-Vsy", "NATGRIDREF"
)

load_ea_wq_folder <- function(folder, site_lookup, var_mapping) {
  files <- list.files(folder, pattern = "^LDFiltered_[0-9]{4}\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No LDFiltered_YYYY.csv files found in ", folder)

  purrr::map_dfr(files, function(f) {
    raw <- read.csv(f, stringsAsFactors = FALSE)
    raw %>%
      filter(!determinand.label %in% EA_NON_WQ_DETERMINANDS) %>%
      transmute(
        datetime = ymd_hms(sample.sampleDateTime, quiet = TRUE),
        date = as.Date(datetime),
        source_site_code = sample.samplingPoint.notation,
        raw_variable = determinand.label,
        raw_value_chr = as.character(result),
        censor_qualifier = resultQualifier.notation,
        censored = resultQualifier.notation %in% c("<", ">"),
        raw_unit_reported = determinand.unit.label
      ) %>%
      mutate(raw_value = suppressWarnings(as.numeric(raw_value_chr))) %>%
      filter(!is.na(raw_value))
  }) %>%
    apply_variable_mapping("raw_variable", "EA_WQ", var_mapping) %>%
    attach_site_id("source_site_code", site_lookup, "EA_notation", "ea_wq_sites") %>%
    transmute(
      source = "EA_WQ", site_id, site_name, source_site_code, date, datetime,
      variable_code, value = raw_value * conversion_factor, unit = standard_unit,
      depth_m = NA_real_, depth_zone = NA_character_, censored,
      raw_variable, raw_value, raw_unit = raw_unit_reported,
      notes = if_else(censored,
                       paste0(notes, if_else(notes == "", "", "; "),
                              "censored result (qualifier: ", censor_qualifier,
                              "), value = reported detection/quantification limit"),
                       notes)
    )
}


# ============================================================================
# 05_depth_aggregation.R
# For depth-resolved sources, collapses each site-date-variable profile into
# three summary rows (per Sven's instruction):
#   - "surface"    : mean of samples with depth_m <= 1
#   - "bottom"     : mean of samples within 1m of the deepest sample taken
#                    that site-date-variable (i.e. bottom 1m of the profile
#                    AS SAMPLED, not of the lake's true max depth - the lake's
#                    true max depth (MXDP) is available in the lake_metadata
#                    table if a stricter definition is later wanted)
#   - "whole_lake" : mean across all sampled depths that site-date-variable
#
# Non-depth-resolved sources (chemistry, Secchi, EA_WQ, LakesTour CHEM) are
# left untouched (depth_zone stays NA, one row per measurement as already
# produced by the loaders).
# ============================================================================

DEPTH_RESOLVED_SOURCES <- c(
  "UKCEH_Monitoring_tempox", "UKCEH_LakesTour_TeOx",
  "UKCEH_UWMN_Maximum", "UKCEH_UWMN_Mean", "UKCEH_UWMN_Minimum"
)

aggregate_depth_profiles <- function(df) {

  depth_rows <- df %>% filter(source %in% DEPTH_RESOLVED_SOURCES)
  other_rows <- df %>% filter(!source %in% DEPTH_RESOLVED_SOURCES)

  grp <- depth_rows %>%
    group_by(source, site_id, site_name, source_site_code, date, datetime,
             variable_code, unit, raw_variable, raw_unit, notes)

  profile_summary <- grp %>%
    summarise(max_depth = max(depth_m, na.rm = TRUE), .groups = "drop")

  depth_with_maxdepth <- depth_rows %>%
    left_join(profile_summary,
              by = c("source", "site_id", "site_name", "source_site_code",
                     "date", "datetime", "variable_code", "unit",
                     "raw_variable", "raw_unit", "notes"))

  make_zone_summary <- function(data, zone_name, filter_expr) {
    data %>%
      filter(!!filter_expr) %>%
      group_by(source, site_id, site_name, source_site_code, date, datetime,
               variable_code, unit, raw_variable, raw_unit, notes) %>%
      summarise(
        value = mean(value, na.rm = TRUE),
        raw_value = mean(raw_value, na.rm = TRUE),
        n_depths_averaged = n(),
        .groups = "drop"
      ) %>%
      mutate(depth_zone = zone_name, depth_m = NA_real_, censored = FALSE)
  }

  surface <- make_zone_summary(depth_with_maxdepth, "surface", quo(depth_m <= 1))
  bottom  <- make_zone_summary(depth_with_maxdepth, "bottom",  quo(depth_m >= max_depth - 1))
  whole   <- make_zone_summary(depth_with_maxdepth, "whole_lake", quo(TRUE))

  depth_summarised <- bind_rows(surface, bottom, whole) %>%
    select(source, site_id, site_name, source_site_code, date, datetime,
           variable_code, value, unit, depth_m, depth_zone, censored,
           raw_variable, raw_value, raw_unit, notes, n_depths_averaged)

  bind_rows(
    other_rows %>% mutate(n_depths_averaged = NA_integer_),
    depth_summarised
  )
}


# ============================================================================
# 06_deduplication.R
# Removes cross-source duplicate observations - the SAME underlying reading
# published in more than one raw database (e.g. one project may have pulled
# in another's readings for a shared site/date). Matching deliberately does
# NOT use "source" as part of the key, per instruction - two rows count as
# the same observation if they agree on:
#   site_id, date, variable_code, unit, depth_zone, censored, and value
# (value compared after rounding to ROUND_DP decimal places, to absorb
# floating-point noise from unit conversion, e.g. mg/L * 1000 -> ug/L)
#
# "source" and "source_site_code" are EXCLUDED from the matching key on
# purpose (that's the whole point - two different databases reporting an
# identical value for the same site/date/variable are exactly the case this
# is meant to catch), but they ARE recorded in the removed-duplicates log so
# you can see which sources were involved in each collision.
#
# Where a group of duplicates is found, ONE row is kept (see keep_source_priority
# below) and the rest are dropped from the combined database and written out
# separately to duplicate_records_removed.csv for review.
# ============================================================================

ROUND_DP <- 6

# deterministic tie-break for which row survives when duplicates are found -
# NOT a claim that one source is more "correct" than another, just a fixed,
# reproducible rule (alphabetical by source name) so re-running the pipeline
# always keeps the same row. Review duplicate_records_removed.csv if a
# different survivor should be picked for a particular case.
deduplicate_records <- function(df) {

  df <- df %>%
    mutate(
      .dedup_row_id = row_number(),
      value_rounded = round(value, ROUND_DP),
      # for MATCHED sites, site_id alone safely identifies the real-world
      # lake regardless of which raw code each source used for it. For
      # UNMATCHED sites (site_id is NA) we have no such confirmed identity -
      # two different unmatched codes (e.g. "CIRA" vs "CIRB") are NOT known
      # to be the same physical location, so folding them into one NA bucket
      # would falsely flag unrelated sub-sites that happen to share a value
      # on the same date as "duplicates". Unmatched rows are therefore only
      # compared against OTHER ROWS WITH THE SAME RAW source_site_code.
      match_key = if_else(!is.na(site_id), paste0("WBID:", site_id),
                           paste0("UNMATCHED:", source_site_code))
    )

  dup_groups <- df %>%
    group_by(match_key, date, variable_code, unit, depth_zone, censored, value_rounded) %>%
    mutate(
      n_in_group = n(), group_id = cur_group_id(),
      n_distinct_sources = n_distinct(source)
    ) %>%
    ungroup() %>%
    # a group only represents a genuine duplicate observation if it spans
    # MORE THAN ONE raw source (that's the actual thing being checked for -
    # the same reading published in more than one database), OR it's a
    # literal repeated row within a single source under the SAME raw site
    # code (an accidental double-entry). A single source reporting the same
    # value on the same date from DIFFERENT site codes (e.g. two distinct EA
    # sampling sub-points on Ullswater both reading suspended solids = 3.0)
    # is not a duplicate - those are deliberately kept as separate sites -
    # so such rows are split back into their own singleton sub-groups.
    mutate(
      dedup_subgroup = if_else(
        n_distinct_sources > 1,
        as.character(group_id),
        paste0(group_id, "::", source, "::", source_site_code)
      )
    ) %>%
    group_by(dedup_subgroup) %>%
    mutate(n_in_group = n(), group_id = cur_group_id()) %>%
    ungroup()

  duplicated_groups <- dup_groups %>% filter(n_in_group > 1)

  if (nrow(duplicated_groups) == 0) {
    return(list(
      deduplicated = df %>% select(-.dedup_row_id, -value_rounded, -match_key),
      removed_log = tibble()
    ))
  }

  # within each duplicate group, keep the first row after sorting by source
  # (alphabetical, deterministic) and log the rest as removed
  duplicated_groups <- duplicated_groups %>%
    group_by(group_id) %>%
    arrange(source, source_site_code, .by_group = TRUE) %>%
    mutate(is_kept = row_number() == 1) %>%
    ungroup()

  kept_from_dups <- duplicated_groups %>% filter(is_kept)
  removed <- duplicated_groups %>% filter(!is_kept)

  # build the review log: one row per REMOVED duplicate, naming which source
  # it was removed in favour of
  removed_log <- removed %>%
    left_join(
      kept_from_dups %>% select(group_id, kept_source = source, kept_source_site_code = source_site_code),
      by = "group_id"
    ) %>%
    transmute(
      site_id, site_name, date, variable_code, value, unit, depth_zone,
      removed_source = source, removed_source_site_code = source_site_code,
      kept_source, kept_source_site_code
    ) %>%
    arrange(site_id, date, variable_code)

  non_duplicated <- dup_groups %>% filter(n_in_group == 1)

  deduplicated <- bind_rows(non_duplicated, kept_from_dups) %>%
    select(-.dedup_row_id, -value_rounded, -match_key, -n_in_group, -group_id,
           -n_distinct_sources, -dedup_subgroup, -is_kept) %>%
    arrange(site_id, date, source, variable_code)

  list(deduplicated = deduplicated, removed_log = removed_log)
}

