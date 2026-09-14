# ==============================================================================
# SDM Training Project — Switzerland (pilot region, ~490 species)
# ------------------------------------------------------------------------------
# Purpose: reproducible practice project recreating a multi-species SDM
# workflow (GBIF occurrences -> environmental predictors -> Maxnet -> AUC).
# This is a training/demo pipeline. It is NOT a recreation of the exact SDMs
# underlying any prior published research — it exists to refresh and
# demonstrate the mechanics of the workflow end-to-end.
#
# Requires a free GBIF account. Credentials are read from environment
# variables (set in a local .Renviron file, which is NOT committed to git —
# see .gitignore and the "Credentials" note near the top of Step 3).
# ==============================================================================


# ------------------------------------------------------------------------------
# SETUP: packages and project folders
# ------------------------------------------------------------------------------
# All paths below are relative to the project root, so the script runs the
# same way on any machine as long as it's run from inside the project folder.

required_packages <- c("rgbif", "dplyr", "terra", "geodata", "maxnet", "pROC")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(rgbif)
library(dplyr)
library(terra)
library(geodata)
library(maxnet)
library(pROC)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("data/predictors", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/evaluation", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# STEP 1: build the species list from GBIF
# ------------------------------------------------------------------------------
# What: query GBIF for how many distinct species (with coordinates) have been
# recorded in Switzerland, then keep the best-recorded ~700+ as a starting list.
# Why: record count is a proxy for "will this species have enough data to
# model" — not a measure of true abundance. This deliberately favours
# easy-to-observe species (birds, butterflies, common plants) over
# cryptic ones (fungi, many insects); acceptable for a training project,
# worth naming as a limitation for anything more consequential.

species_facet <- occ_search(
  country = "CH",
  hasCoordinate = TRUE,
  facet = "scientificName",
  facetLimit = 20000, # set well above the expected true count so it isn't capped
  limit = 0            # metadata/facets only — no raw records pulled
)

species_df <- species_facet$facets$scientificName

species_filtered <- species_df %>%
  filter(count >= 30) %>%          # floor: species need enough national records to be worth trying
  arrange(desc(count))

species_top710 <- head(species_filtered, 710)
nrow(species_top710)


# ------------------------------------------------------------------------------
# STEP 2: resolve species names to GBIF taxon keys
# ------------------------------------------------------------------------------
# What: convert scientific names into GBIF's internal numeric taxon keys.
# Why: occ_download() requires taxon keys, not text names, since keys are
# unambiguous (no synonym/spelling issues) in a way plain strings aren't.

species_names <- species_top710$name
keys_lookup <- name_backbone_checklist(species_names)
keys <- keys_lookup$usageKey

table(keys_lookup$matchType) # sanity check: how many EXACT vs FUZZY vs NONE

keys_clean <- keys[!is.na(keys)] # drop any names GBIF couldn't resolve


# ------------------------------------------------------------------------------
# STEP 3: submit the GBIF occurrence download
# ------------------------------------------------------------------------------
# What: request one server-side GBIF download covering all species keys,
# restricted to Switzerland, coordinate-present, human-observation records
# from 2010 onward.
# Why one download instead of a loop of occ_search() calls: occ_search() is
# rate-limited and capped at 300 records per call — completely impractical
# at hundreds of species and millions of records. occ_download() prepares
# one file server-side and gives a citable DOI, which also doubles as a
# reproducibility reference for this exact dataset.
#
# CREDENTIALS: never hardcode a GBIF username/password/email in this script.
# Create a free account at https://www.gbif.org, then add these three lines
# to a local .Renviron file (Session > Set Working Directory > ... in
# RStudio, or `usethis::edit_r_environ()`), which .gitignore excludes from
# version control:
#
#   GBIF_USER=your_username
#   GBIF_PWD=your_password
#   GBIF_EMAIL=your_email@example.com
#
# Restart R after saving .Renviron so the variables load, then Sys.getenv()
# picks them up automatically below.

download_request <- occ_download(
  pred_in("taxonKey", keys_clean),
  pred("country", "CH"),
  pred("hasCoordinate", TRUE),
  pred("basisOfRecord", "HUMAN_OBSERVATION"),
  pred_gte("year", 2010),
  format = "SIMPLE_CSV",
  user  = Sys.getenv("GBIF_USER"),
  pwd   = Sys.getenv("GBIF_PWD"),
  email = Sys.getenv("GBIF_EMAIL")
)

occ_download_meta(download_request) # check status; wait for "SUCCEEDED"
# occ_download_wait(download_request) # or block until it's ready

# Once SUCCEEDED, download and import (only needs to be done once — the file
# is cached locally after this):
d <- occ_download_get(download_request, path = "data/raw/")
occ_data <- occ_download_import(d)

# To resume a finished download in a fresh R session without resubmitting:
#   d <- occ_download_get("<your-download-key>", path = "data/raw/")
#   occ_data <- occ_download_import(d)

nrow(occ_data)
length(unique(occ_data$species))


# ------------------------------------------------------------------------------
# STEP 4: subsample per species
# ------------------------------------------------------------------------------
# What: cap each species at 1,000 occurrence records (keep all if fewer).
# Why: a handful of very common species (deer, foxes, common butterflies)
# can have 500,000+ records each — far more than Maxnet needs, and heavy
# to process. Capping keeps ample signal per species while cutting total
# row count drastically.

set.seed(123) # fixes the random sample so results are reproducible

occ_subsampled <- occ_data %>%
  group_by(species) %>%
  group_modify(~ slice_sample(.x, n = min(1000, nrow(.x)))) %>%
  ungroup()

# Re-check the 30-record floor: capping can't add records, so a species
# that had e.g. 25 records nationally still has 25 here — drop any that
# fall short.
species_counts <- occ_subsampled %>% group_by(species) %>% summarise(n = n())
occ_subsampled <- occ_subsampled %>%
  filter(species %in% species_counts$species[species_counts$n >= 30])


# ------------------------------------------------------------------------------
# STEP 5: clean coordinates
# ------------------------------------------------------------------------------
# What: remove missing coordinates, the (0,0) GBIF artifact, and exact
# duplicate coordinates within each species; then re-check the 30-record
# floor again (deduplication can expose species that looked fine only
# because the same point was logged many times).
# Why: bad coordinates don't cause errors, they silently distort what the
# model learns — worth catching before modelling, not after.

occ_clean <- occ_subsampled %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude),
    !(decimalLatitude == 0 & decimalLongitude == 0)
  ) %>%
  distinct(species, decimalLatitude, decimalLongitude, .keep_all = TRUE)

species_counts_after_dedup <- occ_clean %>% group_by(species) %>% summarise(n = n())
occ_clean <- occ_clean %>%
  filter(species %in% species_counts_after_dedup$species[species_counts_after_dedup$n >= 30])


# ------------------------------------------------------------------------------
# STEP 6: restrict to the study bounding box
# ------------------------------------------------------------------------------
# What: crop to a small pilot region (~50x55 km, Zurich/Aargau area) rather
# than all of Switzerland, then re-check the 30-record floor a third time —
# local counts inside a small box are often much lower than national counts.
# Why small: keeps raster size and per-species compute fast, which was the
# whole point of choosing a coarse resolution and limited extent. The
# trade-off (fewer species clear the local threshold) is accepted here
# rather than enlarging the box or loosening the floor — see project notes.

bbox <- list(xmin = 8.2, xmax = 8.9, ymin = 47.1, ymax = 47.6)

occ_clean <- occ_clean %>%
  filter(
    decimalLongitude >= bbox$xmin, decimalLongitude <= bbox$xmax,
    decimalLatitude  >= bbox$ymin, decimalLatitude  <= bbox$ymax
  )

species_counts_bbox <- occ_clean %>% group_by(species) %>% summarise(n = n())
occ_clean <- occ_clean %>%
  filter(species %in% species_counts_bbox$species[species_counts_bbox$n >= 30])

length(unique(occ_clean$species)) # final species count after all filtering
nrow(occ_clean)


# ------------------------------------------------------------------------------
# STEP 7: download and prepare environmental predictors
# ------------------------------------------------------------------------------
# What: download CHELSA mean annual temperature (bio1) and annual
# precipitation (bio12) at ~1km resolution, crop to the bounding box, add
# elevation (SRTM via geodata) resampled onto the same grid, then check
# pairwise correlation between all three.
# Why 1km / 2-4 variables: matches the "learn fast, don't wait days" goal —
# finer resolution or many more predictors would multiply compute for little
# benefit in a small, low-relief pilot area.

temp_url   <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio1_1981-2010_V.2.1.tif"
precip_url <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio12_1981-2010_V.2.1.tif"

download.file(temp_url,   destfile = "data/predictors/chelsa_temp.tif",   mode = "wb")
download.file(precip_url, destfile = "data/predictors/chelsa_precip.tif", mode = "wb")

study_extent <- ext(bbox$xmin, bbox$xmax, bbox$ymin, bbox$ymax)

temp_cropped   <- crop(rast("data/predictors/chelsa_temp.tif"),   study_extent)
precip_cropped <- crop(rast("data/predictors/chelsa_precip.tif"), study_extent)

elev_raster <- elevation_3s(
  lon = mean(c(bbox$xmin, bbox$xmax)),
  lat = mean(c(bbox$ymin, bbox$ymax)),
  path = "data/predictors/"
)
elev_cropped <- crop(elev_raster, study_extent)
elev_matched <- resample(elev_cropped, temp_cropped, method = "bilinear") # match CHELSA's coarser grid

names(temp_cropped)   <- "temperature"
names(precip_cropped) <- "precipitation"
names(elev_matched)   <- "elevation"

predictors <- c(temp_cropped, precip_cropped, elev_matched)

# Collinearity check — required before finalising the predictor set.
cor(values(predictors), use = "complete.obs")

# Elevation was dropped here after checking: it correlated at r = -0.98 with
# temperature and r = 0.82 with precipitation in this study area (elevation
# drives both via lapse rate and orographic effect in a small mountainous
# box) — keeping all three would make coefficients unstable and hard to
# interpret. Temperature and precipitation were kept as the more directly
# interpretable, less redundant pair.
predictors_final <- predictors[[c("temperature", "precipitation")]]


# ------------------------------------------------------------------------------
# STEP 8: extract predictor values and generate background points
# ------------------------------------------------------------------------------
# What: look up temperature/precipitation at every occurrence point, and
# draw a background sample representing the general environmental
# conditions across the study area (every raster cell, since the box is
# small enough that "all cells" is both simplest and sufficient here).
# Why background, not pseudo-absence: GBIF never confirms a species was
# absent anywhere. Background points represent "what's generally available
# in the landscape," which Maxnet compares against where the species was
# actually seen — a different, weaker, and more defensible assumption than
# treating unseen locations as confirmed absences.

occ_coords <- occ_clean[, c("decimalLongitude", "decimalLatitude")]
occ_env <- extract(predictors_final, occ_coords)

occ_model_data <- occ_clean %>%
  bind_cols(occ_env) %>%
  filter(!is.na(temperature)) # drops the rare point that falls just outside the raster

set.seed(123)
background_points <- spatSample(
  predictors_final,
  size = ncell(predictors_final), # every cell — the box is small enough that this beats a partial random draw
  method = "random",
  na.rm = TRUE,
  xy = TRUE
)


# ------------------------------------------------------------------------------
# STEP 9: fit Maxnet across every species (robust loop)
# ------------------------------------------------------------------------------
# What: for each species, fit a Maxnet model against the shared background,
# predict back onto the training data, and compute AUC. Wrapped in tryCatch
# so one species' failure doesn't stop the run — failures are logged with
# a status and skipped, not silently dropped.
# Why Maxnet: the standard choice for presence-only citizen-science data —
# it doesn't require true absences, only a presence set and a background set.

fit_species_model <- function(species_name, occ_data, background_data) {

  species_data <- occ_data %>% filter(species == species_name)

  if (nrow(species_data) < 30) {
    return(list(species = species_name, status = "SKIPPED_TOO_FEW", model = NULL, auc = NA))
  }

  occ_env_vals <- species_data %>% select(temperature, precipitation)
  bg_env_vals  <- background_data %>% select(temperature, precipitation)

  env_combined     <- as.data.frame(bind_rows(occ_env_vals, bg_env_vals))
  presence_vector  <- c(rep(1, nrow(occ_env_vals)), rep(0, nrow(bg_env_vals)))

  tryCatch({
    model <- maxnet(p = presence_vector, data = env_combined)
    predicted <- predict(model, env_combined, type = "cloglog")
    model_auc <- as.numeric(auc(presence_vector, as.vector(predicted), quiet = TRUE))
    list(species = species_name, status = "SUCCESS", model = model, auc = model_auc)
  }, error = function(e) {
    list(species = species_name, status = paste("ERROR:", conditionMessage(e)), model = NULL, auc = NA)
  })
}

species_list <- unique(occ_model_data$species)
results_list <- vector("list", length(species_list))

for (i in seq_along(species_list)) {
  results_list[[i]] <- fit_species_model(species_list[i], occ_model_data, background_points)
  if (i %% 25 == 0) cat("Processed", i, "of", length(species_list), "species\n")
}

statuses <- sapply(results_list, function(x) x$status)
table(statuses) # confirm SUCCESS count and inspect any errors

aucs <- sapply(results_list, function(x) x$auc)
summary(aucs)
hist(aucs, breaks = 20, main = "AUC distribution across species", xlab = "AUC")


# ------------------------------------------------------------------------------
# STEP 10: save results
# ------------------------------------------------------------------------------
# What: write a human-readable AUC summary table, and save the full model
# list plus the upstream data objects, so the project survives a session
# restart without re-running steps 1-9.

results_summary <- data.frame(
  species = sapply(results_list, function(x) x$species),
  status  = sapply(results_list, function(x) x$status),
  auc     = sapply(results_list, function(x) x$auc)
) %>% arrange(desc(auc))

write.csv(results_summary, "outputs/evaluation/species_auc_summary.csv", row.names = FALSE)

saveRDS(results_list,       "outputs/models/all_species_models.rds")
saveRDS(occ_model_data,     "data/processed/occ_model_data.rds")
saveRDS(background_points,  "data/processed/background_points.rds")
saveRDS(predictors_final,   "data/predictors/predictors_final.rds")


# ------------------------------------------------------------------------------
# STEP 11: generate suitability maps for a spread of species
# ------------------------------------------------------------------------------
# What: predict each chosen species' fitted model across the full raster
# grid (not just at points) to produce an actual suitability surface, for
# a handful of species spanning the AUC range — turns a model object and a
# single number into a visual, comparable result.

species_to_map <- c(
  results_summary$species[1],                       # highest AUC
  results_summary$species[round(nrow(results_summary) * 0.25)],
  results_summary$species[round(nrow(results_summary) / 2)],  # median
  results_summary$species[round(nrow(results_summary) * 0.75)],
  results_summary$species[nrow(results_summary)]     # lowest AUC
)

species_names_all <- sapply(results_list, function(x) x$species)

par(mfrow = c(2, 3))
for (sp in species_to_map) {
  idx <- which(species_names_all == sp)
  model <- results_list[[idx]]$model
  suitability_map <- predict(predictors_final, model, type = "cloglog", na.rm = TRUE)
  plot(suitability_map, main = paste0(sp, "\nAUC = ", round(results_list[[idx]]$auc, 2)))
}
par(mfrow = c(1, 1))

# ==============================================================================
# END OF PIPELINE
#
# What this project deliberately does NOT claim:
# - This is a training/demo project, not a reproduction of any prior
#   published SDM analysis.
# - Species selection is biased toward well-recorded, easy-to-observe taxa
#   (record count reflects observer effort as much as true distribution).
# - Only two climate predictors were used after dropping elevation for
#   collinearity — a fuller analysis would consider land cover, soil, and
#   finer-resolution or additional bioclimatic variables.
# - Background points (not true absences) were used throughout, standard
#   for presence-only Maxnet workflows but a different assumption than
#   presence/absence methods (GLM, RF) would require.
# ==============================================================================
