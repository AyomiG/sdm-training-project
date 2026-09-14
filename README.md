# SDM Training Project — Switzerland Pilot

A reproducible practice project recreating the general workflow of a
multi-species Species Distribution Modelling (SDM) pipeline: GBIF
occurrence data → environmental predictors → Maxnet models → evaluation →
suitability maps.

## What this is — and isn't

This project was built to refresh and demonstrate practical SDM skills
end-to-end, from raw occurrence data through to evaluated, mapped models
across roughly 490 species in a small pilot region of Switzerland.

It is **not** a reproduction of any prior published research or dataset.
Species selection, predictor choice, and study extent were all chosen for
speed and learning value, not scientific novelty — see *Limitations*
below for what that trades away.

## Workflow overview

1. Build a species list from GBIF (Switzerland, coordinate-present records)
2. Resolve species names to GBIF taxon keys
3. Submit one GBIF occurrence download (not per-species queries)
4. Subsample each species to a manageable record cap
5. Clean coordinates (missing values, duplicates, artifacts)
6. Restrict to a small pilot bounding box
7. Download and align environmental predictors (CHELSA climate + SRTM elevation), check collinearity
8. Extract predictor values at occurrences, generate background points
9. Fit a Maxnet model per species, in a loop that logs and skips failures rather than stopping
10. Save all results to disk
11. Generate suitability maps for a spread of species across the AUC range

Each step is documented inline in `sdm_pipeline.R` with a short **What /
Why** comment explaining the reasoning, not just the code.

## Setup

### 1. Install R packages

The script installs anything missing automatically on first run:
`rgbif`, `dplyr`, `terra`, `geodata`, `maxnet`, `pROC`.

### 2. Create a free GBIF account

Register at [gbif.org](https://www.gbif.org) — needed for `occ_download()`,
which this project uses instead of looping `occ_search()` (see the Step 3
comment in the script for why that matters at this scale).

### 3. Set your credentials locally — never in the script

Create a file named `.Renviron` in the project root (same folder as
`sdm_pipeline.R`) containing:

```
GBIF_USER=your_username
GBIF_PWD=your_password
GBIF_EMAIL=your_email@example.com
```

Restart R after saving this file. The script reads these with
`Sys.getenv()` — your credentials never appear in the code and `.Renviron`
is excluded from version control by `.gitignore`.

### 4. Run the script

Run `sdm_pipeline.R` top to bottom from the project root. The GBIF download
step can take anywhere from a few minutes to over an hour depending on
server load — check status with `occ_download_meta()` before proceeding,
or use `occ_download_wait()` to block until it's ready.

If your session restarts after the download finishes but before you've
imported it, you don't need to resubmit the request — pull the finished
file directly by its download key:

```r
d <- occ_download_get("<your-download-key>", path = "data/raw/")
occ_data <- occ_download_import(d)
```

## Project structure

```
sdm-training-project/
├── README.md
├── sdm_pipeline.R
├── .gitignore
├── data/
│   ├── raw/            # GBIF download archive (gitignored)
│   ├── processed/       # cleaned occurrence data, background points (gitignored)
│   └── predictors/      # CHELSA/SRTM rasters (gitignored)
└── outputs/
    ├── models/          # fitted Maxnet models, one per species (gitignored)
    ├── evaluation/      # species_auc_summary.csv (tracked — small, human-readable)
    └── figures/         # suitability maps (gitignored)
```

Raw data, downloaded rasters, and fitted model objects are excluded from
version control since they're large and fully regenerable by re-running
the script. The AUC summary CSV is small and kept tracked as the one
human-readable results artifact.

## Design choices worth knowing

- **Resolution:** ~1km (CHELSA native), not finer — keeps the small pilot
  area computationally trivial without meaningful ecological loss at this
  scale.
- **Study area:** ~50×55 km box (Zurich/Aargau region), not all of
  Switzerland — chosen for speed; scaling the extent later is a one-line
  change to the bounding box, not a rewrite.
- **Predictors:** started with temperature, precipitation, and elevation;
  elevation was dropped after a collinearity check (r = -0.98 with
  temperature in this study area) — see Step 7 in the script.
- **Background points, not pseudo-absences:** GBIF never confirms true
  absence, so background points represent general landscape conditions
  rather than confirmed "not present" locations — the standard, more
  defensible assumption for a presence-only method like Maxnet.
- **Evaluation metric:** AUC only. TSS/sensitivity/specificity are a
  reasonable extension but weren't computed for the full species set in
  this pass.

## Limitations

- Species selection favours well-recorded, easy-to-observe taxa (birds,
  butterflies, common plants) over cryptic ones (fungi, many insects),
  since it's based on GBIF record count, which reflects observer effort
  as much as true distribution.
- Only two climate predictors are used; a fuller analysis would consider
  land cover, soil, and additional bioclimatic variables.
- Results are specific to this small pilot region and are not intended to
  generalise to all of Switzerland without re-running at the full extent.
