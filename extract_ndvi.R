# ==============================================================================
# Extraction du NDVI Sentinel-2 le plus récent pour un polygone
# Packages : rstac, sf, terra
# ==============================================================================

# === CONFIGURATION ===
chemin_shp <- "C:\\Users\\MP_NDIAYE\\Desktop\\webmap\\Pivot.shp"
sortie_ndvi <- "ndvi_recent.tif"

# === CHARGEMENT DES LIBRAIRIES ===
library(rstac)
library(sf)
library(terra)
library(httr)

# === AUTHENTIFICATION COPERNICUS ===
get_auth <- function() {
  req <- POST(
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
    body = list(
      grant_type = "client_credentials",
      client_id = "cdsepublic",
      client_secret = "J2eY0dBqFpUcUupo9x3Q7b4n6tX8qW1m"
    ),
    encode = "form"
  )
  content(req)$access_token
}

# ==============================================================================
# ÉTAPE 1 : Lire le polygone
# ==============================================================================
cat("Lecture du polygone...\n")

polygone <- st_read(chemin_shp, quiet = TRUE)

if (is.null(polygone) || nrow(polygone) == 0) stop("Shapefile vide!")
if (all(st_is_empty(polygone))) stop("Pas de géométrie valide!")
if (is.na(st_crs(polygone))) stop("CRS manquant!")

polygone <- st_transform(polygone, 4326)
cat("Polygone:", nrow(polygone), "feature(s)\n")
print(st_bbox(polygone))

# ==============================================================================
# ÉTAPE 2 : Recherche Sentinel-2
# ==============================================================================
cat("\nConnexion STAC...\n")

stac_c <- stac("https://catalogue.dataspace.copernicus.eu/stac", force_version = "1.0.0")

bbox_poly <- st_bbox(polygone)
bbox_coords <- as.numeric(bbox_poly)
names(bbox_coords) <- NULL

date_debut <- "2020-01-01T00:00:00Z"
date_fin <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")

query <- stac_c |>
  stac_search(
    collections = "sentinel-2-l2a",
    datetime = paste0(date_debut, "/", date_fin),
    bbox = bbox_coords,
    limit = 100
  ) |>
  post_request() |>
  items_fetch()

cat("Images trouvées:", query$total, "\n")
if (query$total == 0) stop("Aucune image!")

# ==============================================================================
# ÉTAPE 3 : Plus récente avec < 50% nuage
# ==============================================================================
cat("\nTri par date...\n")

items_df <- lapply(query$features, function(x) {
  data.frame(
    id = x$id,
    datetime = x$properties$datetime,
    cloud = x$properties$`eo:cloud_cover`
  )
}) |> do.call(rbind, args = _)

items_df$datetime <- as.POSIXct(items_df$datetime)
items_df <- items_df[order(items_df$datetime, decreasing = TRUE), ]
items_df <- items_df[items_df$cloud < 50, ]

if (nrow(items_df) == 0) items_df <- items_df[1, ]

recent_item <- query$features[[1]]
cat("Image:", items_df$id[1], "\n")
cat("Date:", as.character(items_df$datetime[1]), "\n")
cat("Nuage:", items_df$cloud[1], "%\n")

# ==============================================================================
# ÉTAPE 4 : Bandes B04 et B08
# ==============================================================================
cat("\nRecherche bandes B04 et B08...\n")

assets <- recent_item$assets
bande_red <- bande_nir <- NULL

for (name in names(assets)) {
  if (grepl("B04_10m", name)) bande_red <- assets[[name]]$href
  if (grepl("B08_10m", name)) bande_nir <- assets[[name]]$href
}

# Fallback
if (is.null(bande_red)) {
  bande_red <- grep("B04", names(assets), value = TRUE)[1]
  bande_red <- assets[[bande_red]]$href
}
if (is.null(bande_nir)) {
  bande_nir <- grep("B08", names(assets), value = TRUE)[1]
  bande_nir <- assets[[bande_nir]]$href
}

cat("  B04:", basename(bande_red), "\n")
cat("  B08:", basename(bande_nir), "\n")

# ==============================================================================
# ÉTAPE 5 : Lecture VSI et NDVI
# ==============================================================================
cat("\nLecture VSI...\n")

r_red <- terra::rast(paste0("/vsicurl/", bande_red))
r_nir <- terra::rast(paste0("/vsicurl/", bande_nir))

cat("Dimensions:", dim(r_red), "\n")

cat("Calcul NDVI...\n")
ndvi <- (r_nir - r_red) / (r_nir + r_red)

# ==============================================================================
# ÉTAPE 6 : Recadrage et masque
# ==============================================================================
cat("Recadrage...\n")

vect_poly <- terra::vect(polygone)
ndvi_crop <- terra::crop(ndvi, vect_poly)
ndvi_masked <- terra::mask(ndvi_crop, vect_poly)

# ==============================================================================
# ÉTAPE 7 : Sauvegarde
# ==============================================================================
cat("Sauvegarde...\n")

terra::writeRaster(ndvi_masked, sortie_ndvi, 
                   overwrite = TRUE, filetype = "COG")

# ==============================================================================
# RÉSUMÉ
# ==============================================================================
cat("\n=== OK ===\n")
cat("Fichier:", sortie_ndvi, "\n")
cat("NDVI mean:", round(terra::mean(ndvi_masked, na.rm = TRUE), 3), "\n")
cat("NDVI min:", round(terra::min(ndvi_masked, na.rm = TRUE), 3), "\n")
cat("NDVI max:", round(terra::max(ndvi_masked, na.rm = TRUE), 3), "\n")
