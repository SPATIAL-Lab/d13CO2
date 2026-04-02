

# ----------- write present-day points to a shapefile for GPlates ----------- #
library(sf)

# read in data
comp <- read.csv("Phan/PhanData/assignLatLonSite/PhanCompUpdated_WithComponent.csv")

# Ensure ages in Ma
if (max(comp$age, na.rm = TRUE) > 5e3) comp$age <- comp$age/1e3  # kyr->Ma
if (max(comp$age, na.rm = TRUE) > 1e4) comp$age <- comp$age/1e6  # yr->Ma

# Keep valid rows
ok <- is.finite(comp$assigned.lon) & comp$assigned.lon >= -180 & comp$assigned.lon <= 180 &
  is.finite(comp$assigned.lat) & comp$assigned.lat >=  -90 & comp$assigned.lat <=   90 &
  is.finite(comp$age)          & comp$age >= 0 & comp$age <= 750
comp  <- comp[ok, , drop = FALSE]

# Give every row a stable string id GPlates will carry through as the feature "Name"
comp$id   <- sprintf("row_%07d", seq_len(nrow(comp)))
comp$Name <- comp$id  # GPlates preserves the 'Name' attribute in exports

pts_sf <- st_as_sf(comp,
                   coords = c("assigned.lon", "assigned.lat"),
                   crs = 4326, remove = FALSE)

# Write as ESRI Shapefile (GeoPackage can crash GPlates on big layers)
dir.create("gplates_in", showWarnings = FALSE)
st_write(pts_sf, "gplates_in/points.shp", delete_layer = TRUE, quiet = TRUE)
cat("Wrote gplates_in/points.shp (", nrow(pts_sf), " features)\n", sep = "")



# ----------- Gplates desktop instructions ----------- #

# 1) Load points as a feature collection
# 
# Use File → Open Feature Collection… and open points.geojson/.shp 
# 
# points should appear in the Layers panel (right side). If not: View → Show Layers.
# 
# 2) Make sure the points are “reconstructable”
# 
# In Layers, right-click points layer → Change Layer Type… → choose Reconstructable (Regular Feature Collection).
# 
# If the points don’t follow plates when you scrub time, assign IDs: Tools → Assign plate IDs… (Target =  points; Sources = model’s static polygons/topologies; Time = 0 Ma).
# 
# 3) Enable ONLY points for export
# 
# In the Layers panel, un-check every other layer (coastlines, polygons, topologies, etc.).
# 
# Check only points layer. Export pulls from the set of enabled reconstructable layers.
# 
# 4) Export
# 
# Reconstruction → Export
# 
# Choose Export Time Sequence of Snapshots and set from / to / increment.
# 
# Add Export → Reconstructed Geometries → choose Shapefiles (*.shp) (or GeoJSON). Template: reconstructed_%0.2fMa.
# 
# Back in the parent dialog, click Edit to set file options but layer selection is already handled by step 3.
# 
# Pick a Target directory and Begin Export.
# 
# 5) What you’ll get & how to read it
# 
# With “Export to multiple files” and “Separate output directory per input file/layer” checked, you’ll get a subfolder per input layer and a shapefile for each time step.
# 
# If you ever see files named *_polygon.shp, it means a polygon layer was still enabled—go back and disable it and re-export.


## ---- libs ----
if (!requireNamespace("sf", quietly = TRUE)) stop("Install 'sf'")
library(sf)

## ---- inputs ----
csv_in   <- "Phan/PhanData/assignLatLonSite/PhanCompUpdated_WithComponent.csv"
shp_dir  <- "/Users/dustintharper/Documents/Gplates/gplates_out/points"
csv_out  <- "Phan/PhanData/PhanCompPaleocoord_from_export.csv"

## ---- read CSV and normalize columns ----
comp <- read.csv(csv_in, check.names = FALSE)  # keep spaces in names
stopifnot(all(c("age","assigned lon","assigned lat") %in% names(comp)))

## ages must be Ma
if (max(comp$age, na.rm = TRUE) > 1e4) comp$age <- comp$age/1e6
if (max(comp$age, na.rm = TRUE) >  600) comp$age <- comp$age/1e3

## keep valid rows
ok <- is.finite(comp$age) &
  is.finite(comp[["assigned lon"]]) & comp[["assigned lon"]] >= -180 & comp[["assigned lon"]] <= 180 &
  is.finite(comp[["assigned lat"]]) & comp[["assigned lat"]] >=  -90 & comp[["assigned lat"]] <=  90
comp <- comp[ok, , drop = FALSE]
comp$row_id <- seq_len(nrow(comp))                       # preserve order

## ---- index the exported shapefiles by age ----
shps <- list.files(shp_dir, pattern = "^reconstructed_.*Ma\\.shp$", full.names = TRUE)
if (!length(shps)) stop("No reconstructed_*.shp files found in: ", shp_dir)

## pull ages from filenames like "reconstructed_1.00Ma.shp"
get_age <- function(p) {
  # returns numeric age in Ma (e.g., 1.00)
  x <- sub("^.*reconstructed_([0-9.]+)Ma\\.shp$", "\\1", basename(p))
  suppressWarnings(as.numeric(x))
}
ages_available <- vapply(shps, get_age, numeric(1))
names(shps) <- format(ages_available, nsmall = 2, trim = TRUE)

## function to snap an age to the nearest available file (within tol)
snap_age <- function(a, tol = 0.51) {  # ±0.51 My tolerance around 1-Myr steps
  i <- which.min(abs(ages_available - a))
  if (!length(i) || is.infinite(ages_available[i])) return(NA_real_)
  if (abs(ages_available[i] - a) <= tol) ages_available[i] else NA_real_
}

## compute snapped ages per row
comp$age_snap <- vapply(comp$age, snap_age, numeric(1))
miss_age <- sum(!is.finite(comp$age_snap))
if (miss_age) message("Rows with no nearby export time step: ", miss_age)

## ---- helper: rounded join key on present-day lon/lat ----
key_ll <- function(lon, lat, digits = 5L)
  paste(round(lon, digits), round(lat, digits), sep = "|")

comp$key_ll <- key_ll(comp[["assigned lon"]], comp[["assigned lat"]])

## ---- loop over unique snapped ages and join ----
uages <- sort(unique(comp$age_snap[is.finite(comp$age_snap)]))
out_list <- vector("list", length(uages))

for (j in seq_along(uages)) {
  a  <- uages[j]
  shp <- shps[format(a, nsmall = 2, trim = TRUE)]
  if (!length(shp) || !file.exists(shp)) next
  
  # read this time step
  g  <- st_read(shp, quiet = TRUE)
  xy <- st_coordinates(g)
  
  # build a small lookup by present-day lon/lat (fields from export!)
  if (!all(c("assgnd_ln","assgnd_lt") %in% names(g))) {
    stop("Shapefile is missing 'assgnd_ln'/'assgnd_lt' fields: ", shp)
  }
  look <- data.frame(
    key_ll    = key_ll(g$assgnd_ln, g$assgnd_lt),
    paleo_lon = xy[,1],
    paleo_lat = xy[,2],
    stringsAsFactors = FALSE
  )
  
  # merge back to the rows needing age 'a'
  sel <- comp$age_snap == a
  m   <- merge(comp[sel, c("row_id","key_ll")], look, by = "key_ll", all.x = TRUE, sort = FALSE)
  out_list[[j]] <- m[, c("row_id","paleo_lon","paleo_lat")]
  if (j %% 25 == 0) cat("\r", j, "/", length(uages), " ages", sep = "")
}
cat("\n")

## ---- combine and reattach to comp in original order ----
found <- do.call(rbind, out_list)
comp$paleo_lon <- NA_real_; comp$paleo_lat <- NA_real_
if (!is.null(found) && nrow(found)) {
  comp$paleo_lon[match(found$row_id, comp$row_id)] <- found$paleo_lon
  comp$paleo_lat[match(found$row_id, comp$row_id)] <- found$paleo_lat
}

## ---- report & save ----
n_ok   <- sum(is.finite(comp$paleo_lon) & is.finite(comp$paleo_lat))
n_all  <- nrow(comp)
message("Assigned paleocoords for ", n_ok, " / ", n_all, " rows (",
        sprintf("%.1f%%", 100*n_ok/n_all), ").")

write.csv(comp, csv_out, row.names = FALSE)
message("Wrote: ", csv_out)


