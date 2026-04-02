# packages
library(sf)
library(rgplates)
suppressMessages(sf::sf_use_s2(FALSE))    # avoid antimeridian surprises

## ------------------ USER PATHS (yours) ------------------
rot_file  <- "TectonicModels/Muller2022/1000_0_rotfile.rot"
static_fc <- "TectonicModels/Muller2022/shapes_static_polygons_Merdith_etal.gpml"
topo_dir  <- "TectonicModels/Muller2022/Topologies"
topo_fc   <- list.files(topo_dir, pattern = "\\.gpmlz?$", full.names = TRUE)

stopifnot(file.exists(rot_file),
          file.exists(static_fc),
          length(topo_fc) > 0)

## ------------------ INPUT TABLE -------------------------
comp <- read.csv("Phan/PhanData/assignLatLonSite/PhanCompUpdated_WithComponent.csv")

# Ensure ages are in Ma
if (max(comp$age, na.rm = TRUE) > 5e3) comp$age <- comp$age/1e3  # kyr -> Ma
if (max(comp$age, na.rm = TRUE) > 1e4) comp$age <- comp$age/1e6  # yr  -> Ma

# Keep valid rows and tag order
ok <- is.finite(comp$assigned.lon) & comp$assigned.lon >= -180 & comp$assigned.lon <= 180 &
  is.finite(comp$assigned.lat) & comp$assigned.lat >=  -90 & comp$assigned.lat <=   90 &
  is.finite(comp$age)          & comp$age >= 0 & comp$age <= 750
comp        <- comp[ok, , drop = FALSE]
comp$row_id <- seq_len(nrow(comp))

pts0 <- st_as_sf(comp, coords = c("assigned.lon","assigned.lat"),
                 crs = 4326, remove = FALSE)

## ------------------ DETECT rgplates API -----------------
rc <- get("reconstruct", asNamespace("rgplates"))
argn <- names(formals(rc))

api <- "unknown"
if (all(c("method","rotations","static_polygons","topology_features") %in% argn)) {
  api <- "offline_new"     # reconstruct(x, age, method='offline', rotations=..., static_polygons=..., topology_features=..., anchor_plate_id=...)
} else if (all(c("method","rotation","static_polygons","topologies") %in% argn)) {
  api <- "offline_old"     # reconstruct(x, age, method='offline', rotation=..., static_polygons=..., topologies=..., anchor_plate_id=...)
} else if (all(c("x","age","model") %in% argn)) {
  api <- "online_only"     # reconstruct(x, age, model='PALEOMAP') — no offline args supported
}

message("Detected rgplates::reconstruct() API: ", api)

if (api == "online_only" || api == "unknown") {
  stop(
    "\nYour installed 'rgplates' does not expose an offline reconstruction interface.\n",
    "It only supports 'model=\"PALEOMAP\"' via the web service.\n\n",
    "To run offline with local rotation/static/topology files, install a build of 'rgplates' that supports:\n",
    "  reconstruct(x, age, method='offline', rotations/rotation=..., static_polygons=..., topology_features/topologies=..., ...)\n\n",
    "Quick checks you can run:\n",
    "  - packageVersion('rgplates')\n",
    "  - args(rgplates::reconstruct)\n\n",
    "After upgrading, re-run this script."
  )
}

## ------------------ OFFLINE RECONSTRUCTION -------------------
# Round ages so near-equals share a call (tune digits as needed)
comp$age_key <- round(comp$age, 1)   # try 1.0 if you want fewer unique ages
uages <- sort(unique(comp$age_key))

anchor_id <- 0   # mantle frame (adjust if your model needs a different anchor)

reconstruct_offline <- switch(
  api,
  "offline_new" = function(x, a) {
    rgplates::reconstruct(
      x                 = x,
      age               = a,
      method            = "offline",
      rotations         = rot_file,
      static_polygons   = static_fc,
      topology_features = topo_fc,
      anchor_plate_id   = anchor_id
    )
  },
  "offline_old" = function(x, a) {
    rgplates::reconstruct(
      x                = x,
      age              = a,
      method           = "offline",
      rotation         = rot_file,
      static_polygons  = static_fc,
      topologies       = topo_fc,
      anchor_plate_id  = anchor_id
    )
  }
)

out_list <- vector("list", length(uages))
names(out_list) <- as.character(uages)

cat("Unique rounded ages: ", length(uages), "\n", sep = "")
for (j in seq_along(uages)) {
  a    <- uages[j]
  sel  <- comp$age_key == a
  ptsA <- pts0[sel, , drop = FALSE]
  
  recA <- reconstruct_offline(ptsA, a)
  
  xy <- st_coordinates(st_geometry(recA))
  out_list[[j]] <- data.frame(
    row_id    = recA$row_id,   # carried through from input sf
    paleo_lon = xy[, 1],
    paleo_lat = xy[, 2]
  )
  
  if (j %% 10 == 0) cat("\r", j, "/", length(uages), " ages", sep = "")
}
cat("\nDone.\n")

paleo_df <- do.call(rbind, out_list)

## ------------------ JOIN BACK & SAVE -------------------------
comp_out <- merge(comp, paleo_df, by = "row_id", all.x = TRUE, sort = FALSE)
comp_out <- comp_out[order(comp_out$row_id), ]

out_file <- "Phan/PhanData/PhanCompPaleocoord_local.csv"
write.csv(comp_out, out_file, row.names = FALSE)
cat("Wrote: ", out_file, "\n", sep = "")


####################################################################

# Extract temps from CESM netCDF

comp <- read.csv("Phan/PhanData/PhanCompPaleocoord_local.csv")
ncin <- nc_open("Phan/PhanData/High_Resolution_Climate_Simulation_Dataset_540_Myr.nc")

lon <- ncvar_get(ncin, "lon")
lat <- ncvar_get(ncin, "lat")
month <- ncvar_get(ncin, "month")
simulation <- ncvar_get(ncin, "simulation")
tmp_array <- ncvar_get(ncin, "T")
MAT_array <- apply(tmp_array, c(1, 2, 4), mean)

# Calculate global mean surface temperature (GMST)
lat_rad <- lat * pi / 180
lat_weights <- cos(lat_rad) / sum(cos(lat_rad))  # Normalize weights

n_time <- dim(MAT_array)[3]
GMST_Li22_array <- numeric(n_time)

for (t in 1:n_time) {
  zonal_mean <- apply(MAT_array[, , t], 2, mean, na.rm = TRUE)
  GMST_Li22_array[t] <- sum(zonal_mean * lat_weights, na.rm = TRUE)
}

sim.ages <- seq(from = 540, to = 0, by = -10)
GMST_PhanDA <- read.csv("Phan/PhanData/PhanDA_GMST.csv")

# Initialize vectors
n <- nrow(comp)
MAT <- GMST_Li22 <- numeric(n)
lon.index <- lat.index <- age.index1 <- age.index2 <- rep(NA_integer_, n)
ceiling.age <- floor.age <- lon.closest <- lat.closest <- rep(NA_real_, n)

for (i in seq_len(n)) {
  # Longitude index
  if (!is.na(comp$paleolon[i])) {
    lon_diff <- abs(lon - comp$paleolon[i])
    if (!all(is.na(lon_diff))) {
      lon.index[i] <- which.min(lon_diff)
      lon.closest[i] <- lon[lon.index[i]]
    }
  }
  
  # Latitude index
  if (!is.na(comp$paleolat[i])) {
    lat_diff <- abs(lat - comp$paleolat[i])
    if (!all(is.na(lat_diff))) {
      lat.index[i] <- which.min(lat_diff)
      lat.closest[i] <- lat[lat.index[i]]
    }
  }
  
  # Age indices
  ceiling.age[i] <- ceiling(comp$age[i] / 10) * 10
  floor.age[i] <- floor(comp$age[i] / 10) * 10
  age.index1[i] <- match(ceiling.age[i], sim.ages)
  age.index2[i] <- match(floor.age[i], sim.ages)
  
  # Pull MAT and GMST if all indices are available
  if (!any(is.na(c(lon.index[i], lat.index[i], age.index1[i], age.index2[i])))) {
    MAT[i] <- MAT_array[lon.index[i], lat.index[i], age.index1[i]] * ((ceiling.age[i] - comp$age[i]) / 10) +
      MAT_array[lon.index[i], lat.index[i], age.index2[i]] * ((comp$age[i] - floor.age[i]) / 10)
    
    GMST_Li22[i] <- GMST_Li22_array[age.index1[i]] * ((ceiling.age[i] - comp$age[i]) / 10) +
      GMST_Li22_array[age.index2[i]] * ((comp$age[i] - floor.age[i]) / 10)
  } else {
    MAT[i] <- NA
    GMST_Li22[i] <- NA
  }
}

# Interpolate PhanDA GMST
PhanDA.interp <- approx(GMST_PhanDA$AverageAge, GMST_PhanDA$GMST_50, xout = comp$age)$y
PhanDA.hi <- approx(GMST_PhanDA$AverageAge, GMST_PhanDA$GMST_95, xout = comp$age)$y
PhanDA.lo <- approx(GMST_PhanDA$AverageAge, GMST_PhanDA$GMST_05, xout = comp$age)$y

# Final calculations
temp_offset <- MAT - GMST_Li22
temp_offset_PhanDA <- MAT - PhanDA.interp

# Add to dataframe
comp$MAT <- MAT
comp$GMST_Li22 <- GMST_Li22
comp$GMST_PhanDA <- PhanDA.interp
comp$GMST_PhanDA_hi <- PhanDA.hi
comp$GMST_PhanDA_lo <- PhanDA.lo
comp$temp_offset <- temp_offset
comp$temp_offset_PhanDA <- temp_offset_PhanDA

write.csv(comp, "Phan/PhanData/PhanCompWithTemp_component.csv", row.names = FALSE)




