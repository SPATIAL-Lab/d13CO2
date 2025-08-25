

######################################################################################
# Extracts data from netCDF climate model output for Phanerozoic from Li et al. (2022) 
# Locates nearest grid cell and linearly interpolates temperature from the two closest 
# timeslices 
######################################################################################

# Load libraries 
library(rgplates)
library(ncdf4)

# Reconstruct paleocoordinates 
comp <- read.csv("Phan/PhanData/assignLatLonSite/PhanCompUpdated.csv")
paleocoord <- array(dim = c(length(comp$d13C), 2))

system.time({
  for (i in seq_along(comp$d13C)) {
    paleocoord[i, ] <- reconstruct(x = c(comp$assigned.lon[i], comp$assigned.lat[i]), age = comp$age[i],
                                   model = "MULLER2022") 
  }
})

comp$paleolon <- as.numeric(paleocoord[, 1])
comp$paleolat <- as.numeric(paleocoord[, 2])
write.csv(comp, file = "Phan/PhanData/PhanCompPaleocoord_MU22.csv")

# Extract temps from CESM netCDF
comp <- read.csv("Phan/PhanData/PhanCompPaleocoord_MU22.csv")
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

write.csv(comp, "Phan/PhanData/PhanCompWithTemp_MU22.csv", row.names = FALSE)


