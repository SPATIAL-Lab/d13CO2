

######################################################################################
# Extracts data from netCDF climate model output for Phanerozoic from Li et al. (2022) 
# Locates nearest grid cell and linearly interpolates temperature from the two closest 
# timeslices 
######################################################################################

# Load libraries 
library(rgplates)
library(ncdf4)


######################################################################################
# Reconstruct paleocoordinates 
######################################################################################

comp <- as.data.frame(read.csv(file = "Phan/PhanData/assignLatLonSite/PhanCompUpdated.csv"))

# Object 'comp' should contain modern lon, lat and observation age 

paleocoord <-array(dim = c(length(comp$d13C),2))

system.time({for (i in 1:length(comp$d13C)){
  paleocoord[i,] <- reconstruct(x = c(comp$assigned.lon[i], comp$assigned.lat[i]), age = comp$age[i],
                                model = "TorsvikCocks2017", anchor = 1) 
}})

comp <- transform(comp, paleolon = as.numeric(paleocoord[,1]))
comp <- transform(comp, paleolat = as.numeric(paleocoord[,2]))

write.csv(comp, file = "Phan/PhanData/PhanCompPaleocoord_TC17.csv")

######################################################################################
# Extract temps from CESM netCDF (Li et al., 2022)
######################################################################################

comp <- as.data.frame(read.csv(file = "Phan/PhanData/PhanCompPaleocoord.csv"))

ncpath <- "Phan/PhanData/"
ncname <- "High_Resolution_Climate_Simulation_Dataset_540_Myr"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "T"  
ncin <- nc_open(ncfname)
print(ncin)

# Get longitude and latitude
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
head(lat)

# Get month
month <- ncvar_get(ncin,"month")
month

# Get simulation
simulation <- ncvar_get(ncin,"simulation")
simulation

# Monthly temperature 
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

# Average temperature across all months (MAT)
MAT_array <- apply(tmp_array, c(1,2,4), FUN=mean)

#Global mean surface temperature from Li et al., 2022 sims 

lat_rad <- lat * pi / 180
lat_weights <- cos(lat_rad)  # area proportional to cos(latitude)
# normalize weights to sum to 1 (accounting for number of longitude points)
lat_weights <- lat_weights / sum(lat_weights)

n_time <- dim(MAT_array)[3]
GMST_Li22_array <- numeric(n_time)

for (t in 1:n_time) {
  temp_t <- MAT_array[, , t]  # [lon, lat]
  # mean over longitude (assuming even spacing, i.e. unweighted)
  zonal_mean <- apply(temp_t, 2, mean, na.rm = TRUE)  # length = lat
  # area-weighted mean over latitude
  GMST_Li22_array[t] <- sum(zonal_mean * lat_weights, na.rm = TRUE)
}

sim.ages <- seq(from = 540, to = 0, by = -10)

#Global mean surface temperature from PhanDA
GMST_PhanDA <- read.csv(file="Phan/PhanData/PhanDA_GMST.csv")

# Initialize vectors 
target.index.lon <- vector(length = length(comp$d13C))
lon.closest <- vector(length = length(comp$d13C))
lon.index <- vector(length = length(comp$d13C))
target.index.lat <- vector(length = length(comp$d13C))
lat.closest <- vector(length = length(comp$d13C))
lat.index <- vector(length = length(comp$d13C))
ceiling.age <- vector(length = length(comp$d13C))
floor.age <- vector(length = length(comp$d13C))
age.index1 <- vector(length = length(comp$d13C))
age.index2 <- vector(length = length(comp$d13C))
MAT <- vector(length = length(comp$d13C))
GMST_Li22 <- vector(length = length(comp$d13C))

# For loop extraction using observations' indices 
for (i in 1:length(comp$d13C)){
  
  # Create target index for lon
  target.index.lon[i] <- which(abs(lon - comp$paleolon[i]) == min(abs(lon - comp$paleolon[i])))
  lon.closest[i] <- lon[target.index.lon[i]]
  lon.index[i] <- as.integer(match(lon.closest[i], lon))
  
  
  # TARGET INDEXING MAY REQURE USING 'NEAREST' FUNCTION TO MAKE SURE NA IS NOT BEING SAMPLED
  
  
  # Create target index for lat
  target.index.lat[i] <- which(abs(lat - comp$paleolat[i]) == min(abs(lat - comp$paleolat[i])))
  lat.closest[i] <- lat[target.index.lat[i]]
  lat.index[i] <- as.integer(match(lat.closest[i], lat))
  
  # Create target index for age
  ceiling.age[i] <- (ceiling(comp$age[i]/10))*10
  floor.age[i] <- (floor(comp$age[i]/10))*10
  
  age.index1[i] <- as.integer(match(ceiling.age[i], sim.ages))
  age.index2[i] <- as.integer(match(floor.age[i], sim.ages))
  
  # Pull surface MAT value using indices derived from lon, lat, and age 
  MAT[i] <- MAT_array[lon.index[i],lat.index[i],age.index1[i]]*((ceiling.age[i] - comp$age[i])/10) +
    MAT_array[lon.index[i],lat.index[i],age.index2[i]]*((comp$age[i] - floor.age[i])/10)
  
  # Pull GMST from CESM for compilation ages 
  GMST_Li22[i] <- GMST_Li22_array[age.index1[i]]*((ceiling.age[i] - comp$age[i])/10) +
    GMST_Li22_array[age.index2[i]]*((comp$age[i] - floor.age[i])/10)
}

# Linearly interpolate PhanDA GMST
PhanDA.interp <- approx(GMST_PhanDA$AverageAge, GMST_PhanDA$GMST_50, xout=comp$age, method="linear") 
PhanDA.interp.hi <- approx(GMST_PhanDA$AverageAge, GMST_PhanDA$GMST_95, xout=comp$age, method="linear") 
PhanDA.interp.lo <- approx(GMST_PhanDA$AverageAge, GMST_PhanDA$GMST_05, xout=comp$age, method="linear") 
PhanDA.mean <- PhanDA.interp[["y"]]
PhanDA.hi <- PhanDA.interp.hi[["y"]]
PhanDA.lo <- PhanDA.interp.lo[["y"]]

temp_offset <- MAT - GMST_Li22
temp_offset_PhanDA <- MAT - PhanDA.mean

comp <- transform(comp, MAT = as.numeric(MAT))
comp <- transform(comp, GMST_Li22 = as.numeric(GMST_Li22))
comp <- transform(comp, GMST_PhanDA = as.numeric(PhanDA.mean))
comp <- transform(comp, GMST_PhanDA_hi = as.numeric(PhanDA.hi))
comp <- transform(comp, GMST_PhanDA_lo = as.numeric(PhanDA.lo))
comp <- transform(comp, temp_offset = as.numeric(temp_offset))
comp <- transform(comp, temp_offset_PhanDA = as.numeric(temp_offset_PhanDA))

write.csv(comp, file = "Phan/PhanData/PhanCompWithTemp_TC17.csv")

# comp_bf <- comp[comp$material == 'bf',]
# comp_bulk <- comp[comp$material == 'bulk',]
# comp_brach <- comp[comp$material == 'brach',]
# comp_bf_clean <- comp_bf[,1:9]




