


######## Load data matrix #########

prox.in <- as.data.frame(read.csv(file = "data/processed/PhanCompWithTemp_PALEOMAP.csv"))
prox.in <- cbind(prox.in[,1:7], prox.in[,9:10], prox.in[,21:27],rep(x = 2, times = nrow(prox.in)))
names(prox.in) <- c("age", "d13C", "source", "site", "lat", "lon", "category", 
                    "paleolon","paleolat", "MAT", "GMST_Li22", "GMST_PhanDA", "GMST_PhanDA_hi",
                    "GMST_PhanDA_lo", "temp_offset", "temp_offset_PhanDA", "temp_offset_sd")



######## Culling names #########

# Vector of prefixes you want to keep exactly as-is (add more as needed)
keep_prefixes <- c("165-999", "108-659", "138-851","90-590","184-1146", "7-62", "130-806", "39-354", "207-1258", "DSDP384", "207-1257", "207-1259","171-1052", "DSDP577", "207-1260","Ashford Black Marble Quarry", "Ghent Stream", "Isle of Skye", "Peniche", "Quarry west of Hwy. A-68, NY", "Waterswallow Quarry SK")  # <- add more here

# Safely escape regex metacharacters in the prefixes (e.g., ".", "(", etc.)
escape_regex <- function(x) gsub("([][{}()^$.|*+?\\\\])", "\\\\\\1", x)

# Build a single pattern that matches any of the prefixes at the start of the string
pat <- paste0("^(", paste(escape_regex(keep_prefixes), collapse = "|"), ").*")

# Apply culling: replace whole match with just the captured prefix
all_sites_culled <- sub(pat, "\\1", all_sites)





######## Replacement names #########

replace_string <- function(x, old, new) {
  x[x == old] <- new
  x
}

all_sites_culled <- replace_string(all_sites_culled, old = " Basleo", new = "Basleo")
all_sites_culled <- replace_string(all_sites_culled, old = " Broidi", new = "Broidi")
all_sites_culled <- replace_string(all_sites_culled, old = " Diverses", new = "Diverses")
all_sites_culled <- replace_string(all_sites_culled, old = " Fatu Tuanini", new = "Fatu Tuanini")
all_sites_culled <- replace_string(all_sites_culled, old = " Kachchh Basin", new = "Kachchh Basin")
all_sites_culled <- replace_string(all_sites_culled, old = " Kioeuoko", new = "Kioeuoko")
all_sites_culled <- replace_string(all_sites_culled, old = " Korte 12", new = "Korte 12")
all_sites_culled <- replace_string(all_sites_culled, old = " NE SanVigilio Veneto", new = "NE SanVigilio Veneto")
all_sites_culled <- replace_string(all_sites_culled, old = " Niki-Niki", new = "Niki-Niki")
all_sites_culled <- replace_string(all_sites_culled, old = " Pauseokah", new = "Pauseokah")
all_sites_culled <- replace_string(all_sites_culled, old = " Tubulopo", new = "Tubulopo")
all_sites_culled <- replace_string(all_sites_culled, old = "(Kutor); DorfPetschistschi (Wolgal near Kazan)", new = "Kutor")
all_sites_culled <- replace_string(all_sites_culled, old = "(Kutor); Gauri 1842", new = "Kutor")
all_sites_culled <- replace_string(all_sites_culled, old = "(Kutor); Gauri 1843", new = "Kutor")
all_sites_culled <- replace_string(all_sites_culled, old = "(Kutor); Gauri 1844", new = "Kutor")
all_sites_culled <- replace_string(all_sites_culled, old = "Altai Salair, Russia", new = "Altai-Salair")
all_sites_culled <- replace_string(all_sites_culled, old = "Amti-Atlas", new = "Anti-Atlas")
all_sites_culled <- replace_string(all_sites_culled, old = "Anti Atlas Mts, Morocco", new = "Anti-Atlas")
all_sites_culled <- replace_string(all_sites_culled, old = "Anticosti Island, Canada", new = "Anticosti Island")
all_sites_culled <- replace_string(all_sites_culled, old = "Anticosti Island, Quebec", new = "Anticosti Island")
all_sites_culled <- replace_string(all_sites_culled, old = "1146", new = "184-1146")
all_sites_culled <- replace_string(all_sites_culled, old = "1258", new = "207-1258")
all_sites_culled <- replace_string(all_sites_culled, old = "Beard 1", new = "Beard")
all_sites_culled <- replace_string(all_sites_culled, old = "Beard 2", new = "Beard")
all_sites_culled <- replace_string(all_sites_culled, old = "Blake Nose (DSDP Site 390A", new = "DSDP 390")
all_sites_culled <- replace_string(all_sites_culled, old = "Cantabrian Mts.", new = "Cantabrian Mts, Spain")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP site 288A, Pacific", new = "288")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP site 289, Pacific", new = "289")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP site 305, Pacific", new = "305")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP Site 354.", new = "39-354")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP site 390, Atlantic", new = "390")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP Site 90-590", new = "90-590")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP219", new = "219")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP363", new = "363")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP366A", new = "366")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP384", new = "384")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP390.", new = "390")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP392", new = "392")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP463.", new = "463")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP577", new = "577")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP77B", new = "77")
all_sites_culled <- replace_string(all_sites_culled, old = "Jhra Dome", new = "Jara Dome")
all_sites_culled <- replace_string(all_sites_culled, old = "Meishan", new = "Meishan, China")
all_sites_culled <- replace_string(all_sites_culled, old = "MidPacificDSDP463.", new = "463")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP site 1049C (western North Atlantic)", new = "1049")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP Site 1050", new = "171-1050")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP Site 1259", new = "207-1259")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP Site 366.", new = "366")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP Site 366A.", new = "366")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP Site 540.", new = "540")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP-121/758.", new = "121-758")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP-154/925B.", new = "154-925")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP-198-2010B", new = "198-1210")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP1052", new = "1052")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP1258.", new = "207-1258")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP202-1241.", new = "202-1241")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP758A", new = "121-758")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP865", new = "865")
all_sites_culled <- replace_string(all_sites_culled, old = "ODP926", new = "926")
all_sites_culled <- replace_string(all_sites_culled, old = "Once-a-week Quarry, 1 mile N of Sheldon SK 158 681", new = "Once-a-week Quarry, Sheldon SK")
all_sites_culled <- replace_string(all_sites_culled, old = "Once-a-week Quarry, 1 mile N of Sheldon SK 158 682", new = "Once-a-week Quarry, Sheldon SK")
all_sites_culled <- replace_string(all_sites_culled, old = "Pacific, ODP807.", new = "807")
all_sites_culled <- replace_string(all_sites_culled, old = "Province Palencia; outcrop along Arroyo de Castilleria River approximately 3 km east of Celada de Robeldo", new = "Province Palencia; 3 km east of Celada de Robeldo")
all_sites_culled <- replace_string(all_sites_culled, old = "Site 1049", new = "1049")
all_sites_culled <- replace_string(all_sites_culled, old = "Site 1050", new = "1050")
all_sites_culled <- replace_string(all_sites_culled, old = "Site 305", new = "305")
all_sites_culled <- replace_string(all_sites_culled, old = "Site 463", new = "463")
all_sites_culled <- replace_string(all_sites_culled, old = "Site 806B", new = "806")
all_sites_culled <- replace_string(all_sites_culled, old = "Site 925", new = "925")
all_sites_culled <- replace_string(all_sites_culled, old = "DSDP 390", new = "390")

unique(sort(all_sites_culled))


write.csv(all_sites_culled, file = "data/cache/all_sites_culled.csv")

