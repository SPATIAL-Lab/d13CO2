library(rgplates)

## ---- optional, for withTimeout() ----
if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils")
}
withTimeout <- R.utils::withTimeout

comp <- read.csv("Phan/PhanData/assignLatLonSite/PhanCompUpdated_WithComponent.csv")

## make sure ages are in Ma
if (max(comp$age, na.rm=TRUE) > 2e3) comp$age <- comp$age/1e6 else
  if (max(comp$age, na.rm=TRUE) > 600)  comp$age <- comp$age/1e3

## validate rows
ok <- is.finite(comp$assigned.lon) & comp$assigned.lon >= -180 & comp$assigned.lon <= 180 &
  is.finite(comp$assigned.lat) & comp$assigned.lat >=  -90 & comp$assigned.lat <=   90 &
  is.finite(comp$age)          & comp$age >= 0 & comp$age <= 750
comp <- comp[ok, , drop=FALSE]

## de-duplicate requests: round age so near-equals share a call (tune digits)
age_digits <- 1     # 0.1 Ma; if still too many uniques, try 0.2–0.5
comp$age_key <- round(comp$age, age_digits)

key <- paste(comp$assigned.lon, comp$assigned.lat, comp$age_key)
uniq <- !duplicated(key)
U <- comp[uniq, c("assigned.lon","assigned.lat","age_key")]
rownames(U) <- key[uniq]

## on-disk cache so you can resume without losing progress
cache_file <- "paleocoords_cache.rds"
cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()

## --- knobs you can tweak ---
QPS <- 1.5             # polite rate limit (requests/sec)
DT  <- 1/QPS           # sleep between calls
TRIES <- 4             # per-point tries
TIMEOUT_SEC <- 25      # HARD timeout per call so the loop can't hang
CHK_EVERY <- 50        # checkpoint cache every N new requests
HEARTBEAT_EVERY <- 300 # seconds; print a heartbeat line
## ---------------------------

recon_one <- function(lon, lat, age_ma, tries=TRIES) {
  for (k in seq_len(tries)) {
    res <- try(
      withTimeout(
        reconstruct(c(lon,lat), age=age_ma, model="PALEOMAP"),
        timeout = TIMEOUT_SEC, onTimeout = "error"
      ),
      silent = TRUE
    )
    if (!inherits(res, "try-error") && all(is.finite(res))) return(res)
    ## backoff with a little jitter
    Sys.sleep(0.3*k + runif(1, 0, 0.2))
  }
  c(NA_real_, NA_real_)
}

cat("Unique requests: ", nrow(U), "\n", sep="")
t0 <- Sys.time(); last_beat <- t0
n_done <- 0L
n_new  <- 0L

for (i in seq_len(nrow(U))) {
  k <- rownames(U)[i]
  if (is.null(cache[[k]])) {
    cache[[k]] <- recon_one(U$assigned.lon[i], U$assigned.lat[i], U$age_key[i])
    n_new <- n_new + 1L
    if (n_new %% CHK_EVERY == 0) {
      saveRDS(cache, cache_file)
      cat(sprintf("\rSaved cache at %d/%d (%.1f%%) | elapsed %s",
                  i, nrow(U), 100*i/nrow(U),
                  format(Sys.time()-t0, digits=2)), sep="")
      flush.console()
    }
    Sys.sleep(DT)  # be polite to the service
  }
  n_done <- n_done + 1L
  
  ## heartbeat every HEARTBEAT_EVERY seconds so you can see it's alive
  if (as.numeric(difftime(Sys.time(), last_beat, units="secs")) > HEARTBEAT_EVERY) {
    done_keys <- sum(!vapply(cache, is.null, logical(1)))
    cat(sprintf("\nHeartbeat: %d/%d unique done (%.1f%%). Elapsed %s.\n",
                done_keys, nrow(U), 100*done_keys/nrow(U),
                format(Sys.time()-t0, digits=2)))
    flush.console()
    last_beat <- Sys.time()
  }
}

saveRDS(cache, cache_file); cat("\nDone unique requests.\n")

## map cached results back to all rows
get_xy <- function(lon, lat, age_key) {
  k <- paste(lon, lat, age_key)
  v <- cache[[k]]; if (is.null(v)) c(NA_real_, NA_real_) else v
}
coords <- t(mapply(get_xy, comp$assigned.lon, comp$assigned.lat, comp$age_key))
colnames(coords) <- c("paleo_lon","paleo_lat")

fail <- which(!is.finite(coords[,1]) | !is.finite(coords[,2]))
if (length(fail)) message("Reconstruction failed for ", length(fail), " rows (after retries/timeout).")
