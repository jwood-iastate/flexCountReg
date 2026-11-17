# Prepare data/stream_salamanders.rda from the original Dryad deposit
#   Price, S. J., Muncy, B. L., Bonner, S. J., Drayer, A. N. and Barton, C. D.
#   (2016). Data from: Effects of mountaintop removal mining and valley
#   filling on the occupancy and abundance of stream salamanders. Dryad.
#   doi:10.5061/dryad.5m8f6. License: CC0 1.0.
#
# Run from the package root during development (not during installation):
#   Rscript data-raw/stream_salamanders.R
# Dryad does not allow unauthenticated scripted downloads, so the file was
# downloaded by hand from https://datadryad.org/dataset/doi:10.5061/dryad.5m8f6
# (version of 2016-11-17) and is kept in data-raw/salamanders/. Its MD5
# matches the digest published by Dryad; the SHA-256 is pinned below.
#
# The source is one row per site (23 sites) with, for each of the four
# samples (suffix .1 to .4), the sampling covariates and the counts of seven
# species/life-stage groups (columns c.<group>.<sample>). The reshape gives
# one row per site x sample x group (23 x 4 x 7 = 644 rows). Changes:
#   * site labels normalized ("VF - 4", " R -6" -> "VF-4", "R-6");
#   * Mined (0/1) stored as factor mined with levels "no", "yes";
#   * group codes upper-cased with "-" (c.ec.a -> "EC-A");
#   * columns renamed (Study Site -> site, Cover -> cover, DOP.k -> DOP,
#     Wtemp.k -> Wtemp, DOY.k -> DOY).
# The covariates were already centred and scaled in the source; they are
# kept as they are. No counts are changed.

src <- "data-raw/salamanders/Price et al. JAPPL-2015-count-data.csv"
src_md5 <- "a9f102cb6f6b51da4387b9e64d62d588"
src_sha256 <- paste0("6d4fa18c272add2c6370f9e49c3f24b1",
                     "82e9786e7d527eb6a78c1d9c8fdbd21e")
stopifnot(file.exists(src),
          identical(digest::digest(file = src, algo = "md5"), src_md5),
          identical(digest::digest(file = src, algo = "sha256"), src_sha256))

raw <- utils::read.csv(src, check.names = FALSE, stringsAsFactors = FALSE)
raw <- raw[!is.na(raw$Mined), ]   # drop the blank trailing line
stopifnot(nrow(raw) == 23L, ncol(raw) == 43L)

site <- gsub("\\s+", "", raw[["Study Site"]])
stopifnot(!anyDuplicated(site), all(grepl("^(R|VF)-[0-9]+$", site)))
num <- as.integer(sub(".*-", "", site))
site_levels <- c(paste0("R-", sort(num[startsWith(site, "R")])),
                 paste0("VF-", sort(num[startsWith(site, "VF")])))

groups <- c(gp = "GP", pr = "PR", dm = "DM", ec.a = "EC-A", ec.l = "EC-L",
            des.l = "DES-L", df = "DF")

rows <- list()
for (g in names(groups)) {
  for (k in 1:4) {
    count <- raw[[paste0("c.", g, ".", k)]]
    stopifnot(is.numeric(count), !anyNA(count), all(count >= 0),
              all(count == round(count)))
    rows[[length(rows) + 1L]] <- data.frame(
      site = site,
      mined = ifelse(raw$Mined == 1, "yes", "no"),
      cover = raw$Cover,
      sample = k,
      DOP = raw[[paste0("DOP.", k)]],
      Wtemp = raw[[paste0("Wtemp.", k)]],
      DOY = raw[[paste0("DOY.", k)]],
      spp = groups[[g]],
      count = as.integer(count),
      stringsAsFactors = FALSE
    )
  }
}
stream_salamanders <- do.call(rbind, rows)
stream_salamanders$site <- factor(stream_salamanders$site,
                                  levels = site_levels)
stream_salamanders$mined <- factor(stream_salamanders$mined,
                                   levels = c("no", "yes"))
stream_salamanders$sample <- as.integer(stream_salamanders$sample)
stream_salamanders$spp <- factor(stream_salamanders$spp,
                                 levels = unname(groups))
rownames(stream_salamanders) <- NULL

stopifnot(
  nrow(stream_salamanders) == 644L,
  !anyDuplicated(stream_salamanders[c("site", "sample", "spp")]),
  all(table(stream_salamanders$site) == 28L),
  !anyNA(stream_salamanders)
)

# Cross-check against the independently prepared glmmTMB::Salamanders object
# (if installed), matching rows by site, sample and species
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  e <- new.env()
  utils::data("Salamanders", package = "glmmTMB", envir = e)
  ref <- e$Salamanders
  key <- function(d) paste(d$site, d$sample, d$spp)
  m <- match(key(stream_salamanders), key(ref))
  stopifnot(!anyNA(m), nrow(ref) == 644L)
  for (v in c("count", "cover", "DOP", "Wtemp", "DOY")) {
    stopifnot(isTRUE(all.equal(as.numeric(stream_salamanders[[v]]),
                               as.numeric(ref[[v]][m]))))
  }
  stopifnot(identical(as.character(stream_salamanders$mined),
                      as.character(ref$mined[m])))
  message("Matches glmmTMB::Salamanders (glmmTMB ",
          utils::packageVersion("glmmTMB"), ") on all keys and values.")
}

save(stream_salamanders, file = "data/stream_salamanders.rda",
     compress = "xz", version = 2)
message("data/stream_salamanders.rda SHA-256: ",
        digest::digest(file = "data/stream_salamanders.rda",
                       algo = "sha256"))
