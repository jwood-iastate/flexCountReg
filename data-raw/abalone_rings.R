# Prepare data/abalone_rings.rda from the UCI Machine Learning Repository
# Abalone dataset (Nash et al., 1994; doi:10.24432/C55C7W), licensed under
# CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).
#
# Run from the package root during development (not during installation):
#   Rscript data-raw/abalone_rings.R
# The publisher's archive is kept in data-raw/abalone/abalone.zip. The script
# never accepts a changed archive silently: the SHA-256 hashes of the archive
# and of the extracted data file must match the pinned values below, which
# are also recorded in data-raw/source-manifest.csv.
#
# Changes made to the source data (all recorded in ?abalone_rings):
#   * lower-case snake_case column names (source names in the comment below);
#   * sex stored as a factor with levels F, I, M;
#   * rings stored as integer.
# No rows are removed and no values are rescaled: the continuous
# measurements keep the source scaling (millimetres or grams divided by 200).

archive <- "data-raw/abalone/abalone.zip"
url <- "https://archive.ics.uci.edu/static/public/1/abalone.zip"
archive_sha256 <- paste0("755a6a67c5b266961a3f149ea13be2cf",
                         "b6e6c727e48cee3c83bc0b4526210ee4")
data_sha256 <- paste0("de37cdcdcaaa50c309d514f248f7c230",
                      "2a5f1f88c168905eba23fe2fbc78449f")

if (!file.exists(archive)) {
  dir.create(dirname(archive), recursive = TRUE, showWarnings = FALSE)
  utils::download.file(url, archive, mode = "wb")
}
stopifnot(identical(digest::digest(file = archive, algo = "sha256"),
                    archive_sha256))

members <- utils::unzip(archive, list = TRUE)$Name
member <- members[basename(members) == "abalone.data"]
stopifnot(length(member) == 1L)
tmp <- tempfile()
dir.create(tmp)
utils::unzip(archive, files = member, exdir = tmp)
data_file <- file.path(tmp, member)
stopifnot(identical(digest::digest(file = data_file, algo = "sha256"),
                    data_sha256))

raw <- utils::read.csv(data_file, header = FALSE, stringsAsFactors = FALSE)
unlink(tmp, recursive = TRUE)
stopifnot(nrow(raw) == 4177L, ncol(raw) == 9L)

# Source names: Sex, Length, Diameter, Height, Whole weight, Shucked weight,
# Viscera weight, Shell weight, Rings
names(raw) <- c("sex", "length", "diameter", "height", "whole_weight",
                "shucked_weight", "viscera_weight", "shell_weight", "rings")
stopifnot(
  all(raw$sex %in% c("F", "I", "M")),
  !anyNA(raw),
  all(vapply(raw[2:8], is.numeric, logical(1))),
  all(is.finite(raw$rings)), all(raw$rings >= 0),
  all(raw$rings == floor(raw$rings))
)

abalone_rings <- raw
abalone_rings$sex <- factor(abalone_rings$sex, levels = c("F", "I", "M"))
abalone_rings$rings <- as.integer(abalone_rings$rings)
rownames(abalone_rings) <- NULL

save(abalone_rings, file = "data/abalone_rings.rda", compress = "xz",
     version = 2)
message("data/abalone_rings.rda SHA-256: ",
        digest::digest(file = "data/abalone_rings.rda", algo = "sha256"))
