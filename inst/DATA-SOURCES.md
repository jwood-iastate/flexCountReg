# Sources of bundled datasets

## abalone_rings

- Source: UCI Machine Learning Repository, Abalone dataset,
  doi:10.24432/C55C7W.
- Archive: <https://archive.ics.uci.edu/static/public/1/abalone.zip>,
  retrieved 2026-11-17. SHA-256 of the archive:
  `755a6a67c5b266961a3f149ea13be2cfb6e6c727e48cee3c83bc0b4526210ee4`;
  of `abalone.data`:
  `de37cdcdcaaa50c309d514f248f7c2302a5f1f88c168905eba23fe2fbc78449f`.
- Original study: Nash, W. J., Sellers, T. L., Talbot, S. R., Cawthorn,
  A. J. and Ford, W. B. (1994). The population biology of abalone
  (*Haliotis* species) in Tasmania. I. Blacklip abalone (*H. rubra*) from the
  north coast and islands of Bass Strait. Sea Fisheries Division, Technical
  Report No. 48.
- License: CC BY 4.0.
- Preparation: `data-raw/abalone_rings.R` in the source repository (not
  included in the built package). The manifest
  `data-raw/source-manifest.csv` records the source and output hashes.

## stream_salamanders

- Source: Price, S. J., Muncy, B. L., Bonner, S. J., Drayer, A. N. and
  Barton, C. D. (2016). Data from: Effects of mountaintop removal mining and
  valley filling on the occupancy and abundance of stream salamanders.
  Dryad, doi:10.5061/dryad.5m8f6 (version 2016-11-17).
- File: `Price et al. JAPPL-2015-count-data.csv`, downloaded 2026-11-17.
  MD5 `a9f102cb6f6b51da4387b9e64d62d588` (as published by Dryad); SHA-256
  `6d4fa18c272add2c6370f9e49c3f24b182e9786e7d527eb6a78c1d9c8fdbd21e`.
- Study: Price et al. (2016), Journal of Applied Ecology 53(2), 459-468,
  doi:10.1111/1365-2664.12585.
- License: CC0 1.0.
- Preparation: `data-raw/stream_salamanders.R` (wide-to-long reshape,
  checked against `glmmTMB::Salamanders`, glmmTMB 1.1.14).

## washington_roads

- Source: Washington State Department of Transportation crash, roadway and
  traffic data (2016-2018) obtained through the Highway Safety Information
  System, <https://highways.dot.gov/research/safety/hsis>.
