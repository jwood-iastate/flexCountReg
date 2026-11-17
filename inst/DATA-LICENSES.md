# Licenses of bundled datasets

The flexCountReg package code is licensed under the MIT license (see
`LICENSE`). Some bundled datasets are third-party data distributed under
their own terms, which the MIT license does not replace.

| Dataset | Terms | Attribution |
|---|---|---|
| `washington_roads` | Compiled by the package authors from Highway Safety Information System (HSIS) data | Washington State Department of Transportation via HSIS, <https://highways.dot.gov/research/safety/hsis> |
| `abalone_rings` | Creative Commons Attribution 4.0 International (CC BY 4.0), <https://creativecommons.org/licenses/by/4.0/> | Nash, W., Sellers, T., Talbot, S., Cawthorn, A. and Ford, W. (1994). Abalone (dataset). UCI Machine Learning Repository. doi:10.24432/C55C7W |
| `stream_salamanders` | CC0 1.0 public-domain dedication, <https://creativecommons.org/publicdomain/zero/1.0/> | Price, S. J., Muncy, B. L., Bonner, S. J., Drayer, A. N. and Barton, C. D. (2016). Data from: Effects of mountaintop removal mining and valley filling on the occupancy and abundance of stream salamanders. Dryad. doi:10.5061/dryad.5m8f6 |

Changes made to `abalone_rings`: columns renamed to lower-case snake_case,
`sex` stored as a factor and `rings` as an integer. No rows or values were
changed. `stream_salamanders` was reshaped from one row per site to one row
per site, sample and species group (see `?stream_salamanders`). See also
`DATA-SOURCES.md`.
