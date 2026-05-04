# Changelog

All notable changes to nanoID will be documented in this file.

## [0.2.0] – 2026‑05‑04

### Added

- Automatic selection of the abundance dominance parameter (κ) via grid search and knee‑point detection when κ is not explicitly specified. Automatic selection is now the default behavior.

### Changed

- Removed the option to suppress per‑split singleton consensus sequences (`--ignore_singletons` parameter for `condens`), as retaining only consensus nodes shared across splits is sufficient.
- Lowered the minimum identity threshold for ASV quantification to 0.96 (`--min_identity` parameter for `quant`).

## [0.1.0] – 2026‑05‑01

### Added

- First public release.
