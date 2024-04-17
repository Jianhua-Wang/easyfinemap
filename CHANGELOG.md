# Changelog

## [0.4.6] - 2024-04-17

### Added

### Changed

### Fixed
- fix bugs in cojo cond

## [0.4.5] - 2024-04-17

### Added

### Changed

### Fixed
- fix bugs in polyfun+finemap and polyfun+susie


## [0.4.4] - 2023-11-09

### Added

### Changed
- discard progress bar, because it is not compatible with multiprocessing

### Fixed
- fix bugs in cojo cond


## [0.4.3] - 2023-10-25

### Added

### Changed
- env add data.table
- fix inf beta and se in sumstats

### Fixed



## [0.4.2] - 2023-10-23

### Added

### Changed
- load sumstat after multiprocessing
- aviod error when no intersection ldref with input sumstat
- fix bugs in cojo cond

### Fixed


## [0.4.1] - 2023-10-15

### Added
- docs

### Changed

### Fixed


## [0.4.0] - 2023-10-12

### Added
- docs

### Changed

### Fixed

## [0.3.9] - 2023-09-26

### Added

### Changed

### Fixed
- fix dependency

## [0.3.5] - 2023-06-13

### Added

### Changed
- fix smunger independency
### Fixed



## [0.3.4] - 2023-06-12

### Added

### Changed
- fix bugs in loci by ldblock
### Fixed


## [0.3.3] - 2023-06-12

### Added

### Changed
- fix bugs in loci by ldblock
### Fixed



## [0.3.2] - 2023-05-25

### Added

### Changed
- remove pathos
### Fixed


## [0.3.1] - 2023-05-25

### Added
- locus plot

### Changed

### Fixed


## [0.3.0] - 2023-05-24

### Added
- annotate R2 for locus plot

### Changed

### Fixed


## [0.2.9] - 2023-05-24

### Added


### Changed
- Speed up susie by using fread

### Fixed


## [0.2.8] - 2023-05-23

### Added


### Changed
- Speed up conditional analysis by intersect ldref with input sumstat first

### Fixed


## [0.2.7] - 2023-05-20

### Added


### Changed
- Use tabix to load sumstats when finemap

### Fixed


## [0.2.6] - 2023-05-19

### Added

- Support polyfun for finemap and susie

### Changed
- Use the most significant SNP as lead SNP when COJO fails

### Fixed

## [0.2.5] - 2023-03-30

### Added
- set the most significant SNP as SNP, if there are no SNPs with P-value ≤ threshold.

### Changed



### Fixed

## [0.2.4] - 2023-03-30

### Added
- support susie

### Changed



### Fixed

## [0.2.3] - 2023-03-24

### Added
- suppor using LD blocks as loci boudaries.

### Changed



### Fixed



## [0.2.2] - 2023-03-22

### Added


### Changed
- deleted format function


### Fixed


## [0.2.1] - 2023-01-10

### Added

- Instruction of installation

### Changed


### Fixed
- typo in easyfinemap.py line519

## [0.2.0] - 2023-01-10

### Added

- add CAVIARBF

### Changed


### Fixed

## [0.1.4] - 2023-01-10

### Added

- add PAINTOR

### Changed


### Fixed

## [0.1.3] - 2023-01-09

### Added

- output credible sets

### Changed


### Fixed

## [0.1.2] - 2023-01-07

### Added

- update summary statistics using cojo-cond
- make ld matrix using plink --r2
- fine-mapping tools: FINEMAP

### Changed


### Fixed

## [0.1.1] - 2023-01-07

### Added

- add temp dir decorator
- identify lead SNPs by LD clumping
- identify lead SNPs by conditional analysis

### Changed


### Fixed

## [0.0.5] - 2022-12-26

### Added

- validate GWAS summary statistics

### Changed


### Fixed

## [0.0.4] - 2022-12-26

### Added

- prepare and validate LD reference panel

### Changed


### Fixed

## [0.0.1] - 2022-12-20

### Added

- merge the overlapped independent loci (optional).

### Changed


### Fixed

## [0.0.2] - 2022-12-22

### Added

- identify the independent lead snps by distance only
- expand the independent lead snps to independent loci by given range.

### Changed


### Fixed

## [0.0.3] - 2022-12-22

### Added

- extract LD ref plink bfile and clean it.

### Changed


### Fixed