# Change Log
All notable changes to the *mblt-actigraphy* subdirectory under the *h4085* project will be documented in this file.
This project adheres to standards proposed by [Keep a CHANGELOG](http://keepachangelog.com/).

---

## [Unreleased]
### Changed
- Cleaning up code to adhere to [Hadley's R style guide](http://adv-r.had.co.nz/Style.html)

## [0.1.0] - 2017-07-12
### Added
- 01_ETL.R script reads Actigrapy CSV files; pulls header and scoring information to two separate data frames. 
    - Date-Time calculations and sunset function to stage based on time of year
- 02_Analysis.R script groundwork laid. Some insights and functions include:
    - Daytime activity ratio;
    - Nighttime sleep duration; and 
    - Percent sleep

## [EXAMPLE]
### Added
- New features
### Changed
- Changes in existing functionality
### Deprecated
- Once-stable features removed in upcoming releases
### Removed
- Deprecated features removed in this release
### Fixed
- Bug fixes
