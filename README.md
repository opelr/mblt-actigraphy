# H4085 MBLT Actigraphy

Processing of actigraphy data related for H4085. This project allows all other repositories in the H4085 family to use the same, consistent data structures for actigraphy data.

## Table of Contents

- [Getting Started](#getting-started)
- [Usage](#usage)
- [Development](#development)
- [Credits](#credits)
- [License](#license)

## Getting Started

Below are instructions for getting the project up-and-running on your local machine.

### Prerequisites

Before interacting with this repository, users will need [`R`](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/#Desktop) installed on their local machine.

Users should ensure that the following R packages are installed: `tidyverse`, `Hmisc`, `zoo`, and `openxlsx`.

If you are planning on contributing to the code, you'll need to install [`git`](https://git-scm.com/). For a tutorial on integrating Git and R, check out [Happy Git with R](http://happygitwithr.com/).

### Installation

1. Either `git clone` (_preferred_) or download the repository.
2. Open `mblt-actigraphy.Rproj` in RStudio.

## Usage

### Download Necessary Files

If using the repository specifically for H4085 purposes, you'll need to download the Participant Master List and the Actiware Database from the OHSU Box and the rMEQ and PCL surveys from OHSU REDCap.

The Master List and surveys should be renamed and placed in `/data/raw/`. Actiware files should be exported from the database following the protocol specified in the H4085 project documentation.

### Automated

Open `src/directory_monitory.ps1` in a text editor and change the file paths to match your machine. This script will when check if any untracked or modified exported Actigraphy files exist in that directory, and will create the various serialized R data object files (`.rds`) using the new information and store them in the `/data/rds/` folder.

This script can be set to run on a `chron` job for automated monitoring and data integration.

Finally, the `.rds` files from this repository can be copied into any new project that uses H4085 actigraphy data.

## Development

### Contributing

Before contributing any major changes, developers should checkout a `development` branch. Before merging changes with `master`, these major additions should be reviewed with a colleague who's familiar with the project.

Code style should adhere to the [tidyverse style guide](http://style.tidyverse.org/).

### Running Tests

#### Functional Tests

This project currently has no functional tests.

#### Conceptual Tests

Two test scripts exist in the `/tests/` folder to illustrate how IS/IV work and another for sample size. Run them as you would any `R` script.

#### Style Tests

Developers can use the `styler` package to automatically format their code.

### Versioning

Our group uses [Semantic Versioning](http://semver.org/) for versioning.

## Credits

### Authors

- **Ryan Opel** - *Initial contributions* - ryan.a.opel@gmail.com

### Built With

- [Philips Actiware](http://www.actigraphy.com/solutions/actiware/) - Actigraphy devices and database program

## License

This project is licensed under the GNU General Public License v2.0 - see [LICENSE](LICENSE) for details.
