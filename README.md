# Bootstrap Methods

A R implementation of statistical EM algorithm for parameter estimation and confidence interval construction and bootstrap methods.

## General Information

**EM Algorithm**: Implementation for estimating parameters of exponential distributions with missing data

**Bootstrap Methods**: Statistical resampling techniques for constructing confidence intervals for various estimators

## Repository structure

```
Bootstrap-Methods/
├── .gitignore
├── emAlgoboot.R         
├── README.md             
├── LICENSE              
└── Sujet.pdf            
```

## Dependencies

### R Installation

Install R from [CRAN](https://cran.r-project.org/)

### R Package Dependencies

Install the required R packages by running the following commands in R:

```r
# Install required packages
install.packages(c("dplyr", "boot"))

# Load packages to verify installation
library(dplyr)
library(boot)
```

## Usage

> Rscript emAlgoboot.R

## Results

The script outputs:
- Estimated λ parameter from the EM algorithm
- Confidence intervals for all four estimators using three different methods
- Optimal bootstrap sample sizes for each estimator and interval type


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
