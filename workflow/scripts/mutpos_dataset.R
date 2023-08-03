### PURPOSE OF THIS SCRIPT
## Make a dataset that can be accessed via the arrow package
## in R, to allow for working with larger than memory mutation
## position data.

# Load dependencies ----------------------------------------------------------

library(arrow)
library(dplyr)
library(optparse)
library(MASS)



# Parse parameters -----------------------------------------------------------

### Process arguments

args = commandArgs(trailingOnly = TRUE)


### Read parameters

option_list <- list(
    make_option(c("-m", "--mutpos", type="character"),
                    default = ".",
                    help = "mutpos.csv.gz file path"),
    make_option(c("-o", "--output", type = "character"),
                    default = ".",
                    help = 'output file path'),
    make_option(c("-e", "--echocode", type="logical"),
                    default = "FALSE",
                    help = 'print R code to stdout'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.

options(echo = as.logical(opt$echocode))


# Create dataset -------------------------------------------------------------

ds <- open_dataset(opt$mutpos, format = "csv",
                   schema = schema(
                     sample = string(),
                     rname = string(),
                     gloc = int64(),
                     GF = string(),
                     XF = string(),
                     ai = bool(),
                     tp = string(),
                     trials = int64(),
                     n = int64()
                   ),
                   skip_rows=1)

ds %>%
    select(-ai) %>%
    group_by(sample, tp) %>%
    write_dataset(opt$output, format = "parquet")