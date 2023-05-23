#!/usr/bin/env Rscript
# R script for normalizing TimeLapse data
# July 17, 2017
# Oct 30, 2018 updated to prevent scale file with NaN/Inf entries.
# May 12, 2021 rewrite by Martin to allow stand-alone run
#              processing of multiple file locations at once (in case normalization between several TimeLapse experiments is needed)
#              option to calculate normalization factors for spikeins only
# Simon Lab
# Output a scale file for making tracks.

args = commandArgs(trailingOnly = TRUE)


# Read parameters
    library(optparse)

    option_list <- list(
        make_option(c("-d", "--dirs", type="character"),
                     default = ".",
                     help = "comma separated list of directories with mature.*.txt files"),
        make_option(c("-s", "--spikename", type="character"),
                     default = "",
                     help = 'character string that is common identifier of spikeins gene_name'),
        make_option(c("-e", "--echocode", type="logical"),
                     default = "FALSE",
                     help = 'print R code to stdout'))

    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser) # Load options from command line.

# Set R code printing for debug mode
    options(echo = as.logical(opt$echocode))


# Load libraries:
    library(tidyverse)
    library(edgeR)
    library(MASS)


# Set date

    date <- Sys.Date()
    date

# Define the data frame

    master <- tibble()


# # Load and clean master data

    # Define samples to load
        ## Martin's old code
	# dirs <- opt$dirs %>% str_split(',') %>% unlist()

	# Screw it, going to hardcode directory cause I can
	dirs <- paste0(getwd(), "/results/htseq/")
    print(dirs)

        samplefiles <- list.files(path = dirs,
                                  pattern = 'mature.*.txt',
                                  recursive = FALSE)

	## Martin's old code that doesn't work in my workflow
        # samplenames <- gsub(".dir.*", "", samplefiles)
        samplenames <- paste0(dirs, samplefiles)

    print(samplenames)
	# Actual sample name wildcards
	snames <- gsub(".*mature.|_htseq.*", "", samplefiles)
	#snames <- gsub("htseq*", "", snames)


    print(snames)
    # Loop through the samples to load them:

    for (i in seq_along(samplenames)){

        df <- read_tsv(samplenames[i],
                       col_names = c('gene', 'count'),
                       col_types = 'ci')

        df$sample <- paste0(snames[i])

        master <- rbind(master, df)

    }

    rm(df)

    # glimpse(master)

    m <- master %>%
        pivot_wider(names_from = sample, values_from = count) %>%
        filter(!grepl('(__not_aligned|__alignment_not_unique|__ambiguous|__no_feature|__too_low_aQual)', gene))

    # Filter only for spikeins. gene_name (or other identifier must contain universal character string opt$spikename)
    if (opt$spikename != "") {
        m <- m %>%
                filter(grepl(opt$spikename, gene))

        print(paste0("* Calculating normalization factors from spikeins with name *", opt$spikename))
    }

    # head(m, 20)
    rnames <- m[,1]
    m <- as.matrix(m[,-1])
    storage.mode(m) <- "numeric"
    rownames(m) <- unlist(rnames)
    # head(m, 20)

    print(rnames)

# Correlation analysis
    #cor(m)

# Examine RNA seq by edgeR:

# Add dummy IDs (don't need them yet)
    ers <- rep("X", time = dim(m)[2])
    ers <- factor(ers)

# Filter to only data-containing genes:
    keep <- rowSums(m > 0) > 1
    erde <- DGEList(m[keep,], group = ers)
    erde <- calcNormFactors(erde, method = "upperquartile")
    scale <- erde$samples$norm.factors*erde$samples$lib.size
    # scale

if (sum(is.finite(scale)) == length(scale)){

    x <- mean(scale)/scale

    sdf <- tibble(sample = snames,
                  scale = x)

    print(sdf)

    write.table(sdf, file = 'scale',
                     quote = FALSE,
                     row.names = FALSE,
                     col.names = FALSE)

} else {
    paste('no acceptable scale values')
}

