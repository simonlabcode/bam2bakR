#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Merge feature assignment and mutation counts tables

# Load dependencies ------------------------------------------------------------

library(optparse)
library(data.table)
library(readr)

# Process parameters -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("-g", "--genes", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to genes"),
  make_option(c("-e", "--exons", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exons"),
  make_option(c("-b", "--exonbins", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exonbins"),
  make_option(c("-t", "--transcripts", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to transcripts"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to modified annotation output"),
  make_option(c("-s", "--sample", type = "character"),
              help = "Sample name")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Combine tables --------------------------------------------------------------

sample <- paste0(opt$sample, "_counts")

print(paste0("sample is: ", sample))


### Load mutation counts
muts_file <- list.files(path = "./results/counts/",
                        pattern = sample,
                        full.names = TRUE)[1]

muts <- fread(muts_file)



# merge with gene assignments
if(opt$genes){
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
  genes_file <- list.files("./results/featurecounts_genes/",
                           pattern = sample, full.names = TRUE)[1]
  
  message("genes_file is:")
  print(genes_file)

  genes <- fread(genes_file)
  
  colnames(genes) <- c("qname", "status", "nhits", "GF")
  
  genes <- genes[ nhits > 0 , c("qname", "GF")]
  
  genes[, GF := gsub(",", "+", GF)]

  
  muts <- genes[muts, on = .(qname)]
  
  
}


# merge with exon assignments
if(opt$exons){
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
  exons_file <- list.files("./results/featurecounts_exons/",
                           pattern = sample, full.names = TRUE)[1]
  
  message("exons_file is:")
  print(exons_file)
  
  exons <- fread(exons_file)
  
  colnames(exons) <- c("qname", "status", "nhits", "XF")
  
  exons <- exons[ nhits > 0 , c("qname", "XF")]
  
  exons[, XF := gsub(",", "+", XF)]
  
  muts <- exons[muts, on = .(qname)]
  
  
}


# Merge with exonbin assignments
if(opt$exonbins){
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
  exonbins_file <- list.files("./results/featurecounts_exonbins/",
                           pattern = sample, full.names = TRUE)[1]
  
  exonbins <- fread(exonbins_file)
  
  colnames(exonbins) <- c("qname", "status", "nhits", "exon_bin")
  
  
  exonbins <- exonbins[ nhits > 0 , c("qname", "exon_bin")]
  
  
  exonbins[, exon_bin := gsub(",", "+", exon_bin)]
  
  
  muts <- exonbins[muts, on = .(qname)]
  
  
}


# Merge with transcript assignments
if(opt$transcripts){
  
  sample <- paste0("^", opt$sample, ".s.bam.featureCounts")
  
  transcripts_file <- list.files("./results/featurecounts_transcripts/",
                              pattern = sample, full.names = TRUE)[1]
  
  transcripts <- fread(transcripts_file)
  
  colnames(transcripts) <- c("qname", "status", "nhits", "transcripts")
  
  
  transcripts <- transcripts[ nhits > 0 , c("qname", "transcripts")]
  
  transcripts[, transcripts := gsub(",", "+", transcripts)]
  
  
  muts <- transcripts[muts, on = .(qname)]
  
  
}


# Write to final output
write_csv(muts,
          file = opt$output)
