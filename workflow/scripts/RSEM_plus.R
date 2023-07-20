#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Combine counts.csv and rsem.csv files to perform
## transcript-isoform level analysis

### Load dependencies

library(data.table)
library(tidyverse)
library(optparse)
library(MASS)


# Helper functions that I will use on multiple occasions
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(1+exp(x))

# Likelihood function for mixture model
mixed_lik <- function(pnew, pold, TC, nT, n, logit_fn, p_sd = 1, p_mean = 0){
  logl <- sum(n*(log(inv_logit(logit_fn)*dbinom(TC, size = nT, prob = pnew) + (1-inv_logit(logit_fn))*dbinom(TC, nT, pold) ) ) ) + log(stats::dnorm(logit_fn, mean = p_mean, sd = p_sd))
  return(-logl)
}

### Process arguments

args = commandArgs(trailingOnly = TRUE)


### Read parameters

option_list <- list(
    make_option(c("-c", "--counts", type="character"),
                    default = ".",
                    help = "count.csv file path"),
    make_option(c("-r", "--rsem", type="character"),
                    default = ".",
                    help = 'rsem.csv file path'),
    make_option(c("-o", "--output", type = "character"),
                    default = ".",
                    help = 'output file path'),
    make_option(c("-s", "--sample", type = "character"),
                    default = "",
                    help = "sample name"),
    make_option(c("-e", "--echocode", type="logical"),
                    default = "FALSE",
                    help = 'print R code to stdout'),
    make_option(c("-n", "--pnew", type="double"),
                    default = -1,
                    help = 'mutation rate in new reads'),
    make_option(c("-b", "--pold", type="double"),
                    default = -1,
                    help = 'background mutation rate in old reads'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.

options(echo = as.logical(opt$echocode))

### Check that pnew and pold make sense
if(opts$pnew != -1 & opts$pnew < 0){
  stop("pnew must be > 0 or equal to -1! Set to -1 if you want RSEM+ to estimate pnew for you.")
}

if(opts$pold != -1 & opts$pold < 0){
  stop("pold must be > 0 or equal to -1! Set to -1 if you want RSEM+ to estimate pold for you.")
}


### Load csvs 
counts <- fread(opt$counts)
rsem <- fread(opt$rsem)

### Only keep relevant columns in counts.csv
    # Should eventually allow other mutation types to be analyzed
    # according to mutType argument
counts <- counts[,c("GF", "XF", "EF", "qname", "TC", "nT")]

cB <- counts[!grepl("__", XF), .N, by = .(XF, TC, nT)]

cT <- setDT(inner_join(rsem, counts, by = "qname"))

rm(counts)
rm(rsem)

### Estimate new and old mutation rate

if(opts$pnew == -1 & opts$pold == -1){

  ## USER PROVIDED NEITHER PNEW OR POLD

  # Binomial mixture likelihood
  mixture_lik <- function(param, TC, nT, n){

      logl <- sum(n*log(inv_logit(param[3])*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[2])^TC)*((1 -inv_logit(param[2]))^(nT-TC)) +  (1-inv_logit(param[3]))*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[1])^TC)*((1 - inv_logit(param[1]))^(nT-TC)) ) )

      return(-logl)

  }

  low_ps <- c(-9, -9, -9)
  high_ps <- c(0, 0, 9)

  fit <- stats::optim(par=c(-7,-2,0), mixture_lik, TC = cB$TC, nT = cB$nT,
                                  n = cB$N, method = "L-BFGS-B", 
                                  lower = low_ps, upper = high_ps)


  pnew <- inv_logit(max(fit$par[1:2]))
  pold <- inv_logit(min(fit$par[1:2]))

}else if(opts$pnew == -1 & opts$pold != -1){

  ## USER PROVIDED POLD BUT NOT PNEW

  # Binomial mixture likelihood
  mixture_lik <- function(param, TC, nT, n){

      logl <- sum(n*log(inv_logit(param[2])*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[1])^TC)*((1 -inv_logit(param[1]))^(nT-TC)) +  (1-inv_logit(param[2]))*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(opts$pold)^TC)*((1 - inv_logit(opts$pold))^(nT-TC)) ) )

      return(-logl)

  }

  low_ps <- c(logit(opts$pold), -9)
  high_ps <- c(0, 9)

  fit <- stats::optim(par=c(-2,0), mixture_lik, TC = cB$TC, nT = cB$nT,
                                  n = cB$N, method = "L-BFGS-B", 
                                  lower = low_ps, upper = high_ps)


  pnew <- inv_logit(fit$par[1])
  pold <- opts$pold

}else if(opts$pnew != -1 & opts$pold == -1){

  ## USER PROVIDED PNEW BUT NOT POLD


    # Binomial mixture likelihood
  mixture_lik <- function(param, TC, nT, n){

      logl <- sum(n*log(inv_logit(param[2])*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(opts$pnew)^TC)*((1 -inv_logit(opts$pnew))^(nT-TC)) +  (1-inv_logit(param[2]))*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[1])^TC)*((1 - inv_logit(param[1]))^(nT-TC)) ) )

      return(-logl)

  }

  low_ps <- c(-9, -9)
  high_ps <- c(logit(opts$pnew), 9)

  fit <- stats::optim(par=c(-7,0), mixture_lik, TC = cB$TC, nT = cB$nT,
                                  n = cB$N, method = "L-BFGS-B", 
                                  lower = low_ps, upper = high_ps)


  pnew <- opts$pnew
  pold <- inv_logit(fit$par[1])

}else{

  pnew <- opts$pnew
  pold <- opts$pold

}




### Get priors



# Add mutation rates
cT$pnew <- pnew
cT$pold <- pold

cB$pnew <- pnew
cB$pold <- pold

# Likelihood function for mixture model
mixed_lik <- function(pnew, pold, TC, nT, n, logit_fn, p_sd = 1, p_mean = 0){
  logl <- sum(n*(log(inv_logit(logit_fn)*dbinom(TC, size = nT, prob = pnew) + (1-inv_logit(logit_fn))*dbinom(TC, nT, pold) ) ) ) + log(stats::dnorm(logit_fn, mean = p_mean, sd = p_sd))
  return(-logl)
}

# Estimate GF for prior
Fn_prior <- cB %>% dplyr::ungroup() %>%
  dplyr::group_by(XF, TC, nT, pnew, pold) %>%
  dplyr::summarise(n = sum(N)) %>%
  dplyr::group_by(XF) %>%
  dplyr::summarise(logit_fn_rep = stats::optim(0, mixed_lik, nT = nT, TC = TC, n = n, pnew = pnew, pold = pold, method = "L-BFGS-B", lower = -7, upper = 7)$par, nreads =sum(n), .groups = "keep") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(prior = inv_logit(logit_fn_rep)) %>%
  dplyr::select(XF, prior)


# Add prior info
cT <- setDT(inner_join(cT, Fn_prior, by = "XF"))

# Calculate logit(fn) with analytical Bayesian approach
Fn_est <- cT[,.(fn_est = (sum((pt*dbinom(TC, nT, pnew)*prior)/(dbinom(TC, nT, pnew)*prior + dbinom(TC, nT, pold)*(1-prior)) ))/(sum(pt)),
                nreads = sum(pt)), by = .(XF, TF)]

Fn_est$sample <- opt$sample

write_csv(Fn_est, file = opt$output)

