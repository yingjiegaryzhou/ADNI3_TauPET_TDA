library(BGLR)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)

setwd("/Users/garyzhou/Downloads/Dropbox")
source("EC3D.R")
#source(SavePNGImages.R) if threshold_val not already defined

process_ecf <- function(q_val) {
  ### Set up the Parameters ###
  q <- round(q_val, 2)
  startdir = ""
  dir_base <- paste0("/Users/garyzhou/Downloads/Dropbox/hippo_tau7/Q", q)
  in.dir <- file.path(dir_base, "images")
  out.file <- file.path(dir_base, "MRIECs_hippo.RData")
  img.dir = NULL
  stepsize=100
  rotstep=72
  ### Run The Euler Characteristic Function ###
  ecf = ecf(in.dir = in.dir, out.file = out.file, img.dir = NULL, first.only = FALSE)
}

# Apply the function to each value in threshold_val
lapply(threshold_val, process_ecf)
