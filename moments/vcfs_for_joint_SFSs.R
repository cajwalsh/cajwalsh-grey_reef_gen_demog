
## Rscript --vanilla

setwd("data")
library(radiator)
set.seed(1284)

## Get inds file
strata <- read.table("recall_pops_strata_nogb02.txt", header = T)

## Make id blacklist
idblist <- data.frame(INDIVIDUALS = strata[which(strata$STRATA != "north_gbr" &
                                                   strata$STRATA != "indo"),1])

## Import vcf and create dataset with just indo and north_gbr individuals
ngbr_indo_data <- read_vcf(data = "Walsh_Camblyrhynchos.vcf",
                           strata = strata,
                           blacklist.id = idblist,
                           filter.monomorphic = FALSE,
                           filter.common.markers = FALSE,
                           vcf.stats = FALSE)

## Export as vcf for easySFS to get non-thinned SFS to get sumSFS
## needed in first length recalibration equation in Appendix S1
genomic_converter(data = ngbr_indo_data,
                  filter.monomorphic = FALSE,
                  filter.common.markers = FALSE,
                  output = "vcf",
                  filename = "../../results/moments/ngbr_indo_un-thinned")

## Now filter out sites that are monomorphic in both populations
ngbr_indo_no_monomorphs <- tidy_genomic_data(data = ngbr_indo_data,
                                             filter.monomorphic = TRUE,
                                             filter.common.markers = FALSE)

## (Because initial attempts to do this using version of radiator I was using were
## causing a variable not needed in the analysis to go unfound or something)
source("../otheR_scripts/modified_filter_ld.R")

## Filter for short-distance linkage disequilibrium (i.e. thin),
ngbr_indo_thinned <- ld_filter_attempt(data = ngbr_indo_no_monomorphs,
                                       filter.monomorphic = FALSE,
                                       filter.common.markers = FALSE,
                                       interactive.filter = FALSE,
                                       filter.short.ld = "random")

## Export as a vcf for to get 2D-SFS for joint demographic history modelling
genomic_converter(data = ngbr_indo_thinned,
                  filter.monomorphic = FALSE,
                  filter.common.markers = FALSE,
                  output = "vcf",
                  filename = "../../results/moments/ngbr_indo_thinned")
