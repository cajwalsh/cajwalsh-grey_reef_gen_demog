## Base directory
setwd("~/Desktop/shark_nf/")

## Get inds file
setwd("data")
mystrata <- read.table("recall_pops_strata_nogb02.txt", header = T)

## Load project vcf
library(radiator)
vcf_data <- read_vcf(data = "Walsh_Camblyrhynchos.vcf",
                     strata = mystrata,
                     filter.monomorphic = FALSE,
                     filter.common.markers = FALSE,
                     vcf.stats = FALSE)

## Turn into tidy (should probably be able to do this and previous step in one
## using only this one calling directly to vcf file but it was causing my R
## session to abort session so I do it in two here...maybe not enough RAM)
tidy_data <- tidy_genomic_data(data = vcf_data,
                               strata = mystrata,
                               filter.monomorphic = FALSE,
                               filter.common.markers = FALSE,
                               vcf.stats = FALSE)

## Create chromosome/loci file necessary for NeEstimator to do LDNe method with
## inter-locus comparisons only (by pretending each locus is on its own
## chromosome within using a two column df with matching identical locus and
## chromosome names) across all populations (including all loci)
setwd("../results/NeEstimator")
neestim_chrom_df <- data.frame(CHROM = unique(tidy_data$LOCUS),
                               LOCUS = unique(tidy_data$LOCUS))
write.table(neestim_chrom_df, file = "chrom.txt",
            quote = F, sep = " ", row.names = F, col.names = F)

## For each pop:
setwd("tmp")
for(pop in unique(mystrata$STRATA)) {

idblist <- data.frame(INDIVIDUALS = mystrata[which(mystrata$STRATA != pop),1])

### get correct individuals,
pop_data <- tidy_genomic_data(data = tidy_data,
                              blacklist.id = idblist,
                              filter.monomorphic = FALSE,
                              filter.common.markers = FALSE)

### filter out sites monomorphic in this population,
pop_no_monomorphs <- tidy_genomic_data(data = pop_data,
                                       filter.monomorphic = TRUE,
                                       filter.common.markers = FALSE)

### write genepop file for LDNe in NeEstimator
genomic_converter(data = pop_no_monomorphs,
                  output = "genepop",
                  filename = paste0(as.character(pop),
                                    "_no_monomorphs"))
}
