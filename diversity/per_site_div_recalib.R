### Written by Paolo adapted by Cameron ###

### Recalibrate the sequenced length of each population for per-site
### standardization of nucleotide diversity calculated from SFS

### Rscript --vanilla 

### The following values are taken from the ddRAD_stats.txt file from iPyRAD
### var   sum_va
### 0     1192
### 1     2520
### 2     3355
### 3     3383
### 4     3035
### 5     2284

var0 <- 1192
var1 <- 2520
var2 <- 3355
var3 <- 3383
var4 <- 3035
var5 <- 2284

monomorphic <- var0 ## loci with 0 variable sites
singletons <- var1 ## loci with 1 variable site
polymorphic <- sum(c(var1, var2, var3, var4, var5)) ## sites with at least one variable site
initial_loci <- sum(c(monomorphic+polymorphic))

## Of the polymorphic loci, we kept 7952 loci with 14935 SNPs after quality filtering
filtered_loci <- 7952
filtered_snps <- 14935

### Now we assume that a similar proportion of monomorphic 
### to polymorphic loci would have passed quality filtering

### This would be represented as:
### mono_filtered <- (monomorphic/total)*filtered

### There were also 791 singletons that we filtered out in our replicate-based filtering step
### A similar portion of these were also probably true good monomorphic loci

### Adding these to our mono_filtered would be represented as:
### mono_filtered_final <- mono_filtered + (singletons/total)*791

## Now to rearrange these two equations to go with the equation from Appendix S1:
### filtered loci + (initial monomorphic loci / initial loci) * filtered loci +
###                 (initial singleton loci / initial loci) * singletons removed at technical replicate step

final_loci <- filtered_loci + (monomorphic / initial_loci) * filtered_loci + (singletons / initial_loci) * 791

## The last step to calculate total sequenced length is to multiply
## the final number of loci by their average length (as in Appendix S1)
lengths <- read.table("../../data/seq-length_vcf_loci", header = F) ### list of loci lengths
mean_loc_length <- mean(lengths$V1) ### 63.48745

L_final <- final_loci * mean_loc_length ### 551040


## Import diversity stats table produced by Stats_from_SFS_TD.R
diversity_stats <- read.table("Diversity_stats.txt", header = T)

## Undo per site standardization done already with length of only S passed in SFS
diversity_stats$pi <- diversity_stats$nsites*diversity_stats$pi
diversity_stats$theta <- diversity_stats$nsites*diversity_stats$theta

## Recalibrate each population's total length using the equation from Appendix S1:
## L * (sum of projected SFS[nsites in table] / S[filtered_snps])
diversity_stats$L<-L_final*(diversity_stats$nsites/filtered_snps)

## Standardize per site nucleotide diversity (pi) by recalibrated length
diversity_stats$pi_per_site<-diversity_stats$pi/diversity_stats$L

## Write final diversity table and use it in correlation analyses
write.table(diversity_stats, "div_per_site.txt", row.names = FALSE, quote = FALSE, sep = "\t")
