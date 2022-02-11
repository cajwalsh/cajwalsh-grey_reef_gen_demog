
## Rscript --vanilla
## Recalibration of L for North GBR and Indonesia joint population history modelling
## Notation and variable names correspond to formulae presented in Appendix S1

setwd("results/moments")

unthinned_SFS <- read.table("un-thinned_easySFS/fastsimcoal2/ngbr_indo_un-thinned_jointMAFpop1_0.obs", skip = 1)
sumSFS <- sum(unthinned_SFS)
Spop <- sumSFS-unthinned_SFS[1,1]
Lpop <- (sumSFS/14935)*551040

thinned_SFS <- read.table("thinned_easySFS/fastsimcoal2/ngbr_indo_thinned_jointMAFpop1_0.obs", skip = 1)
Sthin <- sum(thinned_SFS)-thinned_SFS[1,1]
Lthinned <- (Sthin/Spop)*Lpop

print(paste0("sumSFS = ", round(sumSFS)))
print(paste0("Spop = ", round(Spop)))
print(paste0("Lpop = ", round(Lpop)))
print(paste0("Sthin = ", round(Sthin)))
print(paste0("Lpop(thinned) = ", round(Lthinned)))
