
## Rscript --vanilla $pop $absolute_path
## Recalibration of L for all populations analysed using stairway plot v2
## Notation and variable names correspond to formulae presented in Appendix S1

argums <- commandArgs(trailingOnly = FALSE)[7:8]
pop <- argums[1]
basedir <- argums[2]

div_per_site <- read.table(paste0(basedir, "/results/diversity/div_per_site.txt"), header = T)
pop_row <- div_per_site[which(div_per_site$POP == pop),]

Spop <- pop_row$S
Lpop <- pop_row$L

thinned_SFS <- read.table(paste0("thinned_SFS.txt"), skip = 1, header = T)
Sthin <- sum(thinned_SFS)-thinned_SFS[1,1]
Lthinned <- (Sthin/Spop)*Lpop

print(paste0("Spop = ", round(Spop)))
print(paste0("Lpop = ", round(Lpop)))
print(paste0("Sthin = ", round(Sthin)))
print(paste0("Lpop(thinned) = ", round(Lthinned)))
