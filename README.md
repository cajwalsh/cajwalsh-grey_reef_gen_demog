# Code repository for: ***Genomic insights into the historical and contemporary demographics of the grey reef shark***
This repository contains most of the code needed to replicate the analyses from the following paper:

Cameron Walsh, Paolo Momigliano, Germain Boussarie, William Robbins, Lucas Bonnin, Cecile Fauvelot, Jeremy Kiszka, David Mouillot, Laurent Vigliola, and Stéphanie Manel (2022). ***Genomic insights into the historical and contemporary demographics of the grey reef shark***. Heredity 128(4), 225–235 DOI: 10.1038/s41437-022-00514-4

To replicate the analyses presented in this paper, you would first need to download the data from Genbank (Bioproject PRJNA795958) and Zenodo (https://doi.org/10.5281/zenodo.5935133), and follow the calling and filtering processes outlined in the Supplementary Information.

Most analyses in this repository are conducted in R and linked together through bash scripts, but some require [Nextflow](https://www.nextflow.io/). To run the Nextflow scripts used in this study, you will need to be working in a directory containing one of the .nf files and execute "nextflow run FILENAME.nf" on the command line (where FILENAME is the .nf file for each analysis subfolder Nextflow is used in).

To perform the genetic diversity and stairway plot analyses in this repository (with this paper's data or your own), you will first need to run the easySFSpreview.nf script (which uses [isaacovercast/easySFS](https://github.com/isaacovercast/easySFS)).

From this point on, you can check out the subfolders and their descriptions to see how the diversity, *stairway plot*, *moments*, and *NeEstimator* analyses were performed. Any external scripts or software programs needed to perform these analyses are linked in the descriptions of each analysis subfolder. The absolute directory/full file path containing this repository on your local machine must be changed in the first few lines of some scripts.  The analyses contained in each folder will write its results into a subfolder of the same name in the main results folder.
