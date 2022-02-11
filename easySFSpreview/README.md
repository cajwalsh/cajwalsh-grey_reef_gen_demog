This script uses the [easySFS](https://github.com/isaacovercast/easySFS) python script written by Isaac Overcast, which must be downloaded and used as a command in bash (as easySFS.py).

As mentioned on the main page, to run this Nextflow script (prior to doing the genetic diversity and stairway plot analyses), change your working directory to this easySFSpreview folder, and execute "nextflow run easySFSpreview.nf" on the command line.

With the output of easySFSpreview.nf, which is written in the main (initially empty) results folder, make .csv files with the desired projections using the required formatting in its results subfolder. These include:

long rows for unthinned SFSs for diversity analyses called "all_sfs_projs.csv" pop1,pop2,pop3,... \n
xx,yy,zz,... \n

and long columns for thinned SFSs for stairway plot called "sp_sfs_projs.csv"
pop1,xx \n
pop2,yy \n
pop3,zz \n
...

To replicate these analyses, you will need to use the projections we used from Table S2 in the Supplementary Information. For use with other data, see the documentation of [easySFS](https://github.com/isaacovercast/easySFS) for more information about projection.
