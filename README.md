## Scripts belonging to manuscript: 
## The feedback between selection and demography shapes genomic diversity during coevolution
2018.11.23

cas.retel@eawag.ch

github.com/RetelC/TempDynamics_HostVirCoevol

- downsampleBam.sh is a standalone shell script that calls
  samtools view -s to randomly sample reads from a .bam

The files below contain functions, should be available and invoked using source() for analysis
scripts to run: 
- gg_multiplot.R
- round_10e3.R
- fs_syncCalculations.R
- fs_syncFiltering.R
- plotAfs.R
- writeReadSync.R



The four files below are standalone R scripts (invoked using Rscript on the command line)
used to manipulate text files in the .sync format (for explanation of the format, see 
https://sourceforge.net/p/popoolation/wiki/Main/). Some of them depend on functions specified 
below (invoked using R::source()): 
- syncFilterDepth.R  ## substitutes entries by NA if they're outside a given confidence interval
- syncToAaf.R ## substitutes allelic counts at every position for the proportion equal 
    to reference base (= third column .sync)
- syncToDafs_bash.R  ## same, but outputs proportion of most abundant non-reference base
- findVariantAlleles.R  ## from input .aaf file (see syncToAaf.R), returns only rows where 
    absolute difference between minimum and maximum observed frequency is larger than a 
    user-specified value. 



The scripts below are R scripts used to filter relatively small (i.e. a few thousand positions; 
easily managed locally in R) .sync files, containing only information on polymorphisms that show 
variation in allele frequency values through time for sequencing artefacts (see manuscript Methods). 
They were meant to run locally, to have more control of intermediate steps and output. Plotting and 
file-writing wrappers are also included: 
- D00-99_nc64a_aaffiltering.R  ## filters pre-filtered host allele frequency dataset
- D12-99_pbcv1_aaffiltering.R  ## filters pre-filtered virus allele frequency dataset
- D00-99_nc64a_multipos.R	 ## combines genomic datasets of three replicates into one object, 
    checks for overlap, adds annotation information. 
- D12-99_pbcv1_multipos.R  ## same for virus
- D00-99_process_phenotypes.R  ## processes pairwise phenotypic comparisons

