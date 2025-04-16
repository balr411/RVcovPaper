# Pipeline Descriptions

The pipelines and scripts required to reproduce analyses in the paper are available 
in the Scripts/ directory. 

- The directory Scripts/1.mainAnalysis/ contain all of 
the main analysis using UKBB reference panels and the ten UKBB traits. 
- The directory Scripts/2.InPSYghtAnalysis/ contains all of the analyses done 
using the InPSYght reference panel
- The directory Scripts/3.SKAT_weightedBurden/ contains all of the analyses done
using SKAT and weighted burden tests
- The directory Scripts/4.metaAnalysis/ contains all of the analyses done for the 
meta-analysis scenario
- The directory Scripts/5.timeComparison/ conrains all of the analyses done comparing
the one-stage and two-stage approaches.

Note that for each of the Snakefiles in the directories, to run them as-is, they 
must be moved to this directory (../../). More information on Snakemake can be found
at https://snakemake.readthedocs.io/en/stable/. To run this pipelines with the 
corresponding Snakefile moved to this directory, use the following command:

``` bash 
snakemake --profile simpleSlurm 
```
