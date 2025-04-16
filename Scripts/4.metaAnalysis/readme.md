# Pipeline Description 

Note that the files in this directory correspond to two separate Snakemake pipelines. 
In the first pipeline (Snakefile is named Snakefile_data_generation), we generate
the 20 studies of 6,998 individuals. 

In the second Snakemake pipeline (Snakefile is just called Snakefile), we perform 
the analysis using external covariance similar to what is done in the other directories. 

As with the other directories, the two Snakefiles must be placed two levels higher 
in the directory ../../