# Pipeline Description 

Note that the files in this directory correspond to two separate Snakemake pipelines. 
In the first pipeline (Snakefile is named Snakefile_data_generation), we generate
the 20 studies of 6,998 individuals. 

After completing this pipeline, single-variant meta-analysis needs to be completed
in order to get the pooled allele frequencies for each variant. To do this, use 
the script metaAnalysis_summaryFiles.R to create summary files for input into 
raremetal. Then use the script metaAnalysis_singleVariantMetaAnalysis.R to call 
raremetal on each of the created summary files. Note that these pooled allele 
frequencies are needed to create the group files further downstream.

After doing this, due to issues with the {study} wildcard in the next Snakefile, 
you will need to run the scripts metaAnalysis_imp0SummaryFiles.R and metaAnalysis_all0CovFiles.R
in order to output the covfiles necessary to run RAREMETAL.

In the second Snakemake pipeline (Snakefile is just called Snakefile), we perform 
the analysis using external covariance similar to what is done in the other directories. 

As with the other directories, the two Snakefiles must be placed two levels higher 
in the directory ../../