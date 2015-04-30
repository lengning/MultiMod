# MultiMod
Identify genes with multiple modes


Wrapper codes to identify genes with more than one modes. Example to run the code (from command line):

Rscript MultiMod.R SCexample.csv

or

Rscript MultiMod.R SCexample.csv F 0 5 F

The input values:

The 3rd term indicates the name of the input data set. Currently the program takes csv files or tab delimited file. The input file will be treated as a tab delimited file if the suffix is not '.csv'. Rows are genes and columns are samples. Row names and column names are required in the input file.

The 4th term defines whether 0s will be excluded when estimating number of modes. Default is F.

The 5th term defines the lower limit of detection threshold (default is 0). Genes with max expression below this threshold will be removed.

The 6th term defines cutoff of number of expressed cells. Genes having number of non-zeros < cutoff will be exlcuded in the analysis. Default is 5.

The 7th term defines whether normalization is needed.  If T is specified, median-by-ratio normalization will be performed prior to PCA analysis (default is F).

Outputs:

For each gene, number of modes are estimated using Gaussian mixture model (mclust, univariate, unequal variance). 
If normalization is enabled, normalization will be performed first.
Mclust was applied on log(data+1).

XX_T_Multi_mod.csv: contains genes with more than 1 modes. 

XX_Num_clusters.csv: number of clusters for each gene.
