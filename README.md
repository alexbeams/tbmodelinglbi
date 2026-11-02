This repository contains all of the codes needed to produce simulations and analysis in the manuscript "Identifying host and pathogen factors affecting transmission in Mycobacterium tuberculosis phylogenies".

To protect individuals' privacy, sequence names in the metadata and phylogenetic trees have been renamed with numeric labels. 

The metadata file (data_analysis/tblbi_metadata_anonymized.csv) contains columns for covariates corresponding to sampled infections corresponding to tips in the phylogenetic trees (data_analysis/BactDatingLin4_strict_fixed_anonymized.tree and data/analysis_BactDatingLin4_relaxed_fixed_anonymized.tree). The first column of the metadata csv file contains the corresponding tip in the phylogenetic trees, the second column contains the year that the infection was sampled, and all other columns correspond to the covariates summarized in Table S1 of the Electronic Supplementary Material associated with our manuscript. The phylogenetic trees are the same ones described in the Methods section of our manuscript. The branch lengths in the trees represent time as measured in years.

All of the code is run in R version 4.4.3. The file environment.yml contains the settings for the mamba environment we use to run this analysis. The main packages (and their dependencies) you will need are the following:

```{r}
install.packages(TiPS)
install.packages(ape)
install.packages(ggplot2)
install.packages(ggtree)
install.packages(cowplot)
install.packages(tidyr)
install.packages(dplyr)
```

The most important package in the list is TiPS. It provides a way to specify a compartmental model, simulate it using Gillespie/tau leaping, and then use the realized events to simulate phylogenetic trees of sampled infections using a coalescent approach. You can find a tutorial [here](https://cran.r-project.org/web/packages/TiPS/vignettes/TiPS.html) and read about the method in their paper in [Methods in Ecology and Evolution](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14038).

Here is how to reproduce the main figures in the manuscript:

Fig 1, Fig S2, Fig S3 - run getFigMalawi.R

Fig 2, Fig S4, Fig S6 - run nullmodel.R and then nullmodelplot.R.

Fig S5 - run nullmodel_bigsim.R (which calls nullmodel_bigsim_gettrees.R) to simulate all of the trees, then run predictionplots_nullmodel.R to produce ROC plots

Fig 3, Fig S7,  - run varmodel_bigsim.R (which calls varmodel_bigsim_gettrees.R) to simulate the trees, then run predictionplots_varmodel.R 

FigS5, FigS6 - run varmodel.R

Fig4, FigS8, Fig S9 - run pophetmodel_case1_bigsim.R, pophetmodel_case2_bigsim.R (which call pophetmodel_casex_bigsim_gettrees.R), then run predictionplots_pophetmodel.R to produce each of these figures

Fig 5, S10, S11, S12 - run malawi_lbi.R, which is located in the data_analysis folder
