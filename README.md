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

Simulations needed to reproduce figures are included in the sims folder of this repository. There are large sections of code that are currently commented-out to prevent-rerunning simulations, meaning that figures are reproduced from the same simulations that are displayed in our paper. Below, in the instructions to reproduce the main figures in the manuscript, we point to the relevant files where users can un-comment code to rerun new simulations.

Fig 1, Fig S2 - run getFigMalawi.R

Fig S3 - run getFigLBIvsTHD.R

Fig 2, Fig S4, Fig S6 - run nullmodel.R and then nullmodelplot.R.

Fig S5 - run nullmodel_bigsim.R (which calls nullmodel_bigsim_gettrees.R) to simulate all of the trees, then run predictionplots_nullmodel.R to produce ROC plots. In nullmodel_bigsim.R and prediciton_plots_nullmodel.R, there are section that run the simulations which are currently commented out, meaning that figures are reproduced from saved simulations.

Fig 3, data for Tables S4 and S5  - run varmodel_bigsim.R (which calls varmodel_bigsim_gettrees.R) to simulate the trees, then run predictionplots_varmodel.R (and in this last file, there is a section that runs simulations which can be commented out to regenerate figures from saved simulations.)

Fig S7 - run varmodel_phylowave.R

Fig4, FigS8, Fig S9 - run pophetmodel_case1_bigsim.R, pophetmodel_case2_bigsim.R (which call pophetmodel_casex_bigsim_gettrees.R), then run predictionplots_pophetmodel.R to produce each of these figures. In pophetmodel_case1_bigsim.R and popgetmodel_case2_bigsim.R, there are sections of code commented out to prevent re-running simulations, in order to reproduce figures from the previously run simulations used for the manuscript).

Fig 5, S10, S11, S12 - run data_analysis/malawi_lbi.R
