This repostory contains all of the codes used to produce simulations in the manuscript "Identifying host factors and lineage-specific fitness effects in *Mycobacterium tuberculosis* phylogenies using phylogenetic branching statistics".

We are not able to make the Malawi publicly available to protect individuals' privacy, so this repository does not contain code to reproduce Figures 2, S1, or S2. (Figure 1 is a diagram of the models.)

All of the code is run in R. You will need the following packages and their dependencies:

```{r}
install.packages(TiPS)
install.packages(ape)
install.packages(ggplot2)
install.packages(ggtree)
install.packages(cowplot)
install.packages(tidyr)
install.packages(dplyr)
```

The most important package in the list is TiPS. It provides a way to specify a compartmental model, simulate it using Gillespie/tau leaping, and then use the realized events to simulate phylogenetic trees of sampled infections using a coalescent approach. You can find a tutorial [here](https://cran.r-project.org/web/packages/TiPS/vignettes/TiPS.html) and read about the method in their [Methods in Ecology and Evolution](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14038).

Here is how to reproduce the rest of the figures in the manuscript:

Fig 3, Fig S3, Fig S4 - run nullmodel.R. This is a smaller version of the code used to produce Fig4.

Fig 4 - run nullmodel_bigsim.R (which calls nullmodel_bigsim_gettrees.R) to simulate all of the trees, then run predictionplots_nullmodel.R to produce ROC plots

Fig 5 - run varmodel_bigsim.R (which calls varmodel_bigsim_gettrees.R) to simulate the trees, then run predictionplots_varmodel.R 

FigS5, FigS6 - run varmodel.R

Fig6, FigS7, FigS8 - run pophetmodel_case1_bigsim.R, pophetmodel_case2_bigsim.R (which call pophetmodel_casex_bigsim_gettrees.R), then run predictionplots_pophetmodel.R to produce each of these figures
