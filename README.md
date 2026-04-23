This repostory contains all of the codes used to produce simulations in the manuscript "Identifying host factors and lineage-specific fitness effects in \textit{Mycobacterium tuberculosis} phylogenies using phylogenetic branching statistics".

We are not making the Malawi publicly available at this time to protect individuals' privacy.

Here is how to reproduce the rest of the figures in the manuscript:


Fig 3, Fig S3, Fig S4 - nullmodel.R

Fig 4 - run nullmodel_bigsim.R (which calls nullmodel_bigsim_gettrees.R) to simulate all of the trees, then run predictionplots_nullmodel.R to produce ROC plots

Fig 5 - run varmodel_bigsim.R (which calls varmodel_bigsim_gettrees.R) to simulate the trees, then run predictionplots_varmodel.R 

FigS5, FigS6 - run varmodel.R

Fig6, FigS7, FigS8 - run pophetmodel_case1_bigsim.R, pophetmodel_case2_bigsim.R (which call pophetmodel_casex_bigsim_gettrees.R), then run predictionplots_pophetmodel.R to produce each of these figures
