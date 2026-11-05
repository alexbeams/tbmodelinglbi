rm(list=ls())


source('lbi.R')

require(ape)
require(ggtree)
require(scales)
require(ggplot2)
require(cowplot)
require(lubridate)

lin4phy <- read.tree('data_analysis/BactDatingLin4_relaxed_fixed_anonymized.tree')

getmainplot <- function(tree=lin4phy, taulbi=4, tauthd=5, taurels=6, tauclust=6){
	treeheight <- max(node.depth.edgelength(tree))
	x = cophenetic(tree)

	dat <- apply(x, 1, function(z) sum(exp(-z/tauthd)))
	dat <- as.data.frame(dat)
	colnames(dat) <- 'THD'
	dat$label = rownames(dat)
	dat <- dat[,c(2,1)]
	dat$LBI <- lbi(tree,tau=taulbi)[1:length(tree$tip.label)]
	dat$nrelatives <- apply(x,1,function(x) sum(x <= taurels)) - 1
	hc <- hclust(as.dist(x), method='single')
	clusts <- cutree(hc, h=tauclust)
	dat$Cluster <- clusts
	dat[,'Cluster\nSize\n(single)'] <- table(dat$Cluster)[dat$Cluster]
	hc <- hclust(as.dist(x), method='complete')
	clusts <- cutree(hc, h=tauclust)
	dat[,'Cluster\nSize\n(complete)'] <- table(clusts)[clusts]

	# Rescale statistics so they are on the same scale:
	dat$LBI <- dat$LBI - min(dat$LBI)
	dat$LBI <- dat$LBI/max(dat$LBI)

	dat$THD <- dat$THD - min(dat$THD)
	dat$THD <- dat$THD/max(dat$THD)

	dat[,'Cluster\nSize\n(single)'] <- dat[,'Cluster\nSize\n(single)'] - min(dat[,'Cluster\nSize\n(single)'] ) 
	dat[,'Cluster\nSize\n(single)'] <- dat[,'Cluster\nSize\n(single)'] / max(dat[,'Cluster\nSize\n(single)'] ) 

	dat[,'Cluster\nSize\n(complete)'] <- dat[,'Cluster\nSize\n(complete)'] - min(dat[,'Cluster\nSize\n(complete)'] ) 
	dat[,'Cluster\nSize\n(complete)'] <- dat[,'Cluster\nSize\n(complete)'] / max(dat[,'Cluster\nSize\n(complete)'] ) 

	dat$nrelatives <- dat$nrelatives - min(dat$nrelatives)
	dat$nrelatives <- dat$nrelatives/max(dat$nrelatives)

	# Plot the tree:
	plin4 <- ggtree(tree, layout='rectangular', mrsd='2019-01-01') + theme_tree2()

	# tree's actual x- and y-ranges (before the heatmap is appended):
	xrng <- layer_scales(plin4)$x$range$range
	yrng <- layer_scales(plin4)$y$range$range

	brks <- scales::pretty_breaks(n = 8)(xrng)
	brks <- brks[brks >= xrng[1] & brks <= xrng[2]]
	lbls <- as.character(as.integer(format(date_decimal(brks), '%Y')))

	colnames(dat)[4] <- "No. of\nclose\nrelatives"

	# geometry for the hand-drawn axis, positioned just below the tips:
	axis_y  <- yrng[1] - 0.03 * diff(yrng)
	tick_ln <- 0.015 * diff(yrng)

	heatfig <- gheatmap(plin4, dat[,c('LBI','THD',
			'Cluster\nSize\n(single)',
			'Cluster\nSize\n(complete)',
			'No. of\nclose\nrelatives')],
		colnames=T, colnames_position="bottom", hjust=0.5,
		colnames_offset_y=-30,
		color=NA) +
		scale_fill_continuous(name='Value of\nscaled\nstatistic',
	        low='#FEFE62', high='#5D3A9B') +
		theme(
			axis.line.x  = element_blank(),
			axis.ticks.x = element_blank(),
			axis.text.x  = element_blank(),
			axis.title.x = element_blank(),
			plot.margin  = unit(c(1,1,4,1),'cm')
		) +
		# axis line, restricted to the tree's x-range only:
		annotate('segment', x = xrng[1], xend = xrng[2],
		         y = axis_y, yend = axis_y, linewidth = 0.5) +
		# tick marks at each year break:
		annotate('segment', x = brks, xend = brks,
		         y = axis_y, yend = axis_y - tick_ln, linewidth = 0.5) +
		# year labels:
		annotate('text', x = brks, y = axis_y - tick_ln*2.5,
		         label = lbls, size = 3.5, vjust = 1) +
		annotate('text', x = mean(xrng), y = axis_y - tick_ln*5,
		         label = 'Year', size = 4, fontface = 'bold', vjust = 1) +
		coord_cartesian(clip = 'off') +
		ggtitle(paste0('Tuberculosis lineage 4, Blantyre, Malawi \n(height = ',round(treeheight,2), ' yrs)')) +
		theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))

	return(heatfig)
}


lin4treeplot_main <- getmainplot()

ggsave(lin4treeplot_main, file='figures/Fig1.png',
	height=8,width=12,dpi=300)



##############################################
# plot a tree using a strict clock:
##############################################

lin4phy2 <- read.tree('data_analysis/BactDatingLin4_strict_fixed_anonymized.tree')

lin4treeplot_supp <- getmainplot(tree=lin4phy2)

ggsave(lin4treeplot_supp, file='figures/FigS2.png',
	height=8,width=12,dpi=600)



