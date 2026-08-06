require(lubridate)
require(ggnewscale)
require(phytools)

getmainplot <- function(tree,taulbi=4,tauthd=5,taurels=6,tauclust=6,title='title'){

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(time = tree$tip.height,
		label = tree$tip.label)

	# the node labels have the times; extract these:
	m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
	crud2 <- data.frame(time=as.numeric(m), 
			label = names(m))
	# the node labels are super clunky, but we need to keep them to match with the tree
	crud <- rbind(crud,crud2)

	# rearrange columns with labels first:
	crud <- crud[,c(2,1)]

	# calculate LBI for the tips and the nodes:
	crud$lbi <- lbi(tree, tau=taulbi)

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'IH'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'IL'

	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]+0))
	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	p <- ggtree(tree,layout='rectangular', mrsd='2019-01-01') + theme_tree2() 
	p <- p %<+% crud

	p1 <- p + geom_tree(linewidth=0.60) +
		scale_color_manual(name='Host',
		values=c('IL'='gray80','IH'='gray27')) +
		theme(legend.position='none') +
		labs(title=title) + theme(plot.title=element_text(hjust=0.5,face='bold',size=18)) 
	
	# tree's actual x- and y-ranges (before the heatmap is appended):
        xrng <- layer_scales(p1)$x$range$range
        yrng <- layer_scales(p1)$y$range$range

        brks <- scales::pretty_breaks(n = 4)(xrng)
        brks <- brks[brks >= xrng[1] & brks <= xrng[2]]
        lbls <- as.character(as.integer(format(date_decimal(brks), '%Y')))

        # geometry for the hand-drawn axis, positioned just below the tips:
        axis_y  <- yrng[1] - 0.03 * diff(yrng)
        tick_ln <- 0.015 * diff(yrng)

        # calculate tree height:
        treeheight <- max(node.depth.edgelength(tree))

        # use cophenetic distances to calculate THD and No. of close relatives:
        x = cophenetic(tree)

        # calculate statistics for the tree (THD, LBI, No. of close relatives):

        # calculate THD from cophenetic distances:
        dat <- apply(x, 1, function(z) sum(exp(-z/tauthd)))

        # organize into a dataframe:
        dat <- as.data.frame(dat)
        colnames(dat) <- 'THD'
        dat$label = rownames(dat)
        dat <- dat[,c(2,1)]

        # calculate LBI directly from the tree:
        dat$LBI <- lbi(tree,tau=taulbi)[1:length(tree$tip.label)]

	#dat$LBI <- dat$LBI - min(dat$LBI)
	#dat$THD <- dat$THD - min(dat$THD)

        #count the number of close relatives within a threshold distance:
        dat$nrelatives <- apply(x,1,function(x) sum(x <= taurels)) - 1

        # calculate clusters using single linkage:
        hc <- hclust(as.dist(x), method='single')
        clusts <- cutree(hc, h=tauclust)
        dat$Cluster <- clusts

        # calculate cluster size:
        dat[,'Cluster Size (single)'] <- table(dat$Cluster)[dat$Cluster]

        # calculate clusters using complete linkage (identical w/ No. of close relatives?)
        hc <- hclust(as.dist(x), method='complete')
        clusts <- cutree(hc, h=tauclust)
        dat[,'Cluster Size (complete)'] <- table(clusts)[clusts]

	# save the raw values in columns:
	dat$LBIraw <- dat$LBI
	dat$THDraw <- dat$THD
	dat$nrelativesraw <- dat$nrelatives
	dat$clustcompraw <- dat$clustcompraw
	dat$clustsingraw <- dat$clustsingraw


	# standardize statistics for ease of comparison (uncomment to show raw stats):
	dat$LBI <- (dat$LBI-mean(dat$LBI))/sd(dat$LBI)
	dat$THD <- (dat$THD-mean(dat$THD))/sd(dat$THD)
	dat$nrelatives <- (dat$nrelatives-mean(dat$nrelatives))/sd(dat$nrelatives)
	dat[,'Cluster Size (complete)'] <- (dat[,'Cluster Size (complete)']-mean(dat[,'Cluster Size (complete)']))/sd(dat[,'Cluster Size (complete)'])
	dat[,'Cluster Size (single)'] <- (dat[,'Cluster Size (single)']-mean(dat[,'Cluster Size (single)']))/sd(dat[,'Cluster Size (single)'])


	# create a ggtree plot:
        #plin4 <- ggtree(tree, layout='rectangular') 
       	plin4 <- p1 + theme(legend.position='none')
 
	# rename the columns:
        colnames(dat)[4] <- "No. of close relatives"

	statedat <- as.data.frame(crud[,'state'])
	rownames(statedat) <- crud$label
	colnames(statedat) <- 'Host'

	discretefigvar <- gheatmap(plin4, statedat, offset=0, width=0.2, colnames_angle=-45,
			colnames_offset_y = -8, hjust=0.5, font.size=4,
			colnames_position='bottom', color=NA) + 
			scale_fill_discrete(name='Host',breaks=c('IH','IL'), labels=expression(italic(I[h]),italic(I[l])),
				type=c('gray27','gray80'))  + 
			theme(panel.spacing = unit(0,'pt'),legend.title =element_text(size=18),
				legend.text=element_text(size=18)) 

	#discretefigvar2 <- ggplot(discretefigvar, aes(x=c(1),y=label)) +
	#	geom_tile(aes(fill=Host)) + xlab(NULL) + ylab(NULL)
	 
	discretefigvar2 <- discretefigvar + new_scale_fill()

        # use a heatmap to visualize the statistics:
        heatfig <-  gheatmap(discretefigvar2,dat[,c('LBI',
			'THD',
                        'Cluster Size (single)',
                        'Cluster Size (complete)',
                        'No. of close relatives')],
                colnames=T, colnames_position="bottom", hjust=0.0,offset=40,
                colnames_offset_y=-3,colnames_angle=-45)+
                scale_fill_continuous(name='Value of\nstatistic\nat tips\n(standardized)',
                low='#FEFE62',high='#5D3A9B') +
                theme(plot.margin=unit(c(1,1,3,1),'cm')) +
                coord_cartesian(clip = 'off') +
                ggtitle(title) +
                theme(plot.title=element_text(hjust=0.5,size=18,face="bold")) + 
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold')) +
		theme(legend.text = element_text(size=16), legend.key.size = unit(1.0,'cm'),
			legend.title=element_text(size=18)) + 
		# --- suppress the automatic (full-width) axis from theme_tree2() ---
		theme(axis.line.x  = element_blank(),
		      axis.ticks.x = element_blank(),
		      axis.text.x  = element_blank(),
		      axis.title.x = element_blank()) +
		# axis line, restricted to the tree's x-range only:
		annotate('segment', x = xrng[1], xend = xrng[2],
		         y = axis_y, yend = axis_y, linewidth = 0.5) +
		# tick marks at each year break:
		annotate('segment', x = brks, xend = brks,
		         y = axis_y, yend = axis_y - tick_ln, linewidth = 0.5) +
		# year labels:
		annotate('text', x = brks, y = axis_y - tick_ln*2,
		         label = lbls, size = 3.5, angle=-45, vjust = 1) +
		annotate('text', x = mean(xrng), y = axis_y - tick_ln*5,
		         label = 'Year', size = 4, fontface = 'bold', vjust = 1) +
		coord_cartesian(clip = 'off')
	
	#reorder the factor levels in crud$state:
	# If we want to just look at LBI at the tips, we need to just use the first
	# ntests rows of crud:
	crud$state <- factor(crud$state, levels=c('IL','IH'))

	p1.box <- ggplot(crud[1:ntests,]) + geom_boxplot(aes(x=state,y=lbi,fill=state)) +
		scale_fill_manual(name='Host' ,
			values=c('IL'='gray80','IH'='gray27')) + 
		scale_x_discrete(labels = expression(italic(I[l]),italic(I[h]))) +
		theme_classic() + 
		ylab('LBI (raw)') + 
		xlab('Host') + 
		theme(legend.position='none', 
			axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold'))
		

	# make the tree and boxplot figures w/out the title first:
	alnd <- align_plots(heatfig,p1.box,align='v',axis='lr')
	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1, rel_heights=c(2,0.5))



        return(list(fig.p1,dat))
}

figcross <- plot_grid(getmainplot(tree3,title='Cross-sectional sampling,\nsuperspreading',taulbi=0.9,tauthd=2)[[1]],
		getmainplot(tree1,title='Cross-sectional sampling,\nhomogenous infectiousness',taulbi=2,tauthd=5)[[1]],
		nrow=1)

figlong <- plot_grid(
		getmainplot(tree4,title='Longitudinal sampling,\nsuperspreading',taulbi=.6,tauthd=6,taurels=10,tauclust=10)[[1]],
		getmainplot(tree2,title='Longitudinal sampling,\nhomogeneous infectiousness',taulbi=5,tauthd=10,taurels=20,tauclust=20)[[1]],
		nrow=1)

ggsave(figcross, file='figures/nullmodelfigs/Fig3_simtrees.png', dpi=600, height=16, width=12)

ggsave(figlong, file='figures/nullmodelfigs/FigS4_simtrees.png', dpi=600, height=16, width=12)

# Want to visualize the relationship between LBI and superspreader status
dat1 <- getmainplot(tree1,title='')[[2]]
dat2 <- getmainplot(tree2,title='')[[2]]
dat3 <- getmainplot(tree3,title='')[[2]]
dat4 <- getmainplot(tree4,title='')[[2]]


dat1$Infectiousness <- as.numeric(grepl(dat1$label, pattern='IH'))
dat1$Infectiousness <- dat1$Infectiousness + 1
dat1$Infectiousness <- c('Low', 'High')[dat1$Infectiousness]

dat2$Infectiousness <- as.numeric(grepl(dat2$label, pattern='IH'))
dat2$Infectiousness <- dat2$Infectiousness + 1
dat2$Infectiousness <- c('Low', 'High')[dat2$Infectiousness]

dat3$Infectiousness <- as.numeric(grepl(dat3$label, pattern='IH'))
dat3$Infectiousness <- dat3$Infectiousness + 1
dat3$Infectiousness <- c('Low', 'High')[dat3$Infectiousness]

dat4$Infectiousness <- as.numeric(grepl(dat4$label, pattern='IH'))
dat4$Infectiousness <- dat4$Infectiousness + 1
dat4$Infectiousness <- c('Low', 'High')[dat4$Infectiousness]

dat3$Sim <- 3
dat4$Sim <- 4

boxdat <- rbind(dat3,dat4) 

# reorder factors in boxdat:
boxdat$Sim <- factor(boxdat$Sim, levels=c(4,3))
boxdat$Infectiousness <- factor(boxdat$Infectiousness, levels=c('Low','High'))


#pdf(file='figures/nullmodelfigs/boxplots.pdf',width=8,height=8)
#par(mfrow=c(1,1))
#boxplot(LBIraw~Infectiousness*Sim,boxdat,col=c('blue','blue','darkblue','darkblue'),
#	xaxt='n',main='LBI ~ Infectiousness x Sampling scheme',xlab='',
#	ylab='Local Branching Index')
#axis(1, line=2,
#	at = c(1,2,3,4),
#	labels = c(
#		'Low\ninfectiousness+\nlongitudinal\nsampling',
#		'High\ninfectiousness+\nlongitudinal\nsampling',
#		'Low\ninfectiousness+\ncross-sectional\nsampling',
#		'High\ninfectiousness+\ncross-sectional\nsampling'),
#	tick=F, cex=0.3)
##dev.off()

#pdf(file='figures/nullmodelfigs/boxplots.pdf',width=9,height=9)
#par(mfrow=c(1,2))
#boxplot(LBI~Infectiousness,dat3,main='Cross-sectional sampling,\nsuperspreading')
#boxplot(LBI~Infectiousness,dat4,main='Longitudinal sampling,\nsuperspreading')
#dev.off()

gettbls <- function(tree){
	tbls <-	tree$edge.length[ tree$edge[,2] <= length(tree$tip.label) ]
	return(tbls)
}


#dat1 <- getmainplot(tree1)[[2]]
#dat2 <- getmainplot(tree2)[[2]]

# calculate LTT, LBI, and TBL for the empirical tree:
#lin4phi <- read.nexus('/Users/abeams/Documents/projects/TB/lineage_4_delphy/alex/lineage4_delphy.trees')
#lin4tree <- lin4phi[[200]]

lin4tree <- read.tree('/Users/abeams/Documents/projects/tblbi/data/BactDatingLin4_relaxed_fixed.tree')

lin4complete <- lin4tree

# Extract the largest clade from 175 years ago for comparison:
subtrees <- treeSlice(lin4tree, slice = 1317 - 200, orientation = "tipwards", trivial = TRUE)

tip_counts <- sapply(subtrees, function(x) length(x$tip.label))
largest_clade <- subtrees[[which.max(tip_counts)]]

lin4tree <- largest_clade

lin4ltt <- ltt.plot.coords(lin4tree)
lin4lbi <- lbi(lin4tree,4)[1:length(lin4tree$tip.label)]
lin4tbl <- gettbls(lin4tree)


# make a panel displaying LTT, distribution of LBI, and tbl distribution:

pdf(file='figures/nullmodelfigs/Fig3_stats.pdf',height=5,width=15)

par(mfrow=c(1,3),mar=c(5,5,4,2))

# 1. make LTT plots for the trees:
ltt1 <- ltt.plot.coords(tree1)
ltt2 <- ltt.plot.coords(tree2)
ltt3 <- ltt.plot.coords(tree3)
ltt4 <- ltt.plot.coords(tree4)


inds = 1:dim(ltt1)[1]
plot(ltt1[inds,'time'], ltt1[inds,'N'], type = "l",col = "#5D3A9B",lwd = 4,
	xlab = "Time",ylab = "No. of lineages",
	main = "Lineages through time",cex.main=2,cex.axis=2,cex.lab=2,cex=2,
	ylim=c(0,100))
lines(lin4ltt[,'time'],lin4ltt[,'N'],lty=1,lwd=4,col='black')
lines(ltt3[inds,'time'],ltt3[inds,'N'],lty=1,lwd=3,col='#E66100')

legend('topleft',legend=
	c('Homogeneous inf.',
	'Superspreading',
	'Empirical Lineage 4 phy.'),
	lty=c(1,1,1),
	lwd=c(3,3,3),cex=2,col=c('#5D3A9B','#E66100','black'),bty='n')

# 2. Display tbl distributions

plot(density(gettbls(tree1)),lwd=3, xlab='Terminal branch length',
	ylab='Density',main='Terminal Branch Lengths',
	cex.main=2,cex.axis=2,cex.lab=2,cex=2,lty=1,ylim=c(0,0.2),
	col='#5D3A9B',xlim=c(0,30))
lines(density(gettbls(tree3)),lwd=3, xlab='LBI',ylab='Density',	
	main='LBI distribution',lty=1,
	col='#E66100')
lines(density(gettbls(lin4tree)),lwd=3, xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
	col='black')
	
# 3. display LBI distributions:

plot(ecdf(dat1$LBIraw),xlab='LBI',ylab='Cumulative density',main='Local Branching Index',
	cex.main=2,cex.axis=2,cex.lab=2,cex=2,lty=1,ylim=c(0,1.0),xlim=c(0,20),
	verticals=T, do.points=F, lwd=3, col='#5D3A9B')
lines(ecdf(dat3$LBIraw), 
	xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
	verticals=T,do.points=F,col='#E66100',lwd=3)
lines(ecdf(lin4lbi), 
	xlab='LBI',ylab='Density',main='LBI distribution',lty=3,col='black',
	verticals=T, do.points=F,lwd=3)

dev.off()

## making figures that have boxplots for LBI~Group, and superspreader status
## annotated onto the phylogenies

getboxtreefig <- function(tree,title='add a title!',titleadjust=0.45,ptcex=1){

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(time = tree$tip.height,
		label = tree$tip.label)

	# the node labels have the times; extract these:
	m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
	crud2 <- data.frame(time=as.numeric(m), 
			label = names(m))
	# the node labels are super clunky, but we need to keep them to match with the tree
	crud <- rbind(crud,crud2)

	# rearrange columns with labels first:
	crud <- crud[,c(2,1)]

	# calculate LBI for the tips and the nodes:
	crud$lbi20 <- lbi(tree, tau=20)

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'IH'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'IL'

	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]+0))

	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	# add a column for group:
	#crud$group <- NA
	#crud$group[grep(1,crud$state)] <- 'Group 1'
	#crud$group[grep(2,crud$state)] <- 'Group 2'


	# make a ggtree object:
	p <- ggtree(tree)

	# merge the dataframe with LBI onto it:
	p <- p %<+% crud

	# plot it!
	#p1 <- p + geom_tippoint(aes(col=group)) + geom_nodepoint(aes(col=group)) + 
	#	scale_color_discrete(name='Host subpopulation',
	#	type=c('gray80','gray27'))  +
	#	theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face='bold')) +
	#	theme(legend.position='none')	
	#	#theme(legend.text=element_text(size=18),legend.title=element_text(size=16,face='bold')) +
	#	#theme(legend.position='left')

	p1 <- p + aes(col=state) + geom_tree(linewidth=0.60) +
		scale_color_discrete(name='Host',
		values=c('IL'='gray80','IH'='gray27')) +
		theme(legend.position='none') +
		labs(title=title) + theme(plot.title=element_text(hjust=0.5,face='bold',size=18))
		#theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face='bold')) +
		#theme(legend.text=element_text(size=18),legend.title=element_text(size=16,face='bold')) +
		#theme(legend.position='left')

	## make a panel plot

	p1.box <- ggplot(crud) + geom_boxplot(aes(x=state,y=lbi20,fill=state)) +
		scale_fill_manual(name='Host' ,values=c('IL'='gray80','IH'='gray27')) + 
		theme_classic() + 
		ylab('LBI') + 
		xlab('Host') + 
		theme(legend.position='none') +
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold'))
		

	# make the tree and boxplot figures w/out the title first:
	alnd <- align_plots(p1,p1.box,align='v',axis='lr')
	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1, rel_heights=c(2,0.75))

	# make the LBI boxplots an inset:
	fig.p1.inset <- ggdraw(p1) +
		draw_plot(p1.box, 0.05, .6, .2, .35)  

	return(fig.p1)
}


gettauplot <- function(tree,title='add a title!'){

	# want to add in rows for nodes with times and LBIs
	crud <- data.frame(time = tree$tip.height,
		label = tree$tip.label)

	# the node labels have the times; extract these:
	m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
	crud2 <- data.frame(time=as.numeric(m), 
			label = names(m))
	# the node labels are super clunky, but we need to keep them to match with the tree
	crud <- rbind(crud,crud2)

	# rearrange columns with labels first:
	crud <- crud[,c(2,1)]

	# calculate LBI for the tips and the nodes for lots of tau values:

	# tau vals for plotting:
	tauvals <- seq(-6,1.3, length=100)
	tauvals <- 10^tauvals

	for(tau in tauvals) crud[,paste0('lbi',tau)] <- lbi(tree, tau=tau)

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'IH'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'IL'

	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]+0))

	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	#ratios <- apply(crud[,paste0('lbi',tauvals)], 2, function(x) mean(x[crud$state=='IH'])/mean(x[crud$state=='IL'])  )
	ratios <- apply(crud[,paste0('lbi',tauvals)], 2, function(x) summary(aov(x~state,crud))[[1]][["F value"]][1]  )
	#ratios <- apply(crud[,paste0('lbi',tauvals)], 2, function(x) summary(aov(x~state,crud))[[1]][["Pr(>F)"]][1]  )



	plot(tauvals, ratios, xlab=bquote(tau),ylab='F statistic',  
		type='l',lwd=3, ylim=c(0, max(ratios)*1.5),
		cex.axis=1.5,cex.lab=1.5, main=title)
}


pdf(file='figures/nullmodelfigs/FigS3_superspreading.pdf',height=6,width=6)
gettauplot(tree3, title='Superspreading')
dev.off()

pdf(file='figures/nullmodelfigs/FigS3_homogeneous.pdf',height=6,width=6)
gettauplot(tree1, title='Homogeneous infectiousness')
dev.off()


