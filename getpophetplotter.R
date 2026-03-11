getmainplot <- function(tree,taulbi=4,tauthd=5,taurels=6,tauclust=6,title='title'){

	# want to add in rows for nodes with times and LBIs
	#crud <- data.frame(time = tree$tip.height,
	#	label = tree$tip.label)
	crud <- data.frame(label = tree$tip.label)


	# the node labels have the times; extract these:
	#m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
	m<- sapply(tree$node.label, function(z) substr(z, regexpr('string=',z)[1]+11, regexpr('string=',z)[1]+13 ) ) 
	#crud2 <- data.frame(time=as.numeric(m), 
	#		label = names(m))
	crud2 <- data.frame(label = names(m))
	

	# the node labels are super clunky, but we need to keep them to match with the tree
	crud <- rbind(crud,crud2)

	# rearrange columns with labels first:
	#crud <- crud[,c(2,1)]

	# calculate LBI for the tips and the nodes:
	crud$lbi <- lbi(tree, tau=taulbi)

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IL1',crud$label[1:ntests]),'state'] <- 'Group 1'
	crud[grep('IH2',crud$label[1:ntests]),'state'] <- 'Group 2'
	crud[grep('IL2',crud$label[1:ntests]),'state'] <- 'Group 2'


	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+3, regexpr(".[+]=",z)[1]-0))
	nodenms[grep('IH1',nodenms)] <- 'Group 1'
	nodenms[grep('IL1',nodenms)] <- 'Group 1'
	nodenms[grep('IH2',nodenms)] <- 'Group 2'
	nodenms[grep('IL2',nodenms)] <- 'Group 2'


	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	p <- ggtree(tree,layout='rectangular') %<+% crud

	p1 <- p + aes(col=state) + geom_tree(linewidth=0.60) +
		scale_color_manual(name='Host',
		values=c('Group 1'='#E1BE6A','Group 2'='#40B0A6')) +
		theme(legend.position='none') +
		labs(title=title) + theme(plot.title=element_text(hjust=0.5,face='bold',size=18))
		#theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face='bold')) +
		#theme(legend.text=element_text(size=18),legend.title=element_text(size=16,face='bold')) +
		#theme(legend.position='left')


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

	# save the raw values in columns:
	dat$LBIraw <- dat$LBI

	# standardize statistics for ease of comparison (uncomment to show raw stats):
	dat$LBI <- (dat$LBI-mean(dat$LBI))/sd(dat$LBI)

	# create a ggtree plot:
       	plin4 <- p1
 

	lbidat <- as.data.frame(dat[,'LBIraw'])
	rownames(lbidat) <- rownames(dat)
	colnames(lbidat) <- 'LBI'

        # use a heatmap to visualize the statistics:
        heatfig <-  gheatmap(plin4,lbidat,
                colnames=T, colnames_position="bottom", hjust=0.0,
                colnames_offset_y=-3,colnames_angle=-45,width=0.1)+
                scale_fill_continuous(name='Value of\nLBI\nat tips\n(raw)',
                low='#FEFE62',high='#5D3A9B') +
                theme(plot.margin=unit(c(1,1,3,1),'cm')) +
                coord_cartesian(clip = 'off') +
                ggtitle(title) +
                theme(plot.title=element_text(hjust=0.5,size=18,face="bold")) + 
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold')) +
		theme(legend.text = element_text(size=16), legend.key.size = unit(1.0,'cm'),
			legend.title=element_text(size=18)) + 
		guides(
		 color=guide_legend(order = 1, override.aes=list(linewidth=2)),
		 fill = guide_colorbar(order=2)
		)
	
	#reorder the factor levels in crud$state:
	# If we want to just look at LBI at the tips, we need to just use the first
	# ntests rows of crud:
	crud$state <- factor(crud$state, levels=c('Group 1','Group 2'))

	p1.box <- ggplot(crud[1:ntests,]) + geom_boxplot(aes(x=state,y=lbi,fill=state)) +
		scale_fill_manual(name='Host' ,values=c('Group 1'='#E1BE6A','Group 2'='#40B0A6')) + 
		theme_classic() + 
		ylab('LBI (raw)') + 
		xlab('Subpopulation') + 
		theme(legend.position='none') +
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold'))
	

	# make the tree and boxplot figures w/out the title first:
	alnd <- align_plots(heatfig,p1.box,align='v',axis='lr')
	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1, rel_heights=c(2,0.5))



        return(list(fig.p1,dat))
}


