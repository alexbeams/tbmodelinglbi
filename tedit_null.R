rm(list=ls())

source("compile_variant_and_superspreading_model.R")

# 2 useful seeds:
# with Ieq <- 0.01, use 379920
# with Ieq <- 0.001, use 509262

#############
#############

ntests <- 500


time <- c(-125,50)

# when should we start to sample?
tm1 <- 45
tms.lin4 <- runif(ntests,min=tm1,max=time[2])
tms.lin4 <- tms.lin4[order(tms.lin4)]

# what is the total population size?
Npop <- 2.0 * 10^5

gettheta <- function(x,fractrans){
	#of those with active TB, what fraction are superspreaders?
	pH <- 0.1

	# infectiousness multipliers
	f <- fractrans # fraction of transmission from superspreaders

	# parameterize in terms of equilibrium values of S, L, I
	Seq <- 2/3 * Npop
	Ieq <- x * Npop
	Leq <- Npop-Seq-Ieq

	# what proportion of newly infected immediately become infectious?
	# 5% develop active TB in the first two years; use this:
	pfast <- 0.05

	# Average duration of untreated pulmonary TB, historically, is ~ 1- 3 years:
	gamma <- 1/1.5

	# The annual risk of developing active TB is 1e-4 - 2e-4 per year:
	sigma <- 1e-4
	 
	# transmission rates (beta=zeta*c*theta):
	beta <- (gamma*Ieq-sigma*Leq)/pfast/Seq/Ieq

	# determine rho:
	rho <- (beta*Seq*Ieq - gamma*Ieq)/Leq

	# susceptibility 
	zeta <- 1

	# infectiousness
	thetaL <- (1-f)/(1-pH)
	thetaH <- f/pH

	# increased infectiousness of variant:
	deltathetaL <- 0.25 * thetaL
	deltathetaH <- 0.25 * thetaH

	# probability that latent infection mutates into variant strain:
	pmu <- 0.000

	# contact rates within/across groups:
	c <- beta

	theta <- list(
		zeta = zeta,
		thetaH = thetaH,
		thetaL = thetaL,
		deltathetaH = deltathetaH,
		deltathetaL = deltathetaL,	
		c = c,
		gamma = gamma,
		sigma = sigma,
		rho = rho,
		pmu = pmu,
		pH = pH,
		pfast=pfast
	)
	return(theta)
}

# 2 parametrizations: low and high prevalence:
#theta1 <- gettheta(0.001) # low prevalence
#theta2 <- gettheta(0.01)  # high prevalence

theta1 <- gettheta(0.01,0.1)
theta2 <- gettheta(0.01,0.9)


# set the timestep size:
dT <- 1

# specify initial states:
initialStates <- c(Npop,0,1,0,0,0,0)
names(initialStates) <- c('S','L','IH','IL','M','JH','JL') 

out1 <-  sir_simu(
   paramValues = as.list(theta1),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100,
   seed=280361)

out2 <-  sir_simu(
   paramValues = as.list(theta2),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100,
   seed=280361)

# use a smaller dataframe for plotting (just sample rows)
traj1 <- out1$traj
plottraj1 <- traj1[c(1:1500,sample(1:dim(traj1)[1],1000)),]
plottraj1 <- plottraj1[order(plottraj1$Time),]
# transform time back to date for plotting:
plottraj1$date <- as.Date(plottraj1$Time * 365)

traj2 <- out2$traj
plottraj2 <- traj2[c(1:1500,sample(1:dim(traj2)[1],1000)),]
plottraj2 <- plottraj2[order(plottraj2$Time),]
# transform time back to date for plotting:
plottraj2$date <- as.Date(plottraj2$Time * 365)



# plot the trajectories of out1 and out2:

pdf(file='figures/nullmodelfigs/nullmodeltrajs.pdf',height=5,width=9)

par(mfrow=c(1,2))

# low prevalence:
plot(log10(L)~date,plottraj1,type='l',col='#005AB5',main='Homogeneous infectiousness',lwd=2.5,
	ylab=bquote(Log[10](.('No. of infections'))),xlab='Date',lty='dotted', ylim=c(0,6))

lines(log10(IL)~date,plottraj1,type='l',col='#005AB5',main='Active TB',lwd=2.5,
	ylab='No. of infections',xlab='Date',lty='dashed')
lines(log10(IH)~date,plottraj1,type='l',col='#005AA0',lwd=2.5)
legend('top',col=c('#005AB5','#005AA0'),
	lty=c(3,2,1),lwd=2.5,legend=c(bquote(L),bquote(I[1]),bquote(I[2])),
	cex=1.3)

# high prevalence:
plot(log10(L)~date,plottraj2,type='l',col='#005AB5',main='Superspreading',lwd=2.5,
	ylab=bquote(Log[10](.('No. of infections'))),xlab='Date',lty='dotted', ylim=c(0,6))

lines(log10(IL)~date,plottraj2,type='l',col='#005AB5',main='Active TB',lwd=2.5,
	ylab='No. of infections',xlab='Date',lty='dashed')
lines(log10(IH)~date,plottraj2,type='l',col='#005AA0',lwd=2.5)
#legend('topleft',col=c('#005AB5','#005AA0'),
#	lty=c(3,2,1),lwd=2.5,legend=c(bquote(L),bquote(I[1]),bquote(I[2])),
#	cex=1.3,bty='n')
#

dev.off()

# use tms.lin4 to generate simulated trees with the same height (approx.) as the
#	empirical lineage 4 tree

pH <- 0.1

getsimtree <- function(tms,output){
	simulate_tree(
	simuResults=output,
	dates=c(tms),
	deme=c('IH','IL','L','JH','JL','M'),
	sampled=c(
		IH=pH,
		IL=(1-pH),
		JH=0,
		JL=0),
	root = 'IH',
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)

}

#fulltree <- simulate_tree(
#	simuResults=out1,
#	deme=c('IH','IL','L','JH','JL','M'),
#	sampled=c(
#		IH=pH,
#		IL=(1-pH),
#		JH=0,
#		JL=0),
#	root = 'IH',
#	nTrials=50,
#	resampling=FALSE,
#	addInfos = TRUE,
#	isFullTrajectory=TRUE)
#


tree1 <- getsimtree(tms.lin4,out1)

tms2.lin4 = runif(-100,max(tms.lin4),n=ntests)
tree2 <- getsimtree(tms2.lin4,out1)

tree3 <- getsimtree(tms.lin4,out2)
tree4 <- getsimtree(tms2.lin4,out2)

#par(mfrow=c(1,3))

### variant @ 1%
#tree <- tree1
#
#plot(tree, type='fan',show.tip.label=F)
#tips_cols <- ifelse(grepl(x=tree$tip.label,pattern="I"),"#005AB5","#DC3220")
#nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="I"),"#005AB5","#DC3220")
#tiplabels(pch=20,col=tips_cols)
#nodelabels(pch=20,col=nodes_cols)
#
#x <- cophenetic(tree)
#tau <- mean(x)*.01
#dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
#dat <- as.data.frame(dat)
#colnames(dat) <- 'thd'
#dat$label = rownames(dat)
#dat <- dat[,c(2,1)]
#dat$var <- apply(x,1,var)
#dat$min <- apply(x,1,function(x) min(x[x>0]) )
#source('lbi.R')
#
#dat$lbi <- lbi(tree,tau=50)[1:length(tree$tip.label)]
#
## now add in states:
#dat$state <- 0
#dat[grep('I',tree$tip.label),'state'] <- 'I'
#dat[grep('J',tree$tip.label),'state'] <- 'J'


# load the LBI function:
source('lbi.R')


## make plots of each tree
#p1 <- ggtree(tree1,layout='circular')
#p2 <- ggtree(tree2,layout='circular')
#
#p1 <- p1 + labs(title='') + theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))
#p2 <- p2 + labs(title='') + theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))
#
#ptrees <- plot_grid(p1,p2,nrow=1)
#
### mimic the data figure:

# This function will make the heatmap plots:

# This function will make the heatmap plots:
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

	p <- ggtree(tree,layout='rectangular') %<+% crud

	p1 <- p + aes(col=state) + geom_tree(linewidth=0.60) +
		scale_color_manual(name='Host',
		values=c('IL'='gray80','IH'='gray27')) +
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
       	plin4 <- p1
 
	# rename the columns:
        colnames(dat)[4] <- "No. of close relatives"

        # use a heatmap to visualize the statistics:
        heatfig <-  gheatmap(plin4,dat[,c('LBI',
			'THD',
                        'Cluster Size (single)',
                        'Cluster Size (complete)',
                        'No. of close relatives')],
                colnames=T, colnames_position="bottom", hjust=0.0,
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
		guides(color=guide_legend(override.aes=list(linewidth=2)))
	
	#reorder the factor levels in crud$state:
	# If we want to just look at LBI at the tips, we need to just use the first
	# ntests rows of crud:
	crud$state <- factor(crud$state, levels=c('IL','IH'))

	p1.box <- ggplot(crud[1:ntests,]) + geom_boxplot(aes(x=state,y=lbi,fill=state)) +
		scale_fill_manual(name='Host' ,values=c('IL'='gray80','IH'='gray27')) + 
		theme_classic() + 
		ylab('LBI (raw)') + 
		xlab('Host') + 
		theme(legend.position='none') +
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold'))
		

	# make the tree and boxplot figures w/out the title first:
	alnd <- align_plots(heatfig,p1.box,align='v',axis='lr')
	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1, rel_heights=c(2,0.5))



        return(list(fig.p1,dat))
}


figcross <- plot_grid(getmainplot(tree3,title='Cross-sectional sampling,\nsuperspreading',taulbi=1,tauthd=1)[[1]],
		getmainplot(tree1,title='Cross-sectional sampling,\nhomogenous infectiousness')[[1]],
		nrow=1)

figlong <- plot_grid(
		getmainplot(tree4,title='Longitudinal sampling,\nsuperspreading',taulbi=1,tauthd=1)[[1]],
		getmainplot(tree2,title='Longitudinal sampling,\nhomogeneous infectiousness')[[1]],
		nrow=1)

ggsave(figcross, file='figures/nullmodelfigs/nullmodeltrees.png', dpi=600, height=16, width=12)

ggsave(figlong, file='figures/nullmodelfigs/nullmodeltrees_long.png', dpi=600, height=16, width=12)

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
par(mfrow=c(1,1))
boxplot(LBIraw~Infectiousness*Sim,boxdat,col=c('blue','blue','darkblue','darkblue'),
	xaxt='n',main='LBI ~ Infectiousness x Sampling scheme',xlab='',
	ylab='Local Branching Index')
axis(1, line=2,
	at = c(1,2,3,4),
	labels = c(
		'Low\ninfectiousness+\nlongitudinal\nsampling',
		'High\ninfectiousness+\nlongitudinal\nsampling',
		'Low\ninfectiousness+\ncross-sectional\nsampling',
		'High\ninfectiousness+\ncross-sectional\nsampling'),
	tick=F, cex=0.3)
#dev.off()

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
lin4phi <- read.nexus('/Users/abeams/Documents/projects/TB/lineage_4_delphy/alex/lineage4_delphy.trees')
lin4tree <- lin4phi[[200]]

lin4ltt <- ltt.plot.coords(lin4tree)
lin4lbi <- lbi(lin4tree,4)[1:length(lin4tree$tip.label)]
lin4tbl <- gettbls(lin4tree)


# make a panel displaying LTT, distribution of LBI, and tbl distribution:

pdf(file='figures/nullmodelfigs/nullmodelstats.pdf',height=5,width=15)

par(mfrow=c(1,3),mar=c(5,5,4,2))

# 1. make LTT plots for the trees:
ltt1 <- ltt.plot.coords(tree1)
ltt2 <- ltt.plot.coords(tree2)
ltt3 <- ltt.plot.coords(tree3)
ltt4 <- ltt.plot.coords(tree4)


inds = 1:dim(ltt1)[1]
plot(ltt1[inds,'time'], ltt1[inds,'N'], type = "l",col = "#40B0A6",lwd = 1,
	xlab = "Time",ylab = "No. of lineages",
	main = "Lineages through time",cex.main=2,cex.axis=2,cex.lab=2,cex=2,
	ylim=c(0,500))
lines(lin4ltt[,'time'],lin4ltt[,'N'],lty=1,lwd=4,col='#DC3220')
#lines(ltt2[inds,'time'],ltt2[inds,'N'],lty=1,lwd=1,col='#E1BE6A')
lines(ltt3[inds,'time'],ltt3[inds,'N'],lty=1,lwd=3,col='#E1BE6A')
#lines(ltt4[inds,'time'],ltt4[inds,'N'],lty=1,lwd=3,col='#E1BE6A')

legend('left',legend=
	c('Homogeneous inf.',
	'Superspreading',
	'Empirical Lineage 4 phy.'),
	lty=c(1,1,1),
	lwd=c(3,3,3),cex=2,col=c('#40B0A6','#E1BE6A','#DC3220'),bty='n')

# 2. Display tbl distributions

plot(density(gettbls(tree1)),lwd=3, xlab='Terminal branch length',
	ylab='Density',main='Terminal Branch Lengths',
	cex.main=2,cex.axis=2,cex.lab=2,cex=2,lty=1,ylim=c(0,0.4),
	col='#40B0A6',xlim=c(0,30))
#lines(density(gettbls(tree2)),lwd=1, xlab='LBI',ylab='Density',	
#	main='LBI distribution',lty=1,
#	col='#E1BE6A')
lines(density(gettbls(tree3)),lwd=3, xlab='LBI',ylab='Density',	
	main='LBI distribution',lty=1,
	col='#E1BE6A')
#lines(density(gettbls(tree4)),lwd=3, xlab='LBI',ylab='Density',	
#	main='LBI distribution',lty=1,
#	col='#E1BE6A')
lines(density(gettbls(lin4tree)),lwd=3, xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
	col='#DC3220')
	
# 3. display LBI distributions:

plot(ecdf(dat1$LBIraw),xlab='LBI',ylab='Cumulative density',main='Local Branching Index',
	cex.main=2,cex.axis=2,cex.lab=2,cex=2,lty=1,ylim=c(0,1.0),xlim=c(0,30),
	verticals=T, do.points=F, lwd=3, col='#40B0A6')
#lines(ecdf(dat1$LBIraw), 
#	xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
#	verticals=T,do.points=F,col='#40B0A6',lwd=1)
lines(ecdf(dat3$LBIraw), 
	xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
	verticals=T,do.points=F,col='#E1BE6A',lwd=3)
#lines(ecdf(dat4$LBIraw), 
#	xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
#	verticals=T,do.points=F,col='#E1BE6A',lwd=3)
lines(ecdf(lin4lbi), 
	xlab='LBI',ylab='Density',main='LBI distribution',lty=3,col='#DC3220',
	verticals=T, do.points=F,lwd=3)


## 4. Display LBI ~ infectiousness for cross sect'l and long'tdl sampling:
#boxplot(LBI~Infectiousness*Sim,boxdat,col=c('blue','blue','darkblue','darkblue'),
#	xaxt='n',main='LBI ~ Infectiousness x Sampling scheme',xlab='',
#	ylab='Local Branching Index',
#	cex.main=2,cex.axis=2,cex.lab=2)
#legend('topright',col=c('blue','darkblue'),legend=c('Cross-sectional sampling','Longitudinal sampling'),
#	pch=15,cex=2,bty='n')
#axis(1, line=0.25,
#	at = c(1,2,3,4),
#	labels = c(
#		'High',
#		'Low',
#		'High',
#		'Low'),
#	tick=F, cex.axis=2)
#axis(1, line=3.00,
#	at = c(2.5),
#	labels = 'Infectiousness',
#	tick=F, cex.axis=2)
#

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

# create the individual sub-figures:
#fig.1 <- getfig(tree1,'Altered infectiousness',titleadjust=0.40)
#fig.2 <- getfig(tree2,'Altered susceptibility',titleadjust=0.40)
#fig.3 <- getfig(tree3,'Altered infectiousness + preferential mixing',titleadjust=0.30)
#fig.4 <- getfig(tree4,'Altered susceptibility + preferential mixing',titleadjust=0.25)
#
## place them in a combined figure:
#mainfig <- plot_grid(fig.1,fig.2,fig.3,fig.4,byrow=T,nrow=2)




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
	tauvals <- seq(-4,2, length=30)
	tauvals <- 10^tauvals

	for(tau in tauvals) crud[,paste0('lbi',tau)] <- lbi(tree, tau=tau)

	# add in a column for the state of the node/tip:
	crud$state <- NA
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'IH'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'IL'

	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]+0))

	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	ratios <- apply(crud[,paste0('lbi',tauvals)], 2, function(x) mean(x[crud$state=='IH'])/mean(x[crud$state=='IL'])  )

	plot(tauvals, ratios, xlab=bquote(tau),ylab='LBI ratio (IH/IL)',  
		type='l',lwd=3, ylim=c(0, max(ratios)*1.5),
		cex.axis=1.5,cex.lab=1.5, main=title)
}


pdf(file='figures/nullmodelfigs/super.pdf',height=6,width=6)
gettauplot(tree3, title='Superspreading')
dev.off()

pdf(file='figures/nullmodelfigs/homogeneous.pdf',height=6,width=6)
gettauplot(tree1, title='Homogeneous infectiousness')
dev.off()


