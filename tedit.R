rm(list=ls())

source("compile_variant_and_superspreading_model.R")

# 2 useful seeds:
# with Ieq <- 0.01, use 379920
# with Ieq <- 0.001, use 509262

#############
#############

source('loaddates.R')

tms.lin1$hiv <- 'I'
tms.lin2$hiv <- 'I'
tms.lin3$hiv <- 'I'
tms.lin4$hiv <- 'I'

tms.lin1 <- tms.lin1[,1]
tms.lin2 <- tms.lin2[,1]
tms.lin3 <- tms.lin3[,1]
tms.lin4 <- tms.lin4[,1]

# let's simulate a lineage 4 tree for the time period

# this is the julian date of the most recent sampling event:
tmax <- max(tms.lin4)

# we think we want to simulate the TMRCA about 400 years prior
tstart <- tmax-175*365

# time to simulate:
time <- c(tstart,tmax)/365


# let's put the sampling times on the same scale:
tms.lin4 <- tms.lin4/365

# what is the total population size?
Npop <- 2.0 * 10^5

#of those with active TB, what fraction are superspreaders?
pH <- 0.1

# infectiousness multipliers
f <- .9 # fraction of transmission from superspreaders

# parameterize in terms of equilibrium values of S, L, I
Seq <- 2/3 * Npop
Ieq <- 0.0050 * Npop
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
pmu <- 0.005


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


dT <- 1

# specify initial states:
# set phiv% of the population into the HIV state:
initialStates <- c(Npop,0,1,0,0,0,0)
names(initialStates) <- c('S','L','IH','IL','M','JH','JL') 

out <-  sir_simu(
   paramValues = as.list(theta),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100,
   seed=509262)

traj <- out$traj

# calculate the proportion of infections caused by the variant:
traj$p <- with(traj, (JH+JL)/(IH+IL+JH+JL) )


plottraj <- traj[c(1:1500,sample(1:dim(traj)[1],1000)),]
plottraj <- plottraj[order(plottraj$Time),]

# plot the trajectories:

par(mfrow=c(1,3))
# resident:
plot(IL~Time,plottraj,type='l',col='#005AB5',main='Active TB',lwd=2.5)
lines(IH~Time,plottraj,type='l',col='#005AA0',lwd=2.5)
lines(JL~Time,plottraj,type='l',col='#DC3220',lwd=2.5)
#variant:
plot(L~Time,plottraj,type='l',col='#005AB5',main='Latent TB',lwd=2.5)
lines(M~Time,plottraj,type='l',col='#DC3220',lwd=2.5)
#add legend:
legend('right',col=c('#005AB5','#DC3220'),lty=1,lwd=2.5,legend=c('resident','variant'),
	cex=1.3)


# first, let's figure out when the variant reached x% of all infections:
t.01 <- traj[max(which(traj$p <= .01)), 'Time']
t.05 <- traj[max(which(traj$p <= .05)), 'Time']
t.10 <- traj[max(which(traj$p <= .10)), 'Time']

# let's use the same spacing of sampling events, but shift them so the
#	first sampling event coincides with t.x

tms.01 <- tms.lin4
tms.01 <- tms.01 - min(tms.01) + t.01

tms.05 <- tms.lin4
tms.05 <- tms.05 - min(tms.05) + t.05

tms.10 <- tms.lin4
tms.10 <- tms.10 - min(tms.10) + t.10


getsimtree.01 <- function(tms){
	simulate_tree(
	simuResults=out,
	dates=c(tms),
	deme=c('IH','IL','L','JH','JL','M'),
	sampled=c(
		IH=pH*0.99,
		IL=(1-pH)*0.99,
		JH=pH*0.01,
		JL=(1-pH)*0.01),
	root = 'IH',
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)

}

getsimtree.05 <- function(tms){
	simulate_tree(
	simuResults=out,
	dates=c(tms),
	deme=c('IH','IL','L','JH','JL','M'),
	sampled=c(
		IH=pH*0.95,
		IL=(1-pH)*0.95,
		JH=pH*0.05,
		JL=(1-pH)*0.05),
	root = 'IH',
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)

}

getsimtree.10 <- function(tms){
	simulate_tree(
	simuResults=out,
	dates=c(tms),
	deme=c('IH','IL','L','JH','JL','M'),
	sampled=c(
		IH=pH*0.90,
		IL=(1-pH)*0.90,
		JH=pH*0.10,
		JL=(1-pH)*0.10),
	root = 'IH',
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)

}

tree.01 <- getsimtree.01(tms.01)
tree.05 <- getsimtree.05(tms.05)
tree.10 <- getsimtree.10(tms.10)


par(mfrow=c(1,3))

## variant @ 1%

tree <- tree.01
plot(tree, type='fan',show.tip.label=F)
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern="I"),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="I"),"#005AB5","#DC3220")
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)

x <- cophenetic(tree)
tau <- mean(x)*.01
dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
dat <- as.data.frame(dat)
colnames(dat) <- 'thd'
dat$label = rownames(dat)
dat <- dat[,c(2,1)]
dat$var <- apply(x,1,var)
dat$min <- apply(x,1,function(x) min(x[x>0]) )
source('lbi.R')

dat$lbi <- lbi(tree,tau=50)[1:length(tree$tip.label)]

# now add in states:
dat$state <- 0
dat[grep('I',tree$tip.label),'state'] <- 'I'
dat[grep('J',tree$tip.label),'state'] <- 'J'

## variant @ 5%

tree <- tree.05
plot(tree, type='fan',show.tip.label=F)
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern="I"),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="I"),"#005AB5","#DC3220")
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)

## variant @ 10%

tree <- tree.10
plot(tree, type='fan',show.tip.label=F)
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern="I"),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="I"),"#005AB5","#DC3220")
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)


# load the LBI function:
source('lbi.R')

tree1 <- tree.01
tree2 <- tree.05
tree3 <- tree.10

# make plots of each tree
p1 <- ggtree(tree1,layout='circular')
p2 <- ggtree(tree2,layout='circular')
p3 <- ggtree(tree3,layout='circular')

p1 <- p1 + labs(title='') + theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))
p2 <- p2 + labs(title='') + theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))
p3 <- p3 + labs(title='') + theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))

ptrees <- plot_grid(p1,p2,p3,nrow=1)

## mimic the data figure:


getmainplot <- function(tree,taulbi=4,tauthd=5,taurels=6,tauclust=6,title='title'){

        # calculate tree height:
        treeheight <- max(node.depth.edgelength(tree))

        # use cophenetic distances to calculate THD and No. of close relatives:
        x = cophenetic(tree)

        # calculate statistics for the tree (THD, LBI, No. of close relatives):
        #tauthd <- 5 #bandwidth for THD
        #taulbi <- 4 #bandwidth for LBI
        #tauclust <- 12 #threshold distance for No. of Close Relatives

        # calculate THD from cophenetic distances:
        dat <- apply(x, 1, function(z) sum(exp(-z/tauthd)))

        # organize into a dataframe:
        dat <- as.data.frame(dat)
        colnames(dat) <- 'THD'
        dat$label = rownames(dat)
        dat <- dat[,c(2,1)]

        # calculate LBI directly from the tree:
        dat$LBI <- lbi(tree,tau=taulbi)[1:length(tree$tip.label)]

	dat$LBI <- dat$LBI - min(dat$LBI)
	dat$THD <- dat$THD - min(dat$THD)

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

	# create a ggtree plot:
        plin4 <- ggtree(tree, layout='rectangular')

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
                scale_fill_continuous(name='Value of\nstatistic',
                low='#FEFE62',high='#5D3A9B') +
                theme(plot.margin=unit(c(1,1,3,1),'cm')) +
                coord_cartesian(clip = 'off') +
                ggtitle(title) +
                theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))


        return(list(heatfig,dat))
}


fig <- plot_grid(getmainplot(tree1,title='Low progression')[[1]],
		getmainplot(tree2,title='Medium progression')[[1]],
		getmainplot(tree3, title='High progression')[[1]],
		nrow=1)

#ggsave(fig, file='modelAfigs/modelA_trees_v2.png', dpi=600, height=8, width=12)

gettbls <- function(tree){
	tbls <-	tree$edge.length[ tree$edge[,2] <= length(tree$tip.label) ]
	return(tbls)
}


dat1 <- getmainplot(tree1)[[2]]
dat2 <- getmainplot(tree2)[[2]]
dat3 <- getmainplot(tree3)[[2]]

# calculate LTT, LBI, and TBL for the empirical tree:
lin4phi <- read.nexus('/Users/abeams/Documents/projects/TB/lineage_4_delphy/alex/lineage4_delphy.trees')
lin4tree <- lin4phi[[200]]

lin4ltt <- ltt.plot.coords(lin4tree)
lin4lbi <- lbi(lin4tree,4)[1:length(lin4tree$tip.label)]
lin4tbl <- gettbls(lin4tree)


# make a panel displaying LTT, distribution of LBI, and tbl distribution:

#pdf(file='modelAfigs/modelA_stats_v2.pdf',height=6,width=16)

par(mfrow=c(1,3),mar=c(5,5,4,2))
# 1. make LTT plots for the trees:
ltt1 <- ltt.plot.coords(tree1)
ltt2 <- ltt.plot.coords(tree2)
ltt3 <- ltt.plot.coords(tree3)

inds = 1:dim(ltt1)[1]
plot(ltt1[inds,'time'], ltt1[inds,'N'], type = "l",col = "blue",lwd = 1,
	xlab = "Time",ylab = "No. of lineages",
	main = "Lineage-through-time plot",cex.main=2,cex.axis=2,cex.lab=2,cex=2)
lines(ltt2[inds,'time'],ltt2[inds,'N'],lty=1,lwd=3,col='darkblue')
lines(ltt3[inds,'time'],ltt3[inds,'N'],lty=1,lwd=5,col='navy')
lines(lin4ltt[,'time'],lin4ltt[,'N'],lty=1,lwd=4,col='red')

# 2. display LBI distributions:

#plot(density(dat2$LBI),lwd=3, xlab='LBI',ylab='Density',main='LBI distribution',
#	cex.main=2,cex.axis=2,cex.lab=2,cex=2,lty=2,ylim=c(0,0.9),xlim=c(0,30))
#lines(density(dat1$LBI),lwd=3, xlab='LBI',ylab='Density',main='LBI distribution',lty=1)
#lines(density(dat3$LBI),lwd=3, xlab='LBI',ylab='Density',main='LBI distribution',lty=3)
#lines(density(lin4lbi),lwd=3, xlab='LBI',ylab='Density',main='LBI distribution',lty=4,col='red')
#legend('right',legend=c('Low progression','Medium progression','High progression','Empirical tree'),
#	lty=c(1,2,3,4),lwd=3,cex=2,col=c('black','black','black','red'))

plot(ecdf(dat2$LBI),xlab='LBI',ylab='Cumulative density',main='LBI distribution',
	cex.main=2,cex.axis=2,cex.lab=2,cex=2,lty=2,ylim=c(0,1.0),xlim=c(0,30),
	verticals=T, do.points=F, lwd=3, col='darkblue')
lines(ecdf(dat1$LBI), 
	xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
	verticals=T,do.points=F,col='blue',lwd=1)
lines(ecdf(dat3$LBI), 
	xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
	verticals=T, do.points=F,lwd=5, col='navy')
lines(ecdf(lin4lbi), 
	xlab='LBI',ylab='Density',main='LBI distribution',lty=1,col='red',
	verticals=T, do.points=F,lwd=4)
legend('right',legend=c('Low progression','Medium progression','High progression','Empirical tree'),
	lwd=c(1,3,5,4),lty=1,cex=2,col=c('blue','darkblue','navy','red'),bty='n')


# 3. Display tbl distributions

plot(ecdf(gettbls(tree1)),lwd=1, xlab='Terminal branch length',
	ylab='Cumulative density',main='Terminal Branch Length distribution',
	cex.main=2,cex.axis=2,cex.lab=2,cex=2,lty=1,ylim=c(0,1.0),
	verticals=T, do.points=F,col='blue')
lines(ecdf(gettbls(tree2)),lwd=3, xlab='LBI',ylab='Density',	
	main='LBI distribution',lty=1,
	verticals=T, do.points=F,col='darkblue')
lines(ecdf(gettbls(tree3)),lwd=5, xlab='LBI',ylab='Density',
	main='LBI distribution',lty=1,
	verticals=T, do.points=F,col='navy')
lines(ecdf(gettbls(lin4tree)),lwd=4, xlab='LBI',ylab='Density',main='LBI distribution',lty=1,
	col='red',
	verticals=T,do.points=F)







