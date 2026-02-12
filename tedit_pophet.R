rm(list=ls())

library(ggnewscale)

source("compile_pophet_and_superspreading_model.R")

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

gettheta <- function(x1,x2,x3,x4){
	#of those with active TB, what fraction are superspreaders?
	pH <- 0.1

	# infectiousness multipliers
	f <- x1 # fraction of transmission from superspreaders

	# parameterize in terms of equilibrium values of S, L, I
	Seq <- 2/3 * Npop
	Ieq <- 0.01 * Npop
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
	zeta1 <- 1
	zeta2 <- x2

	# infectiousness
	thetaL <- (1-f)/(1-pH)
	thetaH <- f/pH


	# contact rates within/across groups:
	c <- beta
	c11 <- c
	c12 <- c
	c21 <- c
	c22 <- c

	cij <- matrix(c(c11,c12,c21,c22),byrow=T,nrow=2)
	# make the off-diagonals smaller for preferential mixing:
	diag(cij) <- x3 * diag(cij)
	cij[1,2] <- x4 * cij[1,2]
	cij[2,1] <- x4 * cij[2,1]
	#cij <- cij / max(eigen(cij)$values)
	#cij <- cij*c

	c11 <- cij[1,1]	
	c12 <- cij[1,2]	
	c21 <- cij[2,1]	
	c22 <- cij[2,2]	

	theta <- list(
		beta=beta,
		zeta1 = zeta1,
		zeta2 = zeta2,
		thetaH = thetaH,
		thetaL = thetaL,
		c11 = c11,
		c12 = c12,
		c21 = c21,
		c22 = c22,
		gamma = gamma,
		sigma = sigma,
		rho = rho,
		#pmu = pmu,
		pH = pH,
		pfast=pfast
	)
	return(theta)
}

# 6 parametrizations:
theta1 <- gettheta(0.1,5,1,1)   # no superspreading, altsus only
theta2 <- gettheta(0.1,1,2,0.01) # no superspreading, cij only
theta3 <- gettheta(0.1,5,2,0.01) # no superspreading, altsus and cij 

theta4 <- gettheta(0.9,5,1,1)   # superspreading, altsus only
theta5 <- gettheta(0.9,1,2,0.01) # superspreading, cij only
theta6 <- gettheta(0.9,5,2,0.01) # superspreading, altsus and cij 

# set the timestep size:
dT <- 1

# specify initial states:
initialStates <- c(Npop/2,0,1,0,Npop/2,0,0,0)
names(initialStates) <- c('S1','L1','IH1','IL1','S2','L2','IH2','IL2') 

out1 <-  sir_simu(
   paramValues = as.list(theta1),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
#   seed=919755)

out2 <-  sir_simu(
   paramValues = as.list(theta2),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
   #seed=280361)

out3 <-  sir_simu(
   paramValues = as.list(theta3),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
   #seed=280361)

out4 <-  sir_simu(
   paramValues = as.list(theta4),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
   #seed=280361)

out5 <-  sir_simu(
   paramValues = as.list(theta5),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)
   #seed=280361)

out6 <-  sir_simu(
   paramValues = as.list(theta6),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100,
   seed=280361)


# use a smaller dataframe for plotting (just sample rows)
traj1 <- out1$traj
if(dim(traj1)[1]<500){plottraj1<-traj1}else{
plottraj1 <- traj1[c(1:1500,sample(1:dim(traj1)[1],1000)),]
plottraj1 <- plottraj1[order(plottraj1$Time),]}
# transform time back to date for plotting:
plottraj1$date <- as.Date(plottraj1$Time * 365)


traj2 <- out2$traj
plottraj2 <- traj2[c(1:1500,sample(1:dim(traj2)[1],1000)),]
plottraj2 <- plottraj2[order(plottraj2$Time),]
# transform time back to date for plotting:
plottraj2$date <- as.Date(plottraj2$Time * 365)



# plot the trajectories of out1 and out2:

#pdf(file='nullmodelfigs/nullmodeltrajs.pdf',height=5,width=9)

par(mfrow=c(1,2))

# low prevalence:
plot(log10(L1+L2)~date,plottraj1,type='l',col='#005AB5',main='Low prevalence',lwd=2.5,
	ylab=bquote(Log[10](.('No. of infections'))),xlab='Date',lty='dotted', ylim=c(0,6))

lines(log10(IL1+IL2)~date,plottraj1,type='l',col='#005AB5',main='Active TB',lwd=2.5,
	ylab='No. of infections',xlab='Date',lty='dashed')
lines(log10(IH1+IH2)~date,plottraj1,type='l',col='#005AA0',lwd=2.5)
legend('top',col=c('#005AB5','#005AA0'),
	lty=c(3,2,1),lwd=2.5,legend=c(bquote(L),bquote(I[1]),bquote(I[2])),
	cex=1.3)

# high prevalence:
plot(log10(L1+L2)~date,plottraj2,type='l',col='#005AB5',main='High prevalence',lwd=2.5,
	ylab=bquote(Log[10](.('No. of infections'))),xlab='Date',lty='dotted', ylim=c(0,6))

lines(log10(IL1+IL2)~date,plottraj2,type='l',col='#005AB5',main='Active TB',lwd=2.5,
	ylab='No. of infections',xlab='Date',lty='dashed')
lines(log10(IH1+IH2)~date,plottraj2,type='l',col='#005AA0',lwd=2.5)
#legend('topleft',col=c('#005AB5','#005AA0'),
#	lty=c(3,2,1),lwd=2.5,legend=c(bquote(L),bquote(I[1]),bquote(I[2])),
#	cex=1.3,bty='n')
#

#dev.off()

# use tms.lin4 to generate simulated trees with the same height (approx.) as the
#	empirical lineage 4 tree

pH <- 0.1

getsimtree <- function(tms,output){
	simulate_tree(
	simuResults=output,
	dates=c(tms),
	deme=c('IH1','IL1','L1','IH2','IL2','L2'),
	sampled=c(
		IH1=pH*0.5,
		IL1=(1-pH)*0.5,
		IH2=pH*0.5,
		IL2=(1-pH)*0.5),
	root = 'IH1',
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)

}

tree1 <- getsimtree(tms.lin4,out1)
tree2 <- getsimtree(tms.lin4,out2)
tree3 <- getsimtree(tms.lin4,out3)
tree4 <- getsimtree(tms.lin4,out4)
tree5 <- getsimtree(tms.lin4,out5)
tree6 <- getsimtree(tms.lin4,out6)

pdf(file='figures/pophetmodelfigs/trees.pdf',width=6,height=9)
par(mfrow=c(3,2))

tree <- tree1
plot(tree, type='phylogram',show.tip.label=F,main='Altered susceptibility')
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)


tree <- tree4
plot(tree, type='phylogram',show.tip.label=F,main='Altered susceptibility +\n superspreading')
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)

tree <- tree2
plot(tree, type='phylogram',show.tip.label=F,main='Preferential mixing')
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)

tree <- tree5
plot(tree, type='phylogram',show.tip.label=F,main='Preferential mixing +\n superspreading')
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)
legend('topleft',legend=c('Group 1','Group2'),pch=19, col=c('#005AB5','#DC3220'),bty='n')

tree <- tree3
plot(tree, type='phylogram',show.tip.label=F,main='Altered susceptibility +\n preferential mixing')
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)

tree <- tree6
plot(tree, type='phylogram',show.tip.label=F,main='Altered susceptibility +\n preferential mixing +\n superspreading')
tips_cols <- ifelse(grepl(x=tree$tip.label,pattern = c("1_")),"#005AB5","#DC3220")
nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="+IH1+"),"#005AB5",
	ifelse(grepl(x=tree$node.label,pattern="+IL1+"),'#005AB5','#DC3220') )
tiplabels(pch=20,col=tips_cols)
nodelabels(pch=20,col=nodes_cols)

dev.off()

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


# load the LBI function:
source('lbi.R')


## make plots of each tree
p1 <- ggtree(tree1,layout='circular')
p2 <- ggtree(tree2,layout='circular')

p1 <- p1 + labs(title='') + theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))
p2 <- p2 + labs(title='') + theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))

ptrees <- plot_grid(p1,p2,nrow=1)

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

	# Add in a column for group status
	inds1 <- grepl(dat$label, pattern='IL1')
	inds2 <- grepl(dat$label, pattern='IH1')
	inds <- inds1+inds2
	dat$Group <- inds
	dat$Group <- dat$Group + 1
	dat$Group <- as.character(dat$Group) 

        # use a heatmap to visualize the statistics:
	# First, make a heatmap for group status using a discrete colour scale

	grpdat <- dat[,'Group']
	names(grpdat) <- rownames(dat)
	grpdat <- as.data.frame(grpdat)
	colnames(grpdat) <- 'Group'

	heatfiggrp <- gheatmap(plin4, grpdat, offset=0, width=0.2, colnames_angle=-45,
			colnames_offset_y=-30, hjust=0.5, font.size=8,
			colnames_position='bottom', color=NA) +
			scale_fill_discrete(name='Group', breaks=c('1','2'),
				type=c('gray90','grey24')) + 
			theme(panel.spacing = unit(0,'pt'), legend.title=element_text(size=18),
				legend.text=element_text(size=18))

        heatfig <-  gheatmap(heatfiggrp+new_scale_fill(),
			dat[,c('LBI',
			'THD',
                        'Cluster Size (single)',
                        'Cluster Size (complete)',
                        'No. of close relatives')],
               	offset=40, font.size=8, colnames_angle = -45,
		colnames=T, colnames_position="bottom", hjust=0.0,
                colnames_offset_y=-30) +
                scale_fill_continuous(name='Value of\nstatistic',
                low='#FEFE62',high='#5D3A9B') +
                theme(plot.margin=unit(c(1,1,3,1),'cm')) +
                coord_cartesian(clip = 'off') +
                ggtitle(title) +
                theme(plot.title=element_text(hjust=0.5,size=18,face="bold"))


        return(list(heatfig,dat))
}


#fig <- plot_grid(getmainplot(tree1,title='Low prevalence')[[1]],
#		getmainplot(tree2,title='High prevalence')[[1]],
#		nrow=1)

fig1 <- getmainplot(tree1, title='1',taulbi=10)
fig2 <- getmainplot(tree2, title='2',taulbi=10)
fig3 <- getmainplot(tree3, title='3',taulbi=10)
fig4 <- getmainplot(tree4, title='4',taulbi=10)
fig5 <- getmainplot(tree5, title='5',taulbi=10)
fig6 <- getmainplot(tree6, title='6',taulbi=10)

fig <- plot_grid(fig1[[1]], fig2[[1]], fig3[[1]], fig4[[1]], fig5[[1]], fig6[[1]],
	nrow=3,byrow=F)

dat1 <- fig1[[2]]
dat2 <- fig2[[2]]
dat3 <- fig3[[2]]
dat4 <- fig4[[2]]
dat5 <- fig5[[2]]
dat6 <- fig6[[2]]


ggsave(fig, file='figures/pophetmodelfigs/pophetmodeltrees.png', dpi=600, height=15, width=12)

gettbls <- function(tree){
	tbls <-	tree$edge.length[ tree$edge[,2] <= length(tree$tip.label) ]
	return(tbls)
}

pdf(file='figures/pophetmodelfigs/boxplots.pdf',height=9,width=5)
par(mfrow=c(3,2))
boxplot(LBI~Group,dat1, main='Altered susceptibility')
boxplot(LBI~Group,dat4, main='Altered susceptibility +\n superspreading')
boxplot(LBI~Group,dat2, main='Preferential mixing')
boxplot(LBI~Group,dat5, main='Preferential mixing +\n superspreading')
boxplot(LBI~Group,dat3, main='Altered susceptibility +\n Preferential mixing')
boxplot(LBI~Group,dat6, main='Altered susceptibility +\n Preferential mixing +\n superspreading')

dev.off()








