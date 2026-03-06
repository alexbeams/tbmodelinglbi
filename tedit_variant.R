rm(list=ls())

source("compile_variant_and_superspreading_model.R")


#############
#############
ntests <- 500

# time to simulate:
time <- c(-125,50)

# when should we start to sample?
tm1 <- 45
tms.lin4 <- runif(ntests,min=tm1,max=time[2])
tms.lin4 <- tms.lin4[order(tms.lin4)]


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
t.20 <- traj[max(which(traj$p <= .20)), 'Time']
t.30 <- traj[max(which(traj$p <= .30)), 'Time']
t.40 <- traj[max(which(traj$p <= .40)), 'Time']
t.50 <- traj[max(which(traj$p <= .50)), 'Time']
t.60 <- traj[max(which(traj$p <= .60)), 'Time']
t.70 <- traj[max(which(traj$p <= .70)), 'Time']
t.80 <- traj[max(which(traj$p <= .80)), 'Time']
t.90 <- traj[max(which(traj$p <= .90)), 'Time']


# let's use the same spacing of sampling events, but shift them so the
#	first sampling event coincides with t.x

tms.01 <- tms.lin4
tms.01 <- tms.01 - min(tms.01) + t.01

tms.05 <- tms.lin4
tms.05 <- tms.05 - min(tms.05) + t.05

tms.10 <- tms.lin4
tms.10 <- tms.10 - min(tms.10) + t.10

tms.20 <- tms.lin4
tms.20 <- tms.20 - min(tms.20) + t.20

tms.30 <- tms.lin4
tms.30 <- tms.30 - min(tms.30) + t.30

tms.40 <- tms.lin4
tms.40 <- tms.40 - min(tms.40) + t.40

tms.50 <- tms.lin4
tms.50 <- tms.50 - min(tms.50) + t.50

tms.60 <- tms.lin4
tms.60 <- tms.60 - min(tms.60) + t.60

tms.70 <- tms.lin4
tms.70 <- tms.70 - min(tms.70) + t.70

tms.80 <- tms.lin4
tms.80 <- tms.80 - min(tms.80) + t.80

tms.90 <- tms.lin4
tms.90 <- tms.90 - min(tms.90) + t.90

# set the minimum time for longitudinal sampling to commence:
mintime <- -120

crud <- plottraj[plottraj$Time > mintime,]
crud <- crud[sort(sample(1:dim(crud)[1],ntests)),]
crud$pih <- (1-crud$p) * pH
crud$pil <- (1-crud$p) * (1-pH)
crud$pjh <- crud$p * pH
crud$pjl <- crud$p * (1-pH)

crud$sampled <- apply(crud, 1, function(x) c('IH','IL','JH','JL')[which(t(rmultinom(1,1,prob=x[c('pih','pil','pjh','pjl')])) > 0)])

tms.long <- crud[,c('Time','sampled')]
colnames(tms.long) <- c('Date','Comp')

getsimtree.p <- function(tms,p){
	simulate_tree(
	simuResults=out,
	dates=c(tms),
	deme=c('IH','IL','L','JH','JL','M'),
	sampled=c(
		IH=pH*(1-p),
		IL=(1-pH)*(1-p),
		JH=pH*p,
		JL=(1-pH)*p),
	root = 'IH',
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)
}

getsimtree.long <- function(tms){
	simulate_tree(
	simuResults=out,
	dates=tms,
	deme=c('IH','IL','L','JH','JL','M'),
	root = 'IH',
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)
}

# Simulate trees with sampling once variant reaches proportion p:
tree.01 <- getsimtree.p(tms.01,0.01)
tree.05 <- getsimtree.p(tms.05,0.05)
tree.10 <- getsimtree.p(tms.10,0.10)
tree.20 <- getsimtree.p(tms.20,0.20)
tree.30 <- getsimtree.p(tms.30,0.30)
tree.40 <- getsimtree.p(tms.40,0.40)
tree.50 <- getsimtree.p(tms.50,0.50)
tree.60 <- getsimtree.p(tms.60,0.60)
tree.70 <- getsimtree.p(tms.70,0.70)
tree.80 <- getsimtree.p(tms.80,0.80)
tree.90 <- getsimtree.p(tms.90,0.90)

# Sample longitudinally over the whole sim:
tree.long <- getsimtree.long(tms.long)

# Most of this below is now deprecated, but not all of it

#par(mfrow=c(1,3))
#
### variant @ 1%
#
#tree <- tree.01
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
#
### variant @ 5%
#
#tree <- tree.05
#plot(tree, type='fan',show.tip.label=F)
#tips_cols <- ifelse(grepl(x=tree$tip.label,pattern="I"),"#005AB5","#DC3220")
#nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="I"),"#005AB5","#DC3220")
#tiplabels(pch=20,col=tips_cols)
#nodelabels(pch=20,col=nodes_cols)
#
### variant @ 10%
#
#tree <- tree.10
#plot(tree, type='fan',show.tip.label=F)
#tips_cols <- ifelse(grepl(x=tree$tip.label,pattern="I"),"#005AB5","#DC3220")
#nodes_cols <- ifelse(grepl(x=tree$node.label,pattern="I"),"#005AB5","#DC3220")
#tiplabels(pch=20,col=tips_cols)
#nodelabels(pch=20,col=nodes_cols)
#

# load the LBI function:
source('lbi.R')



#p.01 <- ggtree(tree.01,layout='circular')
#p.05 <- ggtree(tree.05,layout='circular')
#p.10 <- ggtree(tree.10,layout='circular')
#p.long <- ggtree(tree.long,layout='circular')
#
#
## calculate data for tree.01:
#tree <- tree.01
#x <- cophenetic(tree)
#tau <- mean(x)*.1
#dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
#dat <- as.data.frame(dat)
#colnames(dat) <- 'thd'
#dat$label = rownames(dat)
#dat <- dat[,c(2,1)]
#dat$var <- apply(x,1,var)
#dat$min <- apply(x,1,function(x) min(x[x>0]) )
#dat$lbi.1000 <- lbi(tree,tau=1000)[1:length(tree$tip.label)]
#dat$lbi.500 <- lbi(tree,tau=500)[1:length(tree$tip.label)]
#dat$lbi.100 <- lbi(tree,tau=100)[1:length(tree$tip.label)]
#dat$lbi.50 <- lbi(tree,tau=50)[1:length(tree$tip.label)]
#dat$lbi.25 <- lbi(tree,tau=25)[1:length(tree$tip.label)]
#dat$lbi.10 <- lbi(tree,tau=10)[1:length(tree$tip.label)]
#dat$lbi.5 <- lbi(tree,tau=5)[1:length(tree$tip.label)]
#dat$lbi.4 <- lbi(tree,tau=4)[1:length(tree$tip.label)]
#dat$lbi.3 <- lbi(tree,tau=3)[1:length(tree$tip.label)]
#dat$lbi.2.5 <- lbi(tree,tau=2.5)[1:length(tree$tip.label)]
#dat$lbi.2 <- lbi(tree,tau=2)[1:length(tree$tip.label)]
#dat$lbi.1 <- lbi(tree,tau=1)[1:length(tree$tip.label)]
## now add in states:
#dat$state <- 0
#dat[grep('I',tree$tip.label),'state'] <- 'I'
#dat[grep('J',tree$tip.label),'state'] <- 'J'
#dat.01 <- dat
#
## merge data onto the ggtree plot:
#p.01 <- p.01 %<+% dat.01
#
#
## calculate data for tree.05:
#tree <- tree.05
#x <- cophenetic(tree)
#tau <- mean(x)*.1
#dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
#dat <- as.data.frame(dat)
#colnames(dat) <- 'thd'
#dat$label = rownames(dat)
#dat <- dat[,c(2,1)]
#dat$var <- apply(x,1,var)
#dat$min <- apply(x,1,function(x) min(x[x>0]) )
#dat$lbi.1000 <- lbi(tree,tau=1000)[1:length(tree$tip.label)]
#dat$lbi.500 <- lbi(tree,tau=500)[1:length(tree$tip.label)]
#dat$lbi.100 <- lbi(tree,tau=100)[1:length(tree$tip.label)]
#dat$lbi.50 <- lbi(tree,tau=50)[1:length(tree$tip.label)]
#dat$lbi.25 <- lbi(tree,tau=25)[1:length(tree$tip.label)]
#dat$lbi.10 <- lbi(tree,tau=10)[1:length(tree$tip.label)]
#dat$lbi.5 <- lbi(tree,tau=5)[1:length(tree$tip.label)]
#dat$lbi.4 <- lbi(tree,tau=4)[1:length(tree$tip.label)]
#dat$lbi.3 <- lbi(tree,tau=3)[1:length(tree$tip.label)]
#dat$lbi.2.5 <- lbi(tree,tau=2.5)[1:length(tree$tip.label)]
#dat$lbi.2 <- lbi(tree,tau=2)[1:length(tree$tip.label)]
#dat$lbi.1 <- lbi(tree,tau=1)[1:length(tree$tip.label)]
#
## now add in states:
#dat$state <- 0
#dat[grep('I',tree$tip.label),'state'] <- 'I'
#dat[grep('J',tree$tip.label),'state'] <- 'J'
#dat.05 <- dat
#p.05 <- p.05 %<+% dat.05
#
## calculate data for tree.10:
#tree <- tree.10
#x <- cophenetic(tree)
#tau <- mean(x)*.1
#dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
#dat <- as.data.frame(dat)
#colnames(dat) <- 'thd'
#dat$label = rownames(dat)
#dat <- dat[,c(2,1)]
#dat$var <- apply(x,1,var)
#dat$min <- apply(x,1,function(x) min(x[x>0]) )
#dat$lbi.1000 <- lbi(tree,tau=1000)[1:length(tree$tip.label)]
#dat$lbi.500 <- lbi(tree,tau=500)[1:length(tree$tip.label)]
#dat$lbi.100 <- lbi(tree,tau=100)[1:length(tree$tip.label)]
#dat$lbi.50 <- lbi(tree,tau=50)[1:length(tree$tip.label)]
#dat$lbi.25 <- lbi(tree,tau=25)[1:length(tree$tip.label)]
#dat$lbi.10 <- lbi(tree,tau=10)[1:length(tree$tip.label)]
#dat$lbi.5 <- lbi(tree,tau=5)[1:length(tree$tip.label)]
#dat$lbi.4 <- lbi(tree,tau=4)[1:length(tree$tip.label)]
#dat$lbi.3 <- lbi(tree,tau=3)[1:length(tree$tip.label)]
#dat$lbi.2.5 <- lbi(tree,tau=2.5)[1:length(tree$tip.label)]
#dat$lbi.2 <- lbi(tree,tau=2)[1:length(tree$tip.label)]
#dat$lbi.1 <- lbi(tree,tau=1)[1:length(tree$tip.label)]
#
## now add in states:
#dat$state <- 0
#dat[grep('I',tree$tip.label),'state'] <- 'I'
#dat[grep('J',tree$tip.label),'state'] <- 'J'
#dat.10 <- dat
#p.10 <- p.10 %<+% dat.10
#
#
## calculate data for tree.long:
#tree <- tree.long
#x <- cophenetic(tree)
#tau <- mean(x)*.1
#dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
#dat <- as.data.frame(dat)
#colnames(dat) <- 'thd'
#dat$label = rownames(dat)
#dat <- dat[,c(2,1)]
#dat$var <- apply(x,1,var)
#dat$min <- apply(x,1,function(x) min(x[x>0]) )
#dat$lbi.1000 <- lbi(tree,tau=1000)[1:length(tree$tip.label)]
#dat$lbi.500 <- lbi(tree,tau=500)[1:length(tree$tip.label)]
#dat$lbi.100 <- lbi(tree,tau=100)[1:length(tree$tip.label)]
#dat$lbi.50 <- lbi(tree,tau=50)[1:length(tree$tip.label)]
#dat$lbi.25 <- lbi(tree,tau=25)[1:length(tree$tip.label)]
#dat$lbi.10 <- lbi(tree,tau=10)[1:length(tree$tip.label)]
#dat$lbi.5 <- lbi(tree,tau=5)[1:length(tree$tip.label)]
#dat$lbi.4 <- lbi(tree,tau=4)[1:length(tree$tip.label)]
#dat$lbi.3 <- lbi(tree,tau=3)[1:length(tree$tip.label)]
#dat$lbi.2.5 <- lbi(tree,tau=2.5)[1:length(tree$tip.label)]
#dat$lbi.2 <- lbi(tree,tau=2)[1:length(tree$tip.label)]
#dat$lbi.1 <- lbi(tree,tau=1)[1:length(tree$tip.label)]
#
#
### Want to generate plots of the trees with edges annotated with LBI,
### the variant clade highlighted, and display LBI~time:
#
## do this for the longitudinal tree:
#tree <- tree.long
#
## want to add in rows for nodes with times and LBIs
#crud <- data.frame(time = tree$tip.height,
#	label = tree$tip.label)
#
## the node labels have the times; extract these:
#m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
#crud2 <- data.frame(time=as.numeric(m), 
#		label = names(m))
## the node labels are super clunky, but we need to keep them to match with the tree
#crud <- rbind(crud,crud2)
#
## rearrange columns with labels first:
#crud <- crud[,c(2,1)]
#
## calculate LBI for the tips and the nodes:
#crud$lbi20 <- lbi(tree, tau=20)
#
## add in a column for the state of the node/tip:
#crud$state <- NA
#crud[grep('I',crud$label[1:ntests]),'state'] <- 'I'
#crud[grep('J',crud$label[1:ntests]),'state'] <- 'J'
#
#nodenms <- sapply(crud[(ntests+1):(ntests+ntests-1),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]-1))
#
#crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms
#
#
## make a ggtree object:
#p <- ggtree(tree)
#
## merge the dataframe with LBI onto it:
#p <- p %<+% crud
#
## plot it!
#p20 <- p + geom_tippoint(aes(col=lbi20)) + geom_nodepoint(aes(col=lbi20))
#
## this one might be best for our purposes:
#
#p20 <- p + geom_tippoint(aes(col=lbi20)) + geom_nodepoint(aes(col=lbi20))
#
#vrns <- crud[1:Ntip(tree),]
#vrns <- vrns[vrns$state=='J','label']
#vrntmrca <- getMRCA(tree, vrns )
#
## double check that this works correctly:
##p + geom_hilight(node=vrntmrca, fill='red') +
##	geom_tippoint(aes(col=state))
## ^this is working correctly, identifying the variant clade. 
#
##ptree <- p + aes(col=lbi20) + 
##	geom_hilight(node=vrntmrca, fill='red',alpha=0.10)
##ptree.long <- ptree
#
#ptree <- p + #aes(col=lbi10) + 
#	geom_tippoint(aes(col=lbi20)) + geom_nodepoint(aes(col=lbi20)) +
#	scale_colour_continuous(name='LBI',
#        low='#FEFE62',high='#5D3A9B') + 
#	geom_hilight(node=vrntmrca, fill='gray',alpha=0.5)
#ptree.long <- ptree + labs(title='Longitudinal observations') + 
#	theme(plot.title=element_text(hjust=0.5,size=16,face='bold'))
#
### make a panel plot
#
## reorder the data so that J's appear at the end (and get plotted last):
#crudI <- crud[which(crud$state=='I'),]
#crudJ <- crud[which(crud$state=='J'),]
#
#crud = rbind(crudI,crudJ)
#
#
#plbi.long <- ggplot(crud) + geom_point(aes(x=time-min(time),y=lbi20,col=state)) +
#	scale_color_discrete(name='Variant',breaks=c('I','J'),
#			type=c('#005AB5','#DC3220'))  + theme_classic() + 
#	ylab('LBI') + xlab('Time (years)')
#
#alnd <- align_plots(ptree.long,plbi.long,align='v',axis='lr')
#fig.long <- plot_grid(alnd[[1]],alnd[[2]],ncol=1)
# 
#
## do this for the 10% tree:
#tree <- tree.10
#
## want to add in rows for nodes with times and LBIs
#crud <- data.frame(time = tree$tip.height,
#	label = tree$tip.label)
#
## the node labels have the times; extract these:
#m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
#crud2 <- data.frame(time=as.numeric(m), 
#		label = names(m))
## the node labels are super clunky, but we need to keep them to match with the tree
#crud <- rbind(crud,crud2)
#
## rearrange columns with labels first:
#crud <- crud[,c(2,1)]
#
## calculate LBI for the tips and the nodes:
#crud$lbi10 <- lbi(tree, tau=10)
#
## add in a column for the state of the node/tip:
#crud$state <- NA
#crud[grep('I',crud$label[1:ntests]),'state'] <- 'I'
#crud[grep('J',crud$label[1:ntests]),'state'] <- 'J'
#
#nodenms <- sapply(crud[(ntests+1):(2*ntests-1),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]-1))
#
#crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms
#
## make a ggtree object:
#p <- ggtree(tree)
#
## merge the dataframe with LBI onto it:
#p <- p %<+% crud
#
#vrns <- crud[1:Ntip(tree),]
#vrns <- vrns[vrns$state=='J','label']
#vrntmrca <- getMRCA(tree, vrns )
#
## double check that this works correctly:
##p + geom_hilight(node=vrntmrca, fill='red') +
##	geom_tippoint(aes(col=state))
## ^this is working correctly, identifying the variant clade. 
#
#ptree <- p + #aes(col=lbi10) + 
#	geom_tippoint(aes(col=lbi10)) + geom_nodepoint(aes(col=lbi10)) +
#	scale_colour_continuous(name='LBI',
#        low='#FEFE62',high='#5D3A9B') + 
#	geom_hilight(node=vrntmrca, fill='gray',alpha=0.5)
#
#ptree.10 <- ptree + labs(title='Cross-sectional observations\nvariant @ 10%') + 
#	theme(plot.title=element_text(hjust=0.5,size=16,face='bold'))
#
## reorder the data so that J's appear at the end (and get plotted last):
#crudI <- crud[which(crud$state=='I'),]
#crudJ <- crud[which(crud$state=='J'),]
#
#crud = rbind(crudI,crudJ)
#
### make a panel plot
#tree_xlim <- layer_scales(ptree.10)$x$get_limits()
#
#plbi.10 <- ggplot(crud) + geom_point(aes(x=time-min(time),y=lbi10,col=state)) +
#	scale_color_discrete(name='Variant',breaks=c('I','J'),
#			type=c('#005AB5','#DC3220'))  + theme_classic() + 
#	ylab('LBI') + xlab('Time (years)')
#
#alnd <- align_plots(ptree.10,plbi.10,align='v',axis='lr')
#
#fig.10 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1) 
#
## do this for the 90% tree:
#tree <- tree.90
#
## want to add in rows for nodes with times and LBIs
#crud <- data.frame(time = tree$tip.height,
#	label = tree$tip.label)
#
## the node labels have the times; extract these:
#m<- sapply(tree$node.label, function(z) substr(z, regexpr('=',z)[1]+1, regexpr(',re',z)[1]-1 ) ) 
#crud2 <- data.frame(time=as.numeric(m), 
#		label = names(m))
## the node labels are super clunky, but we need to keep them to match with the tree
#crud <- rbind(crud,crud2)
#
## rearrange columns with labels first:
#crud <- crud[,c(2,1)]
#
## calculate LBI for the tips and the nodes:
#crud$lbi10 <- lbi(tree, tau=10)
#
## add in a column for the state of the node/tip:
#crud$state <- NA
#crud[grep('I',crud$label[1:ntests]),'state'] <- 'I'
#crud[grep('J',crud$label[1:ntests]),'state'] <- 'J'
#
#nodenms <- sapply(crud[(ntests+1):(2*ntests-1),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]-1))
#
#crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms
#
## make a ggtree object:
#p <- ggtree(tree)
#
## merge the dataframe with LBI onto it:
#p <- p %<+% crud
#
#vrns <- crud[1:Ntip(tree),]
#vrns <- vrns[vrns$state=='J','label']
#vrntmrca <- getMRCA(tree, vrns )
#
## double check that this works correctly:
##p + geom_hilight(node=vrntmrca, fill='red') +
##	geom_tippoint(aes(col=state))
## ^this is working correctly, identifying the variant clade. 
#
#ptree <- p + #aes(col=lbi10) + 
#	geom_tippoint(aes(col=lbi10)) + geom_nodepoint(aes(col=lbi10)) +
#	scale_colour_continuous(name='LBI',
#        low='#FEFE62',high='#5D3A9B') + 
#	geom_hilight(node=vrntmrca, fill='gray',alpha=0.5)
#
#ptree.90 <- ptree + labs(title='Cross-sectional observations\nvariant @ 90%') + 
#	theme(plot.title=element_text(hjust=0.5,size=16,face='bold'))
#
## reorder the data so that J's appear at the end (and get plotted last):
#crudI <- crud[which(crud$state=='I'),]
#crudJ <- crud[which(crud$state=='J'),]
#
#crud = rbind(crudI,crudJ)
#
### make a panel plot
#tree_xlim <- layer_scales(ptree.10)$x$get_limits()
#
#plbi.90 <- ggplot(crud) + geom_point(aes(x=time-min(time),y=lbi10,col=state)) +
#	scale_color_discrete(name='Variant',breaks=c('I','J'),
#			type=c('#005AB5','#DC3220'))  + theme_classic() + 
#	ylab('LBI') + xlab('Time (years)')
#
#alnd <- align_plots(ptree.90,plbi.90,align='v',axis='lr')
#
#fig.90 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1) 
#
#
#fig.cross <- plot_grid(fig.10,fig.90,ncol=2)
#fig.var <- plot_grid(fig.10,fig.90,fig.long,ncol=3)
#
#ggsave(fig.var, file='figures/varmodelfigs/phylowave.pdf',height=7,width=21)
#
#
## now add in states:
#dat$state <- 0
#dat[grep('I',tree$tip.label),'state'] <- 'I'
#dat[grep('J',tree$tip.label),'state'] <- 'J'
#dat.long <- dat
#p.long <- p.long %<+% dat.long
#
#rm(dat)
#rm(x)
#rm(tau)
#
#
#p.10r <- ggtree(tree.10, layout='rectangular') +
#	geom_tree(size=0.1,color='gray60')
#heatdat.10r <- dat.10
#
## standardize LBI values for plotting:
#z.lbi.1 <- heatdat.10r$lbi.1
#z.lbi.1 <- z.lbi.1-mean(z.lbi.1)
#z.lbi.1 <- z.lbi.1/sd(z.lbi.1)
#
#z.lbi.2<- heatdat.10r$lbi.2
#z.lbi.2 <- z.lbi.2-mean(z.lbi.2)
#z.lbi.2 <- z.lbi.2/sd(z.lbi.2)
#
#z.lbi.2.5<- heatdat.10r$lbi.2.5
#z.lbi.2.5 <- z.lbi.2.5-mean(z.lbi.2.5)
#z.lbi.2.5 <- z.lbi.2.5/sd(z.lbi.2.5)
#
#z.lbi.3<- heatdat.10r$lbi.3
#z.lbi.3 <- z.lbi.3-mean(z.lbi.3)
#z.lbi.3 <- z.lbi.3/sd(z.lbi.3)
#
#z.lbi.4<- heatdat.10r$lbi.4
#z.lbi.4 <- z.lbi.4-mean(z.lbi.4)
#z.lbi.4 <- z.lbi.4/sd(z.lbi.4)
#
#z.lbi.5<- heatdat.10r$lbi.5
#z.lbi.5 <- z.lbi.5-mean(z.lbi.5)
#z.lbi.5 <- z.lbi.5/sd(z.lbi.5)
#
#z.lbi.10 <- heatdat.10r$lbi.10
#z.lbi.10 <- z.lbi.10-mean(z.lbi.10)
#z.lbi.10 <- z.lbi.10/sd(z.lbi.10)
#
#z.lbi.25 <- heatdat.10r$lbi.25
#z.lbi.25 <- z.lbi.25-mean(z.lbi.25)
#z.lbi.25 <- z.lbi.25/sd(z.lbi.25)
#
#z.lbi.50 <- heatdat.10r$lbi.50
#z.lbi.50 <- z.lbi.50-mean(z.lbi.50)
#z.lbi.50 <- z.lbi.50/sd(z.lbi.50)
#
#z.lbi.100 <- heatdat.10r$lbi.100
#z.lbi.100 <- z.lbi.100-mean(z.lbi.100)
#z.lbi.100 <- z.lbi.100/sd(z.lbi.100)
#
#z.lbi.1000 <- heatdat.10r$lbi.1000
#z.lbi.1000 <- z.lbi.1000-mean(z.lbi.1000)
#z.lbi.1000 <- z.lbi.1000/sd(z.lbi.1000)
#
#heatdat.10r$z.lbi.1 <- z.lbi.1
#heatdat.10r$z.lbi.2 <- z.lbi.2
#heatdat.10r$z.lbi.2.5 <- z.lbi.2.5
#heatdat.10r$z.lbi.3 <- z.lbi.3
#heatdat.10r$z.lbi.4 <- z.lbi.4
#heatdat.10r$z.lbi.5 <- z.lbi.5
#heatdat.10r$z.lbi.10 <- z.lbi.10
#heatdat.10r$z.lbi.25 <- z.lbi.25
#heatdat.10r$z.lbi.50 <- z.lbi.50
#heatdat.10r$z.lbi.100 <- z.lbi.100
#heatdat.10r$z.lbi.1000 <- z.lbi.1000
#
#p.10r <- p.10r %<+% heatdat.10r
#
#
##new_labels <- c(expression(tau~"=1"),'10','50','100','1000' )
#
##newlabs <- c('LBI_1','LBI_5','LBI_10','LBI_25','LBI_50','LBI_100','LBI_1000')
#newlabs <- c('LBI_1','LBI_2.5','LBI_5','LBI_10','LBI_25')
#
##colnames(heatdat.10r)[colnames(heatdat.10r) %in% c(paste0('z.lbi.',c(1,5,10,25,50,100,1000)))] <- newlabs
#colnames(heatdat.10r)[colnames(heatdat.10r) %in% c(paste0('z.lbi.',c(1,2.5,5,10,25)))] <- newlabs
#
#vardat <- heatdat.10r[,'state']
#names(vardat) <- rownames(heatdat.10r)
#vardat <- as.data.frame(vardat)
#colnames(vardat) <- 'Variant'
#
## first, we make the heatmap for variant status, which is discrete:
#heatfigvar <- gheatmap(p.10r, vardat, offset=0, width=0.2, colnames_angle=-45,
#		colnames_offset_y = -30, hjust=0.5, font.size=8,
#		colnames_position='bottom', color=NA) + 
#		scale_fill_discrete(name='Variant',breaks=c('I','J'),
#			type=c('#005AB5','#DC3220'))  + 
#		theme(panel.spacing = unit(0,'pt'),legend.title =element_text(size=18),
#			legend.text=element_text(size=18))
#
##heatfigvar <- ggplot(heatdat.10r, aes(x=c(1),y=label)) +
##	geom_tile(aes(fill=state)) + xlab(NULL) + ylab(NULL)
# 
#library(ggnewscale)
#p.10rheatvar <- heatfigvar + new_scale_fill()
#
#
##now, we make another one for the LBI values: 
#heatfig <-  gheatmap(p.10rheatvar,heatdat.10r[,newlabs],
#       	offset=40, font.size=8, colnames_angle = -45,
#	colnames=T, colnames_position="bottom", hjust=0.5, 
#	colnames_offset_y=-30, color=NA)+ 
#	scale_fill_continuous(name='Value of\nstandardized\nLBI',
#        low='#FEFE62',high='#5D3A9B') + 
#	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
#	coord_cartesian(clip = 'off') +
#	ggtitle('Variant at 10%') +
#	theme(panel.spacing=unit(0,'pt'),
#		legend.title=element_text(size=18),
#		legend.text=element_text(size=18), 
#		plot.title=element_text(hjust=0.6,size=24,face='bold') )
#
#ggsave(heatfig, file='figures/varmodelfigs/variant_and_lbi.pdf',height=15,width=15)
#ggsave(heatfig, file='figures/varmodelfigs/variant_and_lbi.png',dpi=600,height=15,width=15)
#                                                                               
#
### Let's make the same LBI heatmap style of plot with a tree sampled longitudinally:
#
#p.longr <- ggtree(tree.long, layout='rectangular') +
#	geom_tree(size=0.1,color='gray60')
#heatdat.longr <- dat.long
#
## standardize LBI values for plotting:
#z.lbi.1 <- heatdat.longr$lbi.1
#z.lbi.1 <- z.lbi.1-mean(z.lbi.1)
#z.lbi.1 <- z.lbi.1/sd(z.lbi.1)
#
#z.lbi.2<- heatdat.longr$lbi.2
#z.lbi.2 <- z.lbi.2-mean(z.lbi.2)
#z.lbi.2 <- z.lbi.2/sd(z.lbi.2)
#
#z.lbi.2.5<- heatdat.longr$lbi.2.5
#z.lbi.2.5 <- z.lbi.2.5-mean(z.lbi.2.5)
#z.lbi.2.5 <- z.lbi.2.5/sd(z.lbi.2.5)
#
#z.lbi.3<- heatdat.longr$lbi.3
#z.lbi.3 <- z.lbi.3-mean(z.lbi.3)
#z.lbi.3 <- z.lbi.3/sd(z.lbi.3)
#
#z.lbi.4<- heatdat.longr$lbi.4
#z.lbi.4 <- z.lbi.4-mean(z.lbi.4)
#z.lbi.4 <- z.lbi.4/sd(z.lbi.4)
#
#z.lbi.5<- heatdat.longr$lbi.5
#z.lbi.5 <- z.lbi.5-mean(z.lbi.5)
#z.lbi.5 <- z.lbi.5/sd(z.lbi.5)
#
#z.lbi.10 <- heatdat.longr$lbi.10
#z.lbi.10 <- z.lbi.10-mean(z.lbi.10)
#z.lbi.10 <- z.lbi.10/sd(z.lbi.10)
#
#z.lbi.25 <- heatdat.longr$lbi.25
#z.lbi.25 <- z.lbi.25-mean(z.lbi.25)
#z.lbi.25 <- z.lbi.25/sd(z.lbi.25)
#
#z.lbi.50 <- heatdat.longr$lbi.50
#z.lbi.50 <- z.lbi.50-mean(z.lbi.50)
#z.lbi.50 <- z.lbi.50/sd(z.lbi.50)
#
#z.lbi.100 <- heatdat.longr$lbi.100
#z.lbi.100 <- z.lbi.100-mean(z.lbi.100)
#z.lbi.100 <- z.lbi.100/sd(z.lbi.100)
#
#z.lbi.500 <- heatdat.longr$lbi.500
#z.lbi.500 <- z.lbi.500-mean(z.lbi.500)
#z.lbi.500 <- z.lbi.500/sd(z.lbi.500)
#
#z.lbi.1000 <- heatdat.longr$lbi.1000
#z.lbi.1000 <- z.lbi.1000-mean(z.lbi.1000)
#z.lbi.1000 <- z.lbi.1000/sd(z.lbi.1000)
#
#heatdat.longr$z.lbi.1 <- z.lbi.1
#heatdat.longr$z.lbi.2 <- z.lbi.2
#heatdat.longr$z.lbi.2.5 <- z.lbi.2.5
#heatdat.longr$z.lbi.3 <- z.lbi.3
#heatdat.longr$z.lbi.4 <- z.lbi.4
#heatdat.longr$z.lbi.5 <- z.lbi.5
#heatdat.longr$z.lbi.10 <- z.lbi.10
#heatdat.longr$z.lbi.25 <- z.lbi.25
#heatdat.longr$z.lbi.50 <- z.lbi.50
#heatdat.longr$z.lbi.100 <- z.lbi.100
#heatdat.longr$z.lbi.500 <- z.lbi.500
#heatdat.longr$z.lbi.1000 <- z.lbi.1000
#
#p.longr <- p.longr %<+% heatdat.longr
#
##newlabs <- c('LBI_1','LBI_2.5','LBI_5','LBI_10','LBI_25')
#newlabs <- c('LBI_5','LBI_10','LBI_50','LBI_100','LBI_500','LBI_1000')
#
##colnames(heatdat.longr)[colnames(heatdat.longr) %in% c(paste0('z.lbi.',c(1,2.5,5,10,25)))] <- newlabs
#colnames(heatdat.longr)[colnames(heatdat.longr) %in% c(paste0('z.lbi.',c(5,10,50,100,500,1000)))] <- newlabs
#
#vardat <- heatdat.longr[,'state']
#names(vardat) <- rownames(heatdat.longr)
#vardat <- as.data.frame(vardat)
#colnames(vardat) <- 'Variant'
#
## first, we make the heatmap for variant status, which is discrete:
#heatfigvarlong <- gheatmap(p.longr, vardat, offset=0, width=0.2, colnames_angle=-45,
#		colnames_offset_y = -30, hjust=0.5, font.size=8,
#		colnames_position='bottom', color=NA) + 
#		scale_fill_discrete(name='Variant',breaks=c('I','J'),
#			type=c('#005AB5','#DC3220'))  + 
#		theme(panel.spacing = unit(0,'pt'),legend.title =element_text(size=18),
#			legend.text=element_text(size=18))
#
# 
#library(ggnewscale)
#p.longrheatvar <- heatfigvarlong + new_scale_fill()
#
#
##now, we make another one for the LBI values: 
#heatfiglong <-  gheatmap(p.longrheatvar,heatdat.longr[,newlabs],
#       	offset=40, font.size=8, colnames_angle = -45,
#	colnames=T, colnames_position="bottom", hjust=0.5, 
#	colnames_offset_y=-30, color=NA)+ 
#	scale_fill_continuous(name='Value of\nstandardized\nLBI',
#        low='#FEFE62',high='#5D3A9B') + 
#	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
#	coord_cartesian(clip = 'off') +
#	ggtitle('Longitudinal sampling') +
#	theme(panel.spacing=unit(0,'pt'),
#		legend.title=element_text(size=18),
#		legend.text=element_text(size=18), 
#		plot.title=element_text(hjust=0.6,size=24,face='bold') )
#
#
#ggsave(heatfiglong, file='figures/varmodelfigs/variant_and_lbi_longitudinal.pdf',height=15,width=15)
#ggsave(heatfiglong, file='figures/varmodelfigs/variant_and_lbi_longitudinal.png',dpi=600,height=15,width=15)
#                                                                               
#
## Look at cross-sectional schemes at other variant levels (1% and 5%):
#
#
#p.01r <- ggtree(tree.01, layout='rectangular')+
#	geom_tree(size=0.1, color='gray60')
#heatdat.01r <- dat.01
#
#p.05r <- ggtree(tree.05, layout='rectangular') + 
#	geom_tree(size=0.1, color='gray60')
#heatdat.05r <- dat.05
#
#p.10r <- ggtree(tree.10, layout='rectangular')+
#	geom_tree(size=0.1, color='gray60')
#heatdat.10r <- dat.10
#
#newlabs <- c('LBI_2.5')
#colnames(heatdat.01r)[colnames(heatdat.01r) %in% c('lbi.2.5')] <- newlabs
#colnames(heatdat.05r)[colnames(heatdat.05r) %in% c('lbi.2.5')] <- newlabs
#colnames(heatdat.10r)[colnames(heatdat.10r) %in% c('lbi.2.5')] <- newlabs
#
#
### 1%:
## make the variant heatmaps:
#
#vardat.01r <- heatdat.01r[,'state']
#names(vardat.01r) <- rownames(heatdat.01r)
#vardat.01r <- as.data.frame(vardat.01r)
#colnames(vardat.01r) <- 'Variant'
#
## first, we make the heatmap for variant status, which is discrete:
#heatfigvar.01r <- gheatmap(p.01r, vardat.01r, offset=08, width=0.2, color=NA,
#		colnames_offset_y = -20, hjust=0.5,
#		font.size=8,
#		colnames_position='bottom') + 
#		scale_fill_discrete(name='Variant',breaks=c('I','J'),
#			type=c('#005AB5','#DC3220'))+
#		theme(panel.spacing=unit(0,'pt'))
#
#
#
## now add the LBI stuff
#p.01rheatvar <- heatfigvar.01r + new_scale_fill()
#
#lbidat.01r <- heatdat.01r[,newlabs]
#names(lbidat.01r) <- rownames(heatdat.01r)
#lbidat.01r <- as.data.frame(lbidat.01r)
#colnames(lbidat.01r) <- newlabs
#
##now, we make another one for the LBI values: 
#heatfig.01r <-  gheatmap(p.01rheatvar,lbidat.01r, color=NA,
#       	offset=60, width=0.2,
#	font.size=8,
#	colnames=T, colnames_position="bottom", 
#	hjust=0.5, colnames_offset_y=-20)+ 
#	scale_fill_continuous(name='LBI',
#        low='#FEFE62',high='#5D3A9B') + 
#	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
#	coord_cartesian(clip = 'off') + 
#	theme(panel.spacing = unit(0,'pt')) +
#	ggtitle('Variant at 1%') +
#	theme(plot.title=element_text(hjust=0.6,size=24,face='bold')) +
#	theme(legend.title=element_text(size=18),
#	legend.text=element_text(size=18))+ theme(legend.position="none") 
#
#
#
### 5%:
## make the variant heatmaps:
#
#vardat.05r <- heatdat.05r[,'state']
#names(vardat.05r) <- rownames(heatdat.05r)
#vardat.05r <- as.data.frame(vardat.05r)
#colnames(vardat.05r) <- 'Variant'
#
## first, we make the heatmap for variant status, which is discrete:
#heatfigvar.05r <- gheatmap(p.05r, vardat.05r, offset=08, width=0.2, color=NA,
#		colnames_offset_y = -20, hjust=0.5,
#		font.size=8,
#		colnames_position='bottom') + 
#		scale_fill_discrete(name='Variant',breaks=c('I','J'),
#			type=c('#005AB5','#DC3220')) +
#		theme(panel.spacing = unit(0,'pt'))
#
#
## now add the LBI stuff
#p.05rheatvar <- heatfigvar.05r + new_scale_fill()
#
#lbidat.05r <- heatdat.05r[,newlabs]
#names(lbidat.05r) <- rownames(heatdat.05r)
#lbidat.05r <- as.data.frame(lbidat.05r)
#colnames(lbidat.05r) <- newlabs
#
##now, we make another one for the LBI values: 
#heatfig.05r <-  gheatmap(p.05rheatvar,lbidat.05r, color=NA,
#       	offset=80, width=0.2,
#	font.size=8,
#	colnames=T, colnames_position="bottom", 
#	hjust=0.5, colnames_offset_y=-20)+ 
#	scale_fill_continuous(name='LBI',
#        low='#FEFE62',high='#5D3A9B') + 
#	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
#	coord_cartesian(clip = 'off') +
#	theme(panel.spacing = unit(0,'pt')) +
#	ggtitle('Variant at 5%') +
#	theme(plot.title=element_text(hjust=0.6,size=24,face='bold')) +
#	theme(legend.title=element_text(size=18),
#	legend.text=element_text(size=18)) 
#
#
#
#
### 10%:
## make the variant heatmaps:
#
#vardat.10r <- heatdat.10r[,'state']
#names(vardat.10r) <- rownames(heatdat.10r)
#vardat.10r <- as.data.frame(vardat.10r)
#colnames(vardat.10r) <- 'Variant'
#
## first, we make the heatmap for variant status, which is discrete:
#heatfigvar.10r <- gheatmap(p.10r, vardat.10r, offset=08, width=0.2, color=NA,
#		font.size=8,
#		colnames_offset_y = -20, hjust=0.5,
#		colnames_position='bottom') + 
#		scale_fill_discrete(name='Variant',breaks=c('I','J'),
#			type=c('#005AB5','#DC3220')) + 
#		theme(panel.spacing=unit(0,'null')) +
#		ggtitle('Variant at 10%')
#
#
#
## now add the LBI stuff
#p.10rheatvar <- heatfigvar.10r + new_scale_fill()
#
#lbidat.10r <- heatdat.10r[,newlabs]
#names(lbidat.10r) <- rownames(heatdat.10r)
#lbidat.10r <- as.data.frame(lbidat.10r)
#colnames(lbidat.10r) <- newlabs
#
##now, we make another one for the LBI values: 
#heatfig.10r <-  gheatmap(p.10rheatvar,lbidat.10r, color=NA,
#       	offset=60, width=0.2,
#	colnames=T, colnames_position="bottom", 
#	font.size=8,
#	hjust=0.5, colnames_offset_y=-20)+ 
#	scale_fill_continuous(name='LBI',
#        low='#FEFE62',high='#5D3A9B') + 
#	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
#	coord_cartesian(clip = 'off') + 
#	theme(panel.spacing=unit(0,'null')) +
#	theme(plot.title=element_text(hjust=0.5,size=24,face='bold')) +
#	theme(legend.title=element_text(size=18),
#	legend.text=element_text(size=18)) 
#
#
#
#
#varfig <- plot_grid(heatfig.01r, heatfig.05r, nrow=1)
#
#ggsave(varfig, file='figures/varmodelfigs/variant_sweeps.pdf', width=15, height=12)
#ggsave(varfig, file='figures/varmodelfigs/variant_sweeps.png', dpi=600, width=15, height=12)


# Make a figure with the same basic format as the null model figures, but with just LBI in the heatmap:

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
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'I'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'I'
	crud[grep('JH',crud$label[1:ntests]),'state'] <- 'J'
	crud[grep('JL',crud$label[1:ntests]),'state'] <- 'J'


	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]-1))
	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	p <- ggtree(tree,layout='rectangular') %<+% crud

	p1 <- p + aes(col=state) + geom_tree(linewidth=0.60) +
		scale_color_manual(name='Host',
		values=c('I'='#1A85FF','J'='#D41159')) +
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
		guides(color=guide_legend(override.aes=list(linewidth=2)))
	
	#reorder the factor levels in crud$state:
	# If we want to just look at LBI at the tips, we need to just use the first
	# ntests rows of crud:
	crud$state <- factor(crud$state, levels=c('I','J'))

	p1.box <- ggplot(crud[1:ntests,]) + geom_boxplot(aes(x=state,y=lbi,fill=state)) +
		scale_fill_manual(name='Host' ,values=c('I'='#1A85FF','J'='#D41159')) + 
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


getmainplotlong <- function(tree,taulbi=4,tauthd=5,taurels=6,tauclust=6,title='title'){

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
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'I'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'I'
	crud[grep('JH',crud$label[1:ntests]),'state'] <- 'J'
	crud[grep('JL',crud$label[1:ntests]),'state'] <- 'J'


	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]-1))
	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	p <- ggtree(tree,layout='rectangular') %<+% crud

	p1 <- p + aes(col=state) + geom_tree(linewidth=0.60) +
		scale_color_manual(name='Host',
		values=c('I'='#1A85FF','J'='#D41159')) +
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
		guides(color=guide_legend(override.aes=list(linewidth=2)))
	
	#reorder the factor levels in crud$state:
	# If we want to just look at LBI at the tips, we need to just use the first
	# ntests rows of crud:
	crud$state <- factor(crud$state, levels=c('I','J'))

	p1.points <- ggplot(crud[1:ntests,]) + geom_point(aes(x=time,y=lbi,color=state)) +
		scale_color_manual(name='Host' ,values=c('I'='#1A85FF','J'='#D41159')) + 
		theme_classic() + 
		ylab('LBI (raw)') + 
		xlab('Time') + 
		theme(legend.position='none') +
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold'))
	

	# make the tree and pointsplot figures w/out the title first:
	alnd <- align_plots(plin4,p1.points,align='v',axis='lr')
	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1, rel_heights=c(2,1))



        return(list(fig.p1,dat))
}


fig.10 <- getmainplot(tree.10, title='Variant @ 10%', taulbi=1)
fig.50 <- getmainplot(tree.50, title='Variant @ 50%', taulbi=15)
fig.90 <- getmainplot(tree.90, title='Variant @ 90%', taulbi=25)

figcross <- plot_grid(fig.10[[1]],fig.50[[1]],fig.90[[1]], nrow=1)

ggsave(figcross, file='figures/varmodelfigs/varmodeltrees.png', dpi=600, height=16, width=24)

fig.long <- getmainplotlong(tree.long, title='Longitudinal observations', taulbi=20)[[1]]
ggsave(fig.long, file='figures/varmodelfigs/varmodeltrees_long.png', dpi=600, height=16, width=24)



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
	crud[grep('IH',crud$label[1:ntests]),'state'] <- 'I'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- 'I'
	crud[grep('JH',crud$label[1:ntests]),'state'] <- 'J'
	crud[grep('JL',crud$label[1:ntests]),'state'] <- 'J'


	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]+0))

	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	ratios <- apply(crud[,paste0('lbi',tauvals)], 2, function(x) mean(x[crud$state=='J'])/mean(x[crud$state=='I'])  )

	plot(tauvals, ratios, xlab=bquote(tau),ylab='LBI ratio (J/I)',  
		type='l',lwd=3, ylim=c(0, max(ratios)*1.5),
		cex.axis=1.3,cex.lab=1.3, main=title)
}


pdf(file='figures/varmodelfigs/var10.pdf',height=6,width=6)
gettauplot(tree.10, title='Variant @ 10%')
dev.off()

pdf(file='figures/varmodelfigs/var50.pdf',height=6,width=6)
gettauplot(tree.50, title='Variant @ 50%')
dev.off()

pdf(file='figures/varmodelfigs/var90.pdf',height=6,width=6)
gettauplot(tree.90, title='Variant @ 90%')
dev.off()




