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



p.01 <- ggtree(tree.01,layout='circular')
p.05 <- ggtree(tree.05,layout='circular')
p.10 <- ggtree(tree.10,layout='circular')

# calculate data for tree.01:
tree <- tree.01
x <- cophenetic(tree)
tau <- mean(x)*.1
dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
dat <- as.data.frame(dat)
colnames(dat) <- 'thd'
dat$label = rownames(dat)
dat <- dat[,c(2,1)]
dat$var <- apply(x,1,var)
dat$min <- apply(x,1,function(x) min(x[x>0]) )
dat$lbi.1000 <- lbi(tree,tau=1000)[1:length(tree$tip.label)]
dat$lbi.500 <- lbi(tree,tau=500)[1:length(tree$tip.label)]
dat$lbi.100 <- lbi(tree,tau=100)[1:length(tree$tip.label)]
dat$lbi.50 <- lbi(tree,tau=50)[1:length(tree$tip.label)]
dat$lbi.10 <- lbi(tree,tau=10)[1:length(tree$tip.label)]
dat$lbi.1 <- lbi(tree,tau=1)[1:length(tree$tip.label)]
# now add in states:
dat$state <- 0
dat[grep('I',tree$tip.label),'state'] <- 'I'
dat[grep('J',tree$tip.label),'state'] <- 'J'
dat.01 <- dat

# merge data onto the ggtree plot:
p.01 <- p.01 %<+% dat.01


# calculate data for tree.05:
tree <- tree.05
x <- cophenetic(tree)
tau <- mean(x)*.1
dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
dat <- as.data.frame(dat)
colnames(dat) <- 'thd'
dat$label = rownames(dat)
dat <- dat[,c(2,1)]
dat$var <- apply(x,1,var)
dat$min <- apply(x,1,function(x) min(x[x>0]) )
dat$lbi.1000 <- lbi(tree,tau=1000)[1:length(tree$tip.label)]
dat$lbi.500 <- lbi(tree,tau=500)[1:length(tree$tip.label)]
dat$lbi.100 <- lbi(tree,tau=100)[1:length(tree$tip.label)]
dat$lbi.50 <- lbi(tree,tau=50)[1:length(tree$tip.label)]
dat$lbi.10 <- lbi(tree,tau=10)[1:length(tree$tip.label)]
dat$lbi.1 <- lbi(tree,tau=1)[1:length(tree$tip.label)]
# now add in states:
dat$state <- 0
dat[grep('I',tree$tip.label),'state'] <- 'I'
dat[grep('J',tree$tip.label),'state'] <- 'J'
dat.05 <- dat
p.05 <- p.05 %<+% dat.05

# calculate data for tree.10:
tree <- tree.10
x <- cophenetic(tree)
tau <- mean(x)*.1
dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
dat <- as.data.frame(dat)
colnames(dat) <- 'thd'
dat$label = rownames(dat)
dat <- dat[,c(2,1)]
dat$var <- apply(x,1,var)
dat$min <- apply(x,1,function(x) min(x[x>0]) )
dat$lbi.1000 <- lbi(tree,tau=1000)[1:length(tree$tip.label)]
dat$lbi.500 <- lbi(tree,tau=500)[1:length(tree$tip.label)]
dat$lbi.100 <- lbi(tree,tau=100)[1:length(tree$tip.label)]
dat$lbi.50 <- lbi(tree,tau=50)[1:length(tree$tip.label)]
dat$lbi.10 <- lbi(tree,tau=10)[1:length(tree$tip.label)]
dat$lbi.1 <- lbi(tree,tau=1)[1:length(tree$tip.label)]
# now add in states:
dat$state <- 0
dat[grep('I',tree$tip.label),'state'] <- 'I'
dat[grep('J',tree$tip.label),'state'] <- 'J'
dat.10 <- dat
p.10 <- p.10 %<+% dat.10

rm(dat)
rm(x)
rm(tau)


p.10r <- ggtree(tree.10, layout='rectangular') +
	geom_tree(size=0.1,color='gray60')
heatdat.10r <- dat.10

# standardize LBI values for plotting:
z.lbi.1 <- heatdat.10r$lbi.1
z.lbi.1 <- z.lbi.1-mean(z.lbi.1)
z.lbi.1 <- z.lbi.1/sd(z.lbi.1)

z.lbi.10 <- heatdat.10r$lbi.10
z.lbi.10 <- z.lbi.10-mean(z.lbi.10)
z.lbi.10 <- z.lbi.10/sd(z.lbi.10)

z.lbi.50 <- heatdat.10r$lbi.50
z.lbi.50 <- z.lbi.50-mean(z.lbi.50)
z.lbi.50 <- z.lbi.50/sd(z.lbi.50)

z.lbi.100 <- heatdat.10r$lbi.100
z.lbi.100 <- z.lbi.100-mean(z.lbi.100)
z.lbi.100 <- z.lbi.100/sd(z.lbi.100)

z.lbi.1000 <- heatdat.10r$lbi.1000
z.lbi.1000 <- z.lbi.1000-mean(z.lbi.1000)
z.lbi.1000 <- z.lbi.1000/sd(z.lbi.1000)

heatdat.10r$z.lbi.1 <- z.lbi.1
heatdat.10r$z.lbi.10 <- z.lbi.10
heatdat.10r$z.lbi.50 <- z.lbi.50
heatdat.10r$z.lbi.100 <- z.lbi.100
heatdat.10r$z.lbi.1000 <- z.lbi.1000

p.10r <- p.10r %<+% heatdat.10r


#new_labels <- c(expression(tau~"=1"),'10','50','100','1000' )

newlabs <- c('LBI_1','LBI_10','LBI_50','LBI_100','LBI_1000')

colnames(heatdat.10r)[colnames(heatdat.10r) %in% c(paste0('z.lbi.',c(1,10,50,100,1000)))] <- newlabs

vardat <- heatdat.10r[,'state']
names(vardat) <- rownames(heatdat.10r)
vardat <- as.data.frame(vardat)
colnames(vardat) <- 'Variant'

# first, we make the heatmap for variant status, which is discrete:
heatfigvar <- gheatmap(p.10r, vardat, offset=0, width=0.2, colnames_angle=-45,
		colnames_offset_y = -30, hjust=0.5, font.size=8,
		colnames_position='bottom', color=NA) + 
		scale_fill_discrete(name='Variant',breaks=c('I','J'),
			type=c('#005AB5','#DC3220'))  + 
		theme(panel.spacing = unit(0,'pt'),legend.title =element_text(size=18),
			legend.text=element_text(size=18))

#heatfigvar <- ggplot(heatdat.10r, aes(x=c(1),y=label)) +
#	geom_tile(aes(fill=state)) + xlab(NULL) + ylab(NULL)
 
library(ggnewscale)
p.10rheatvar <- heatfigvar + new_scale_fill()


#now, we make another one for the LBI values: 
heatfig <-  gheatmap(p.10rheatvar,heatdat.10r[,newlabs],
       	offset=40, font.size=8, colnames_angle = -45,
	colnames=T, colnames_position="bottom", hjust=0.5, 
	colnames_offset_y=-30, color=NA)+ 
	scale_fill_continuous(name='Value of\nstandardized\nLBI',
        low='#FEFE62',high='#5D3A9B') + 
	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
	coord_cartesian(clip = 'off') +
	ggtitle('Variant at 10%') +
	theme(panel.spacing=unit(0,'pt'),
		legend.title=element_text(size=18),
		legend.text=element_text(size=18), 
		plot.title=element_text(hjust=0.6,size=24,face='bold') )

ggsave(heatfig, file='varmodelfigs/variant_and_lbi.pdf',height=15,width=15)
ggsave(heatfig, file='varmodelfigs/variant_and_lbi.png',dpi=600,height=15,width=15)
                                                                               



# make figure 9:

p.01r <- ggtree(tree.01, layout='rectangular')+
	geom_tree(size=0.1, color='gray60')
heatdat.01r <- dat.01

p.05r <- ggtree(tree.05, layout='rectangular') + 
	geom_tree(size=0.1, color='gray60')
heatdat.05r <- dat.05

p.10r <- ggtree(tree.10, layout='rectangular')+
	geom_tree(size=0.1, color='gray60')
heatdat.10r <- dat.10

newlabs <- c('LBI_50')
colnames(heatdat.01r)[colnames(heatdat.01r) %in% c('lbi.50')] <- newlabs
colnames(heatdat.05r)[colnames(heatdat.05r) %in% c('lbi.50')] <- newlabs
colnames(heatdat.10r)[colnames(heatdat.10r) %in% c('lbi.50')] <- newlabs


## 1%:
# make the variant heatmaps:

vardat.01r <- heatdat.01r[,'state']
names(vardat.01r) <- rownames(heatdat.01r)
vardat.01r <- as.data.frame(vardat.01r)
colnames(vardat.01r) <- 'Variant'

# first, we make the heatmap for variant status, which is discrete:
heatfigvar.01r <- gheatmap(p.01r, vardat.01r, offset=08, width=0.2, color=NA,
		colnames_offset_y = -20, hjust=0.5,
		font.size=8,
		colnames_position='bottom') + 
		scale_fill_discrete(name='Variant',breaks=c('I','J'),
			type=c('#005AB5','#DC3220'))+
		theme(panel.spacing=unit(0,'pt'))



# now add the LBI stuff
p.01rheatvar <- heatfigvar.01r + new_scale_fill()

lbidat.01r <- heatdat.01r[,newlabs]
names(lbidat.01r) <- rownames(heatdat.01r)
lbidat.01r <- as.data.frame(lbidat.01r)
colnames(lbidat.01r) <- newlabs

#now, we make another one for the LBI values: 
heatfig.01r <-  gheatmap(p.01rheatvar,lbidat.01r, color=NA,
       	offset=60, width=0.2,
	font.size=8,
	colnames=T, colnames_position="bottom", 
	hjust=0.5, colnames_offset_y=-20)+ 
	scale_fill_continuous(name='LBI',
        low='#FEFE62',high='#5D3A9B') + 
	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
	coord_cartesian(clip = 'off') + 
	theme(panel.spacing = unit(0,'pt')) +
	ggtitle('Variant at 1%') +
	theme(plot.title=element_text(hjust=0.6,size=24,face='bold')) +
	theme(legend.title=element_text(size=18),
	legend.text=element_text(size=18))+ theme(legend.position="none") 



## 5%:
# make the variant heatmaps:

vardat.05r <- heatdat.05r[,'state']
names(vardat.05r) <- rownames(heatdat.05r)
vardat.05r <- as.data.frame(vardat.05r)
colnames(vardat.05r) <- 'Variant'

# first, we make the heatmap for variant status, which is discrete:
heatfigvar.05r <- gheatmap(p.05r, vardat.05r, offset=08, width=0.2, color=NA,
		colnames_offset_y = -20, hjust=0.5,
		font.size=8,
		colnames_position='bottom') + 
		scale_fill_discrete(name='Variant',breaks=c('I','J'),
			type=c('#005AB5','#DC3220')) +
		theme(panel.spacing = unit(0,'pt'))


# now add the LBI stuff
p.05rheatvar <- heatfigvar.05r + new_scale_fill()

lbidat.05r <- heatdat.05r[,newlabs]
names(lbidat.05r) <- rownames(heatdat.05r)
lbidat.05r <- as.data.frame(lbidat.05r)
colnames(lbidat.05r) <- newlabs

#now, we make another one for the LBI values: 
heatfig.05r <-  gheatmap(p.05rheatvar,lbidat.05r, color=NA,
       	offset=80, width=0.2,
	font.size=8,
	colnames=T, colnames_position="bottom", 
	hjust=0.5, colnames_offset_y=-20)+ 
	scale_fill_continuous(name='LBI',
        low='#FEFE62',high='#5D3A9B') + 
	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
	coord_cartesian(clip = 'off') +
	theme(panel.spacing = unit(0,'pt')) +
	ggtitle('Variant at 5%') +
	theme(plot.title=element_text(hjust=0.6,size=24,face='bold')) +
	theme(legend.title=element_text(size=18),
	legend.text=element_text(size=18)) 




## 10%:
# make the variant heatmaps:

vardat.10r <- heatdat.10r[,'state']
names(vardat.10r) <- rownames(heatdat.10r)
vardat.10r <- as.data.frame(vardat.10r)
colnames(vardat.10r) <- 'Variant'

# first, we make the heatmap for variant status, which is discrete:
heatfigvar.10r <- gheatmap(p.10r, vardat.10r, offset=08, width=0.2, color=NA,
		font.size=8,
		colnames_offset_y = -20, hjust=0.5,
		colnames_position='bottom') + 
		scale_fill_discrete(name='Variant',breaks=c('I','J'),
			type=c('#005AB5','#DC3220')) + 
		theme(panel.spacing=unit(0,'null')) +
		ggtitle('Variant at 10%')



# now add the LBI stuff
p.10rheatvar <- heatfigvar.10r + new_scale_fill()

lbidat.10r <- heatdat.10r[,newlabs]
names(lbidat.10r) <- rownames(heatdat.10r)
lbidat.10r <- as.data.frame(lbidat.10r)
colnames(lbidat.10r) <- newlabs

#now, we make another one for the LBI values: 
heatfig.10r <-  gheatmap(p.10rheatvar,lbidat.10r, color=NA,
       	offset=60, width=0.2,
	colnames=T, colnames_position="bottom", 
	font.size=8,
	hjust=0.5, colnames_offset_y=-20)+ 
	scale_fill_continuous(name='LBI',
        low='#FEFE62',high='#5D3A9B') + 
	theme(plot.margin=unit(c(1,1,2,1),'cm')) + 
	coord_cartesian(clip = 'off') + 
	theme(panel.spacing=unit(0,'null')) +
	theme(plot.title=element_text(hjust=0.5,size=24,face='bold')) +
	theme(legend.title=element_text(size=18),
	legend.text=element_text(size=18)) 




varfig <- plot_grid(heatfig.01r, heatfig.05r, nrow=1)

ggsave(varfig, file='varmodelfigs/variant_sweeps.pdf', width=15, height=12)
ggsave(varfig, file='varmodelfigs/variant_sweeps.png', dpi=600, width=15, height=12)




