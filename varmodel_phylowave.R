rm(list=ls())

source("compile_variant_and_superspreading_model.R")
source('lbi.R')

require(ggnewscale)


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


# Sample longitudinally over the whole sim:
tree.long <- getsimtree.long(tms.long)

# Make a figure with the same basic format as the null model figures:

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
	crud[grep('IH',crud$label[1:ntests]),'state'] <- '1'
	crud[grep('IL',crud$label[1:ntests]),'state'] <- '1'
	crud[grep('JH',crud$label[1:ntests]),'state'] <- '2'
	crud[grep('JL',crud$label[1:ntests]),'state'] <- '2'


	nodenms <- sapply(crud[(ntests+1):(ntests+tree$Nnode),'label'], function(z) substr(z, regexpr("S+",z)[1]+2, regexpr(".[+]=",z)[1]-1))
	nodenms[nodenms=='I'] <- 1
	nodenms[nodenms=='J'] <- 2
	crud[(Ntip(tree)+1):(tree$Nnode + Ntip(tree)),'state'] <- nodenms

	p <- ggtree(tree,layout='rectangular') %<+% crud

	p1 <- p + aes(col=state) + geom_tree(linewidth=0.60) +
		scale_color_manual(name='Subtype',
		values=c('1'='#1A85FF','2'='#D41159')) +
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
	crud$state <- factor(crud$state, levels=c('1','2'))

	p1.points <- ggplot(crud[1:ntests,]) + geom_point(aes(x=time,y=lbi,color=state)) +
		scale_color_manual(name='Subtype' ,values=c('1'='#1A85FF','2'='#D41159')) + 
		scale_x_continuous(labels = function(x) round(x + 1970)) +
		theme_classic() + 
		ylab('LBI (raw)') + 
		xlab('Year') + 
		theme(legend.position='none') +
		theme(axis.text=element_text(size=18), 
			axis.title=element_text(size=18, face='bold'))
	

	# make the tree and pointsplot figures w/out the title first:
	alnd <- align_plots(plin4,p1.points,align='v',axis='lr')
	fig.p1 <- plot_grid(alnd[[1]],alnd[[2]],ncol=1, rel_heights=c(2,1))



        return(list(fig.p1,dat))
}

fig.long <- getmainplotlong(tree.long, title='Longitudinal observations', taulbi=20)[[1]]
ggsave(fig.long, file='figures/varmodelfigs/FigS6.png', dpi=300, height=8, width=12)


