rm(list=ls())

# This simulates the pophet model under Case 1: The two host sub-populations
# are equal in size, but one group is more susceptible (and therefore over-represented
# among infections that get sampled).

# The key parameters that get varied will be 
# 1. y = L2/(L1+L2) = I2/(I1+I2), the level of over-representation
#    of the Group 2 individuals among infections. We care about 0.5 < y < 1:

# 2. w, the parameter describing the level of mixing between the two groups: 
#    0 < w < 1, with w=1 homogeneous mixing and w=0 no cross-group transmission

library(ggnewscale)

#compile the model:
source("compile_pophet_and_superspreading_model.R")

# load in the local branching index (LBI) function
source('lbi.R')

#############
#############

ntests <- 500

# Time to simulate:
time <- c(-125,50)

# when should we start to sample?
tm1 <- 45
tms.lin4 <- runif(ntests,min=tm1,max=time[2])
tms.lin4 <- tms.lin4[order(tms.lin4)]


# what is the total population size?
Npop <- 2.0 * 10^5

# specify the level of preferential mixing, w, and the representation of 
# Group 2 individuals among infections (0.5<y<1 modeling over-representation of Group 2):
# and the level of superspreading functioning independently of all of this (fH,pH):

gettheta <- function(y,w,fH=0.9,pH=0.1){
	
	#of those with active TB, what fraction are superspreaders?
	#pH <- 0.1
	# infectiousness multipliers
	#fH <- 0.9 # fraction of transmission from superspreaders

	N <- Npop

	# parameterize in terms of equilibrium values of S, L, I
	Seq <- 2/3 * Npop
	Ieq <- 0.01 * Npop
	Leq <- Npop-Seq-Ieq

	# what proportion of newly infected immediately become infectious?
	# 5% develop active TB in the first two years; use this:
	pfast <- 0.05
	# This variable corresponds to "p" in our analysis from Maple that gives us zeta, etc.
	p <- pfast

	# Average duration of untreated pulmonary TB, historically, is ~ 1- 3 years:
	gamma <- 1/1.5
	# The solutions fro Maple use this:	
	gam <- gamma

	# The annual risk of developing active TB is 1e-4 - 2e-4 per year:
	sigma <- 1e-4
	 
	# transmission rates (beta=zeta*c*theta):
	beta <- (gamma*Ieq-sigma*Leq)/pfast/Seq/Ieq

	# determine rho:
	#rho <- (beta*Seq*Ieq - gamma*Ieq)/Leq
	rho <- ((1-p)*gamma*Ieq - sigma*Leq)/(p*Leq)

	# What is the infectiousness of group 2 relative to group 1?
	theta <- 1

	# What is the susceptibility of Group 2, given y=L2/(L1+L2) and a value of w, where w <1
	# describes reduced cross-group contact and w =1 is homogeneous mixing between the groups?
	zeta <- -y*((2*Ieq + 2*Leq)*y + N - 2*Ieq - 2*Leq)/(((-2*Ieq - 2*Leq)*y + N)*(((theta + 1)*w - theta)*y - w))

	# What is the baseline transmission rate, beta, in the perfect mixing case, given
	# y=L2/(L1+L2)?
	beta0 <- -2*(Ieq*gam*y - Leq*sigma*y - Ieq*gam + Leq*sigma)/(Ieq*(2*Ieq*y + 2*Leq*y - 2*Ieq - 2*Leq + N)*p*(theta*y - y + 1)) 

	# If betaij = w * beta0, so that 0 < w < 1 models assortative mixing, what is betaii?
	betaii <- -2*(-1 + (1 + (w - 1)*theta)*y)*(Ieq*gam - Leq*sigma)/(Ieq*(1 + (theta - 1)*y)*((2*Ieq + 2*Leq)*y + N - 2*Ieq - 2*Leq)*p)

	# Baseline superspreading parameters: infectiousness
	thetaL <- (1-fH)/(1-pH)
	thetaH <- fH/pH

	# specify the transmission rates:
	c11 <- betaii
	c12 <- w * beta0
	c21 <- w * beta0
	c22 <- betaii
	
	theta <- list(
		beta=beta,
		zeta1 = 1,
		zeta2 = zeta,
		theta1 = 1,
		theta2 = theta,
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

# create a function to simulate trees given tms and output:
getsimtree <- function(tms,output){
	simulate_tree(
	simuResults=output,
	dates=tms,
	deme=c('IH1','IL1','L1','IH2','IL2','L2'),
	root = 'IH2', # set this to match the initial condition used to produce simuResults=output
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)
}


# 1. No preferential mixing, heightened susceptibility for Group 2 (they are 80% of all infs):
theta1 <- gettheta(y=0.8,w=1,pH=0.1,fH=0.9)

# 2. Strong preferential mixing, "...":
theta2 <- gettheta(y=0.8,w=.1,pH=0.1,fH=0.9)


# 1. Both groups are equal size in the overall population:
initialStates <- c(Npop * 0.5,0,0,0,Npop * 0.5,0,1,0)
names(initialStates) <- c('S1','L1','IH1','IL1','S2','L2','IH2','IL2') 

# Folder to save outputs:
case <- 'case1/'

# First, simulate using theta1 (homogeneous mixing):
numouts <- 30
siminds <- 1:numouts

# sub-folder to save outputs:
thetafolder <- 'theta1/'

for(outind in siminds){
 outpath <- paste0('sims/pophetmodel/case1/theta1/',outind,'.txt')

 # Run the sim, save the seed!
 out <-  sir_simu(
    paramValues = as.list(theta1),
    initialStates = initialStates,
    tau = .0001,
    times = time,
    method = "mixed",
    verbose = TRUE,
    nTrials = 100)
 # Save the seed:
 writeLines(as.character(out$seed), paste0('sims/pophetmodel/case1/theta1/seed',outind,'.txt'))

 # Simulate a tree :
 source('tedit_pophet_case1_bigsim_gettrees.R')
} 


# Second, simulate using theta2 (preferential mixing):
numouts <- 30
siminds <- 1:numouts

# sub-folder to save outputs:
thetafolder <- 'theta2/'


for(outind in siminds){
 outpath <- paste0('sims/pophetmodel/case1/theta2/',outind,'.txt')

 # Run the sim, save the seed!
 out <-  sir_simu(
    paramValues = as.list(theta2),
    initialStates = initialStates,
    tau = .0001,
    times = time,
    method = "mixed",
    verbose = TRUE,
    nTrials = 100)
 # Save the seed:
 writeLines(as.character(out$seed), paste0('sims/pophetmodel/case1/theta2/seed',outind,'.txt'))

 # Simulate a tree :
 source('tedit_pophet_case1_bigsim_gettrees.R')
} 






## plotting function:

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
		guides(color=guide_legend(override.aes=list(linewidth=2)))
	
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


