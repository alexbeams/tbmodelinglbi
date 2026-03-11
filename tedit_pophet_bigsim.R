rm(list=ls())

#compile the model:
source("compile_pophet_and_superspreading_model.R")

# load in the local branching index (LBI) function
source('lbi.R')

#############
#############

time <- c(0,200)

# what is the total population size?
Npop <- 2.0 * 10^5

# specify the level of preferential mixing (mii, mij), the different
# levels of susceptibility across the groups (zeta1, zeta2), the levels
# of infectiousness across the two grups (theta1, theta2), and the 
# level of superspreading functioning independently of all of this (fH,pH):

gettheta <- function(mii,mij,zeta1,zeta2,theta1,theta2,fH=0.9,pH=0.1){
	
	#of those with active TB, what fraction are superspreaders?
	#pH <- 0.1
	# infectiousness multipliers
	#fH <- 0.9 # fraction of transmission from superspreaders

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
	#zeta1 <- zeta1
	#zeta2 <- zeta2

	# infectiousness
	thetaL <- (1-fH)/(1-pH)
	thetaH <- fH/pH


	# contact rates within/across groups:
	# make c12 and c22 larger to model increased
	# infectiousness from those groups indep of superspreading: 
	c <- beta
	c11 <- beta
	c12 <- beta
	c21 <- beta
	c22 <- beta

	cij <- matrix(c(c11,c12,c21,c22),byrow=T,nrow=2)
	# make the off-diagonals smaller for preferential mixing:
	diag(cij) <- mii * diag(cij)
	cij[1,2] <- mij * cij[1,2]
	cij[2,1] <- mij * cij[2,1]
	#cij <- cij / max(eigen(cij)$values)
	#cij <- cij*c

	c11 <- cij[1,1]	
	c12 <- cij[1,2]	
	c21 <- cij[2,1]	
	c22 <- cij[2,2]	

	# can make c12 and c22 larger to model increased
	# infectiousness from group 2: 
	c11 <- c11 * theta1
	c21 <- c21 * theta1
	c12 <- c12 * theta2
	c22 <- c22 * theta2

	theta <- list(
		beta=beta,
		zeta1 = zeta1,
		zeta2 = zeta2,
		theta1 = theta1,
		theta2 = theta2,
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


## We will simulate with a large difference in susceptibility, 
## and look at a range of preferential mixing parameters

# 1. No preferential mixing
theta1 <- gettheta(mii=1,mij=1,zeta1=1/4,zeta2=4,
	theta1=1,theta2=1,pH=0.1,fH=0.9)

# 2. Some preferential mixing
theta2 <- gettheta(mii=5,mij=.5,zeta1=1/4,zeta2=4,
	theta1=1,theta2=1,pH=0.1,fH=0.9)

# 3. Strong preferential mixing
theta3 <- gettheta(mii=5,mij=.1,zeta1=1/4,zeta2=4,
	theta1=1,theta2=1,pH=0.1,fH=0.9)

# 4. Extreme preferential mixing
theta4 <- gettheta(mii=5,mij=.01,zeta1=1/4,zeta2=4,
	theta1=1,theta2=1,pH=0.1,fH=0.9)


# Set the up the simulation time and the initial conditions:
# set the timestep size:
dT <- 1

# specify initial conditions -- we need to be careful about this now to make
# sure that we have the correct proportion of superspreaders for the
# various cases if we change proportions above:

# Specify initial states
initialStates <- c(Npop * 0.8,0,0,0,Npop * 0.2,0,1,0)
names(initialStates) <- c('S1','L1','IH1','IL1','S2','L2','IH2','IL2') 


# Create a function to simulate trees. The save feature in simulate_tree
# seems buggy, but write.tree works just fine for our purposes

getsimtree <- function(tms,output){
	simulate_tree(
	simuResults=output,
	dates=tms,
	deme=c('IH1','IL1','L1','IH2','IL2','L2'),
	root = 'IH2', # set this to match the initial condition 
	nTrials=50,
	resampling=FALSE,
	addInfos = TRUE)
}

out1 <-  sir_simu(
   paramValues = as.list(theta1),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)

out2 <-  sir_simu(
   paramValues = as.list(theta2),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)

out3 <-  sir_simu(
   paramValues = as.list(theta3),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)

out4 <-  sir_simu(
   paramValues = as.list(theta4),
   initialStates = initialStates,
   tau = .0001,
   times = time,
   method = "mixed",
   verbose = TRUE,
   nTrials = 100)


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

traj3 <- out3$traj
plottraj3 <- traj3[c(1:1500,sample(1:dim(traj3)[1],1000)),]
plottraj3 <- plottraj3[order(plottraj3$Time),]
# transform time back to date for plotting:
plottraj3$date <- as.Date(plottraj3$Time * 365)

traj4 <- out4$traj
plottraj4 <- traj4[c(1:1500,sample(1:dim(traj4)[1],1000)),]
plottraj4 <- plottraj4[order(plottraj4$Time),]
# transform time back to date for plotting:
plottraj4$date <- as.Date(plottraj4$Time * 365)



# for cross-sectional sampling similar to lineage 4:
tmstart <- max(time)-5
tmend <- max(time) 
ntests <- 500


x1 <- plottraj1[plottraj1$Time < tmend & plottraj1$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]
x2 <- plottraj2[plottraj2$Time < tmend & plottraj2$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]
x3 <- plottraj3[plottraj3$Time < tmend & plottraj3$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]
x4 <- plottraj4[plottraj4$Time < tmend & plottraj4$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]

y1 <- colSums(x1)/sum(colSums(x1))
y2 <- colSums(x2)/sum(colSums(x2))
y3 <- colSums(x3)/sum(colSums(x3))
y4 <- colSums(x4)/sum(colSums(x4))

# randomly sample states according to the probabilities with which they were observed in the simulation:
s1 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y1[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]
s2 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y2[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]
s3 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y3[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]
s4 <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y4[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]

# create data frames with the tms.lin4 sampling times and the sampled state:
tms.1 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s1) 
tms.2 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s2) 
tms.3 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s3) 
tms.4 <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s4) 

tms.1 <- tms.1[order(tms.1$Date),]
tms.2 <- tms.2[order(tms.2$Date),]
tms.3 <- tms.3[order(tms.3$Date),]
tms.4 <- tms.4[order(tms.4$Date),]


# Make sure that pop2frac here matches the initial conditions above
tree1 <- getsimtree(tms.1,out1)
tree2 <- getsimtree(tms.2,out2)
tree3 <- getsimtree(tms.3,out3)
tree4 <- getsimtree(tms.4,out4)


