rm(list=ls())

source("compile_variant_and_superspreading_model.R")

#############
#############

# How long to simulate? 
time <- c(-125,50)

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


theta <- gettheta(0.01,0.9)

# set the timestep size:
dT <- 1

# specify initial states:
initialStates <- c(Npop,0,1,0,0,0,0)
names(initialStates) <- c('S','L','IH','IL','M','JH','JL') 

# the save functionality in simualte_tree seems buggy, but
# the write.tree function seems to work just fine.
getsimtree <- function(tms,output,pH){
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

# How many simulations do we want? 
numouts <- 30


# produce the trajectory from which to obtain a tree - 
#	where should we save it?

for(outind in 1:numouts){
 outpath <- paste0('sims/nullmodel/out',outind,'.txt')

 # run the sim. These files are huge (~500 Mb), so we'll just
 # save the seeds used to run these
 out <-  sir_simu(
    paramValues = as.list(theta),
    initialStates = initialStates,
    tau = .0001,
    times = time,
    method = "mixed",
    verbose = TRUE,
    nTrials = 100
    )
 #Save the seed:
 writeLines(as.character(out$seed), paste0('sims/nullmodel/seed',outind,'.txt'))

 # Simulate trees for different combinations of ntests and Nyr.

 # How many tests, and how long to test?
 ntestsvec <- c(500,2000)
 Nyrvec <- c(5,10,15,20)


 # Simulate a tree for each combo of ntests and Nyr:
 for(ntests in ntestsvec){
  for(Nyr in Nyrvec){
   # Where should we save our simulated tree?
   treepath <- paste0('sims/nullmodel/',ntests,'/',Nyr,'yr/tree',outind,'.nwk')
   # produce a simulated tree:
   source('tedit_null_bigsim_gettrees.R')
  }
 }

 # close loop, run another trajectory
}

