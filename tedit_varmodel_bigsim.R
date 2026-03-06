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

#########
######### Set up the epidemiological parameters:
#########

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
pmu <- 0.05


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

# set the timestep size:
dT <- 1

# specify initial states:
initialStates <- c(Npop,0,1,0,0,0,0)
names(initialStates) <- c('S','L','IH','IL','M','JH','JL') 

# Create a function to simulae a tree from a realization of
# the compartmental model - it depends on tms, output, and p:

# the save functionality in simulate_tree seems buggy,
# but the write.tree function seems to work just fine.

getsimtree.p <- function(tms,output,p){
	simulate_tree(
	simuResults=output,
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


# How many simulations do we want?
numouts <- 30


# produce the trajectory from which to obtain trees - where to save?

# Not all of the these will succeed at generating trees -- for now, just
# going in to manually change siminds to re-run the small number that fail
# under a different random seed until we get 30 simulations 

#siminds <- 1:numouts
#siminds <- c(1,4,7,11,12,17,18,19,21,22,23,25,26,27)
#siminds <- c(1,4,7,11,12,18,22,23,26,27)
siminds <- 12

for(outind in siminds){
 outpath <- paste0('sims/varmodel/out',outind,'.txt')

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
 writeLines(as.character(out$seed), paste0('sims/varmodel/seed',outind,'.txt'))

 # Simulate trees for different variant frequencies from 10% -- 50%.
 source('tedit_varmodel_bigsim_gettrees.R')

 # close loop, run another trajectory
}


