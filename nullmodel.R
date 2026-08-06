rm(list=ls())

# This code runs simulations of the null model (A) to generate Figure 3, 
# Figure S3, and Figure S4.

source("compile_variant_and_superspreading_model.R")

#############
#############

ntests <- 84

time <- c(-150,50)
#time <- c(-750,50)

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

# Produce trees under homogeneous transmission first:
tree1 <- getsimtree(tms.lin4,out1)

tms2.lin4 = runif(-50,max(tms.lin4),n=ntests)
tree2 <- getsimtree(tms2.lin4,out1)

# Superspreading trees next:
tree3 <- getsimtree(tms.lin4,out2)
tree4 <- getsimtree(tms2.lin4,out2)


# load the LBI function:
source('lbi.R')



