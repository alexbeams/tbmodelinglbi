# when should we start to sample?
tm1 <- 45
tms.lin4 <- runif(ntests,min=tm1,max=time[2])
tms.lin4 <- tms.lin4[order(tms.lin4)]
# ^^ this gets overwritten below; we are sampling over a 5-yr period

# Calculate the sampling times for the trajectory.
# Start by looking at variant frequency over time:

# Extract the trajectory record:
traj <- out$traj

# calculate the proportion of infections caused by the variant:
traj$p <- with(traj, (JH+JL)/(IH+IL+JH+JL) )

## If the variant doesn't show up and take over, just toss out the simulation
## and try again, otherwise run the code to simulate trees at different variant frequencies:

if(max(traj$p) > 0.6){

	# it's a huge object, so we need to subsample it:
	plottraj <- traj[c(1:1500,sample(1:dim(traj)[1],1000)),]
	plottraj <- plottraj[order(plottraj$Time),]


	# first, let's figure out when the variant reached x% of all infections:
	t.10 <- traj[max(which(traj$p <= .10)), 'Time']
	t.20 <- traj[max(which(traj$p <= .20)), 'Time']
	t.30 <- traj[max(which(traj$p <= .30)), 'Time']
	t.40 <- traj[max(which(traj$p <= .40)), 'Time']
	t.50 <- traj[max(which(traj$p <= .50)), 'Time']

	# let's use the same spacing of sampling events, but shift them so the
	#	first sampling event coincides with t.x

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

	# set the minimum time for longitudinal sampling to commence:
	mintime <- -120

	# For longitudinal sampling, we would want to account for
	# changing variant frequency over time, but over a 5-yr period
	# it isn't necessary -- can uncomment to change if needed 

	#crud <- plottraj[plottraj$Time > mintime,]
	#crud <- crud[sort(sample(1:dim(crud)[1],ntests)),]
	#crud$pih <- (1-crud$p) * pH
	#crud$pil <- (1-crud$p) * (1-pH)
	#crud$pjh <- crud$p * pH
	#crud$pjl <- crud$p * (1-pH)
	#
	#crud$sampled <- apply(crud, 1, function(x) c('IH','IL','JH','JL')[which(t(rmultinom(1,1,prob=x[c('pih','pil','pjh','pjl')])) > 0)])
	#
	#tms.long <- crud[,c('Time','sampled')]
	#colnames(tms.long) <- c('Date','Comp')

	# Simulate trees with sampling once variant reaches proportion p:
	tree.10 <- getsimtree.p(tms.10,out,0.10)
	write.tree(tree.10, file= paste0('sims/varmodel/10/tree',outind,'.nwk'))

	tree.20 <- getsimtree.p(tms.20,out,0.20)
	write.tree(tree.20, file= paste0('sims/varmodel/20/tree',outind,'.nwk'))

	tree.30 <- getsimtree.p(tms.30,out,0.30)
	write.tree(tree.30, file= paste0('sims/varmodel/30/tree',outind,'.nwk'))

	tree.40 <- getsimtree.p(tms.40,out,0.40)
	write.tree(tree.40, file= paste0('sims/varmodel/40/tree',outind,'.nwk'))

	tree.50 <- getsimtree.p(tms.50,out,0.50)
	write.tree(tree.50, file= paste0('sims/varmodel/50/tree',outind,'.nwk'))

}

