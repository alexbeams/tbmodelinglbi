# when should we start to sample?
tm1 <- 45
tms.lin4 <- runif(ntests,min=tm1,max=time[2])
tms.lin4 <- tms.lin4[order(tms.lin4)]

# use a smaller dataframe for calculating sampling probabilities (just sample rows):
traj <- out$traj
if(dim(traj)[1]<500){plottraj<-traj}else{
plottraj <- traj[c(1:1500,sample(1:dim(traj)[1],1000)),]
plottraj <- plottraj[order(plottraj$Time),]}
# transform time back to date for plotting:
plottraj$date <- as.Date(plottraj$Time * 365)

# for cross-sectional sampling similar to lineage 4:
tmstart <- 45 #min(tms.lin4)
tmend <- 50 #max(tms.lin4)


x <- plottraj[plottraj$Time < tmend & plottraj$Time > tmstart,c('Time','IH1','IL1','IH2','IL2')]

y <- colSums(x)/sum(colSums(x))

# randomly sample states according to the probabilities with which they were observed in the simulation:
s <- c('IH1','IL1','IH2','IL2')[t(rmultinom(ntests,1,prob=y[c('IH1','IL1','IH2','IL2')])) %*% c(1,2,3,4)]

# create data frames with the tms.lin4 sampling times and the sampled state:
tms <- data.frame(Date=runif(ntests,tmstart,tmend),Comp=s) 
tms <- tms[order(tms$Date),]



# Make sure that pop2frac here matches the initial conditions above
tree <- getsimtree(tms,out)

write.tree(tree, file=paste0('sims/pophetmodel/',case, thetafolder,'/tree',outind,'.nwk') )

