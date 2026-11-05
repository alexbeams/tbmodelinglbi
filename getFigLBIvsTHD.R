rm(list=ls())


## load in the lineage 4 tree and examine the distribution of LBI, THD, and No. of close relatives
##	within a threshold patristic distance

## Set thresholds of 12 years for the different clustering methods, and choose bandwidths of LBI and
## THD that correlate most strongly with these

source('lbi.R')

require(ape)
require(ggtree)
require(scales)

# Load in the lineage 4 tree:
lin4phy <- read.tree('data_analysis/BactDatingLin4_relaxed_fixed_anonymized.tree')

# calculate lbi:
lbivals <- lbi(lin4phy, tau = 1) 
names(lbivals) <- lin4phy$tip.label
lbivals <- lbivals - min(lbivals)
lbivals <- lbivals/max(lbivals)

tiplbis <- lbivals[1:513]
nodelbis <- lbivals[514:1025]
tipcols <- alpha('purple', alpha=tiplbis)
nodecols <- alpha('purple', alpha=nodelbis)


# make some plots showing that THD and LBI are highly correlated:

thdvslbi <- function(tau1,tau2){
	x <- cophenetic(lin4phy)
	tau <- mean(x)*tau1
	dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
	dat <- as.data.frame(dat)
	colnames(dat) <- 'thd'
	dat$label = rownames(dat)
	dat <- dat[,c(2,1)]

	dat$var <- apply(x,1,var)
	dat$min <- apply(x,1,function(x) min(x[x>0]) )

	dat$lbi <- lbi(lin4phy,tau=tau2)[1:length(lin4phy$tip.label)]
	
	w <- cor(dat$lbi,dat$thd)

	plot(lbi~thd,dat,xlab='THD',ylab='LBI',cex.main=2,cex.lab=2,
		cex.axis=2, main=bquote('Correlation'==.(round(w,2)) ) )
	return(w)
}

#this has a high correlation. Can we do any better?
thdvslbi(4/mean(cophenetic(lin4phy)),5)

obj <- function(logtau2) thdvslbi(.035,exp(logtau2)) 

opt <- optim(log(5), obj, method='Brent', control=list(fnscale=-1),
	lower=-5,upper=10)

tau2opt <- exp(opt$par)

obj2 <- function(logtau1) thdvslbi(exp(logtau1), tau2opt)

opt2 <- optim(log(.035), obj2, method='Brent', control=list(fnscale=-1),
	lower=-5,upper=10)

tau1opt <- exp(opt2$par)

snpvslbi <- function(tau1,tau2){
	x <- cophenetic(lin4phy)
	tau <- mean(x)*tau1
	dat <- apply(x, 1, function(z) sum(exp(-z/tau)))
	dat <- as.data.frame(dat)
	colnames(dat) <- 'thd'
	dat$label = rownames(dat)
	dat <- dat[,c(2,1)]

	# what is the variance in patristic distances to all other tips?
	dat$var <- apply(x,1,var)

	# what is the distance to the closest relative?
	dat$min <- apply(x,1,function(x) min(x[x>0]) )

	# how many other tips within tau1 units of patristic distance?
	dat$clust <- apply(x,1,function(x) sum(x < tau1)) - 1

	# calculate LBI for each tip:
	dat$lbi <- lbi(lin4phy,tau=tau2)[1:length(lin4phy$tip.label)]
	
	w <- cor(dat$lbi,dat$clust)

	plot(lbi~clust,dat,xlab='No. of close relatives',ylab='LBI',cex.main=2,cex.lab=2,
		cex.axis=2, main=bquote('Correlation'==.(round(w,2))) )
	return(w)
}


tree <- lin4phy

taulbi <- 3
tauthd <- 5
taurels <- 12
tauclust <- 12

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

#count the number of close relatives within a threshold distance:
dat$nrelatives <- apply(x,1,function(x) sum(x <= taurels)) - 1

# calculate clusters using single linkage:
hc <- hclust(as.dist(x), method='single')
clusts <- cutree(hc, h=tauclust)
dat$Cluster <- clusts

# calculate cluster size:
dat[,'Cluster Size (single)'] <- table(dat$Cluster)[dat$Cluster]

# calculate clusters using complete linkage (identical w/ No. of close relatives?)
hc <- hclust(as.dist(x), method='complete')
clusts <- cutree(hc, h=tauclust)
dat[,'Cluster Size (complete)'] <- table(clusts)[clusts]

# create a ggtree plot:
plin4 <- ggtree(tree, layout='rectangular')

# rename the columns:
colnames(dat)[4] <- "No. of close relatives"


pdf(file='figures/FigS3.pdf',width=6,height=15)
par(mfrow=c(5,2),mar=c(5,5,4,2))


# 1. LBI ~ THD

plot(LBI~THD,dat,xlab='THD',ylab='LBI',cex.main=2,cex.lab=2,
		cex.axis=2, main=bquote('Correlation'==.(round(cor(dat$LBI,dat$THD),2))) )

# 2. LBI ~ Cluster size (single)
x <- dat[,'Cluster Size (single)']
x <- c(unlist(x))

plot(x,dat[,'LBI'],
	xlab='Cluster Size (single)',
	ylab='LBI',cex.main=2,cex.lab=2,
	cex.axis=2,
	main=bquote('Correlation'==.(round(cor(dat$LBI,x),2))) )

# 3. LBI ~ Cluste size (complete)
x <- dat[,'Cluster Size (complete)']
x <- c(unlist(x))
x <- jitter(x)

plot(x,dat[,'LBI'],
	xlab='Cluster Size (complete)',
	ylab='LBI',cex.main=2,cex.lab=2,
	cex.axis=2,
	main=bquote('Correlation'==.(round(cor(dat$LBI,x),2))) )

# 4. LBI ~ No. of close relatives

plot(dat[,'No. of close relatives'],dat[,'LBI'],xlab='No. of close relatives',ylab='LBI',cex.main=2,cex.lab=2,
		cex.axis=2, main=bquote('Correlation'==.(round(cor(dat$LBI,dat[,'No. of close relatives']),2))) )

# 5. THD ~ Cluster Size (single)

x <- dat[,'Cluster Size (single)']
x <- c(unlist(x))
x <- jitter(x)

plot(x,dat[,'THD'],
	xlab='Cluster Size (single)',
	ylab='THD',cex.main=2,cex.lab=2,
	cex.axis=2,
	main=bquote('Correlation'==.(round(cor(dat$THD,x),2))) )

# 6. THD ~ Cluster Size (complete)

x <- dat[,'Cluster Size (complete)']
x <- c(unlist(x))
x <- jitter(x)

plot(x,dat[,'THD'],
	xlab='Cluster Size (complete)',
	ylab='THD',cex.main=2,cex.lab=2,
	cex.axis=2,
	main=bquote('Correlation'==.(round(cor(dat$THD,x),2))) )

# 7. THD ~ No. of close relatives

plot(dat[,'No. of close relatives'],dat[,'THD'],xlab='No. of close relatives',ylab='THD',cex.main=2,cex.lab=2,
		cex.axis=2, main=bquote('Correlation'==.(round(cor(dat$THD,dat[,'No. of close relatives']),2))) )

# 8 . Cluster size (single) ~ Cluster Size (complete)

x <- dat[,'Cluster Size (single)']
x <- c(unlist(x))
x <- jitter(x)

y <- dat[,'Cluster Size (complete)']
y <- c(unlist(y))
y <- jitter(y)

plot(y,x,
	ylab='Cluster Size (single)',
	xlab='Cluster Size (complete)',cex.main=2,cex.lab=2,
	cex.axis=2,
	main=bquote('Correlation'==.(round(cor(y,x),2))) )

# 9. Cluster size (single) ~ No. of close relatives

x <- dat[,'Cluster Size (complete)']
x <- c(unlist(x))
x <- jitter(x)

plot(dat[,'No. of close relatives'],x,
	ylab='Cluster Size (complete)',
	xlab='No. of close relatives',cex.main=2,cex.lab=2,
	cex.axis=2,
	main=bquote('Correlation'==.(round(cor(dat$THD,x),2))) )

# 10. Cluster size (single) ~ No. of close relatives

x <- dat[,'Cluster Size (single)']
x <- c(unlist(x))
x <- jitter(x)

plot(dat[,'No. of close relatives'],x,xlab='No. of close relatives',ylab='Cluster Size (single)',cex.main=2,cex.lab=2,
		cex.axis=2, main=bquote('Correlation'==.(round(cor(x,dat[,'No. of close relatives']),2))) )

dev.off()



