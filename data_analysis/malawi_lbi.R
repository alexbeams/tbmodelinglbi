rm(list=ls())

# Calculate sequence densities using patristic distances weighted by decaying exponentials
#	(similar to THD and LBI)
 
library(ape)
library(ggplot2)
library(ggtree)
#library(Biostrings)
library(RColorBrewer)
library(gridExtra)
library(lme4)

# load the lbi function
source('lbi.R')

# read in the tree:
tree <- ape::read.tree('BactDatingLin4_relaxed_fixed_anonymized.tree')
tree <- multi2di(tree)

# load in the metadata and the lineage/drug resistance data
dat <- read.csv('tblbi_metadata_anonymized.csv')
rownames(dat) <- dat$Sequence_name

# create some convenient names for groups of variables

# Symptoms
symnms <- names(dat)[4:9]

# Socioeconomic indicators:
sesnms <- c('fridge','carmotobike','bed','radio','phone','haselectricity')

# TB characteristics:
tbnms <- c('Drug','x17tbtype')

# HIV infection
hivnms <- c('x09iptpst','x10iptdiag','x20hiv','x21art')

# Demographics
demonms <- c('x11sex')

# Smoking history
smokenms <- c('smoke')

getF <- function(tauval,var){
	# calculate LBI for each covariate
	#var = symnms[1]
	#tauval <- 1
	lbivals <- lbi(tree, tauval)
	lbivals <- lbivals[1:Ntip(tree)]
	names(lbivals) <- tree$tip.label
	lbivals <- data.frame(Sequence_name = tree$tip.label, lbi=lbivals)
	crud <- merge(dat,lbivals, by="Sequence_name")
	fstat <- summary(aov(crud$lbi~crud[,which(colnames(crud)==var)]))[[1]]$F[1]
	return(fstat)
}

#tauvals <- 10^seq(-5,2,length=30)
#Fvals <- sapply(tauvals, getF)
#plot(tauvals,Fvals)

gettauopt <- function(var){
	tauvals <- 10^seq(-5,3,length=60)
	Fvals <- sapply(tauvals, function(x) getF(x,var))
	tauopt <- tauvals[which(Fvals==max(Fvals))]
	return(tauopt)
}

# get the optimal bandwidths for the variables

symtaus <- sapply(symnms, gettauopt)
sestaus <- sapply(sesnms, gettauopt)
tbtaus <- sapply(tbnms, gettauopt)
hivtaus <- sapply(hivnms, gettauopt)
demotaus <- sapply(demonms, gettauopt)
smoketaus <- sapply(smokenms, gettauopt)

taudat <- c(symtaus,sestaus,tbtaus,hivtaus,demotaus,smoketaus)
taudat <- data.frame(var=names(taudat),tau=taudat)

plotvar <- function(var,title){
	tauval <- taudat[var,2]
	lbivals <- lbi(tree, tauval)
	lbivals <- lbivals[1:Ntip(tree)]
	names(lbivals) <- tree$tip.label
	lbivals <- data.frame(Sequence_name = tree$tip.label, lbi=lbivals)
	crud <- merge(dat,lbivals, by="Sequence_name")
	# modify the names in 'smoke' for legibility:
	crud$smoke[crud$smoke=='Currently -- regularly'] <- 'Currently\n(reg.)'
	crud$smoke[crud$smoke=="Currently – occasionally"] <- 'Currently\n(occ.)'
	crud$smoke[crud$smoke=="In the past"] <- 'In\nthe past'
	boxplot(lbi~crud[,which(colnames(crud)==var)],crud,main=title,xlab='',ylab='LBI',
		las=0)
}

# male biological sex seems to have slightly higher LBI, as we might expect from our model.
# HIV status shows nothing, also as we might expect given the model's results

# Is the effect of biological sex statistically significant? (or any others?)
getpval <- function(var){
	tauval <- taudat[var,2]
	lbivals <- lbi(tree, tauval)
	lbivals <- lbivals[1:Ntip(tree)]
	names(lbivals) <- tree$tip.label
	lbivals <- data.frame(Sequence_name = tree$tip.label, lbi=lbivals)
	crud <- merge(dat,lbivals, by="Sequence_name")
	pval <- summary(aov(crud$lbi~crud[,which(colnames(crud)==var)]))[[1]]$Pr[1]
	return(pval)
}

pvals <- sapply(taudat$var, getpval)

# Let's try using the Wilcoxon rank-sum test instead of aov to see if results change 
getkruskalpval <- function(var){
	tauval <- taudat[var,2]
	lbivals <- lbi(tree, tauval)
	lbivals <- lbivals[1:Ntip(tree)]
	names(lbivals) <- tree$tip.label
	lbivals <- data.frame(Sequence_name = tree$tip.label, lbi=lbivals)
	crud <- merge(dat,lbivals, by="Sequence_name")
	pval <- kruskal.test(crud$lbi~crud[,which(colnames(crud)==var)])$p.value
	return(pval)
}

kruskalpvals <- sapply(taudat$var, getkruskalpval)

pvals = sort(pvals)
kruskalpvals = sort(kruskalpvals)

# Let's organize these into some tables by the different covariate categories

pvals.syms <- data.frame(anova=pvals[symnms], kruskal=kruskalpvals[symnms])
pvals.ses <- data.frame(anova=pvals[sesnms], kruskal=kruskalpvals[sesnms])
pvals.tb <- data.frame(anova=pvals[tbnms], kruskal=kruskalpvals[tbnms])
pvals.hiv <- data.frame(anova=pvals[hivnms], kruskal=kruskalpvals[hivnms])
pvals.demo <- data.frame(anova=pvals[demonms], kruskal=kruskalpvals[demonms])
pvals.smoke <- data.frame(anova=pvals[smokenms], kruskal=kruskalpvals[smokenms])

aovvars.signif <- pvals[pvals < .05]
kruskvars.signif <- kruskalpvals[kruskalpvals < .05]

unionnms <- union(names(aovvars.signif),names(kruskvars.signif))

# reorder these manually
unionnms <- c('Drug','x11sex','smoke','sympbreath','phone','bed','fridge','carmotobike')

pdf(file='malawi_lbi_figure1.pdf',height=3,width=8)
par(mfrow=c(1,3))
plotvar('Drug','Drug resistance')
plotvar('x11sex','Biological sex')
plotvar('x20hiv','HIV infection') #14 in the 'Unknown' category; difference isn't statistically significant
dev.off()

pdf(file='malawi_lbi_figure_sympbreath.pdf',height=4,width=4)
plotvar('sympbreath','Reported difficulty breathing')
dev.off()

pdf(file='malawi_lbi_figure_socio.pdf',height=6,width=6)
par(mfrow=c(2,2))
plotvar('phone','Has a phone')
plotvar('bed','Has a bed')
plotvar('carmotobike','Has a car/motorcycle')
plotvar('fridge','Has a refrigerator')
#for(nm in c('phone','bed','carmotobike','fridge')){plotvar(nm)}
dev.off()

pdf(file='malawi_lbi_figure_smoke.pdf',height=4,width=6)
par(mfrow=c(1,1))
plotvar('smoke','Smoking history')
dev.off()

pdf(file='malawi_lbi_figure_combined.pdf', height=9, width=10)

# 12-column layout
# Row 1: 3 panels of width 4 each
# Row 2: 4 panels of width 3 each
# Row 3: smoking wider (e.g. cols 2-11, width 10 of 12)

layout(matrix(c(
   1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3,   # Row 1: 3 equal panels
   4, 4, 4,  5, 5, 5,  6, 6, 6,  7, 7, 7,  # Row 2: 4 equal panels
   0, 0, 0, 8, 8, 8, 8, 8, 8, 0, 0, 0      # Row 3: smoking, centered & wider
), nrow=3, byrow=TRUE))

par(mar=c(4, 4, 3, 1))

# Group 1
plotvar('Drug',        'Drug resistance')
plotvar('x11sex',      'Biological sex')
plotvar('sympbreath',  'Reported difficulty breathing')

# Group 2
plotvar('phone',       'Has a phone')
plotvar('bed',         'Has a bed')
plotvar('carmotobike', 'Has a car/motorcycle')
plotvar('fridge',      'Has a refrigerator')

# Group 3
plotvar('smoke', 'Smoking history')

dev.off()
