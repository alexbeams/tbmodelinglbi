require(ape)
require(xtable)

treenms10 <- list.files('10')
treenms20 <- list.files('20')
treenms30 <- list.files('30')
treenms40 <- list.files('40')
treenms50 <- list.files('50')

treelist10 <- lapply(treenms10, function(x) read.tree(paste0('10/',x)))
treelist20 <- lapply(treenms20, function(x) read.tree(paste0('20/',x)))
treelist30 <- lapply(treenms30, function(x) read.tree(paste0('30/',x)))
treelist40 <- lapply(treenms40, function(x) read.tree(paste0('40/',x)))
treelist50 <- lapply(treenms50, function(x) read.tree(paste0('50/',x)))


colless10 <- sapply(treelist10, collessI)
colless20 <- sapply(treelist20, collessI)
colless30 <- sapply(treelist30, collessI)
colless40 <- sapply(treelist40, collessI)
colless50 <- sapply(treelist50, collessI)


sackin10 <- sapply(treelist10, sackinI)
sackin20 <- sapply(treelist20, sackinI)
sackin30 <- sapply(treelist30, sackinI)
sackin40 <- sapply(treelist40, sackinI)
sackin50 <- sapply(treelist50, sackinI)

imbalancedat <- data.frame(VariantFrequency = c(10,20,30,40,50), 
	Colless = c(mean(colless10),mean(colless20),mean(colless30),mean(colless40),mean(colless50)),
	Sackin = c(mean(sackin10),mean(sackin20),mean(sackin30),mean(sackin40),mean(sackin50))

)

xtable(imbalancedat)