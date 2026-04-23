# do N-yr sampling with ntests and save the tree

# set sampling times:
tms <- runif(ntests,min=time[2]-Nyr,max=time[2])
tms <- tms[order(tms)]

tree <- getsimtree(tms,out, pH=0.1)
write.tree(tree, file=treepath)


