rm(list=ls())

require(TiPS)
require(ape)
require(ggplot2)
require(ggtree)
require(cowplot)
require(tidyr)
require(dplyr)
require(parallel)


# model with population heterogeneity and superspreaders:
 
reactions <- c(
	"S1 [(1-pfast)*zeta1*S1*c11*thetaH*IH1] -> L1",
	"S1 [(1-pfast)*zeta1*S1*c11*thetaL*IL1] -> L1",
	"S1 [(1-pfast)*zeta1*S1*c12*thetaH*IH2] -> L1",
	"S1 [(1-pfast)*zeta1*S1*c12*thetaL*IL2] -> L1",
	"S1 [pfast*pH*zeta1*S1*c11*thetaH*IH1] -> IH1",
	"S1 [pfast*pH*zeta1*S1*c11*thetaL*IL1] -> IH1",
	"S1 [pfast*pH*zeta1*S1*c12*thetaH*IH2] -> IH1",
	"S1 [pfast*pH*zeta1*S1*c12*thetaL*IL2] -> IH1",
	"S1 [pfast*(1-pH)*zeta1*S1*c11*thetaH*IH1] -> IL1",
	"S1 [pfast*(1-pH)*zeta1*S1*c11*thetaL*IL1] -> IL1",
	"S1 [pfast*(1-pH)*zeta1*S1*c12*thetaH*IH2] -> IL1",
	"S1 [pfast*(1-pH)*zeta1*S1*c12*thetaL*IL2] -> IL1",
	"L1 [pH*sigma*L1] -> IH1",
	"L1 [(1-pH)*sigma*L1] -> IL1",
	"IH1 [gamma*IH1] -> S1",
	"IL1 [gamma*IL1] -> S1",
	"L1 [rho*L1] -> S1",
	"S2 [(1-pfast)*zeta2*S2*c21*thetaH*IH1] -> L2",
	"S2 [(1-pfast)*zeta2*S2*c21*thetaL*IL1] -> L2",
	"S2 [(1-pfast)*zeta2*S2*c22*thetaH*IH2] -> L2",
	"S2 [(1-pfast)*zeta2*S2*c22*thetaL*IL2] -> L2",
	"S2 [pfast*pH*zeta2*S2*c21*thetaH*IH1] -> IH2",
	"S2 [pfast*pH*zeta2*S2*c21*thetaL*IL1] -> IH2",
	"S2 [pfast*pH*zeta2*S2*c22*thetaH*IH2] -> IH2",
	"S2 [pfast*pH*zeta2*S2*c22*thetaL*IL2] -> IH2",
	"S2 [pfast*(1-pH)*zeta2*S2*c21*thetaH*IH1] -> IL2",
	"S2 [pfast*(1-pH)*zeta2*S2*c21*thetaL*IL1] -> IL2",
	"S2 [pfast*(1-pH)*zeta2*S2*c22*thetaH*IH2] -> IL2",
	"S2 [pfast*(1-pH)*zeta2*S2*c22*thetaL*IL2] -> IL2",
	"L2 [pH*sigma*L2] -> IH2",
	"L2 [(1-pH)*sigma*L2] -> IL2",
	"IH2 [gamma*IH2] -> S2",
	"IL2 [gamma*IL2] -> S2",
	"L2 [rho*L2] -> S2"
)

# compile the model (switch out as needed)
sir_simu <- build_simulator(reactions)

