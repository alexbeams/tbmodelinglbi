rm(list=ls())

require(TiPS)
require(ape)
require(ggplot2)
require(ggtree)
require(cowplot)
require(tidyr)
require(dplyr)


# model with a variant and superspreaders:
 
reactions <- c(
	"S [(1-pfast)*zeta*S*c*thetaH*IH] -> L",
	"S [(1-pfast)*zeta*S*c*thetaL*IL] -> L",
	"S [pfast*pH*zeta*S*c*thetaH*IH] -> IH",
	"S [pfast*pH*zeta*S*c*thetaL*IL] -> IH",
	"S [pfast*(1-pH)*zeta*S*c*thetaH*IH] -> IL",
	"S [pfast*(1-pH)*zeta*S*c*thetaL*IL] -> IL",
	"L [pH*(1-pmu)*sigma*L] -> IH",
	"L [(1-pH)*(1-pmu)*sigma*L] -> IL",
	"IH [gamma*IH] -> S",
	"IL [gamma*IL] -> S",
	"L [rho*L] -> S",
	"L [pH*pmu*sigma*L] -> JH",
	"L [(1-pH)*pmu*sigma*L] -> JL",
	"S [(1-pfast)*zeta*S*c*(thetaH+deltathetaH)*JH] -> M",
	"S [(1-pfast)*zeta*S*c*(thetaL+deltathetaL)*JL] -> M",
	"S [pfast*pH*zeta*S*c*(thetaH+deltathetaH)*JH] -> JH",
	"S [pfast*pH*zeta*S*c*(thetaL+deltathetaL)*JL] -> JH",
	"S [pfast*(1-pH)*zeta*S*c*(thetaH+deltathetaH)*JH] -> JL",
	"S [pfast*(1-pH)*zeta*S*c*(thetaL+deltathetaL)*JL] -> JL",
	"M [pH*sigma*M] -> JH",
	"M [(1-pH)*sigma*M] -> JL",
	"M [rho*M] -> S",
	"JH [gamma*JH] -> S",
	"JL [gamma*JL] -> S"
)

# compile the model (switch out as needed)
sir_simu <- build_simulator(reactions)

