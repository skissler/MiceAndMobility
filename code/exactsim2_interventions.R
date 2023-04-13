# =============================================================================
# Import
# =============================================================================

library(tidyverse) 
library(purrr)
library(deSolve) 

source('~/DropboxHarvard/Projects/EpiFuncs/code/infdurabm.R')

# =============================================================================
# Define infectiousness kernels
# =============================================================================

N <- 100
beta <- 1/2
gamma <- 1/5

# SIR (stepwise) kernels: -----------------------------------------------------
pkernel_sir <- function(tau,pars){
	with(as.list(pars), {
		if(tau>0 & tau<=tstar){
			return(beta)
		} else {
			return(0)
		}
	})
}

ckernel_sir <- function(tau,pars){
	with(as.list(pars),{
		if(tau<=0){
			return(0)
		} else if(tau>0 & tau<=tstar){
			return(beta*tau)
		} else {
			return(beta*tstar)
		}
		})
}

pkernel_sir_int <- function(tau,pars){
	with(as.list(pars), {
		if(tau>0 & tau<=tstar & tau<=log(2)/gamma){
			return(beta)
		} else {
			return(0)
		}
	})
}

ckernel_sir_int <- function(tau,pars){
	with(as.list(pars),{
		if(tau<=0){
			return(0)
		} else if(tau>0 & tau<=tstar & tau<=log(2)/gamma){
			return(beta*tau)
		} else {
			return(beta*min(tstar,log(2)/gamma))
		}
		})
}

# Exponential kernels: --------------------------------------------------------
pkernel_exp <- function(tau,pars){
	with(as.list(pars),{

		if(tau>0){
			return(beta*exp(-gamma*tau))
		} else {
			return(0)
		}

		})
}

ckernel_exp <- function(tau,pars){
	with(as.list(pars),{

		if(tau>0){
			return(beta/gamma*(1-exp(-gamma*tau)))
		} else {
			return(0)
		}

		})
}


pkernel_exp_int <- function(tau,pars){
	with(as.list(pars),{

		if(tau>0 & tau<=(log(2)/gamma)){
			return(beta*exp(-gamma*tau))
		} else {
			return(0)
		}

		})
}

ckernel_exp_int <- function(tau,pars){
	with(as.list(pars),{

		if(tau<=0){
			return(0)
		} else if(tau>0 & tau<=(log(2)/gamma)){
			return(beta/gamma*(1-exp(-gamma*tau)))
		} else {
			return(0.5*beta/gamma)
		}

		})
}


# "Delta" kernels: ------------------------------------------------------------
pkernel_del <- function(tau,pars){
	with(as.list(pars),{

		if(tau<=tstar){
			return(0)
		} else if(tau>tstar & tau<(tstar+kwidth)){
			return(beta/gamma/kwidth)
		} else {
			return(0)
		}

		})
}

ckernel_del <- function(tau,pars){
	with(as.list(pars),{

		if(tau<=tstar){
			return(0)
		} else if(tau>tstar & tau<(tstar+kwidth)){
			return(beta/gamma/kwidth*(tau-tstar))
		} else {
			return(beta/gamma)
		}

		})
}

pkernel_del_int <- function(tau,pars){
	with(as.list(pars),{

		if(tau<=tstar){
			return(0)
		} else if(tau>tstar & tau<=(tstar+kwidth) & tau<=log(2)/gamma){
			return(beta/gamma/kwidth)
		} else {
			return(0)
		}

		})
}

ckernel_del_int <- function(tau,pars){
	with(as.list(pars),{

		if(tau<=tstar){
			return(0)
		} else if(tau>tstar & tau<=(tstar+kwidth) & tau<=log(2)/gamma){
			return(beta/gamma/kwidth*(tau-tstar))
		} else {
			return(beta/gamma/kwidth*(min(tstar+kwidth, log(2)/gamma)-tstar))
		}

		})
}

# =============================================================================
# Define functions for repeated simulations
# =============================================================================

parsgen_sir <- function(N,beta,gamma){
	parslist <- as.list(rexp(N,gamma)) %>% 
			map(~ c(tstar=., beta=beta, gamma=gamma))
	return(parslist)
}

parsgen_exp <- function(N,beta,gamma){
	parslist <- as.list(rep(beta, N)) %>% 
		map(~ c(beta=., gamma=gamma))
	return(parslist)
}

parsgen_del <- function(N,beta,gamma){
	parslist <- as.list(rexp(N,gamma)) %>% 
		map(~ c(tstar=., beta=beta, gamma=gamma, kwidth=0.25))
	return(parslist)
}

repsim <- function(pkernel,ckernel,parsgen,nsims,N,beta,gamma,invckernel=NA,initinf=1,initstep=1,maxits=1000,quiet=TRUE){

	# Run simulations
	simoutlist <- list()
	print(paste0("Loop started: ",Sys.time()))
		for(simnum in 1:nsims){
		parslist <- parsgen(N,beta,gamma)
		simout <- infdurabm(pkernel, ckernel, parslist, invckernel=invckernel, initinf=initinf, initstep=initstep, maxits=maxits, quiet=quiet)
		simoutlist[[simnum]] <- simout
		if(simnum%%10==0){
			print(paste0("Simulation number ",simnum," completed"))
		}
	}
	print(paste0("Loop ended: ",Sys.time()))

	# Tidy the output
	simoutdf <- simoutlist %>% 
		imap(~ mutate(.x, sim=.y)) %>% 
		bind_rows() %>% 
		mutate(counter=case_when(tinf<Inf~1, TRUE~0)) %>% 
		group_by(sim) %>% 
		arrange(sim, tinf) %>% 
		group_by(sim) %>% 
		mutate(cuminf=cumsum(counter)) %>% 
		select(-counter) %>% 
		split(.$sim) %>% 
		map(~ bind_rows(., tibble(id=-1, tinf=Inf, whoinf=-1, sim=min(.$sim), cuminf=last(.$cuminf)))) %>% 
		bind_rows()

	return(simoutdf)
}

# =============================================================================
# Generate SIR comparator
# =============================================================================

sir <- function(t,state,parameters){
	with(as.list(c(state,parameters)), {

		dS <- -beta*I*S
		dI <- beta*I*S - gamma*I
		dcumI <- beta*I*S
		dR <- gamma*I

		list(c(dS, dI, dcumI, dR))

		})
}

parameters <- c(beta=beta, gamma=gamma) 
state <- c(S=1-1/N, I=1/N, cumI=1/N, R=0)
times <- seq(from=0, to=100, by=0.01)

sirout <- ode(y = state, times = times, func = sir, parms = parameters) %>% 
	data.frame() %>% 
	as_tibble()

# =============================================================================
# Do repeated simulations
# =============================================================================

simoutdf_sir <- repsim(pkernel_sir, ckernel_sir, parsgen_sir, nsims=100, N=N, beta=beta, gamma=gamma)

simoutdf_sir_int <- repsim(pkernel_sir_int, ckernel_sir_int, parsgen_sir, nsims=100, N=N, beta=beta, gamma=gamma)

simoutdf_exp <- repsim(pkernel_exp, ckernel_exp, parsgen_exp, nsims=100, N=N, beta=beta, gamma=gamma)

simoutdf_exp_int <- repsim(pkernel_exp_int, ckernel_exp_int, parsgen_exp, nsims=100, N=N, beta=beta, gamma=gamma)

simoutdf_del <- repsim(pkernel_del, ckernel_del, parsgen_del, nsims=100, N=N, beta=beta, gamma=gamma)

simoutdf_del_int <- repsim(pkernel_del_int, ckernel_del_int, parsgen_del, nsims=100, N=N, beta=beta, gamma=gamma)

# save(simoutdf_sir, file="output/simoutdf_sir.RData")
# save(simoutdf_exp, file="output/simoutdf_exp.RData")
# save(simoutdf_del, file="output/simoutdf_del.RData")

# =============================================================================
# Basic plots 
# =============================================================================

# Simulated trajectories against the theoretical SIR: -------------------------

plot_simcurves <- function(simoutdf, sirout=NA){
	out <- simoutdf %>% 
		ggplot() + 
			geom_step(aes(x=tinf, y=cuminf, group=sim), col="grey", alpha=0.2) + 
			theme_classic() + 
			labs(x="Time", y="Cumulative infections") + 
			theme(text=element_text(size=9))

	if(is.data.frame(sirout)){
		out <- out + geom_line(data=sirout, aes(x=time, y=cumI*N), col="blue", size=1)
	}

	return(out)
}

figsimcurves_sir <- plot_simcurves(simoutdf_sir, sirout)
figsimcurves_sir_int <- plot_simcurves(simoutdf_sir_int, sirout)
figsimcurves_exp <- plot_simcurves(simoutdf_exp, sirout)
figsimcurves_exp_int <- plot_simcurves(simoutdf_exp_int, sirout)
figsimcurves_del <- plot_simcurves(simoutdf_del, sirout)
figsimcurves_del_int <- plot_simcurves(simoutdf_del_int, sirout)

# ggsave(figsimcurves_sir, file="figures/simcurves_sir.pdf", width=3.5, height=3.5/1.6)
# ggsave(figsimcurves_exp, file="figures/simcurves_exp.pdf", width=3.5, height=3.5/1.6)
# ggsave(figsimcurves_del, file="figures/simcurves_del.pdf", width=3.5, height=3.5/1.6)

# a_i: -----------------------------------------------------------------------

parslist <- parsgen_sir(50,beta,gamma)
fig_aisir <- plotai(pkernel_sir, parslist, tmin=0.01, tmax=30)
fig_aisir_int <- plotai(pkernel_sir_int, parslist, tmin=0.01, tmax=30)

parslist <- parsgen_exp(50,beta,gamma)
fig_aiexp <- plotai(pkernel_exp, parslist, tmin=0.01, tmax=30)
fig_aiexp_int <- plotai(pkernel_exp_int, parslist, tmin=0.01, tmax=30)

parslist <- parsgen_del(50,beta,gamma)
fig_aidel <- plotai(pkernel_del, parslist, tmin=0.01, tmax=30)
fig_aidel_int <- plotai(pkernel_del_int, parslist, tmin=0.01, tmax=30)

# ggsave(fig_aisir, file="figures/aisir.pdf", width=3.5, height=3.5/1.6)
# ggsave(fig_aiexp, file="figures/aiexp.pdf", width=3.5, height=3.5/1.6)
# ggsave(fig_aidel, file="figures/aidel.pdf", width=3.5, height=3.5/1.6)

# A: ------------------------------------------------------------------------

cf <- function(x){beta*exp(-gamma*x)}

parslist <- parsgen_sir(500,beta,gamma)
fig_Asir <- plotA(pkernel_sir, parslist, tmax=30, compfun=cf)
fig_Asir_int <- plotA(pkernel_sir_int, parslist, tmax=30, compfun=cf)

parslist <- parsgen_exp(500,beta,gamma)
fig_Aexp <- plotA(pkernel_exp, parslist, tmax=30, compfun=cf)
fig_Aexp_int <- plotA(pkernel_exp_int, parslist, tmax=30, compfun=cf)

parslist <- parsgen_del(500,beta,gamma)
fig_Adel <- plotA(pkernel_del, parslist, tmax=30, compfun=cf)
fig_Adel_int <- plotA(pkernel_del_int, parslist, tmax=30, compfun=cf)

# ggsave(fig_Asir, file="figures/Asir.pdf", width=3.5, height=3.5/1.6)
# ggsave(fig_Aexp, file="figures/Aexp.pdf", width=3.5, height=3.5/1.6)
# ggsave(fig_Adel, file="figures/Adel.pdf", width=3.5, height=3.5/1.6)

# Final size: -----------------------------------------------------------------

plotFSdists <- function(simoutdflist, modelnames=1:1000, binwidth=5, plotdensity=FALSE){

	cuminfdf <- simoutdflist[[1]] %>% 
		group_by(sim) %>% 
		summarise(cuminf=max(cuminf)) %>% 
		mutate(model=modelnames[1])

	if(length(simoutdflist)>1){

		for(indexA in 2:length(simoutdflist)){
			cuminfdf <- bind_rows(cuminfdf, (simoutdflist[[indexA]] %>% 
				group_by(sim) %>% 
				summarise(cuminf=max(cuminf)) %>% 
				mutate(model=modelnames[indexA])))
		}

	}


	out <- cuminfdf %>% 
		mutate(model=factor(model, levels=modelnames)) %>% 
		ggplot(aes(x=cuminf)) + 
			geom_histogram(aes(y=..density..), binwidth=binwidth, fill="white", col="darkgrey") + 
			facet_wrap(~model) + 
			theme_classic() + 
			labs(x="Cumulative infections", y="Proportion of simulations") + 
			theme(text=element_text(size=9))

	if(plotdensity){
		out <- out + geom_density(adjust=0.5)
	}

	return(out)

}

fig_FSdists <- plotFSdists(list(simoutdf_sir, simoutdf_exp, simoutdf_del, simoutdf_sir_int, simoutdf_exp_int, simoutdf_del_int), modelnames=c("SIR","Exponential","Delta function","SIR int.","Exponential int.","Delta function int."), binwidth=2)
# ggsave(fig_FSdists, file="figures/FSdists.pdf", width=3.5, height=3.5/1.6)

# Time to 100 infections: -----------------------------------------------------

plottimetox <- function(simoutdflist, threshold, modelnames=1:1000, binwidth=5, plotdensity=FALSE){

	ttxdf <- simoutdflist[[1]] %>% 
		group_by(sim) %>% 
		filter(cuminf >= threshold) %>% 
		arrange(tinf) %>% 
		slice(1) %>% 	
		mutate(model=modelnames[1]) 

	if(length(simoutdflist)>1){

		for(indexA in 2:length(simoutdflist)){
			ttxdf <- bind_rows(ttxdf, (simoutdflist[[indexA]] %>% 
				group_by(sim) %>% 
				filter(cuminf >= threshold) %>% 
				arrange(tinf) %>% 
				slice(1) %>% 	
				mutate(model=modelnames[indexA])))
		}

	}

	out <- ttxdf %>% 
		mutate(model=factor(model, levels=modelnames)) %>% 
		ggplot(aes(x=tinf)) + 
			geom_histogram(binwidth=binwidth, fill="white", col="darkgrey") + 
			facet_wrap(~model) + 
			theme_classic() + 
			labs(x=paste0("Time to ",threshold," infections"), y="Number of simulations") + 
			theme(text=element_text(size=9))


	return(out)

}

fig_timeto50 <- plottimetox(list(simoutdf_sir, simoutdf_exp, simoutdf_del, simoutdf_sir_int, simoutdf_exp_int, simoutdf_del_int), threshold=50, modelnames=c("SIR","Exponential","Delta function","SIR int.","Exponential int.","Delta function int."), binwidth=2)
# ggsave(fig_timeto50, file="figures/timeto50.pdf", width=3.5, height=3.5/1.6)


# =============================================================================
# Raw Gillespie comparison
# =============================================================================

# # Compare with a gillespie algorithm: 
# source('~/DropboxHarvard/Projects/EpiFuncs/code/gillespie.R')

# states <- c(S=199, I=1, R=0)
# parms <- c(beta=beta, gamma=gamma, N=sum(states))
# rates <- c(
# 	"S -> I"="beta*S*I/N", 
# 	"I -> R"="gamma*I")

# gilloutlist <- list()
# for(simnum in 1:200){
# 	gillout <- gillespie(states,parms,rates,maxits=10000)
# 	gilloutlist[[simnum]] <- gillout 
# 	if(simnum%%10==0){
# 		print(paste0("Simulation number ",simnum," completed"))
# 	}	
# }

# gilloutdf <- gilloutlist %>% 
# 	imap(~ mutate(.x, sim=.y)) %>% 
# 	map(~ bind_rows(., tibble(t=Inf, sim=min(.$sim), S=last(.$S), I=last(.$I), R=last(.$R)))) %>% 
# 	bind_rows()

# figgillcurves <- gilloutdf %>% 
# 	ggplot() + 
# 		geom_step(aes(x=t, y=R, group=sim), col="grey", alpha=0.2) + 
# 		geom_line(data=sirout, aes(x=time, y=R*N), col="blue", size=1) + 
# 		theme_classic() 


# Key points for the talk: 
# You need much larger sample sizes in the delta scenario to overcome the stochasticity than you do in the sir or exponential scenarios 
# Uncertainty in time of infection onsets can spread a delta into something much broader. 
# Stochastic dynamics depend massively on individual infectiousness profiles 
# Interventions also depend on these profiles 
# 