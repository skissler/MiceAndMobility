# =============================================================================
# Initial definitions
# =============================================================================

library(tidyverse) 
library(purrr)

N <- 1000
beta <- 1/3
gamma <- 1/5

tauistar <- rexp(N, gamma)
# vi <- beta*tauistar

# ai <- function(tau, beta, tauistar){
# 	if(tau<0){stop("tau must be positive")}
# 	if(tau<=tauistar){
# 		out <- beta
# 	} else {
# 		out <- 0
# 	}
# 	return(out)
# }

# Fi <- function(tau, beta, tauistar){
# 	if(tau<0){stop("tau must be positive")}
# 	if(tau<=tauistar){
# 		out <- beta*tau
# 	} else {
# 		out <- beta*tauistar
# 	}
# 	return(out)
# }


# Fhati <- function(tau, tauistar){
# 	if(tau<0){stop("tau must be positive")}
# 	if(tau<=tauistar){
# 		out <- tau/tauistar
# 	} else {
# 		out <- 1
# 	}
# 	return(out)
# }

# Fhatinvi <- function(x, tauistar){

# 	out <- x*tauistar
# 	return(out)

# }

# =============================================================================
# Simulate
# =============================================================================

# source('code/exactsim.R'); tinfvec[tinfvec<Inf]

tmax <- 100

tinfvec <- rep(Inf,N)

tinfvec[1] <- 0

t <- 0
infcounter <- 2 

while(t < tmax){
	infremaining <- tauistar - (t - tinfvec)
	delta <- min(infremaining[infremaining>0.000001]) # if Inf, stop sim
	if(delta==Inf){
		break()
	}

	# How many people are currently contributing to infectiousness? 
	# That's people who have been infected and not yet recovered: 
	ninf <- sum((tinfvec < Inf) & (tauistar > t-tinfvec))

	# Calculate the total force over the time period from these people: 
	Ftot <- ninf*beta*delta
	pinf <- 1-exp(-Ftot/N)

	newinf <- rbinom(1,N-sum(tinfvec<Inf),pinf)

	# Draw times for the new infections: (add clause where if there are no infections we skip this and set t <- t+delta) 
	if(newinf > 0){
		t <- t + min(runif(newinf))*delta
		tinfvec[infcounter] <- t
		infcounter <- infcounter + 1	
	} else {
		t <- t+delta
	}
}


tinfdf <- tibble(tinf=tinfvec, taui=tauistar) %>% 
	mutate(trec=tinf+taui)

plottimes <- seq(from=0,to=tmax,by=1)
plotinfs <- lapply(plottimes, function(x){
	return(nrow(tinfdf %>% filter(tinf<=x & trec>x)))
	}) %>% unlist()

fig_icurve <- tibble(t=plottimes,inf=plotinfs) %>% 
	ggplot(aes(x=t, y=inf)) +
		geom_line() + 
		theme_classic()


# =============================================================================
# Test
# =============================================================================


# # Initialize time 
# t <- 0

# while(t<tmax){

# 	# Calculate total force 
# 	tau <- t-tinfvec
# 	tau[tau>tauistar] <- tauistar[tau>tauistar]
# 	Fi <- ((beta*(tau+tstep)) - (beta*tau))[ivec==1]
# 	Ftot <- sum(Fi)

# 	# Draw new infections 
# 	pinf <- 1-exp(-Ftot/N)
# 	ninf <- rbinom(1,N-sum(ivec),pinf)

# 	# Determine the first infection
# 	if(ninf>0){


# 	} else {
# 		t <- t+tstep
# 	}


# }


# # =============================================================================

# tinfvec <- rep(Inf,N)
# tinfvec[sample(N, 50)] <- runif(50)*10
# ivec <- rep(0,N)
# ivec[tinfvec<Inf]<-1

# t <- 11
# ftot <- function(t){
# 	out <- sum((beta*pmin(tauistar, t-tinfvec))[ivec==1])/N
# 	return(out)
# }

# # The full cumulative infectiousness profile:
# tvec <- seq(from=10, to=100, by=0.1) 
# yvec <- unlist(lapply(tvec, ftot))
# tibble(t=tvec,y=yvec) %>% ggplot(aes(x=t,y=y)) + geom_line() + theme_classic()

# # The truncated infectiousness CDF: 
# tvec <- seq(from=10, to=40, by=0.1) 
# yvec <- (unlist(lapply(tvec, ftot))-ftot(10))/(ftot(40)-ftot(10))
# tibble(t=tvec,y=yvec) %>% ggplot(aes(x=t,y=y)) + geom_line() + theme_classic()
