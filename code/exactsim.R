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

tstep <- 1
tmax <- 100

ivec <- rep(0,N)
tinfvec <- rep(Inf,N)

ivec[sample(N,1)] <- 1
tinfvec[ivec==1] <- 0

# calc_total_force(t, ivec, tinfvec) 

# draw_new_infs(ivec, force) 

# draw_inf_timing(tinf, ivec, tinfvec) 


# =============================================================================
# Test
# =============================================================================

# Initialize time 
t <- 0

while(t<tmax){

	# Calculate total force 
	tau <- t-tinfvec
	tau[tau>tauistar] <- tauistar[tau>tauistar]
	Fi <- ((beta*(tau+tstep)) - (beta*tau))[ivec==1]
	Ftot <- sum(Fi)

	# Draw new infections 
	pinf <- 1-exp(-Ftot/N)
	ninf <- rbinom(1,N-sum(ivec),pinf)

	# Determine the first infection
	if(ninf>0){


	} else {
		t <- t+tstep
	}


}


# =============================================================================

tinfvec <- rep(Inf,N)
tinfvec[sample(N, 50)] <- runif(50)*10
ivec <- rep(0,N)
ivec[tinfvec<Inf]<-1

t <- 11
ftot <- function(t){
	out <- sum((beta*pmin(tauistar, t-tinfvec))[ivec==1])/N
	return(out)
}

# The full cumulative infectiousness profile:
tvec <- seq(from=10, to=100, by=0.1) 
yvec <- unlist(lapply(tvec, ftot))
tibble(t=tvec,y=yvec) %>% ggplot(aes(x=t,y=y)) + geom_line() + theme_classic()

# The truncated infectiousness CDF: 
tvec <- seq(from=10, to=40, by=0.1) 
yvec <- (unlist(lapply(tvec, ftot))-ftot(10))/(ftot(40)-ftot(10))
tibble(t=tvec,y=yvec) %>% ggplot(aes(x=t,y=y)) + geom_line() + theme_classic()
