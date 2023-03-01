# =============================================================================
# Initial definitions
# =============================================================================

library(tidyverse) 
library(purrr)

N <- 1000
beta <- 1/3
gamma <- 1/5

tauistar <- rexp(N, gamma)

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
