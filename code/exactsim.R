# =============================================================================
# Initial definitions
# =============================================================================

# Import packages:
library(tidyverse) 
library(purrr)


# Set epidemiological parameters:
N <- 1000
beta <- 1/3
gamma <- 1/5



# =============================================================================
# Simulate
# =============================================================================

tinflist <- list()
for(its in 1:500){

# Draw infectious periods for each person: 
tauistar <- rexp(N, gamma)

# Set a max number of days for simulation
tmax <- 100

# Initialize vector to track infection times 
tinfvec <- rep(Inf,N)
# Initialize the first infection
tinfvec[1] <- 0

# Start at time t = 0, and the next infection will go into position 2
t <- 0
infcounter <- 2 

while(t < tmax){
	# How much time is remaining in each person's infectious period? 
	infremaining <- tauistar - (t - tinfvec)
	# Find the minimum remaining infectiousness (adjust for machine precision)
	delta <- min(infremaining[infremaining>0.000001]) 
	# If delta = Inf, no one is infectious anymore 
	if(delta==Inf){ 
		break()
	}

	# Count how many people are currently infectious:
	ninf <- sum((tinfvec < Inf) & (tauistar > t-tinfvec))

	# Calculate the total force over the time period from these people: 
	Ftot <- ninf*beta*delta

	# Calculate the probability of inection for each susceptible person:
	pinf <- 1-exp(-Ftot/N)

	# Draw the number of new infections
	newinf <- rbinom(1,N-sum(tinfvec<Inf),pinf)

	# If there were new infections, draw their times and take the min: 
	if(newinf > 0){
		t <- t + min(runif(newinf))*delta
		tinfvec[infcounter] <- t
		infcounter <- infcounter + 1	
	} else {
		t <- t+delta
	}
}
tinflist[[its]] <- tibble(tinf=tinfvec,tauistar=tauistar)
}

# Collect infection times into a table: 
tinfdf <- tibble(tinf=tinfvec, taui=tauistar) %>% 
	mutate(trec=tinf+taui)
plottimes <- seq(from=0,to=tmax,by=1)
plotinfs <- lapply(plottimes, function(x){
	return(nrow(tinfdf %>% filter(tinf<=x & trec>x)))
	}) %>% unlist()

# Plot number of people infectious over time: 
fig_icurve <- tibble(t=plottimes,inf=plotinfs) %>% 
	ggplot(aes(x=t, y=inf)) +
		geom_line() + 
		theme_classic()


geti <- function(tinfdf, plottimes){

	out <- lapply(plottimes, function(x){
		return(nrow(tinfdf %>% filter(tinf<=x & trec>x)))
		}) %>% 
		unlist()
	return(out)

}

# Plot all the curves in the list: 
plottimes <- seq(from=0,to=tmax,by=1)
temp <- tinflist %>% 
	map(~ mutate(., trec=tinf+tauistar)) %>% 
	map(~ geti(., plottimes)) %>% 
	map(~ tibble(t=plottimes, inf=.)) %>% 
	bind_rows(.id="iteration")

fig_temp <- temp %>% 
	ggplot(aes(x=t, y=inf, group=factor(iteration))) + 
		geom_line(alpha=0.2) + 
		theme_classic() + 
		theme(legend.position="none")



