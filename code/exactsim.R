# =============================================================================
# Initial definitions
# =============================================================================

# Import packages:
library(tidyverse) 
library(purrr)
library(deSolve)

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


sir <- function(t,state,parameters){
	with(as.list(c(state,parameters)), {

		dS <- -beta*I*S
		dI <- beta*I*S - gamma*I
		dR <- gamma*I

		list(c(dS, dI, dR))

		})
}

parameters <- c(beta=beta, gamma=gamma) 
state <- c(S=1-1/N, I=1/N, R=0)
times <- seq(from=0, to=100, by=0.01)

sirout <- ode(y = state, times = times, func = sir, parms = parameters) %>% 
	data.frame() %>% 
	as_tibble()

tempmean <- temp %>% 
	group_by(iteration) %>% 
	mutate(peak = max(inf)) %>% 
	filter(peak>20) %>% 
	group_by(t) %>% 
	summarise(inf=mean(inf))


fig_temp + 
	geom_line(data=sirout, aes(x=time, y=I*N), col="lightgrey", group=1, size=1, alpha=1) + 
	geom_line(data=tempmean, aes(x=t, y=inf), col="lightgrey", alpha=1, group=1, size=1, lty="dashed")



