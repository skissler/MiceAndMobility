# Import packages
library(tidyverse) 
library(animation) 
library(gganimate)
source('code/utils.R')

# Define model parameters
domainwidth <- 10  # Width of the domain of movement
sigma <- 1         # Standard deviation of movement distance 
mu <- 1            # Rate of movements 
k <- 1             # Force of infection when distance = 0
phi <- 1           # Rate of decay of the distance kernel 


tinf <- lapply(1:1000, function(x){
	episim_exp(k=k, phi=phi, mu=mu, sigma=sigma, domainwidth=domainwidth) %>% 
	filter(IsInf==1) %>% 
	pull(t)}) %>% unlist()

# makevideo(eventlist, tstep=0.1, file="figures/episim.gif")

ggplot(data=tibble(t=tinf), aes(x=t)) + 
	geom_histogram(aes(y=..density..),binwidth=10,fill="white",col="darkgray") +
	labs(x="Time of infection", y="Density") + 
	theme_classic()