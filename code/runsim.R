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


episim_exp(k=k, phi=phi, mu=mu, sigma=sigma, domainwidth=domainwidth)


# makevideo(eventlist, tstep=0.1, file="figures/episim.gif")

