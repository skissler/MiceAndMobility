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

# Set initial conditions 
t <- 0
Sposx <- runif(1)*domainwidth
Sposy <- runif(1)*domainwidth
Iposx <- runif(1)*domainwidth
Iposy <- runif(1)*domainwidth
IsInf <- 0
eventlist <- tibble(t=t, 
	Sposx=Sposx, Sposy=Sposy, 
	Iposx=Iposx, Iposy=Iposy, 
	IsInf=IsInf)

# Run simulation
while(IsInf<1){

	# Calculate key parameters 
	d <- sqrt((Sposx-Iposx)^2 + (Sposy-Iposy)^2)
	lambda <- k*exp(-phi*d)
	cumrate <- 2*mu + lambda

	# Simulate the time of the next event
	t <- t+rexp(1, rate=cumrate)

	# Draw the event 
	eventdraw <- runif(1)
	if(eventdraw < mu/cumrate){ # S moves
		Sposx <- update_pos(Sposx, sigma, domainwidth)
		Sposy <- update_pos(Sposy, sigma, domainwidth)
	} else if(eventdraw < 2*mu/cumrate){ # I moves
		Iposx <- update_pos(Iposx, sigma, domainwidth)
		Iposy <- update_pos(Iposy, sigma, domainwidth)
	} else { # S gets infected
		IsInf <- 1
	}

	# Update the event list
	eventlist <- bind_rows(eventlist, tibble(t=t, 
		Sposx=Sposx, Sposy=Sposy, 
		Iposx=Iposx, Iposy=Iposy, 
		IsInf=IsInf))

}

# makevideo(eventlist, tstep=0.1, file="figures/episim.gif")