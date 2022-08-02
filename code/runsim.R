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
scl <- 1           # Scale parameter for the distance kernel

# Plot the kernels: 
dvals <- seq(from=0, to=domainwidth, by=0.1)
kerneldf <- tibble(d=dvals) %>% 
	mutate(lexp=k*exp(-scl*d)) %>%
	mutate(lpow=k/(1+d^scl)) %>% 
	mutate(lstp=ifelse(d<=scl,k,0))
fig_kernels <- kerneldf %>% 
	pivot_longer(-d) %>% 
	ggplot(aes(x=d, y=value, col=name)) + 
		geom_line() + 
		scale_color_manual(values=c("lexp"="blue","lpow"="red","lstp"="black"), labels=c("lexp"="Exponential","lpow"="Power","lstp"="Step")) + 
		theme_classic() + 
		labs(x="Distance", y="Force of infection (λ)") + 
		theme(legend.title=element_blank())
ggsave(fig_kernels, file="notes/kernels.png", width=5, height=3)


# Run 1000 sims and extract the time and position of infection:
tinfdf_exp <- lapply(1:1000, function(x){
	episim(k=k, scale=scl, mu=mu, sigma=sigma, domainwidth=domainwidth, kernel="exp") %>% 
	filter(IsInf==1)}) %>% bind_rows()

tinfdf_pow <- lapply(1:1000, function(x){
	episim(k=k, scale=scl, mu=mu, sigma=sigma, domainwidth=domainwidth, kernel="pow") %>% 
	filter(IsInf==1)}) %>% bind_rows()

tinfdf_stp <- lapply(1:1000, function(x){
	episim(k=k, scale=scl, mu=mu, sigma=sigma, domainwidth=domainwidth, kernel="stp") %>% 
	filter(IsInf==1)}) %>% bind_rows()

# makevideo(eventlist, tstep=0.1, file="figures/episim.gif")

fig_tinf_exp <- ggplot(data=tinfdf_exp, aes(x=t)) + 
	geom_histogram(aes(y=..density..),binwidth=10,fill="white",col="darkgray") +
	labs(title="Exponential kernel", x="Time of infection", y="Density") + 
	theme_classic()

fig_tinf_pow <- ggplot(data=tinfdf_pow, aes(x=t)) + 
	geom_histogram(aes(y=..density..),binwidth=10,fill="white",col="darkgray") +
	labs(title="Power kernel", x="Time of infection", y="Density") + 
	theme_classic()

fig_tinf_stp <- ggplot(data=tinfdf_stp, aes(x=t)) + 
	geom_histogram(aes(y=..density..),binwidth=10,fill="white",col="darkgray") +
	labs(title="Step kernel", x="Time of infection", y="Density") + 
	theme_classic()

ggsave(fig_tinf_exp, file="notes/tinf_exp.png", width=5, height=3)
ggsave(fig_tinf_pow, file="notes/tinf_pow.png", width=5, height=3)
ggsave(fig_tinf_stp, file="notes/tinf_stp.png", width=5, height=3)



fig_dinf_exp <- ggplot(data=tinfdf_exp, aes(x=sqrt((Sposx-Iposx)^2+(Sposy-Iposy)^2))) + 
	geom_histogram(aes(y=..density..),binwidth=0.25,fill="white",col="darkgray") +
	scale_x_continuous(limits=c(0,domainwidth)) + 
	labs(title="Exponential kernel", x="Distance at time of infection", y="Density") + 
	theme_classic()

fig_dinf_pow <- ggplot(data=tinfdf_pow, aes(x=sqrt((Sposx-Iposx)^2+(Sposy-Iposy)^2))) + 
	geom_histogram(aes(y=..density..),binwidth=0.25,fill="white",col="darkgray") +
	scale_x_continuous(limits=c(0,domainwidth)) + 
	labs(title="Power kernel", x="Distance at time of infection", y="Density") + 
	theme_classic()

fig_dinf_stp <- ggplot(data=tinfdf_stp, aes(x=sqrt((Sposx-Iposx)^2+(Sposy-Iposy)^2))) + 
	geom_histogram(aes(y=..density..),binwidth=0.25,fill="white",col="darkgray") +
	scale_x_continuous(limits=c(0,domainwidth)) + 
	labs(title="Step kernel", x="Distance at time of infection", y="Density") + 
	theme_classic()

ggsave(fig_dinf_exp, file="notes/dinf_exp.png", width=5, height=3)
ggsave(fig_dinf_pow, file="notes/dinf_pow.png", width=5, height=3)
ggsave(fig_dinf_stp, file="notes/dinf_stp.png", width=5, height=3)
