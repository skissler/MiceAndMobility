# =============================================================================
# Import and define 
# =============================================================================

library(tidyverse)
library(purrr)

N <- 50

beta <- 1/2 
gamma <- 1/5
ninf <- rpois(N, beta/gamma)

tstar <- rexp(N, gamma)

kwidth <- 0.25

# =============================================================================
# Define the true kernels and related quantities 
# =============================================================================

inftimes <- pmap(list(ninf, tstar), function(x,y){return(c(ninf=x,tstar=y))}) %>% 
	map(~ .[["tstar"]] + kwidth*runif(.[["ninf"]]))

truekernels <- as.list(tstar) %>% 
	imap(~ tibble(id=.y, x=c(0, .x, .x, .x+kwidth, .x+kwidth, Inf), y=c(0, 0, beta/gamma/kwidth, beta/gamma/kwidth, 0, 0))) %>% 
	bind_rows()

# Person-wise mean
means <- map(inftimes, mean)

# Overall per-person SD (gotten by subtracting each mean )
sd <- inftimes %>% 
	map(~ .-mean(.)) %>% 
	unlist() %>% 
	sd()

# =============================================================================
# Plot various inferences of the infectiousness kernel 
# =============================================================================

xvals <- seq(from=0, to=30, by=0.01)

infkernels_indivmean <- means %>% 
	imap(~ tibble(id=.y, x=xvals, y=dnorm(xvals, .x, sd))) %>% 
	bind_rows() %>% 
	mutate(y=y*beta/gamma)

infkernels_sharedmean <- inftimes %>% 
	unlist() %>% 
	(function(x){return(tibble(id=0,x=xvals,mean=mean(x),sd=sd(x)))}) %>% 
	# mutate(y=N*beta/gamma*dnorm(x,mean,sd)) %>% 
	mutate(y=N*beta/gamma*dexp(x,1/mean)) %>% 
	select(id,x,y)

get_noisy_onsets <- function(inftimes, epsilon){

	N <- length(inftimes)

	# get onsets for the index cases 

	inftimes_noisy <- inftimes %>% 
		# Noisy observation of the child cases: 
		map(~ unlist(lapply(., rnorm, n=1, sd=epsilon))) %>% 
		# Noisy observation of the parent cases: 
		map(~ .+rnorm(length(.),0,epsilon))

	return(inftimes_noisy)

}

inftimes_noisy <- get_noisy_onsets(inftimes, 0.75)

infkernels_indivmean_noisy <- map(inftimes_noisy, mean) %>% 
	imap(~ tibble(id=.y, x=xvals, y=dnorm(xvals, .x, sd))) %>% 
	bind_rows() %>% 
	mutate(y=y*beta/gamma)


fig_truekernel <- truekernels %>% 
	ggplot(aes(x=x, y=y, group=factor(id))) + 
		geom_line(col="black", alpha=0.5) + 
		theme_classic() + 
		scale_y_continuous(limits=c(0, beta/gamma/kwidth*1.5)) + 
		labs(x="Time since infection", y="Infectiousness") + 
		theme(text=element_text(size=9))

ggsave(fig_truekernel, file="figures/truekernel.pdf", width=3.5, height=3.5/1.6)		
	

fig_kernelcomp <- truekernels %>% 
	ggplot(aes(x=x, y=y, group=factor(id))) + 
		geom_line() + 
		geom_line(data=infkernels_indivmean, col="gray", alpha=0.6) + 
		geom_line(data=infkernels_indivmean_noisy, col="red", alpha=0.6) + 
		geom_line(data=infkernels_sharedmean, col="dodgerblue", alpha=0.4, size=1) + 
		theme_classic() 

fig_inftimehist <- tibble(tinf=unlist(inftimes_noisy)) %>% 
	ggplot(aes(x=abs(tinf))) + 
		geom_histogram(binwidth=1, fill="white", col="black") + 
		theme_classic() + 
		labs(x="Observed generation time", y="Count") + 
		theme(text=element_text(size=9))

ggsave(fig_inftimehist, file="figures/inftimehist.pdf", width=3.5, height=3.5/1.6)		


fig_tinf_noisy <- inftimes_noisy %>% 
	imap(~ tibble(id=.y, tinf=.x)) %>% 
	bind_rows() %>% 
	group_by(id) %>% 
	mutate(meantinf=mean(tinf)) %>% 
	arrange(meantinf, id, tinf) %>% 
	(function(x){
		orderedids <- x %>% 
			group_by(id) %>% 
			slice(1) %>% 
			ungroup() %>% 
			arrange(meantinf) %>% 
			mutate(newid=1:n()) %>% 
			select(id, newid) 
		out <- x %>% 
			left_join(orderedids, by="id") 
		return(out)
		}) %>% 
	ggplot(aes(x=abs(tinf), y=newid)) + 
		geom_point() + 
		theme_classic() + 
		labs(x="Observed infection time", y="Index case ID") + 
		theme(text=element_text(size=9))

ggsave(fig_tinf_noisy, file="figures/tinf_noisy.pdf", width=3.5, height=3.5/1.6)	



