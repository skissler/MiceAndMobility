library(tidyverse) 

N <- 1000

# tau will be a time-localization parameter
# eta will be the individual R0

Finf_delta <- function(t, tinf, tau, eta){
	if(t>=(tinf+tau)){
		return(eta)
	} else {
		return(0)
	}
}
Finf_delta <- Vectorize(Finf_delta)

Finf_delta_v <- function(t, tinfvec, tauvec, etavec){
	tibble(tinf=tinfvec, tau=tauvec, eta=etavec)
}

tinfvec <- rep(Inf, N)
tinfvec[ceiling(runif(1)*N)] <- 0
tauvec <- rexp(N, 1/5)
etavec <- rep(2, N)

t <- 0
tstep <- 1
tmax <- 20
while(t<tmax){
	
	# print(paste0("told: ",t))
	# print(paste0("tnew: ",t+tstep))
	# print(paste0("tinf + tau: ",(tinfvec+tauvec)[(tinfvec+tauvec)<Inf]))

	rate <- (tibble(tnew=t+tstep, tinf=tinfvec, tau=tauvec, eta=etavec) %>% 
		mutate(F=Finf_delta(tnew,tinf,tau,eta)) %>% 
		pull(F) %>% 
		sum()) - 
	(tibble(told=t, tinf=tinfvec, tau=tauvec, eta=etavec) %>% 
		mutate(F=Finf_delta(told,tinf,tau,eta)) %>% 
		pull(F) %>% 
		sum())
	# print(paste0("rate: ",rate))

	ninf <- rpois(1,rate)
	# print(paste0("ninf: ",ninf))

	if(ninf>0){
		draw <- min(runif(ninf))
		# print(paste0("draw: ",draw))

		tspikevec <- tinfvec + tauvec
		qualifyingspikes <- sort(tspikevec[tspikevec>=t & tspikevec<(t+tstep)])
		# print(paste0("qualifyingspikes: ",qualifyingspikes))
		tinf <- qualifyingspikes[ceiling(draw*length(qualifyingspikes))]
		# print(paste0("tinf: ",tinf))
		newid <- sample(((1:N)[tinfvec>(2*tmax)]),1)
		# print(paste0("newid: ",newid))
		tinfvec[newid] <- tinf
		t <- tinf
		if(max(tinfvec)<Inf){break()}
	} else {
		t <- t+tstep
	}

	# print("-------------------------")


}


getunder <- function(t){
	return(sum(as.numeric(tinfvec<t)))
}
getunder <- Vectorize(getunder)

fig_cuminf <- tibble(t=seq(from=0,to=tmax,by=0.1)) %>% 
	mutate(cuminf=getunder(t)) %>% 
	ggplot(aes(x=t, y=cuminf)) + 
		geom_point(size=0.2) + 
		geom_line() + 
		scale_y_continuous(limits=c(0,N)) + 
		theme_classic()


