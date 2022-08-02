library(tidyverse) 
library(animation)
library(gganimate)

# Update an agent's position in 1D using a N(0,sigma)
update_pos <- function(x, sigma, domainwidth){
	accept <- 0
	while(accept<1){
		xprop <- rnorm(1,x,sigma)
		if(xprop>0 & xprop<domainwidth){accept <- 1}
	}
	return(xprop)
}

# Simulate a two-agent epidemic
#   k: kernel intercept
#   scale: the kernel scale paramter, either d*, alpha, or phi (see below)
#   mu: movement rate (units 1/time) 
#   sigma: sd of movement distance 
#   domainwidth: width of the square domain 
#   kernel: specification of transmission kernel, one of: 
#     stp: k if d ≤ d*, 0 otherwise.
#     pow: k/(1 + d^alpha)
#     exp: k e^(-phi*d)
episim <- function(k=1, scale=1, mu=1, sigma=1, domainwidth=10, kernel="exp"){

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

  while(IsInf<1){

  # Calculate key parameters 
  d <- sqrt((Sposx-Iposx)^2 + (Sposy-Iposy)^2)
  if(kernel=="exp"){
    lambda <- k*exp(-scale*d)  
  } else if(kernel=="pow"){
    lambda <- k/(1 + d^scale)
  } else if(kernel=="stp"){
    lambda <- ifelse(d <= scale, k, 0)
  } else {
    stop("Unspecified kernel")
  }
  
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

  return(eventlist)

}


# episim_exp <- function(k=1, phi=1, mu=1, sigma=1, domainwidth=10){

#   # Set initial conditions 
#   t <- 0
#   Sposx <- runif(1)*domainwidth
#   Sposy <- runif(1)*domainwidth
#   Iposx <- runif(1)*domainwidth
#   Iposy <- runif(1)*domainwidth
#   IsInf <- 0
#   eventlist <- tibble(t=t, 
#     Sposx=Sposx, Sposy=Sposy, 
#     Iposx=Iposx, Iposy=Iposy, 
#     IsInf=IsInf)

#   while(IsInf<1){

#   # Calculate key parameters 
#   d <- sqrt((Sposx-Iposx)^2 + (Sposy-Iposy)^2)
#   lambda <- k*exp(-phi*d)
#   cumrate <- 2*mu + lambda

#   # Simulate the time of the next event
#   t <- t+rexp(1, rate=cumrate)

#   # Draw the event 
#   eventdraw <- runif(1)
#   if(eventdraw < mu/cumrate){ # S moves
#     Sposx <- update_pos(Sposx, sigma, domainwidth)
#     Sposy <- update_pos(Sposy, sigma, domainwidth)
#   } else if(eventdraw < 2*mu/cumrate){ # I moves
#     Iposx <- update_pos(Iposx, sigma, domainwidth)
#     Iposy <- update_pos(Iposy, sigma, domainwidth)
#   } else { # S gets infected
#     IsInf <- 1
#   }

#   # Update the event list
#   eventlist <- bind_rows(eventlist, tibble(t=t, 
#     Sposx=Sposx, Sposy=Sposy, 
#     Iposx=Iposx, Iposy=Iposy, 
#     IsInf=IsInf))
  
#   }

#   return(eventlist)

# }

makevideo <- function(eventlist, tstep, file="anim.gif"){
  eventlist_mod <- eventlist %>% 
    mutate(t=ceiling(t*1/tstep)*tstep) %>% 
    group_by(t) %>% 
    slice(n())

  eventlist_reg <- tibble(t=seq(from=0, to=(max(eventlist$t)+tstep), by=tstep)) %>% 
    left_join(eventlist_mod, by="t") %>% 
    fill(Sposx, Sposy, Iposx, Iposy, IsInf)

  posplot <- ggplot(data=eventlist_reg) + 
    geom_point(aes(x=Sposx, y=Sposy), col="blue", size=4) + 
    geom_point(aes(x=Iposx, y=Iposy), col="red", size=4) + 
    scale_x_continuous(limits=c(0,domainwidth)) + 
    scale_y_continuous(limits=c(0,domainwidth)) + 
    labs(x="x",y="y") + 
    theme_classic()
    
  posplot.animation <- posplot + 
    transition_time(t) # + 
    # shadow_wake(wake_length=0.1)

  animate(posplot.animation, height = 500, width = 500, fps = 20, duration = 20,
        end_pause = 60, res = 100)

  anim_save(paste0(file))

}



















# graph1 <- gapminder %>%
#   ggplot(aes(x=gdpPercap, y=lifeExp, color=continent, size=pop)) +
#   geom_point(alpha = 0.7, stroke = 0) +
#   theme_minimal() +
#   scale_x_log10() +
#   labs(title = "Life Expectancy vs GDP",
#        x = "Income",
#        y = "Life Expectancy",
#        color = "Continent",
#        caption = "Source: Gapminder") +
#   theme(axis.title = element_text(),
#         legend.text=element_text(size=10)) +
#   scale_color_brewer(palette = "Set2")

#  graph1.animation = graph1 +
#   transition_time(year) +
#   labs(subtitle = "Year: {frame_time}") +
#   shadow_wake(wake_length = 0.1)

#  animate(graph1.animation, height = 500, width = 800, fps = 20, duration = 10,
#         end_pause = 60, res = 100)

#  anim_save("gapminder graph.gif")