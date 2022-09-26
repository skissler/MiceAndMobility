library(tidyverse) 

# Set the radius of the space 
radius <- 1

circle_df <- tibble(x=seq(from=-radius, to=radius, by=radius/100)) %>% 
	mutate(ytop=sqrt(radius^2-x^2)) %>% 
	mutate(ybot=-sqrt(radius^2-x^2))

# Set the number of individuals
nindiv <- 250
ninf0 <- 1

# Set the movement rate 
mu <- 1

# Define the infectiousness kernel 
k <- 1             # Force of infection when distance = 0
scl <- 1           # Scale parameter for the distance kernel

# Draw initial positions for everyone
rvals <- radius*sqrt(runif(nindiv))
thetavals <- 2*pi*runif(nindiv)
pos_df <- tibble(x=rvals*cos(thetavals), y=rvals*sin(thetavals), inf=0)
pos_df <- pos_df %>% mutate(id=1:n())
# Randomly draw initial infected
pos_df$inf[sample(1:nrow(pos_df),ninf0)] <- 1


# Calculate the force on each susceptible person
force_df <- pos_df %>% 
	split(.$inf) %>% 
	map(~ select(., -inf)) %>% 
	(function(x){full_join(x[[1]],x[[2]], by=character(), suffix=c("_S","_I"))}) %>% 
	# mutate(d=sqrt((x_S-x_I)^2 + (y_S-y_I)^2))
	mutate(force=k) %>% 
	group_by(id_S) %>% 
	summarise(force=sum(force)) %>% 
	select(id=id_S, force)





tempfig <- pos_df %>% 
	ggplot(aes(x=x, y=y, col=factor(inf))) + 
		geom_rect(xmin=-1.1*radius, xmax=1.1*radius, ymin=-1.1*radius, ymax=1.1*radius, col="#303841", fill="#303841") + 
		geom_point(size=1) + 
		geom_point(data=filter(pos_df, inf==1), size=3) + 
		scale_color_manual(values=c("0"="#A6ACBA", "1"="#DC343B"), guide="none") + 
		geom_line(data=circle_df, aes(x=x, y=ytop), col="#A6ACBA") + 
		geom_line(data=circle_df, aes(x=x, y=ybot), col="#A6ACBA") + 
		theme_void() 

# Determine the next event 
cumrate <- k






# Append the event to event_df 


