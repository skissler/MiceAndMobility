library(tidyverse) 
library(animation)
library(gganimate)

update_pos <- function(x, sigma, domainwidth){
	accept <- 0
	while(accept<1){
		xprop <- rnorm(1, x,sigma)
		if(xprop>0 & xprop<domainwidth){accept <- 1}
	}
	return(xprop)
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