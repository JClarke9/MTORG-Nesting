# Loading libraries -------------------------------------------------------

library(wesanderson)
library(patchwork)

windowsFonts(my_font = windowsFont("Gandi Sans"))

# Data import -------------------------------------------------------------


tot.rich <- read.csv("working/TotalRichness.csv",
                     row.names=1)

tot.simp <- read.csv("working/SimpsonTotal.csv",
                     row.names=1)


tot.rich$Year <- as.factor(tot.rich$Year)
tot.rich$Treat <- as.factor(tot.rich$Treat)

tot.simp$Year <- as.factor(tot.simp$Year)
tot.simp$Treat <- as.factor(tot.simp$Treat)

# Creating violin plots ---------------------------------------------------


rich.viol <- ggplot(tot.rich, 
                    aes(fill=Year,
                        x=Treat, 
                        y=Richness)) +
  geom_violin(trim=FALSE,
              draw_quantiles=c(.25, 
                               .75),
              colour= "black",
              size=1.0) +
  scale_fill_manual(values=c("#798E87", "#CCC591")) +
  theme(plot.title = element_text(family="my_font",
                                  hjust=.5), 
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="transparent",                     # make the interior background transparent
                                        colour = NA), 
        plot.background = element_rect(fill="transparent",                      # make the outer background transparent
                                       colour=NA), 
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.x = element_text(family="my_font",
                                   size=13, 
                                   colour = "black"),                             # color the axis text
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.ticks = element_line(colour = "black"),
        text=element_text(family="my_font",
                          size=13,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  labs(title = "Avian Nesting Richness", 
       x = NULL, 
       y = "Species Richness") +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", 
                            "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0,20))

rich.viol

simp.viol <- ggplot(tot.simp, 
                    aes(fill=Year,
                        x=Treat, 
                        y=simpson)) +
  geom_violin(trim=FALSE,
              draw_quantiles=c(.25, 
                               .75),
              colour= "black",
              size=1.0) +
  scale_fill_manual(values=c("#798E87", "#CCC591")) +
  theme(plot.title = element_text(family="my_font",
                                  hjust=.5), 
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="transparent",                     # make the interior background transparent
                                        colour = NA), 
        plot.background = element_rect(fill="transparent",                      # make the outer background transparent
                                       colour=NA), 
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.x = element_text(family="my_font",
                                   size=13, 
                                   colour = "black"),                             # color the axis text
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.ticks = element_line(colour = "black"),
        text=element_text(family="my_font",
                          size=13,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "right") +                                             # remove the legend
  labs(title = "Avian Nesting Diversity", 
       x = NULL, 
       y = "Simpson's Diversity") +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", 
                            "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0,1))

simp.viol

# Combining violin plots --------------------------------------------------

birds.viol <- (rich.viol/simp.viol)

birds.viol <- birds.viol + plot_annotation(theme=theme(panel.grid.major = element_blank(), # remove the vertical grid lines
                                                       panel.grid.minor = element_blank(), # remove the horizontal grid lines
                                                       panel.background = element_rect(fill=NA, # make the interior background transparent
                                                                                       colour = NA), 
                                                       plot.background = element_rect(fill =NA, # make the outer background transparent
                                                                                      colour = NA), 
                                                       axis.line = element_line(colour = "black"), # color the x and y axis
                                                       axis.text = element_text(colour = "black"), # color the axis text
                                                       axis.ticks = element_line(colour = "black"),
                                                       text=element_text(colour = "black"),# change the color of the axis titles
                                                       legend.position = "none"))

birds.viol

ggsave(birds.viol, 
       filename = "outputs/figs/AvianDiversityViolin.png",  
       dpi = "print", 
       bg = "transparent",
       width = 6,
       height = 7)
```