# Loading libraries -------------------------------------------------------

library(wesanderson)
library(cowplot)
library(tidyverse)

windowsFonts(my_font = windowsFont("Arial"))

# Data import -------------------------------------------------------------

tot.rich <- read.csv("working/TotalRichness.csv")

tot.simp <- read.csv("working/SimpsonTotal.csv")

tot.rich$Year <- as.factor(tot.rich$Year)
tot.rich$Treat <- as.factor(tot.rich$Treat)

tot.simp$Year <- as.factor(tot.simp$Year)
tot.simp$Treat <- as.factor(tot.simp$Treat)

# Creating violin plots ---------------------------------------------------


my_theme <- theme(plot.title = element_text(family = "my_font",
                                              hjust = .5), 
                  panel.grid.major = element_blank(),                           # remove the vertical grid lines
                  panel.grid.minor = element_blank(),                           # remove the horizontal grid lines
                  panel.background = element_rect(fill = "transparent",           # make the interior background transparent
                                                  colour = NA),
                  plot.background = element_rect(fill = "transparent",            # make the outer background transparent
                                                 colour = NA),
                  axis.line = element_line(colour = "black"),                   # color the x and y axis
                  axis.text.x = element_text(family = "my_font",
                                             size = 20, 
                                             colour = "black"),                 # color the axis text
                  axis.text.y = element_text(colour = "black",
                                             size = 20),
                  axis.ticks = element_line(colour = "black"),
                  text = element_text(family = "my_font",
                                    size = 20,                                    # change the size of the axis titles
                                    colour = "black"),                          # change the color of the axis titles
                  legend.position = "right")

rich.viol <- ggplot(tot.rich,
                    aes(fill=Year,
                        x=Treat, 
                        y=Richness)) +
  geom_violin(trim=FALSE,
              draw_quantiles=c(.25, 
                               .75),
              colour= "black",
              linewidth=.5) +
  scale_fill_manual(values=c('#A2A4A2',
                             '#717F5B',
                             '#D4A634')) +
  my_theme +
  labs(title = "Avian Nesting Richness", 
       x = NULL, 
       y = "Species Richness") +
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", 
                            "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest"))

rich.viol

# old colors - "#798E87", "#CCC591"

simp.viol <- ggplot(tot.simp, 
                    aes(fill=Year,
                        x=Treat, 
                        y=simpson)) +
  geom_violin(trim=FALSE,
              draw_quantiles=c(.25, 
                               .75),
              colour= "black",
              linewidth=.5) +
  scale_fill_manual(values=c('#A2A4A2',
                             '#717F5B',
                             '#D4A634')) +
  my_theme +
  labs(title = "Avian Nesting Diversity", 
       x = NULL, 
       y = "Simpson's Diversity") +
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", 
                            "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest"))

simp.viol

# Combining violin plots --------------------------------------------------


birds.viol <- plot_grid(rich.viol, 
                        simp.viol, 
                        nrow = 2, 
                        ncol = 1)

birds.viol

ggsave(birds.viol, 
       filename = "outputs/figs/AvianDiversityViolin.png",  
       dpi = 300,
       bg = "white",
       height = 12,
       width = 15,
       units = "in")
