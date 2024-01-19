# Loading libraries -------------------------------------------------------

library(wesanderson)
library(cowplot)
library(tidyverse)

windowsFonts(my_font = windowsFont("Arial"))

# Data import -------------------------------------------------------------

tot.rich <- read.csv("working/TotalRichness.csv")

tot.simp <- read.csv("working/SimpsonTotal.csv")

tot.rich$Year <- factor(tot.rich$Year,
                        levels = c("2021", "2022", "2023"))
tot.rich$Treat <- factor(tot.rich$Treat,
                         levels = c("Rest", "Moderate", "Full", "Heavy"))

tot.simp$Year <- factor(tot.simp$Year,
                        levels = c("2021", "2022", "2023"))
tot.simp$Treat <- factor(tot.simp$Treat,
                         levels = c("Rest", "Moderate", "Full", "Heavy"))


# Creating Boxplots -------------------------------------------------------


(box.rich <- ggplot(tot.rich,                                                   # select the data to graph
                    aes(fill = Year,
                        x = Treat,                                              # define the x axis
                        y = Richness)) +
   geom_boxplot(colour = "#273046",                                             # create black outlines around the bar plot
                size = 1,
                notch = FALSE, 
                outlier.color = NULL,
                outlier.size = NA,
                outlier.shape = NA,
                fatten = 1) +
   scale_fill_manual(values = c("#798E87", "#CCC591", "brown")) +               # select the color for group 2
   theme(plot.title = element_text(family = "my_font",                          # select the font for the title
                                   size = 30,
                                   hjust = .5),
         panel.grid.major = element_blank(),                                    # remove the vertical grid lines
         panel.grid.minor = element_blank(),                                    # remove the horizontal grid lines
         panel.background = element_rect(fill = NA,                             # make the interior background transparent
                                         colour = NA),                          # remove any other colors
         plot.background = element_rect(fill = NA,                              # make the outer background transparent
                                        colour = NA),                           # remove any other colors
         axis.line = element_line(colour = "black"),                            # color the x and y axis
         axis.text.y = element_text(size = 20, colour = "black"),               # color the axis text
         axis.text.x = element_text(size = 20, colour = "black"),               # color the axis text
         axis.ticks = element_line(colour = "black"),                           # change the colors of the axis tick marks
         text=element_text(size = 24,                                           # change the size of the axis titles
                           colour = "black"),                                   # change the color of the axis titles
         legend.position = "none") +                                           # remove the legend
   labs(title = "Breeding Bird Richness", 
        x = NULL, 
        y = "Species Richness") +                                               # changing axis titles
   coord_cartesian(xlim = NULL, 
                   ylim = c(0,14)))


(box.simp <- ggplot(tot.simp,                                                  # select the data to graph
                    aes(fill=Year,
                        x=Treat,                                                # define the x axis
                        y=simpson)) +
    geom_boxplot(colour= "#273046",                                                     # create black outlines around the bar plot
                 size=1,
                 notch = FALSE, 
                 outlier.color = NULL,
                 outlier.size = NA,
                 outlier.shape = NA,
                 fatten=1) +
    scale_fill_manual(values=c("#798E87", "#CCC591", "brown")) +                                   # select the color for group 2
    theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                    size=30,
                                    hjust=.5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill=NA,                     # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill=NA,                      # make the outer background transparent
                                         colour=NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_text(size=20, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size=20, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size=24,                                              # change the size of the axis titles
                            colour = "black"),                                    # change the color of the axis titles
          legend.position = "right") +                                             # remove the legend
    labs(title = "Breeding Bird Diversity", 
         x = NULL, 
         y = "Simpson's Diversity") +                                               # changing axis titles
    coord_cartesian(xlim=NULL, ylim=c(0.3,0.9)))


# Combining Plots ---------------------------------------------------------


(birds.div <- plot_grid(box.rich, box.simp))
