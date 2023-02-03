# Loading libraries -------------------------------------------------------

library(vegan)
library(tidyverse)
library(wesanderson)

# Data import -------------------------------------------------------------

birds21 <- read.csv("~/Git/NDSU/Avian Community Analysis/WorkingData/birds21.csv", 
                    header = TRUE, 
                    row.names =1)
birds22 <- read.csv("~/Git/NDSU/Avian Community Analysis/WorkingData/birds22.csv", 
                    header = TRUE, 
                    row.names =1)

windowsFonts(my_font = windowsFont("Gandhi Sans"))                          # downloading a font to be used for my ordination

# Data wrangling ----------------------------------------------------------

species21 <- aggregate(birds21[2:25],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds21$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the simpson diversity for each of the
# grazing intensities for 2021
simpson.group21 <- diversity(species21[2:25],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="simpson",                                           # using the simpson index
                             MARGIN = 1)                                                # not positive what this does

summary(simpson.group21)                                                        # show the output of the simpson diversity index
simpson.group21 <- as.data.frame(simpson.group21)                               # assign the simpson diversity index to a data frame
simpson.group21$Treat <- species21$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
simpson.group21                                                                 # show the output for the simpson diversity

# I recalculated simpson diversity for each paddock to
# determine whether diversity was significantly different
# between grazing intensities

simpson21 <- diversity(birds21[2:25],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       index="simpson",                                           # using the simpson index
                       MARGIN = 1)                                                # not positive what this does

summary(simpson21)                                                              # show the output of the simpson diversity index
simpson21 <- as.data.frame(simpson21)                                           # assign the simpson diversity index to a data frame
simpson21$Treat <- birds21$cTreat                                               # create a treatment column in the new data frame to use for an ANOVA
simpson21                                                                       # show the output for the simpson diversity


aov.simp21 <- aov(simpson21~Treat,                                              # select the variables to compare with an ANOVA
                  simpson21)                                                      # select the data frame
summary(aov.simp21)                                                             # show the output of the ANOVA

# Calculate species diversity ---------------------------------------------

species22 <- aggregate(birds22[2:28],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds22$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the simpson diversity for each of the
# grazing intensities for 2022
simpson.group22 <- diversity(species22[2:28],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="simpson",                                           # using the simpson index
                             MARGIN = 1)                                                # not positive what this does

summary(simpson.group22)                                                        # show the output of the simpson diversity index
simpson.group22 <- as.data.frame(simpson.group22)                               # assign the simpson diversity index to a data frame
simpson.group22$Treat <- species22$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
simpson.group22                                                                 # show the output for the simpson diversity

# I recalculated simpson diversity for each paddock to
# determine whether diversity was significantly different
# between grazing intensities

simpson22 <- diversity(birds22[2:28],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       index="simpson",                                           # using the simpson index
                       MARGIN = 1)                                                # not positive what this does

summary(simpson22)                                                              # show the output of the simpson diversity index
simpson22 <- as.data.frame(simpson22)                                           # assign the simpson diversity index to a data frame
simpson22$Treat <- birds22$cTreat                                               # create a treatment column in the new data frame to use for an ANOVA
simpson22                                                                       # show the output for the simpson diversity


aov.simp22 <- aov(simpson22~Treat,                                              # select the variables to compare with an ANOVA
                  simpson22)                                                      # select the data frame
summary(aov.simp22)                                                             # show the output of the ANOVA

# Creating diversity data frame -------------------------------------------

simp21 <- rename(simpson21, 
                 "simpson"="simpson21")
simp21$Year <- 2021
simp22 <- rename(simpson22, 
                 "simpson"="simpson22")
simp22$Year <- 2022
simp.total <- rbind(simp21, simp22)

simp.total$Year <- as.factor(simp.total$Year)

simp.aov <- aov(simpson~ Treat/Year, simp.total)
summary(simp.aov)

TukeyHSD(simp.aov)

simp.div <- full_join(simpson.group21, 
                      simpson.group22, 
                      by = "Treat") |> 
  relocate("Treat", 
           .before="simpson.group21") |> 
  rename("Grazing Intensity" = "Treat", 
         "2021 Diversity" = "simpson.group21", 
         "2022 Diversity" = "simpson.group22")

write.csv(simp.total, "outputs/figs/SimpsonTotal.csv")

# Creating diversity plot -------------------------------------------------


simp.box <- ggplot(simp.total,                                                  # select the data to graph
                   aes(fill=Year,
                       x=Treat,                                                # define the x axis
                       y=simpson)) +
  geom_boxplot(colour= "white",                                            # create black outlines around the bar plot
               size=2.5) +  
  geom_boxplot(colour= "#273046",                                                     # create black outlines around the bar plot
               size=2,
               notch = FALSE, 
               outlier.color = NULL,
               outlier.size = NA,
               outlier.shape = NA,
               fatten=1) +
  scale_fill_manual(values=c("#798E87", "#CCC591")) +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=40,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill=NA,                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill=NA,                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "white"),                             # color the x and y axis
        axis.text.y = element_text(size=28, colour = "white"),                    # color the axis text
        axis.text.x = element_text(size=28, colour = "white"),
        axis.ticks = element_line(colour = "white"),                            # change the colors of the axis tick marks
        text=element_text(size=32,                                              # change the size of the axis titles
                          colour = "white"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  labs(title = "Avian Nesting Diversity", 
       x = NULL, 
       y = "Simpson's Diversity") +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", 
                            "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0.3,0.9))

simp.box

ggsave(simp.box, 
       filename = "outputs/figs/SimpsonsBox.png",  
       dpi = "print", 
       bg = "transparent",
       width = 15,
       height = 10)
```