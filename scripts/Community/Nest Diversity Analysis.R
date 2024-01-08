# Loading libraries -------------------------------------------------------

library(vegan)
library(tidyverse)
library(wesanderson)

# Data import -------------------------------------------------------------

birds21 <- read.csv("working/richness21.csv", 
                    header = TRUE, 
                    row.names =1)

birds22 <- read.csv("working/richness22.csv", 
                    header = TRUE, 
                    row.names =1)

birds23 <- read.csv("working/richness23.csv", 
                    header = TRUE, 
                    row.names =1)

windowsFonts(my_font = windowsFont("Gandhi Sans"))                          # downloading a font to be used for my ordination

# Data wrangling ----------------------------------------------------------

species21 <- aggregate(birds21[3:23],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds21$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the simpson diversity for each of the
# grazing intensities for 2021
simpson.group21 <- diversity(species21[2:22],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="simpson",                                           # using the simpson index
                             MARGIN = 1)                                                # not positive what this does

summary(simpson.group21)                                                        # show the output of the simpson diversity index
simpson.group21 <- as.data.frame(simpson.group21)                               # assign the simpson diversity index to a data frame
simpson.group21$Treat <- species21$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
simpson.group21                                                                 # show the output for the simpson diversity

# I recalculated simpson diversity for each paddock to
# determine whether diversity was significantly different
# between grazing intensities

simpson21 <- diversity(birds21[3:23],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       index="simpson",                                           # using the simpson index
                       MARGIN = 1)                                                # not positive what this does

summary(simpson21)                                                              # show the output of the simpson diversity index
simpson21 <- as.data.frame(simpson21)                                           # assign the simpson diversity index to a data frame
simpson21$Treat <- birds21$cTreat                                               # create a treatment column in the new data frame to use for an ANOVA
simpson21$Paddock <- birds21$Paddock
simpson21                                                                       # show the output for the simpson diversity

simpson21$Treat <- factor(simpson21$Treat,
                          levels = c("Rest",
                                     "Moderate",
                                     "Full",
                                     "Heavy"))

simpson21$Paddock <- factor(simpson21$Paddock,
                            levels = c("1", "2", "3", "4",
                                       "5", "6", "7", "8",
                                       "9", "10", "11", "12",
                                       "13", "14", "15", "16"))

aov.simp21 <- aov(simpson21~Treat,                                              # select the variables to compare with an ANOVA
                  simpson21)
summary(aov.simp21) 

aov.simp21 <- glmmTMB(simpson21~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      simpson21,
                      family = Gamma(link = "log"))
summary(aov.simp21)                                                             # show the output of the ANOVA

simulationOutput <- simulateResiduals(aov.simp21)
plot(simulationOutput)
diagnose(aov.simp21)

# Calculate species diversity ---------------------------------------------

species22 <- aggregate(birds22[3:23],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds22$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the simpson diversity for each of the
# grazing intensities for 2022
simpson.group22 <- diversity(species22[2:22],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="simpson",                                           # using the simpson index
                             MARGIN = 1)                                                # not positive what this does

summary(simpson.group22)                                                        # show the output of the simpson diversity index
simpson.group22 <- as.data.frame(simpson.group22)                               # assign the simpson diversity index to a data frame
simpson.group22$Treat <- species22$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
simpson.group22                                                                 # show the output for the simpson diversity

# I recalculated simpson diversity for each paddock to
# determine whether diversity was significantly different
# between grazing intensities

simpson22 <- diversity(birds22[3:23],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       index="simpson",                                           # using the simpson index
                       MARGIN = 1)                                                # not positive what this does

summary(simpson22)                                                              # show the output of the simpson diversity index
simpson22 <- as.data.frame(simpson22)                                           # assign the simpson diversity index to a data frame
simpson22$Treat <- birds22$cTreat                                               # create a treatment column in the new data frame to use for an ANOVA
simpson22$Paddock <- birds22$Paddock
simpson22                                                                       # show the output for the simpson diversity

simpson22$Treat <- factor(simpson22$Treat,
                          levels = c("Rest",
                                     "Moderate",
                                     "Full",
                                     "Heavy"))

simpson22$Paddock <- factor(simpson22$Paddock,
                            levels = c("1", "2", "3", "4",
                                       "5", "6", "7", "8",
                                       "9", "10", "11", "12",
                                       "13", "14", "15", "16"))

aov.simp22 <- aov(simpson22~Treat,                                              # select the variables to compare with an ANOVA
                  simpson22)
summary(aov.simp22) 

aov.simp22 <- glmmTMB(simpson22~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      simpson22,
                      family = Gamma(link = "log"))
summary(aov.simp22)                                                             # show the output of the ANOVA

simulationOutput <- simulateResiduals(aov.simp22)
plot(simulationOutput)
diagnose(aov.simp22)

# Calculate species diversity ---------------------------------------------

species23 <- aggregate(birds23[3:24],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds23$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the simpson diversity for each of the
# grazing intensities for 2023
simpson.group23 <- diversity(species23[2:23],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="simpson",                                           # using the simpson index
                             MARGIN = 1)                                                # not positive what this does

summary(simpson.group23)                                                        # show the output of the simpson diversity index
simpson.group23 <- as.data.frame(simpson.group23)                               # assign the simpson diversity index to a data frame
simpson.group23$Treat <- species23$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
simpson.group23                                                                 # show the output for the simpson diversity

# I recalculated simpson diversity for each paddock to
# determine whether diversity was significantly different
# between grazing intensities

simpson23 <- diversity(birds23[3:24],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       index="simpson",                                           # using the simpson index
                       MARGIN = 1)                                                # not positive what this does

summary(simpson23)                                                              # show the output of the simpson diversity index
simpson23 <- as.data.frame(simpson23)                                           # assign the simpson diversity index to a data frame
simpson23$Treat <- birds23$cTreat                                               # create a treatment column in the new data frame to use for an ANOVA
simpson23$Paddock <- birds23$Paddock
simpson23                                                                       # show the output for the simpson diversity

simpson23$Treat <- factor(simpson23$Treat,
                          levels = c("Rest",
                                     "Moderate",
                                     "Full",
                                     "Heavy"))

simpson23$Paddock <- factor(simpson23$Paddock,
                            levels = c("1", "2", "3", "4",
                                       "5", "6", "7", "8",
                                       "9", "10", "11", "12",
                                       "13", "14", "15", "16"))

aov.simp23 <- aov(simpson23~Treat,                                              # select the variables to compare with an ANOVA
                  simpson23)
summary(aov.simp23) 

aov.simp23 <- glmmTMB(simpson23~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      simpson23,
                      family = Gamma(link = "log"))
summary(aov.simp23)                                                             # show the output of the ANOVA

simulationOutput <- simulateResiduals(aov.simp23)
plot(simulationOutput)
diagnose(aov.simp23)

# Creating diversity data frame -------------------------------------------

simp21 <- rename(simpson21, 
                 "simpson"="simpson21")
simp21$Year <- 2021

simp22 <- rename(simpson22, 
                 "simpson"="simpson22")
simp22$Year <- 2022

simp23 <- rename(simpson23, 
                 "simpson"="simpson23")
simp23$Year <- 2023

simp.total <- rbind(simp21, 
                    simp22, 
                    simp23)

simp.total$Year <- factor(simp.total$Year,
                          levels = c("2021",
                                     "2022",
                                     "2023"))

str(simp.total)

qqnorm(simp.total$simpson)
qqline(simp.total$simpson)
hist(simp.total$simpson)

simp.aov <- glmmTMB(simpson ~ Treat + (1|Year) + (1|Paddock), 
                simp.total,
                family = Gamma(link = "log"))
summary(simp.aov)

simp.aov <- glmmTMB(simpson ~ Treat + Year + (1|Paddock), 
                    simp.total,
                    family = Gamma(link = "log"))

summary(simp.aov)
diagnose(simp.aov)

simulationOutput <- simulateResiduals(simp.aov)
plot(simulationOutput)
diagnose(simp.aov)

simp.div <- full_join(simpson.group21, 
                      simpson.group22,
                      by = "Treat") |> 
  full_join(simpson.group23,
            by = "Treat") |> 
  relocate("Treat", 
           .before="simpson.group21") |> 
  rename("Grazing Intensity" = "Treat", 
         "2021 Diversity" = "simpson.group21", 
         "2022 Diversity" = "simpson.group22",
         "2023 Diversity" = "simpson.group23")

write.csv(simp.total, "outputs/figs/SimpsonTotal.csv")


# Creating diversity plot -------------------------------------------------


(simp.box <- ggplot(simp.total,                                                  # select the data to graph
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
         legend.position = "none") +                                             # remove the legend
   labs(title = "Breeding Bird Diversity", 
        x = NULL, 
        y = "Simpson's Diversity") +                                               # changing axis titles
   scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                    labels=c("Rest"="Rest", "Moderate"="Moderate", 
                             "Full"="Full", "Heavy"="Heavy"),
                    limits=c("Heavy", "Full", "Moderate", "Rest")) +
   coord_cartesian(xlim=NULL, ylim=c(0.3,0.9))

ggsave(simp.box, 
       filename = "outputs/figs/SimpsonsBox.png",  
       dpi = "retina", 
       bg = "white",
       width = 6,
       height = 6.5)
