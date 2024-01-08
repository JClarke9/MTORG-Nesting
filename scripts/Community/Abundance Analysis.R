# Loading libraries -------------------------------------------------------


library(vegan)
library(tidyverse)
library(RColorBrewer)
library(emmeans)
library(DHARMa)
library(glmmTMB)
library(cowplot)
library(wesanderson)


# Data import -------------------------------------------------------------


totals21 <- read.csv("working/totals21.csv", 
                    row.names = 1)                             # read in the data set

totals22 <- read.csv("working/totals22.csv", 
                    row.names = 1)                             # read in the data set

totals23 <- read.csv("working/totals23.csv", 
                    row.names = 1)                             # read in the data set

windowsFonts(my_font = windowsFont("Gandhi Sans"))


# Manipulating data -------------------------------------------------------

totals21.bar <- totals21 |> 
  filter(Group != "GEN") |> 
  group_by(Paddock,
           cTreat,
           Group) |> 
  summarize(Abundance = sum(Abundance)) |> 
  ungroup()

totals21.bar <- complete(totals21.bar,
                         nesting(Paddock,
                                 cTreat),
                         Group,
                         fill=list(Abundance=0))

totals22.bar <- totals22 |> 
  filter(Group != "GEN") |> 
  group_by(Paddock,
           cTreat,
           Group) |> 
  summarize(Abundance = sum(Abundance)) |> 
  ungroup()

totals22.bar <- complete(totals22.bar,
                         nesting(Paddock,
                                 cTreat),
                         Group,
                         fill=list(Abundance=0))

totals23.bar <- totals23 |> 
  filter(Group != "GEN") |> 
  group_by(Paddock,
           cTreat,
           Group) |> 
  summarize(Abundance = sum(Abundance)) |> 
  ungroup()

totals23.bar <- complete(totals23.bar,
                         nesting(Paddock,
                                 cTreat),
                         Group,
                         fill=list(Abundance=0))

totals21.bar$Year <- 2021
totals22.bar$Year <- 2022
totals23.bar$Year <- 2023

totals.bar <- rbind(totals21.bar,
                    totals22.bar,
                    totals23.bar) |> 
  filter(Group != "WET") |> 
  group_by(Paddock,
           cTreat,
           Group,
           Year) |> 
  summarize(Abundance = sum(Abundance)) |> 
  ungroup()

totals.bar$cTreat <- factor(totals.bar$cTreat,
                            levels = c("Rest",
                                       "Moderate",
                                       "Full",
                                       "Heavy"))


totals.bar$Paddock <- factor(totals.bar$Paddock,
                             levels = c("1", "2", "3", "4",
                                        "5", "6", "7", "8",
                                        "9", "10", "11", "12",
                                        "13", "14", "15", "16"))

totals.bar$Year <- factor(totals.bar$Year,
                          levels = c("2021",
                                     "2022",
                                     "2023"))

FACtotals.bar <- filter(totals.bar, 
                        Group == "FAC")
OBLtotals.bar <- filter(totals.bar, 
                        Group == "OBL")
FAC21.bar <- filter(totals.bar, 
                    Group == "FAC" & Year == 2021)
OBL21.bar <- filter(totals.bar, 
                    Group == "OBL" & Year == 2021)
FAC22.bar <- filter(totals.bar, 
                    Group == "FAC" & Year == 2022)
OBL22.bar <- filter(totals.bar, 
                    Group == "OBL" & Year == 2022)
FAC23.bar <- filter(totals.bar, 
                    Group == "FAC" & Year == 2023)
OBL23.bar <- filter(totals.bar, 
                    Group == "OBL" & Year == 2023)

# Testing for differences between TRT -------------------------------------

obl.model <- glmmTMB(Abundance ~ cTreat + (1|Paddock) + (1|Year),               # test whether average robel measurements were different in each grazing intensity
                     data = OBLtotals.bar,                                      # select the data
                     family = poisson(link="log"))                              # use family Poisson because it is count data

summary(obl.model)

diagnose(obl.model)

simulationOutput <- simulateResiduals(obl.model)
plot(simulationOutput)
testZeroInflation(simulationOutput)

# This runs the obligate with years broken out

obl.model21 <- glmmTMB(Abundance ~ cTreat + (1|Paddock),                        # test whether average robel measurements were different in each grazing intensity
                       data = OBL21.bar,                                        # select the data
                       family = poisson(link = "log"))                          # use family Poisson because it is count data

summary(obl.model21)

diagnose(obl.model21)

simulationOutput <- simulateResiduals(obl.model21)
plot(simulationOutput)
testZeroInflation(simulationOutput) 

emmeans(obl.model21,                                                            # select the model
        pairwise ~ cTreat)                                                      # compare each treatment to each other


obl.model22 <- glmmTMB(Abundance ~ cTreat + (1|Paddock),                            # test whether average robel measurements were different in each grazing intensity
                       data = OBL22.bar,                                              # select the data
                       family = poisson(link = "log"))                                  # use family Poisson because it is count data

summary(obl.model22)

diagnose(obl.model22)

simulationOutput <- simulateResiduals(obl.model22)
plot(simulationOutput)
testZeroInflation(simulationOutput) 


obl.model23 <- glmmTMB(Abundance ~ cTreat + (1|Paddock),                            # test whether average robel measurements were different in each grazing intensity
                       data = OBL23.bar,                                              # select the data
                       family = poisson(link = "log"))                                  # use family Poisson because it is count data

summary(obl.model23)

diagnose(obl.model23)

simulationOutput <- simulateResiduals(obl.model23)
plot(simulationOutput)
testZeroInflation(simulationOutput) 

# This runs the facultative models with are broken down by year

fac.model <- glmmTMB(Abundance ~ cTreat + (1|Paddock) + (1|Year),              # test whether average robel measurements were different in each grazing intensity
                     data = FACtotals.bar,                                      # select the data
                     family = poisson(link="log"))                              # use family Poisson because it is count data

summary(fac.model)

diagnose(fac.model)

simulationOutput <- simulateResiduals(fac.model)
plot(simulationOutput)
testZeroInflation(simulationOutput) 

emmeans(fac.model,                                                              # select the model
        pairwise ~ cTreat)

# This runs the facigate with years broken out

fac.model21 <- glmmTMB(Abundance ~ cTreat + (1|Paddock),                            # test whether average robel measurements were different in each grazing intensity
                       data = FAC21.bar,                                              # select the data
                       family = poisson(link = "log"))                                  # use family Poisson because it is count data

summary(fac.model21)

diagnose(fac.model21)

simulationOutput <- simulateResiduals(fac.model21)
plot(simulationOutput)
testZeroInflation(simulationOutput) 

emmeans(fac.model21,                                                            # select the model
        pairwise ~ cTreat)                                                      # compare each treatment to each other

fac.model22 <- glmmTMB(Abundance ~ cTreat + (1|Paddock),                            # test whether average robel measurements were different in each grazing intensity
                       data = FAC22.bar,                                              # select the data
                       family = poisson(link = "log"))                                  # use family Poisson because it is count data

summary(fac.model22)

diagnose(fac.model22)

simulationOutput <- simulateResiduals(fac.model22)
plot(simulationOutput)
testZeroInflation(simulationOutput) 


fac.model23 <- glmmTMB(Abundance ~ cTreat + (1|Paddock),                            # test whether average robel measurements were different in each grazing intensity
                       data = FAC23.bar,                                              # select the data
                       family = poisson(link = "log"))                                  # use family Poisson because it is count data

summary(fac.model23)

diagnose(fac.model23)

simulationOutput <- simulateResiduals(fac.model23)
plot(simulationOutput)
testZeroInflation(simulationOutput) 

emmeans(fac.model23,
        pairwise ~ cTreat)

# Creating OBG and FAC bar plots ------------------------------------------

(obl.plot21 <- ggplot(OBL21.bar,                                                  # select the data to graph
                     aes(x=cTreat,                                                # define the x axis
                         y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod4") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="white",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="white",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.y = element_text(size=20, colour = "black"),                    # color the axis text
        axis.text.x = element_blank(),
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = "2021 Average Avian Abundance", 
       x = NULL, 
       y = "Mean OBL Abundance") +                                               # changing axis titles
  coord_cartesian(xlim=NULL, ylim=c(0,10)))

(fac.plot21 <- ggplot(FAC21.bar,                                                  # select the data to graph
                      aes(x=cTreat,                                                # define the x axis
                          y=Abundance)) +
    geom_boxplot(aes(fill=Group),
                 colour= "black",                                                       # create black outlines around the bar plot
                 size=.5,
                 notch = FALSE, 
                 outlier.color = "black",
                 outlier.size = 3,
                 outlier.shape = NA) +                                                      # not sure what this does
    scale_fill_manual(values="goldenrod3") +                                   # select the color for group 2
    theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                    size=30,
                                    hjust=.5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill="white",                     # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill="white",                      # make the outer background transparent
                                         colour=NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text = element_text(size=20, colour = "black"),                    # color the axis text
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size=20,                                              # change the size of the axis titles
                            colour = "black"),                                    # change the color of the axis titles
          legend.position = "none") +                                             # remove the legend
    geom_point(aes(fill=Group),
               position = position_jitterdodge(jitter.width = 0, 
                                               dodge.width = .75)) +
    labs(title = NULL, x=NULL, y = "Mean FAC Abundance") +
    coord_cartesian(xlim=NULL, ylim=c(0,50)))

(obl.plot22 <- ggplot(OBL22.bar,                                                  # select the data to graph
                      aes(x=cTreat,                                                # define the x axis
                          y=Abundance)) +
    geom_boxplot(aes(fill=Group),
                 colour= "black",                                                       # create black outlines around the bar plot
                 size=.5,
                 notch = FALSE, 
                 outlier.color = "black",
                 outlier.size = 3,
                 outlier.shape = NA) +                                                      # not sure what this does
    scale_fill_manual(values="goldenrod4") +                                   # select the color for group 2
    theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                    size=30,
                                    hjust=.5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill="white",                           # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill="white",                            # make the outer background transparent
                                         colour=NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_blank(),                                          # color the axis text
          axis.text.x = element_blank(),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_blank(),                                                   # change the color of the axis titles
          legend.position = "none") +                                             # remove the legend
    geom_point(aes(fill=Group),
               position = position_jitterdodge(jitter.width = 0, 
                                               dodge.width = .75)) +
    labs(title = "2022 Average Avian Abundance", 
         x = NULL, 
         y = NULL) +
    coord_cartesian(xlim=NULL, ylim=c(0,10)))

(fac.plot22 <- ggplot(FAC22.bar,                                                  # select the data to graph
                     aes(x=cTreat,                                                # define the x axis
                         y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod3") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="white",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="white",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = NULL, x=NULL, y = NULL) +
  coord_cartesian(xlim=NULL, ylim=c(0,50)))

(obl.plot23 <- ggplot(OBL23.bar,                                                  # select the data to graph
                      aes(x=cTreat,                                                # define the x axis
                          y=Abundance)) +
    geom_boxplot(aes(fill=Group),
                 colour= "black",                                                       # create black outlines around the bar plot
                 size=.5,
                 notch = FALSE, 
                 outlier.color = "black",
                 outlier.size = 3,
                 outlier.shape = NA) +                                                      # not sure what this does
    scale_fill_manual(values="goldenrod4") +                                   # select the color for group 2
    theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                    size=30,
                                    hjust=.5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill="white",                           # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill="white",                            # make the outer background transparent
                                         colour=NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_blank(),                                          # color the axis text
          axis.text.x = element_blank(),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_blank(),                                                   # change the color of the axis titles
          legend.position = "none") +                                             # remove the legend
    geom_point(aes(fill=Group),
               position = position_jitterdodge(jitter.width = 0, 
                                               dodge.width = .75)) +
    labs(title = "2023 Average Avian Abundance", 
         x = NULL, 
         y = NULL) +
    coord_cartesian(xlim=NULL, ylim=c(0,10)))

(fac.plot23 <- ggplot(FAC23.bar,                                                  # select the data to graph
                     aes(x=cTreat,                                                # define the x axis
                         y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod3") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="white",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="white",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = NULL, x=NULL, y = NULL) +
  coord_cartesian(xlim=NULL, ylim=c(0,50)))


# Combing OBL and FAC bar plots -------------------------------------------


(birds21.bar <- (obl.plot21/fac.plot21))

(birds22.bar <- (obl.plot22/fac.plot22))

(birds23.bar <- (obl.plot23/fac.plot23))

(birds.bar <- (birds21.bar|birds22.bar|birds23.bar))

# Creating totals bar plots -----------------------------------------------

(obl.plot <- ggplot(OBLtotals.bar,                                                  # select the data to graph
                   aes(x=cTreat,                                                # define the x axis
                       y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod4") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="white",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="white",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = "Mean Avian Abundance", x=NULL, y = "Mean OBL Abundance") +
  coord_cartesian(xlim=NULL, ylim=c(0,20)))

(fac.plot <- ggplot(FACtotals.bar,                                                  # select the data to graph
                   aes(x=cTreat,                                                # define the x axis
                       y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod3") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="white",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="white",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = NULL, x=NULL, y = "Mean FAC Abundance") +
  coord_cartesian(xlim=NULL, ylim=c(0,100)))

(totalbirds.bar <- (obl.plot/fac.plot))

ggsave(birds.bar, 
       filename = "outputs/figs/BirdsBarYear.png",  
       dpi = "print", 
       bg = "black",
       width = 15,
       height = 10)

ggsave(totalbirds.bar, 
       filename = "outputs/figs/BirdsBar.png",  
       dpi = "print", 
       bg = "black",
       width = 15,
       height = 10)
