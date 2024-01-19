# Loading libraries -------------------------------------------------------

library(vegan)
library(tidyverse)
library(wesanderson)
library(glmmTMB)
library(DHARMa)
library(emmeans)

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
                       by=list(birds21$cTreat),                                 # select values to group by
                       FUN=sum)                                                 # add together abundances for each species

# I wanted to get the shannon diversity for each of the
# grazing intensities for 2021
shannon21.group <- diversity(species21[2:22],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="shannon",                                   # using the shannon index
                             MARGIN = 1)                                        # not positive what this does

summary(shannon21.group)                                                        # show the output of the shannon diversity index
shannon21.group <- as.data.frame(shannon21.group)                               # assign the shannon diversity index to a data frame
shannon21.group$Treat <- species21$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
shannon21.group                                                                 # show the output for the shannon diversity

# I recalculated shannon diversity for each paddock to
# determine whether diversity was significantly different
# between grazing intensities

shannon21 <- diversity(birds21[3:23],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       index="shannon",                                           # using the shannon index
                       MARGIN = 1)                                                # not positive what this does

summary(shannon21)                                                              # show the output of the shannon diversity index
shannon21 <- as.data.frame(shannon21)                                           # assign the shannon diversity index to a data frame
shannon21$Treat <- birds21$cTreat                                               # create a treatment column in the new data frame to use for an ANOVA
shannon21$Paddock <- birds21$Paddock
shannon21                                                                       # show the output for the shannon diversity

shannon21$Treat <- factor(shannon21$Treat,
                          levels = c("Rest",
                                     "Moderate",
                                     "Full",
                                     "Heavy"))

shannon21$Paddock <- factor(shannon21$Paddock,
                            levels = c("1", "2", "3", "4",
                                       "5", "6", "7", "8",
                                       "9", "10", "11", "12",
                                       "13", "14", "15", "16"))

aov.shan21 <- glmmTMB(shannon21~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      shannon21,
                      family = Gamma(link = "log"))

# remove random effect
diagnose(aov.shan21)

summary(aov.shan21)

simulationOutput <- simulateResiduals(aov.shan21)
plot(simulationOutput)
testZeroInflation(simulationOutput)

# Calculate species diversity ---------------------------------------------

species22 <- aggregate(birds22[3:23],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds22$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the shannon diversity for each of the
# grazing intensities for 2022
shannon22.group <- diversity(species22[2:22],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="shannon",                                           # using the shannon index
                             MARGIN = 1)                                                # not positive what this does

summary(shannon22.group)                                                        # show the output of the shannon diversity index
shannon22.group <- as.data.frame(shannon22.group)                               # assign the shannon diversity index to a data frame
shannon22.group$Treat <- species22$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
shannon22.group                                                                 # show the output for the shannon diversity

# I recalculated shannon diversity for each paddock to
# determine whether diversity was significantly different
# between grazing intensities

shannon22 <- diversity(birds22[3:23],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       index="shannon",                                           # using the shannon index
                       MARGIN = 1)                                                # not positive what this does

summary(shannon22)                                                              # show the output of the shannon diversity index
shannon22 <- as.data.frame(shannon22)                                           # assign the shannon diversity index to a data frame
shannon22$Treat <- birds22$cTreat                                               # create a treatment column in the new data frame to use for an ANOVA
shannon22$Paddock <- birds22$Paddock
shannon22                                                                       # show the output for the shannon diversity

shannon22$Treat <- factor(shannon22$Treat,
                          levels = c("Rest",
                                     "Moderate",
                                     "Full",
                                     "Heavy"))

shannon22$Paddock <- factor(shannon22$Paddock,
                            levels = c("1", "2", "3", "4",
                                       "5", "6", "7", "8",
                                       "9", "10", "11", "12",
                                       "13", "14", "15", "16"))

aov.shan22 <- glmmTMB(shannon22~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      data = shannon22,
                      family = Gamma(link = "log"))

# remove random effect
diagnose(aov.shan22)

summary(aov.shan22)

simulationOutput <- simulateResiduals(aov.shan22)
plot(simulationOutput)
testZeroInflation(simulationOutput)

# Calculate species diversity ---------------------------------------------

species23 <- aggregate(birds23[3:24],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds23$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the shannon diversity for each of the
# grazing intensities for 2023
shannon23.group <- diversity(species23[2:23],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="shannon",                                           # using the shannon index
                             MARGIN = 1)                                                # not positive what this does

summary(shannon23.group)                                                        # show the output of the shannon diversity index
shannon23.group <- as.data.frame(shannon23.group)                               # assign the shannon diversity index to a data frame
shannon23.group$Treat <- species23$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
shannon23.group                                                                 # show the output for the shannon diversity

# I recalculated shannon diversity for each paddock to
# determine whether diversity was significantly different
# between grazing intensities

shannon23 <- diversity(birds23[3:24],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       index="shannon",                                           # using the shannon index
                       MARGIN = 1)                                                # not positive what this does

summary(shannon23)                                                              # show the output of the shannon diversity index
shannon23 <- as.data.frame(shannon23)                                           # assign the shannon diversity index to a data frame
shannon23$Treat <- birds23$cTreat                                               # create a treatment column in the new data frame to use for an ANOVA
shannon23$Paddock <- birds23$Paddock
shannon23                                                                       # show the output for the shannon diversity

shannon23$Treat <- factor(shannon23$Treat,
                          levels = c("Rest",
                                     "Moderate",
                                     "Full",
                                     "Heavy"))

shannon23$Paddock <- factor(shannon23$Paddock,
                            levels = c("1", "2", "3", "4",
                                       "5", "6", "7", "8",
                                       "9", "10", "11", "12",
                                       "13", "14", "15", "16"))

aov.shan23 <- glmmTMB(shannon23~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      data = shannon23,
                      family = Gamma(link = "log"))

# remove random effect
diagnose(aov.shan23)

summary(aov.shan23)

simulationOutput <- simulateResiduals(aov.shan23)
plot(simulationOutput)
testZeroInflation(simulationOutput)

# Creating diversity data frame -------------------------------------------

shan21 <- rename(shannon21, 
                 "shannon"="shannon21")
shan21$Year <- 2021

shan22 <- rename(shannon22, 
                 "shannon"="shannon22")
shan22$Year <- 2022

shan23 <- rename(shannon23, 
                 "shannon"="shannon23")
shan23$Year <- 2023

shan.total <- rbind(shan21, 
                    shan22, 
                    shan23)

shan.total$Year <- factor(shan.total$Year,
                          levels = c("2021",
                                     "2022",
                                     "2023"))

str(shan.total)

qqnorm(shan.total$shannon)
qqline(shan.total$shannon)
hist(shan.total$shannon)



# Diversity analysis by year ----------------------------------------------


aov.shan21 <- glmmTMB(shannon ~ Treat, 
                      data = shan21)

# remove random effect
diagnose(aov.shan21)

simulationOutput <- simulateResiduals(aov.shan21, plot = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

summary(aov.shan21)

emmeans(aov.shan21,
        pairwise~Treat,
        adjust = "tukey")


aov.shan22 <- glmmTMB(shannon ~ Treat + (1|Paddock), 
                      data = shan22,
                      family = gaussian(link = "identity"))

# remove random effect
diagnose(aov.shan22)

simulationOutput <- simulateResiduals(aov.shan22, plot = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

summary(aov.shan22)

emmeans(aov.shan22,
        pairwise~Treat,
        adjust = "tukey")


aov.shan23 <- glmmTMB(shannon ~ Treat + (1|Paddock), 
                      data = shan23)

# remove random effect
diagnose(aov.shan23)

simulationOutput <- simulateResiduals(aov.shan23, plot = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

summary(aov.shan23)

emmeans(aov.shan23,
        pairwise~Treat,
        adjust = "tukey")


# Total diversity analysis ------------------------------------------------


aov.shan <- glmmTMB(shannon ~ Treat + Year + (1|Paddock), 
                    data = shan.total)

# remove random effect
diagnose(aov.shan)

simulationOutput <- simulateResiduals(aov.shan, plot = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

summary(aov.shan)

emmeans(aov.shan,
        pairwise~Treat,
        type = "Tukey")

emmeans(aov.shan,
        pairwise~Year,
        type = "Tukey")

shan.div <- full_join(shannon21.group, 
                      shannon22.group,
                      by = "Treat") |> 
  full_join(shannon23.group,
            by = "Treat") |> 
  relocate("Treat", 
           .before="shannon21.group") |> 
  rename("Grazing Intensity" = "Treat", 
         "2021 Diversity" = "shannon21.group", 
         "2022 Diversity" = "shannon22.group",
         "2023 Diversity" = "shannon23.group")

write_csv(shan21, "working/shannonTotal21.csv")
write_csv(shan22, "working/shannonTotal22.csv")
write_csv(shan23, "working/shannonTotal23.csv")
write_csv(shan.total, "working/shannonTotal.csv")

