# Loading libraries -------------------------------------------------------

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

# I wanted to get the simpson diversity for each of the
# grazing intensities for 2021
simpson21.group <- diversity(species21[2:22],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="simpson",                                   # using the simpson index
                             MARGIN = 1)                                        # not positive what this does

summary(simpson21.group)                                                        # show the output of the simpson diversity index
simpson21.group <- as.data.frame(simpson21.group)                               # assign the simpson diversity index to a data frame
simpson21.group$Treat <- species21$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
simpson21.group                                                                 # show the output for the simpson diversity

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

aov.simp21 <- glmmTMB(simpson21~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      simpson21,
                      family = Gamma(link = "log"))

# remove random effect
diagnose(aov.simp21)

summary(aov.simp21)

simulationOutput <- simulateResiduals(aov.simp21)
plot(simulationOutput)
testZeroInflation(simulationOutput)

# Calculate species diversity ---------------------------------------------

species22 <- aggregate(birds22[3:23],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds22$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the simpson diversity for each of the
# grazing intensities for 2022
simpson22.group <- diversity(species22[2:22],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="simpson",                                           # using the simpson index
                             MARGIN = 1)                                                # not positive what this does

summary(simpson22.group)                                                        # show the output of the simpson diversity index
simpson22.group <- as.data.frame(simpson22.group)                               # assign the simpson diversity index to a data frame
simpson22.group$Treat <- species22$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
simpson22.group                                                                 # show the output for the simpson diversity

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

aov.simp22 <- glmmTMB(simpson22~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      data = simpson22,
                      family = Gamma(link = "log"))

# remove random effect
diagnose(aov.simp22)

summary(aov.simp22)

simulationOutput <- simulateResiduals(aov.simp22)
plot(simulationOutput)
testZeroInflation(simulationOutput)

# Calculate species diversity ---------------------------------------------

species23 <- aggregate(birds23[3:24],                                           # select the data frame for the analysis (I excluded the row with treatments) 
                       by=list(birds23$cTreat),                               # select values to group by
                       FUN=sum)                                               # add together abundances for each species

# I wanted to get the simpson diversity for each of the
# grazing intensities for 2023
simpson23.group <- diversity(species23[2:23],                                   # select the data frame for the analysis (I excluded the row with treatments) 
                             index="simpson",                                           # using the simpson index
                             MARGIN = 1)                                                # not positive what this does

summary(simpson23.group)                                                        # show the output of the simpson diversity index
simpson23.group <- as.data.frame(simpson23.group)                               # assign the simpson diversity index to a data frame
simpson23.group$Treat <- species23$Group.1                                      # create a treatment column in the new data frame to use for an ANOVA
simpson23.group                                                                 # show the output for the simpson diversity

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

aov.simp23 <- glmmTMB(simpson23~Treat + (1|Paddock),                                              # select the variables to compare with an ANOVA
                      data = simpson23,
                      family = Gamma(link = "log"))

# remove random effect
diagnose(aov.simp23)

summary(aov.simp23)

simulationOutput <- simulateResiduals(aov.simp23)
plot(simulationOutput)
testZeroInflation(simulationOutput)

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



# Diversity analysis by year ----------------------------------------------


aov.simp21 <- glmmTMB(simpson ~ Treat + (1|Paddock), 
                      data = simp21,
                      family = gaussian(link = "log"))

# remove random effect
diagnose(aov.simp21)

simulationOutput <- simulateResiduals(aov.simp21, plot = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

summary(aov.simp21)

emmeans(aov.simp21,
        pairwise~Treat,
        adjust = "bonferroni")


aov.simp22 <- glmmTMB(simpson ~ Treat + (1|Paddock), 
                      data = simp22,
                      family = Gamma(link = "log"))

# remove random effect
diagnose(aov.simp22)

simulationOutput <- simulateResiduals(aov.simp22, plot = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

summary(aov.simp22)

emmeans(aov.simp22,
        pairwise~Treat,
        adjust = "tukey")


aov.simp23 <- glmmTMB(simpson ~ Treat + (1|Paddock), 
                      data = simp23,
                      family = gaussian(link = "identity"))

# remove random effect
diagnose(aov.simp23)

simulationOutput <- simulateResiduals(aov.simp23, plot = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

summary(aov.simp23)

emmeans(aov.simp23,
        pairwise~Treat,
        adjust = "tukey")


# Total diversity analysis ------------------------------------------------


aov.simp <- glmmTMB(simpson ~ Treat + Year + (1|Paddock), 
                    data = simp.total)

# remove random effect
diagnose(aov.simp)

simulationOutput <- simulateResiduals(aov.simp, plot = T)
testZeroInflation(simulationOutput)
testDispersion(simulationOutput)

summary(aov.simp)

emmeans(aov.simp,
        pairwise~Treat,
        type = "Tukey")

emmeans(aov.simp,
        pairwise~Year,
        type = "Tukey")

simp.div <- full_join(simpson21.group, 
                      simpson22.group,
                      by = "Treat") |> 
  full_join(simpson23.group,
            by = "Treat") |> 
  relocate("Treat", 
           .before="simpson21.group") |> 
  rename("Grazing Intensity" = "Treat", 
         "2021 Diversity" = "simpson21.group", 
         "2022 Diversity" = "simpson22.group",
         "2023 Diversity" = "simpson23.group")

write_csv(simp21, "working/SimpsonTotal21.csv")
write_csv(simp22, "working/SimpsonTotal22.csv")
write_csv(simp23, "working/SimpsonTotal23.csv")
write_csv(simp.total, "working/SimpsonTotal.csv")
