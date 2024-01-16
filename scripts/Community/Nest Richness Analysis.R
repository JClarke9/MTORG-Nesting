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


birds21 <- read.csv("working/birds21.csv")                             # read in the data set

birds22 <- read.csv("working/birds22.csv")                             # read in the data set

birds23 <- read.csv("working/birds23.csv")                             # read in the data set

windowsFonts(my_font = windowsFont("Gandhi Sans"))


# Breaking down site by species matrices ----------------------------------


birds21$cTreat <- factor(birds21$cTreat,
                         levels = c("Rest",
                                    "Moderate",
                                    "Full",
                                    "Heavy"))

birds22$cTreat <- factor(birds22$cTreat,
                         levels = c("Rest",
                                    "Moderate",
                                    "Full",
                                    "Heavy"))

birds23$cTreat <- factor(birds23$cTreat,
                         levels = c("Rest",
                                    "Moderate",
                                    "Full",
                                    "Heavy"))

birds21.obl <- select(birds21, c(cTreat, STGR, UPSA, SAVS, GRSP, CCLO, DICK, BOBO, WEME))
birds21.fac <- select(birds21, c(cTreat, GADW, AMWI, MALL, BWTE, NSHO, NOPI, WILL, MODO, 
                                 EAKI, CCSP, RWBL, BRBL, KILL,))
birds21.wet <- select(birds21, c(cTreat, WIPH, WISN, YHBL))

birds22.obl <- select(birds22, c(cTreat, STGR, UPSA, SAVS, GRSP, CCLO, DICK, BOBO, WEME))
birds22.fac <- select(birds22, c(cTreat, GADW, MALL, BWTE, NSHO, NOPI, GWTE, MODO, 
                                 EAKI, CCSP, RWBL, BRBL, KILL, COGR))
birds22.wet <- select(birds22, c(cTreat, WIPH, YHBL, PBGR, SORA, AMCO, MAWR))

birds23.obl <- select(birds23, c(cTreat, STGR, UPSA, GRSP, CCLO, BOBO, WEME))
birds23.fac <- select(birds23, c(cTreat, GADW, AMWI, MALL, BWTE, NSHO, NOPI, GWTE, WILL, MODO, 
                                 EAKI, CCSP, RWBL, BRBL, KILL, AMBI, LESC))
birds23.wet <- select(birds23, c(cTreat, WISN, YHBL, PBGR, SORA, AMCO, MAWR))


# Total species richness calculations -------------------------------------


rich21.tot <- specnumber(birds21[2:25],
                         birds21$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich21.tot) <- c("treat", 
                          "rich21.tot")

rich22.tot <- specnumber(birds22[2:28],
                         birds22$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich22.tot) <- c("treat", 
                          "rich22.tot")

rich23.tot <- specnumber(birds23[2:29],
                         birds23$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich23.tot) <- c("treat", 
                          "rich23.tot")

rich21.tot
rich22.tot
rich23.tot


# Obligate species richness -----------------------------------------------


rich21.obl <- specnumber(birds21.obl[2:8],
                         birds21.obl$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich21.obl) <- c("treat", 
                          "rich21.obl")

rich22.obl <- specnumber(birds22.obl[2:8],
                         birds22.obl$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich22.obl) <- c("treat", 
                          "rich22.obl")

rich23.obl <- specnumber(birds23.obl[2:7],
                         birds23.obl$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich23.obl) <- c("treat", 
                          "rich23.obl")

rich21.obl
rich22.obl
rich23.obl


# Facultative species richness --------------------------------------------


rich21.fac <- specnumber(birds21.fac[2:14],
                         birds21.fac$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich21.fac) <- c("treat", 
                          "rich21.fac")

rich22.fac <- specnumber(birds22.fac[2:14],
                         birds22.fac$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich22.fac) <- c("treat", 
                          "rich22.fac")

rich23.fac <- specnumber(birds23.fac[2:17],
                         birds23.fac$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich23.fac) <- c("treat", 
                          "rich23.fac")

rich21.fac
rich22.fac
rich23.fac


# Wetland species richness ------------------------------------------------


rich21.wet <- specnumber(birds21.wet[2:4],
                         birds21.wet$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich21.wet) <- c("treat", 
                          "rich21.wet")

rich22.wet <- specnumber(birds22.wet[2:7],
                         birds22.wet$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich22.wet) <- c("treat", 
                          "rich22.wet")

rich23.wet <- specnumber(birds23.wet[2:7],
                         birds23.wet$cTreat,
                         MARGIN = 1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(rich23.wet) <- c("treat", 
                          "rich23.wet")

rich21.wet
rich22.wet
rich23.wet


# Creating Species Richness Dataframe -------------------------------------


rich21 <- full_join(rich21.tot, 
                    rich21.obl, 
                    by = "treat") |> 
  full_join(rich21.fac, 
            by = "treat") |> 
  full_join(rich21.wet, 
            by = "treat")

rich22 <- full_join(rich22.tot, 
                    rich22.obl, 
                    by = "treat") |> 
  full_join(rich22.fac, 
            by = "treat") |> 
  full_join(rich22.wet, 
            by = "treat")

rich23 <- full_join(rich23.tot, 
                    rich23.obl, 
                    by = "treat") |> 
  full_join(rich23.fac, 
            by = "treat") |> 
  full_join(rich23.wet, 
            by = "treat")

rich <- full_join(rich21.tot, 
                  rich22.tot,
                  by = "treat") |> 
  full_join(rich23.tot,
            by = "treat") |> 
  full_join(rich21.obl, 
            by = "treat") |> 
  full_join(rich22.obl, 
            by = "treat") |> 
  full_join(rich23.obl,
            by = "treat") |> 
  full_join(rich21.fac, 
            by = "treat") |> 
  full_join(rich22.fac, 
            by = "treat") |> 
  full_join(rich23.fac,
            by = "treat") |> 
  full_join(rich21.wet, 
            by = "treat") |> 
  full_join(rich22.wet, 
            by = "treat") |> 
  full_join(rich23.wet,
            by = "treat")

rich <- rich |> 
  rename("Grazing Intensity" = "treat",
         "Total Richness 2021" = "rich21.tot",
         "Total Richness 2022" = "rich22.tot",
         "Total Richness 2023" = "rich23.tot",
         "OBL Richness 2021" = "rich21.obl",
         "OBL Richness 2022" = "rich22.obl",
         "OBL Richness 2023" = "rich23.obl",
         "FAC Richness 2021" = "rich21.fac",
         "FAC Richness 2022" = "rich22.fac",
         "FAC Richness 2023" = "rich23.fac",
         "WET Richness 2021" = "rich21.wet",
         "WET Richness 2022" = "rich22.wet",
         "WET Richness 2023" = "rich23.wet")

write.csv(rich, file = "working/richness.csv")


# Manipulating data -------------------------------------------------------


padd21.obl <- birds21.obl |> 
  rownames_to_column(var = "Paddock")

padd22.obl <- birds22.obl |> 
  rownames_to_column(var = "Paddock")

padd23.obl <- birds23.obl |> 
  rownames_to_column(var = "Paddock")

padd21.obl$Paddock <- factor(padd21.obl$Paddock,
                             levels = c("1", "2", "3", "4",
                                        "5", "6", "7", "8",
                                        "9", "10", "11", "12",
                                        "13", "14", "15", "16"))

padd22.obl$Paddock <- factor(padd22.obl$Paddock,
                             levels = c("1", "2", "3", "4",
                                        "5", "6", "7", "8",
                                        "9", "10", "11", "12",
                                        "13", "14", "15", "16"))

padd23.obl$Paddock <- factor(padd23.obl$Paddock,
                             levels = c("1", "2", "3", "4",
                                        "5", "6", "7", "8",
                                        "9", "10", "11", "12",
                                        "13", "14", "15", "16"))


# Calculating obligate species richness for each paddock ------------------


rich21.oblp <- data.frame("Richness" = specnumber(padd21.obl[3:10],
                                                  padd21.obl$Paddock,
                                                  MARGIN = 1),
                          Treat = padd21.obl$cTreat,
                          Year = 2021,
                          Paddock = padd21.obl$Paddock)

rich22.oblp <- data.frame("Richness" = specnumber(padd22.obl[3:10],
                                                  padd22.obl$Paddock,
                                                  MARGIN = 1),
                          Treat = padd22.obl$cTreat,
                          Year = 2022,
                          Paddock = padd22.obl$Paddock)

rich23.oblp <- data.frame("Richness" = specnumber(padd23.obl[3:8],
                                                  padd23.obl$Paddock,
                                                  MARGIN = 1),
                          Treat = padd23.obl$cTreat,
                          Year = 2023,
                          Paddock = padd23.obl$Paddock)


# Testing obligate richness -----------------------------------------------


aov21.obl <- glmmTMB(Richness ~ Treat, 
                     data = rich21.oblp,
                     family = poisson(link = "log"))

diagnose(aov21.obl)
# remove random effect

summary(aov21.obl)


simulationOutput <- simulateResiduals(aov21.obl)
plot(simulationOutput)
testZeroInflation(simulationOutput)

aov22.obl <- glmmTMB(Richness ~ Treat, 
                     data = rich22.oblp,
                     family = poisson(link = "log"))

diagnose(aov22.obl)
# remove random effect

summary(aov22.obl)

simulationOutput <- simulateResiduals(aov22.obl)
plot(simulationOutput)
testZeroInflation(simulationOutput)


aov23.obl <- glmmTMB(Richness ~ Treat, 
                     data = rich23.oblp,
                     family = poisson(link = "log"))

diagnose(aov23.obl)
# remove random effect

summary(aov23.obl)

simulationOutput <- simulateResiduals(aov23.obl)
plot(simulationOutput)
testZeroInflation(simulationOutput)


# Manipulating data -------------------------------------------------------


padd21.fac <- birds21.fac |> 
  rownames_to_column(var = "Paddock")

padd22.fac <- birds22.fac |> 
  rownames_to_column(var = "Paddock")

padd23.fac <- birds23.fac |> 
  rownames_to_column(var = "Paddock")

padd21.fac$Paddock <- factor(padd21.fac$Paddock,
                             levels = c("1", "2", "3", "4",
                                        "5", "6", "7", "8",
                                        "9", "10", "11", "12",
                                        "13", "14", "15", "16"))

padd22.fac$Paddock <- factor(padd22.fac$Paddock,
                             levels = c("1", "2", "3", "4",
                                        "5", "6", "7", "8",
                                        "9", "10", "11", "12",
                                        "13", "14", "15", "16"))

padd23.fac$Paddock <- factor(padd23.fac$Paddock,
                             levels = c("1", "2", "3", "4",
                                        "5", "6", "7", "8",
                                        "9", "10", "11", "12",
                                        "13", "14", "15", "16"))


# Calculating facultative species richness for each paddock ---------------


rich21.facp <- data.frame("Richness" = specnumber(padd21.fac[3:15],
                                                  padd21.fac$Paddock,
                                                  MARGIN = 1),
                          Treat = padd21.fac$cTreat,
                          Year = 2021,
                          Paddock = padd21.fac$Paddock)

rich22.facp <- data.frame("Richness" = specnumber(padd22.fac[3:15],
                                                  padd22.fac$Paddock,
                                                  MARGIN = 1),
                          Treat = padd22.fac$cTreat,
                          Year = 2022,
                          Paddock = padd22.fac$Paddock)

rich23.facp <- data.frame("Richness" = specnumber(padd23.fac[3:18],
                                                  padd23.fac$Paddock,
                                                  MARGIN = 1),
                          Treat = padd23.fac$cTreat,
                          Year = 2023,
                          Paddock = padd23.fac$Paddock)


# Testing facultative richness --------------------------------------------


aov21.fac <- glmmTMB(Richness ~ Treat + (1|Paddock), 
                     data = rich21.facp,
                     family = poisson(link = "log"))

diagnose(aov21.fac)
# remove random effect

summary(aov21.fac)

simulationOutput <- simulateResiduals(aov21.fac)
plot(simulationOutput)
testZeroInflation(simulationOutput)

# remove random effect

aov22.fac <- glmmTMB(Richness ~ Treat + (1|Paddock), 
                     data = rich22.facp,
                     family = poisson(link = "log"))
summary(aov22.fac)

simulationOutput <- simulateResiduals(aov22.fac)
plot(simulationOutput)
testZeroInflation(simulationOutput)

diagnose(aov22.fac)
# remove random effect

aov23.fac <- glmmTMB(Richness ~ Treat + (1|Paddock), 
                     data = rich23.facp,
                     family = poisson(link = "log"))
summary(aov23.fac)

simulationOutput <- simulateResiduals(aov23.fac)
plot(simulationOutput)
testZeroInflation(simulationOutput)

diagnose(aov23.fac)
# remove random effect


# Combining the data frames -----------------------------------------------


rich.oblp <- cbind(rich21.oblp,
                   rich22.oblp,
                   rich23.oblp) |> 
  ungroup()

rich.facp <- rbind(rich21.facp,
                   rich22.facp,
                   rich23.facp) |> 
  ungroup()

aov.obl <- glmmTMB(Richness ~ Treat,
                   data = rich.oblp,
                   family = poisson(link = "log"))

# remove random effects
diagnose(aov.obl)

summary(aov.obl)

simulationOutput <- simulateResiduals(aov.obl)
plot(simulationOutput)
testZeroInflation(simulationOutput)

aov.fac <- glmmTMB(Richness ~ Treat,
                   data = rich.facp,
                   family = poisson(link = "log"))

# remove random effects
diagnose(aov.fac)

summary(aov.fac)

simulationOutput <- simulateResiduals(aov.fac)
plot(simulationOutput)
testZeroInflation(simulationOutput)


# Manipulating data -------------------------------------------------------

padd21 <- full_join(padd21.obl, padd21.fac)
padd22 <- full_join(padd22.obl, padd22.fac)
padd23 <- full_join(padd23.obl, padd23.fac)

# Calculating the total richness in each paddock --------------------------


rich21.padd <- data.frame("Richness" = specnumber(padd21[3:23],
                                                  padd21$Paddock,
                                                  MARGIN = 1),
                          Treat = padd21$cTreat,
                          Year = 2021,
                          Paddock = padd21$Paddock)

rich22.padd <- data.frame("Richness" = specnumber(padd22[3:23],
                                                  padd22$Paddock,
                                                  MARGIN = 1),
                          Treat = padd22$cTreat,
                          Year = 2022,
                          Paddock = padd22$Paddock)

rich23.padd <- data.frame("Richness" = specnumber(padd23[3:24],
                                                 padd23$Paddock,
                                                 MARGIN = 1),
                         Treat = padd23$cTreat,
                         Year = 2023,
                         Paddock = padd23$Paddock)

rich.tot <- rbind(rich21.padd,
                  rich22.padd,
                  rich23.padd) |> 
  ungroup()

rich.tot$Year <- factor(rich.tot$Year,
                        levels = c("2021",
                                   "2022",
                                   "2023"))


# Testing total richness --------------------------------------------------


aov.tot <- glmmTMB(Richness ~ Treat, 
                   data = rich.tot,
                   family = poisson(link = "log"))

diagnose(aov.tot)

summary(aov.tot)

simulationOutput <- simulateResiduals(aov.tot)
plot(simulationOutput)
testZeroInflation(simulationOutput)

write.csv(rich.tot, "working/TotalRichness.csv")
