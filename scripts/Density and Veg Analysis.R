# Load libraries ----------------------------------------------------------


library(tidyverse)


# Data import -------------------------------------------------------------


raw <- read.csv("working/VegAdj.csv")

density <- read.csv("working/Birds_Treatment_Density.csv")


# Data wrangling ----------------------------------------------------------


raw$Date <- as.Date(raw$Date,                                                   # read the Date column as a date
                    "%m/%d/%y")    

MTORG.veg <- filter(raw,                                               # select the data frame
                    Treatment=="MTORG")

MTORG.veg$a.Robel <- MTORG.veg |> 
  select(R1:R4) |> 
  rowMeans()

MTORG.veg$Transect.ID <- factor(MTORG.veg$Transect.ID, levels = c(97:144))

MTORG.veg$Pasture <- factor(MTORG.veg$Pasture, levels = c("NE", "SE",
                                                          "SW", "NW"))

MTORG.veg$SubPatch <- factor(MTORG.veg$SubPatch, levels = c("Rest", "Moderate", 
                                                            "Full", "Heavy"))

MTORG.veg$Patch <- factor(MTORG.veg$Patch, levels = c(1:16))

MTORG.veg$TotalVegCover <- MTORG.veg |> 
  select(KBG:Woody) |> 
  rowSums(na.rm=TRUE)

MTORG.veg <- MTORG.veg |> 
  mutate(across(KBG:Woody, ~ .x/TotalVegCover * 100))

MTORG.str <- MTORG.veg |> 
  group_by(Year,
           Pasture,
           SubPatch) |> 
  summarize(Litter.Depth=mean(Litter.Depth),
            a.Robel = mean(a.Robel)) |> 
  ungroup()

MTORG.str <- rename(MTORG.str,
                    "Replicate" = "Pasture",
                    "cTreat" = "SubPatch")


# Cleaning up density dataframe ------------------------------------------------------------------------------

MTORG.density <- data.frame()

for (Spec in unique(density$Species)) {
  temp <- filter(density, Species == Spec)
  temp <- complete(temp,
                   Year, Replicate, cTreat)
  temp <- full_join(temp, MTORG.str,
                    by = c("Replicate", "Year", "cTreat"))
  MTORG.density <- bind_rows(temp, MTORG.density)
  rm(list = "temp")
}

MTORG.density$Year <- factor(MTORG.density$Year,
                             levels = c("2021", "2022", "2023"))

MTORG.density$cTreat <- factor(MTORG.density$cTreat,
                               levels = c("Rest", "Moderate", "Full", "Heavy"))

MTORG.density$Replicate <- factor(MTORG.density$Replicate,
                                  levels = c("NE", "SE", "SW", "NW"))


# Testing litter impacts on density --------------------------------------------------------------------------



library(glmmTMB)
library(emmeans)
library(DHARMa)


BRBL.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "BRBL"),
                    ziformula = ~Litter.Depth,
                    family = ziGamma(link = "log"))

BRBLoutput <- simulateResiduals(BRBL.lit,
                                n = 999,
                                plot = T)

diagnose(BRBL.lit)
testResiduals(BRBLoutput)
testZeroInflation(BRBLoutput)
testQuantiles(BRBLoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(BRBL.lit)


WEME.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "WEME"),
                    ziformula = ~Litter.Depth,
                    family = ziGamma(link = "log"))

WEMEoutput <- simulateResiduals(WEME.lit,
                                n = 999,
                                plot = T)

diagnose(WEME.lit)
testResiduals(WEMEoutput)
testZeroInflation(WEMEoutput)
testQuantiles(WEMEoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(WEME.lit)


CCSP.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "CCSP"),
                    ziformula = ~Litter.Depth,
                    family = ziGamma(link = "log"))

CCSPoutput <- simulateResiduals(CCSP.lit,
                                n = 999,
                                plot = T)

diagnose(CCSP.lit)
testResiduals(CCSPoutput)
testZeroInflation(CCSPoutput)
testQuantiles(CCSPoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(CCSP.lit)

# The conditional model was significant for CCSP

MODO.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "MODO"),
                    ziformula = ~Litter.Depth,
                    family = ziGamma(link = "log"))

MODOoutput <- simulateResiduals(MODO.lit,
                                n = 999,
                                plot = T)

diagnose(MODO.lit)
testResiduals(MODOoutput)
testZeroInflation(MODOoutput)
testQuantiles(MODOoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(MODO.lit)


RWBL.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "RWBL"),
                    ziformula = ~Litter.Depth,
                    family = ziGamma(link = "log"))

RWBLoutput <- simulateResiduals(RWBL.lit,
                                n = 999,
                                plot = T)

diagnose(RWBL.lit)
testResiduals(RWBLoutput)
testZeroInflation(RWBLoutput)
testQuantiles(RWBLoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(RWBL.lit)


NOPI.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "NOPI"),
                    ziformula = ~Litter.Depth,
                    family = ziGamma(link = "log"))

NOPIoutput <- simulateResiduals(NOPI.lit,
                                n = 999,
                                plot = T)

diagnose(NOPI.lit)
testResiduals(NOPIoutput)
testZeroInflation(NOPIoutput)
testQuantiles(NOPIoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(NOPI.lit)


GADW.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "GADW"),
                    ziformula = ~Litter.Depth,
                    family = ziGamma(link = "log"))

GADWoutput <- simulateResiduals(GADW.lit,
                                n = 999,
                                plot = T)

diagnose(GADW.lit)
testResiduals(GADWoutput)
testZeroInflation(GADWoutput)
testQuantiles(GADWoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(GADW.lit)

# the zero-inflated model was significant for GADW

BWTE.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "BWTE"),
                    ziformula = ~Litter.Depth,
                    family = ziGamma(link = "log"))

BWTEoutput <- simulateResiduals(BWTE.lit,
                                n = 999,
                                plot = T)

diagnose(BWTE.lit)
testResiduals(BWTEoutput)
testZeroInflation(BWTEoutput)
testQuantiles(BWTEoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(BWTE.lit)


# Testing VOR impacts on density -----------------------------------------------------------------------------


BRBL.vor <- glmmTMB(estimate~a.Robel + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "BRBL"),
                    ziformula = ~.,
                    family = ziGamma(link = "log"))

BRBLoutput <- simulateResiduals(BRBL.vor,
                                n = 999,
                                plot = T)

diagnose(BRBL.vor)
testResiduals(BRBLoutput)
testZeroInflation(BRBLoutput)
testQuantiles(BRBLoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(BRBL.vor)


WEME.vor <- glmmTMB(estimate~a.Robel + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "WEME"),
                    ziformula = ~.,
                    family = ziGamma(link = "log"))

WEMEoutput <- simulateResiduals(WEME.vor,
                                n = 999,
                                plot = T)

diagnose(WEME.vor)
testResiduals(WEMEoutput)
testZeroInflation(WEMEoutput)
testQuantiles(WEMEoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(WEME.vor)

# This model does not seem to work
CCSP.vor <- glmmTMB(estimate~a.Robel + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "CCSP"),
                    ziformula = ~a.Robel,
                    family = ziGamma(link = "log"),
                    control = glmmTMBControl())

CCSPoutput <- simulateResiduals(CCSP.vor,
                                n = 999,
                                plot = T)

diagnose(CCSP.vor)
testResiduals(CCSPoutput)
testZeroInflation(CCSPoutput)
testQuantiles(CCSPoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(CCSP.vor)

# The conditional model was significant for CCSP

MODO.vor <- glmmTMB(estimate~a.Robel + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "MODO"),
                    ziformula = ~.,
                    family = ziGamma(link = "log"))

MODOoutput <- simulateResiduals(MODO.vor,
                                n = 999,
                                plot = T)

diagnose(MODO.vor)
testResiduals(MODOoutput)
testZeroInflation(MODOoutput)
testQuantiles(MODOoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(MODO.vor)


RWBL.vor <- glmmTMB(estimate~a.Robel + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "RWBL"),
                    ziformula = ~a.Robel,
                    family = ziGamma(link = "log"))

RWBLoutput <- simulateResiduals(RWBL.vor,
                                n = 999,
                                plot = T)

diagnose(RWBL.vor)
testResiduals(RWBLoutput)
testZeroInflation(RWBLoutput)
testQuantiles(RWBLoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(RWBL.vor)


NOPI.vor <- glmmTMB(estimate~a.Robel + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "NOPI"),
                    ziformula = ~a.Robel,
                    family = ziGamma(link = "log"))

NOPIoutput <- simulateResiduals(NOPI.vor,
                                n = 999,
                                plot = T)

diagnose(NOPI.vor)
testResiduals(NOPIoutput)
testZeroInflation(NOPIoutput)
testQuantiles(NOPIoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(NOPI.vor)


GADW.vor <- glmmTMB(estimate~a.Robel + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "GADW"),
                    ziformula = ~a.Robel,
                    family = ziGamma(link = "log"))

GADWoutput <- simulateResiduals(GADW.vor,
                                n = 999,
                                plot = T)

diagnose(GADW.vor)
testResiduals(GADWoutput)
testZeroInflation(GADWoutput)
testQuantiles(GADWoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(GADW.vor)

# the zero-inflated model was significant for GADW

BWTE.vor <- glmmTMB(estimate~a.Robel + (1|Year/Replicate),
                    data = filter(MTORG.density, Species == "BWTE"),
                    ziformula = ~a.Robel,
                    family = ziGamma(link = "log"))

BWTEoutput <- simulateResiduals(BWTE.vor,
                                n = 999,
                                plot = T)

diagnose(BWTE.vor)
testResiduals(BWTEoutput)
testZeroInflation(BWTEoutput)
testQuantiles(BWTEoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(BWTE.vor)
