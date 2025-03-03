# Load libraries ----------------------------------------------------------


library(tidyverse)
library(glmmTMB)


# Data import -------------------------------------------------------------


raw <- read.csv("working/MTORG_str.csv") |> 
  mutate(Replicate = factor(Replicate, levels = c("NE", "SE",
                                                  "SW", "NW")),
         Intensity = factor(Intensity, levels = c("Rest", "Moderate", 
                                                  "Full", "Heavy")),
         Year = factor(Year, levels = c("2021", "2022", 
                                        "2023", "2024")),
         Patch = factor(Patch, levels = c(1:16)))

density <- read.csv("working/Birds_Treatment_Density.csv") |> 
  mutate(Replicate = factor(Replicate, levels = c("NE", "SE",
                                                  "SW", "NW")),
         cTreat = factor(cTreat, levels = c("Rest", "Moderate", 
                                            "Full", "Heavy")),
         Year = factor(Year, levels = c("2021", "2022", 
                                        "2023", "2024")))


# Cleaning up density dataframe ------------------------------------------------------------------------------


MTORG.density <- data.frame()

for (Spec in unique(density$Species)) {
  temp <- filter(density, Species == Spec)
  temp <- complete(temp,
                   Year, Replicate, cTreat)
  temp <- full_join(temp, raw,
                    by = c("Replicate", "Year", "cTreat" = "Intensity"))
  MTORG.density <- bind_rows(temp, MTORG.density)
  rm(list = "temp")
}

MTORG.density <- MTORG.density |> 
  mutate(Year = factor(Year,
                       levels = c("2021", "2022", 
                                  "2023", "2024")),
         cTreat = factor(cTreat,
                         levels = c("Rest", "Moderate", 
                                    "Full", "Heavy")),
         Replicate = factor(Replicate,
                            levels = c("NE", "SE", 
                                       "SW", "NW")))


# Testing litter impacts on density --------------------------------------------------------------------------



library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)


BRBL.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "BRBL"),
                    family = t_family(link = "identity"))

BRBLoutput <- simulateResiduals(BRBL.lit,
                                n = 999,
                                plot = T)

diagnose(BRBL.lit)
testResiduals(BRBLoutput)
testZeroInflation(BRBLoutput)
testQuantiles(BRBLoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(BRBL.lit, type = 2)


WEME.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "WEME"),
                    family = t_family(link = "identity"))

WEMEoutput <- simulateResiduals(WEME.lit,
                                n = 999,
                                plot = T)

diagnose(WEME.lit)
testResiduals(WEMEoutput)
testZeroInflation(WEMEoutput)
testQuantiles(WEMEoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(WEME.lit, type = 2)


CCSP.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "CCSP"),
                    family = t_family(link = "identity"))

CCSPoutput <- simulateResiduals(CCSP.lit,
                                n = 999,
                                plot = T)

diagnose(CCSP.lit)
testResiduals(CCSPoutput)
testZeroInflation(CCSPoutput)
testQuantiles(CCSPoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(CCSP.lit, type = 2)

# The conditional model was significant for CCSP

MODO.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Patch),
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

Anova(MODO.lit, type = 2)


RWBL.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "RWBL"),
                    family = t_family(link = "identity"))

RWBLoutput <- simulateResiduals(RWBL.lit,
                                n = 999,
                                plot = T)

diagnose(RWBL.lit)
testResiduals(RWBLoutput)
testZeroInflation(RWBLoutput)
testQuantiles(RWBLoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(RWBL.lit, type = 2)


NOPI.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "NOPI"),
                    family = t_family(link = "identity"))

NOPIoutput <- simulateResiduals(NOPI.lit,
                                n = 999,
                                plot = T)

diagnose(NOPI.lit)
testResiduals(NOPIoutput)
testZeroInflation(NOPIoutput)
testQuantiles(NOPIoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(NOPI.lit, type = 2)


GADW.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "GADW"),
                    family = t_family(link = "identity"))

GADWoutput <- simulateResiduals(GADW.lit,
                                n = 999,
                                plot = T)

diagnose(GADW.lit)
testResiduals(GADWoutput)
testZeroInflation(GADWoutput)
testQuantiles(GADWoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(GADW.lit, type = 2)

# the zero-inflated model was significant for GADW
# 
# BWTE.lit <- glmmTMB(estimate~Litter.Depth + (1|Year/Patch),
#                     data = filter(MTORG.density, Species == "BWTE"),
#                     ziformula = ~Litter.Depth,
#                     family = ziGamma(link = "log"))

BWTEoutput <- simulateResiduals(BWTE.lit,
                                n = 999,
                                plot = T)

diagnose(BWTE.lit)
testResiduals(BWTEoutput)
testZeroInflation(BWTEoutput)
testQuantiles(BWTEoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(BWTE.lit)


# Testing VOR impacts on density -----------------------------------------------------------------------------


BRBL.vor <- glmmTMB(estimate~a.Robel + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "BRBL"),
                    family = t_family(link = "identity"))

BRBLoutput <- simulateResiduals(BRBL.vor,
                                n = 999,
                                plot = T)

diagnose(BRBL.vor)
testResiduals(BRBLoutput)
testZeroInflation(BRBLoutput)
testQuantiles(BRBLoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(BRBL.vor, type = 2)


WEME.vor <- glmmTMB(estimate~a.Robel + (1|Year/Patch),
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

Anova(WEME.vor, type = 2)

# This model does not seem to work
CCSP.vor <- glmmTMB(estimate~a.Robel + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "CCSP"),
                    family = t_family(link = "identity"))

CCSPoutput <- simulateResiduals(CCSP.vor,
                                n = 999,
                                plot = T)

diagnose(CCSP.vor)
testResiduals(CCSPoutput)
testZeroInflation(CCSPoutput)
testQuantiles(CCSPoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(CCSP.vor, type = 2)

# The conditional model was significant for CCSP

# MODO.vor <- glmmTMB(estimate~a.Robel + (1|Year/Patch),
#                     data = filter(MTORG.density, Species == "MODO"),
#                     ziformula = ~.,
#                     family = ziGamma(link = "identity"))

MODOoutput <- simulateResiduals(MODO.vor,
                                n = 999,
                                plot = T)

diagnose(MODO.vor)
testResiduals(MODOoutput)
testZeroInflation(MODOoutput)
testQuantiles(MODOoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(MODO.vor, type = 2)


RWBL.vor <- glmmTMB(estimate~a.Robel + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "RWBL"),
                    family = t_family(link = "identity"))

RWBLoutput <- simulateResiduals(RWBL.vor,
                                n = 999,
                                plot = T)

diagnose(RWBL.vor)
testResiduals(RWBLoutput)
testZeroInflation(RWBLoutput)
testQuantiles(RWBLoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(RWBL.vor, type = 2)


NOPI.vor <- glmmTMB(estimate~a.Robel + (1|Year/Patch),
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

Anova(NOPI.vor, type = 2)


GADW.vor <- glmmTMB(estimate~a.Robel + (1|Year/Patch),
                    data = filter(MTORG.density, Species == "GADW"),
                    ziformula = ~.,
                    family = ziGamma(link = "log"))

GADWoutput <- simulateResiduals(GADW.vor,
                                n = 999,
                                plot = T)

diagnose(GADW.vor)
testResiduals(GADWoutput)
testZeroInflation(GADWoutput)
testQuantiles(GADWoutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

Anova(GADW.vor, type = 2)

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

Anova(BWTE.vor, type = 2)