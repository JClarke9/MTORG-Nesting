# Load libraries ----------------------------------------------------------


library(lubridate)
library(vegan)
library(tidyverse)

# Data import -------------------------------------------------------------

raw <- read.csv("working/LogExpnesting.csv")                                       # read in the data set

raw <- raw |> 
  dplyr::select(-Open) |> 
  mutate(Year = factor(Year, levels = c("2021", "2022", "2023", "2024")),
         id = as.factor(id),
         Spec = as.factor(Spec),
         Fate = as.factor(Fate),
         Fate2 = as.factor(Fate2),
         Stage = factor(Stage, levels = c("Laying", "Incubating", "Nestling")),
         Paddock = factor(Paddock, levels = c(1:16)),
         Replicate = as.factor(Replicate),
         Treatment = as.factor(Treatment),
         cTreat = factor(cTreat, levels = c("0", "39", "49", "68")),
         pTreat = factor(pTreat, levels = c("0", "39", "49", "68")))


##################################################################################
##red-winged blackbird nest survival code from Shew et al. 2018 (Journal of Applied Ecology) 
#################################################################################

## Logistic exposure link function from (https://rpubs.com/bbolker/logregexp)

library(MASS)
logexp <- function(exposure = 1)
{
  linkfun <- function(mu) qlogis(mu^(1/exposure))
  ## FIXME: is there some trick we can play here to allow
  ##   evaluation in the context of the 'data' argument?
  linkinv <- function(eta)  plogis(eta)^exposure
  logit_mu_eta <- function(eta) {
    ifelse(abs(eta)>30,.Machine$double.eps,
           exp(eta)/(1+exp(eta))^2)
    ## OR .Call(stats:::C_logit_mu_eta, eta, PACKAGE = "stats")
  }
  mu.eta <- function(eta) {       
    exposure * plogis(eta)^(exposure-1) *
      logit_mu_eta(eta)
  }
  valideta <- function(eta) TRUE
  link <- paste("logexp(", deparse(substitute(exposure)), ")",
                sep="")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}

## Import data from proper folder using R-studio and CSV file
##install and load “lme4”, “MuMIn”, and “AICcmodavg” packages before proceeding

library(lme4)
library(MuMIn)
library(AICcmodavg)

WEME <- raw |> 
  filter(Spec == "WEME" & Stage != "Laying") |> 
  na.omit() |> 
  ungroup()

attributes(WEME)$na.action <- NULL

WEME$Stage <- recode(WEME$Stage,
                     "Incubating" = "0",
                     "Nestling" = "1")
WEME$Stage <- as.numeric(WEME$Stage)



## STEP 1
R1 <- glmer(Fate ~ (1|id),
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
R2 <- glmer(Fate ~ (1|id) + (1|Year),
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
R3 <- glmer(Fate ~ (1|Replicate/id),
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
R4 <- glmer(Fate ~ (1|Replicate/id) + (1|Year),
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)

## AIC Table

Cand.Rand.mods <- list(
  "id-R1" = R1,
  "id+Year-R2" = R2,
  "Pasture/id-R3" = R3, 
  "Pasture/id+Year-R4" = R4)

aictab(cand.set = Cand.Rand.mods,
       second.ord = TRUE)

## summary of best model to investigate beta’s and SE

R1
summary(R1)
coeffs(R1)
confint(R1,
        method = c("Wald"))

Cp(R1)
Cp(R2)
Cp(R3)
Cp(R4)



## STEP 2
T1 <- glmer(Fate ~ (1|id) + DateChecked,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
T2 <- glmer(Fate ~ (1|id) + Year,
            family=binomial(logexp(exposure = WEME$Expos)),
            data = WEME)

## AIC Table

Cand.Exp.mods <- list(
  "Null" = R1,
  "Linear" = T1,
  "Year" = T2)

aictab(cand.set = Cand.Exp.mods,
       second.ord = TRUE)

## summary of best model to investigate beta’s and SE

R1
summary(R1)
coeffs(R1)
confint(R1,
        method = c("Wald"))

Cp(R1)



##STEP 3
S1 <- glmer(Fate ~ (1|id) + Stage,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)

## AIC Table

Cand.Temp.mods <- list(
  "Null" = R1,
  "stage" = S1)

aictab(cand.set = Cand.Temp.mods,
       second.ord = TRUE)

## summary of best model to investigate beta’s and SE

summary(S1)
S1

coeffs(S1)
confint(S1,
        method = c("Wald"))

Cp(S1)



##STEP 4
B1 <- glmer(Fate ~ (1|id) + Stage + BHCONum,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)

## AIC Table

Cand.Temp.mods <- list(
  "stage" = S1,
  "BHCO" = B1)

aictab(cand.set = Cand.Temp.mods,
       second.ord = TRUE)


##### summary of best model to investigate beta’s and SE

summary(S1)
S1

coeffs(S1)
confint(S1,
        method = c("Wald"))

Cp(S1)



##STEP 5
G1 <- glmer(Fate ~ (1|id) + Stage + grazed,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
G2 <- glmer(Fate ~ (1|id) + Stage + pDoD,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)

## AIC Table

Cand.Temp.mods <- list(
  "stage" = S1,
  "grazed" = G1,
  "pDoD" = G2)

aictab(cand.set = Cand.Temp.mods,
       second.ord = TRUE)


##### summary of best model to investigate beta’s and SE

summary(S1)
S1

coeffs(S1)
confint(S1,
        method = c("Wald"))

Cp(S1)
## R-squared calculations
r.squaredLR(S1,R1)

##Marginal and Conditional R-squared calculations from Nakagawa, S. & Schielzeth H. (2013)

##fitted model = T2
##random effect name = “id”
## Calculation of the variance in fitted values
VarF <- var(as.vector(fixef(S1) %*% t(model.matrix(S1))))

## R2GLMM(m) - marginal R2GLMM
VarF/(VarF + VarCorr(S1)$id[1] + pi^2/3)

## R2GLMM(c) - conditional R2GLMM for full model
(VarF + VarCorr(S1)$id[1])/(VarF + VarCorr(S1)$id[1] + pi^2/3)

## Proportion Change in Variance (PCV) for fitted model (rand.ef1)
100*(1 - (VarCorr(S1)$id[1])/(VarCorr(R1)$id[1]))



##STEP 6

V1 <- glmer(Fate ~ (1|id) + Stage + KBG,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
V2 <- glmer(Fate ~ (1|id) + Stage + SmoothB,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
V3 <- glmer(Fate ~ (1|id) + Stage + Litter,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
V4 <- glmer(Fate ~ (1|id) + Stage + Bare,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
V5 <- glmer(Fate ~ (1|id) + Stage + Forb,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
V6 <- glmer(Fate ~ (1|id) + Stage + Grasslike,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
V7 <- glmer(Fate ~ (1|id) + Stage + Woody,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
V8 <- glmer(Fate ~ (1|id) + Stage + LitterD,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)
V9 <- glmer(Fate ~ (1|id) + Stage + VOR,
            family = binomial(logexp(exposure = WEME$Expos)),
            data = WEME)

##AIC Table
Cand.Nest.mods <- list(
  "stage" = S1,
  "KBG" = V1,
  "SmoothB" = V2,
  "litter" = V3,
  "bare" = V4,
  "forb" = V5,
  "grass" = V6,
  "woody" = V7,
  "litdep" = V8,
  "vor" = V9)

aictab(cand.set = Cand.Nest.mods, second.ord = TRUE)


### summary of best model to investigate beta’s and SE

summary(S1)
S1

coeffs(S1)
confint(S1 ,method = c("Wald"))

Cp(S1)


## R-squared calculations
r.squaredLR(S1,R1)
