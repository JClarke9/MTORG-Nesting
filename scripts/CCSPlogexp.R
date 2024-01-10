# Load libraries ----------------------------------------------------------


library(lubridate)
library(vegan)
library(tidyverse)

# Data import -------------------------------------------------------------

raw <- read.csv("working/Community.csv",
                row.names = 1)                                       # read in the data set

raw$Litter.Depth <- as.numeric(raw$Litter.Depth)
raw$Veg.Height <- as.numeric(raw$Veg.Height)
raw$AgeFound <- as.numeric(raw$AgeFound)

raw$id <- factor(raw$id)
raw$Spec <- factor(raw$Spec)

raw$cTreat <- factor(raw$cTreat,
                     c("Rest", "Moderate", "Full", "Heavy"))
raw$pTreat <- factor(raw$pTreat,
                     c("Rest", "Moderate", "Full", "Heavy"))
raw$Stage <- factor(raw$Stage,
                    c("Laying", "Incubating", "Nestling"))
raw$Paddock <- factor(raw$Paddock,
                      c(1:16))
raw$Replicate <- factor(raw$Replicate,
                        c("NE", "SE", "SW", "NW"))
raw$Year <- factor(raw$Year,
                   levels = c("2021", "2022", "2023"))

str(raw)

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

CCSP <- raw |> 
  filter(Spec=="CCSP" & Fate2 != "Unknown") |> 
  ungroup()

MISSING <- is.na(CCSP$KBG)
MISSING <- is.na(CCSP$Smooth.Brome)
MISSING <- is.na(CCSP$Litter)
MISSING <- is.na(CCSP$Bare)
MISSING <- is.na(CCSP$Forb)
MISSING <- is.na(CCSP$Grasslike)
MISSING <- is.na(CCSP$Woody)
MISSING <- is.na(CCSP$Litter.Depth)
MISSING <- is.na(CCSP$Veg.Height)
MISSING <- is.na(CCSP$VOR)

sum(MISSING)

CCSP <- subset(CCSP, subset = !MISSING)

CCSP <- CCSP |> 
  filter(CCSP$Expos != 0)

CCSP$trials <- 1

## STEP 1a

R1 <-glmer(Fate/trials ~ (1|id), 
           family = binomial(logexp(exposure = CCSP$Expos)),
           data = CCSP)

R2 <-glmer(Fate/trials ~ ( 1|Paddock/id), 
           family = binomial(logexp(exposure = CCSP$Expos)) ,
           data = CCSP)

R3<-glmer(Fate/trials ~ (1|Paddock/id) + (1|Year), 
          family = binomial(logexp(exposure = CCSP$Expos)),
          data = CCSP)

R4<-glmer(Fate/trials ~ (1|Replicate/Paddock/id), 
          family = binomial(logexp(exposure = CCSP$Expos)),
          data = CCSP)

R5 <-glmer(Fate/trials ~  (1|Replicate/Paddock/id) +(1|Year), 
           family = binomial(logexp(exposure = CCSP$Expos)),
           data = CCSP)

## AIC Table

Cand.Rand.mods <- list(
  "id-R1" = R1,
  "Paddock/id-R2" = R2, 
  "Paddock/id+Year-R3" = R3, 
  "Replicate/Paddock/id-R4" = R4,
  "Replicate/Paddock/id+Year-R5" = R5)

aictab(cand.set = Cand.Rand.mods,
       second.ord = TRUE)

## summary of best model to investigate beta’s and SE

R3
summary(R3)
coeffs(R3)
confint(R3,
        method = c("Wald"))

R5
summary(R5)
coeffs(R5)
confint(R5,
        method = c("Wald"))

Cp(R1)
Cp(R2)
Cp(R3)
Cp(R4)

##STEP 2

T1 <- glmer(Fate/trials ~ (1|Paddock/id) + (1|Year) + Date,
            family = binomial(logexp(exposure = CCSP$Expos)),
            data = CCSP)

T2 <- glmer(Fate/trials ~ (1|Paddock/id) + (1|Year) + Date + I(Date^2),
            family = binomial(logexp(exposure = CCSP$Expos)),
            data = CCSP)

T3 <- glmer(Fate/trials ~ (1|Paddock/id) + (1|Year) + Stage,
            family = binomial(logexp(exposure = CCSP$Expos)),
            data = CCSP)

## AIC Table

Cand.Temp.mods <- list(
  "id-Null-R3" = R3, 
  "Date-T1" = T1, 
  "Quad Date-T2" = T2,
  "Stage-T3" = T3)

aictab(cand.set = Cand.Temp.mods,
       second.ord = TRUE)


##### summary of best model to investigate beta’s and SE

summary(T3)
T3

coeffs(T3)
confint(T3,
        method = c("Wald"))

Cp(R1)
Cp(T1)
Cp(T2)
Cp(T3)
Cp(T4)
Cp(T5)

## R-squared calculations
r.squaredLR(R1,R1)
r.squaredLR(T1,R1)
r.squaredLR(T2,R1)
r.squaredLR(T3,R1)
r.squaredLR(T4,R1)
r.squaredLR(T5,R1)

##Marginal and Conditional R-squared calculations from Nakagawa, S. & Schielzeth H. (2013)

##fitted model = T2
##random effect name = “id”
## Calculation of the variance in fitted values
VarF <- var(as.vector(fixef(T2) %*% t(model.matrix(T2))))

## R2GLMM(m) - marginal R2GLMM
VarF/(VarF + VarCorr(T2)$id[1] + pi^2/3)

## R2GLMM(c) - conditional R2GLMM for full model
(VarF + VarCorr(T2)$id[1])/(VarF + VarCorr(T2)$id[1] + pi^2/3)

## Proportion Change in Variance (PCV) for fitted model (rand.ef1)
100*(1 - (VarCorr(T2)$id[1])/(VarCorr(R1)$id[1]))

##STEP 3a

N1 <- glmer(Fate/trials ~ (1|id) + Stage + KBG,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N2 <- glmer(Fate/trials ~ (1|id) + Stage + Smooth.Brome,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N3 <- glmer(Fate/trials ~ (1|id) + Stage + Litter,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N4 <- glmer(Fate/trials ~ (1|id) + Stage + Bare,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N5 <- glmer(Fate/trials ~ (1|id) + Stage + Forb,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N6 <- glmer(Fate/trials ~ (1|id) + Stage + Grasslike,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N7 <- glmer(Fate/trials ~ (1|id) + Stage + Woody,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N8 <- glmer(Fate/trials ~ (1|id) + Stage + Litter.Depth,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N9 <- glmer(Fate/trials ~ (1|id) + Stage + Veg.Height,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N10<- glmer(Fate/trials ~ (1|id) + Stage + aRobel,
            family=binomial(logexp(exposure=CCSP$Expos)),
            data=CCSP)
N11 <- glmer(Fate/trials ~ (1|id) + Stage + Litter.Depth + aRobel,
             family=binomial(logexp(exposure=CCSP$Expos)),
             data=CCSP)

##AIC Table
Cand.Nest.mods <- list(
  "id–Null" = R1,
  "Stage–T2" = T2,
  "KBG-N1" = N1,
  "Smooth Brome-N2" = N2,
  "Litter-N3" = N3,
  "Bare-N4" = N4,
  "Forb-N5" = N5,
  "Grasslike-N6" = N6,
  "Woody-N7" = N7,
  "Litter Depth-N8" = N8,
  "Veg Height-N9" = N9,
  "aRobel-N10" = N10,
  "Litter Depth+aRobel-N11" = N11)

aictab(cand.set = Cand.Nest.mods, second.ord = TRUE)


### summary of best model to investigate beta’s and SE

summary(N4)
N4

coeffs(N4)
confint(N4 ,method = c("Wald"))

Cp(R2)
Cp(T3)
Cp(N1)
Cp(N2)
Cp(N3)
Cp(N4)
Cp(N5)
Cp(N6)
Cp(N7)
Cp(N8)
Cp(N9)
Cp(N10)
Cp(N11)


## R-squared calculations
r.squaredLR(R2,R1)
r.squaredLR(T3,R1)
r.squaredLR(N1,R1)
r.squaredLR(N2,R1)
r.squaredLR(N3,R1)
r.squaredLR(N4,R1)
r.squaredLR(N5,R1)
r.squaredLR(N6,R1)
r.squaredLR(N7,R1)
r.squaredLR(N8,R1)
r.squaredLR(N9,R1)
r.squaredLR(N10,R1)
r.squaredLR(N11,R1)

##fitted model = T3
##randome effect name = “Pasture”, “id”
## Calculation of the variance in fitted values
VarF <- var(as.vector(fixef(N3) %*% t(model.matrix(N3))))


## R2GLMM(m) - marginal R2GLMM
VarF/(VarF + VarCorr(N3)$id[1] + VarCorr(N3)$Pasture[1] +  pi^2/3)

## R2GLMM(c) - conditional R2GLMM for full model
(VarF + VarCorr(N3)$id[1] + VarCorr(N3)$Pasture[1])/(VarF + VarCorr(N3)$id[1] + VarCorr(N3)$Pasture[1] + pi^2/3)

## Proportion Change in Variance (PCV) for fitted model (rand.ef1)
100*(1 - (VarCorr(N3)$id[1])/(VarCorr(R2)$id[1]))

## Proportion Change in Variance (PCV) for fitted model (rand.ef2)
100*(1 - (VarCorr(N3)$Pasture[1])/(VarCorr(R2)$Pasture[1]))