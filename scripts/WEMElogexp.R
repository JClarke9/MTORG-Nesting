# Load libraries ----------------------------------------------------------


library(lubridate)
library(vegan)
library(tidyverse)

# Data import -------------------------------------------------------------

raw <- read.csv("~/Git/NDSU/Nesting.csv")                                       # read in the data set

raw$Date <- as.Date(raw$Date,                                                   # read the Date column as a date
                    "%m/%d/%y") |>                                            # show the format of the existing dates
  yday()                                                                      # calculate the Date day

raw$FirstFound <- as.Date(raw$FirstFound,                                                   # read the Date column as a date
                          "%m/%d/%y") |>                                            # show the format of the existing dates
  yday()                                                                      # calculate the Date day

raw$aRobel <- rowMeans(raw[,28:31])                                             # create a new column with the robel readings averaged

raw$Litter.Depth <- as.numeric(raw$Litter.Depth)
raw$Veg.Height <- as.numeric(raw$Veg.Height)
raw$cTreat <- as.factor(raw$cTreat)
raw$pTreat <- as.factor(raw$pTreat)
raw$Stage <- as.factor(raw$Stage)
raw$Pasture <- as.factor(raw$Pasture)
raw$Patch <- as.factor(raw$Patch)
raw$Survive <- as.factor(raw$Survive)
raw$Year <- as.factor(raw$Year)
raw$id <- as.factor(raw$id)
raw$Spec <- as.factor(raw$Spec)
raw$AgeFound <- as.numeric(raw$AgeFound)

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
  filter(Spec=="WEME") |> 
  na.omit() |> 
  ungroup()

attributes(WEME)$na.action <- NULL

WEME$trials <- 1

WEME[,c(6,19:31,39)] <- scale(WEME[,c(6,19:31,39)], center = TRUE, scale = TRUE)

WEME$Stage <- recode(WEME$Stage,
                     "Incubating" = "0",
                     "Laying" = "1",
                     "Nestling" = "2")
WEME$Stage <- as.numeric(WEME$Stage)
## STEP 1a

R1 <- glmer(Fate/trials ~ (1|id),
            family = binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
R2 <- glmer(Fate/trials ~ (1|id) + (1|Year),
            family = binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
R3 <- glmer(Fate/trials ~ (1|Pasture/id),
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
R4 <- glmer(Fate/trials ~ (1|Pasture/id) + (1|Year),
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)

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

## STEP 1b

E1 <- glmer(Fate/trials ~ (1|id) + cTreat,
            family = binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
S1 <- glmer(Fate/trials ~ (1|id) + Stage,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
B1 <- glmer(Fate/trials ~ (1|id) + BHCOpres,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
T1 <- glmer(Fate/trials ~ (1|id) + Date,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
T4 <- glmer(Fate/trials ~ (1|id) + Date + Date^2,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)

## AIC Table

Cand.Exp.mods <- list(
  "id-Null" = R1,
  "id+cTreat-E1" = E1,
  "id+Stage=S1" = S1,
  "id+BHCOpres=B1" = B1,
  "id+Date=T1" = T1,
  "id+Date2=T4" = T4)

aictab(cand.set = Cand.Exp.mods,
       second.ord = TRUE)

## summary of best model to investigate beta’s and SE

R1
summary(R1)
coeffs(R1)
confint(R1,
        method = c("Wald"))

Cp(R1)
Cp(E1)

##STEP 2

T1 <- glmer(Fate/trials ~ (1|id) + Date,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
T4 <- glmer(Fate/trials ~ (1|id) + Date + Date^2,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
T5 <- glmer(Fate/trials ~ (1|id) + AgeFound,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)

## AIC Table

Cand.Temp.mods <- list(
  "id-Null" = R1, 
  "Date-T1" = T1, 
  "Stage-T2" = T2, 
  "Date+Stage-T3" = T3,
  "Quad Date-T4" = T4,
  "Age-T5" = T5)

aictab(cand.set = Cand.Temp.mods,
       second.ord = TRUE)


##### summary of best model to investigate beta’s and SE

summary(T2)
T2

coeffs(T2)
confint(T2,
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
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N2 <- glmer(Fate/trials ~ (1|id) + Stage + Smooth.Brome,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N3 <- glmer(Fate/trials ~ (1|id) + Stage + Litter,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N4 <- glmer(Fate/trials ~ (1|id) + Stage + Bare,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N5 <- glmer(Fate/trials ~ (1|id) + Stage + Forb,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N6 <- glmer(Fate/trials ~ (1|id) + Stage + Grasslike,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N7 <- glmer(Fate/trials ~ (1|id) + Stage + Woody,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N8 <- glmer(Fate/trials ~ (1|id) + Stage + Litter.Depth,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N9 <- glmer(Fate/trials ~ (1|id) + Stage + Veg.Height,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N10<- glmer(Fate/trials ~ (1|id) + Stage + aRobel,
            family=binomial(logexp(exposure=WEME$Expos)),
            data=WEME)
N11 <- glmer(Fate/trials ~ (1|id) + Stage + Litter.Depth + aRobel,
             family=binomial(logexp(exposure=WEME$Expos)),
             data=WEME)

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
# Log Exposure Models -----------------------------------------------------


