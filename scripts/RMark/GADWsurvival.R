# Loading libraries -------------------------------------------------------

library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)
library(MuMIn)
source("scripts/Functions/RMark_Stage_Code.R")

windowsFonts(my_font = windowsFont("Gandhi Sans"))

# Data import -------------------------------------------------------------

nest <- read.csv("working/RMarknesting.csv", 
                 row.names=1)

# Subsetting data ---------------------------------------------------------

GADW.surv <- filter(nest, 
                    Spec=="GADW")                                         # select out only GADW nest

test <- filter(GADW.surv,
               is.na(KBG) |
                 is.na(SmoothB) |
                 is.na(Litter) |
                 is.na(Bare) |
                 is.na(Forb) |
                 is.na(Grasslike) |
                 is.na(Woody) |
                 is.na(LitterD) |
                 is.na(Veg.Height) |
                 is.na(VOR) |
                 is.na(cTreat))

MISSING <- is.na(GADW.surv$AgeFound)

sum(MISSING)

GADW.surv <- subset(GADW.surv, 
                    subset = !MISSING)

GADW.surv$Year <- factor(GADW.surv$Year,
                         levels = c("2021", "2022", "2023"))

str(GADW.surv)

# Creating stage variable -------------------------------------------------

x <- create.stage.var(GADW.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(GADW.surv$LastChecked)), 
                      12)

GADW.surv <- bind_cols(GADW.surv, x)

rm(list = ls()[!ls() %in% c("GADW.surv")])

# Daily survival rate models ----------------------------------------------

GADW.pr <- process.data(GADW.surv,
                        nocc=max(GADW.surv$LastChecked),
                        groups = c("Year"),
                        model="Nest")

# Temporal candidate model set
GADW1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  GADW.model.list = create.model.list("Nest")
  GADW1.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

# Results of candidate model set
GADW1.results <- GADW1.run()
GADW1.results

coef(GADW1.results$S.year)
confint(GADW1.results$S.year, level = 0.85)


# Biological candidate model set
GADW2.run <- function()
{
  # 1. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Year + NestAge)
  
  GADW.model.list = create.model.list("Nest")
  GADW2.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
GADW2.results <- GADW2.run()
GADW2.results

coef(GADW2.results$S.age)

confint(GADW2.results$S.age, level = 0.85)


# Grazing candidate model set
GADW3.run <- function()
{
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Year + NestAge)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + NestAge + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Year + NestAge + grazep)
  
  # 4. DSR varies with the previous Year + NestAges grazing intensity
  S.pTreat = list(formula = ~1 + Year + NestAge + pTreat)
  
  # 4. DSR varies with the previous Year + NestAges grazing intensity
  S.grazedpTreat = list(formula = ~1 + Year + NestAge + grazed + pTreat)
  
  GADW.model.list = create.model.list("Nest")
  GADW3.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
GADW3.results <- GADW3.run()
GADW3.results


# Vegetation candidate model set
GADW4.run <- function()
{
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Year + NestAge)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Year + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Year + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Year + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Year + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Year + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Year + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Year + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + VOR)
  
  GADW.model.list = create.model.list("Nest")
  GADW4.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
GADW4.results <- GADW4.run()
GADW4.results

coef(GADW4.results$S.forb)
confint(GADW4.results$S.forb, level = 0.85)


# Plotting Beta Coefficients ----------------------------------------------


GADW.mod <- mark(GADW.surv, 
                 nocc=max(GADW.surv$LastChecked), 
                 model = "Nest", 
                 model.parameters = list(S = list(formula =  ~1 + LitterD)))


# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all=TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
#  files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)
