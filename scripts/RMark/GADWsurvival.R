# Loading libraries -------------------------------------------------------

library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)

windowsFonts(my_font = windowsFont("Gandhi Sans"))

# Data import -------------------------------------------------------------

nest <- read.csv("working/RMarknesting.csv", 
                 row.names=1)

# Subsetting data ---------------------------------------------------------

GADW.surv <- filter(nest, 
                    Spec=="GADW")                                         # select out only GADW nests

GADW.surv$AgeDay1 <- GADW.surv$AgeFound - GADW.surv$FirstFound + 1
GADW.surv$Year <- as.factor(GADW.surv$Year)
GADW.surv$cTreat <- as.factor(GADW.surv$cTreat)

# Creating stage variable -------------------------------------------------

create.stage.var=function(data,agevar.name,stagevar.name,time.intervals,cutoff)
{
  nocc=length(time.intervals)
  age.mat=matrix(data[,agevar.name],nrow=dim(data)[1],ncol=nocc-1)
  age.mat=t(t(age.mat)+cumsum(c(0,time.intervals[1:(nocc-2)])))
  stage.mat=t(apply(age.mat,1,function(x) as.numeric(x<=cutoff)))
  stage.mat=data.frame(stage.mat)
  names(stage.mat)=paste(stagevar.name,1:(nocc-1),sep="")
  return(stage.mat)
}

x <- create.stage.var(GADW.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(GADW.surv$LastChecked)), 
                      23)

GADW.surv <- bind_cols(GADW.surv, x)

rm(list = ls()[!ls() %in% c("GADW.surv")])

# Daily survival rate models ----------------------------------------------

GADW.pr <- process.data(GADW.surv,
                        nocc=max(GADW.surv$LastChecked), 
                        groups = c("cTreat",
                                   "Year"),
                        model="Nest")

# candidate models for the null vs. treatment model
GADW1.run <- function()
{
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  # 3. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies With year interacting with treatment (treatment identity varies by year)
  S.treatyear = list(formula = ~1 + cTreat:Year)
  
  GADW.model.list = create.model.list("Nest")
  GADW1.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
GADW1.results <- GADW1.run()
GADW1.results

GADW1.results$S.Dot$results$beta
GADW1.results$S.year$results$beta

#candidate model set for time trends
GADW2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  GADW.model.list = create.model.list("Nest")
  GADW2.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
GADW2.results <- GADW2.run()
GADW2.results

GADW2.results$S.Dot$results$beta
GADW2.results$S.quad$results$beta
GADW2.results$S.julian$results$beta
GADW2.results$S.year$results$beta

#candidate model set for grazing
GADW3.run <- function()
{
  # 1. DSR varies with BHCO number
  S.grazed = list(formula = ~1 + grazed)
  
  S.grazep = list(formula = ~1 + grazep)
  
  # 2. DSR varies with quadratic effect of date
  S.Dot = list(formula = ~1)
  GADW.model.list = create.model.list("Nest")
  GADW3.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
GADW3.results <- GADW3.run()
GADW3.results

GADW3.results$S.Dot$results$beta
GADW3.results$S.grazed$results$beta
GADW3.results$S.grazep$results$beta

#candidate model set for parasitism and nest stage
GADW4.run <- function()
{
  # 1. DSR varies with KBG
  S.kbg = list(formula =  ~1 + KBG)
  
  # 2. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.smooth.brome = list(formula = ~1 + SmoothB)
  
  # 3. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Litter)
  
  # 4. DSR varies with Bare
  S.bare = list(formula =  ~1 + Bare)
  
  # 5. DSR varies with Forb
  S.forb = list(formula =  ~1 + Forb)
  
  # 6. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Grasslike)
  
  # 7. DSR varies with Woody
  S.woody = list(formula =  ~1 + Woody)
  
  # 8. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + LitterD)
  
  # 9. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Veg.Height)
  
  # 10. DSR varies with VOR
  S.vor = list(formula =  ~1 + VOR)
  
  # 11. DSR varies with litter and Veg Height
  S.litheight = list(formula = ~1 + Litter + Veg.Height)
  
  # 12. DSR varies with woody and litter depth
  S.woodylitdep = list(formula = ~1 + Woody + LitterD)
  
  # 13. DSR varies with KBG and litter depth
  S.kbglitdep = list(formula = ~1 + KBG + LitterD)
  
  # 14. DSR varies with KBG and Veg.Height
  S.kbgheight = list(formula = ~1 + KBG + Veg.Height)
  
  # 15. DSR varies with woody and Veg Height
  S.woodyheight = list(formula = ~1 + Woody + Veg.Height)
  
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  GADW.model.list = create.model.list("Nest")
  GADW4.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

#results of candidate model set 
#all results
GADW4.results <- GADW4.run()
GADW4.results

GADW4.results$S.litdep$results$beta
GADW4.results$S.bare$results$beta
GADW4.results$S.forb$results$beta
GADW4.results$S.height$results$beta
GADW4.results$S.kbglitdep$results$beta
GADW4.results$S.woodylitdep$results$beta


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
