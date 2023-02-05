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

NOPI.surv <- filter(nest, 
                    Spec=="NOPI")                                         # select out only NOPI nests

NOPI.surv$AgeDay1 <- NOPI.surv$AgeFound - NOPI.surv$FirstFound + 1
NOPI.surv$Year <- as.factor(NOPI.surv$Year)
NOPI.surv$cTreat <- as.factor(NOPI.surv$cTreat)

NOPI.surv$AgeDay1 <- NOPI.surv$AgeFound - NOPI.surv$FirstFound + 1
NOPI.surv$Year <- as.factor(NOPI.surv$Year)
NOPI.surv$cTreat <- as.factor(NOPI.surv$cTreat)

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

x <- create.stage.var(NOPI.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(NOPI.surv$LastChecked)), 
                      23)

NOPI.surv <- bind_cols(NOPI.surv, x)

rm(list = ls()[!ls() %in% "NOPI.surv"])

# Daily survival rate models ----------------------------------------------

NOPI.pr <- process.data(NOPI.surv,
                        nocc=max(NOPI.surv$LastChecked),
                        groups = c("cTreat", 
                                   "Year"),
                        model="Nest")

# candidate models for the null vs. treatment model
NOPI1.run <- function()
{
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  # 3. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies With year interacting with treatment (treatment identity varies by year)
  S.treatyear = list(formula = ~1 + cTreat:Year)
  
  NOPI.model.list = create.model.list("Nest")
  NOPI1.results = mark.wrapper(NOPI.model.list,
                               data = NOPI.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
NOPI1.results <- NOPI1.run()
NOPI1.results

NOPI1.results$S.Dot$results$beta
NOPI1.results$S.year$results$beta

#candidate model set for time trends
NOPI2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  NOPI.model.list = create.model.list("Nest")
  NOPI2.results = mark.wrapper(NOPI.model.list,
                               data = NOPI.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
NOPI2.results <- NOPI2.run()
NOPI2.results

NOPI2.results$S.Dot$results$beta
NOPI2.results$S.time$results$beta
NOPI2.results$S.quad$results$beta

#candidate model set for grazing
NOPI3.run <- function()
{
  # 1. DSR varies with BHCO number
  S.grazed = list(formula = ~1 + grazed)
  
  S.grazep = list(formula = ~1 + grazep)
  
  # 2. DSR varies with quadratic effect of date
  S.Dot = list(formula = ~1)
  
  NOPI.model.list = create.model.list("Nest")
  NOPI3.results = mark.wrapper(NOPI.model.list,
                               data = NOPI.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
NOPI3.results <- NOPI3.run()
NOPI3.results

#candidate model set for parasitism and nest stage
NOPI4.run <- function()
{
  # 1. DSR varies with KBG
  S.kbg = list(formula = ~1 + KBG)
  
  # 2. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.smooth.brome = list(formula = ~1 + SmoothB)
  
  # 3. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula = ~1 + Litter)
  
  # 4. DSR varies with Bare
  S.bare = list(formula = ~1 + Bare)
  
  # 5. DSR varies with Forb
  S.forb = list(formula = ~1 + Forb)
  
  # 6. DSR varies with Grasslike (correlated with KBG)
  S.grass = list(formula = ~1 + Grasslike)
  
  # 7. DSR varies with Woody
  S.woody = list(formula = ~1 + Woody)
  
  # 8. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula = ~1 + LitterD)
  
  # 9. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula = ~1 + Veg.Height)
  
  # 10. DSR varies with VOR
  S.vor = list(formula = ~1 + VOR)
  
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
  
  NOPI.model.list = create.model.list("Nest")
  NOPI4.results = mark.wrapper(NOPI.model.list,
                               data = NOPI.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

#results of candidate model set 
#all results
NOPI4.results <- NOPI4.run()
NOPI4.results

NOPI4.results$S.lit$results$beta
NOPI4.results$S.litheight$results$beta

# Plotting beta coefficients ----------------------------------------------

NOPI.mod <- mark(NOPI.surv, 
                 nocc=max(NOPI.surv$LastChecked), 
                 model = "Nest", 
                 model.parameters = list(S = list(formula =  ~1 + Litter)))

NOPI4.beta <- NOPI4.results$S.lit$results$beta
NOPI4.beta

NOPI.top <- rownames_to_column(NOPI4.beta, 
                               "Variable")
NOPI.top <- rename(NOPI.top, 
                   "Coefficient"="estimate")

NOPI.top[,1] <- gsub("S:", "", NOPI.top[,1])

str(NOPI.top)

NOPI.top$SEup <- NOPI.top$Coefficient + NOPI.top$se
NOPI.top$SElow <- NOPI.top$Coefficient - NOPI.top$se

NOPI.plot <- ggplot(NOPI.top, aes(x = Variable, 
                                  y = Coefficient)) +
  geom_hline(yintercept = 0,
             colour = gray(1/2), 
             lty = 2) +
  geom_point(aes(x = Variable,
                 y = Coefficient),
             size = 5) +
  geom_errorbar(aes(x = Variable,
                    ymin = lcl,
                    ymax = ucl),
                width = .5,
                size = 1) +
  ggtitle("NOPI Beta Coefficient") +
  theme(panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="white",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="white",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=12, 
                                 colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=12,                                              # change the size of the axis titles
                          colour = "black")) +                                    # change the color of the axis titles
  coord_flip()

# Creating predictive plots -----------------------------------------------

plotdata <- NOPI4.results$S.lit

minlitter <- min(NOPI.surv$Litter)
maxlitter <- max(NOPI.surv$Litter)
litter.values = seq(from = minlitter, 
                    to = maxlitter, 
                    length = 100)

litter.pred <- covariate.predictions(plotdata,
                                     data=data.frame(Litter=litter.values),
                                     indices = 1)

litter.pred$estimates$Group <- "NOPI"

NOPIlitter.plot <- ggplot(litter.pred$estimates, 
                          aes(x = covdata, 
                              y = estimate,
                              fill = Group)) +
  geom_line(size = 1.5,
            aes(colour = Group)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  scale_colour_manual(values = c('black')) +
  scale_fill_manual(values = c('black')) +
  xlab("Litter Cover") +
  ylab("Estimated DSR") +
  ylim(.3, 1) + 
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=16,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill=NA,                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill=NA,                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.y = element_text(size=12, colour = "black"),                    # color the axis text
        axis.text.x = element_text(size=12, colour = "black"),
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=12,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  labs(title = "Northern Pintail")

NOPI.plot
NOPIlitter.plot

ggsave(NOPI.plot,
       filename = "outputs/figs/NOPIbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(NOPIlitter.plot,
       filename = "outputs/figs/NOPIlitter.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all=TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
#  files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)