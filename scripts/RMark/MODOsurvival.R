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

MODO.surv <- filter(nest, 
                    Spec=="MODO")                     # select out only MODO nests

MODO.surv$AgeDay1 <- MODO.surv$AgeFound - MODO.surv$FirstFound + 1
MODO.surv$Year <- as.factor(MODO.surv$Year)
MODO.surv$cTreat <- as.factor(MODO.surv$cTreat)

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

x <- create.stage.var(MODO.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(MODO.surv$LastChecked)), 
                      14)

MODO.surv <- bind_cols(MODO.surv, x)

rm(list = ls()[!ls() %in% "MODO.surv"])

# Daily survival rate models ----------------------------------------------

MODO.pr <- process.data(MODO.surv,
                        nocc=max(MODO.surv$LastChecked), 
                        groups = c("cTreat",
                                   "Year"),
                        model="Nest")

# candidate models for the null vs. treatment model
MODO1.run <- function()
{
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  # 3. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies With year interacting with treatment (treatment identity varies by year)
  S.treatyear = list(formula = ~1 + cTreat:Year)
  
  MODO.model.list = create.model.list("Nest")
  MODO1.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
MODO1.results <- MODO1.run()
MODO1.results

MODO1.results$S.Dot$results$beta
MODO1.results$S.year$results$beta

#candidate model set for time trends
MODO2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies across stages
  S.stage = list(formula = ~1 + Incub)
  
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  MODO.model.list = create.model.list("Nest")
  MODO2.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
MODO2.results <- MODO2.run()
MODO2.results

MODO2.results$S.quad$results$beta
MODO2.results$S.Dot$results$beta

#candidate model set for grazing
MODO3.run <- function()
{
  S.grazed = list(formula = ~1 + Time + I(Time^2) + grazed)
  
  S.grazep = list(formula = ~1 + Time + I(Time^2) + grazep)
  
  # 2. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  MODO.model.list = create.model.list("Nest")
  MODO3.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
MODO3.results <- MODO3.run()
MODO3.results

MODO3.results$S.quad$results$beta
MODO3.results$S.grazed$results$beta
MODO3.results$S.grazep$results$beta

#candidate model set for parasitism and nest BHCONum
MODO4.run <- function()
{
  # 1. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Time + I(Time^2) + KBG)
  
  # 2. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.smooth.brome = list(formula = ~1 + Time + I(Time^2) + SmoothB)
  
  # 3. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Time + I(Time^2) + Litter)
  
  # 4. DSR varies with Bare
  S.bare = list(formula =  ~1 + Time + I(Time^2) + Bare)
  
  # 5. DSR varies with Forb
  S.forb = list(formula =  ~1 + Time + I(Time^2) + Forb)
  
  # 6. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Time + I(Time^2) + Grasslike)
  
  # 7. DSR varies with Woody
  S.woody = list(formula =  ~1 + Time + I(Time^2) + Woody)
  
  # 8. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Time + I(Time^2) + LitterD)
  
  # 9. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Time + I(Time^2) + Veg.Height)
  
  # 10. DSR varies with VOR
  S.vor = list(formula =  ~1 + Time + I(Time^2) + VOR)
  
  # 11. DSR varies with litter and Veg Height
  S.litheight = list(formula = ~1 + Time + I(Time^2) + Litter + Veg.Height)
  
  # 12. DSR varies with woody and litter depth
  S.woodylitdep = list(formula = ~1 + Time + I(Time^2) + Woody + LitterD)
  
  # 13. DSR varies with KBG and litter depth
  S.kbglitdep = list(formula = ~1 + Time + I(Time^2) + KBG + LitterD)
  
  # 14. DSR varies with KBG and Veg.Height
  S.kbgheight = list(formula = ~1 + Time + I(Time^2) + KBG + Veg.Height)
  
  # 15. DSR varies with woody and Veg Height
  S.woodyheight = list(formula = ~1 + Time + I(Time^2) + Woody + Veg.Height)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  MODO.model.list = create.model.list("Nest")
  MODO4.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

#results of candidate model set 
#all results
MODO4.results <- MODO4.run()
MODO4.results
names(MODO4.results)

MODO4.results$S.kbglitdep$results$beta
MODO4.results$S.kbg$results$beta
MODO4.results$S.height$results$beta

MODO.real <- as.data.frame(MODO4.results$S.kbglitdep$results$real)

mean(MODO.real$estimate)

# Plotting beta coefficients ----------------------------------------------

MODO.mod <- mark(MODO.surv, 
                 nocc=max(MODO.surv$LastChecked), 
                 model = "Nest", 
                 model.parameters = list(S = list(formula = ~1 + Time + I(Time^2) + KBG + LitterD)))

MODO4.beta <- MODO4.results$S.kbglitdep$results$beta
MODO4.beta

MODO.top <- rownames_to_column(MODO4.beta, 
                               "Variable")
MODO.top <- rename(MODO.top, 
                   "Coefficient"="estimate")
str(MODO.top)

MODO.top[,1] <- gsub("S:", "", MODO.top[,1])


MODO.top$SEup <- MODO.top$Coefficient + MODO.top$se
MODO.top$SElow <- MODO.top$Coefficient - MODO.top$se

MODO.plot <- ggplot(MODO.top, aes(x = Variable,
                                  y = Coefficient)) +
  geom_hline(yintercept = 0,
             colour = gray(1/2), lty = 2) +
  geom_point(aes(x = Variable,
                 y = Coefficient),
             size = 5) +
  geom_errorbar(aes(x = Variable,
                    ymin = lcl,
                    ymax = ucl),
                width = .5,
                size = 1) +
  ggtitle("MODO Beta Coefficient") +
  theme(panel.grid.major = element_blank(),                   # remove the vertical grid lines
        panel.grid.minor = element_blank(),                   # remove the horizontal grid lines
        panel.background = element_rect(fill="white",           # make the interior background transparent
                                        colour = NA),              # remove any other colors
        plot.background = element_rect(fill="white",           # make the outer background transparent
                                       colour=NA),               # remove any other colors
        axis.line = element_line(colour = "black"),               # color the x and y axis
        axis.text = element_text(size=12, colour = "black"),          # color the axis text
        axis.ticks = element_line(colour = "black"),              # change the colors of the axis tick marks
        text=element_text(size=12,                       # change the size of the axis titles
                          colour = "black")) +                  # change the color of the axis titles
  coord_flip()

# Creating predictive plots -----------------------------------------------

MODO.ddl <- make.design.data(MODO.pr)
MODO.ddl <- as.data.frame(MODO.ddl)

plotdata <- MODO4.results$S.kbglitdep

minkbg <- min(MODO.surv$KBG)
maxkbg <- max(MODO.surv$KBG)
kbg.values = seq(from = minkbg, 
                 to = maxkbg, 
                 length = 100)

minlitdep <- min(MODO.surv$LitterD)
maxlitdep <- max(MODO.surv$LitterD)
litdep.values = seq(from = minlitdep, 
                    to = maxlitdep, 
                    length = 100)

kbg.pred <- covariate.predictions(plotdata,
                                  data=data.frame(KBG=kbg.values,
                                                  LitterD=mean(litdep.values)),
                                  indices = c(6, 21, 41, 61))

Day5 <- which(kbg.pred$estimates$par.index == 6)
Day20 <- which(kbg.pred$estimates$par.index == 21)
Day40 <- which(kbg.pred$estimates$par.index == 41)
Day60 <- which(kbg.pred$estimates$par.index == 61)

kbg.pred$estimates$Date <- NA
kbg.pred$estimates$Date[Day5] <- "Early"
kbg.pred$estimates$Date[Day20] <- "Mid1"
kbg.pred$estimates$Date[Day40] <- "Mid2"
kbg.pred$estimates$Date[Day60] <- "Late"

MODOkbg.plot <- ggplot(transform(kbg.pred$estimates,
                                 Date = factor(Date, levels=c("Early", "Mid1", "Mid2" , "Late"))), 
                       aes(x = KBG, 
                           y = estimate,
                           groups = Date)) +
  geom_line(size = 1.5,
            aes(linetype = Date)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  scale_colour_manual(values = c('black')) +
  scale_fill_manual(values = c('black')) +
  xlab("KBG Cover") +
  ylab("Estimated DSR") +
  ylim(0, 1) + 
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
        legend.background = element_rect(fill=NA),
        legend.position = c(.1, .1),
        legend.box = "horizontal") +
  labs(title = "Mourning Dove")

litdep.pred <- covariate.predictions(plotdata,
                                     data=data.frame(LitterD=litdep.values,
                                                     KBG=mean(kbg.values)),
                                     indices = c(6, 21, 41, 61))

Day5 <- which(litdep.pred$estimates$par.index == 6)
Day20 <- which(litdep.pred$estimates$par.index == 21)
Day40 <- which(litdep.pred$estimates$par.index == 41)
Day60 <- which(litdep.pred$estimates$par.index == 61)

litdep.pred$estimates$Date <- NA
litdep.pred$estimates$Date[Day5] <- "Early"
litdep.pred$estimates$Date[Day20] <- "Mid1"
litdep.pred$estimates$Date[Day40] <- "Mid2"
litdep.pred$estimates$Date[Day60] <- "Late"

MODOlitdep.plot <- ggplot(transform(litdep.pred$estimates,
                                    Date = factor(Date, levels=c("Early", "Mid1", "Mid2", "Late"))),
                          aes(x = LitterD, 
                              y = estimate,
                              groups = Date)) +
  geom_line(size = 1.5,
            aes(linetype = Date)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  scale_colour_manual(values = c('black')) +
  scale_fill_manual(values = c('black')) +
  xlab("Litter Depth (mm)") +
  ylab("Estimated DSR") +
  ylim(0, 1) + 
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
        legend.background = element_rect(fill=NA),
        legend.position = c(.1, .1),
        legend.box = "horizontal") +
  labs(title = "Mourning Dove")

julian.pred <- covariate.predictions(plotdata,
                                     data=data.frame(LitterD=mean(litdep.values),
                                                     KBG=mean(kbg.values)),
                                     indices = c(1:71))

x <- which(julian.pred$estimates$par.index == c(1:71))

julian.pred$estimates$Group[x] <- "MODO"

julian.pred$estimates$Julian <- c(1:71)

MODOjulian.plot <- ggplot(julian.pred$estimates, 
                          aes(x = Julian, 
                              y = estimate,
                              fill = Group)) +
  geom_line(size = 1.5,
            aes(colour = Group)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  scale_colour_manual(values = c('black')) +
  scale_fill_manual(values = c('black')) +
  xlab("Julian Day") +
  ylab("Estimated DSR") +
  ylim(0, 1) + 
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
        legend.background = element_rect(fill=NA),
        legend.position = "none",
        legend.box = "horizontal") +                                             # remove the legend
  labs(title = "Mourning Dove")

MODO.plot
MODOkbg.plot
MODOlitdep.plot
MODOjulian.plot

ggsave(MODO.plot,
       filename = "outputs/figs/MODObeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(MODOkbg.plot,
       filename = "outputs/figs/MODOkbg.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(MODOlitdep.plot,
       filename = "outputs/figs/MODOlitdep.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(MODOjulian.plot,
       filename = "outputs/figs/MODOjulian.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

# If you want to clean up the mark*.inp, .vcv, .res and .out
# and .tmp files created by RMark in the working directory,
# execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all=TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
# files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)