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

MODO.surv <- filter(nest, 
                    Spec=="MODO")                                         # select out only MODO nest

test <- filter(MODO.surv,
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

MISSING <- is.na(MODO.surv$AgeFound)

sum(MISSING)

MODO.surv <- subset(MODO.surv, 
                    subset = !MISSING)

MODO.surv$Year <- factor(MODO.surv$Year,
                         levels = c("2021", "2022", "2023"))

str(MODO.surv)

# Creating stage variable -------------------------------------------------

x <- create.stage.var(MODO.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(MODO.surv$LastChecked)), 
                      12)

MODO.surv <- bind_cols(MODO.surv, x)

rm(list = ls()[!ls() %in% c("MODO.surv")])

# Daily survival rate models ----------------------------------------------

MODO.pr <- process.data(MODO.surv,
                        nocc=max(MODO.surv$LastChecked),
                        groups = c("Year"),
                        model="Nest")

# Temporal candidate model set
MODO1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  MODO.model.list = create.model.list("Nest")
  MODO1.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

# Results of candidate model set
MODO1.results <- MODO1.run()
MODO1.results

coef(MODO1.results$S.quad)
confint(MODO1.results$S.quad, level = 0.85)


# Biological candidate model set
MODO2.run <- function()
{
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Time + I(Time^2) + NestAge)
  
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Time + I(Time^2) + Incub)
  
  MODO.model.list = create.model.list("Nest")
  MODO2.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MODO2.results <- MODO2.run()
MODO2.results

coef(MODO2.results$S.stage)
confint(MODO2.results$S.stage, level = 0.85)


# Grazing candidate model set
MODO3.run <- function()
{
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Time + I(Time^2) + Incub)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Time + I(Time^2) + Incub + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Time + I(Time^2) + Incub + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pTreat = list(formula = ~1 + Time + I(Time^2) + Incub + pTreat)
  
  # 4. DSR varies with the previous years grazing intensity
  S.grazedpTreat = list(formula = ~1 + Time + I(Time^2) + Incub + grazed + pTreat)
  
  MODO.model.list = create.model.list("Nest")
  MODO3.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MODO3.results <- MODO3.run()
MODO3.results

# Vegetation candidate model set
MODO4.run <- function()
{
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Time + I(Time^2) + Incub)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Time + I(Time^2) + Incub + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Time + I(Time^2) + Incub + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Time + I(Time^2) + Incub + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Time + I(Time^2) + Incub + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Time + I(Time^2) + Incub + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Time + I(Time^2) + Incub + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Time + I(Time^2) + Incub + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Time + I(Time^2) + Incub + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Time + I(Time^2) + Incub + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Time + I(Time^2) + Incub + VOR)
  
  MODO.model.list = create.model.list("Nest")
  MODO4.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MODO4.results <- MODO4.run()
MODO4.results

coef(MODO4.results$S.kbg)
confint(MODO4.results$S.kbg, level = 0.85)

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