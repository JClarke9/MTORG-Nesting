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
                               delete = TRUE)
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
  
  MODO.model.list = create.model.list("Nest")
  MODO3.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MODO3.results <- MODO3.run()
MODO3.results

coef(MODO3.results$S.stage)
confint(MODO3.results$S.stage, level = 0.85)


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


MODO.real <- as.data.frame(MODO4.results$S.kbg$results$real) |> 
  summarise(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))


# Plotting beta coefficients ----------------------------------------------


MODO.beta <- coef(MODO4.results$S.kbg) |>
  cbind(confint(MODO4.results$S.kbg, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

MODO.beta$Variable <- gsub("S:", "", MODO.beta$Variable)

str(MODO.beta)

(MODO.plot <- ggplot(MODO.beta[4:5,], 
                     aes(x = Variable,
                         y = Coefficient)) +
    geom_hline(yintercept = 0,
               colour = gray(1/2), 
               lty = 2) +
    geom_point(aes(x = Variable,
                   y = Coefficient),
               size = 4) +
    geom_errorbar(aes(x = Variable,
                      ymin = lcl,
                      ymax = ucl),
                  width = .5,
                  linewidth = 1) +
    theme(plot.title = element_text(family = "my_font",
                                    hjust = .5,
                                    size = 20,
                                    vjust = 1,
                                    colour = "black"),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = NA,                     # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = NA,                      # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text = element_text(size=12, 
                                   colour = "black"),                    # color the axis text
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size=12,                                              # change the size of the axis titles
                            colour = "black")) +                                    # change the color of the axis titles
    labs(title = "Clay-colored Sparrow",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


MODO.ddl <- make.design.data(MODO.pr) |> 
  as.data.frame()

plotdata <- MODO4.results$S.kbg

kbg.values <- seq(from = min(MODO.surv$KBG), 
                 to = max(MODO.surv$KBG), 
                 length = 100)

incub.values <- c(0,1)

time.pred <- covariate.predictions(plotdata,
                                   data = data.frame(KBG = mean(kbg.values),
                                                     Incub = 0),
                                   indices = c(1:72))

(MODOtime.plot <- ggplot(time.pred$estimates, 
                        aes(x = par.index, 
                            y = estimate)) +
    geom_line(linewidth = 1.5,
              aes(color = "model.index")) +
    scale_colour_manual(values = c('#D4A634')) +
    scale_fill_manual(values = c('#D4A634')) +
    theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                    size=16,
                                    hjust=.5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill=NA,                                # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill=NA,                                 # make the outer background transparent
                                         colour=NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_text(size=12, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size=12, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size=12,                                              # change the size of the axis titles
                            colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill=NA),
          legend.position = "none") +
    labs(title = "Mourning Dove",
         color = "Year",
         x = "Julian Day",
         y = "Daily Survival Rate"))


stage.pred <- covariate.predictions(plotdata,
                                    data = data.frame(KBG = mean(kbg.values),
                                                      Incub = incub.values,
                                                      Time = 30),
                                    indices = c(1:72))

(MODOstage.plot <- ggplot(stage.pred$estimates, 
                         aes(x = par.index, 
                             y = estimate,
                             groups = factor(Incub,
                                             levels = c(0,1)))) +
    geom_line(linewidth = 1.5,
              aes(color = factor(Incub,
                                 levels = c(0,1)))) +
    scale_colour_manual(values = c('#717F5B',
                                   '#D4A634')) +
    scale_fill_manual(values = c('#717F5B',
                                 '#D4A634')) +
    theme(plot.title = element_text(family = "my_font",                             # select the font for the title
                                    size = 16,
                                    hjust = .5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = NA,                                # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = NA,                                 # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_text(size = 12, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size = 12,                                              # change the size of the axis titles
                            colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = "none") +
    labs(title = "Mourning Dove",
         color = "Year",
         x = "Stage",
         y = "Daily Survival Rate"))


kbg.pred <- covariate.predictions(plotdata,
                                  data = data.frame(KBG = kbg.values,
                                                    Incub = 0),
                                  indices = c(1, 18, 36, 54, 72))

Day1 <- which(kbg.pred$estimates$par.index == 1)
Day18 <- which(kbg.pred$estimates$par.index == 18)
Day36 <- which(kbg.pred$estimates$par.index == 36)
Day54 <- which(kbg.pred$estimates$par.index == 54)
Day72 <- which(kbg.pred$estimates$par.index == 72)

kbg.pred$estimates$Day <- NA
kbg.pred$estimates$Day[Day1] <- "Day 1"
kbg.pred$estimates$Day[Day18] <- "Day 18"
kbg.pred$estimates$Day[Day36] <- "Day 36"
kbg.pred$estimates$Day[Day54] <- "Day 54"
kbg.pred$estimates$Day[Day72] <- "Day 72"

(MODOkbg.plot <- ggplot(kbg.pred$estimates, 
                         aes(x = KBG, 
                             y = estimate)) +
    geom_line(linewidth = 1.5,
              aes(color = "model.index")) +
    scale_colour_manual(values = c('#D4A634')) +
    scale_fill_manual(values = c('#D4A634')) +
    theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                    size=16,
                                    hjust=.5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill=NA,                                # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill=NA,                                 # make the outer background transparent
                                         colour=NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_text(size=12, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size=12, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size=12,                                              # change the size of the axis titles
                            colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill=NA),
          legend.position = "none") +
    facet_grid(~Day) +
    labs(title = "Mourning Dove",
         color = "Year",
         x = "Kentucky Bluegrass (Percent Cover)",
         y = "Daily Survival Rate"))


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