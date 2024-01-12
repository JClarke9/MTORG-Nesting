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


WEME.surv <- filter(nest, 
                    Spec=="WEME")                                         # select out only WEME nest

test <- filter(WEME.surv,
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

MISSING <- is.na(WEME.surv$AgeFound)

sum(MISSING)

WEME.surv <- subset(WEME.surv, 
                    subset = !MISSING)

WEME.surv$Year <- factor(WEME.surv$Year,
                         levels = c("2021", "2022", "2023"))

str(WEME.surv)


# Creating stage variable -------------------------------------------------


x <- create.stage.var(WEME.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(WEME.surv$LastChecked)), 
                      12)

WEME.surv <- bind_cols(WEME.surv, x)

rm(list = ls()[!ls() %in% c("WEME.surv")])


# Daily survival rate models ----------------------------------------------


WEME.pr <- process.data(WEME.surv,
                        nocc=max(WEME.surv$LastChecked),
                        groups = c("Year"),
                        model="Nest")

# Temporal candidate model set
WEME1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  WEME.model.list = create.model.list("Nest")
  WEME1.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME1.results <- WEME1.run()
WEME1.results

coef(WEME1.results$S.null)


# Biological candidate model set
WEME2.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + BHCONum)
  
  # 2. DSR varies with BHCO number
  S.bhcop = list(formula = ~1 + BHCOPres)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + NestAge)
  
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Incub)
  
  WEME.model.list = create.model.list("Nest")
  WEME2.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME2.results <- WEME2.run()
WEME2.results

coef(WEME2.results$S.null)


# Grazing candidate model set
WEME3.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pTreat = list(formula = ~1 + pTreat)
  
  WEME.model.list = create.model.list("Nest")
  WEME3.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME3.results <- WEME3.run()
WEME3.results

coef(WEME3.results$S.grazed)
confint(WEME3.results$S.grazed, level = 0.85)


# Vegetation candidate model set
WEME4.run <- function()
{
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + grazed)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + grazed + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + grazed + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + grazed + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + grazed + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + grazed + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + grazed + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + grazed + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + grazed + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + grazed + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + grazed + VOR)
  
  WEME.model.list = create.model.list("Nest")
  WEME4.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME4.results <- WEME4.run()
WEME4.results

coef(WEME4.results$S.litdep)
confint(WEME4.results$S.litdep, level = 0.85)


# Plotting beta coefficients ----------------------------------------------


WEME.beta <- coef(WEME4.results$S.litdep) |>
  cbind(confint(WEME4.results$S.litdep, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

WEME.beta$Variable <- gsub("S:", "", WEME.beta$Variable)

str(WEME.beta)

(WEME.plot <- ggplot(WEME.beta[2:3,], 
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
    labs(title = "Western Meadowlark",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


WEME.ddl <- make.design.data(WEME.pr) |> 
  as.data.frame()

plotdata <- WEME4.results$S.litdep

Grazedvalues <- seq(from = min(WEME.surv$grazed),
                    to = max(WEME.surv$grazed),
                    length = 100)

LitDvalues <- seq(from = min(WEME.surv$LitterD),
                  to = max(WEME.surv$LitterD),
                  length = 100)


Grazed.pred <- covariate.predictions(plotdata,
                                     data = data.frame(grazed = Grazedvalues,
                                                       LitterD = mean(LitDvalues)),
                                     indices = 1)

(WEMEgrazed.plot <- ggplot(Grazed.pred$estimates, 
                           aes(x = grazed, 
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
    labs(title = "Western Meadowlark",
         color = "Year",
         x = "Days Grazed",
         y = "Daily Survival Rate"))


LitterD.pred <- covariate.predictions(plotdata,
                                     data = data.frame(grazed = mean(Grazedvalues),
                                                       LitterD = LitDvalues),
                                     indices = 1)

(WEMElitd.plot <- ggplot(LitterD.pred$estimates, 
                         aes(x = LitterD, 
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
    labs(title = "Western Meadowlark",
         color = "Year",
         x = "Litter Depth (mm)",
         y = "Daily Survival Rate"))


ggsave(WEME.plot,
       filename = "outputs/figs/WEMEbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(WEMEgrazed.plot,
       filename = "outputs/figs/WEMEgrazed.png",
       dpi = "print",
       height = 3,
       width = 5)

ggsave(WEMElitd.plot,
       filename = "outputs/figs/WEMElitterD.png",
       dpi = "print",
       height = 3,
       width = 5)


# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all=TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
#  files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)
