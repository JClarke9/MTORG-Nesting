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

BWTE.surv <- filter(nest, 
                    Spec=="BWTE")                                         # select out only BWTE nest

test <- filter(BWTE.surv,
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

MISSING <- is.na(BWTE.surv$AgeFound)

sum(MISSING)

BWTE.surv <- subset(BWTE.surv, 
                    subset = !MISSING)

BWTE.surv$Year <- factor(BWTE.surv$Year,
                         levels = c("2021", "2022", "2023"))

str(BWTE.surv)

# Creating stage variable -------------------------------------------------

x <- create.stage.var(BWTE.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(BWTE.surv$LastChecked)), 
                      12)

BWTE.surv <- bind_cols(BWTE.surv, x)

rm(list = ls()[!ls() %in% c("BWTE.surv")])

# Daily survival rate models ----------------------------------------------

BWTE.pr <- process.data(BWTE.surv,
                        nocc=max(BWTE.surv$LastChecked),
                        groups = c("Year"),
                        model="Nest")

# Temporal candidate model set
BWTE1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  BWTE.model.list = create.model.list("Nest")
  BWTE1.results = mark.wrapper(BWTE.model.list,
                               data = BWTE.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BWTE1.results <- BWTE1.run()
BWTE1.results

coef(BWTE1.results$S.year)
confint(BWTE1.results$S.year, level = 0.85)


# Biological candidate model set
BWTE2.run <- function()
{
  # 1. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Year + NestAge)
  
  BWTE.model.list = create.model.list("Nest")
  BWTE2.results = mark.wrapper(BWTE.model.list,
                               data = BWTE.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BWTE2.results <- BWTE2.run()
BWTE2.results

coef(BWTE2.results$S.age)
confint(BWTE2.results$S.age, level = 0.85)


# Grazing candidate model set
BWTE3.run <- function()
{
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Year + NestAge)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + NestAge + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Year + NestAge + grazep)
  
  # 4. DSR varies with the previous Year + NestAges grazing intensity
  S.pTreat = list(formula = ~1 + Year + NestAge + pTreat)
  
  BWTE.model.list = create.model.list("Nest")
  BWTE3.results = mark.wrapper(BWTE.model.list,
                               data = BWTE.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BWTE3.results <- BWTE3.run()
BWTE3.results

coef(BWTE3.results$S.grazep)
confint(BWTE3.results$S.grazep, level = 0.85)

# Vegetation candidate model set
BWTE4.run <- function()
{
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Year + NestAge + grazep)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Year + NestAge + grazep + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Year + NestAge + grazep + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Year + NestAge + grazep + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Year + NestAge + grazep + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Year + NestAge + grazep + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Year + NestAge + grazep + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Year + NestAge + grazep + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + NestAge + grazep + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + NestAge + grazep + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + NestAge + grazep + VOR)
  
  BWTE.model.list = create.model.list("Nest")
  BWTE4.results = mark.wrapper(BWTE.model.list,
                               data = BWTE.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BWTE4.results <- BWTE4.run()
BWTE4.results

coef(BWTE4.results$S.height)
coef(BWTE4.results$S.litdep)
coef(BWTE4.results$S.vor)

confint(BWTE4.results$S.height, level = 0.85)
confint(BWTE4.results$S.litdep, level = 0.85)
confint(BWTE4.results$S.vor, level = 0.85)

BWTE.avg <- model.avg(BWTE4.results$S.height,
                      BWTE4.results$S.litdep,
                      BWTE4.results$S.vor)

summary(BWTE.avg)
coef(BWTE.avg)
confint(BWTE.avg, level = 0.85)

# Vegetation candidate model set
BWTE.run <- function()
{
    # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + NestAge + grazep + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + NestAge + grazep + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + NestAge + grazep + VOR)
  
  BWTE.model.list = create.model.list("Nest")
  BWTE.results = mark.wrapper(BWTE.model.list,
                              data = BWTE.pr,
                              adjust = FALSE,
                              delete = TRUE)
}

BWTE.results <- BWTE.run()
BWTE.results

BWTE.dsr <- model.average(BWTE.results)

BWTE.dsr <- as.data.frame(BWTE.dsr) |> 
  rename("Year" = "par.index")

BWTE.dsr$Year <- case_match(BWTE.dsr$Year,
                            1 ~ "2021",
                            67 ~ "2022",
                            133 ~ "2023")

BWTE.dsr <- fill(BWTE.dsr, Year, .direction = "down")

BWTE.dsr <- BWTE.dsr |> 
  group_by(Year) |> 
  summarise(estimate = mean(estimate))


# Plotting beta coefficients ----------------------------------------------

BWTE.avg

(BWTE.beta <- as.data.frame(t(BWTE.avg$coefficients)))

BWTE.beta <- as.data.frame(t(BWTE.avg$coefficients)) |> 
  cbind(confint(BWTE.avg, level = 0.85, full = T)) |> 
  select(full, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "full",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

BWTE.beta$Variable <- case_match(BWTE.beta$Variable,
                                 "S((Intercept))" ~ "Intercept",
                                 "S(Year2022)" ~ "Year2022",
                                 "S(Year2023)" ~ "Year2023",
                                 "S(NestAge)" ~ "NestAge",
                                 "S(grazep)" ~ "Grazing_Presence",
                                 "S(Veg.Height)" ~ "Veg_Height",
                                 "S(LitterD)" ~ "Litter_Depth",
                                 "S(VOR)" ~ "VOR")

str(BWTE.beta)

(BWTE.plot <- ggplot(BWTE.beta[4:8,], 
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
    labs(title = "Blue-winged Teal",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


BWTE.ddl <- make.design.data(BWTE.pr) |> 
  as.data.frame()

LitDvalues <- seq(from = min(BWTE.surv$LitterD),
                 to = max(BWTE.surv$LitterD),
                 length = 100)

grazep.values <- c(0,1)

VegHvalues <- seq(from = min(BWTE.surv$Veg.Height),
                 to = max(BWTE.surv$Veg.Height),
                 length = 100)

VORvalues <- seq(from = min(BWTE.surv$VOR),
                 to = max(BWTE.surv$VOR),
                 length = 100)


LitD.pred <- covariate.predictions(BWTE.results,
                                   data = data.frame(LitterD = LitDvalues,
                                                     grazep = 0,
                                                     Veg.Height = mean(VegHvalues),
                                                     VOR = mean(VORvalues)),
                                   indices = c(2, 15, 25,
                                               68, 82, 91,
                                               134, 147, 157))

D1Y2021 <- which(LitD.pred$estimates$par.index == 2)
D14Y2021 <- which(LitD.pred$estimates$par.index == 15)
D24Y2021 <- which(LitD.pred$estimates$par.index == 25)
D1Y2022 <- which(LitD.pred$estimates$par.index == 68)
D14Y2022 <- which(LitD.pred$estimates$par.index == 82)
D24Y2022 <- which(LitD.pred$estimates$par.index == 91)
D1Y2023 <- which(LitD.pred$estimates$par.index == 134)
D14Y2023 <- which(LitD.pred$estimates$par.index == 147)
D24Y2023 <- which(LitD.pred$estimates$par.index == 157)

LitD.pred$estimates$Year <- NA
LitD.pred$estimates$Year[D1Y2021] <- "2021"
LitD.pred$estimates$Year[D14Y2021] <- "2021"
LitD.pred$estimates$Year[D24Y2021] <- "2021"
LitD.pred$estimates$Year[D1Y2022] <- "2022"
LitD.pred$estimates$Year[D14Y2022] <- "2022"
LitD.pred$estimates$Year[D24Y2022] <- "2022"
LitD.pred$estimates$Year[D1Y2023] <- "2023"
LitD.pred$estimates$Year[D14Y2023] <- "2023"
LitD.pred$estimates$Year[D24Y2023] <- "2023"

LitD.pred$estimates$Day <- NA
LitD.pred$estimates$Day[D1Y2021] <- "Day1"
LitD.pred$estimates$Day[D14Y2021] <- "Day14"
LitD.pred$estimates$Day[D24Y2021] <- "Day24"
LitD.pred$estimates$Day[D1Y2022] <- "Day1"
LitD.pred$estimates$Day[D14Y2022] <- "Day14"
LitD.pred$estimates$Day[D24Y2022] <- "Day24"
LitD.pred$estimates$Day[D1Y2023] <- "Day1"
LitD.pred$estimates$Day[D14Y2023] <- "Day14"
LitD.pred$estimates$Day[D24Y2023] <- "Day24"

(BWTElitd.plot <- ggplot(transform(LitD.pred$estimates,
                                   Year = factor(Year, levels=c("2021", "2022", "2023"))), 
                         aes(x = LitterD, 
                             y = estimate,
                             groups = Year,
                             fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#A2A4A2',
                                   '#717F5B',
                                   '#D4A634')) +
    scale_fill_manual(values = c('black',
                                 '#717F5B',
                                 '#D4A634')) +
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
          legend.position = c(.85, .1),
          legend.box = "none") +
    facet_grid(~Day) +
    labs(title = "Blue-winged Teal",
         color = "Year",
         x = "Litter Depth (mm)",
         y = "Daily Survival Rate"))


grazep.pred <- covariate.predictions(BWTE.results,
                                     data = data.frame(LitterD = mean(LitDvalues),
                                                       grazep = factor(grazep.values,
                                                                       levels = c(0,1)),
                                                       Veg.Height = mean(VegHvalues),
                                                       VOR = mean(VORvalues)),
                                     indices = c(2, 15, 25,
                                                 68, 82, 91,
                                                 134, 147, 157))

D1Y2021 <- which(grazep.pred$estimates$par.index == 2)
D14Y2021 <- which(grazep.pred$estimates$par.index == 15)
D24Y2021 <- which(grazep.pred$estimates$par.index == 25)
D1Y2022 <- which(grazep.pred$estimates$par.index == 68)
D14Y2022 <- which(grazep.pred$estimates$par.index == 82)
D24Y2022 <- which(grazep.pred$estimates$par.index == 91)
D1Y2023 <- which(grazep.pred$estimates$par.index == 134)
D14Y2023 <- which(grazep.pred$estimates$par.index == 147)
D24Y2023 <- which(grazep.pred$estimates$par.index == 157)

grazep.pred$estimates$Year <- NA
grazep.pred$estimates$Year[D1Y2021] <- "2021"
grazep.pred$estimates$Year[D14Y2021] <- "2021"
grazep.pred$estimates$Year[D24Y2021] <- "2021"
grazep.pred$estimates$Year[D1Y2022] <- "2022"
grazep.pred$estimates$Year[D14Y2022] <- "2022"
grazep.pred$estimates$Year[D24Y2022] <- "2022"
grazep.pred$estimates$Year[D1Y2023] <- "2023"
grazep.pred$estimates$Year[D14Y2023] <- "2023"
grazep.pred$estimates$Year[D24Y2023] <- "2023"

grazep.pred$estimates$Day <- NA
grazep.pred$estimates$Day[D1Y2021] <- "Day1"
grazep.pred$estimates$Day[D14Y2021] <- "Day14"
grazep.pred$estimates$Day[D24Y2021] <- "Day24"
grazep.pred$estimates$Day[D1Y2022] <- "Day1"
grazep.pred$estimates$Day[D14Y2022] <- "Day14"
grazep.pred$estimates$Day[D24Y2022] <- "Day24"
grazep.pred$estimates$Day[D1Y2023] <- "Day1"
grazep.pred$estimates$Day[D14Y2023] <- "Day14"
grazep.pred$estimates$Day[D24Y2023] <- "Day24"

(BWTEgrazep.plot <- ggplot(transform(grazep.pred$estimates,
                                   Year = factor(Year, levels=c("2021", "2022", "2023"))), 
                         aes(x = grazep, 
                             y = estimate,
                             groups = Year,
                             fill = Year)) +
    geom_point(size = 4,
               aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#A2A4A2',
                                   '#717F5B',
                                   '#D4A634')) +
    scale_fill_manual(values = c('black',
                                 '#717F5B',
                                 '#D4A634')) +
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
          legend.position = c(.85, .1),
          legend.box = "none") +
    facet_grid(~Day) +
    labs(title = "Blue-winged Teal",
         color = "Year",
         x = "Vegetation Height (mm)",
         y = "Daily Survival Rate"))


VegH.pred <- covariate.predictions(BWTE.results,
                                   data = data.frame(LitterD = mean(LitDvalues),
                                                     grazep = 0,
                                                     Veg.Height = VegHvalues,
                                                     VOR = mean(VORvalues)),
                                   indices = c(2, 15, 25,
                                               68, 82, 91,
                                               134, 147, 157))

D1Y2021 <- which(VegH.pred$estimates$par.index == 2)
D14Y2021 <- which(VegH.pred$estimates$par.index == 15)
D24Y2021 <- which(VegH.pred$estimates$par.index == 25)
D1Y2022 <- which(VegH.pred$estimates$par.index == 68)
D14Y2022 <- which(VegH.pred$estimates$par.index == 82)
D24Y2022 <- which(VegH.pred$estimates$par.index == 91)
D1Y2023 <- which(VegH.pred$estimates$par.index == 134)
D14Y2023 <- which(VegH.pred$estimates$par.index == 147)
D24Y2023 <- which(VegH.pred$estimates$par.index == 157)

VegH.pred$estimates$Year <- NA
VegH.pred$estimates$Year[D1Y2021] <- "2021"
VegH.pred$estimates$Year[D14Y2021] <- "2021"
VegH.pred$estimates$Year[D24Y2021] <- "2021"
VegH.pred$estimates$Year[D1Y2022] <- "2022"
VegH.pred$estimates$Year[D14Y2022] <- "2022"
VegH.pred$estimates$Year[D24Y2022] <- "2022"
VegH.pred$estimates$Year[D1Y2023] <- "2023"
VegH.pred$estimates$Year[D14Y2023] <- "2023"
VegH.pred$estimates$Year[D24Y2023] <- "2023"

VegH.pred$estimates$Day <- NA
VegH.pred$estimates$Day[D1Y2021] <- "Day1"
VegH.pred$estimates$Day[D14Y2021] <- "Day14"
VegH.pred$estimates$Day[D24Y2021] <- "Day24"
VegH.pred$estimates$Day[D1Y2022] <- "Day1"
VegH.pred$estimates$Day[D14Y2022] <- "Day14"
VegH.pred$estimates$Day[D24Y2022] <- "Day24"
VegH.pred$estimates$Day[D1Y2023] <- "Day1"
VegH.pred$estimates$Day[D14Y2023] <- "Day14"
VegH.pred$estimates$Day[D24Y2023] <- "Day24"

(BWTEVegH.plot <- ggplot(transform(VegH.pred$estimates,
                                   Year = factor(Year, levels=c("2021", "2022", "2023"))), 
                         aes(x = Veg.Height, 
                             y = estimate,
                             groups = Year,
                             fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#A2A4A2',
                                   '#717F5B',
                                   '#D4A634')) +
    scale_fill_manual(values = c('black',
                                 '#717F5B',
                                 '#D4A634')) +
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
          legend.position = c(.85, .1),
          legend.box = "none") +
    facet_grid(~Day) +
    labs(title = "Blue-winged Teal",
         color = "Year",
         x = "Vegetation Height (mm)",
         y = "Daily Survival Rate"))


VOR.pred <- covariate.predictions(BWTE.results,
                                  data = data.frame(LitterD = mean(LitDvalues),
                                                    grazep = 0,
                                                    Veg.Height = mean(VegHvalues),
                                                    VOR = VORvalues),
                                  indices = c(2, 15, 25,
                                              68, 82, 91,
                                              134, 147, 157))

D1Y2021 <- which(VOR.pred$estimates$par.index == 2)
D14Y2021 <- which(VOR.pred$estimates$par.index == 15)
D24Y2021 <- which(VOR.pred$estimates$par.index == 25)
D1Y2022 <- which(VOR.pred$estimates$par.index == 68)
D14Y2022 <- which(VOR.pred$estimates$par.index == 82)
D24Y2022 <- which(VOR.pred$estimates$par.index == 91)
D1Y2023 <- which(VOR.pred$estimates$par.index == 134)
D14Y2023 <- which(VOR.pred$estimates$par.index == 147)
D24Y2023 <- which(VOR.pred$estimates$par.index == 157)

VOR.pred$estimates$Year <- NA
VOR.pred$estimates$Year[D1Y2021] <- "2021"
VOR.pred$estimates$Year[D14Y2021] <- "2021"
VOR.pred$estimates$Year[D24Y2021] <- "2021"
VOR.pred$estimates$Year[D1Y2022] <- "2022"
VOR.pred$estimates$Year[D14Y2022] <- "2022"
VOR.pred$estimates$Year[D24Y2022] <- "2022"
VOR.pred$estimates$Year[D1Y2023] <- "2023"
VOR.pred$estimates$Year[D14Y2023] <- "2023"
VOR.pred$estimates$Year[D24Y2023] <- "2023"

VOR.pred$estimates$Day <- NA
VOR.pred$estimates$Day[D1Y2021] <- "Day1"
VOR.pred$estimates$Day[D14Y2021] <- "Day14"
VOR.pred$estimates$Day[D24Y2021] <- "Day24"
VOR.pred$estimates$Day[D1Y2022] <- "Day1"
VOR.pred$estimates$Day[D14Y2022] <- "Day14"
VOR.pred$estimates$Day[D24Y2022] <- "Day24"
VOR.pred$estimates$Day[D1Y2023] <- "Day1"
VOR.pred$estimates$Day[D14Y2023] <- "Day14"
VOR.pred$estimates$Day[D24Y2023] <- "Day24"

(BWTEVOR.plot <- ggplot(transform(VOR.pred$estimates,
                                   Year = factor(Year, levels=c("2021", "2022", "2023"))), 
                         aes(x = VOR, 
                             y = estimate,
                             groups = Year,
                             fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#A2A4A2',
                                   '#717F5B',
                                   '#D4A634')) +
    scale_fill_manual(values = c('black',
                                 '#717F5B',
                                 '#D4A634')) +
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
          legend.position = c(.85, .1),
          legend.box = "none") +
    facet_grid(~Day) +
    labs(title = "Blue-winged Teal",
         color = "Year",
         x = "Visual Obstruction Reading (dm)",
         y = "Daily Survival Rate"))


AGE.pred <- covariate.predictions(BWTE.results,
                                  data = data.frame(LitterD = mean(LitDvalues),
                                                    grazep = 0,
                                                    Veg.Height = mean(VegHvalues),
                                                    VOR = mean(VORvalues)),
                                  indices = c(2:25,
                                              68:91,
                                              134:157))

D1Y2021 <- which(AGE.pred$estimates$par.index == 2)
D1Y2022 <- which(AGE.pred$estimates$par.index == 68)
D1Y2023 <- which(AGE.pred$estimates$par.index == 134)

AGE.pred$estimates$Year <- NA
AGE.pred$estimates$Year[D1Y2021] <- "2021"
AGE.pred$estimates$Year[D1Y2022] <- "2022"
AGE.pred$estimates$Year[D1Y2023] <- "2023"

AGE.pred$estimates <- fill(AGE.pred$estimates, Year, .direction = "down")

AGE.pred$estimates$Day <- c(1:24)

(BWTEage.plot <- ggplot(transform(AGE.pred$estimates,
                                  Year = factor(Year, levels=c("2021", "2022", "2023"))), 
                        aes(x = Day, 
                            y = estimate,
                            groups = Year,
                            fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#A2A4A2',
                                   '#717F5B',
                                   '#D4A634')) +
    scale_fill_manual(values = c('black',
                                 '#717F5B',
                                 '#D4A634')) +
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
          legend.position = c(.85, .1),
          legend.box = "none") +
    labs(title = "Blue-winged Teal",
         color = "Year",
         x = "Nest Age",
         y = "Daily Survival Rate"))


ggsave(BWTE.plot,
       filename = "~/Git/NDSU/RMARK/Figures/BWTEbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BWTElitdep.plot,
       filename = "outputs/figs/BWTElitdep.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BWTEgrazep.plot,
       filename = "outputs/figs/BWTEgrazep.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BWTEvegh.plot,
       filename = "outputs/figs/BWTEvegh.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BWTEvor.plot,
       filename = "outputs/figs/BWTEvor.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BWTEage.plot,
       filename = "outputs/figs/BWTEage.png",
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
