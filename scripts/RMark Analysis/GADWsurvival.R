# Loading libraries -------------------------------------------------------


library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)
library(MuMIn)
source("scripts/Functions/RMark_Stage_Code.R")

windowsFonts(my_font = windowsFont("Gandhi Sans"))


# Data import -------------------------------------------------------------


nest <- read.csv("working/RMarknesting.csv")


# Subsetting data ---------------------------------------------------------


GADW.surv <- filter(nest, Spec=="GADW" & Stage != "Laying")                                         # select out only GADW nest

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
                         levels = c("2021", "2022", "2023", "2024"))

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
                        nocc = max(GADW.surv$LastChecked),
                        groups = c("Year"),
                        model = "Nest")



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
                               delete = TRUE)
}

# Results of candidate model set
GADW1.results <- GADW1.run()
GADW1.results

coef(GADW1.results$S.year)
confint(GADW1.results$S.year, level = 0.85)



# Biological candidate model set
GADW2.run <- function()
{
  # 4. DSR varies with year
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
  # 1. DSR varies with nest age
  S.age = list(formula = ~1 + Year + NestAge)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + NestAge + grazed)
  
  # 4. DSR varies with the previous Year + NestAges grazing intensity
  S.pDoD = list(formula = ~1 + Year + NestAge + pDoD)
  
  GADW.model.list = create.model.list("Nest")
  GADW3.results = mark.wrapper(GADW.model.list,
                               data = GADW.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
GADW3.results <- GADW3.run()
GADW3.results

coef(GADW3.results$S.grazed)
confint(GADW3.results$S.grazed, level = 0.85)

coef(GADW3.results$S.pDoD)
confint(GADW3.results$S.pDoD, level = 0.85)


# Vegetation candidate model set
GADW4.run <- function()
{
  # 1. DSR varies with grazing
  S.grazing = list(formula = ~1 + Year + NestAge + grazed + pDoD)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Year + NestAge + grazed + pDoD + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Year + NestAge + grazed + pDoD + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Year + NestAge + grazed + pDoD + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Year + NestAge + grazed + pDoD + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Year + NestAge + grazed + pDoD + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Year + NestAge + grazed + pDoD + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Year + NestAge + grazed + pDoD + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + NestAge + grazed + pDoD + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + NestAge + grazed + pDoD + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + NestAge + grazed + pDoD + VOR)
  
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

GADW4.results$S.forb$results$real |>
  summarize(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))

GADW.real <- as.data.frame(GADW4.results$S.forb$results$real) |> 
  rownames_to_column(var = "Group") |> 
  mutate(Year = case_when(
    grepl("2021", Group) ~ "2021",
    grepl("2022", Group) ~ "2022",
    grepl("2023", Group) ~ "2023",
    grepl("2024", Group) ~ "2024")) |> 
  select(Year, estimate, se, lcl, ucl)

(GADW.dsr <- GADW.real |> 
    group_by(Year) |> 
    summarise(estimate = mean(estimate),
              se = mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl)))


# Plotting beta coefficients ----------------------------------------------


GADW.beta <- coef(GADW4.results$S.forb) |>
  cbind(confint(GADW4.results$S.forb, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

GADW.beta$Variable <- gsub("S:", "", GADW.beta$Variable)

str(GADW.beta)

(GADW.plot <- ggplot(GADW.beta[c(6:8),], 
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
          axis.text = element_text(size = 12, 
                                   colour = "black"),                    # color the axis text
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black")) +                                    # change the color of the axis titles
    labs(title = "Gadwall",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


GADW.ddl <- make.design.data(GADW.pr) |> 
  as.data.frame()

plotdata <- GADW4.results$S.forb

Forbvalues <- seq(from = min(GADW.surv$Forb),
                  to = max(GADW.surv$Forb),
                  length = 100)

Grazedvalues <- seq(from = min(GADW.surv$grazed),
                    to = max(GADW.surv$grazed),
                    length = 100)

pDoDvalues <- seq(from = min(GADW.surv$pDoD),
                  to = max(GADW.surv$pDoD),
                  length = 100)


filter(GADW.ddl, S.Year == 2021 & S.age %in% c(1:27))
filter(GADW.ddl, S.Year == 2022 & S.age %in% c(1:27))
filter(GADW.ddl, S.Year == 2023 & S.age %in% c(1:27))
filter(GADW.ddl, S.Year == 2024 & S.age %in% c(1:27))


AGE.pred <- covariate.predictions(plotdata,
                                  data = data.frame(Forb = mean(Forbvalues),
                                                    grazed = mean(Grazedvalues),
                                                    pDoD = mean(pDoDvalues)),
                                  indices = c(2:28,
                                              77:103,
                                              152:178,
                                              227:253))

D1Y2021 <- which(AGE.pred$estimates$par.index == 2)
D1Y2022 <- which(AGE.pred$estimates$par.index == 77)
D1Y2023 <- which(AGE.pred$estimates$par.index == 152)
D1Y2024 <- which(AGE.pred$estimates$par.index == 227)

AGE.pred$estimates$Year <- NA
AGE.pred$estimates$Year[D1Y2021] <- "2021"
AGE.pred$estimates$Year[D1Y2022] <- "2022"
AGE.pred$estimates$Year[D1Y2023] <- "2023"
AGE.pred$estimates$Year[D1Y2024] <- "2024"

AGE.pred$estimates <- fill(AGE.pred$estimates, Year, .direction = "down")

AGE.pred$estimates$Day <- c(1:27)

(GADWage.plot <- ggplot(transform(AGE.pred$estimates,
                                  Year = factor(Year, levels = c("2021", "2022", "2023", "2024"))), 
                        aes(x = Day, 
                            y = estimate,
                            groups = Year,
                            fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c("#A2A4A2",
                                   "#717F5B",
                                   "#D4A634",
                                   "goldenrod4")) +
    scale_fill_manual(values = c("#A2A4A2",
                                 "#717F5B",
                                 "#D4A634",
                                 "goldenrod4")) +
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
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = c(.85, .1),
          legend.box = "horizontal") +
    labs(title = "Nest Survival",
         color = "Year",
         x = "Nest Age",
         y = "Daily Survival Rate"))



Forb.pred <- covariate.predictions(plotdata,
                                   data = data.frame(Forb = Forbvalues,
                                                     grazed = mean(Grazedvalues),
                                                     pDoD = mean(pDoDvalues)),
                                   indices = c(2, 16, 28,
                                               77, 91, 103,
                                               152, 166, 178,
                                               227, 241, 253))

D1Y2021 <- which(Forb.pred$estimates$par.index == 2)
D15Y2021 <- which(Forb.pred$estimates$par.index == 16)
D27Y2021 <- which(Forb.pred$estimates$par.index == 28)
D1Y2022 <- which(Forb.pred$estimates$par.index == 77)
D15Y2022 <- which(Forb.pred$estimates$par.index == 91)
D27Y2022 <- which(Forb.pred$estimates$par.index == 103)
D1Y2023 <- which(Forb.pred$estimates$par.index == 152)
D15Y2023 <- which(Forb.pred$estimates$par.index == 166)
D27Y2023 <- which(Forb.pred$estimates$par.index == 178)
D1Y2024 <- which(Forb.pred$estimates$par.index == 227)
D15Y2024 <- which(Forb.pred$estimates$par.index == 241)
D27Y2024 <- which(Forb.pred$estimates$par.index == 253)

Forb.pred$estimates$Year <- NA
Forb.pred$estimates$Year[D1Y2021] <- "2021"
Forb.pred$estimates$Year[D15Y2021] <- "2021"
Forb.pred$estimates$Year[D27Y2021] <- "2021"
Forb.pred$estimates$Year[D1Y2022] <- "2022"
Forb.pred$estimates$Year[D15Y2022] <- "2022"
Forb.pred$estimates$Year[D27Y2022] <- "2022"
Forb.pred$estimates$Year[D1Y2023] <- "2023"
Forb.pred$estimates$Year[D15Y2023] <- "2023"
Forb.pred$estimates$Year[D27Y2023] <- "2023"
Forb.pred$estimates$Year[D1Y2024] <- "2024"
Forb.pred$estimates$Year[D15Y2024] <- "2024"
Forb.pred$estimates$Year[D27Y2024] <- "2024"

Forb.pred$estimates$Day <- NA
Forb.pred$estimates$Day[D1Y2021] <- "Day1"
Forb.pred$estimates$Day[D15Y2021] <- "Day15"
Forb.pred$estimates$Day[D27Y2021] <- "Day27"
Forb.pred$estimates$Day[D1Y2022] <- "Day1"
Forb.pred$estimates$Day[D15Y2022] <- "Day15"
Forb.pred$estimates$Day[D27Y2022] <- "Day27"
Forb.pred$estimates$Day[D1Y2023] <- "Day1"
Forb.pred$estimates$Day[D15Y2023] <- "Day15"
Forb.pred$estimates$Day[D27Y2023] <- "Day27"
Forb.pred$estimates$Day[D1Y2024] <- "Day1"
Forb.pred$estimates$Day[D15Y2024] <- "Day15"
Forb.pred$estimates$Day[D27Y2024] <- "Day27"

(GADWforb.plot <- ggplot(transform(Forb.pred$estimates,
                                   Year = factor(Year, levels = c("2021", "2022", "2023", "2024"))), 
                         aes(x = Forb, 
                             y = estimate,
                             groups = Year,
                             fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c("#A2A4A2",
                                   "#717F5B",
                                   "#D4A634",
                                   "goldenrod4")) +
    scale_fill_manual(values = c("#A2A4A2",
                                 "#717F5B",
                                 "#D4A634",
                                 "goldenrod4")) +
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
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill=NA),
          legend.position = c(.85, .1),
          legend.box = "horizontal") +
    facet_grid(~Day) +
    labs(title = "Nest Survival",
         color = "Year",
         x = "Forb (Percent Cover)",
         y = "Daily Survival Rate"))


grazed.pred <- covariate.predictions(plotdata,
                                     data = data.frame(Forb = mean(Forbvalues),
                                                       grazed = Grazedvalues,
                                                       pDoD = mean(pDoDvalues)),
                                     indices = c(2, 16, 28,
                                                 77, 91, 103,
                                                 152, 166, 178,
                                                 227, 241, 253))

D1Y2021 <- which(grazed.pred$estimates$par.index == 2)
D15Y2021 <- which(grazed.pred$estimates$par.index == 16)
D27Y2021 <- which(grazed.pred$estimates$par.index == 28)
D1Y2022 <- which(grazed.pred$estimates$par.index == 77)
D15Y2022 <- which(grazed.pred$estimates$par.index == 91)
D27Y2022 <- which(grazed.pred$estimates$par.index == 103)
D1Y2023 <- which(grazed.pred$estimates$par.index == 152)
D15Y2023 <- which(grazed.pred$estimates$par.index == 166)
D27Y2023 <- which(grazed.pred$estimates$par.index == 178)
D1Y2024 <- which(grazed.pred$estimates$par.index == 227)
D15Y2024 <- which(grazed.pred$estimates$par.index == 241)
D27Y2024 <- which(grazed.pred$estimates$par.index == 253)

grazed.pred$estimates$Year <- NA
grazed.pred$estimates$Year[D1Y2021] <- "2021"
grazed.pred$estimates$Year[D15Y2021] <- "2021"
grazed.pred$estimates$Year[D27Y2021] <- "2021"
grazed.pred$estimates$Year[D1Y2022] <- "2022"
grazed.pred$estimates$Year[D15Y2022] <- "2022"
grazed.pred$estimates$Year[D27Y2022] <- "2022"
grazed.pred$estimates$Year[D1Y2023] <- "2023"
grazed.pred$estimates$Year[D15Y2023] <- "2023"
grazed.pred$estimates$Year[D27Y2023] <- "2023"
grazed.pred$estimates$Year[D1Y2024] <- "2024"
grazed.pred$estimates$Year[D15Y2024] <- "2024"
grazed.pred$estimates$Year[D27Y2024] <- "2024"

grazed.pred$estimates$Day <- NA
grazed.pred$estimates$Day[D1Y2021] <- "Day1"
grazed.pred$estimates$Day[D15Y2021] <- "Day15"
grazed.pred$estimates$Day[D27Y2021] <- "Day27"
grazed.pred$estimates$Day[D1Y2022] <- "Day1"
grazed.pred$estimates$Day[D15Y2022] <- "Day15"
grazed.pred$estimates$Day[D27Y2022] <- "Day27"
grazed.pred$estimates$Day[D1Y2023] <- "Day1"
grazed.pred$estimates$Day[D15Y2023] <- "Day15"
grazed.pred$estimates$Day[D27Y2023] <- "Day27"
grazed.pred$estimates$Day[D1Y2024] <- "Day1"
grazed.pred$estimates$Day[D15Y2024] <- "Day15"
grazed.pred$estimates$Day[D27Y2024] <- "Day27"

(GADWgrazed.plot <- ggplot(transform(grazed.pred$estimates,
                                     Year = factor(Year, levels = c("2021", "2022", "2023", "2024"))), 
                           aes(x = grazed, 
                               y = estimate,
                               groups = Year,
                               fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c("#A2A4A2",
                                   "#717F5B",
                                   "#D4A634",
                                   "goldenrod4")) +
    scale_fill_manual(values = c("#A2A4A2",
                                 "#717F5B",
                                 "#D4A634",
                                 "goldenrod4")) +
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
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill=NA),
          legend.position = c(.85, .1),
          legend.box = "horizontal") +
    facet_grid(~Day) +
    labs(title = "Nest Survival",
         color = "Year",
         x = "Days Grazed During Nesting",
         y = "Daily Survival Rate"))


pDoD.pred <- covariate.predictions(plotdata,
                                   data = data.frame(Forb = mean(Forbvalues),
                                                     grazed = mean(Grazedvalues),
                                                     pDoD = pDoDvalues),
                                   indices = c(2, 16, 28,
                                               77, 91, 103,
                                               152, 166, 178,
                                               227, 241, 253))

D1Y2021 <- which(pDoD.pred$estimates$par.index == 2)
D15Y2021 <- which(pDoD.pred$estimates$par.index == 16)
D27Y2021 <- which(pDoD.pred$estimates$par.index == 28)
D1Y2022 <- which(pDoD.pred$estimates$par.index == 77)
D15Y2022 <- which(pDoD.pred$estimates$par.index == 91)
D27Y2022 <- which(pDoD.pred$estimates$par.index == 103)
D1Y2023 <- which(pDoD.pred$estimates$par.index == 152)
D15Y2023 <- which(pDoD.pred$estimates$par.index == 166)
D27Y2023 <- which(pDoD.pred$estimates$par.index == 178)
D1Y2024 <- which(pDoD.pred$estimates$par.index == 227)
D15Y2024 <- which(pDoD.pred$estimates$par.index == 241)
D27Y2024 <- which(pDoD.pred$estimates$par.index == 253)

pDoD.pred$estimates$Year <- NA
pDoD.pred$estimates$Year[D1Y2021] <- "2021"
pDoD.pred$estimates$Year[D15Y2021] <- "2021"
pDoD.pred$estimates$Year[D27Y2021] <- "2021"
pDoD.pred$estimates$Year[D1Y2022] <- "2022"
pDoD.pred$estimates$Year[D15Y2022] <- "2022"
pDoD.pred$estimates$Year[D27Y2022] <- "2022"
pDoD.pred$estimates$Year[D1Y2023] <- "2023"
pDoD.pred$estimates$Year[D15Y2023] <- "2023"
pDoD.pred$estimates$Year[D27Y2023] <- "2023"
pDoD.pred$estimates$Year[D1Y2024] <- "2024"
pDoD.pred$estimates$Year[D15Y2024] <- "2024"
pDoD.pred$estimates$Year[D27Y2024] <- "2024"

pDoD.pred$estimates$Day <- NA
pDoD.pred$estimates$Day[D1Y2021] <- "Day1"
pDoD.pred$estimates$Day[D15Y2021] <- "Day15"
pDoD.pred$estimates$Day[D27Y2021] <- "Day27"
pDoD.pred$estimates$Day[D1Y2022] <- "Day1"
pDoD.pred$estimates$Day[D15Y2022] <- "Day15"
pDoD.pred$estimates$Day[D27Y2022] <- "Day27"
pDoD.pred$estimates$Day[D1Y2023] <- "Day1"
pDoD.pred$estimates$Day[D15Y2023] <- "Day15"
pDoD.pred$estimates$Day[D27Y2023] <- "Day27"
pDoD.pred$estimates$Day[D1Y2024] <- "Day1"
pDoD.pred$estimates$Day[D15Y2024] <- "Day15"
pDoD.pred$estimates$Day[D27Y2024] <- "Day27"

(GADWpDoD.plot <- ggplot(transform(pDoD.pred$estimates,
                                   Year = factor(Year, levels = c("2021", "2022", "2023", "2024"))), 
                         aes(x = pDoD, 
                             y = estimate,
                             groups = Year,
                             fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c("#A2A4A2",
                                   "#717F5B",
                                   "#D4A634",
                                   "goldenrod4")) +
    scale_fill_manual(values = c("#A2A4A2",
                                 "#717F5B",
                                 "#D4A634",
                                 "goldenrod4")) +
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
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill=NA),
          legend.position = c(.85, .1),
          legend.box = "horizontal") +
    facet_grid(~Day) +
    labs(title = "Nest Survival",
         color = "Year",
         x = "Previous Year Degree of Disappearance",
         y = "Daily Survival Rate"))


ggsave(GADW.plot,
       filename = "outputs/figs/betaGADW.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(GADWage.plot,
       filename = "outputs/figs/ageGADW.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(GADWforb.plot,
       filename = "outputs/figs/forbGADW.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(GADWgrazed.plot,
       filename = "outputs/figs/grazedGADW.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(GADWpDoD.plot,
       filename = "outputs/figs/pDoDGADW.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)


# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute "rm(list = ls(all = TRUE))" - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all = TRUE))
# Then, execute "cleanup(ask = FALSE)" to delete orphaned output
#  files from MARK. Execute "?cleanup" to learn more
cleanup(ask = FALSE)

 