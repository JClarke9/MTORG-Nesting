# Loading libraries -------------------------------------------------------


library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)
library(MuMIn)
library(cowplot)
source("scripts/Functions/RMark_Stage_Code.R")

windowsFonts(my_font = windowsFont("Gandhi Sans"))


# Data import -------------------------------------------------------------


nest <- read.csv("working/RMarknesting.csv")


# Subsetting data ---------------------------------------------------------


MODO.surv <- filter(nest, 
                    Spec == "MODO" & Stage != "Laying")                                         # select out only MODO nest

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
                         levels = c("2021", "2022", "2023", "2024"))

MODO.surv$Nestling <- factor(MODO.surv$Nestling,
                             level = c("0", "1"))

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
                        nocc = max(MODO.surv$LastChecked),
                        groups = c("Year",
                                   "Nestling"),
                        model = "Nest")



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

# only 2023 doesn't overlap zero
coef(MODO1.results$S.year)
confint(MODO1.results$S.year, level = 0.85)



# Nest stage/age candidate model set
MODO2.run <- function()
{
  # 1. DSR varies with year and quadratic
  S.yearQuad = list(formula = ~1 + Year + Time + I(Time^2))
  
  # 2. DSR varies with nest stage
  S.stage = list(formula = ~1 + Year + Time + I(Time^2) + Nestling)
  
  # 3. DSR varies with nest age
  S.age = list(formula = ~1 + Year + Time + I(Time^2) + NestAge)
  
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
MODO4.run <- function()
{
  # 1. DSR varies with nest stage
  S.stage = list(formula = ~1 + Year + Time + I(Time^2) + Nestling)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + grazed)
  
  # 3. DSR varies with the previous Year + Nestlings grazing intensity
  S.pDoD = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + pDoD)
  
  MODO.model.list = create.model.list("Nest")
  MODO4.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MODO4.results <- MODO4.run()
MODO4.results

coef(MODO4.results$S.stage)
confint(MODO4.results$S.stage, level = 0.85)



# Vegetation candidate model set
MODO5.run <- function()
{
  # 1. DSR varies with nest stage
  S.stage = list(formula = ~1 + Year + Time + I(Time^2) + Nestling)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula = ~1 + Year + Time + I(Time^2) + Nestling + VOR)
  
  MODO.model.list = create.model.list("Nest")
  MODO5.results = mark.wrapper(MODO.model.list,
                               data = MODO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MODO5.results <- MODO5.run()
MODO5.results

coef(MODO5.results$S.kbg)
confint(MODO5.results$S.kbg, level = 0.85)



MODO5.results$S.kbg$results$real |> 
  summarize(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))

(MODO.real <- as.data.frame(MODO5.results$S.kbg$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023",
      grepl("2024", Group) ~ "2024"),
      Stage = case_when(
        grepl("20210", Group) ~ "Incubating",
        grepl("20230", Group) ~ "Incubating",
        grepl("20240", Group) ~ "Incubating",
        grepl("20211", Group) ~ "Nestling",
        grepl("20221", Group) ~ "Nestling",
        grepl("20231", Group) ~ "Nestling",
        grepl("20241", Group) ~ "Nestling")) |> 
    select(Year, Stage, estimate, se, lcl, ucl))

(MODO.year <- MODO.real |> 
    group_by(Year) |> 
    summarize(mean = mean(estimate)))

(MODO.stage <- MODO.real |> 
    group_by(Stage) |> 
    summarize(mean = mean(estimate)))


# Plotting beta coefficients ----------------------------------------------


MODO.beta <- coef(MODO5.results$S.kbg) |>
  cbind(confint(MODO5.results$S.kbg, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

MODO.beta$Variable <- gsub("S:", "", MODO.beta$Variable)

str(MODO.beta)

(MODO.plot <- ggplot(MODO.beta[c(8),], 
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
    labs(title = "Mourning Dove",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


MODO.ddl <- make.design.data(MODO.pr) |> 
  as.data.frame()

plotdata <- MODO5.results$S.kbg

kbg.values <- seq(from = min(MODO.surv$KBG), 
                  to = max(MODO.surv$KBG), 
                  length = 100)


max(MODO.ddl$S.Age)

filter(MODO.ddl, S.Nestling == 0 & S.time %in% c(1:74) & S.Year == 2021)
filter(MODO.ddl, S.Nestling == 0 & S.time %in% c(1:74) & S.Year == 2022)
filter(MODO.ddl, S.Nestling == 0 & S.time %in% c(1:74) & S.Year == 2023)
filter(MODO.ddl, S.Nestling == 0 & S.time %in% c(1:74) & S.Year == 2024)
filter(MODO.ddl, S.Nestling == 1 & S.time %in% c(1:74) & S.Year == 2021)
filter(MODO.ddl, S.Nestling == 1 & S.time %in% c(1:74) & S.Year == 2022)
filter(MODO.ddl, S.Nestling == 1 & S.time %in% c(1:74) & S.Year == 2023)
filter(MODO.ddl, S.Nestling == 1 & S.time %in% c(1:74) & S.Year == 2024)

time.pred <- covariate.predictions(plotdata,
                                   data = data.frame(KBG = mean(kbg.values)),
                                   indices = c(1:74,
                                               76:149,
                                               151:224,
                                               226:299,
                                               301:374,
                                               376:449,
                                               451:524,
                                               526:599))

inc2021 <- which(time.pred$estimates$par.index %in% c(1:74))
inc2022 <- which(time.pred$estimates$par.index %in% c(76:149))
inc2023 <- which(time.pred$estimates$par.index %in% c(151:224))
inc2024 <- which(time.pred$estimates$par.index %in% c(226:299))
nst2021 <- which(time.pred$estimates$par.index %in% c(301:374))
nst2022 <- which(time.pred$estimates$par.index %in% c(376:449))
nst2023 <- which(time.pred$estimates$par.index %in% c(451:524))
nst2024 <- which(time.pred$estimates$par.index %in% c(526:599))

time.pred$estimates$Year <- NA
time.pred$estimates$Year[inc2021] <- "2021"
time.pred$estimates$Year[inc2022] <- "2022"
time.pred$estimates$Year[inc2023] <- "2023"
time.pred$estimates$Year[inc2024] <- "2024"
time.pred$estimates$Year[nst2021] <- "2021"
time.pred$estimates$Year[nst2022] <- "2022"
time.pred$estimates$Year[nst2023] <- "2023"
time.pred$estimates$Year[nst2024] <- "2024"

time.pred$estimates$Stage <- NA
time.pred$estimates$Stage[inc2021] <- "Incubating"
time.pred$estimates$Stage[inc2022] <- "Incubating"
time.pred$estimates$Stage[inc2023] <- "Incubating"
time.pred$estimates$Stage[inc2024] <- "Incubating"
time.pred$estimates$Stage[nst2021] <- "Nestling"
time.pred$estimates$Stage[nst2022] <- "Nestling"
time.pred$estimates$Stage[nst2023] <- "Nestling"
time.pred$estimates$Stage[nst2024] <- "Nestling"

time.pred$estimates$Day <- c(1:74)

(MODOtime.plot <- ggplot(transform(time.pred$estimates,
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
    facet_grid(~Stage) +
    labs(title = "Nest Survival",
         x = "time",
         y = "Daily Survival Rate"))


stage.pred <- covariate.predictions(plotdata,
                                    data = data.frame(KBG = mean(kbg.values)),
                                    indices = c(1, 18, 36, 54, 74,
                                                76, 93, 111, 129, 149,
                                                151, 168, 186, 204, 224,
                                                226, 243, 261, 279, 229,
                                                301, 318, 336, 354, 374,
                                                376, 393, 411, 429, 449,
                                                451, 468, 486, 504, 524,
                                                526, 543, 561, 579, 599))

inc2021 <- which(stage.pred$estimates$par.index == 1)
inc2022 <- which(stage.pred$estimates$par.index == 76)
inc2023 <- which(stage.pred$estimates$par.index == 151)
inc2024 <- which(stage.pred$estimates$par.index == 226)
nst2021 <- which(stage.pred$estimates$par.index == 301)
nst2022 <- which(stage.pred$estimates$par.index == 376)
nst2023 <- which(stage.pred$estimates$par.index == 451)
nst2024 <- which(stage.pred$estimates$par.index == 526)

stage.pred$estimates$Year <- NA
stage.pred$estimates$Year[inc2021] <- "2021"
stage.pred$estimates$Year[inc2022] <- "2022"
stage.pred$estimates$Year[inc2023] <- "2023"
stage.pred$estimates$Year[inc2024] <- "2024"
stage.pred$estimates$Year[nst2021] <- "2021"
stage.pred$estimates$Year[nst2022] <- "2022"
stage.pred$estimates$Year[nst2023] <- "2023"
stage.pred$estimates$Year[nst2024] <- "2024"

stage.pred$estimates$Stage <- NA
stage.pred$estimates$Stage[inc2021] <- "Incubating"
stage.pred$estimates$Stage[inc2022] <- "Incubating"
stage.pred$estimates$Stage[inc2023] <- "Incubating"
stage.pred$estimates$Stage[inc2024] <- "Incubating"
stage.pred$estimates$Stage[nst2021] <- "Nestling"
stage.pred$estimates$Stage[nst2022] <- "Nestling"
stage.pred$estimates$Stage[nst2023] <- "Nestling"
stage.pred$estimates$Stage[nst2024] <- "Nestling"

stage.pred$estimates <- fill(stage.pred$estimates, Year, .direction = "down")
stage.pred$estimates <- fill(stage.pred$estimates, Stage, .direction = "down")

stage.pred$estimates$Day <- c(1, 18, 36, 54, 74)

(MODOstage.plot <- ggplot(transform(stage.pred$estimates,
                                    Day = factor(Day, levels = c("1", "18", "36",
                                                                 "54", "74"))), 
                          aes(x = Stage, 
                              y = estimate,
                              groups = Year,
                              fill = Year)) +
    geom_point(size = 4,
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
    facet_wrap(~Day) +
    labs(title = "Nest Survival",
         x = "Stage",
         y = "Daily Survival Rate"))


kbg.pred <- covariate.predictions(plotdata,
                                  data = data.frame(KBG = kbg.values),
                                  indices = c(1, 18, 36, 54, 74,
                                              76, 93, 111, 129, 149,
                                              151, 168, 186, 204, 224,
                                              226, 243, 261, 279, 229,
                                              301, 318, 336, 354, 374,
                                              376, 393, 411, 429, 449,
                                              451, 468, 486, 504, 524,
                                              526, 543, 561, 579, 599))

inc2021 <- which(kbg.pred$estimates$par.index == 1)
inc2022 <- which(kbg.pred$estimates$par.index == 76)
inc2023 <- which(kbg.pred$estimates$par.index == 151)
inc2024 <- which(kbg.pred$estimates$par.index == 226)
nst2021 <- which(kbg.pred$estimates$par.index == 301)
nst2022 <- which(kbg.pred$estimates$par.index == 376)
nst2023 <- which(kbg.pred$estimates$par.index == 451)
nst2024 <- which(kbg.pred$estimates$par.index == 526)

kbg.pred$estimates$Year <- NA
kbg.pred$estimates$Year[inc2021] <- "2021"
kbg.pred$estimates$Year[inc2022] <- "2022"
kbg.pred$estimates$Year[inc2023] <- "2023"
kbg.pred$estimates$Year[inc2024] <- "2024"
kbg.pred$estimates$Year[nst2021] <- "2021"
kbg.pred$estimates$Year[nst2022] <- "2022"
kbg.pred$estimates$Year[nst2023] <- "2023"
kbg.pred$estimates$Year[nst2024] <- "2024"

kbg.pred$estimates$Stage <- NA
kbg.pred$estimates$Stage[inc2021] <- "Incubating"
kbg.pred$estimates$Stage[inc2022] <- "Incubating"
kbg.pred$estimates$Stage[inc2023] <- "Incubating"
kbg.pred$estimates$Stage[inc2024] <- "Incubating"
kbg.pred$estimates$Stage[nst2021] <- "Nestling"
kbg.pred$estimates$Stage[nst2022] <- "Nestling"
kbg.pred$estimates$Stage[nst2023] <- "Nestling"
kbg.pred$estimates$Stage[nst2024] <- "Nestling"

kbg.pred$estimates <- fill(kbg.pred$estimates, Year, .direction = "down")
kbg.pred$estimates <- fill(kbg.pred$estimates, Stage, .direction = "down")

kbg.pred$estimates$Day <- c(1, 18, 36, 54, 74)

(MODOkbg.plot <- ggplot(transform(kbg.pred$estimates,
                                  Year = factor(Year, levels = c("2021", "2022", "2023", "2024"))), 
                        aes(x = covdata, 
                            y = estimate,
                            groups = Year,
                            fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
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
          legend.position = c(.95, .15),
          legend.box = "horizontal") +
    facet_grid(~Stage + Day) +
    labs(title = "Nest Survival",
         color = "Day",
         x = "Kentucky Bluegrass (Percent Cover)",
         y = "Daily Survival Rate"))


ggsave(MODO.plot,
       filename = "outputs/figs/betaMODO.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(MODOtime.plot,
       filename = "outputs/figs/timeMODO.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(MODOstage.plot,
       filename = "outputs/figs/stageMODO.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(MODOkbg.plot,
       filename = "outputs/figs/kbgMODO.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)


# If you want to clean up the mark*.inp, .vcv, .res and .out
# and .tmp files created by RMark in the working directory,
# execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all = TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
# files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)

