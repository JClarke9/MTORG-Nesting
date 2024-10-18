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


RWBL.surv <- filter(nest, 
                    Spec == "RWBL" & Stage != "Laying")                                         # select out only RWBL nest

test <- filter(RWBL.surv,
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

MISSING <- is.na(RWBL.surv$AgeFound)

sum(MISSING)

RWBL.surv <- subset(RWBL.surv, 
                    subset = !MISSING)

RWBL.surv$Year <- factor(RWBL.surv$Year,
                         levels = c("2021", "2022", "2023", "2024"))

RWBL.surv$Nestling <- factor(RWBL.surv$Nestling,
                             level = c("0", "1"))

str(RWBL.surv)


# Creating stage variable -------------------------------------------------


x <- create.stage.var(RWBL.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(RWBL.surv$LastChecked)), 
                      12)

RWBL.surv <- bind_cols(RWBL.surv, x)

rm(list = ls()[!ls() %in% c("RWBL.surv")])


# Daily survival rate models ----------------------------------------------


RWBL.pr <- process.data(RWBL.surv,
                        nocc = max(RWBL.surv$LastChecked),
                        groups = c("Year",
                                   "Nestling"),
                        model = "Nest")



# Temporal candidate model set
RWBL1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL1.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL1.results <- RWBL1.run()
RWBL1.results

# both of the quadratic effects overlap zero
coef(RWBL1.results$S.quad)
confint(RWBL1.results$S.quad, level = 0.85)

coef(RWBL1.results$S.time)
confint(RWBL1.results$S.time, level = 0.85)


# Nest stage/age candidate model set
RWBL2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 2. DSR varies with nest stage
  S.stage = list(formula = ~1 + Time + Nestling)
  
  # 3. DSR varies with nest age
  S.age = list(formula = ~1 + Time + NestAge)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL2.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL2.results <- RWBL2.run()
RWBL2.results

coef(RWBL2.results$S.stage)
confint(RWBL2.results$S.stage, level = 0.85)



# Biological candidate model set
RWBL3.run <- function()
{
  # 1. DSR varies with year
  S.stage = list(formula = ~1 + Time + Nestling)
  
  # 2. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Time + Nestling + BHCONum)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL3.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL3.results <- RWBL3.run()
RWBL3.results

coef(RWBL3.results$S.bhcon)
confint(RWBL3.results$S.bhcon, level = 0.85)



# Grazing candidate model set
RWBL4.run <- function()
{
  # 1. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Time + Nestling + BHCONum)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Time + Nestling + BHCONum + grazed)
  
  # 3. DSR varies with the previous Year + Nestlings grazing intensity
  S.pDoD = list(formula = ~1 + Time + Nestling + BHCONum + pDoD)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL4.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL4.results <- RWBL4.run()
RWBL4.results

coef(RWBL4.results$S.bhcon)
confint(RWBL4.results$S.bhcon, level = 0.85)


# Vegetation candidate model set
RWBL5.run <- function()
{
  # 1. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Time + Nestling + BHCONum)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula = ~1 + Time + Nestling + BHCONum + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Time + Nestling + BHCONum + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG) - NA values
  # S.lit = list(formula = ~1 + Time + Nestling + BHCONum + Litter)
  
  # 5. DSR varies with Bare - NA values
  # S.bare = list(formula = ~1 + Time + Nestling + BHCONum + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula = ~1 + Time + Nestling + BHCONum + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula = ~1 + Time + Nestling + BHCONum + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula = ~1 + Time + Nestling + BHCONum + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR) - NA values
  # S.litdep = list(formula = ~1 + Time + Nestling + BHCONum + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula = ~1 + Time + Nestling + BHCONum + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula = ~1 + Time + Nestling + BHCONum + VOR)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL5.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL5.results <- RWBL5.run()
RWBL5.results

coef(RWBL5.results$S.bhcon)
confint(RWBL5.results$S.bhcon, level = 0.85)



RWBL5.results$S.bhcon$results$real |> 
  summarize(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))

(RWBL.real <- as.data.frame(RWBL5.results$S.bhcon$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Stage = case_when(
      grepl("20210", Group) ~ "Incubating",
      grepl("20230", Group) ~ "Incubating",
      grepl("20240", Group) ~ "Incubating",
      grepl("20211", Group) ~ "Nestling",
      grepl("20221", Group) ~ "Nestling",
      grepl("20231", Group) ~ "Nestling",
      grepl("20241", Group) ~ "Nestling")) |> 
    select(Stage, estimate, se, lcl, ucl) |> 
    group_by(Stage) |> 
    summarize(mean = mean(estimate)))


# Plotting beta coefficients ----------------------------------------------


RWBL.beta <- coef(RWBL4.results$S.bhcon) |>
  cbind(confint(RWBL4.results$S.bhcon, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

RWBL.beta$Variable <- gsub("S:", "", RWBL.beta$Variable)

str(RWBL.beta)

(RWBL.plot <- ggplot(RWBL.beta[4,], 
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
    labs(title = "Red-winged Blackbird",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


RWBL.ddl <- make.design.data(RWBL.pr) |> 
  as.data.frame()

plotdata <- RWBL5.results$S.bhcon


max(RWBL.ddl$S.Age)

filter(RWBL.ddl, 
       S.Nestling == 0 & S.time == 1 & S.Year == 2021| 
         S.Nestling == 0 & S.time == 18  & S.Year == 2021| 
         S.Nestling == 0 & S.time == 36  & S.Year == 2021| 
         S.Nestling == 0 & S.time == 54  & S.Year == 2021| 
         S.Nestling == 0 & S.time == 75 & S.Year == 2021|
         S.Nestling == 1 & S.time == 1 & S.Year == 2021| 
         S.Nestling == 1 & S.time == 18  & S.Year == 2021| 
         S.Nestling == 1 & S.time == 36  & S.Year == 2021| 
         S.Nestling == 1 & S.time == 54  & S.Year == 2021| 
         S.Nestling == 1 & S.time == 75 & S.Year == 2021)

time.pred <- covariate.predictions(plotdata,
                                   indices = c(1:75, 305:379))

inc <- which(time.pred$estimates$par.index == 1)
nst <- which(time.pred$estimates$par.index == 305)

time.pred$estimates$Stage <- NA
time.pred$estimates$Stage[inc] <- "Incubating"
time.pred$estimates$Stage[nst] <- "Nestling"

time.pred$estimates <- fill(time.pred$estimates, Stage, .direction = "down")

time.pred$estimates$Day <- 1:75

(RWBLtime.plot <- ggplot(time.pred$estimates, 
                         aes(x = Day, 
                             y = estimate)) +
    geom_line(linewidth = 1.5,
              aes(color = "model.index")) +
    scale_colour_manual(values = c('#D4A634')) +
    scale_fill_manual(values = c('#D4A634')) +
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
          axis.text.y = element_text(size=12, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size=12, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = "none") +
    facet_grid(~Stage) +
    labs(title = "Nest Survival",
         x = "Julian Day",
         y = "Daily Survival Rate"))


stage.pred <- covariate.predictions(plotdata,
                                    indices = c(1, 18, 36, 54, 75,
                                                305, 322, 340, 358, 379))

inc <- which(stage.pred$estimates$par.index == 1)
nst <- which(stage.pred$estimates$par.index == 305)

stage.pred$estimates$Stage <- NA
stage.pred$estimates$Stage[inc] <- "Incubating"
stage.pred$estimates$Stage[nst] <- "Nestling"

stage.pred$estimates <- fill(stage.pred$estimates, Stage, .direction = "down")

stage.pred$estimates$Day <- c(1, 18, 36, 54, 75)

(RWBLstage.plot <- ggplot(transform(stage.pred$estimates,
                                    Day = factor(Day, levels = c("1", "18", "36",
                                                                 "54", "75"))), 
                          aes(x = Stage, 
                              y = estimate,
                              groups = Day,
                              fill = Day)) +
    geom_point(size = 4,
               aes(color = Day)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#A2A4A2',
                                   '#717F5B',
                                   '#D4A634',
                                   'red',
                                   'black')) +
    scale_fill_manual(values = c('#A2A4A2',
                                 '#717F5B',
                                 '#D4A634',
                                 'red',
                                 'black')) +
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
         x = "Stage",
         y = "Daily Survival Rate"))


ggsave(RWBL.plot,
       filename = "outputs/figs/betaRWBL.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(RWBLtime.plot,
       filename = "outputs/figs/timeRWBL.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(RWBLstage.plot,
       filename = "outputs/figs/stageRWBL.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)


# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all = TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
#  files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)
