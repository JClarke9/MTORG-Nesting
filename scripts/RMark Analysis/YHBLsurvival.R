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


YHBL.surv <- filter(nest, 
                    Spec == "YHBL" & Stage != "Laying")                                         # select out only YHBL nest

test <- filter(YHBL.surv,
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

MISSING <- is.na(YHBL.surv$AgeFound)

sum(MISSING)

YHBL.surv <- subset(YHBL.surv, 
                    subset = !MISSING)

YHBL.surv$Year <- factor(YHBL.surv$Year,
                         levels = c("2021", "2022", "2023", "2024"))

YHBL.surv$Nestling <- factor(YHBL.surv$Nestling,
                             level = c("0", "1"))

str(YHBL.surv)


# Creating stage variable -------------------------------------------------


x <- create.stage.var(YHBL.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(YHBL.surv$LastChecked)), 
                      12)

YHBL.surv <- bind_cols(YHBL.surv, x)

rm(list = ls()[!ls() %in% c("YHBL.surv")])


# Daily survival rate models ----------------------------------------------


YHBL.pr <- process.data(YHBL.surv,
                        nocc = max(YHBL.surv$LastChecked),
                        groups = c("Year",
                                   "Nestling"),
                        model = "Nest")

# Temporal candidate model set
# I didn't include year because most years we didn't have many nests
YHBL1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  YHBL.model.list = create.model.list("Nest")
  YHBL1.results = mark.wrapper(YHBL.model.list,
                               data = YHBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
YHBL1.results <- YHBL1.run()
YHBL1.results

coef(YHBL1.results$S.time)
confint(YHBL1.results$S.time, level = 0.85)

coef(YHBL1.results$S.quad)
confint(YHBL1.results$S.quad, level = 0.85)



# Biological candidate model set
YHBL2.run <- function()
{
  # 3. DSR varies with quadratic effect of date
  S.time = list(formula = ~1 + Time)
  
  # 2. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Time + BHCONum)
  
  # 2. DSR varies with BHCO number
  S.bhcop = list(formula = ~1 + Time + BHCOPres)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Time + NestAge)
  
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Time + Nestling)
  
  YHBL.model.list = create.model.list("Nest")
  YHBL2.results = mark.wrapper(YHBL.model.list,
                               data = YHBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
YHBL2.results <- YHBL2.run()
YHBL2.results

coef(YHBL2.results$S.stage)
confint(YHBL2.results$S.stage, level = 0.85)


# Grazing candidate model set
YHBL3.run <- function()
{
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Time + I(Time^2) + Nestling)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Time + I(Time^2) + Nestling + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Time + I(Time^2) + Nestling + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pDoD = list(formula = ~1 + Time + I(Time^2) + Nestling + pDoD)
  
  YHBL.model.list = create.model.list("Nest")
  YHBL3.results = mark.wrapper(YHBL.model.list,
                               data = YHBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
YHBL3.results <- YHBL3.run()
YHBL3.results

coef(YHBL2.results$S.stage)
confint(YHBL2.results$S.stage, level = 0.85)


# Vegetation candidate model set
YHBL4.run <- function()
{
  # 2. DSR varies with the number of days a nest experienced grazing
  S.stage = list(formula = ~1 + Time + I(Time^2) + Nestling)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Time + I(Time^2) + Nestling + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Time + I(Time^2) + Nestling + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  # S.lit = list(formula =  ~1 + Time + I(Time^2) + Nestling + Litter)
  
  # 5. DSR varies with Bare
  #S.bare = list(formula =  ~1 + Time + I(Time^2) + Nestling + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Time + I(Time^2) + Nestling + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Time + I(Time^2) + Nestling + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Time + I(Time^2) + Nestling + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  #S.litdep = list(formula =  ~1 + Time + I(Time^2) + Nestling + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Time + I(Time^2) + Nestling + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Time + I(Time^2) + Nestling + VOR)
  
  YHBL.model.list = create.model.list("Nest")
  YHBL4.results = mark.wrapper(YHBL.model.list,
                               data = YHBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
YHBL4.results <- YHBL4.run()
YHBL4.results

coef(YHBL4.results$S.stage)
confint(YHBL4.results$S.stage, level = 0.85)

YHBL4.results$S.stage$results$real |> 
  summarize(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))

(YHBL.real <- as.data.frame(YHBL4.results$S.stage$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Stage = case_when(
      grepl("20210", Group) ~ "Incubating",
      grepl("20220", Group) ~ "Incubating",
      grepl("20230", Group) ~ "Incubating",
      grepl("20240", Group) ~ "Incubating",
      grepl("20211", Group) ~ "Nestling",
      grepl("20221", Group) ~ "Nestling",
      grepl("20231", Group) ~ "Nestling",
      grepl("20241", Group) ~ "Nestling")) |> 
    group_by(Stage) |> 
    summarize(estimate = mean(estimate), 
              se = mean(se), 
              lcl = mean(lcl), 
              ucl = mean(ucl)))


# Plotting beta coefficients ----------------------------------------------


YHBL.beta <- coef(YHBL4.results$S.grazed) |>
  cbind(confint(YHBL4.results$S.grazed, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

YHBL.beta$Variable <- gsub("S:", "", YHBL.beta$Variable)

str(YHBL.beta)

(YHBL.plot <- ggplot(YHBL.beta[4:5,], 
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


YHBL.ddl <- make.design.data(YHBL.pr) |> 
  as.data.frame()

max(YHBL.ddl$S.Age)

filter(YHBL.ddl, 
       S.Nestling == 0 & S.time == 1 & S.Year == 2021| 
         S.Nestling == 0 & S.time == 18  & S.Year == 2021| 
         S.Nestling == 0 & S.time == 36  & S.Year == 2021| 
         S.Nestling == 0 & S.time == 54  & S.Year == 2021| 
         S.Nestling == 0 & S.time == 72 & S.Year == 2021|
         S.Nestling == 1 & S.time == 1 & S.Year == 2021| 
         S.Nestling == 1 & S.time == 18  & S.Year == 2021| 
         S.Nestling == 1 & S.time == 36  & S.Year == 2021| 
         S.Nestling == 1 & S.time == 54  & S.Year == 2021| 
         S.Nestling == 1 & S.time == 72 & S.Year == 2021)

plotdata <- YHBL4.results$S.grazed

grazed.values <- seq(from = min(YHBL.surv$grazed), 
                     to = max(YHBL.surv$grazed), 
                     length = 100)


time.pred <- covariate.predictions(plotdata,
                                   data = data.frame(grazed = mean(grazed.values)),
                                   indices = c(1:73, 220:292))

inc <- which(time.pred$estimates$par.index == 1)
nst <- which(time.pred$estimates$par.index == 220)

time.pred$estimates$Stage <- NA
time.pred$estimates$Stage[inc] <- "Incubating"
time.pred$estimates$Stage[nst] <- "Nestling"

time.pred$estimates <- fill(time.pred$estimates, Stage, .direction = "down")

time.pred$estimates$Day <- 1:73

(YHBLtime.plot <- ggplot(time.pred$estimates, 
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
    labs(title = "Red-winged Blackbird",
         x = "Julian Day",
         y = "Daily Survival Rate"))


stage.pred <- covariate.predictions(plotdata,
                                    data = data.frame(grazed = mean(grazed.values)),
                                    indices = c(1, 18, 36, 54, 72,
                                                220, 237, 255, 273, 291))

inc <- which(stage.pred$estimates$par.index == 1)
nst <- which(stage.pred$estimates$par.index == 220)

stage.pred$estimates$Stage <- NA
stage.pred$estimates$Stage[inc] <- "Incubating"
stage.pred$estimates$Stage[nst] <- "Nestling"

stage.pred$estimates <- fill(stage.pred$estimates, Stage, .direction = "down")

stage.pred$estimates$Day <- c(1, 18, 36, 54, 72)

(YHBLstage.plot <- ggplot(transform(stage.pred$estimates,
                                    Day = factor(Day, levels = c("1", "18", "36",
                                                                 "54", "72"))), 
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
    labs(title = "Red-winged Blackbird",
         x = "Stage",
         y = "Daily Survival Rate"))


grazed.pred <- covariate.predictions(plotdata,
                                     data = data.frame(grazed = grazed.values),
                                     indices = c(1, 18, 36, 54, 72,
                                                 220, 237, 255, 273, 291))

inc1 <- which(grazed.pred$estimates$par.index == 1)
inc18 <- which(grazed.pred$estimates$par.index == 18)
inc36 <- which(grazed.pred$estimates$par.index == 36)
inc54 <- which(grazed.pred$estimates$par.index == 54)
inc72 <- which(grazed.pred$estimates$par.index == 72)
nst1 <- which(grazed.pred$estimates$par.index == 220)
nst18 <- which(grazed.pred$estimates$par.index == 237)
nst36 <- which(grazed.pred$estimates$par.index == 255)
nst54 <- which(grazed.pred$estimates$par.index == 273)
nst72 <- which(grazed.pred$estimates$par.index == 291)

grazed.pred$estimates$Group <- NA
grazed.pred$estimates$Group[inc1] <- "Incubating1"
grazed.pred$estimates$Group[inc18] <- "Incubating18"
grazed.pred$estimates$Group[inc36] <- "Incubating36"
grazed.pred$estimates$Group[inc54] <- "Incubating54"
grazed.pred$estimates$Group[inc72] <- "Incubating72"
grazed.pred$estimates$Group[nst1] <- "Nestling1"
grazed.pred$estimates$Group[nst18] <- "Nestling18"
grazed.pred$estimates$Group[nst36] <- "Nestling36"
grazed.pred$estimates$Group[nst54] <- "Nestling54"
grazed.pred$estimates$Group[nst72] <- "Nestling72"

grazed.pred$estimates$Stage <- NA
grazed.pred$estimates$Stage[inc1] <- "Incubating"
grazed.pred$estimates$Stage[inc18] <- "Incubating"
grazed.pred$estimates$Stage[inc36] <- "Incubating"
grazed.pred$estimates$Stage[inc54] <- "Incubating"
grazed.pred$estimates$Stage[inc72] <- "Incubating"
grazed.pred$estimates$Stage[nst1] <- "Nestling"
grazed.pred$estimates$Stage[nst18] <- "Nestling"
grazed.pred$estimates$Stage[nst36] <- "Nestling"
grazed.pred$estimates$Stage[nst54] <- "Nestling"
grazed.pred$estimates$Stage[nst72] <- "Nestling"

grazed.pred$estimates$Day <- NA
grazed.pred$estimates$Day[inc1] <- "1"
grazed.pred$estimates$Day[inc18] <- "18"
grazed.pred$estimates$Day[inc36] <- "36"
grazed.pred$estimates$Day[inc54] <- "54"
grazed.pred$estimates$Day[inc72] <- "72"
grazed.pred$estimates$Day[nst1] <- "1"
grazed.pred$estimates$Day[nst18] <- "18"
grazed.pred$estimates$Day[nst36] <- "36"
grazed.pred$estimates$Day[nst54] <- "54"
grazed.pred$estimates$Day[nst72] <- "72"

(YHBLgrazed.plot <- ggplot(transform(grazed.pred$estimates,
                                     Day = factor(Day, levels = c("1", "18", "36",
                                                                  "54", "72"))), 
                           aes(x = covdata, 
                               y = estimate,
                               groups = Group)) +
    geom_line(linewidth = 1.5,
              aes(color = Day)) +
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
          legend.position = c(.95, .15),
          legend.box = "horizontal") +
    facet_grid(~Stage) +
    labs(title = "Red-winged Blackbird",
         color = "Day",
         x = "Days Grazed",
         y = "Daily Survival Rate"))


ggsave(YHBL.plot,
       filename = "outputs/figs/betaYHBL.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(YHBLtime.plot,
       filename = "outputs/figs/timeYHBL.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(YHBLstage.plot,
       filename = "outputs/figs/stageYHBL.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(YHBLgrazed.plot,
       filename = "outputs/figs/grazedYHBL.png",
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
