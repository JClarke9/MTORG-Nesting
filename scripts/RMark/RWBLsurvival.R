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
                         levels = c("2021", "2022", "2023"))

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

coef(RWBL1.results$S.quad)
confint(RWBL1.results$S.quad, level = 0.85)



# Biological candidate model set
RWBL2.run <- function()
{
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 2. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Time + I(Time^2) + BHCONum)
  
  # 2. DSR varies with BHCO number
  S.bhcop = list(formula = ~1 + Time + I(Time^2) + BHCOPres)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Time + I(Time^2) + NestAge)
  
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Time + I(Time^2) + Nestling)
  
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


# Grazing candidate model set
RWBL3.run <- function()
{
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Time + I(Time^2) + Nestling)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Time + I(Time^2) + Nestling + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Time + I(Time^2) + Nestling + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pTreat = list(formula = ~1 + Time + I(Time^2) + Nestling + pTreat)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL3.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL3.results <- RWBL3.run()
RWBL3.results

coef(RWBL3.results$S.grazed)
coef(RWBL3.results$S.grazep)

confint(RWBL3.results$S.grazed, level = 0.85)
confint(RWBL3.results$S.grazep, level = 0.85)


# Vegetation candidate model set
RWBL4.run <- function()
{
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Time + I(Time^2) + Nestling + grazed)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Time + I(Time^2) + Nestling + grazed + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + Litter)
  
  # 5. DSR varies with Bare
  #S.bare = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  #S.litdep = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Time + I(Time^2) + Nestling + grazed + VOR)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL4.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL4.results <- RWBL4.run()
RWBL4.results

coef(RWBL4.results$S.grazed)
confint(RWBL4.results$S.grazed, level = 0.85)


(RWBL.real <- as.data.frame(RWBL4.results$S.grazed$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Stage = case_when(
      grepl("20210", Group) ~ "Incubating",
      grepl("20220", Group) ~ "Incubating",
      grepl("20230", Group) ~ "Incubating",
      grepl("20211", Group) ~ "Nestling",
      grepl("20221", Group) ~ "Nestling",
      grepl("20231", Group) ~ "Nestling")) |> 
    group_by(Stage) |> 
    summarize(estimate = mean(estimate), 
              se = mean(se), 
              lcl = mean(lcl), 
              ucl = mean(ucl)))


# Plotting beta coefficients ----------------------------------------------


RWBL.beta <- coef(RWBL4.results$S.grazed) |>
  cbind(confint(RWBL4.results$S.grazed, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

RWBL.beta$Variable <- gsub("S:", "", RWBL.beta$Variable)

str(RWBL.beta)

(RWBL.plot <- ggplot(RWBL.beta[4:5,], 
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

max(RWBL.ddl$S.Age)

filter(RWBL.ddl, 
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

plotdata <- RWBL4.results$S.grazed

grazed.values <- seq(from = min(RWBL.surv$grazed), 
                     to = max(RWBL.surv$grazed), 
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

(RWBLstage.plot <- ggplot(transform(stage.pred$estimates,
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

(RWBLgrazed.plot <- ggplot(transform(grazed.pred$estimates,
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

ggsave(RWBLgrazed.plot,
       filename = "outputs/figs/grazedRWBL.png",
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
