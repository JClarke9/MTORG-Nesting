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


CCSP.surv <- filter(nest, 
                    Spec == "CCSP" & Stage != "Laying")

test <- filter(CCSP.surv,
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

MISSING <- is.na(CCSP.surv$AgeFound)

sum(MISSING)

CCSP.surv <- subset(CCSP.surv, 
                    subset = !MISSING)

CCSP.surv$Year <- factor(CCSP.surv$Year,
                         levels = c("2021", "2022", "2023"))

CCSP.surv$Nestling <- factor(CCSP.surv$Nestling,
                             level = c("0", "1"))

str(CCSP.surv)


# Creating stage variable -------------------------------------------------


x <- create.stage.var(CCSP.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(CCSP.surv$LastChecked)), 
                      12)

CCSP.surv <- bind_cols(CCSP.surv, x)

rm(list = ls()[!ls() %in% c("CCSP.surv")])


# Daily survival rate models ----------------------------------------------


CCSP.pr <- process.data(CCSP.surv,
                        nocc = max(CCSP.surv$LastChecked),
                        groups = c("Year",
                                   "Nestling"),
                        model = "Nest")

# Temporal candidate model set
CCSP1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP1.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP1.results <- CCSP1.run()
CCSP1.results

coef(CCSP1.results$S.year)
confint(CCSP1.results$S.year, level = 0.85)


# Biological candidate model set
CCSP2.run <- function()
{
  # 1. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 2. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Year + BHCONum)
  
  # 2. DSR varies with BHCO number
  S.bhcop = list(formula = ~1 + Year + BHCOPres)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Year + NestAge)
  
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Year + Nestling)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP2.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP2.results <- CCSP2.run()
CCSP2.results

coef(CCSP2.results$S.stage)
confint(CCSP2.results$S.stage, level = 0.85)


# Grazing candidate model set
CCSP3.run <- function()
{
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Year + Nestling)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + Nestling + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Year + Nestling + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pDoD = list(formula = ~1 + Year + Nestling + pDoD)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP3.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP3.results <- CCSP3.run()
CCSP3.results

coef(CCSP3.results$S.stage)
confint(CCSP3.results$S.stage, level = 0.85)


# Vegetation candidate model set
CCSP4.run <- function()
{
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Year + Nestling)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Year + Nestling + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Year + Nestling + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Year + Nestling + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Year + Nestling + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Year + Nestling + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Year + Nestling + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Year + Nestling + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + Nestling + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + Nestling + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + Nestling + VOR)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP4.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP4.results <- CCSP4.run()
CCSP4.results

coef(CCSP4.results$S.stage)
confint(CCSP4.results$S.stage, level = 0.85)


(CCSP.real <- as.data.frame(CCSP4.results$S.stage$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023"),
      Stage = case_when(
        grepl("20210", Group) ~ "Incubating",
        grepl("20220", Group) ~ "Incubating",
        grepl("20230", Group) ~ "Incubating",
        grepl("20211", Group) ~ "Nestling",
        grepl("20221", Group) ~ "Nestling",
        grepl("20231", Group) ~ "Nestling")) |> 
    select(Year, Stage, estimate, se, lcl, ucl))


# Plotting beta coefficients ----------------------------------------------


CCSP.beta <- coef(CCSP4.results$S.stage) |>
  cbind(confint(CCSP4.results$S.stage, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

CCSP.beta$Variable <- gsub("S:", "", CCSP.beta$Variable)

str(CCSP.beta)

(CCSP.plot <- ggplot(CCSP.beta[4,], 
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
    labs(title = "Clay-colored Sparrow",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


CCSP.ddl <- make.design.data(CCSP.pr) |> 
  as.data.frame()

filter(CCSP.ddl, S.age == 0)

plotdata <- CCSP4.results$S.stage

stage.pred <- covariate.predictions(plotdata,
                                    indices = c(1, 79, 157,
                                                235, 313, 391))

inc2021 <- which(stage.pred$estimates$par.index == 1)
inc2022 <- which(stage.pred$estimates$par.index == 79)
inc2023 <- which(stage.pred$estimates$par.index == 157)
nst2021 <- which(stage.pred$estimates$par.index == 235)
nst2022 <- which(stage.pred$estimates$par.index == 313)
nst2023 <- which(stage.pred$estimates$par.index == 391)

stage.pred$estimates$Year <- NA
stage.pred$estimates$Year[inc2021] <- "2021"
stage.pred$estimates$Year[inc2022] <- "2022"
stage.pred$estimates$Year[inc2023] <- "2023"
stage.pred$estimates$Year[nst2021] <- "2021"
stage.pred$estimates$Year[nst2022] <- "2022"
stage.pred$estimates$Year[nst2023] <- "2023"

stage.pred$estimates$Stage <- NA
stage.pred$estimates$Stage[inc2021] <- "Incubating"
stage.pred$estimates$Stage[inc2022] <- "Incubating"
stage.pred$estimates$Stage[inc2023] <- "Incubating"
stage.pred$estimates$Stage[nst2021] <- "Nestling"
stage.pred$estimates$Stage[nst2022] <- "Nestling"
stage.pred$estimates$Stage[nst2023] <- "Nestling"

(CCSPstage.plot <- ggplot(transform(stage.pred$estimates,
                                    Year = factor(Year, levels = c("2021", "2022", "2023"))), 
                          aes(x = Stage, 
                              y = estimate,
                              groups = Year,
                              fill = Year)) +
    geom_point(size = 4,
               aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#A2A4A2',
                                   '#717F5B',
                                   '#D4A634')) +
    scale_fill_manual(values = c('#A2A4A2',
                                 '#717F5B',
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
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = c(.85, .1),
          legend.box = "horizontal") +
    labs(title = "Clay-colored Sparrow",
         color = "Year",
         x = "Stage",
         y = "Daily Survival Rate"))


ggsave(CCSP.plot,
       filename = "outputs/figs/betaCCSP.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(CCSPstage.plot,
       filename = "outputs/figs/stageCCSP.png",
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
