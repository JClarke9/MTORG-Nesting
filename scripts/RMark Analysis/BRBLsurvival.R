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


BRBL.surv <- filter(nest, 
                    Spec == "BRBL" & Stage != "Laying")

test <- filter(BRBL.surv,
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

MISSING <- is.na(BRBL.surv$AgeFound)

sum(MISSING)

BRBL.surv <- subset(BRBL.surv, 
                    subset = !MISSING)

BRBL.surv$Year <- factor(BRBL.surv$Year,
                         levels = c("2021", "2022", "2023"))

BRBL.surv$Nestling <- factor(BRBL.surv$Nestling,
                             level = c("0", "1"))

str(BRBL.surv)


# Creating stage variable -------------------------------------------------


x <- create.stage.var(BRBL.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(BRBL.surv$LastChecked)), 
                      12)

BRBL.surv <- bind_cols(BRBL.surv, x)

rm(list = ls()[!ls() %in% c("BRBL.surv")])


# Daily survival rate models ----------------------------------------------


BRBL.pr <- process.data(BRBL.surv,
                        nocc = max(BRBL.surv$LastChecked),
                        groups = c("Year",
                                   "Nestling"),
                        model = "Nest")

# Temporal candidate model set
BRBL1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL1.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BRBL1.results <- BRBL1.run()
BRBL1.results

coef(BRBL1.results$S.year)
confint(BRBL1.results$S.year, level = 0.85)


# Biological candidate model set
BRBL2.run <- function()
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
  
  BRBL.model.list = create.model.list("Nest")
  BRBL2.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BRBL2.results <- BRBL2.run()
BRBL2.results

coef(BRBL2.results$S.stage)
confint(BRBL2.results$S.stage, level = 0.85)


# Grazing candidate model set
BRBL3.run <- function()
{
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Year + Nestling)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + Nestling + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Year + Nestling + grazep)
  
  # 4. DSR varies with the previous Year + Nestlings grazing intensity
  S.pDoD = list(formula = ~1 + Year + Nestling + pDoD)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL3.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BRBL3.results <- BRBL3.run()
BRBL3.results

coef(BRBL3.results$S.stage)
confint(BRBL3.results$S.stage, level = 0.85)


# Vegetation candidate model set
BRBL4.run <- function()
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
  
  BRBL.model.list = create.model.list("Nest")
  BRBL4.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BRBL4.results <- BRBL4.run()
BRBL4.results

coef(BRBL4.results$S.lit)
confint(BRBL4.results$S.lit, level = 0.85)

BRBL4.results$S.lit$results$real |> 
  summarize(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))

(BRBL.real <- as.data.frame(BRBL4.results$S.lit$results$real) |> 
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

(BRBL.year <- BRBL.real |> 
    group_by(Year) |> 
    summarize(mean = mean(estimate)))

(BRBL.stage <- BRBL.real |> 
    group_by(Stage) |> 
    summarize(mean = mean(estimate)))


# Plotting Beta Coefficients ----------------------------------------------


BRBL.beta <- coef(BRBL4.results$S.lit) |>
  cbind(confint(BRBL4.results$S.lit, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

BRBL.beta$Variable <- gsub("S:", "", BRBL.beta$Variable)

str(BRBL.beta)

(BRBL.plot <- ggplot(BRBL.beta[4:5,], 
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
    labs(title = "BRBL",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating Predictive Plots -----------------------------------------------


BRBL.ddl <- make.design.data(BRBL.pr) |> 
  as.data.frame()

filter(BRBL.ddl, S.age == 0)

plotdata <- BRBL4.results$S.lit

lit.values <- seq(from = min(BRBL.surv$Litter), 
                  to = max(BRBL.surv$Litter), 
                  length = 100)


lit.pred <- covariate.predictions(plotdata,
                                  data = data.frame(Litter = lit.values),
                                  indices = c(1, 77,
                                              153, 229, 305))

inc2021 <- which(lit.pred$estimates$par.index == 1)
inc2023 <- which(lit.pred$estimates$par.index == 77)
nst2021 <- which(lit.pred$estimates$par.index == 153)
nst2022 <- which(lit.pred$estimates$par.index == 229)
nst2023 <- which(lit.pred$estimates$par.index == 305)

lit.pred$estimates$Year <- NA
lit.pred$estimates$Year[inc2021] <- "2021"
lit.pred$estimates$Year[inc2023] <- "2023"
lit.pred$estimates$Year[nst2021] <- "2021"
lit.pred$estimates$Year[nst2022] <- "2022"
lit.pred$estimates$Year[nst2023] <- "2023"

lit.pred$estimates$Stage <- NA
lit.pred$estimates$Stage[inc2021] <- "Incubating"
lit.pred$estimates$Stage[inc2023] <- "Incubating"
lit.pred$estimates$Stage[nst2021] <- "Nestling"
lit.pred$estimates$Stage[nst2022] <- "Nestling"
lit.pred$estimates$Stage[nst2023] <- "Nestling"

(BRBLlit.plot <- ggplot(transform(lit.pred$estimates,
                                  Year = factor(Year, levels = c("2021", "2022", "2023"))), 
                        aes(x = covdata, 
                            y = estimate,
                            groups = Year,
                            fill = Year)) +
    geom_line(linewidth = 1.5,
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
    facet_grid(~Stage) + 
    labs(title = "Brewer's Blackbird",
         color = "Year",
         x = "Litter (Percent Cover)",
         y = "Daily Survival Rate"))


ggsave(BRBL.plot,
       filename = "outputs/figs/betaBRBL.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BRBLlit.plot,
       filename = "outputs/figs/litBRBL.png",
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
