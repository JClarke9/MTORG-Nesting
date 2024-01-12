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


BRBL.surv <- filter(nest, 
                    Spec=="BRBL")                                         # select out only BRBL nest

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
                        nocc=max(BRBL.surv$LastChecked),
                        groups = c("Year"),
                        model="Nest")

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
  S.stage = list(formula = ~1 + Year + Incub)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL2.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BRBL2.results <- BRBL2.run()
BRBL2.results

coef(BRBL2.results$S.year)
confint(BRBL2.results$S.year, level = 0.85)


# Grazing candidate model set
BRBL3.run <- function()
{
  # 1. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Year + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pTreat = list(formula = ~1 + Year + pTreat)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL3.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BRBL3.results <- BRBL3.run()
BRBL3.results

coef(BRBL3.results$S.year)
confint(BRBL3.results$S.year, level = 0.85)


# Vegetation candidate model set
BRBL4.run <- function()
{
  # 1. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Year + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Year + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Year + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Year + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Year + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Year + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Year + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + VOR)
  
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


BRBL.real <- as.data.frame(BRBL4.results$S.lit$results$real) |>
  rownames_to_column(var = "Group") |> 
  mutate(Year = case_when(
    grepl("2021", Group) ~ "2021",
    grepl("2022", Group) ~ "2022",
    grepl("2023", Group) ~ "2023")) |> 
  select(Year, estimate, se, lcl, ucl)


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

(BRBL.plot <- ggplot(BRBL.beta[4,], 
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
    labs(title = "BRBL",
         x = "Percent Cover",
         y = expression("Beta " (beta))))


# Creating Predictive Plots -----------------------------------------------


plotdata <- BRBL4.results$S.lit

BRBL.ddl <- make.design.data(BRBL.pr)
BRBL.ddl <- as.data.frame(BRBL.ddl)

lit.values <- seq(from = min(BRBL.surv$Litter), 
                 to = max(BRBL.surv$Litter), 
                 length = 100)


lit.pred <- covariate.predictions(plotdata,
                                  data = data.frame(Litter = lit.values),
                                  indices = c(1, 88, 175))

`2021` <- which(lit.pred$estimates$par.index == 1)
`2022` <- which(lit.pred$estimates$par.index == 88)
`2023` <- which(lit.pred$estimates$par.index == 175)

lit.pred$estimates$Year <- NA
lit.pred$estimates$Year[`2021`] <- "2021"
lit.pred$estimates$Year[`2022`] <- "2022"
lit.pred$estimates$Year[`2023`] <- "2023"

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
          legend.box = "horizontal") +
    labs(title = "Brewer's Blackbird",
         color = "Year",
         x = "Litter (Percent Cover)",
         y = "Daily Survival Rate"))


ggsave(BRBL.plot,
       filename = "outputs/figs/BRBLbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BRBLlit.plot,
       filename = "outputs/figs/BRBLlit.png",
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
