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


RWBL.surv <- filter(nest, 
                    Spec=="RWBL")                                         # select out only RWBL nest

# I didn't filter for bare or Litter Depth because these were often
# over water so we didn't measure that
test <- filter(RWBL.surv,
               is.na(KBG) |
                 is.na(SmoothB) |
                 is.na(Litter) |
                 is.na(Forb) |
                 is.na(Grasslike) |
                 is.na(Woody) |
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
                        nocc=max(RWBL.surv$LastChecked),
                        groups = c("Year"),
                        model="Nest")

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

coef(RWBL1.results$S.time)
coef(RWBL1.results$S.quad)

confint(RWBL1.results$S.time, level = 0.85)
confint(RWBL1.results$S.quad, level = 0.85)



# Biological candidate model set
RWBL2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 2. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Time + BHCONum)
  
  # 2. DSR varies with BHCO number
  S.bhcop = list(formula = ~1 + Time + BHCOPres)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + Time + NestAge)
  
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Time + Incub)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL2.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL2.results <- RWBL2.run()
RWBL2.results

coef(RWBL2.results$S.time)
confint(RWBL2.results$S.time, level = 0.85)


# Grazing candidate model set
RWBL3.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Time + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Time + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pTreat = list(formula = ~1 + Time + pTreat)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL3.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL3.results <- RWBL3.run()
RWBL3.results

coef(RWBL3.results$S.time)
confint(RWBL3.results$S.time, level = 0.85)


# Vegetation candidate model set
RWBL4.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Time + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Time + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Time + Litter)
  
  # 5. DSR varies with Bare
  #S.bare = list(formula =  ~1 + Time + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Time + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Time + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Time + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  #S.litdep = list(formula =  ~1 + Time + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Time + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Time + VOR)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL4.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
RWBL4.results <- RWBL4.run()
RWBL4.results

coef(RWBL4.results$S.time)
confint(RWBL4.results$S.time, level = 0.85)

RWBL.real <- as.data.frame(RWBL4.results$S.kbg$results$real) |> 
  summarise(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))


# Plotting beta coefficients ----------------------------------------------


RWBL.beta <- coef(RWBL4.results$S.time) |>
  cbind(confint(RWBL4.results$S.time, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

RWBL.beta$Variable <- gsub("S:", "", RWBL.beta$Variable)

str(RWBL.beta)

(RWBL.plot <- ggplot(RWBL.beta[2,], 
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
    labs(title = "Red-winged Blackbird",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


RWBL.ddl <- make.design.data(RWBL.pr) |> 
  as.data.frame()

plotdata <- RWBL4.results$S.time

Time.pred <- covariate.predictions(plotdata,
                                   indices = c(1:73))

Time.pred$estimates$group <- 1

(RWBLtime.plot <- ggplot(Time.pred$estimates, 
                         aes(x = index, 
                             y = estimate)) +
    geom_line(linewidth = 1.5,
              aes(color = "group")) +
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
    labs(title = "Red-winged Blackbird",
         color = "Year",
         x = "Julian Day",
         y = "Daily Survival Rate"))


ggsave(RWBL.plot,
       filename = "outputs/figs/RWBLbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(RWBLtime.plot,
       filename = "outputs/figs/RWBLtime.png",
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
