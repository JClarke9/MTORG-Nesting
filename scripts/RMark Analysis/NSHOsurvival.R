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


NSHO.surv <- filter(nest, Spec=="NSHO" & Stage != "Laying")                                         # select out only NSHO nest

test <- filter(NSHO.surv,
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

MISSING <- is.na(NSHO.surv$AgeFound)

sum(MISSING)

NSHO.surv <- subset(NSHO.surv, 
                    subset = !MISSING)

NSHO.surv$Year <- factor(NSHO.surv$Year,
                         levels = c("2021", "2022", "2023", "2024"))

str(NSHO.surv)

NSHO.surv$Nestling <- factor(NSHO.surv$Nestling,
                             levels = c("0", "1"))


# Creating stage variable -------------------------------------------------


x <- create.stage.var(NSHO.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(NSHO.surv$LastChecked)), 
                      12)

NSHO.surv <- bind_cols(NSHO.surv, x)

rm(list = ls()[!ls() %in% c("NSHO.surv")])


# Daily survival rate models ----------------------------------------------


NSHO.pr <- process.data(NSHO.surv,
                        nocc = max(NSHO.surv$LastChecked),
                        groups = c("Year"),
                        model = "Nest")



# Grazing candidate model set
NSHO1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + grazed)
  
  # 4. DSR varies with the previous  + NestAges grazing intensity
  S.pDoD = list(formula = ~1 + pDoD)
  
  NSHO.model.list = create.model.list("Nest")
  NSHO1.results = mark.wrapper(NSHO.model.list,
                               data = NSHO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
NSHO1.results <- NSHO1.run()
NSHO1.results

coef(NSHO1.results$S.null)



# Temporal candidate model set
NSHO2.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  NSHO.model.list = create.model.list("Nest")
  NSHO2.results = mark.wrapper(NSHO.model.list,
                               data = NSHO.pr,
                               adjust = FALSE,
                               delete =TRUE)
}

# Results of candidate model set
NSHO2.results <- NSHO2.run()
NSHO2.results

coef(NSHO2.results$S.null)


# Biological candidate model set
NSHO3.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + NestAge)
  
  NSHO.model.list = create.model.list("Nest")
  NSHO3.results = mark.wrapper(NSHO.model.list,
                               data = NSHO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
NSHO3.results <- NSHO3.run()
NSHO3.results

coef(NSHO3.results$S.null)



# Vegetation candidate model set
NSHO4.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + VOR)
  
  NSHO.model.list = create.model.list("Nest")
  NSHO4.results = mark.wrapper(NSHO.model.list,
                               data = NSHO.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
NSHO4.results <- NSHO4.run()
NSHO4.results

coef(NSHO4.results$S.brome)
confint(NSHO4.results$S.brome, level = 0.85)

(NSHO.dsr <- as.data.frame(NSHO4.results$S.brome$results$real))


# Plotting beta coefficients ----------------------------------------------


NSHO.beta <- coef(NSHO4.results$S.brome) |>
  cbind(confint(NSHO4.results$S.brome, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

NSHO.beta$Variable <- gsub("S:", "", NSHO.beta$Variable)

str(NSHO.beta)

(NSHO.plot <- ggplot(NSHO.beta[2,], 
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
    labs(title = "Northern Shoveler",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


NSHO.ddl <- make.design.data(NSHO.pr) |> 
  as.data.frame()

BROMEvalues <- seq(from = min(NSHO.surv$SmoothB),
                   to = max(NSHO.surv$SmoothB),
                   length = 100)


BROME.pred <- covariate.predictions(NSHO4.results$S.brome,
                                    data = data.frame(SmoothB = BROMEvalues),
                                    indices = 1)

(NSHObrome.plot <- ggplot(BROME.pred$estimates, 
                        aes(x = covdata, 
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
          axis.text.y = element_text(size = 12, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = "none") +
    labs(title = "Northern Shoveler",
         color = "Year",
         x = "Smooth Brome",
         y = "Daily Survival Rate"))


ggsave(NSHO.plot,
       filename = "outputs/figs/betaNSHO.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(NSHObrome.plot,
       filename = "outputs/figs/bromeNSHO.png",
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
