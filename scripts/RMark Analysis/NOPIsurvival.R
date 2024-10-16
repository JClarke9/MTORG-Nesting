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


NOPI.surv <- filter(nest, Spec=="NOPI" & Stage != "Laying")                                         # select out only NOPI nest

test <- filter(NOPI.surv,
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

MISSING <- is.na(NOPI.surv$AgeFound)

sum(MISSING)

NOPI.surv <- subset(NOPI.surv, 
                    subset = !MISSING)

NOPI.surv$Year <- factor(NOPI.surv$Year,
                         levels = c("2021", "2022", "2023", "2024"))

str(NOPI.surv)

NOPI.surv$Nestling <- factor(NOPI.surv$Nestling,
                             levels = c("0", "1"))


# Creating stage variable -------------------------------------------------


x <- create.stage.var(NOPI.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(NOPI.surv$LastChecked)), 
                      12)

NOPI.surv <- bind_cols(NOPI.surv, x)

rm(list = ls()[!ls() %in% c("NOPI.surv")])


# Daily survival rate models ----------------------------------------------


NOPI.pr <- process.data(NOPI.surv,
                        nocc = max(NOPI.surv$LastChecked),
                        groups = c("Year"),
                        model = "Nest")



# Temporal candidate model set
NOPI1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  NOPI.model.list = create.model.list("Nest")
  NOPI1.results = mark.wrapper(NOPI.model.list,
                               data = NOPI.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
NOPI1.results <- NOPI1.run()
NOPI1.results

coef(NOPI1.results$S.null)



# Biological candidate model set
NOPI2.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + NestAge)
  
  NOPI.model.list = create.model.list("Nest")
  NOPI2.results = mark.wrapper(NOPI.model.list,
                               data = NOPI.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
NOPI2.results <- NOPI2.run()
NOPI2.results

coef(NOPI2.results$S.null)


# Grazing candidate model set
NOPI3.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + grazed)
  
  # 4. DSR varies with the previous Year + NestAges grazing intensity
  S.pDoD = list(formula = ~1 + pDoD)
  
  NOPI.model.list = create.model.list("Nest")
  NOPI3.results = mark.wrapper(NOPI.model.list,
                               data = NOPI.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
NOPI3.results <- NOPI3.run()
NOPI3.results

coef(NOPI3.results$S.null)



# Vegetation candidate model set
NOPI4.run <- function()
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
  
  NOPI.model.list = create.model.list("Nest")
  NOPI4.results = mark.wrapper(NOPI.model.list,
                               data = NOPI.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
NOPI4.results <- NOPI4.run()
NOPI4.results

coef(NOPI4.results$S.lit)
confint(NOPI4.results$S.lit, level = 0.85)

coef(NOPI4.results$S.vor)
confint(NOPI4.results$S.vor, level = 0.85)

(NOPI.dsrLIT <- as.data.frame(NOPI4.results$S.lit$results$real))
(NOPI.dsrVOR <- as.data.frame(NOPI4.results$S.vor$results$real))

NOPI5.results <- mark(NOPI.surv,
                      nocc = max(NOPI.surv$LastChecked),
                      model = "Nest",
                      groups = "Year",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula = ~1 + Litter + VOR)))

coef(NOPI.final)
confint(NOPI.final, level = 0.85)



# Plotting beta coefficients ----------------------------------------------


NOPI.beta <- coef(NOPI5.results) |>
  cbind(confint(NOPI5.results, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

NOPI.beta$Variable <- gsub("S:", "", NOPI.beta$Variable)

str(NOPI.beta)

(NOPI.plot <- ggplot(NOPI.beta[2:3,], 
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
    labs(title = "Northern Pintail",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


NOPI.ddl <- make.design.data(NOPI.pr) |> 
  as.data.frame()

LITvalues <- seq(from = min(NOPI.surv$Litter),
                 to = max(NOPI.surv$Litter),
                 length = 100)

VORvalues <- seq(from = min(NOPI.surv$VOR),
                 to = max(NOPI.surv$VOR),
                 length = 100)


VOR.pred <- covariate.predictions(NOPI5.results,
                                  data = data.frame(VOR = VORvalues,
                                                    Litter = mean(LITvalues)),
                                  indices = 1)

(NOPIvor.plot <- ggplot(VOR.pred$estimates, 
                        aes(x = VOR, 
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
    labs(title = "Northern Pintail",
         color = "Year",
         x = "Visual Obstruction (dm)",
         y = "Daily Survival Rate"))


LIT.pred <- covariate.predictions(NOPI5.results,
                                  data = data.frame(VOR = mean(VORvalues),
                                                    Litter = LITvalues),
                                  indices = 1)

(NOPIlit.plot <- ggplot(LIT.pred$estimates, 
                        aes(x = Litter, 
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
    labs(title = "Northern Pintail",
         color = "Year",
         x = "Litter (Percent Cover)",
         y = "Daily Survival Rate"))

ggsave(NOPI.plot,
       filename = "outputs/figs/betaNOPI.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(NOPIvor.plot,
       filename = "outputs/figs/vorNOPI.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(NOPIlit.plot,
       filename = "outputs/figs/litNOPI.png",
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
