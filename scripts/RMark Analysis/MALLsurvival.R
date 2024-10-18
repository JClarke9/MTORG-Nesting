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


MALL.surv <- filter(nest, Spec=="MALL" & Stage != "Laying")                                         # select out only MALL nest

test <- filter(MALL.surv,
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

MISSING <- is.na(MALL.surv$AgeFound)

sum(MISSING)

MALL.surv <- subset(MALL.surv, 
                    subset = !MISSING)

MALL.surv$Year <- factor(MALL.surv$Year,
                         levels = c("2021", "2022", "2023", "2024"))

str(MALL.surv)

MALL.surv$Nestling <- factor(MALL.surv$Nestling,
                             levels = c("0", "1"))


# Creating stage variable -------------------------------------------------


x <- create.stage.var(MALL.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(MALL.surv$LastChecked)), 
                      12)

MALL.surv <- bind_cols(MALL.surv, x)

rm(list = ls()[!ls() %in% c("MALL.surv")])


# Daily survival rate models ----------------------------------------------


MALL.pr <- process.data(MALL.surv,
                        nocc = max(MALL.surv$LastChecked),
                        groups = c("Year"),
                        model = "Nest")



# Temporal candidate model set
MALL1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  MALL.model.list = create.model.list("Nest")
  MALL1.results = mark.wrapper(MALL.model.list,
                               data = MALL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MALL1.results <- MALL1.run()
MALL1.results

coef(MALL1.results$S.null)



# Biological candidate model set
MALL2.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + NestAge)
  
  MALL.model.list = create.model.list("Nest")
  MALL2.results = mark.wrapper(MALL.model.list,
                               data = MALL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MALL2.results <- MALL2.run()
MALL2.results

coef(MALL2.results$S.null)


# Grazing candidate model set
MALL3.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + grazed)
  
  # 4. DSR varies with the previous Year + NestAges grazing intensity
  S.pDoD = list(formula = ~1 + pDoD)
  
  MALL.model.list = create.model.list("Nest")
  MALL3.results = mark.wrapper(MALL.model.list,
                               data = MALL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MALL3.results <- MALL3.run()
MALL3.results

coef(MALL3.results$S.null)



# Vegetation candidate model set
MALL4.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Litter)
  
  # 5. DSR varies with Bare - most of the bare was 0
  # S.bare = list(formula =  ~1 + Bare)
  
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
  
  MALL.model.list = create.model.list("Nest")
  MALL4.results = mark.wrapper(MALL.model.list,
                               data = MALL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
MALL4.results <- MALL4.run()
MALL4.results

coef(MALL4.results$S.lit)
confint(MALL4.results$S.lit, level = 0.85)

(MALL.dsr <- as.data.frame(MALL4.results$S.lit$results$real))


# Plotting beta coefficients ----------------------------------------------


MALL.beta <- coef(MALL4.results$S.lit) |>
  cbind(confint(MALL4.results$S.lit, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

MALL.beta$Variable <- gsub("S:", "", MALL.beta$Variable)

str(MALL.beta)

(MALL.plot <- ggplot(MALL.beta[2,], 
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
    labs(title = "Mallard",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


MALL.ddl <- make.design.data(MALL.pr) |> 
  as.data.frame()

Litter.values <- seq(from = min(MALL.surv$Litter),
                 to = max(MALL.surv$Litter),
                 length = 100)


lit.pred <- covariate.predictions(MALL4.results$S.lit,
                                  data = data.frame(Litter = Litter.values),
                                  indices = 1)

(MALLlit.plot <- ggplot(lit.pred$estimates, 
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
    labs(title = "Nest Survival",
         color = "Year",
         x = "Litter (Percent Cover)",
         y = "Daily Survival Rate"))


ggsave(MALL.plot,
       filename = "outputs/figs/betaMALL.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(MALLlit.plot,
       filename = "outputs/figs/litMALL.png",
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
