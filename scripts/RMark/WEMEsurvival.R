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

WEME.surv <- filter(nest, 
                    Spec=="WEME")                                         # select out only WEME nest

test <- filter(WEME.surv,
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

MISSING <- is.na(WEME.surv$AgeFound)

sum(MISSING)

WEME.surv <- subset(WEME.surv, 
                    subset = !MISSING)

WEME.surv$Year <- factor(WEME.surv$Year,
                         levels = c("2021", "2022", "2023"))

str(WEME.surv)

# Creating stage variable -------------------------------------------------

x <- create.stage.var(WEME.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(WEME.surv$LastChecked)), 
                      12)

WEME.surv <- bind_cols(WEME.surv, x)

rm(list = ls()[!ls() %in% c("WEME.surv")])

# Daily survival rate models ----------------------------------------------

WEME.pr <- process.data(WEME.surv,
                        nocc=max(WEME.surv$LastChecked),
                        groups = c("Year"),
                        model="Nest")

# Temporal candidate model set
WEME1.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  WEME.model.list = create.model.list("Nest")
  WEME1.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME1.results <- WEME1.run()
WEME1.results

coef(WEME1.results$S.null)

# Biological candidate model set
WEME2.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + BHCONum)
  
  # 2. DSR varies with BHCO number
  S.bhcop = list(formula = ~1 + BHCOPres)
  
  # 4. DSR varies with nest age
  S.age = list(formula = ~1 + NestAge)
  
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Incub)
  
  WEME.model.list = create.model.list("Nest")
  WEME2.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME2.results <- WEME2.run()
WEME2.results

coef(WEME2.results$S.null)


# Grazing candidate model set
WEME3.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pTreat = list(formula = ~1 + pTreat)
  
  WEME.model.list = create.model.list("Nest")
  WEME3.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME3.results <- WEME3.run()
WEME3.results

coef(WEME3.results$S.grazed)
confint(WEME3.results$S.grazed)

# Vegetation candidate model set
WEME4.run <- function()
{
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + grazed)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + grazed + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + grazed + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + grazed + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + grazed + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + grazed + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + grazed + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + grazed + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + grazed + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + grazed + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + grazed + VOR)
  
  WEME.model.list = create.model.list("Nest")
  WEME4.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME4.results <- WEME4.run()
WEME4.results

coef(WEME4.results$S.litdep)
confint(WEME4.results$S.litdep)


# Plotting beta coefficients ----------------------------------------------


WEME4.beta <- WEME4.results$S.bare$results$beta
WEME4.beta

WEME.top <- rownames_to_column(WEME4.beta, 
                               "Variable")
WEME.top <- rename(WEME.top, 
                   "Coefficient"="estimate")

WEME.top[,1] <- gsub("S:", "", WEME.top[,1])
str(WEME.top)

WEME.top$SEup <- WEME.top$Coefficient + WEME.top$se
WEME.top$SElow <- WEME.top$Coefficient - WEME.top$se

WEME.plot <- ggplot(WEME.top, aes(x = Variable,
                                  y = Coefficient)) +
  geom_hline(yintercept = 0,
             colour = gray(1/2), 
             lty = 2) +
  geom_point(aes(x = Variable,
                 y = Coefficient),
             size = 5) +
  geom_errorbar(aes(x = Variable,
                    ymin = lcl,
                    ymax = ucl),
                width = .5,
                size = 1) +
  ggtitle("WEME Beta Coefficient") +
  theme(panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="white",                           # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="white",                            # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=12, 
                                 colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=12,                                              # change the size of the axis titles
                          colour = "black")) +                                  # change the color of the axis titles
  coord_flip()

WEME.plot


# Creating predictive plots -----------------------------------------------


WEME.pr <- process.data(WEME.surv,
                        nocc=max(WEME.surv$LastChecked),
                        model="Nest")

WEME.ddl <- make.design.data(WEME.pr)
WEME.ddl <- as.data.frame(WEME.ddl)

WEME.height = seq(from = min(WEME.surv$Veg.Height), 
                  to = max(WEME.surv$Veg.Height), 
                  length = 100)

WEME.predH <- covariate.predictions(WEME.mod,
                                    data=data.frame(Veg.Height = WEME.height),
                                    indices = 1)

WEME.predH$estimates$Group <- "WEME"

WEMEheight.plot <- ggplot(WEME.predH$estimates,
                          aes(x = covdata,
                              y = estimate,
                              fill = Group)) +
  geom_line(size = 1.5,
            aes(colour = Group)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  scale_colour_manual(values = c('black')) +
  scale_fill_manual(values = c('black')) +
  ylim(.65, 1) + 
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=14,
                                  hjust=.5),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=12, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=12,                                              # change the size of the axis titles
                          colour = "black",
                          family = "my_font"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  labs(title = "Vegetation Height",
       x = "Centimeters",
       y = "Daily Survival Rate")


WEME.KBG = seq(from = min(WEME.surv$KBG),
               to = max(WEME.surv$KBG),
               length = 100)

WEME.predKBG <- covariate.predictions(WEME.mod,
                                      data=data.frame(KBG = WEME.KBG),
                                      indices = 1)

WEME.predKBG$estimates$Group <- "WEME"

WEMEkbg.plot <- ggplot(WEME.predKBG$estimates,
                       aes(x = covdata,
                           y = estimate,
                           fill = Group)) +
  geom_line(size = 1.5,
            aes(colour = Group)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  scale_colour_manual(values = c('black')) +
  scale_fill_manual(values = c('black')) +
  ylim(.65, 1) + 
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=14,
                                  hjust=.5),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=12, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=12,                                              # change the size of the axis titles
                          colour = "black",
                          family = "my_font"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  labs(title = "Kentucky Bluegrass",
       x = "Percent Cover",
       y = "Daily Survival Rate")

max(WEME.predH$estimates$estimate)
min(WEME.predH$estimates$estimate)

max(WEME.predKBG$estimates$estimate)
min(WEME.predKBG$estimates$estimate)

WEME.plot
WEMEheight.plot
WEMEkbg.plot

WEMEpred.plot <- plot_grid(WEMEheight.plot, 
                           WEMEkbg.plot + theme(axis.text.y = element_blank()) + labs(y = NULL))

title <- ggdraw() +
  draw_label("Western Meadowlark Daily Nest Survival",
             fontface = "bold",
             fontfamily = "my_font",
             size = 14) +
  theme(plot.margin = margin(0,0,7))

WEMEveg.plot <- plot_grid(title,
                           WEMEpred.plot,
                           ncol = 1,
                          rel_heights = c(0.1, 1))
WEMEveg.plot

ggsave(WEME.plot,
       filename = "outputs/figs/WEMEbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(WEMEveg.plot,
       filename = "outputs/figs/WEMEveg.png",
       dpi = "print",
       height = 3,
       width = 5)

# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all=TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
#  files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)
