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

CCSP.surv <- filter(nest, 
                    Spec=="CCSP")                                         # select out only CCSP nest

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
                        nocc=max(CCSP.surv$LastChecked),
                        groups = c("Year"),
                        model="Nest")

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
  S.stage = list(formula = ~1 + Year + Incub)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP2.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP2.results <- CCSP2.run()
CCSP2.results

coef(CCSP2.results$S.bhcon)
confint(CCSP2.results$S.bhcon, level = 0.85)


# Grazing candidate model set
CCSP3.run <- function()
{
  # 1. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Year + BHCONum)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + BHCONum + grazed)
  
  # 3. DSR varies with the number of days a nest experienced grazing
  S.grazep = list(formula = ~1 + Year + BHCONum + grazep)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pTreat = list(formula = ~1 + Year + BHCONum + pTreat)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP3.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP3.results <- CCSP3.run()
CCSP3.results

coef(CCSP2.results$S.bhcon)
confint(CCSP2.results$S.bhcon, level = 0.85)

# Vegetation candidate model set
CCSP4.run <- function()
{
  # 1. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Year + BHCONum)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Year + BHCONum + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Year + BHCONum + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Year + BHCONum + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Year + BHCONum + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Year + BHCONum + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Year + BHCONum + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Year + BHCONum + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + BHCONum + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + BHCONum + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + BHCONum + VOR)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP4.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP4.results <- CCSP4.run()
CCSP4.results

summary(CCSP.avg)
confint(CCSP.avg, level = 0.85)


CCSP4.results$S.bare$results$beta
CCSP4.results$S.brome$results$beta
CCSP4.results$S.vor$results$beta
CCSP4.results$S.kbg$results$beta

CCSP4.results$S.height$results$real

CCSP.avg <- model.avg(CCSP4.results,
                      rank = AIC)

(CCSP4.avg <- model.average(CCSP4.results,
                            parameter = "S",
                            vcv = TRUE))

CCSP.dsr <- rename(CCSP4.avg,
                   "Year" = "par.index")

CCSP.dsr$Year <- case_match(CCSP.dsr$Year,
                            1 ~ "2021",
                            349 ~ "2022",
                            697 ~ "2023")

CCSP.dsr


# Plotting beta coefficients ----------------------------------------------


CCSP4.beta <- CCSP4.results$S.vor$results$beta
CCSP4.beta

CCSP.top <- rownames_to_column(CCSP4.beta, 
                               "Variable")

CCSP.top <- rename(CCSP.top, 
                   "Coefficient"="estimate")

CCSP.top[,1] <- gsub("S:", "", CCSP.top[,1])
str(CCSP.top)

CCSP.top$SEup <- CCSP.top$Coefficient + CCSP.top$se
CCSP.top$SElow <- CCSP.top$Coefficient - CCSP.top$se

(CCSP.plot <- ggplot(CCSP.top, 
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
    ggtitle("CCSP Beta Coefficient") +
    theme(panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill="white",                     # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill="white",                      # make the outer background transparent
                                         colour=NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text = element_text(size=12, 
                                   colour = "black"),                    # color the axis text
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size=12,                                              # change the size of the axis titles
                            colour = "black")) +                                    # change the color of the axis titles
    coord_flip())

# Creating predictive plots -----------------------------------------------

CCSP.ddl <- make.design.data(CCSP.pr)
CCSP.ddl <- as.data.frame(CCSP.ddl)

plotdata <- CCSP4.results$S.vor

VORvalues = seq(from = quantile(CCSP.surv$VOR, 0.05),
                to = quantile(CCSP.surv$VOR, 0.95),
                length = 100)

BHCONum.values = seq(from = quantile(CCSP.surv$BHCONum, 0.05),
                     to = max(CCSP.surv$BHCONum, 0.95),
                     length = 100)

vor.pred <- covariate.predictions(plotdata,
                                  data = data.frame(VOR = VORvalues,
                                                    BHCONum = mean(BHCONum.values)),
                                  indices = c(1, 88, 175, 262,
                                              349, 436, 523, 610,
                                              697, 784, 871, 958))

Rest.21 <- which(vor.pred$estimates$par.index == 1)
Moderate.21 <- which(vor.pred$estimates$par.index == 88)
Full.21 <- which(vor.pred$estimates$par.index == 175)
Heavy.21 <- which(vor.pred$estimates$par.index == 262)
Rest.22 <- which(vor.pred$estimates$par.index == 349)
Moderate.22 <- which(vor.pred$estimates$par.index == 436)
Full.22 <- which(vor.pred$estimates$par.index == 523)
Heavy.22 <- which(vor.pred$estimates$par.index == 610)
Rest.23 <- which(vor.pred$estimates$par.index == 697)
Moderate.23 <- which(vor.pred$estimates$par.index == 784)
Full.23 <- which(vor.pred$estimates$par.index == 871)
Heavy.23 <- which(vor.pred$estimates$par.index == 958)

vor.pred$estimates$Group <- NA
vor.pred$estimates$Group[Rest.21] <- "Rest 2021"
vor.pred$estimates$Group[Moderate.21] <- "Moderate 2021"
vor.pred$estimates$Group[Full.21] <- "Full 2021"
vor.pred$estimates$Group[Heavy.21] <- "Heavy 2022"
vor.pred$estimates$Group[Rest.22] <- "Rest 2022"
vor.pred$estimates$Group[Moderate.22] <- "Moderate 2022"
vor.pred$estimates$Group[Full.22] <- "Full 2022"
vor.pred$estimates$Group[Heavy.22] <- "Heavy 2022"
vor.pred$estimates$Group[Rest.23] <- "Rest 2023"
vor.pred$estimates$Group[Moderate.23] <- "Moderate 2023"
vor.pred$estimates$Group[Full.23] <- "Full 2023"
vor.pred$estimates$Group[Heavy.23] <- "Heavy 2023"

vor.pred$estimates$Treat <- NA
vor.pred$estimates$Treat[Rest.21] <- "Rest"
vor.pred$estimates$Treat[Moderate.21] <- "Moderate"
vor.pred$estimates$Treat[Full.21] <- "Full"
vor.pred$estimates$Treat[Heavy.21] <- "Heavy"
vor.pred$estimates$Treat[Rest.22] <- "Rest"
vor.pred$estimates$Treat[Moderate.22] <- "Moderate"
vor.pred$estimates$Treat[Full.22] <- "Full"
vor.pred$estimates$Treat[Heavy.22] <- "Heavy"
vor.pred$estimates$Treat[Rest.23] <- "Rest"
vor.pred$estimates$Treat[Moderate.23] <- "Moderate"
vor.pred$estimates$Treat[Full.23] <- "Full"
vor.pred$estimates$Treat[Heavy.23] <- "Heavy"

vor.pred$estimates$Year <- NA
vor.pred$estimates$Year[Rest.21] <- "2021"
vor.pred$estimates$Year[Moderate.21] <- "2021"
vor.pred$estimates$Year[Full.21] <- "2021"
vor.pred$estimates$Year[Heavy.21] <- "2021"
vor.pred$estimates$Year[Rest.22] <- "2022"
vor.pred$estimates$Year[Moderate.22] <- "2022"
vor.pred$estimates$Year[Full.22] <- "2022"
vor.pred$estimates$Year[Heavy.22] <- "2022"
vor.pred$estimates$Year[Rest.23] <- "2023"
vor.pred$estimates$Year[Moderate.23] <- "2023"
vor.pred$estimates$Year[Full.23] <- "2023"
vor.pred$estimates$Year[Heavy.23] <- "2023"
head(vor.pred$estimates)

(CCSPvor.plot <- ggplot(transform(vor.pred$estimates,
                                  Year = factor(Year, levels=c("2021", "2022", "2023"))), 
                        aes(x = VOR, 
                            y = estimate,
                            groups = Year,
                            fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#A2A4A2',
                                   '#717F5B',
                                   '#D4A634')) +
    scale_fill_manual(values = c('black',
                                 '#717F5B',
                                 '#D4A634')) +
    xlab("VOR (dm)") +
    ylab("Daily Survival Rate") +
    ylim(0.7, 1) + 
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
    labs(title = "Clay-colored Sparrow",
         color = "Year"))

bhco.pred <- covariate.predictions(plotdata,
                                   data=data.frame(VOR = mean(VORvalues),
                                                   BHCONum = BHCONum.values),
                                   indices = c(1, 88, 175, 262,
                                               349, 436, 523, 610,
                                               697, 784, 871, 958))

Rest.21 <- which(bhco.pred$estimates$par.index == 1)
Moderate.21 <- which(bhco.pred$estimates$par.index == 88)
Full.21 <- which(bhco.pred$estimates$par.index == 175)
Heavy.21 <- which(bhco.pred$estimates$par.index == 262)
Rest.22 <- which(bhco.pred$estimates$par.index == 349)
Moderate.22 <- which(bhco.pred$estimates$par.index == 436)
Full.22 <- which(bhco.pred$estimates$par.index == 523)
Heavy.22 <- which(bhco.pred$estimates$par.index == 610)
Rest.23 <- which(bhco.pred$estimates$par.index == 697)
Moderate.23 <- which(bhco.pred$estimates$par.index == 784)
Full.23 <- which(bhco.pred$estimates$par.index == 871)
Heavy.23 <- which(bhco.pred$estimates$par.index == 958)

bhco.pred$estimates$Group <- NA
bhco.pred$estimates$Group[Rest.21] <- "Rest 2021"
bhco.pred$estimates$Group[Moderate.21] <- "Moderate 2021"
bhco.pred$estimates$Group[Full.21] <- "Full 2021"
bhco.pred$estimates$Group[Heavy.21] <- "Heavy 2022"
bhco.pred$estimates$Group[Rest.22] <- "Rest 2022"
bhco.pred$estimates$Group[Moderate.22] <- "Moderate 2022"
bhco.pred$estimates$Group[Full.22] <- "Full 2022"
bhco.pred$estimates$Group[Heavy.22] <- "Heavy 2022"
bhco.pred$estimates$Group[Rest.23] <- "Rest 2023"
bhco.pred$estimates$Group[Moderate.23] <- "Moderate 2023"
bhco.pred$estimates$Group[Full.23] <- "Full 2023"
bhco.pred$estimates$Group[Heavy.23] <- "Heavy 2023"

bhco.pred$estimates$Treat <- NA
bhco.pred$estimates$Treat[Rest.21] <- "Rest"
bhco.pred$estimates$Treat[Moderate.21] <- "Moderate"
bhco.pred$estimates$Treat[Full.21] <- "Full"
bhco.pred$estimates$Treat[Heavy.21] <- "Heavy"
bhco.pred$estimates$Treat[Rest.22] <- "Rest"
bhco.pred$estimates$Treat[Moderate.22] <- "Moderate"
bhco.pred$estimates$Treat[Full.22] <- "Full"
bhco.pred$estimates$Treat[Heavy.22] <- "Heavy"
bhco.pred$estimates$Treat[Rest.23] <- "Rest"
bhco.pred$estimates$Treat[Moderate.23] <- "Moderate"
bhco.pred$estimates$Treat[Full.23] <- "Full"
bhco.pred$estimates$Treat[Heavy.23] <- "Heavy"

bhco.pred$estimates$Year <- NA
bhco.pred$estimates$Year[Rest.21] <- "2021"
bhco.pred$estimates$Year[Moderate.21] <- "2021"
bhco.pred$estimates$Year[Full.21] <- "2021"
bhco.pred$estimates$Year[Heavy.21] <- "2021"
bhco.pred$estimates$Year[Rest.22] <- "2022"
bhco.pred$estimates$Year[Moderate.22] <- "2022"
bhco.pred$estimates$Year[Full.22] <- "2022"
bhco.pred$estimates$Year[Heavy.22] <- "2022"
bhco.pred$estimates$Year[Rest.23] <- "2023"
bhco.pred$estimates$Year[Moderate.23] <- "2023"
bhco.pred$estimates$Year[Full.23] <- "2023"
bhco.pred$estimates$Year[Heavy.23] <- "2023"
head(bhco.pred$estimates)

(CCSPbhco.plot <- ggplot(transform(bhco.pred$estimates,
                                   Year = factor(Year, levels=c("2021", "2022", "2023"))), 
                         aes(x = BHCONum, 
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
    xlab("BHCO Eggs") +
    ylab("Daily Survival Rate") +
    ylim(0.7, 1) + 
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
    labs(title = "Clay-colored Sparrow",
         color = "Year"))

CCSP.plot
CCSPvor.plot
CCSPbhco.plot

ggsave(CCSP.plot,
       filename = "outputs/figs/CCSPbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(CCSPvor.plot,
       filename = "outputs/figs/CCSPvor.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(CCSPbhco.plot,
       filename = "outputs/figs/CCSPbhco.png",
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
