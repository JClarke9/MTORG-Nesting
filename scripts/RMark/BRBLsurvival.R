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
                               delete = FALSE)
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

coef(BRBL2.results$S.age)
coef(BRBL2.results$S.year)
coef(BRBL2.results$S.stage)

confint(BRBL2.results$S.age, level = 0.85)
confint(BRBL2.results$S.year, level = 0.85)
confint(BRBL2.results$S.stage, level = 0.85)


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
  
  # 4. DSR varies with the previous years grazing intensity
  S.grazedpTreat = list(formula = ~1 + Year + grazed + pTreat)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL3.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
BRBL3.results <- BRBL3.run()
BRBL3.results

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

BRBL.real <- as.data.frame(BRBL4.results$S.vor$results$real)
BRBL.real <- rownames_to_column(BRBL.real, var = "Group")
BRBL.real[,1] <- gsub("S g0", "", BRBL.real[,1])

BRBL.real <- BRBL.real |> 
  mutate(Year = case_when(
    startsWith(Group, "2021") ~ "2021",
    startsWith(Group, "2022") ~ "2022"
  ))

BRBL.avgDSR <- BRBL.real |> 
  group_by(Year) |> 
  summarize(estimate = mean(estimate))

# Plotting Beta Coefficients ----------------------------------------------

BRBL4.results$S.vor$results$real

BRBL.mod <- mark(BRBL.surv, nocc=max(BRBL.surv$LastChecked), 
                 model = "Nest", 
                 groups = "Year", 
                 model.parameters = list(S = list(formula =  ~1 + Year + Time + VOR)))

BRBL4.beta <- BRBL4.results$S.vor$results$beta
BRBL4.beta

BRBL.top <- rownames_to_column(BRBL4.beta, 
                               "Variable")
BRBL.top <- rename(BRBL.top, 
                   "Coefficient"="estimate")
str(BRBL.top)

BRBL.top[,1] <- gsub("S:", "", BRBL.top[,1])

BRBL.top$SEup <- BRBL.top$Coefficient + BRBL.top$se
BRBL.top$SElow <- BRBL.top$Coefficient - BRBL.top$se

BRBL.plot <- ggplot(BRBL.top, aes(x = Variable,
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
  ggtitle("BRBL Beta Coefficient") +
  theme(panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="white",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="white",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=12, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=12,                                              # change the size of the axis titles
                          colour = "black")) +                                    # change the color of the axis titles
  coord_flip()

# Creating Predictive Plots -----------------------------------------------

plotdata <- BRBL4.results$S.vor

BRBL.ddl <- make.design.data(BRBL.pr)
BRBL.ddl <- as.data.frame(BRBL.ddl)

minvor <- min(BRBL.surv$VOR)
maxvor <- max(BRBL.surv$VOR)
VOR.values = seq(from = minvor, 
                 to = maxvor, 
                 length = 100)

vor.pred <- covariate.predictions(plotdata,
                                  data=data.frame(VOR=VOR.values),
                                  indices = c(6, 41, 71,309, 344, 374))

Day5.21 <- which(vor.pred$estimates$par.index == 6)
Day40.21 <- which(vor.pred$estimates$par.index == 41)
Day70.21 <- which(vor.pred$estimates$par.index == 71)
Day5.22 <- which(vor.pred$estimates$par.index == 309)
Day40.22 <- which(vor.pred$estimates$par.index == 344)
Day70.22 <- which(vor.pred$estimates$par.index == 374)

vor.pred$estimates$Group <- NA
vor.pred$estimates$Group[Day5.21] <- "Early 2021"
vor.pred$estimates$Group[Day40.21] <- "Mid 2021"
vor.pred$estimates$Group[Day70.21] <- "Late 2021"
vor.pred$estimates$Group[Day5.22] <- "Early 2022"
vor.pred$estimates$Group[Day40.22] <- "Mid 2022"
vor.pred$estimates$Group[Day70.22] <- "Late 2022"
head(vor.pred$estimates)

vor.pred$estimates$Date <- NA
vor.pred$estimates$Date[Day5.21] <- "Early"
vor.pred$estimates$Date[Day40.21] <- "Mid"
vor.pred$estimates$Date[Day70.21] <- "Late"
vor.pred$estimates$Date[Day5.22] <- "Early"
vor.pred$estimates$Date[Day40.22] <- "Mid"
vor.pred$estimates$Date[Day70.22] <- "Late"
head(vor.pred$estimates)

vor.pred$estimates$Year <- NA
vor.pred$estimates$Year[Day5.21] <- "2021"
vor.pred$estimates$Year[Day40.21] <- "2021"
vor.pred$estimates$Year[Day70.21] <- "2021"
vor.pred$estimates$Year[Day5.22] <- "2022"
vor.pred$estimates$Year[Day40.22] <- "2022"
vor.pred$estimates$Year[Day70.22] <- "2022"
head(vor.pred$estimates)

BRBLvor.plot <- ggplot(transform(vor.pred$estimates,
                                 Year = factor(Year, levels=c("2021", "2022")),
                                 Date = factor(Date, levels=c("Early", "Mid", "Late"))),
                       aes(x = covdata, 
                           y = estimate,
                           groups = Group,
                           fill = Year)) +
  geom_line(size = 1.5,
            aes(linetype = Date,
                colour = Year)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.05)  +
  facet_grid(Date~.) +
  scale_colour_manual(values = c('#56B4E9',
                                          '#D55E00')) +
                                            scale_fill_manual(values = c('#56B4E9',
                                                                                  '#D55E00')) +
                                                                                    xlab("VOR") +
  ylab("Estimated DSR") +
  xlim(1 , 16) +
  ylim(0, 1) +
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
  labs(title = "Brewer's Blackbird")

BRBLvor.plot

julian.pred <- covariate.predictions(plotdata,
                                     data=data.frame(VOR=mean(VOR.values)),
                                     indices=c(1:76, 305:380))

Year21 <- which(julian.pred$estimates$par.index == c(1:76))
Year22 <- which(julian.pred$estimates$par.index == c(305:380))

julian.pred$estimates$Year[Year21] <- "2021"
julian.pred$estimates$Year[Year22] <- "2022"

julian.pred$estimates$Julian[Year21] <- c(1:76)
julian.pred$estimates$Julian[Year22] <- c(1:76)

BRBLjulian.plot <- ggplot(transform(julian.pred$estimates,
                                    Year = factor(Year, levels=c("2021", "2022"))), 
                          aes(x = Julian,
                              y = estimate,
                              fill = Year)) +
  geom_line(size = 1.5,
            aes(colour=Year)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.2) +
  scale_colour_manual(values = c('#56B4E9',
                                          '#D55E00')) +
                                            scale_fill_manual(values = c('#56B4E9',
                                                                                  '#D55E00')) +
                                                                                    xlab("Julian Day") +
  ylab("Estimated DSR") +
  ylim(0, 1.0) +
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
        legend.position = c(.1, .1),
        legend.box = "horizontal") +
  labs(title = "Brewer's Blackbird")

BRBL.plot
BRBLvor.plot
BRBLjulian.plot

ggsave(BRBL.plot,
       filename = "outputs/figs/BRBLbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BRBLvor.plot,
       filename = "outputs/figs/BRBLvor.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BRBLjulian.plot,
       filename = "outputs/figs/BRBLjulian.png",
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
