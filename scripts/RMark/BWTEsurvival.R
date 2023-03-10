# Loading libraries -------------------------------------------------------

library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)

windowsFonts(my_font = windowsFont("Gandhi Sans"))

# Data import -------------------------------------------------------------

nest <- read.csv("working/RMarknesting.csv", 
                 row.names=1)

# Subsetting data ---------------------------------------------------------

BWTE.surv <- dplyr::filter(nest,
                           Spec=="BWTE")                                         # select out only BWTE nests

BWTE.surv$AgeDay1 <- BWTE.surv$AgeFound - BWTE.surv$FirstFound + 1
BWTE.surv$Year <- as.factor(BWTE.surv$Year)
BWTE.surv$cTreat <- as.factor(BWTE.surv$cTreat)

# Creating nest stage variable --------------------------------------------

create.stage.var=function(data,agevar.name,stagevar.name,time.intervals,cutoff)
{
  nocc=length(time.intervals)
  age.mat=matrix(data[,agevar.name],nrow=dim(data)[1],ncol=nocc-1)
  age.mat=t(t(age.mat)+cumsum(c(0,time.intervals[1:(nocc-2)])))
  stage.mat=t(apply(age.mat,1,function(x) as.numeric(x<=cutoff)))
  stage.mat=data.frame(stage.mat)
  names(stage.mat)=paste(stagevar.name,1:(nocc-1),sep="")
  return(stage.mat)
}

x <- create.stage.var(BWTE.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(BWTE.surv$LastChecked)), 
                      23)

BWTE.surv <- bind_cols(BWTE.surv, x)

rm(list = ls()[!ls() %in% c("BWTE.surv")])

# Daily survival rate models ----------------------------------------------

BWTE.pr <- process.data(BWTE.surv,
                        nocc=max(BWTE.surv$LastChecked), 
                        groups = c("cTreat",
                                   "Year"),
                        model="Nest")

# candidate models for the null vs. treatment model
BWTE1.run <- function()
{
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  # 3. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies With year interacting with treatment (treatment identity varies by year)
  S.treatyear = list(formula = ~1 + cTreat:Year)
  
  BWTE.model.list = create.model.list("Nest")
  BWTE1.results = mark.wrapper(BWTE.model.list,
                               data = BWTE.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
BWTE1.results <- BWTE1.run()
BWTE1.results

BWTE1.results$S.cTreat$results$beta
BWTE1.results$S.treatyear$results$beta

#candidate model set for time trends
BWTE2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + cTreat + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + cTreat + Time + I(Time^2))
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  BWTE.model.list = create.model.list("Nest")
  BWTE2.results = mark.wrapper(BWTE.model.list,
                               data = BWTE.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
BWTE2.results <- BWTE2.run()
BWTE2.results

BWTE2.results$S.time$results$beta
BWTE2.results$S.quad$results$beta

#candidate model set for grazing
BWTE3.run <- function()
{
  # 1. DSR varies with BHCO number
  S.grazed = list(formula = ~1 + cTreat + Time + grazed)
  
  S.grazep = list(formula = ~1 + cTreat + Time + grazep)
  
  # 2. DSR varies with quadratic effect of date
  S.time = list(formula = ~1 + cTreat + Time)
  BWTE.model.list = create.model.list("Nest")
  BWTE3.results = mark.wrapper(BWTE.model.list,
                               data = BWTE.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
BWTE3.results <- BWTE3.run()
BWTE3.results

BWTE3.results$S.time$results$beta
BWTE3.results$S.grazed$results$beta
BWTE3.results$S.grazep$results$beta

#candidate model set for parasitism and nest stage
BWTE4.run <- function()
{
  # 1. DSR varies with KBG
  S.kbg = list(formula =  ~1 + cTreat + Time + KBG)
  
  # 2. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.smooth.brome = list(formula = ~1 + cTreat + Time + SmoothB)
  
  # 3. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + cTreat + Time + Litter)
  
  # 4. DSR varies with Bare
  S.bare = list(formula =  ~1 + cTreat + Time + Bare)
  
  # 5. DSR varies with Forb
  S.forb = list(formula =  ~1 + cTreat + Time + Forb)
  
  # 6. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + cTreat + Time + Grasslike)
  
  # 7. DSR varies with Woody
  S.woody = list(formula =  ~1 + cTreat + Time + Woody)
  
  # 8. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + cTreat + Time + LitterD)
  
  # 9. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + cTreat + Time + Veg.Height)
  
  # 10. DSR varies with VOR
  S.vor = list(formula =  ~1 + cTreat + Time + VOR)
  
  # 11. DSR varies with litter and Veg Height
  S.litheight = list(formula = ~1 + cTreat + Time + Litter + Veg.Height)
  
  # 12. DSR varies with woody and litter depth
  S.woodylitdep = list(formula = ~1 + cTreat + Time + Woody + LitterD)
  
  # 13. DSR varies with KBG and litter depth
  S.kbglitdep = list(formula = ~1 + cTreat + Time + KBG + LitterD)
  
  # 14. DSR varies with KBG and Veg.Height
  S.kbgheight = list(formula = ~1 + cTreat + Time + KBG + Veg.Height)
  
  # 15. DSR varies with woody and Veg Height
  S.woodyheight = list(formula = ~1 + cTreat + Time + Woody + Veg.Height)
  
  # 1. DSR varies with BHCO number
  S.time = list(formula =  ~1 + cTreat + Time)
  
  BWTE.model.list = create.model.list("Nest")
  BWTE4.results = mark.wrapper(BWTE.model.list,
                               data = BWTE.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

#results of candidate model set 
#all results
BWTE4.results <- BWTE4.run()
BWTE4.results

BWTE4.results$S.litdep$results$beta
BWTE4.results$S.woodylitdep$results$beta
BWTE4.results$S.kbglitdep$results$beta

BWTE4.results$S.bare$results$real

BWTE.real <- as.data.frame(BWTE4.results$S.bare$results$real)
BWTE.real <- rownames_to_column(BWTE.real, var = "Group")
BWTE.real[,1] <- gsub("S g", "", BWTE.real[,1])

BWTE.real <- BWTE.real |> 
  mutate(Treat = case_when(
    startsWith(Group, "0") ~ "Rest",
    startsWith(Group, "39") ~ "Moderate",
    startsWith(Group, "49") ~ "Full",
    startsWith(Group, "68") ~ "Heavy"
  ))

BWTE.avgDSR <- BWTE.real |> 
  group_by(Treat) |> 
  summarize(estimate = mean(estimate))

# Plotting beta coefficients ----------------------------------------------

BWTE.mod <- mark(BWTE.surv, 
                 nocc=max(BWTE.surv$LastChecked), 
                 model = "Nest", 
                 groups = "cTreat", 
                 model.parameters = list(S = list(formula =  ~1 + cTreat + Time + LitterD)))

BWTE.mod$results$beta

BWTE4.beta <- BWTE4.results$S.litdep$results$beta
BWTE4.beta

BWTE.top <- rownames_to_column(BWTE4.beta, 
                               "Variable")
BWTE.top <- rename(BWTE.top, 
                   "Coefficient"="estimate")

BWTE.top[,1] <- gsub("S:", "", BWTE.top[,1])

str(BWTE.top)

BWTE.top$SEup <- BWTE.top$Coefficient + BWTE.top$se
BWTE.top$SElow <- BWTE.top$Coefficient - BWTE.top$se

BWTE.plot <- ggplot(BWTE.top, aes(x = Variable,
                                  y = Coefficient)) +
  geom_hline(yintercept = 0,
             colour = gray(1/2), lty = 2) +
  geom_point(aes(x = Variable,
                 y = Coefficient),
             size = 5) +
  geom_errorbar(aes(x = Variable,
                    ymin = lcl,
                    ymax = ucl),
                width = .5,
                linewidth = 1) +
  ggtitle("BWTE Beta Coefficient") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white",
                                        colour = NA),
        plot.background = element_rect(fill="white",
                                       colour=NA),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=12, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        text=element_text(size=12,
                          colour = "black")) +
  coord_flip()

# Creating predictive plots -----------------------------------------------

plotdata <- BWTE4.results$S.litdep

BWTE.ddl <- make.design.data(BWTE.pr)
BWTE.ddl <- as.data.frame(BWTE.ddl)

minlitdep <- min(BWTE.surv$LitterD)
maxlitdep <- max(BWTE.surv$LitterD)
LitterD.values = seq(from = minlitdep, 
                     to = maxlitdep, 
                     length = 100)

litdep.pred <- covariate.predictions(plotdata,
                                     data=data.frame(LitterD=LitterD.values),
                                     indices = c(6, 31, 58, 69, 94, 121, 132, 157, 184, 195, 220, 247))

rows39.5 <- which(litdep.pred$estimates$par.index == 6)
rows39.30 <- which(litdep.pred$estimates$par.index == 31)
rows39.57 <- which(litdep.pred$estimates$par.index == 58)
rows49.5 <- which(litdep.pred$estimates$par.index == 69)
rows49.30 <- which(litdep.pred$estimates$par.index == 94)
rows49.57 <- which(litdep.pred$estimates$par.index == 121)
rows68.5 <- which(litdep.pred$estimates$par.index == 132)
rows68.30 <- which(litdep.pred$estimates$par.index == 157)
rows68.57 <- which(litdep.pred$estimates$par.index == 184)
rows0.5 <- which(litdep.pred$estimates$par.index == 195)
rows0.30 <- which(litdep.pred$estimates$par.index == 220)
rows0.57 <- which(litdep.pred$estimates$par.index == 247)

litdep.pred$estimates$Group <- NA
litdep.pred$estimates$Group[rows39.5] <- "Early Moderate"
litdep.pred$estimates$Group[rows39.30] <- "Mid Moderate"
litdep.pred$estimates$Group[rows39.57] <- "Late Moderate"
litdep.pred$estimates$Group[rows49.5] <- "Early Full"
litdep.pred$estimates$Group[rows49.30] <- "Mid Full"
litdep.pred$estimates$Group[rows49.57] <- "Late Full"
litdep.pred$estimates$Group[rows68.5] <- "Early Heavy"
litdep.pred$estimates$Group[rows68.30] <- "Mid Heavy"
litdep.pred$estimates$Group[rows68.57] <- "Late Heavy"
litdep.pred$estimates$Group[rows0.5] <- "Early Rest"
litdep.pred$estimates$Group[rows0.30] <- "Mid Rest"
litdep.pred$estimates$Group[rows0.57] <- "Late Rest"
head(litdep.pred$estimates)

litdep.pred$estimates$Date <- NA
litdep.pred$estimates$Date[rows39.5] <- "Early"
litdep.pred$estimates$Date[rows39.30] <- "Mid"
litdep.pred$estimates$Date[rows39.57] <- "Late"
litdep.pred$estimates$Date[rows49.5] <- "Early"
litdep.pred$estimates$Date[rows49.30] <- "Mid"
litdep.pred$estimates$Date[rows49.57] <- "Late"
litdep.pred$estimates$Date[rows68.5] <- "Early"
litdep.pred$estimates$Date[rows68.30] <- "Mid"
litdep.pred$estimates$Date[rows68.57] <- "Late"
litdep.pred$estimates$Date[rows0.5] <- "Early"
litdep.pred$estimates$Date[rows0.30] <- "Mid"
litdep.pred$estimates$Date[rows0.57] <- "Late"
head(litdep.pred$estimates)

litdep.pred$estimates$Treat <- NA
litdep.pred$estimates$Treat[rows39.5] <- "Moderate"
litdep.pred$estimates$Treat[rows39.30] <- "Moderate"
litdep.pred$estimates$Treat[rows39.57] <- "Moderate"
litdep.pred$estimates$Treat[rows49.5] <- "Full"
litdep.pred$estimates$Treat[rows49.30] <- "Full"
litdep.pred$estimates$Treat[rows49.57] <- "Full"
litdep.pred$estimates$Treat[rows68.5] <- "Heavy"
litdep.pred$estimates$Treat[rows68.30] <- "Heavy"
litdep.pred$estimates$Treat[rows68.57] <- "Heavy"
litdep.pred$estimates$Treat[rows0.5] <- "Rest"
litdep.pred$estimates$Treat[rows0.30] <- "Rest"
litdep.pred$estimates$Treat[rows0.57] <- "Rest"
head(litdep.pred$estimates)

litdep.pred.rest <- litdep.pred[c(6,31,58)]

BWTElitdep.plot <- ggplot(transform(litdep.pred$estimates,
                                    Date = factor(Date,
                                                  levels=c("Early", "Mid", "Late")),
                                    Treat = factor(Treat,
                                                   levels=c("Rest", "Moderate", "Full", "Heavy"))), 
                          aes(x = covdata, 
                              y = estimate,
                              groups = Group,
                              fill = Treat)) +
  geom_line(size = 1.5,
            aes(colour = Treat)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.05) +
  facet_grid(Date~.) +
  scale_colour_manual(values = c('#A2A4A2',
                                          'lightgoldenrod2',
                                          '#D4A634',
                                          '#717F5B')) +
                                            scale_fill_manual(values = c('#A2A4A2',
                                                                                  'lightgoldenrod2',
                                                                                  '#D4A634',
                                                                                  '#717F5B')) +
                                                                                    xlab("Litter Depth (mm)") +
  ylab("Estimated DSR") +
  ylim(0, 1.0) +
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=16,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill=NA,                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill=NA,                      # make the outer background transparent
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
  labs(title = "Blue-winged Teal")

julian.pred <- covariate.predictions(plotdata,
                                     data=data.frame(LitterD=mean(LitterD.values)),
                                     indices = c(1:252))

rows0 <- which(julian.pred$estimates$par.index == c(190:252))
rows39 <- which(julian.pred$estimates$par.index == c(1:63))
rows49 <- which(julian.pred$estimates$par.index == c(64:126))
rows68 <- which(julian.pred$estimates$par.index == c(127:189))

julian.pred$estimates$Treat <- NA
julian.pred$estimates$Treat[rows0] <- "Rest"
julian.pred$estimates$Treat[rows39] <- "Moderate"
julian.pred$estimates$Treat[rows49] <- "Full"
julian.pred$estimates$Treat[rows68] <- "Heavy"
head(julian.pred$estimates)

julian.pred$estimates$Julian[rows0] <- c(1:63)
julian.pred$estimates$Julian[rows39] <- c(1:63)
julian.pred$estimates$Julian[rows49] <- c(1:63)
julian.pred$estimates$Julian[rows68] <- c(1:63)

BWTEjulian.plot <- ggplot(transform(julian.pred$estimates,
                                    Treat = factor(Treat,
                                                   levels=c("Rest", "Moderate", "Full", "Heavy"))), 
                          aes(x = Julian,
                              y = estimate,
                              fill = Treat)) +
  geom_line(size = 1.5,
            aes(colour = Treat)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  scale_colour_manual(values = c('#A2A4A2',
                                          'lightgoldenrod2',
                                          '#D4A634',
                                          '#717F5B')) +
                                            scale_fill_manual(values = c('#A2A4A2',
                                                                                  'lightgoldenrod2',
                                                                                  '#D4A634',
                                                                                  '#717F5B')) +
                                                                                    xlab("Julian Day") +
  ylab("Estimated DSR") +
  ylim(0, 1.0) +
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=16,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill=NA,                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill=NA,                      # make the outer background transparent
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
  labs(title = "Blue-winged Teal")

BWTEjulian.plot

BWTE.plot
BWTElitdep.plot
BWTEjulian.plot

ggsave(BWTE.plot,
       filename = "~/Git/NDSU/RMARK/Figures/BWTEbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BWTElitdep.plot,
       filename = "outputs/figs/BWTElitdep.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(BWTEjulian.plot,
       filename = "outputs/figs/BWTEjulian.png",
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
