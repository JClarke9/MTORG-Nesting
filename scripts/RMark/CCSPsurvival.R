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

CCSP.surv <- filter(nest, 
                    Spec=="CCSP")                                         # select out only CCSP nests

CCSP.surv$AgeDay1 <- CCSP.surv$AgeFound - CCSP.surv$FirstFound + 1
CCSP.surv$Year <- as.factor(CCSP.surv$Year)
CCSP.surv$cTreat <- as.factor(CCSP.surv$cTreat)

# Creating stage variable -------------------------------------------------

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
                        groups = c("cTreat", 
                                   "Year"),
                        model="Nest")

# candidate models for the null vs. treatment model
CCSP1.run <- function()
{
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  # 3. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies With year interacting with treatment (treatment identity varies by year)
  S.treatyear = list(formula = ~1 + cTreat:Year)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP1.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
CCSP1.results <- CCSP1.run()
CCSP1.results

CCSP1.results$S.year$results$beta
CCSP1.results$S.treatyear$results$beta

#candidate model set for time trends
CCSP2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Year + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Year+ Time + I(Time^2))
  
  # 4. DSR varies across stages
  S.stage = list(formula = ~1 + Year + Incub)
  
  # 2. DSR varies with current treatment
  S.year = list(formula = ~1 + Year)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP2.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
CCSP2.results <- CCSP2.run()
CCSP2.results

CCSP2.results$S.time$results$beta
CCSP2.results$S.stage$results$beta
CCSP2.results$S.year$results$beta
CCSP2.results$S.quad$results$beta

#candidate model set for parasitism and nest stage
CCSP3.run <- function()
{
  # 1. DSR varies with BHCO number
  S.bhco = list(formula = ~1 + Year + Time + BHCONum)
  
  S.grazed = list(formula = ~1 + Year + Time + grazed)
  
  S.bhcop = list(formula = ~1 + Year + Time + BHCOpres)
  
  S.grazep = list(formula = ~1 + Year + Time + grazep)
  
  # 2. DSR varies with quadratic effect of date
  S.time = list(formula = ~1 + Year + Time)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP3.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
CCSP3.results <- CCSP3.run()
CCSP3.results

CCSP3.results$S.bhco$results$beta
CCSP3.results$S.time$results$beta

#candidate model set for vegetation data
CCSP4.run <- function()
{
  # 1. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Year + Time + BHCOpres + KBG)
  
  # 2. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.smooth.brome = list(formula = ~1 + Year + Time + BHCOpres + SmoothB)
  
  # 3. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Year + Time + BHCOpres + Litter)
  
  # 4. DSR varies with Bare
  S.bare = list(formula =  ~1 + Year + Time + BHCOpres + Bare)
  
  # 5. DSR varies with Forb
  S.forb = list(formula =  ~1 + Year + Time + BHCOpres + Forb)
  
  # 6. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Year + Time + BHCOpres + Grasslike)
  
  # 7. DSR varies with Woody
  S.woody = list(formula =  ~1 + Year + Time + BHCOpres + Woody)
  
  # 8. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + Time + BHCOpres + LitterD)
  
  # 9. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + Time + BHCOpres + Veg.Height)
  
  # 10. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + Time + BHCOpres + VOR)
  
  # 11. DSR varies with litter and Veg Height
  S.litheight = list(formula = ~1 + Year + Time + BHCOpres + Litter + Veg.Height)
  
  # 12. DSR varies with woody and litter depth
  S.woodylitdep = list(formula = ~1 + Year + Time + BHCOpres + Woody + LitterD)
  
  # 13. DSR varies with KBG and litter depth
  S.kbglitdep = list(formula = ~1 + Year + Time + BHCOpres + KBG + LitterD)
  
  # 14. DSR varies with KBG and Veg.Height
  S.kbgheight = list(formula = ~1 + Year + Time + BHCOpres + KBG + Veg.Height)
  
  # 15. DSR varies with woody and Veg Height
  S.woodyheight = list(formula = ~1 + Year + Time + BHCOpres + Woody + Veg.Height)
  
  # 1. DSR varies with BHCO number
  S.bhcop = list(formula =  ~1 + Year + Time + BHCOpres)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP4.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
CCSP4.results <- CCSP4.run()
CCSP4.results

CCSP4.results$S.height$results$beta
CCSP4.results$S.litheight$results$beta
CCSP4.results$S.woodyheight$results$beta
CCSP4.results$S.kbgheight$results$beta

CCSP.real <- as.data.frame(CCSP4.results$S.height$results$real)
CCSP.real <- rownames_to_column(CCSP.real, var = "Group")
CCSP.real[,1] <- gsub("S g0", "", CCSP.real[,1])

CCSP.real <- CCSP.real |> 
  mutate(Year = case_when(
    startsWith(Group, "2021") ~ "2021",
    startsWith(Group, "2022") ~ "2022"
  ))

CCSP.avgDSR <- CCSP.real |> 
  group_by(Year) |> 
  summarize(estimate = mean(estimate))

# Plotting beta coefficients ----------------------------------------------

CCSP.mod <- mark(CCSP.surv, 
                 nocc=max(CCSP.surv$LastChecked), 
                 model = "Nest", 
                 groups = c("Year"), 
                 model.parameters = list(S = list(formula =  ~1 + Year + Time + BHCOpres + Veg.Height)))

CCSP4.results$S.height$results$real
CCSP.mod$results$beta

CCSP4.beta <- CCSP4.results$S.height$results$beta
CCSP4.beta

CCSP.top <- rownames_to_column(CCSP4.beta, 
                               "Variable")
CCSP.top <- rename(CCSP.top, 
                   "Coefficient"="estimate")

CCSP.top[,1] <- gsub("S:", "", CCSP.top[,1])
str(CCSP.top)

CCSP.top$SEup <- CCSP.top$Coefficient + CCSP.top$se
CCSP.top$SElow <- CCSP.top$Coefficient - CCSP.top$se

CCSP.plot <- ggplot(CCSP.top, 
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
                size = 1) +
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
  coord_flip()

# Creating predictive plots -----------------------------------------------

CCSP.ddl <- make.design.data(CCSP.pr)
CCSP.ddl <- as.data.frame(CCSP.ddl)

plotdata <- CCSP.mod <- CCSP4.results$S.height

minheight <- min(CCSP.surv$Veg.Height)
maxheight <- max(CCSP.surv$Veg.Height)
Veg.Heightvalues = seq(from = minheight,
                       to = maxheight,
                       length = 100)

minbhco <- min(CCSP.surv$BHCONum)
maxbhco <- max(CCSP.surv$BHCONum)
BHCONum.values = seq(from = minbhco,
                     to = maxbhco,
                     length = 100)

height.pred <- covariate.predictions(plotdata,
                                     data=data.frame(Veg.Height=Veg.Heightvalues,
                                                     BHCONum=mean(BHCONum.values)),
                                     indices = c(6, 41, 81, 346, 381, 421))

Day5.21 <- which(height.pred$estimates$par.index == 6)
Day40.21 <- which(height.pred$estimates$par.index == 41)
Day80.21 <- which(height.pred$estimates$par.index == 81)
Day5.22 <- which(height.pred$estimates$par.index == 346)
Day40.22 <- which(height.pred$estimates$par.index == 381)
Day80.22 <- which(height.pred$estimates$par.index == 421)

height.pred$estimates$Group <- NA
height.pred$estimates$Group[Day5.21] <- "Early 2021"
height.pred$estimates$Group[Day40.21] <- "Mid 2021"
height.pred$estimates$Group[Day80.21] <- "Late 2021"
height.pred$estimates$Group[Day5.22] <- "Early 2022"
height.pred$estimates$Group[Day40.22] <- "Mid 2022"
height.pred$estimates$Group[Day80.22] <- "Late 2022"
head(height.pred$estimates)

height.pred$estimates$Date <- NA
height.pred$estimates$Date[Day5.21] <- "Early"
height.pred$estimates$Date[Day40.21] <- "Mid"
height.pred$estimates$Date[Day80.21] <- "Late"
height.pred$estimates$Date[Day5.22] <- "Early"
height.pred$estimates$Date[Day40.22] <- "Mid"
height.pred$estimates$Date[Day80.22] <- "Late"
head(height.pred$estimates)

height.pred$estimates$Year <- NA
height.pred$estimates$Year[Day5.21] <- "2021"
height.pred$estimates$Year[Day40.21] <- "2021"
height.pred$estimates$Year[Day80.21] <- "2021"
height.pred$estimates$Year[Day5.22] <- "2022"
height.pred$estimates$Year[Day40.22] <- "2022"
height.pred$estimates$Year[Day80.22] <- "2022"
head(height.pred$estimates)

CCSPheight.plot <- ggplot(transform(height.pred$estimates,
                                    Year = factor(Year, levels=c("2021", "2022")),
                                    Date = factor(Date, levels=c("Early", "Mid", "Late"))), 
                          aes(x = Veg.Height, 
                              y = estimate,
                              groups = Group,
                              fill = Year)) +
  geom_line(size = 1.5,
            aes(linetype = Date,
                colour = Year)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  facet_grid(Date~.) +
  scale_colour_manual(values = c('#56B4E9',
                                          '#D55E00')) +
                                            scale_fill_manual(values = c('#56B4E9',
                                                                                  '#D55E00')) +
                                                                                    xlab("Veg Height (mm)") +
  ylab("Estimated DSR") +
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
  labs(title = "Clay-colored Sparrow")

CCSPheight.plot

bhco.pred <- covariate.predictions(plotdata,
                                   data=data.frame(Veg.Height=mean(Veg.Heightvalues),
                                                   BHCONum=BHCONum.values),
                                   indices = c(6, 41, 81, 346, 481, 421))

Day5.21 <- which(bhco.pred$estimates$par.index == 6)
Day40.21 <- which(bhco.pred$estimates$par.index == 41)
Day80.21 <- which(bhco.pred$estimates$par.index == 81)
Day5.22 <- which(bhco.pred$estimates$par.index == 346)
Day40.22 <- which(bhco.pred$estimates$par.index == 481)
Day80.22 <- which(bhco.pred$estimates$par.index == 421)

bhco.pred$estimates$Group <- NA
bhco.pred$estimates$Group[Day5.21] <- "Early 2021"
bhco.pred$estimates$Group[Day40.21] <- "Mid 2021"
bhco.pred$estimates$Group[Day80.21] <- "Late 2021"
bhco.pred$estimates$Group[Day5.22] <- "Early 2022"
bhco.pred$estimates$Group[Day40.22] <- "Mid 2022"
bhco.pred$estimates$Group[Day80.22] <- "Late 2022"
head(height.pred$estimates)

bhco.pred$estimates$Date <- NA
bhco.pred$estimates$Date[Day5.21] <- "Early"
bhco.pred$estimates$Date[Day40.21] <- "Mid"
bhco.pred$estimates$Date[Day80.21] <- "Late"
bhco.pred$estimates$Date[Day5.22] <- "Early"
bhco.pred$estimates$Date[Day40.22] <- "Mid"
bhco.pred$estimates$Date[Day80.22] <- "Late"
head(bhco.pred$estimates)

bhco.pred$estimates$Year <- NA
bhco.pred$estimates$Year[Day5.21] <- "2021"
bhco.pred$estimates$Year[Day40.21] <- "2021"
bhco.pred$estimates$Year[Day80.21] <- "2021"
bhco.pred$estimates$Year[Day5.22] <- "2022"
bhco.pred$estimates$Year[Day40.22] <- "2022"
bhco.pred$estimates$Year[Day80.22] <- "2022"
head(bhco.pred$estimates)

CCSPbhco.plot <- ggplot(transform(bhco.pred$estimates,
                                  Year = factor(Year, levels=c("2021", "2022")),
                                  Date = factor(Date, levels=c("Early", "Mid", "Late"))),
                        aes(x = BHCONum,
                            y = estimate,
                            groups = Group,
                            fill = Year)) +
  geom_line(size = 1.5,
            aes(linetype = Date,
                colour=Year)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.1) +
  facet_grid(Date~.) +
  scale_colour_manual(values = c('#56B4E9',
                                          '#D55E00')) +
                                            scale_fill_manual(values = c('#56B4E9',
                                                                                  '#D55E00')) +
                                                                                    xlab("Julian Day") +
  ylab("Estimated DSR") +
  ylim(.6, 1.0) +
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
        legend.position = c(.15, .1),
        legend.box = "horizontal") +
  labs(title = "Clay-colored Sparrow")

julian.pred <- covariate.predictions(plotdata,
                                     data=data.frame(Veg.Height=mean(Veg.Heightvalues),
                                                     BHCONum=mean(BHCONum.values)),
                                     indices=c(1:85, 341:425))

Year21 <- which(julian.pred$estimates$par.index == c(1:85))
Year22 <- which(julian.pred$estimates$par.index == c(341:425))

julian.pred$estimates$Year[Year21] <- "2021"
julian.pred$estimates$Year[Year22] <- "2022"

julian.pred$estimates$Julian[Year21] <- c(1:85)
julian.pred$estimates$Julian[Year22] <- c(1:85)

CCSPjulian.plot <- ggplot(transform(julian.pred$estimates,
                                    Year = factor(Year, levels=c("2021", "2022"))),
                          aes(x = Julian,
                              y = estimate,
                              fill = Year)) +
  geom_line(size = 1.5,
            aes(colour=Year)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.3) +
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
  labs(title = "Clay-colored Sparrow")

CCSP.plot
CCSPheight.plot
CCSPbhco.plot
CCSPjulian.plot

ggsave(CCSP.plot,
       filename = "outputs/figs/CCSPbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(CCSPheight.plot,
       filename = "outputs/figs/CCSPheight.png",
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

ggsave(CCSPjulian.plot,
       filename = "outputs/figs/CCSPjulian.png",
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