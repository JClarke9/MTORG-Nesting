# Loading Libraries -------------------------------------------------------

library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)

windowsFonts(my_font = windowsFont("Gandhi Sans"))

# Data import -------------------------------------------------------------

nest <- read.csv("working/RMarknesting.csv", 
                 row.names=1)

# Subsetting Data ---------------------------------------------------------

BRBL.surv <- filter(nest, 
                    Spec=="BRBL")                                         # select out only BRBL nests

BRBL.surv$AgeDay1 <- BRBL.surv$AgeFound - BRBL.surv$FirstFound + 1
BRBL.surv$Year <- as.factor(BRBL.surv$Year)
BRBL.surv$cTreat <- as.factor(BRBL.surv$cTreat)

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

x <- create.stage.var(BRBL.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(BRBL.surv$LastChecked)), 
                      12)

BRBL.surv <- bind_cols(BRBL.surv, x)

#remove everything but the final dataframe
rm(list = ls()[!ls() %in%  c("BRBL.surv")])

# Daily survival rate models ----------------------------------------------

BRBL.pr <- process.data(BRBL.surv,
                        nocc=max(BRBL.surv$LastChecked),
                        groups = c("cTreat", 
                                   "Year"),
                        model="Nest")

# candidate models for the null vs. treatment model
BRBL1.run <- function()
{
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  # 3. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies With year interacting with treatment (treatment identity varies by year)
  S.treatyear = list(formula = ~1 + cTreat:Year)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL1.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
BRBL1.results <- BRBL1.run()
BRBL1.results

BRBL1.results$S.year$results$beta
BRBL1.results$S.Dot$results$beta

#candidate model set for time trends
BRBL2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Year + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Year+ Time + I(Time^2))
  
  # 4. DSR varies across stages
  S.stage = list(formula = ~1 + Year + Incub)
  
  # 2. DSR varies with current treatment
  S.year = list(formula = ~1 + Year)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL2.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
BRBL2.results <- BRBL2.run()
BRBL2.results

BRBL2.results$S.time$results$beta
BRBL2.results$S.year$results$beta
BRBL2.results$S.quad$results$beta

#candidate model set for parasitism and grazing
BRBL3.run <- function()
{
  # 1. DSR varies with BHCO number
  S.bhco = list(formula = ~1 + Year + Time + BHCONum)
  
  S.grazed = list(formula = ~1 + Year + Time + grazed)
  
  S.bhcop = list(formula = ~1 + Year + Time + BHCOpres)
  
  S.grazep = list(formula = ~1 + Year + Time + grazep)
  
  # 2. DSR varies with quadratic effect of date
  S.time = list(formula = ~1 + Year + Time)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL3.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
BRBL3.results <- BRBL3.run()
BRBL3.results

BRBL3.results$S.time$results$beta
BRBL3.results$S.bhco$results$beta

#candidate model set for parasitism and nest stage
BRBL4.run <- function()
{
  # 1. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Year + Time + KBG)
  
  # 2. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.smooth.brome = list(formula = ~1 + Year + Time + SmoothB)
  
  # 3. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Year + Time + Litter)
  
  # 4. DSR varies with Bare
  S.bare = list(formula =  ~1 + Year + Time + Bare)
  
  # 5. DSR varies with Forb
  S.forb = list(formula =  ~1 + Year + Time + Forb)
  
  # 6. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Year + Time + Grasslike)
  
  # 7. DSR varies with Woody
  S.woody = list(formula =  ~1 + Year + Time + Woody)
  
  # 8. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Year + Time + LitterD)
  
  # 9. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Year + Time + Veg.Height)
  
  # 10. DSR varies with VOR
  S.vor = list(formula =  ~1 + Year + Time + VOR)
  
  # 11. DSR varies with litter and Veg Height
  S.litheight = list(formula = ~1 + Year + Time + Litter + Veg.Height)
  
  # 12. DSR varies with woody and litter depth
  S.woodylitdep = list(formula = ~1 + Year + Time + Woody + LitterD)
  
  # 13. DSR varies with KBG and litter depth
  S.kbglitdep = list(formula = ~1 + Year + Time + KBG + LitterD)
  
  # 14. DSR varies with KBG and Veg.Height
  S.kbgheight = list(formula = ~1 + Year + Time + KBG + Veg.Height)
  
  # 15. DSR varies with woody and Veg Height
  S.woodyheight = list(formula = ~1 + Year + Time + Woody + Veg.Height)
  
  # 2. DSR varies with quadratic effect of date
  S.time = list(formula = ~1 + Year + Time)
  
  BRBL.model.list = create.model.list("Nest")
  BRBL4.results = mark.wrapper(BRBL.model.list,
                               data = BRBL.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

#results of candidate model set 
#all results
BRBL4.results <- BRBL4.run()
BRBL4.results
names(BRBL4.results)

BRBL4.results$S.vor$results$beta
BRBL4.results$S.forb$results$beta

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
