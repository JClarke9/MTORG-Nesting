# Loading libraries -------------------------------------------------------

library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)

windowsFonts(my_font = windowsFont("Gandhi Sans"))

# Data  import ------------------------------------------------------------

nest <- read.csv("working/RMarknesting.csv", 
                 row.names=1)

# Subsetting data ---------------------------------------------------------

RWBL.surv <- filter(nest, 
                    Spec=="RWBL")                                         # select out only RWBL nests

RWBL.surv$AgeDay1 <- RWBL.surv$AgeFound - RWBL.surv$FirstFound + 1
RWBL.surv$Year <- as.factor(RWBL.surv$Year)
RWBL.surv$cTreat <- as.factor(RWBL.surv$cTreat)

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

x <- create.stage.var(RWBL.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(RWBL.surv$LastChecked)), 
                      12)

RWBL.surv <- bind_cols(RWBL.surv, x)

rm(list = ls()[!ls() %in% "RWBL.surv"])

# Daily survival rate models ----------------------------------------------

RWBL.pr <- process.data(RWBL.surv,
                        nocc=max(RWBL.surv$LastChecked),
                        groups = c("cTreat", "Year"),
                        model="Nest")

# candidate models for the null vs. treatment model
RWBL1.run <- function()
{
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  # 3. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies With year interacting with treatment (treatment identity varies by year)
  S.treatyear = list(formula = ~1 + cTreat:Year)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL1.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
RWBL1.results <- RWBL1.run()
RWBL1.results

RWBL1.results$S.cTreat$results$beta
RWBL1.results$S.treatyear$results$beta

#candidate model set for time trends
RWBL2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + cTreat + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + cTreat + Time + I(Time^2))
  
  # 4. DSR varies across stages
  S.stage = list(formula = ~1 + cTreat + Incub)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL2.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
RWBL2.results <- RWBL2.run()
RWBL2.results

RWBL2.results$S.time$results$beta
RWBL2.results$S.cTreat$results$beta
RWBL2.results$S.quad$results$beta

#candidate model set for parasitism and grazing
RWBL3.run <- function()
{
  # 1. DSR varies with BHCO number
  S.bhco = list(formula = ~1 + cTreat + Time + BHCONum)
  
  S.grazed = list(formula = ~1 + cTreat + Time + grazed)
  
  S.bhcop = list(formula = ~1 + cTreat + Time + BHCOpres)
  
  S.grazep = list(formula = ~1 + cTreat + Time + grazep)
  
  # 2. DSR varies with quadratic effect of date
  S.time = list(formula = ~1 + cTreat + Time)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL3.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
RWBL3.results <- RWBL3.run()
RWBL3.results

RWBL3.results$S.time$results$beta
RWBL3.results$S.bhcop$results$beta
RWBL3.results$S.grazed$results$beta
RWBL3.results$S.bhco$results$beta
RWBL3.results$S.grazep$results$beta

#candidate model set for parasitism and nest stage
RWBL4.run <- function()
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
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~cTreat)
  
  RWBL.model.list = create.model.list("Nest")
  RWBL4.results = mark.wrapper(RWBL.model.list,
                               data = RWBL.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

#results of candidate model set 
#all results
RWBL4.results <- RWBL4.run()
RWBL4.results
names(RWBL4.results)

RWBL4.results$S.bare$results$beta
RWBL4.results$S.kbg$results$beta
RWBL4.results$S.litdep$results$beta

RWBL4.results$S.bare$results$real

RWBL.real <- as.data.frame(RWBL4.results$S.bare$results$real)
RWBL.real <- rownames_to_column(RWBL.real, var = "Group")
RWBL.real[,1] <- gsub("S g", "", RWBL.real[,1])

RWBL.real <- RWBL.real |> 
  mutate(Treat = case_when(
    startsWith(Group, "0") ~ "Rest",
    startsWith(Group, "39") ~ "Moderate",
    startsWith(Group, "49") ~ "Full",
    startsWith(Group, "68") ~ "Heavy"
  ))

RWBL.real <- RWBL.real |> 
  mutate(Year = case_when(
    grepl("2021", Group) ~ "2021",
    grepl("2022", Group) ~ "2022"
  ))

RWBL.avgDSR <- RWBL.real |> 
  group_by(Treat,
           Year) |> 
  summarize(estimate = mean(estimate))

# Plotting beta coefficients ----------------------------------------------

RWBL.mod <- mark(RWBL.surv, 
                 nocc=max(RWBL.surv$LastChecked), 
                 model = "Nest", 
                 groups ="cTreat", 
                 model.parameters = list(S = list(formula =  ~1 + cTreat + Time + Bare)))

RWBL4.beta <- RWBL4.results$S.bare$results$beta
RWBL4.beta

RWBL.top <- rownames_to_column(RWBL4.beta, 
                               "Variable")
RWBL.top <- rename(RWBL.top, 
                   "Coefficient"="estimate")

RWBL.top[,1] <- gsub("S:", "", RWBL.top[,1])
str(RWBL.top)

RWBL.top$SEup <- RWBL.top$Coefficient + RWBL.top$se
RWBL.top$SElow <- RWBL.top$Coefficient - RWBL.top$se

RWBL.plot <- ggplot(RWBL.top, aes(x = Variable, 
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
                size = 1) +
  ggtitle("RWBL Beta Coefficient") +
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

RWBL.plot

# Creating predictive plots -----------------------------------------------

plotdata <- RWBL4.results$S.bare

RWBL.ddl <- make.design.data(RWBL.pr)
RWBL.ddl <- as.data.frame(RWBL.ddl)

minbare <- min(RWBL.surv$Bare)
maxbare <- max(RWBL.surv$Bare)
bare.values = seq(from = minbare, 
                  to = maxbare, 
                  length = 100)

bare.pred <- covariate.predictions(plotdata,
                                   data=data.frame(Bare=bare.values),
                                   indices = c(6, 31, 66, 77, 102, 137, 148, 173, 208, 361, 386, 421))

Day5.39 <- which(bare.pred$estimates$par.index == 6)
Day30.39 <- which(bare.pred$estimates$par.index == 31)
Day65.39 <- which(bare.pred$estimates$par.index == 66)
Day5.49 <- which(bare.pred$estimates$par.index == 77)
Day30.49 <- which(bare.pred$estimates$par.index == 102)
Day65.49 <- which(bare.pred$estimates$par.index == 137)
Day5.0 <- which(bare.pred$estimates$par.index == 148)
Day30.0 <- which(bare.pred$estimates$par.index == 173)
Day65.0 <- which(bare.pred$estimates$par.index == 208)
Day5.68 <- which(bare.pred$estimates$par.index == 361)
Day30.68 <- which(bare.pred$estimates$par.index == 386)
Day65.68 <- which(bare.pred$estimates$par.index == 421)

bare.pred$estimates$Group <- NA
bare.pred$estimates$Group[Day5.39] <- "Early Moderate"
bare.pred$estimates$Group[Day30.39] <- "Mid Moderate"
bare.pred$estimates$Group[Day65.39] <- "Late Moderate"
bare.pred$estimates$Group[Day5.49] <- "Early Full"
bare.pred$estimates$Group[Day30.49] <- "Mid Full"
bare.pred$estimates$Group[Day65.49] <- "Late Full"
bare.pred$estimates$Group[Day5.0] <- "Early Rest"
bare.pred$estimates$Group[Day30.0] <- "Mid Rest"
bare.pred$estimates$Group[Day65.0] <- "Late Rest"
bare.pred$estimates$Group[Day5.68] <- "Early Heavy"
bare.pred$estimates$Group[Day30.68] <- "Mid Heavy"
bare.pred$estimates$Group[Day65.68] <- "Late Heavy"
head(bare.pred$estimates)

bare.pred$estimates$Date <- NA
bare.pred$estimates$Date[Day5.39] <- "Early"
bare.pred$estimates$Date[Day30.39] <- "Mid"
bare.pred$estimates$Date[Day65.39] <- "Late"
bare.pred$estimates$Date[Day5.49] <- "Early"
bare.pred$estimates$Date[Day30.49] <- "Mid"
bare.pred$estimates$Date[Day65.49] <- "Late"
bare.pred$estimates$Date[Day5.0] <- "Early"
bare.pred$estimates$Date[Day30.0] <- "Mid"
bare.pred$estimates$Date[Day65.0] <- "Late"
bare.pred$estimates$Date[Day5.68] <- "Early"
bare.pred$estimates$Date[Day30.68] <- "Mid"
bare.pred$estimates$Date[Day65.68] <- "Late"
head(bare.pred$estimates)

bare.pred$estimates$Treat <- NA
bare.pred$estimates$Treat[Day5.39] <- "Moderate"
bare.pred$estimates$Treat[Day30.39] <- "Moderate"
bare.pred$estimates$Treat[Day65.39] <- "Moderate"
bare.pred$estimates$Treat[Day5.49] <- "Full"
bare.pred$estimates$Treat[Day30.49] <- "Full"
bare.pred$estimates$Treat[Day65.49] <- "Full"
bare.pred$estimates$Treat[Day5.0] <- "Rest"
bare.pred$estimates$Treat[Day30.0] <- "Rest"
bare.pred$estimates$Treat[Day65.0] <- "Rest"
bare.pred$estimates$Treat[Day5.68] <- "Heavy"
bare.pred$estimates$Treat[Day30.68] <- "Heavy"
bare.pred$estimates$Treat[Day65.68] <- "Heavy"
head(bare.pred$estimates)

RWBLbare.plot <- ggplot(transform(bare.pred$estimates,
                                  Treat = factor(Treat, levels=c("Rest", "Moderate", "Full", "Heavy")),
                                  Date = factor(Date, levels=c("Early", "Mid", "Late"))), 
                        aes(x = covdata, 
                            y = estimate,
                            groups = Group,
                            fill = Treat)) +
  geom_line(size = 1.5,
            aes(linetype = Date,
                colour = Treat)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.05)  +
  facet_grid(Date~.) +
  scale_colour_manual(values = c('#A2A4A2',
                                          'lightgoldenrod2',
                                          '#D4A634',
                                          '#717F5B')) +
                                            scale_fill_manual(values = c('#A2A4A2',
                                                                                  'lightgoldenrod2',
                                                                                  '#D4A634',
                                                                                  '#717F5B')) +
                                                                                    xlab("Bare Cover") +
  ylab("Estimated DSR") +
  ylim(0, 1) + 
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
        legend.position = c(.85, .1),
        legend.box = "horizontal") +                                            # remove the legend
  labs(title = "Red-winged Blackbird")

julian.pred <- covariate.predictions(plotdata,
                                     data=data.frame(Bare=mean(bare.values)),
                                     indices = c(1:71, 72:142, 143:213, 356:426))

Moderate <- which(julian.pred$estimates$par.index == c(1:71))
Full <- which(julian.pred$estimates$par.index == c(72:142))
Rest <- which(julian.pred$estimates$par.index == c(143:213))
Heavy <- which(julian.pred$estimates$par.index == c(356:426))

julian.pred$estimates$Treat[Moderate] <- "Moderate"
julian.pred$estimates$Treat[Full] <- "Full"
julian.pred$estimates$Treat[Rest] <- "Rest"
julian.pred$estimates$Treat[Heavy] <- "Heavy"
head(julian.pred$estimates)

julian.pred$estimates$Julian[Moderate] <- c(1:71)
julian.pred$estimates$Julian[Full] <- c(1:71)
julian.pred$estimates$Julian[Rest] <- c(1:71)
julian.pred$estimates$Julian[Heavy] <- c(1:71)
head(julian.pred$estimates)

RWBLjulian.plot <- ggplot(transform(julian.pred$estimates,
                                    Treat = factor(Treat, levels=c("Rest", "Moderate", "Full", "Heavy"))), 
                          aes(x = Julian,
                              y = estimate,
                              fill = Treat)) +
  geom_line(size = 1.5,
            aes(colour=Treat)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.3) +
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
        legend.position = c(.85, .1),
        legend.box = "horizontal") +
  labs(title = "Red-winged Blackbird")

RWBL.plot
RWBLbare.plot
RWBLjulian.plot

ggsave(RWBL.plot,
       filename = "outputs/figs/RWBLbeta.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(RWBLbare.plot,
       filename = "outputs/figs/RWBLbare.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(RWBLjulian.plot,
       filename = "outputs/figs/RWBLjulian.png",
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
