# Loading libraries -------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(RMark)
library(cowplot)


# Data import -------------------------------------------------------------


nest <- read.csv("working/RMarknesting.csv", 
                 row.names=1)


windowsFonts(my_font = windowsFont("Times New Roman"))


# Subsetting data ---------------------------------------------------------


WEME.surv <- filter(nest, 
                    Spec=="WEME")                                         # select out only WEME nests

WEME.surv$AgeDay1 <- WEME.surv$AgeFound - WEME.surv$FirstFound + 1
WEME.surv$Year <- as.factor(WEME.surv$Year)
WEME.surv$cTreat <- as.factor(WEME.surv$cTreat)


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

x <- create.stage.var(WEME.surv, 
                      "AgeDay1", 
                      "Incub", 
                      rep(1,max(WEME.surv$LastChecked)), 
                      14)

WEME.surv <- bind_cols(WEME.surv, x)

rm(list = ls()[!ls() %in% "WEME.surv"])


# Daily survival rate models ----------------------------------------------


WEME.pr <- process.data(WEME.surv,
                        nocc=max(WEME.surv$LastChecked),
                        groups = c("cTreat",
                                   "Year"),
                        model="Nest")

# candidate models for the null vs. treatment model
WEME1.run <- function()
{
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  # 2. DSR varies with current treatment
  S.cTreat = list(formula = ~1 + cTreat)
  
  # 3. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 4. DSR varies With year interacting with treatment (treatment identity varies by year)
  S.treatyear = list(formula = ~1 + cTreat:Year)
  
  WEME.model.list = create.model.list("Nest")
  WEME1.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
WEME1.results <- WEME1.run()
WEME1.results

WEME1.results$S.Dot$results$beta
WEME1.results$S.year$results$beta

#candidate model set for time trends
WEME2.run <- function()
{
  # 1. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies across stages
  S.stage = list(formula = ~1 + Incub)
  
  # 1. a model for constant daily survival rate
  S.Dot = list(formula =  ~1)
  
  WEME.model.list = create.model.list("Nest")
  WEME2.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

#results of candidate model set 
#all results
WEME2.results <- WEME2.run()
WEME2.results

WEME2.results$S.Dot$results$beta
WEME2.results$S.time$results$beta
WEME2.results$S.stage$results$beta

#candidate model set for parasitism and grazing
WEME3.run <- function()
{
  # 1. DSR varies with BHCO number
  S.bhco = list(formula = ~1 + BHCONum)
  
  S.grazed = list(formula = ~1 + grazed)
  
  S.bhcop = list(formula = ~1 + BHCOpres)
  
  S.grazep = list(formula = ~1 + grazep)
  
  S.Dot = list(formula = ~1)
  
  WEME.model.list = create.model.list("Nest")
  WEME3.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

#results of candidate model set 
#all results
WEME3.results <- WEME3.run()
WEME3.results

WEME3.results$S.grazed$results$beta
WEME3.results$S.Dot$results$beta
WEME3.results$S.bhcop$results$beta
WEME3.results$S.grazep$results$beta

#candidate model set for parasitism and nest stage
WEME4.run <- function()
{
  # 1. DSR varies with KBG
  S.kbg = list(formula =  ~1 + grazed + KBG)
  
  # 2. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.smooth.brome = list(formula = ~1 + grazed + SmoothB)
  
  # 3. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + grazed + Litter)
  
  # 4. DSR varies with Bare
  S.bare = list(formula =  ~1 + grazed + Bare)
  
  # 5. DSR varies with Forb
  S.forb = list(formula =  ~1 + grazed + Forb)
  
  # 6. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + grazed + Grasslike)
  
  # 7. DSR varies with Woody
  S.woody = list(formula =  ~1 + grazed + Woody)
  
  # 8. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + grazed + LitterD)
  
  # 9. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + grazed + Veg.Height)
  
  # 1 + grazed0. DSR varies with VOR
  S.vor = list(formula =  ~1 + grazed + VOR)
  
  # 1 + grazed1 + grazed. DSR varies with litter and Veg Height
  S.litheight = list(formula = ~1 + grazed + Litter + Veg.Height)
  
  # 1 + grazed2. DSR varies with woody and litter depth
  S.woodylitdep = list(formula = ~1 + grazed + Woody + LitterD)
  
  # 1 + grazed3. DSR varies with KBG and litter depth
  S.kbglitdep = list(formula = ~1 + grazed + KBG + LitterD)
  
  # 1 + grazed4. DSR varies with KBG and Veg.Height
  S.kbgheight = list(formula = ~1 + grazed + KBG + Veg.Height)
  
  # 1 + grazed5. DSR varies with woody and Veg Height
  S.woodyheight = list(formula = ~1 + grazed + Woody + Veg.Height)
  
  # 1. a model for constant daily survival rate
  S.grazed = list(formula =  ~1 + grazed)
  
  WEME.model.list = create.model.list("Nest")
  WEME4.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = FALSE)
}

#results of candidate model set 
#all results
WEME4.results <- WEME4.run()
WEME4.results

WEME4.results$S.kbgheight$results$beta
WEME4.results$S.bare$results$beta
WEME4.results$S.height$results$beta
WEME4.results$S.forb$results$beta
WEME4.results$S.grazed$results$beta

WEME4.results$S.bare$results$real


# Plotting beta coefficients ----------------------------------------------


WEME.mod <- mark(WEME.surv, 
                 nocc=max(WEME.surv$LastChecked), 
                 model = "Nest", 
                 model.parameters = list(S = list(formula =  ~1 + grazed + KBG + Veg.Height)))

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
