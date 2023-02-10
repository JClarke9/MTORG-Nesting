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


nest$Year <- as.factor(nest$Year)
nest$cTreat <- as.factor(nest$cTreat)

WEME.surv <- filter(nest, 
                    Spec=="WEME")

BRBL.surv <- filter(nest, 
                    Spec=="BRBL")

RWBL.surv <- filter(nest, 
                    Spec=="RWBL")

CCSP.surv <- filter(nest, 
                    Spec=="CCSP")

BWTE.surv <- filter(nest,
                    Spec=="BWTE")

GADW.surv <- filter(nest, 
                    Spec=="GADW")

NOPI.surv <- filter(nest, 
                    Spec=="NOPI")

MODO.surv <- filter(nest, 
                    Spec=="MODO")


WEME.surv$AgeDay1 <- WEME.surv$AgeFound - WEME.surv$FirstFound + 1
BRBL.surv$AgeDay1 <- BRBL.surv$AgeFound - BRBL.surv$FirstFound + 1
RWBL.surv$AgeDay1 <- RWBL.surv$AgeFound - RWBL.surv$FirstFound + 1
CCSP.surv$AgeDay1 <- CCSP.surv$AgeFound - CCSP.surv$FirstFound + 1
BWTE.surv$AgeDay1 <- BWTE.surv$AgeFound - BWTE.surv$FirstFound + 1
GADW.surv$AgeDay1 <- GADW.surv$AgeFound - GADW.surv$FirstFound + 1
NOPI.surv$AgeDay1 <- NOPI.surv$AgeFound - NOPI.surv$FirstFound + 1
MODO.surv$AgeDay1 <- MODO.surv$AgeFound - MODO.surv$FirstFound + 1

#remove everything but the final dataframe
rm(list = ls()[!ls() %in%  c("WEME.surv", 
                             "BRBL.surv", 
                             "RWBL.surv", 
                             "CCSP.surv", 
                             "BWTE.surv", 
                             "GADW.surv", 
                             "NOPI.surv", 
                             "MODO.surv")])


# Daily survival rate models ----------------------------------------------


WEME.mod <- mark(WEME.surv, 
                 nocc=max(WEME.surv$LastChecked), 
                 model = "Nest", 
                 model.parameters = list(S = list(formula =  ~1 + grazed + KBG + Veg.Height)))


BRBL.mod <- mark(BRBL.surv, nocc=max(BRBL.surv$LastChecked), 
                 model = "Nest", 
                 groups = "Year", 
                 model.parameters = list(S = list(formula =  ~1 + Year + Time + VOR)))


RWBL.mod <- mark(RWBL.surv, 
                 nocc=max(RWBL.surv$LastChecked), 
                 model = "Nest", 
                 groups ="cTreat", 
                 model.parameters = list(S = list(formula =  ~1 + cTreat + Time + Bare)))


CCSP.mod <- mark(CCSP.surv, 
                 nocc=max(CCSP.surv$LastChecked), 
                 model = "Nest", 
                 groups = c("Year"), 
                 model.parameters = list(S = list(formula =  ~1 + Year + Time + BHCOpres + Veg.Height)))


BWTE.mod <- mark(BWTE.surv, 
                 nocc=max(BWTE.surv$LastChecked), 
                 model = "Nest", 
                 groups = "cTreat", 
                 model.parameters = list(S = list(formula =  ~1 + cTreat + Time + LitterD)))


GADW.mod <- mark(GADW.surv, 
                 nocc=max(GADW.surv$LastChecked), 
                 model = "Nest", 
                 model.parameters = list(S = list(formula =  ~1 + LitterD)))


NOPI.mod <- mark(NOPI.surv, 
                 nocc=max(NOPI.surv$LastChecked), 
                 model = "Nest", 
                 model.parameters = list(S = list(formula =  ~1 + Litter)))


MODO.mod <- mark(MODO.surv, 
                 nocc=max(MODO.surv$LastChecked), 
                 model = "Nest", 
                 model.parameters = list(S = list(formula = ~1 + Time + I(Time^2) + KBG + LitterD)))

WEME.mod$results$beta
BRBL.mod$results$beta
RWBL.mod$results$beta
CCSP.mod$results$beta
BWTE.mod$results$beta
GADW.mod$results$beta
NOPI.mod$results$beta
MODO.mod$results$beta


# Plotting Vegetation Height Pred -----------------------------------------

WEME.pr <- process.data(WEME.surv,
                        nocc=max(WEME.surv$LastChecked),
                        model="Nest")

WEME.ddl <- make.design.data(WEME.pr)
WEME.ddl <- as.data.frame(WEME.ddl)

WEME.minH <- min(WEME.surv$Veg.Height)
WEME.maxH <- max(WEME.surv$Veg.Height)
WEME.height = seq(from = WEME.minH, 
                  to = WEME.maxH, 
                  length = 100)

WEME.predH <- covariate.predictions(WEME.mod,
                                    data=data.frame(Veg.Height = WEME.height),
                                    indices = 1)

WEME.predH$estimates$Group <- "WEME"


CCSP.pr <- process.data(CCSP.surv,
                        nocc=max(CCSP.surv$LastChecked),
                        group = "Year",
                        model="Nest")

CCSP.ddl <- make.design.data(CCSP.pr)
CCSP.ddl <- as.data.frame(CCSP.ddl)

CCSP.minH <- min(CCSP.surv$Veg.Height)
CCSP.maxH <- max(CCSP.surv$Veg.Height)
CCSP.height = seq(from = CCSP.minH, 
                  to = CCSP.maxH, 
                  length = 100)

CCSP.predH <- covariate.predictions(CCSP.mod,
                                    data=data.frame(Veg.Height = CCSP.height),
                                    indices = 29)

CCSP.predH$estimates$Group <- "CCSP"

predH <- list(WEME = as.data.frame(WEME.predH),
              CCSP = as.data.frame(CCSP.predH))

predH.plot <- predH |> 
  map_df(~bind_rows(.x, .id = "source"))


height.plot <- ggplot(predH.plot,
                      aes(x = estimates.covdata,
                          y = estimates.estimate,
                          groups = estimates.Group,
                          fill = estimates.Group)) +
  geom_line(size = 1.5,
            aes(color = estimates.Group)) +
  geom_ribbon(aes(ymin = estimates.lcl, 
                  ymax = estimates.ucl), 
              alpha = 0.20,
              show.legend = FALSE)  +
  scale_color_manual(values = c("#B0A8A5", "#F6D056")) +
  scale_fill_manual(values = c("#B0A8A5", "#F6D056")) +
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
        text=element_text(size=16,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.background = element_rect(fill=NA),
        legend.position = c(.9, .1),
        legend.box = "horizontal") +
  labs(title = "Vegetation Height",
       x = "Vegetation Height (cm)",
       y = "Estimated DSR",
       color = "Species",
       fill = NULL) +
  xlim(90, 1000) +
  ylim(0.6, 1)

height.plot


BRBL.pr <- process.data(BRBL.surv,
                        nocc=max(BRBL.surv$LastChecked),
                        group = "Year",
                        model="Nest")

BRBL.ddl <- make.design.data(BRBL.pr)
BRBL.ddl <- as.data.frame(BRBL.ddl)

BRBL.minV <- min(BRBL.surv$VOR)
BRBL.maxV <- max(BRBL.surv$VOR)
BRBL.vor = seq(from = BRBL.minV, 
               to = BRBL.maxV, 
               length = 100)

BRBL.predV <- covariate.predictions(BRBL.mod,
                                    data=data.frame(VOR = BRBL.vor),
                                    indices = 29)

BRBL.predV$estimates$Group <- "BRBL"

VOR.plot <- ggplot(BRBL.predV$estimates,
                      aes(x = covdata,
                          y = estimate,
                          groups = Group,
                          fill = Group)) +
  geom_line(size = 1.5,
            aes(color = Group)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.20,
              show.legend = FALSE)  +
  scale_color_manual(values = "#060202") +
  scale_fill_manual(values = "#060202") +
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
        text=element_text(size=16,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.background = element_rect(fill=NA),
        legend.position = c(.9, .1),
        legend.box = "horizontal") +
  labs(title = "Visual Obstruction Reading",
       x = "VOR",
       y = "Estimated DSR",
       color = "Species") +
  ylim(0.25, 1)

VOR.plot


BWTE.pr <- process.data(BWTE.surv,
                        nocc=max(BWTE.surv$LastChecked),
                        model="Nest")

BWTE.ddl <- make.design.data(BWTE.pr)
BWTE.ddl <- as.data.frame(BWTE.ddl)

BWTE.minL <- min(BWTE.surv$LitterD)
BWTE.maxL <- max(BWTE.surv$LitterD)
BWTE.litdep = seq(from = BWTE.minL, 
                  to = BWTE.maxL, 
                  length = 100)

BWTE.predLD <- covariate.predictions(BWTE.mod,
                                    data=data.frame(LitterD = BWTE.litdep),
                                    indices = 1)

BWTE.predLD$estimates$Group <- "BWTE"

GADW.pr <- process.data(GADW.surv,
                        nocc=max(GADW.surv$LastChecked),
                        model="Nest")

GADW.ddl <- make.design.data(GADW.pr)
GADW.ddl <- as.data.frame(GADW.ddl)

GADW.minL <- min(GADW.surv$LitterD)
GADW.maxL <- max(GADW.surv$LitterD)
GADW.litdep = seq(from = GADW.minL, 
                  to = GADW.maxL, 
                  length = 100)

GADW.predLD <- covariate.predictions(GADW.mod,
                                    data=data.frame(LitterD = GADW.litdep),
                                    indices = 1)

GADW.predLD$estimates$Group <- "GADW"

MODO.pr <- process.data(MODO.surv,
                        nocc=max(MODO.surv$LastChecked),
                        model="Nest")

MODO.ddl <- make.design.data(MODO.pr)
MODO.ddl <- as.data.frame(MODO.ddl)

MODO.minL <- min(MODO.surv$LitterD)
MODO.maxL <- max(MODO.surv$LitterD)
MODO.litdep = seq(from = MODO.minL, 
                  to = MODO.maxL, 
                  length = 100)

MODO.predLD <- covariate.predictions(MODO.mod,
                                    data=data.frame(LitterD = MODO.litdep),
                                    indices = 1)

MODO.predLD$estimates$Group <- "MODO"

predLD <- list(BWTE = as.data.frame(BWTE.predLD),
              GADW = as.data.frame(GADW.predLD),
              MODO = as.data.frame(MODO.predLD))

predLD.plot <- predLD |> 
  map_df(~bind_rows(.x, .id = "source"))


litdep.plot <- ggplot(predLD.plot,
                      aes(x = estimates.covdata,
                          y = estimates.estimate,
                          groups = estimates.Group,
                          fill = estimates.Group)) +
  geom_line(size = 1.5,
            aes(color = estimates.Group)) +
  geom_ribbon(aes(ymin = estimates.lcl, 
                  ymax = estimates.ucl), 
              alpha = 0.20,
              show.legend = FALSE)  +
  scale_color_manual(values = c("#C9D5E0","#999395", "#DAB27F")) +
  scale_fill_manual(values = c("#C9D5E0","#999395", "#DAB27F")) +
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
        text=element_text(size=16,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.background = element_rect(fill=NA),
        legend.position = c(.9, .1),
        legend.box = "horizontal") +
  labs(title = "Litter Depth",
       x = "Litter Depth (mm)",
       y = "Estimated DSR",
       color = "Species") +
  xlim(0, 100) +
  ylim(.85, 1)

litdep.plot


NOPI.pr <- process.data(NOPI.surv,
                        nocc=max(NOPI.surv$LastChecked),
                        group = "Year",
                        model="Nest")

NOPI.ddl <- make.design.data(NOPI.pr)
NOPI.ddl <- as.data.frame(NOPI.ddl)

NOPI.minL <- min(NOPI.surv$Litter)
NOPI.maxL <- max(NOPI.surv$Litter)
NOPI.lit = seq(from = NOPI.minL, 
               to = NOPI.maxL, 
               length = 100)

NOPI.predL <- covariate.predictions(NOPI.mod,
                                    data=data.frame(Litter = NOPI.lit),
                                    indices = 1)

NOPI.predL$estimates$Group <- "NOPI"

lit.plot <- ggplot(NOPI.predL$estimates,
                   aes(x = covdata,
                       y = estimate,
                       groups = Group,
                       fill = Group)) +
  geom_line(size = 1.5,
            aes(color = Group)) +
  geom_ribbon(aes(ymin = lcl, 
                  ymax = ucl), 
              alpha = 0.20,
              show.legend = FALSE)  +
  scale_color_manual(values = "#5E3028")  +
  scale_fill_manual(values = "#5E3028") +
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
        text=element_text(size=16,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.background = element_rect(fill=NA),
        legend.position = c(.9, .1),
        legend.box = "horizontal") +
  labs(title = "Litter Cover",
       x = "Litter Cover (%)",
       y = "Estimated DSR",
       color = "Species") +
  ylim(.4, 1)

lit.plot

DSR_pred <- list(CCSP.predH,
                 WEME.predH,
                 BRBL.predV,
                 GADW.predLD,
                 BWTE.predLD,
                 MODO.predLD,
                 NOPI.predL)

DSR_range <- data.frame()

for(i in 1:length(DSR_pred)) {
  z <- data.frame(spec = unique(DSR_pred[[i]]$estimate$Group),
                  max = max(DSR_pred[[i]]$estimates$estimate))
  y <- data.frame(spec = unique(DSR_pred[[i]]$estimate$Group),
                  min = min(DSR_pred[[i]]$estimates$estimate))
  
  x <- full_join(z, y)
  
  DSR_range <- bind_rows(DSR_range, x)
}

DSR_range <- data.frame(max = max(CCSP.predH$estimates$estimate,
                                  WEME.predH$estimates$estimate,
                                  BRBL.predV$estimates$estimate,
                                  GADW.predLD$estimates$estimate,
                                  BWTE.predLD$estimates$estimate,
                                  MODO.predLD$estimates$estimate,
                                  NOPI.predL$estimates$estimate),
                        min = min(CCSP.predH$estimates$estimate,
                                  WEME.predH$estimates$estimate,
                                  BRBL.predV$estimates$estimate,
                                  GADW.predLD$estimates$estimate,
                                  BWTE.predLD$estimates$estimate,
                                  MODO.predLD$estimates$estimate,
                                  NOPI.predL$estimates$estimate))

ggsave(height.plot,
       filename = "outputs/figs/heightDSR_pred.png",
       dpi = "retina",
       bg = "white",
       height = 6.5,
       width = 6.5)

ggsave(VOR.plot,
       filename = "outputs/figs/vorDSR_pred.png",
       dpi = "retina",
       bg = "white",
       height = 6.5,
       width = 6.5)

ggsave(litdep.plot,
       filename = "outputs/figs/litdepDSR_pred.png",
       dpi = "retina",
       bg = "white",
       height = 6.5,
       width = 6.5)

ggsave(lit.plot,
       filename = "outputs/figs/litDSR_pred.png",
       dpi = "retina",
       bg = "white",
       height = 6.5,
       width = 6.5)

# If you want to clean up the mark*.inp, .vcv, .res and .out
# and .tmp files created by RMark in the working directory,
# execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all=TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
# files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)
