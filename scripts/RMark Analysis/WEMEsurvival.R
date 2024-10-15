# Loading libraries -------------------------------------------------------


library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)
library(MuMIn)
source('scripts/Functions/RMark_Stage_Code.R')

windowsFonts(my_font = windowsFont('Gandhi Sans'))


# Data import -------------------------------------------------------------


nest <- read.csv('working/RMarknesting.csv')


# Subsetting data ---------------------------------------------------------


WEME.surv <- filter(nest, 
                    Spec == 'WEME' & Stage != 'Laying')                                         # select out only WEME nest

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
                         levels = c('2021', '2022', '2023', '2024'))

WEME.surv$Nestling <- factor(WEME.surv$Nestling,
                             level = c('0', '1'))

str(WEME.surv)


# Creating stage variable -------------------------------------------------


x <- create.stage.var(WEME.surv, 
                      'AgeDay1', 
                      'Incub', 
                      rep(1,max(WEME.surv$LastChecked)), 
                      12)

WEME.surv <- bind_cols(WEME.surv, x)

rm(list = ls()[!ls() %in% c('WEME.surv')])


# Daily survival rate models ----------------------------------------------


WEME.pr <- process.data(WEME.surv,
                        nocc = max(WEME.surv$LastChecked),
                        groups = c('Year',
                                   'Nestling'),
                        model = 'Nest')



# Grazing candidate model set
WEME1.run <- function()
{
  # 5. DSR varies with nest age
  S.null = list(formula = ~1)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + grazed)
  
  # 4. DSR varies with the previous years grazing intensity
  S.pDoD = list(formula = ~1 + pDoD)
  
  WEME.model.list = create.model.list('Nest')
  WEME1.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME1.results <- WEME1.run()
WEME1.results

coef(WEME1.results$S.null)
confint(WEME1.results$S.null, level = 0.85)



# Temporal candidate model set
WEME2.run <- function()
{
  # 1. DSR varies with time
  S.null = list(formula = ~1)
  
  # 2. DSR varies with time
  S.time = list(formula = ~1 + Time)
  
  # 3. DSR varies with quadratic effect of date
  S.quad = list(formula = ~1 + Time + I(Time^2))
  
  # 4. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  WEME.model.list = create.model.list('Nest')
  WEME2.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME2.results <- WEME2.run()
WEME2.results

coef(WEME2.results$S.null)



# Biological candidate model set
WEME3.run <- function()
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
  S.stage = list(formula = ~1 + Nestling)
  
  WEME.model.list = create.model.list('Nest')
  WEME3.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME3.results <- WEME3.run()
WEME3.results

coef(WEME3.results$S.stage)
confint(WEME3.results$S.stage, level = 0.85)



# Vegetation candidate model set
WEME4.run <- function()
{
  # 5. DSR varies with nest age
  S.stage = list(formula = ~1 + Nestling)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula =  ~1 + Nestling + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Nestling + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula =  ~1 + Nestling + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula =  ~1 + Nestling + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula =  ~1 + Nestling + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula =  ~1 + Nestling + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula =  ~1 + Nestling + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula =  ~1 + Nestling + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula =  ~1 + Nestling + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula =  ~1 + Nestling + VOR)
  
  WEME.model.list = create.model.list('Nest')
  WEME4.results = mark.wrapper(WEME.model.list,
                               data = WEME.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
WEME4.results <- WEME4.run()
WEME4.results

coef(WEME4.results$S.stage)
confint(WEME4.results$S.stage, level = 0.85)

WEME4.results$S.stage$results$real |> 
  summarize(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))

(WEME.real <- as.data.frame(WEME4.results$S.stage$results$real) |> 
    rownames_to_column(var = 'Group') |> 
    mutate(Stage = case_when(
      grepl('20210', Group) ~ 'Incubating',
      grepl('20220', Group) ~ 'Incubating',
      grepl('20230', Group) ~ 'Incubating',
      grepl('20240', Group) ~ 'Incubating',
      grepl('20211', Group) ~ 'Nestling',
      grepl('20221', Group) ~ 'Nestling',
      grepl('20231', Group) ~ 'Nestling',
      grepl('20241', Group) ~ 'Nestling')) |> 
    group_by(Stage) |> 
    summarize(estimate = mean(estimate), 
              se = mean(se), 
              lcl = mean(lcl), 
              ucl = mean(ucl)))


# Plotting beta coefficients ----------------------------------------------


WEME.beta <- coef(WEME4.results$S.stage) |>
  cbind(confint(WEME4.results$S.stage, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = 'Variable') |> 
  rename(c('Coefficient' = 'estimate',
           'lcl' = '7.5 %',
           'ucl' = '92.5 %'))

WEME.beta$Variable <- gsub('S:', '', WEME.beta$Variable)

str(WEME.beta)

(WEME.plot <- ggplot(WEME.beta[2,], 
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
    theme(plot.title = element_text(family = 'my_font',
                                    hjust = .5,
                                    size = 20,
                                    vjust = 1,
                                    colour = 'black'),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = NA,                     # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = NA,                      # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          axis.line = element_line(colour = 'black'),                             # color the x and y axis
          axis.text = element_text(size = 12, 
                                   colour = 'black'),                    # color the axis text
          axis.ticks = element_line(colour = 'black'),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = 'black')) +                                    # change the color of the axis titles
    labs(title = 'Western Meadowlark',
         x = NULL,
         y = expression('Beta ' (beta))))


# Creating predictive plots -----------------------------------------------


WEME.ddl <- make.design.data(WEME.pr) |> 
  as.data.frame()

plotdata <- WEME4.results$S.stage


filter(WEME.ddl, 
       S.Nestling == 0 & S.time == 1 & S.Year == 2021| 
         S.Nestling == 1 & S.time == 1 & S.Year == 2021)

stage.pred <- covariate.predictions(plotdata,
                                    indices = c(1, 305))

inc <- which(stage.pred$estimates$par.index == 1)
nst <- which(stage.pred$estimates$par.index == 305)

stage.pred$estimates$Stage <- NA
stage.pred$estimates$Stage[inc] <- 'Incubating'
stage.pred$estimates$Stage[nst] <- 'Nestling'

stage.pred$estimates$Day <- 1

(WEMEstage.plot <- ggplot(transform(stage.pred$estimates,
                                    Day = factor(Day, levels = '1')), 
                          aes(x = Stage, 
                              y = estimate,
                              groups = Day,
                              fill = Day)) +
    geom_point(size = 4,
               aes(color = Day)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c('#D4A634')) +
    scale_fill_manual(values = c('#D4A634')) +
    theme(plot.title = element_text(family = 'my_font',                             # select the font for the title
                                    size = 16,
                                    hjust = .5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = NA,                                # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = NA,                                 # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          axis.line = element_line(colour = 'black'),                             # color the x and y axis
          axis.text.y = element_text(size = 12, colour = 'black'),                    # color the axis text
          axis.text.x = element_text(size = 12, colour = 'black'),
          axis.ticks = element_line(colour = 'black'),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = 'black'),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = 'none') +
    labs(title = 'Western Meadowlark',
         x = 'Stage',
         y = 'Daily Survival Rate'))


ggsave(WEME.plot,
       filename = 'outputs/figs/betaWEME.png',
       dpi = 'print',
       bg = 'white',
       height = 6,
       width = 6)

ggsave(WEMEstage.plot,
       filename = 'outputs/figs/stageWEME.png',
       dpi = 'print',
       bg = 'white',
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
