# Loading libraries -------------------------------------------------------


library(ggplot2)
library(vegan)
library(tidyverse)
library(RMark)
library(MuMIn)
library(cowplot)
source("scripts/Functions/RMark_Stage_Code.R")

windowsFonts(my_font = windowsFont("Gandhi Sans"))


# Data import -------------------------------------------------------------


nest <- read.csv("working/RMarknesting.csv")


# Subsetting data ---------------------------------------------------------


CCSP.surv <- filter(nest, 
                    Spec == "CCSP" & Stage != "Laying")

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
                         levels = c("2021", "2022", "2023", "2024"))

CCSP.surv$Nestling <- factor(CCSP.surv$Nestling,
                             level = c("0", "1"))

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
                        nocc = max(CCSP.surv$LastChecked),
                        groups = c("Year",
                                   "Nestling"),
                        model = "Nest")



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



# Nest stage/age candidate model set
CCSP2.run <- function()
{
  # 1. DSR varies with year
  S.year = list(formula = ~1 + Year)
  
  # 2. DSR varies with nest stage
  S.stage = list(formula = ~1 + Year + Nestling)
  
  # 3. DSR varies with nest age
  S.age = list(formula = ~1 + Year + NestAge)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP2.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP2.results <- CCSP2.run()
CCSP2.results

coef(CCSP2.results$S.stage)
confint(CCSP2.results$S.stage, level = 0.85)



# Biological candidate model set
CCSP3.run <- function()
{
  # 1. DSR varies with nest stage and age
  S.stage = list(formula = ~1 + Year + Nestling)
  
  # 2. DSR varies with BHCO number
  S.bhcon = list(formula = ~1 + Year + Nestling + BHCONum)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP3.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP3.results <- CCSP3.run()
CCSP3.results

coef(CCSP3.results$S.stage)
confint(CCSP3.results$S.stage, level = 0.85)



# Grazing candidate model set
CCSP4.run <- function()
{
  # 1. DSR varies with nest stage and age
  S.stage = list(formula = ~1 + Year + Nestling)
  
  # 2. DSR varies with the number of days a nest experienced grazing
  S.grazed = list(formula = ~1 + Year + Nestling + grazed)
  
  # 3. DSR varies with the previous Year + Nestlings grazing intensity
  S.pDoD = list(formula = ~1 + Year + Nestling + pDoD)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP4.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP4.results <- CCSP4.run()
CCSP4.results

coef(CCSP4.results$S.stage)
confint(CCSP4.results$S.stage, level = 0.85)



# Vegetation candidate model set
CCSP5.run <- function()
{
  # 1. DSR varies with nest stage and age
  S.stage = list(formula = ~1 + Year + Nestling)
  
  # 2. DSR varies with KBG
  S.kbg = list(formula = ~1 + Year + Nestling + KBG)
  
  # 3. DSR varies with Smooth Brome (correlated with KBG and Litter Depth)
  S.brome = list(formula = ~1 + Year + Nestling + SmoothB)
  
  # 4. DSR varies with Litter (correlated with KBG)
  S.lit = list(formula = ~1 + Year + Nestling + Litter)
  
  # 5. DSR varies with Bare
  S.bare = list(formula = ~1 + Year + Nestling + Bare)
  
  # 6. DSR varies with Forb
  S.forb = list(formula = ~1 + Year + Nestling + Forb)
  
  # 7. DSR varies with Grasslike  (correlated with KBG)
  S.grass = list(formula = ~1 + Year + Nestling + Grasslike)
  
  # 8. DSR varies with Woody
  S.woody = list(formula = ~1 + Year + Nestling + Woody)
  
  # 9. DSR varies with Litter Depth (correlated with VOR)
  S.litdep = list(formula = ~1 + Year + Nestling + LitterD)
  
  # 10. DSR varies with Veg Height (correlated with VOR)
  S.height = list(formula = ~1 + Year + Nestling + Veg.Height)
  
  # 11. DSR varies with VOR
  S.vor = list(formula = ~1 + Year + Nestling + VOR)
  
  CCSP.model.list = create.model.list("Nest")
  CCSP5.results = mark.wrapper(CCSP.model.list,
                               data = CCSP.pr,
                               adjust = FALSE,
                               delete = TRUE)
}

# Results of candidate model set
CCSP5.results <- CCSP5.run()
CCSP5.results

coef(CCSP5.results$S.bare)
confint(CCSP5.results$S.bare, level = 0.85)

coef(CCSP5.results$S.litdep)
confint(CCSP5.results$S.litdep, level = 0.85)


CCSP6.results <- mark(CCSP.surv, 
                      nocc = max(CCSP.surv$LastChecked), 
                      model = "Nest", 
                      groups = c("Year",
                                 "Nestling"),
                      adjust = FALSE,
                      delete = TRUE, 
                      model.parameters = list(S = list(formula = ~1 + Year + Nestling + LitterD + Bare)))


CCSP6.results$results$real |> 
  summarize(estimate = mean(estimate),
            se = mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))

(CCSP.real <- as.data.frame(CCSP6.results$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023",
      grepl("2024", Group) ~ "2024"),
      Stage = case_when(
        grepl("20210", Group) ~ "Incubating",
        grepl("20220", Group) ~ "Incubating",
        grepl("20230", Group) ~ "Incubating",
        grepl("20240", Group) ~ "Incubating",
        grepl("20211", Group) ~ "Nestling",
        grepl("20221", Group) ~ "Nestling",
        grepl("20231", Group) ~ "Nestling",
        grepl("20241", Group) ~ "Nestling")) |> 
    select(Year, Stage, estimate, se, lcl, ucl))

(CCSP.year <- CCSP.real |> 
    group_by(Year) |> 
    summarize(mean = mean(estimate)))

(CCSP.stage <- CCSP.real |> 
    group_by(Stage) |> 
    summarize(mean = mean(estimate)))


# Plotting beta coefficients ----------------------------------------------

# I need to model average my estimates


CCSP.beta <- coef(CCSP6.results) |>
  cbind(confint(CCSP6.results, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %"))

CCSP.beta$Variable <- gsub("S:", "", CCSP.beta$Variable)

str(CCSP.beta)

(CCSP.plot <- ggplot(CCSP.beta[6:7,], 
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
    theme(plot.title = element_text(family = "my_font",
                                    hjust = .5,
                                    size = 20,
                                    vjust = 1,
                                    colour = "black"),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = NA,                     # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = NA,                      # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text = element_text(size = 12, 
                                   colour = "black"),                    # color the axis text
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black")) +                                    # change the color of the axis titles
    labs(title = "Clay-colored Sparrow",
         x = NULL,
         y = expression("Beta " (beta))))


# Creating predictive plots -----------------------------------------------


CCSP.ddl <- make.design.data(CCSP.pr) |> 
  as.data.frame()

plotdata <- CCSP6.results

barevalues <- seq(from = min(CCSP.surv$Bare),
                  to = max(CCSP.surv$Bare),
                  length = 100)

litdvalues <- seq(from = min(CCSP.surv$LitterD),
                  to = max(CCSP.surv$LitterD),
                  length = 100)



filter(CCSP.ddl, S.age == 0)

stage.pred <- covariate.predictions(plotdata,
                                    data = data.frame(Bare = mean(barevalues),
                                                      LitterD = mean(litdvalues)),
                                    indices = c(1, 89, 177, 265, 
                                                353, 441, 529, 617))

inc2021 <- which(stage.pred$estimates$par.index == 1)
inc2022 <- which(stage.pred$estimates$par.index == 89)
inc2023 <- which(stage.pred$estimates$par.index == 177)
inc2024 <- which(stage.pred$estimates$par.index == 265)
nst2021 <- which(stage.pred$estimates$par.index == 353)
nst2022 <- which(stage.pred$estimates$par.index == 441)
nst2023 <- which(stage.pred$estimates$par.index == 529)
nst2024 <- which(stage.pred$estimates$par.index == 617)

stage.pred$estimates$Year <- NA
stage.pred$estimates$Year[inc2021] <- "2021"
stage.pred$estimates$Year[inc2022] <- "2022"
stage.pred$estimates$Year[inc2023] <- "2023"
stage.pred$estimates$Year[inc2024] <- "2024"
stage.pred$estimates$Year[nst2021] <- "2021"
stage.pred$estimates$Year[nst2022] <- "2022"
stage.pred$estimates$Year[nst2023] <- "2023"
stage.pred$estimates$Year[nst2024] <- "2024"

stage.pred$estimates$Stage <- NA
stage.pred$estimates$Stage[inc2021] <- "Incubating"
stage.pred$estimates$Stage[inc2022] <- "Incubating"
stage.pred$estimates$Stage[inc2023] <- "Incubating"
stage.pred$estimates$Stage[inc2024] <- "Incubating"
stage.pred$estimates$Stage[nst2021] <- "Nestling"
stage.pred$estimates$Stage[nst2022] <- "Nestling"
stage.pred$estimates$Stage[nst2023] <- "Nestling"
stage.pred$estimates$Stage[nst2024] <- "Nestling"

(CCSPstage.plot <- ggplot(transform(stage.pred$estimates,
                                    Year = factor(Year, levels = c("2021", "2022", "2023", "2024"))), 
                          aes(x = Stage, 
                              y = estimate,
                              groups = Year,
                              fill = Year)) +
    geom_point(size = 4,
               aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c("#A2A4A2",
                                   "#717F5B",
                                   "#D4A634",
                                   "goldenrod4")) +
    scale_fill_manual(values = c("#A2A4A2",
                                 "#717F5B",
                                 "#D4A634",
                                 "goldenrod4")) +
    theme(plot.title = element_text(family = "my_font",                             # select the font for the title
                                    size = 16,
                                    hjust = .5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = NA,                                # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = NA,                                 # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_text(size = 12, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = c(.85, .1),
          legend.box = "horizontal") +
    labs(title = "Nest Survival",
         color = "Year",
         x = "Stage",
         y = "Daily Survival Rate"))


bare.pred <- covariate.predictions(plotdata,
                                   data = data.frame(Bare = barevalues,
                                                     LitterD = mean(litdvalues)),
                                   indices = c(1, 89, 177, 265, 
                                               353, 441, 529, 617))

inc2021 <- which(bare.pred$estimates$par.index == 1)
inc2022 <- which(bare.pred$estimates$par.index == 89)
inc2023 <- which(bare.pred$estimates$par.index == 177)
inc2024 <- which(bare.pred$estimates$par.index == 265)
nst2021 <- which(bare.pred$estimates$par.index == 353)
nst2022 <- which(bare.pred$estimates$par.index == 441)
nst2023 <- which(bare.pred$estimates$par.index == 529)
nst2024 <- which(bare.pred$estimates$par.index == 617)

bare.pred$estimates$Year <- NA
bare.pred$estimates$Year[inc2021] <- "2021"
bare.pred$estimates$Year[inc2022] <- "2022"
bare.pred$estimates$Year[inc2023] <- "2023"
bare.pred$estimates$Year[inc2024] <- "2024"
bare.pred$estimates$Year[nst2021] <- "2021"
bare.pred$estimates$Year[nst2022] <- "2022"
bare.pred$estimates$Year[nst2023] <- "2023"
bare.pred$estimates$Year[nst2024] <- "2024"

bare.pred$estimates$Stage <- NA
bare.pred$estimates$Stage[inc2021] <- "Incubating"
bare.pred$estimates$Stage[inc2022] <- "Incubating"
bare.pred$estimates$Stage[inc2023] <- "Incubating"
bare.pred$estimates$Stage[inc2024] <- "Incubating"
bare.pred$estimates$Stage[nst2021] <- "Nestling"
bare.pred$estimates$Stage[nst2022] <- "Nestling"
bare.pred$estimates$Stage[nst2023] <- "Nestling"
bare.pred$estimates$Stage[nst2024] <- "Nestling"

(CCSPbare.plot <- ggplot(transform(bare.pred$estimates,
                                    Year = factor(Year, levels = c("2021", "2022", "2023", "2024"))), 
                          aes(x = Bare, 
                              y = estimate,
                              groups = Year,
                              fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c("#A2A4A2",
                                   "#717F5B",
                                   "#D4A634",
                                   "goldenrod4")) +
    scale_fill_manual(values = c("#A2A4A2",
                                 "#717F5B",
                                 "#D4A634",
                                 "goldenrod4")) +
    theme(plot.title = element_text(family = "my_font",                             # select the font for the title
                                    size = 16,
                                    hjust = .5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = NA,                                # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = NA,                                 # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_text(size = 12, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = c(.85, .1),
          legend.box = "horizontal") +
    facet_grid(~Stage) + 
    labs(title = "Nest Survival",
         color = "Year",
         x = "Bare Ground (Percent Cover)",
         y = "Daily Survival Rate"))


litd.pred <- covariate.predictions(plotdata,
                                   data = data.frame(Bare = mean(barevalues),
                                                     LitterD = litdvalues),
                                   indices = c(1, 89, 177, 265, 
                                               353, 441, 529, 617))

inc2021 <- which(litd.pred$estimates$par.index == 1)
inc2022 <- which(litd.pred$estimates$par.index == 89)
inc2023 <- which(litd.pred$estimates$par.index == 177)
inc2024 <- which(litd.pred$estimates$par.index == 265)
nst2021 <- which(litd.pred$estimates$par.index == 353)
nst2022 <- which(litd.pred$estimates$par.index == 441)
nst2023 <- which(litd.pred$estimates$par.index == 529)
nst2024 <- which(litd.pred$estimates$par.index == 617)

litd.pred$estimates$Year <- NA
litd.pred$estimates$Year[inc2021] <- "2021"
litd.pred$estimates$Year[inc2022] <- "2022"
litd.pred$estimates$Year[inc2023] <- "2023"
litd.pred$estimates$Year[inc2024] <- "2024"
litd.pred$estimates$Year[nst2021] <- "2021"
litd.pred$estimates$Year[nst2022] <- "2022"
litd.pred$estimates$Year[nst2023] <- "2023"
litd.pred$estimates$Year[nst2024] <- "2024"

litd.pred$estimates$Stage <- NA
litd.pred$estimates$Stage[inc2021] <- "Incubating"
litd.pred$estimates$Stage[inc2022] <- "Incubating"
litd.pred$estimates$Stage[inc2023] <- "Incubating"
litd.pred$estimates$Stage[inc2024] <- "Incubating"
litd.pred$estimates$Stage[nst2021] <- "Nestling"
litd.pred$estimates$Stage[nst2022] <- "Nestling"
litd.pred$estimates$Stage[nst2023] <- "Nestling"
litd.pred$estimates$Stage[nst2024] <- "Nestling"

(CCSPlitd.plot <- ggplot(transform(litd.pred$estimates,
                                   Year = factor(Year, levels = c("2021", "2022", "2023", "2024"))), 
                         aes(x = LitterD, 
                             y = estimate,
                             groups = Year,
                             fill = Year)) +
    geom_line(linewidth = 1.5,
              aes(color = Year)) +
    scale_linetype_manual(values = c(1, 3, 2)) +
    scale_colour_manual(values = c("#A2A4A2",
                                   "#717F5B",
                                   "#D4A634",
                                   "goldenrod4")) +
    scale_fill_manual(values = c("#A2A4A2",
                                 "#717F5B",
                                 "#D4A634",
                                 "goldenrod4")) +
    theme(plot.title = element_text(family = "my_font",                             # select the font for the title
                                    size = 16,
                                    hjust = .5),
          panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = NA,                                # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = NA,                                 # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text.y = element_text(size = 12, colour = "black"),                    # color the axis text
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text = element_text(size = 12,                                              # change the size of the axis titles
                              colour = "black"),                                    # change the color of the axis titles
          legend.background = element_rect(fill = NA),
          legend.position = c(.85, .1),
          legend.box = "horizontal") +
    facet_grid(~Stage) + 
    labs(title = "Nest Survival",
         color = "Year",
         x = "Litter Depth (mm)",
         y = "Daily Survival Rate"))


ggsave(CCSP.plot,
       filename = "outputs/figs/betaCCSP.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(CCSPstage.plot,
       filename = "outputs/figs/stageCCSP.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(CCSPbare.plot,
       filename = "outputs/figs/bareCCSP.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)

ggsave(CCSPlitd.plot,
       filename = "outputs/figs/litdCCSP.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 6)


# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute "rm(list = ls(all = TRUE))" - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all = TRUE))
# Then, execute "cleanup(ask = FALSE)" to delete orphaned output
#  files from MARK. Execute "?cleanup" to learn more
cleanup(ask = FALSE)
