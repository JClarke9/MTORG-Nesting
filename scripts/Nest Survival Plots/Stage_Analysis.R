# Load libraries ----------------------------------------------------------


library(tidyverse)
library(RMark)
library(MuMIn)
library(ggpattern)


# Read in Data ------------------------------------------------------------


nest <- read.csv('working/RMarknesting.csv') |> 
  mutate(Nestling = factor(Nestling, levels = c('0', '1')),
         Year = factor(Year, levels = c('2021', '2022', '2023', '2024')))


# Data Wrangling ----------------------------------------------------------


# This loop creates a new data frame for each species and removes any
# dataframes from the environment that aren't over 30 observations.


for (i in unique(nest$Spec)) {
  
  #Create a dataframe name based on the cleaned species name with '.surv' appended
  df_name <- paste0(i, '.surv')
  
  # Filter the data based on the condition for the current species
  current_df <- filter(nest, Spec == i & Stage != 'Laying')
  
  # Remove NA values from the dataframe
  MISSING <- is.na(current_df$AgeFound)
  
  sum(MISSING)
  
  current_df <- subset(current_df, 
                       subset = !MISSING)
  
  # Assign the cleaned and filtered dataframe to the dynamically generated name
  assign(df_name, current_df)
}


# List all objects in the environment
all_objects <- ls()

# Identify dataframes ending with '.surv' and having more than 30 observations
surv_dataframes_to_keep <- all_objects[sapply(all_objects, function(df) {
  inherits(get(df), 'data.frame') && grepl('\\.surv$', df) && nrow(get(df)) > 30
})]

# Remove all objects except those meeting the criteria
objects_to_remove <- setdiff(all_objects, surv_dataframes_to_keep)

if (length(objects_to_remove) > 0) {
  rm(list = objects_to_remove)
}


# Run the top models ------------------------------------------------------


WEME.stage <- mark(WEME.surv, 
                   nocc = max(WEME.surv$LastChecked), 
                   model = 'Nest',
                   groups = 'Nestling',
                   adjust = FALSE,
                   delete = TRUE, 
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))

BRBL.stage <- mark(BRBL.surv, 
                   nocc = max(BRBL.surv$LastChecked), 
                   model = 'Nest',
                   groups = 'Nestling',
                   adjust = FALSE,
                   delete = TRUE, 
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))

CCSP.stage <- mark(CCSP.surv, 
                   nocc = max(CCSP.surv$LastChecked), 
                   model = 'Nest',
                   groups = 'Nestling',
                   adjust = FALSE,
                   delete = TRUE, 
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))

MODO.stage <- mark(MODO.surv, 
                   nocc = max(MODO.surv$LastChecked), 
                   model = 'Nest',
                   groups = 'Nestling',
                   adjust = FALSE,
                   delete = TRUE, 
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))

RWBL.stage <- mark(RWBL.surv, 
                   nocc = max(RWBL.surv$LastChecked), 
                   model = 'Nest',
                   groups = 'Nestling',
                   adjust = FALSE,
                   delete = TRUE, 
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))

YHBL.stage <- mark(YHBL.surv, 
                   nocc = max(YHBL.surv$LastChecked), 
                   model = 'Nest',
                   groups = 'Nestling',
                   adjust = FALSE,
                   delete = TRUE, 
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))


# Pull out Beta coefficients ----------------------------------------------


WEME.Sbeta <- coef(WEME.stage) |>
  cbind(confint(WEME.stage, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = 'Variable') |> 
  rename(c('Coefficient' = 'estimate',
           'lcl' = '7.5 %',
           'ucl' = '92.5 %')) |> 
  mutate(Species = 'WEME')

BRBL.Sbeta <- coef(BRBL.stage) |>
  cbind(confint(BRBL.stage, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = 'Variable') |> 
  rename(c('Coefficient' = 'estimate',
           'lcl' = '7.5 %',
           'ucl' = '92.5 %')) |> 
  mutate(Species = 'BRBL')

CCSP.Sbeta <- coef(CCSP.stage) |>
  cbind(confint(CCSP.stage, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = 'Variable') |> 
  rename(c('Coefficient' = 'estimate',
           'lcl' = '7.5 %',
           'ucl' = '92.5 %')) |> 
  mutate(Species = 'CCSP')

MODO.Sbeta <- coef(MODO.stage) |>
  cbind(confint(MODO.stage, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = 'Variable') |> 
  rename(c('Coefficient' = 'estimate',
           'lcl' = '7.5 %',
           'ucl' = '92.5 %')) |> 
  mutate(Species = 'MODO')

RWBL.Sbeta <- coef(RWBL.stage) |>
  cbind(confint(RWBL.stage, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = 'Variable') |> 
  rename(c('Coefficient' = 'estimate',
           'lcl' = '7.5 %',
           'ucl' = '92.5 %')) |> 
  mutate(Species = 'RWBL')

YHBL.Sbeta <- coef(YHBL.stage) |>
  cbind(confint(YHBL.stage, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = 'Variable') |> 
  rename(c('Coefficient' = 'estimate',
           'lcl' = '7.5 %',
           'ucl' = '92.5 %')) |> 
  mutate(Species = 'YHBL')


stage.beta <- bind_rows(WEME.Sbeta,
                        BRBL.Sbeta,
                        CCSP.Sbeta, 
                        MODO.Sbeta,
                        RWBL.Sbeta,
                        YHBL.Sbeta)

unique(stage.beta$Variable)

stage.beta$Variable <- case_match(stage.beta$Variable,
                                  'S:(Intercept)' ~ 'Intercept',
                                  'S:Nestling1' ~ 'Nestling')

stage.beta$Variable <- factor(stage.beta$Variable,
                              levels = c('Intercept', 'Nestling'))


# Plotting Nest Stage Effect Size -----------------------------------------


windowsFonts(my_font = windowsFont('Gandi Sans'))                          # downloading a font to be used for my ordination


(Sbeta.plot <- ggplot(stage.beta, 
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
                                    hjust = 0.5,
                                    size = 20,
                                    vjust = 1,
                                    colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA,
                                          colour = NA),
          plot.background = element_rect(fill = NA,
                                         colour = NA),
          axis.line = element_line(colour = 'black'),
          axis.text = element_text(size = 12, 
                                   colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          text = element_text(size = 12,
                              colour = 'black')) +
    facet_wrap(~Species) +
    labs(title = 'Stage Effects',
         x = NULL,
         y = expression('Beta ' (beta))))


# Calculating DSR ---------------------------------------------------------


(WEME.real <- as.data.frame(WEME.stage$results$real) |> 
   rownames_to_column(var = 'Group') |> 
   mutate(Stage = case_when(
     grepl('S g0 a0 t1', Group) ~ 'Incubating',
     grepl('S g1 a0 t1', Group) ~ 'Nestling')) |> 
   group_by(Stage) |> 
   summarize(estimate = mean(estimate), 
             se = mean(se), 
             lcl = mean(lcl), 
             ucl = mean(ucl),
             Species = 'WEME'))

(BRBL.real <- as.data.frame(BRBL.stage$results$real) |> 
    rownames_to_column(var = 'Group') |> 
    mutate(Stage = case_when(
      grepl('S g0 a0 t1', Group) ~ 'Incubating',
      grepl('S g1 a0 t1', Group) ~ 'Nestling')) |> 
    group_by(Stage) |> 
    summarize(estimate = mean(estimate), 
              se = mean(se), 
              lcl = mean(lcl), 
              ucl = mean(ucl),
              Species = 'BRBL'))

(CCSP.real <- as.data.frame(CCSP.stage$results$real) |> 
    rownames_to_column(var = 'Group') |> 
    mutate(Stage = case_when(
      grepl('S g0 a0 t1', Group) ~ 'Incubating',
      grepl('S g1 a0 t1', Group) ~ 'Nestling')) |> 
    group_by(Stage) |> 
    summarize(estimate = mean(estimate), 
              se = mean(se), 
              lcl = mean(lcl), 
              ucl = mean(ucl),
              Species = 'CCSP'))

(MODO.real <- as.data.frame(MODO.stage$results$real) |> 
    rownames_to_column(var = 'Group') |> 
    mutate(Stage = case_when(
      grepl('S g0 a0 t1', Group) ~ 'Incubating',
      grepl('S g1 a0 t1', Group) ~ 'Nestling')) |> 
    group_by(Stage) |> 
    summarize(estimate = mean(estimate), 
              se = mean(se), 
              lcl = mean(lcl), 
              ucl = mean(ucl),
              Species = 'MODO'))

(RWBL.real <- as.data.frame(RWBL.stage$results$real) |> 
    rownames_to_column(var = 'Group') |> 
    mutate(Stage = case_when(
      grepl('S g0 a0 t1', Group) ~ 'Incubating',
      grepl('S g1 a0 t1', Group) ~ 'Nestling')) |> 
    group_by(Stage) |> 
    summarize(estimate = mean(estimate), 
              se = mean(se), 
              lcl = mean(lcl), 
              ucl = mean(ucl),
              Species = 'RWBL'))

(YHBL.real <- as.data.frame(YHBL.stage$results$real) |> 
    rownames_to_column(var = 'Group') |> 
    mutate(Stage = case_when(
      grepl('S g0 a0 t1', Group) ~ 'Incubating',
      grepl('S g1 a0 t1', Group) ~ 'Nestling')) |> 
    group_by(Stage) |> 
    summarize(estimate = mean(estimate), 
              se = mean(se), 
              lcl = mean(lcl), 
              ucl = mean(ucl),
              Species = 'YHBL'))


real.dsr <- bind_rows(WEME.real,
                      BRBL.real,
                      CCSP.real, 
                      MODO.real,
                      RWBL.real,
                      YHBL.real)


# Plotting Stage DSR ------------------------------------------------------


windowsFonts(my_font = windowsFont('Gandi Sans'))                          # downloading a font to be used for my ordination


(Sdsr.plot <- ggplot(real.dsr, 
                     aes(x = Stage,
                         y = estimate)) +
    geom_point(aes(x = Stage,
                   y = estimate),
               size = 4) +
    geom_errorbar(aes(x = Stage,
                      ymin = lcl,
                      ymax = ucl),
                  width = .5,
                  linewidth = 1) +
    theme(plot.title = element_text(family = 'my_font',
                                    hjust = 0.5,
                                    size = 20,
                                    vjust = 1,
                                    colour = 'black'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA,
                                          colour = NA),
          plot.background = element_rect(fill = NA,
                                         colour = NA),
          axis.line = element_line(colour = 'black'),
          axis.text = element_text(size = 12, 
                                   colour = 'black'),
          axis.ticks = element_line(colour = 'black'),
          text = element_text(size = 12,
                              colour = 'black')) +
    facet_wrap(~Species) +
    labs(title = 'DSR Across Stages',
         x = NULL,
         y = 'Daily Survival Rate'))

ggsave('outputs/figs/Beta_stage.png',
       Sbeta.plot,
       bg = 'white',
       dpi = 600)

ggsave('outputs/figs/DSR_stage.png',
       Sdsr.plot,
       bg = 'white',
       dpi = 600)
