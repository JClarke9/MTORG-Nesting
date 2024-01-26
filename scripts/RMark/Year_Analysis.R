# Load libraries ----------------------------------------------------------


library(tidyverse)
library(RMark)
library(MuMIn)
library(ggpattern)


# Read in Data ------------------------------------------------------------


nest <- read.csv("working/RMarknesting.csv")

windowsFonts(my_font = windowsFont("Gandi Sans"))                          # downloading a font to be used for my ordination


# Data Wrangling ----------------------------------------------------------


nest$Nestling <- factor(nest$Nestling,
                        levels = c("0", "1"))

nest$Year <- factor(nest$Year,
                    levels = c("2021", "2022", "2023"))

# This loop creates a new data frame for each species and removes any
# dataframes from the environment that aren't over 30 observations.


for (i in unique(nest$Spec)) {
  
  #Create a dataframe name based on the cleaned species name with ".surv" appended
  df_name <- paste0(i, ".surv")
  
  # Filter the data based on the condition for the current species
  current_df <- filter(nest, Spec == i & Stage != "Laying")
  
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

# Identify dataframes ending with ".surv" and having more than 30 observations
surv_dataframes_to_keep <- all_objects[sapply(all_objects, function(df) {
  inherits(get(df), "data.frame") && grepl("\\.surv$", df) && nrow(get(df)) > 30
})]

# Remove all objects except those meeting the criteria
objects_to_remove <- setdiff(all_objects, surv_dataframes_to_keep)
if (length(objects_to_remove) > 0) {
  rm(list = objects_to_remove)
}


# Run the top models ------------------------------------------------------


WEME.year <- mark(WEME.surv, 
                  nocc = max(WEME.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE, 
                  model.parameters = list(S = list(formula =  ~1 + Year)))

BRBL.year <- mark(BRBL.surv, 
                  nocc = max(BRBL.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE, 
                  model.parameters = list(S = list(formula =  ~1 + Year)))

CCSP.year <- mark(CCSP.surv, 
                  nocc = max(CCSP.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE, 
                  model.parameters = list(S = list(formula =  ~1 + Year)))

MODO.year <- mark(MODO.surv, 
                  nocc = max(MODO.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE, 
                  model.parameters = list(S = list(formula =  ~1 + Year)))

RWBL.year <- mark(RWBL.surv, 
                  nocc = max(RWBL.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE, 
                  model.parameters = list(S = list(formula =  ~1 + Year)))

GADW.year <- mark(GADW.surv, 
                  nocc = max(GADW.surv$LastChecked), 
                  model = "Nest", 
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE, 
                  model.parameters = list(S = list(formula =  ~1 + Year)))

NOPI.year <- mark(NOPI.surv, 
                  nocc = max(NOPI.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE, 
                  model.parameters = list(S = list(formula =  ~1 + Year)))

BWTE.year <- mark(BWTE.surv, 
                  nocc = max(BWTE.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE, 
                  model.parameters = list(S = list(formula =  ~1 + Year)))


# Pull out Beta coefficients ----------------------------------------------


WEME.Ybeta <- coef(WEME.year) |>
  cbind(confint(WEME.year, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "WEME")

BRBL.Ybeta <- coef(BRBL.year) |>
  cbind(confint(BRBL.year, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "BRBL")

CCSP.Ybeta <- coef(CCSP.year) |>
  cbind(confint(CCSP.year, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "CCSP")

MODO.Ybeta <- coef(MODO.year) |>
  cbind(confint(MODO.year, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "MODO")

RWBL.Ybeta <- coef(RWBL.year) |>
  cbind(confint(RWBL.year, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "RWBL")

GADW.Ybeta <- coef(GADW.year) |>
  cbind(confint(GADW.year, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "GADW")

NOPI.Ybeta <- coef(NOPI.year) |>
  cbind(confint(NOPI.year, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "NOPI")

BWTE.Ybeta <- coef(BWTE.year) |>
  cbind(confint(BWTE.year, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "BWTE")

year.beta <- bind_rows(WEME.Ybeta,
                       BRBL.Ybeta,
                       CCSP.Ybeta, 
                       MODO.Ybeta,
                       RWBL.Ybeta,
                       GADW.Ybeta,
                       NOPI.Ybeta,
                       BWTE.Ybeta)

unique(year.beta$Variable)

year.beta$Variable <- case_match(year.beta$Variable,
                                 "S:(Intercept)" ~ "Intercept",
                                 "S:Year2022" ~ "Year 2022",
                                 "S:Year2023" ~ "Year 2023")

year.beta$Variable <- factor(year.beta$Variable,
                             levels = c("Intercept", "Year 2022", "Year 2023"))

(Ybeta.plot <- ggplot(year.beta, 
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
                                    hjust = 0.5,
                                    size = 20,
                                    vjust = 1,
                                    colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA,
                                          colour = NA),
          plot.background = element_rect(fill = NA,
                                         colour = NA),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 12, 
                                   colour = "black"),
          axis.ticks = element_line(colour = "black"),
          text = element_text(size = 12,
                              colour = "black")) +
    facet_wrap(~Species) +
    labs(title = "Year Effects",
         x = NULL,
         y = expression("Beta " (beta))))


# Calculating DSR ---------------------------------------------------------


(WEME.real <- as.data.frame(WEME.year$results$real) |> 
   rownames_to_column(var = "Group") |> 
   mutate(Year = case_when(
     grepl("2021", Group) ~ "2021",
     grepl("2022", Group) ~ "2022",
     grepl("2023", Group) ~ "2023")) |> 
   select(Year, estimate, se, lcl, ucl) |> 
   group_by(Year) |> 
   summarise(estimate = mean(estimate),
             se = mean(se),
             lcl = mean(lcl),
             ucl = mean(ucl),
             Species = "WEME"))

(BRBL.real <- as.data.frame(BRBL.year$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023")) |> 
    select(Year, estimate, se, lcl, ucl) |> 
    group_by(Year) |> 
    summarise(estimate = mean(estimate),
              se = mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl),
              Species = "BRBL"))

(CCSP.real <- as.data.frame(CCSP.year$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023")) |> 
    select(Year, estimate, se, lcl, ucl) |> 
    group_by(Year) |> 
    summarise(estimate = mean(estimate),
              se = mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl),
              Species = "CCSP"))

(MODO.real <- as.data.frame(MODO.year$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023")) |> 
    select(Year, estimate, se, lcl, ucl) |> 
    group_by(Year) |> 
    summarise(estimate = mean(estimate),
              se = mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl),
              Species = "MODO"))

(RWBL.real <- as.data.frame(RWBL.year$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023")) |> 
    select(Year, estimate, se, lcl, ucl) |> 
    group_by(Year) |> 
    summarise(estimate = mean(estimate),
              se = mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl),
              Species = "RWBL"))

(GADW.real <- as.data.frame(GADW.year$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023")) |> 
    select(Year, estimate, se, lcl, ucl) |> 
    group_by(Year) |> 
    summarise(estimate = mean(estimate),
              se = mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl),
              Species = "GADW"))

(NOPI.real <- as.data.frame(NOPI.year$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023")) |> 
    select(Year, estimate, se, lcl, ucl) |> 
    group_by(Year) |> 
    summarise(estimate = mean(estimate),
              se = mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl),
              Species = "NOPI"))

(BWTE.real <- as.data.frame(BWTE.year$results$real) |> 
    rownames_to_column(var = "Group") |> 
    mutate(Year = case_when(
      grepl("2021", Group) ~ "2021",
      grepl("2022", Group) ~ "2022",
      grepl("2023", Group) ~ "2023")) |> 
    select(Year, estimate, se, lcl, ucl) |> 
    group_by(Year) |> 
    summarise(estimate = mean(estimate),
              se = mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl),
              Species = "BWTE"))

real.dsr <- bind_rows(WEME.real,
                      BRBL.real,
                      CCSP.real, 
                      MODO.real,
                      RWBL.real,
                      GADW.real,
                      NOPI.real,
                      BWTE.real)

(Ydsr.plot <- ggplot(real.dsr, 
                     aes(x = Year,
                         y = estimate)) +
    geom_point(aes(x = Year,
                   y = estimate),
               size = 4) +
    geom_errorbar(aes(x = Year,
                      ymin = lcl,
                      ymax = ucl),
                  width = .5,
                  linewidth = 1) +
    theme(plot.title = element_text(family = "my_font",
                                    hjust = 0.5,
                                    size = 20,
                                    vjust = 1,
                                    colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA,
                                          colour = NA),
          plot.background = element_rect(fill = NA,
                                         colour = NA),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 12, 
                                   colour = "black"),
          axis.ticks = element_line(colour = "black"),
          text = element_text(size = 12,
                              colour = "black")) +
    facet_wrap(~Species) +
    labs(title = "DSR Across Years",
         x = NULL,
         y = "Daily Survival Rate"))

ggsave("outputs/figs/Beta_year.png",
       Ybeta.plot,
       bg = "white",
       dpi = 600)

ggsave("outputs/figs/DSR_year.png",
       Ydsr.plot,
       bg = "white",
       dpi = 600)
