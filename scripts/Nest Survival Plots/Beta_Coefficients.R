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
                    levels = c("2021", "2022", "2023", "2024"))

# This loop creates a new data frame for each species and removes any
# dataframes from the environment that aren"t over 30 observations.


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


WEME.trt <- mark(WEME.surv, 
                 nocc = max(WEME.surv$LastChecked), 
                 model = "Nest", 
                 groups = "Nestling",
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Nestling)))

BRBL.trt <- mark(BRBL.surv, 
                 nocc = max(BRBL.surv$LastChecked), 
                 model = "Nest", 
                 groups = c("Year",
                            "Nestling"),
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Year + Nestling + Litter)))

CCSP.trt <- mark(CCSP.surv, 
                 nocc = max(CCSP.surv$LastChecked), 
                 model = "Nest", 
                 groups = c("Year",
                            "Nestling"),
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Year + Nestling + LitterD + Bare)))

MODO.trt <- mark(MODO.surv, 
                 nocc = max(MODO.surv$LastChecked), 
                 model = "Nest", 
                 groups = c("Year",
                            "Nestling"),
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Year + Time + I(Time^2) + Nestling + KBG)))

RWBL.trt <- mark(RWBL.surv, 
                 nocc = max(RWBL.surv$LastChecked), 
                 model = "Nest", 
                 groups = "Nestling",
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Time + Nestling + BHCONum)))

GADW.trt <- mark(GADW.surv, 
                 nocc = max(GADW.surv$LastChecked), 
                 model = "Nest", 
                 groups = "Year",
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Year + NestAge + grazed + pDoD + Forb)))

NOPI.trt <- mark(NOPI.surv, 
                 nocc = max(NOPI.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Litter + VOR)))

BWTE.trt <- mark(BWTE.surv, 
                 nocc = max(BWTE.surv$LastChecked), 
                 model = "Nest",
                 groups = "Year",
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Year + NestAge + Veg.Height)))

MALL.trt <- mark(MALL.surv, 
                 nocc = max(MALL.surv$LastChecked), 
                 model = "Nest",
                 groups = "Year",
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + Litter)))

NSHO.trt <- mark(NSHO.surv, 
                 nocc = max(NSHO.surv$LastChecked), 
                 model = "Nest",
                 groups = "Year",
                 adjust = FALSE,
                 delete = TRUE, 
                 model.parameters = list(S = list(formula =  ~1 + SmoothB)))


# Pull out Beta coefficients ----------------------------------------------


WEME.beta <- coef(WEME.trt) |>
  cbind(confint(WEME.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "WEME")

BRBL.beta <- coef(BRBL.trt) |>
  cbind(confint(BRBL.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "BRBL")

CCSP.beta <- coef(CCSP.trt) |>
  cbind(confint(CCSP.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "CCSP")

MODO.beta <- coef(MODO.trt) |>
  cbind(confint(MODO.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "MODO")

RWBL.beta <- coef(RWBL.trt) |>
  cbind(confint(RWBL.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "RWBL")

GADW.beta <- coef(GADW.trt) |>
  cbind(confint(GADW.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "GADW")

NOPI.beta <- coef(NOPI.trt) |>
  cbind(confint(NOPI.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "NOPI")

BWTE.beta <- coef(BWTE.trt) |>
  cbind(confint(BWTE.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "BWTE")

MALL.beta <- coef(MALL.trt) |>
  cbind(confint(MALL.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "MALL")

NSHO.beta <- coef(NSHO.trt) |>
  cbind(confint(NSHO.trt, level = 0.85)) |> 
  select(estimate, `7.5 %`, `92.5 %`) |> 
  rownames_to_column(var = "Variable") |> 
  rename(c("Coefficient" = "estimate",
           "lcl" = "7.5 %",
           "ucl" = "92.5 %")) |> 
  mutate(Species = "NSHO")

beta <- bind_rows(WEME.beta,
                  BRBL.beta,
                  CCSP.beta, 
                  MODO.beta,
                  RWBL.beta,
                  GADW.beta,
                  NOPI.beta,
                  BWTE.beta,
                  MALL.beta,
                  NSHO.beta)

unique(beta$Variable)

beta$Variable <- case_match(beta$Variable,
                            "S:(Intercept)" ~ "Intercept",
                            "S:Nestling1" ~ "Nestling",
                            "S:Year2022" ~ "Year 2022",
                            "S:Year2023" ~ "Year 2023",
                            "S:Year2024" ~ "Year 2024",
                            "S:Litter" ~ "Litter Cover",
                            "S:LitterD" ~ "Litter Depth",
                            "S:Bare" ~ "Bare Ground Cover",
                            "S:Time" ~ "Time",
                            "S:I(Time^2)" ~ "Time^2",
                            "S:KBG" ~ "Kentucky Bluegrass Cover",
                            "S:NestAge" ~ "Nest Age",
                            "S:grazed" ~ "Days Grazed",
                            "S:pDoD" ~ "Prior Grazing Intensity",
                            "S:Forb" ~ "Forb Cover",
                            "S:VOR" ~ "Veg Density",
                            "S:Veg.Height" ~ "Veg Height",
                            "S:SmoothB" ~ "Smooth Brome Cover")

beta$Variable <- factor(beta$Variable,
                        levels = c("Intercept", "Year 2022", "Year 2023", "Year 20224", "Time",
                                   "Time^2", "Nest Age", "Nestling", "Prior Grazing Intensity", "Days Grazed",
                                   "Kentucky Bluegrass Cover", "Forb Cover", "Smooth Brome Cover", "Litter Cover",
                                   "Bare Ground Cover", "Litter Depth", "Veg Height", "Veg Density"))


# Plot the data -----------------------------------------------------------


(betaF.plot <- ggplot(beta, 
                     aes(x = Variable,
                         y = Coefficient,
                         fill = factor(Species,
                                       levels = c("WEME", "BRBL", "MODO", "CCSP", "RWBL",
                                                  "GADW", "BWTE", "NOPI", "MALL", "NSHO")))) +
    geom_hline(yintercept = 0,
               colour = gray(1/2), 
               lty = 2) +
    geom_bar(position = "dodge",
             stat = "identity",
             colour = "black",
             width = 0.7) +
    scale_fill_manual(values = c("#D4A634", 
                                 "#D4A634", 
                                 "#D4A634", 
                                 "#D4A634", 
                                 "#D4A634",
                                 "#A2A4A2", 
                                 "#A2A4A2", 
                                 "#A2A4A2", 
                                 "#A2A4A2", 
                                 "#A2A4A2")) +
    geom_errorbar(aes(ymin = lcl,
                      ymax = ucl),
                  position = position_dodge(0.7),
                  width = 0.25,
                  colour = "black") +
    theme(plot.title = element_text(family = "my_font",
                                    hjust = 0.5,
                                    size = 40,
                                    vjust = 1,
                                    colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA,
                                          colour = NA),
          plot.background = element_rect(fill = NA,
                                         colour = NA),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 30, 
                                   colour = "black"),
          axis.ticks = element_line(colour = "black"),
          text = element_text(size = 30,
                              colour = "black")) +
    labs(title = "Top Model Effect Sizes",
         x = NULL,
         fill = "Species",
         y = expression("Beta " (beta))))


beta_f <- filter(beta, Variable != "Intercept" & Variable != "Year 2022" & Variable != "Year 2023" &
                   Variable != "Time" & Variable != "Time^2" & Variable != "Nestling" & Variable != "Nest Age")

beta_f$Species <- factor(beta_f$Species,
                         levels = c("BRBL", "MODO", "CCSP",
                                    "GADW", "BWTE", "NOPI", "MALL", "NSHO"))

beta_f$Type <- ifelse(beta_f$Species == "GADW", "FAC",
                      ifelse(beta_f$Species == "BWTE", "FAC",
                             ifelse(beta_f$Species == "NOPI", "FAC",
                                    ifelse(beta_f$Species == "MALL", "FAC",
                                           ifelse(beta_f$Species == "NSHO", "FAC",
                                                  ifelse(beta_f$Species == "RWBL", "FAC",
                                                         ifelse(beta_f$Species == "MODO", "FAC",
                                                                ifelse(beta_f$Species == "BRBL", "OBL",
                                                                       ifelse(beta_f$Species == "CCSP", "FAC",
                                                                              ifelse(beta_f$Species == "WEME", "OBL",
                                                                                     NA))))))))))

beta_f$Type <- factor(beta_f$Type,
                      levels = c("FAC", "OBL"))

beta_f$Group <- ifelse(beta_f$Variable %in% c("Litter Cover", "Litter Depth", "Bare Ground Cover", "Veg Height", "Veg Density"),
                       "Structure",
                       ifelse(beta_f$Variable %in% c("Days Grazed", "Prior Grazing Intensity"),
                              "Grazing",
                              ifelse(beta_f$Variable %in% c("Kentucky Bluegrass Cover", "Forb Cover", "Smooth Brome Cover"),
                                                            "Composition",
                                                            NA)))

(beta.plotG <- ggplot(beta_f[beta_f$Group == "Grazing",], 
                     aes(x = Variable,
                         y = Coefficient,
                         fill = Species)) +
    geom_hline(yintercept = 0,
               colour = "black", 
               lty = 1,) +
    geom_bar(position = "dodge",
             stat = "identity",
             colour = "black",
             width = 0.7) +
    geom_errorbar(aes(ymin = lcl,
                      ymax = ucl),
                  position = position_dodge(0.7),
                  width = 0.5,
                  linewidth = 0.7,
                  colour = "black") +
    scale_fill_manual(values = c("GADW" = "#613323")) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",
                                          color = NA), 
          plot.background = element_rect(fill = "transparent",
                                         color = NA),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(family = "my_font",
                                    hjust = .5,
                                    vjust = 1,
                                    size = 40,
                                    color = "black"),
          axis.title.x = element_text(family = "my_font",
                                      size = 28,
                                      color = "black",
                                      vjust = 1),
          axis.title.y = element_text(family = "my_font",
                                      size = 28,
                                      color = "black",
                                      angle = 90,
                                      vjust = 1),
          axis.text.x = element_text(family = "my_font",
                                     size = 28,
                                     color = "black"),
          axis.text.y = element_text(family = "my_font",
                                     size = 28,
                                     color = "black"),
          legend.background = element_blank(),
          legend.position = "none") +
    coord_cartesian(ylim = c(-0.025, 0.10)) +
    labs(title = "Grazing Impacts",
         x = NULL,
         fill = "Type",
         y = expression("Beta " (beta))))

ggsave(beta.plotG,
       filename = "outputs/figs/betaGrazing.png",
       bg = "transparent",
       dpi = 600,
       height = 6.34,
       width = 8.78)

(beta.plotC <- ggplot(beta_f[beta_f$Group == "Composition",], 
                     aes(x = Variable,
                         y = Coefficient,
                         fill = Species)) +
    geom_hline(yintercept = 0,
               colour = "black", 
               lty = 1,) +
    geom_bar(position = "dodge",
             stat = "identity",
             colour = "black",
             width = 0.7) +
    geom_errorbar(aes(ymin = lcl,
                      ymax = ucl),
                  position = position_dodge(0.7),
                  width = 0.5,
                  linewidth = 0.7,
                  colour = "black") +
    scale_fill_manual(values = c("MODO" = "#C8A696", 
                                 "GADW" = "#613323",
                                 "NSHO" = "#AF4A05")) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",
                                          color = NA), 
          plot.background = element_rect(fill = "transparent",
                                         color = NA),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(family = "my_font",
                                    hjust = .5,
                                    vjust = 1,
                                    size = 40,
                                    color = "black"),
          axis.title.x = element_text(family = "my_font",
                                      size = 28,
                                      color = "black",
                                      vjust = 1),
          axis.title.y = element_text(family = "my_font",
                                      size = 28,
                                      color = "black",
                                      angle = 90,
                                      vjust = 1),
          axis.text.x = element_text(family = "my_font",
                                     size = 28,
                                     color = "black"),
          axis.text.y = element_text(family = "my_font",
                                     size = 28,
                                     color = "black"),
          legend.background = element_blank(),
          legend.position = "none") +
    coord_cartesian(ylim = c(-0.025, 0.10)) +
    labs(title = "Vegetation Composition Impacts",
         x = NULL,
         fill = "Type",
         y = NULL))

  ggsave(beta.plotC,
       filename = "outputs/figs/betaComposition.png",
       bg = "transparent",
       dpi = 600,
       height = 6.34,
       width = 13.13)

(beta.plotS <- ggplot(beta_f[beta_f$Group == "Structure",], 
                     aes(x = Variable,
                         y = Coefficient,
                         fill = Species)) +
    geom_hline(yintercept = 0,
               colour = "black", 
               lty = 1,) +
    geom_bar(position = "dodge",
             stat = "identity",
             colour = "black",
             width = 0.7) +
    geom_errorbar(aes(ymin = lcl,
                      ymax = ucl),
                  position = position_dodge(0.7),
                  width = 0.5,
                  linewidth = 0.7,
                  colour = "black") +
    scale_fill_manual(values = c("BRBL" = "#364C59", 
                                 "NOPI" = "#736357",
                                 "MALL" = "#0E7B54",
                                 "CCSP" = "#ADA6B2",
                                 "BWTE" = "#9BAFD0")) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",
                                          color = NA), 
          plot.background = element_rect(fill = "transparent",
                                         color = NA),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(family = "my_font",
                                    hjust = .5,
                                    vjust = 1,
                                    size = 40,
                                    color = "black"),
          axis.title.x = element_text(family = "my_font",
                                      size = 28,
                                      color = "black",
                                      vjust = 1),
          axis.title.y = element_text(family = "my_font",
                                      size = 28,
                                      color = "black",
                                      angle = 90,
                                      vjust = 1),
          axis.text.x = element_text(family = "my_font",
                                     size = 28,
                                     color = "black"),
          axis.text.y = element_text(family = "my_font",
                                     size = 28,
                                     color = "black"),
          legend.background = element_blank(),
          legend.position = "none") +
    labs(title = "Vegetation Structure Impacts",
         x = NULL,
         fill = "Type",
         y = expression("Beta " (beta))))

ggsave(beta.plotS,
       filename = "outputs/figs/betaStructure.png",
       bg = "transparent",
       dpi = 600,
       height = 9.14,
       width = 21.88)

# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute "rm(list = ls(all = TRUE))" - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all = TRUE))
# Then, execute "cleanup(ask = FALSE)" to delete orphaned output
#  files from MARK. Execute "?cleanup" to learn more
cleanup(ask = FALSE)
