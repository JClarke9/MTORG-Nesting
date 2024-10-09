# Load libraries ----------------------------------------------------------


library(tidyverse)
library(RMark)


# Read in Data ------------------------------------------------------------


nest <- read.csv("working/RMarknesting.csv")

windowsFonts(my_font = windowsFont("Gandhi Sans"))


# Data Wrangling ----------------------------------------------------------


nest$Nestling <- factor(nest$Nestling,
                        levels = c("0", "1"))

nest$Year <- factor(nest$Year,
                    levels = c("2021", "2022", "2023"))

nest$cTreat <- factor(nest$cTreat,
                      levels = c("0", "39", "49", "68"))

nest$Replicate <- factor(nest$Replicate,
                         levels = c("NE", "SE", "SW", "NW"))

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


# Constant Survival Models ------------------------------------------------


WEME.constant <- mark(WEME.surv, 
                      nocc = max(WEME.surv$LastChecked), 
                      model = "Nest",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula =  ~1)))

WEME.year <- mark(WEME.surv, 
                  nocc = max(WEME.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + Year)))

WEME.stage <- mark(WEME.surv, 
                   nocc = max(WEME.surv$LastChecked), 
                   model = "Nest",
                   groups = "Nestling",
                   adjust = FALSE,
                   delete = TRUE,
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))


BRBL.constant <- mark(BRBL.surv, 
                      nocc = max(BRBL.surv$LastChecked), 
                      model = "Nest",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula =  ~1)))

BRBL.year <- mark(BRBL.surv, 
                  nocc = max(BRBL.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + Year)))

BRBL.stage <- mark(BRBL.surv, 
                   nocc = max(BRBL.surv$LastChecked), 
                   model = "Nest",
                   groups = "Nestling",
                   adjust = FALSE,
                   delete = TRUE,
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))


CCSP.constant <- mark(CCSP.surv, 
                      nocc = max(CCSP.surv$LastChecked), 
                      model = "Nest",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula =  ~1)))

CCSP.year <- mark(CCSP.surv, 
                  nocc = max(CCSP.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + Year)))

CCSP.stage <- mark(CCSP.surv, 
                   nocc = max(CCSP.surv$LastChecked), 
                   model = "Nest",
                   groups = "Nestling",
                   adjust = FALSE,
                   delete = TRUE,
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))


MODO.constant <- mark(MODO.surv, 
                      nocc = max(MODO.surv$LastChecked), 
                      model = "Nest",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula =  ~1)))

MODO.year <- mark(MODO.surv, 
                  nocc = max(MODO.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + Year)))

MODO.stage <- mark(MODO.surv, 
                   nocc = max(MODO.surv$LastChecked), 
                   model = "Nest",
                   groups = "Nestling",
                   adjust = FALSE,
                   delete = TRUE,
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))


RWBL.constant <- mark(RWBL.surv, 
                      nocc = max(RWBL.surv$LastChecked), 
                      model = "Nest",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula =  ~1)))

RWBL.year <- mark(RWBL.surv, 
                  nocc = max(RWBL.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + Year)))

RWBL.stage <- mark(RWBL.surv, 
                   nocc = max(RWBL.surv$LastChecked), 
                   model = "Nest",
                   groups = "Nestling",
                   adjust = FALSE,
                   delete = TRUE,
                   model.parameters = list(S = list(formula =  ~1 + Nestling)))


GADW.constant <- mark(GADW.surv, 
                      nocc = max(GADW.surv$LastChecked), 
                      model = "Nest",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula =  ~1)))

GADW.year <- mark(GADW.surv, 
                  nocc = max(GADW.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + Year)))


BWTE.constant <- mark(BWTE.surv, 
                      nocc = max(BWTE.surv$LastChecked), 
                      model = "Nest",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula =  ~1)))

BWTE.year <- mark(BWTE.surv, 
                  nocc = max(BWTE.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + Year)))


NOPI.constant <- mark(NOPI.surv, 
                      nocc = max(NOPI.surv$LastChecked), 
                      model = "Nest",
                      adjust = FALSE,
                      delete = TRUE,
                      model.parameters = list(S = list(formula =  ~1)))

NOPI.year <- mark(NOPI.surv, 
                  nocc = max(NOPI.surv$LastChecked), 
                  model = "Nest",
                  groups = "Year",
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + Year)))


# Creating the table ------------------------------------------------------


# List all objects in the environment
all_objects <- ls()

# Identify dataframes ending with ".trt"
dataframes_to_keep <- all_objects[sapply(all_objects, function(df) {
  inherits(get(df), "mark") && grepl("\\.constant$", df) | grepl("\\.year$", df) | 
    grepl("\\.stage$", df)})]

real <- data.frame()

for(df_name in dataframes_to_keep) {
  
  spec <- if (endsWith(df_name, ".constant")) {
    gsub(".constant", "", df_name)
  } else if (endsWith(df_name, ".year")) {
    gsub(".year", "", df_name)
  } else if(endsWith(df_name, ".stage")) {
    gsub(".stage", "", df_name)
  } else(stop("Model does not end in '.constant', '.year', or '.stage'."))
  
  df <- get(df_name)
  
  spec.real <- if (endsWith(df_name, ".constant")) {
    as.data.frame(df$results$real) |> 
      rownames_to_column(var = "Group") |> 
      mutate(Variable = "Constant",
             Species = spec) |> 
      select(Variable, Species, estimate, se, lcl, ucl)
  } else if(endsWith(df_name, ".year")) {
    as.data.frame(df$results$real) |> 
      rownames_to_column(var = "Group") |> 
      mutate(Variable = case_when(
        grepl("g2021", Group) ~ "2021",
        grepl("g2022", Group) ~ "2022",
        grepl("g2023", Group) ~ "2023"),
        Species = spec) |> 
      select(Variable, Species, estimate, se, lcl, ucl)
  } else if(endsWith(df_name, ".stage")) {
    as.data.frame(df$results$real) |> 
      rownames_to_column(var = "Group") |> 
      mutate(Variable = case_when(
        grepl("g0", Group) ~ "Incubating",
        grepl("g1", Group) ~ "Nestling"),
        Species = spec) |> 
      select(Variable, Species, estimate, se, lcl, ucl)}
  else {stop("Model does not end in '.constant', '.year', or '.stage'.")}
  
  real <- rbind(real, spec.real)
}

wide_real <- pivot_wider(real[,c(1:3)],
                         names_from = Variable,
                         values_from = estimate)
