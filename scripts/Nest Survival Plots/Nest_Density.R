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
                    levels = c("2021", "2022", "2023", "2024"))

nest$cTreat <- factor(nest$cTreat,
                      levels = c("0", "39", "49", "68"))

nest$Replicate <- factor(nest$Replicate,
                         levels = c("NE", "SE", "SW", "NW"))

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


# Run constant survival models --------------------------------------------


WEME.trt <- mark(WEME.surv, 
                 nocc = max(WEME.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

BRBL.trt <- mark(BRBL.surv, 
                 nocc = max(BRBL.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

CCSP.trt <- mark(CCSP.surv, 
                 nocc = max(CCSP.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

MODO.trt <- mark(MODO.surv, 
                 nocc = max(MODO.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

RWBL.trt <- mark(RWBL.surv, 
                 nocc = max(RWBL.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

YHBL.trt <- mark(YHBL.surv, 
                 nocc = max(YHBL.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

GADW.trt <- mark(GADW.surv, 
                 nocc = max(GADW.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

BWTE.trt <- mark(BWTE.surv, 
                 nocc = max(BWTE.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

NOPI.trt <- mark(NOPI.surv, 
                 nocc = max(NOPI.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

MALL.trt <- mark(MALL.surv, 
                 nocc = max(MALL.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))

NSHO.trt <- mark(NSHO.surv, 
                 nocc = max(NSHO.surv$LastChecked), 
                 model = "Nest",
                 groups = c("Year", "Replicate"),
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1 + Year + Replicate)))


# Calculate densities -----------------------------------------------------


# This loop calculates the densities and confidence intervals for all
# of the models listed above

birds.trtY <- NULL  # Initialize the combined dataframe

# Main loop for multiple species
all_objects <- ls()
species_list <- all_objects[grepl("\\.surv$", all_objects)]
species_list <- gsub("\\.surv$", "", species_list)

for (species in species_list) {
  surv_obj <- paste0(species, ".surv")
  trt_obj <- paste0(species, ".trt")
  
  if (exists(surv_obj) && exists(trt_obj)) {
    
    # Data preparation within the loop
    real_dataframe <- as.data.frame(get(trt_obj)$results$real) |>
      rownames_to_column(var = "Group") |> 
      mutate(Group = str_replace_all(Group, "S g| a0 t1", "")) |> 
      mutate(Replicate = case_when(
        endsWith(Group, "NE") ~ "NE",
        endsWith(Group, "SE") ~ "SE",
        endsWith(Group, "SW") ~ "SW",
        endsWith(Group, "NW") ~ "NW"),
        Year = case_when(
          startsWith(Group, "2021") ~ "2021",
          startsWith(Group, "2022") ~ "2022",
          startsWith(Group, "2023") ~ "2023",
          startsWith(Group, "2024") ~ "2024"))
    
    real_dataframe <- complete(real_dataframe,
                               Replicate,
                               Year,
                               fill = list(estimate = 0,
                                           se = 0,
                                           lcl = 0,
                                           ucl = 0))
    
    real_dataframe$Year <- factor(real_dataframe$Year,
                                  level = c("2021", "2022", "2023", "2024"))
    
    real_dataframe$Replicate <- factor(real_dataframe$Replicate,
                                       levels = c("NE", "SE", "SW", "NW"))
    
    # Summarize the ".real" dataframe
    summarized_real <- real_dataframe |> 
      group_by(Replicate,
               Year) |> 
      summarize(estimate = mean(estimate),
                se = mean(se),
                lcl = mean(lcl),
                ucl = mean(ucl))
    
    group <- summarized_real |>
      mutate(
        estimate = case_when(
          Year == 2021 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2022 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2023 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2024 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2021 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2022 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2023 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2024 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2021 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2022 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2023 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2024 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2021 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2022 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2023 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2024 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW"]) * 64),
          estimate == 0 ~ 0),
        lcl = case_when(
          Year == 2021 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2022 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2023 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2024 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2021 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2022 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2023 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2024 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2021 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2022 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2023 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2024 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2021 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2022 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2023 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2024 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW"]) * 64),
          estimate == 0 ~ 0),
        ucl = case_when(
          Year == 2021 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2022 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2023 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2024 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2021 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2022 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2023 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2024 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2021 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2022 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2023 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2024 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2021 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2022 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2023 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2024 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW"]) * 64),
          estimate == 0 ~ 0),)
    
    # Combine with the existing birds.trt dataframe
    birds.trtY <- rbind(birds.trtY, group |> mutate(Species = species))
  } else {
    warning(paste0("Objects ", surv_obj, " and/or ", trt_obj, " not found in the environment. Skipping species ", species))
  }
}

(birdY.density <- birds.trtY |> 
    group_by(Species) |> 
    summarize(birds_ha = mean(estimate),
              se =  mean(se),
              lcl = mean(lcl),
              ucl = mean(ucl)))


# Plot densities ----------------------------------------------------------


(density.plotY <- ggplot(birds.trtY, 
                         aes(x = Species, 
                             y = estimate)) +
   stat_summary(geom = "point",
                fun = "mean",
                size = 5) +
   stat_summary(geom = "errorbar",
                fun.data = mean_cl_boot,
                width = 0.2) +
   # geom_text(aes(label = round(estimate, 2)), 
   #           vjust = -0.5, 
   #           size = 3, 
   #           position = position_dodge(width = 0.75)) +
   scale_y_continuous(breaks = c(0.0, 0.25, 0.50, 0.75, 1.25)) +
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
                             colour = "black")) +                                  # change the color of the axis titles
   labs(title = "Average Avian Densities",
        x = NULL, 
        y = "Nests Per Ha"))

ggsave(density.plotY,
       filename = "outputs/figs/AvianDensity_Year.png",
       dpi = "print",
       bg = "white",
       height = 13,
       width = 22.5)


# Run constant survival models --------------------------------------------


WEME.trtG <- mark(WEME.surv, 
                  nocc = max(WEME.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

BRBL.trtG <- mark(BRBL.surv, 
                  nocc = max(BRBL.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

CCSP.trtG <- mark(CCSP.surv, 
                  nocc = max(CCSP.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

MODO.trtG <- mark(MODO.surv, 
                  nocc = max(MODO.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

RWBL.trtG <- mark(RWBL.surv, 
                  nocc = max(RWBL.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

YHBL.trtG <- mark(YHBL.surv, 
                  nocc = max(YHBL.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

GADW.trtG <- mark(GADW.surv, 
                  nocc = max(GADW.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

BWTE.trtG <- mark(BWTE.surv, 
                  nocc = max(BWTE.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

NOPI.trtG <- mark(NOPI.surv, 
                  nocc = max(NOPI.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

MALL.trtG <- mark(MALL.surv, 
                  nocc = max(MALL.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))

NSHO.trtG <- mark(NSHO.surv, 
                  nocc = max(NSHO.surv$LastChecked), 
                  model = "Nest",
                  groups = c("cTreat", "Year", "Replicate"),
                  adjust = FALSE,
                  delete = TRUE,
                  model.parameters = list(S = list(formula =  ~1 + cTreat + Year + Replicate)))


# Calculate densities -----------------------------------------------------


# This loop calculates the densities and confidence intervals for all
# of the models listed above

birds.trtT <- NULL  # Initialize the combined dataframe

# Main loop for multiple species
all_objects <- ls()
species_list <- all_objects[grepl("\\.surv$", all_objects)]
species_list <- gsub("\\.surv$", "", species_list)

for (species in species_list) {
  surv_obj <- paste0(species, ".surv")
  trt_obj <- paste0(species, ".trtG")
  
  if (exists(surv_obj) && exists(trt_obj)) {
    # Data preparation within the loop
    real_dataframe <- as.data.frame(get(trt_obj)$results$real) |>
      rownames_to_column(var = "Group") |> 
      mutate(Group = str_replace_all(Group, "S g| a0 t1", "")) |> 
      mutate(
        cTreat = case_when(
          startsWith(Group, "0") ~ "Rest",
          startsWith(Group, "39") ~ "Moderate",
          startsWith(Group, "49") ~ "Full",
          startsWith(Group, "68") ~ "Heavy"),
        Year = case_when(
          grepl("2021", Group) ~ "2021",
          grepl("2022", Group) ~ "2022",
          grepl("2023", Group) ~ "2023",
          grepl("2024", Group) ~ "2024"),
        Replicate = case_when(
          endsWith(Group, "NE") ~ "NE",
          endsWith(Group, "SE") ~ "SE",
          endsWith(Group, "SW") ~ "SW",
          endsWith(Group, "NW") ~ "NW"))
    
    real_dataframe$cTreat <- factor(real_dataframe$cTreat,
                                    levels = c("Rest", "Moderate", "Full", "Heavy"))
    
    real_dataframe$Year <- factor(real_dataframe$Year,
                                  level = c("2021", "2022", "2023", "2024"))
    
    real_dataframe$Replicate <- factor(real_dataframe$Replicate,
                                       level = c("NE", "SE", "SW", "NW"))
    
    real_dataframe <- complete(real_dataframe,
                               nesting(Year, cTreat),
                               Replicate,
                               fill = list(estimate = 0,
                                           se = 0,
                                           lcl = 0,
                                           ucl = 0))
    
    real_dataframe <- unite(real_dataframe,
                            Group,
                            c(Replicate, Year, cTreat),
                            sep = "_",
                            remove = F)
    
    # Summarize the ".real" dataframe
    summarized_real <- real_dataframe |> 
      group_by(Group, 
               cTreat,
               Replicate,
               Year) |> 
      summarize(estimate = mean(estimate),
                se = mean(se),
                lcl = mean(lcl),
                ucl = mean(ucl))
    
    group <- summarized_real |>
      mutate(
        estimate = case_when(
          Year == 2021 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          
          estimate == 0 ~ 0),
        lcl = case_when(
          Year == 2021 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          estimate == 0 ~ 0),
        ucl = case_when(
          Year == 2021 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2022 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2023 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "68"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "0"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "39"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "49"]) * 64),
          Year == 2024 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2024" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "68"]) * 64),
          estimate == 0 ~ 0))
    
    # Combine with the existing birds.trt dataframe
    birds.trtT <- rbind(birds.trtT, group |> mutate(Species = species))
  } else {
    warning(paste0("Objects ", surv_obj, " and/or ", trt_obj, " not found in the environment. Skipping species ", species))
  }
}

birds.trtT <- birds.trtT |> 
  ungroup()

birdT.density <- birds.trtT |> 
  group_by(cTreat,
           Species) |> 
  summarize(birds_ha = mean(estimate),
            se =  mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))


# Plot densities ----------------------------------------------------------


(density.plotT <- ggplot(birds.trtT[birds.trtT$Species %in% c("CCSP", "RWBL", "NOPI", "GADW"),], 
                         aes(x = factor(Species,
                                        levels = c("CCSP", "RWBL", "NOPI", "GADW")), 
                             y = estimate,
                             color = factor(cTreat,
                                            levels = c("Rest", "Moderate",
                                                       "Full", "Heavy")))) +
   stat_summary(geom = "point",
                fun = "mean",
                size = 1,
                position = position_dodge(.5)) +
   stat_summary(geom = "errorbar",
                fun.data = mean_cl_boot,
                linewidth = 0.5,
                position = position_dodge(.5)) +
   scale_color_manual(values=c("#A2A4A2", "lightgoldenrod2", "#D4A634", "#717F5B")) +
   scale_x_discrete(labels = c("CCSP" = "Clay-colored \n Sparrow", 
                               "RWBL" = "Red-winged \n Blackbird", 
                               "NOPI" = "Northern \n Pintail", 
                               "GADW" = "Gadwall")) +
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(), 
         plot.background = element_blank(),
         axis.line = element_line(color = "black"),
         axis.ticks = element_line(colour = "black"),
         plot.title = element_text(family = "my_font",
                                   hjust = .5,
                                   vjust = 1,
                                   size = 16,
                                   color = "black"),
         axis.title.x = element_blank(),
         axis.title.y = element_text(family = "my_font",
                                     size = 14,
                                     color = "black"),
         axis.text.x = element_text(family = "my_font",
                                    size = 10, 
                                    colour = "black",
                                    margin = margin(0, 0, -10, 0, unit = "pt")),
         axis.text.y = element_text(family = "my_font",
                                    size = 10,
                                    colour = "black"),
         legend.text = element_text(family = "my_font",
                                    size = 10,
                                    colour = "black"),
         legend.title = element_text(family = "my_font",
                                     size = 10,
                                     colour = "black"),
         panel.spacing = unit(20, "lines"),
         legend.position = "bottom",
         legend.key.height =  unit(0.01, "lines")) +
   labs(title = "Nesting Densities by Species",
        x = NULL, 
        y = "Nests Per Ha",
        color = "Grazing Intensity"))

ggsave(density.plotT,
       filename = "outputs/figs/AvianDensity_Treat.png",
       bg = "white",
       dpi = "print",
       height = 3,
       width = 5.58)

write_csv(birds.trtT, "working/Birds_Treatment_Density.csv")


# Running Kruskal-Wallace Tests ------------------------------------------------------------------------------


library(dunn.test)

dunn.test(birds.trtT[birds.trtT$Species == "BRBL", ]$estimate, 
          birds.trtT[birds.trtT$Species == "BRBL", ]$cTreat,
          method = "bonferroni")

dunn.test(birds.trtT[birds.trtT$Species == "WEME", ]$estimate, 
          birds.trtT[birds.trtT$Species == "WEME", ]$cTreat,
          method = "bonferroni")

# significant difference between rest and full
dunn.test(birds.trtT[birds.trtT$Species == "CCSP", ]$estimate, 
          birds.trtT[birds.trtT$Species == "CCSP", ]$cTreat,
          method = "bonferroni")

dunn.test(birds.trtT[birds.trtT$Species == "MODO", ]$estimate, 
          birds.trtT[birds.trtT$Species == "MODO", ]$cTreat,
          method = "bonferroni")

# significant difference between rest and heavy
dunn.test(birds.trtT[birds.trtT$Species == "RWBL", ]$estimate, 
          birds.trtT[birds.trtT$Species == "RWBL", ]$cTreat,
          method = "bonferroni")

dunn.test(birds.trtT[birds.trtT$Species == "YHBL", ]$estimate, 
          birds.trtT[birds.trtT$Species == "YHBL", ]$cTreat,
          method = "bonferroni")

# significant difference between moderate and heavy and rest and moderate
dunn.test(birds.trtT[birds.trtT$Species == "NOPI", ]$estimate, 
          birds.trtT[birds.trtT$Species == "NOPI", ]$cTreat,
          method = "bonferroni")

# significant difference between moderate and heavy
dunn.test(birds.trtT[birds.trtT$Species == "GADW", ]$estimate, 
          birds.trtT[birds.trtT$Species == "GADW", ]$cTreat,
          method = "bonferroni")

dunn.test(birds.trtT[birds.trtT$Species == "BWTE", ]$estimate, 
          birds.trtT[birds.trtT$Species == "BWTE", ]$cTreat,
          method = "bonferroni")

dunn.test(birds.trtT[birds.trtT$Species == "MALL", ]$estimate, 
          birds.trtT[birds.trtT$Species == "MALL", ]$cTreat,
          method = "bonferroni")

dunn.test(birds.trtT[birds.trtT$Species == "NSHO", ]$estimate, 
          birds.trtT[birds.trtT$Species == "NSHO", ]$cTreat,
          method = "bonferroni")


# Modeling differences in densities --------------------------------------------------------------------------


library(glmmTMB)
library(DHARMa)
library(emmeans)


hist(birds.trtT$estimate[birds.trtT$Species == "BRBL"])
hist(birds.trtT$estimate[birds.trtT$Species == "BRBL" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "BRBL" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "BRBL" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "BRBL" & birds.trtT$cTreat == "Heavy"])

aov.BRBL <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "BRBL"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.BRBL,
                                      n = 999,
                                      plot = T)

diagnose(aov.BRBL)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.BRBL)
emmeans(aov.BRBL, 
        pairwise ~ cTreat,
        adjust = "tukey")



hist(birds.trtT$estimate[birds.trtT$Species == "WEME"])
hist(birds.trtT$estimate[birds.trtT$Species == "WEME" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "WEME" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "WEME" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "WEME" & birds.trtT$cTreat == "Heavy"])

aov.WEME <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "WEME"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.WEME,
                                      n = 999,
                                      plot = T)

diagnose(aov.WEME)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.WEME)
emmeans(aov.WEME,
        pairwise~cTreat,
        adjust ="tukey")



hist(birds.trtT$estimate[birds.trtT$Species == "CCSP"])
hist(birds.trtT$estimate[birds.trtT$Species == "CCSP" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "CCSP" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "CCSP" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "CCSP" & birds.trtT$cTreat == "Heavy"])

aov.CCSP <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "CCSP"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.CCSP,
                                      n = 999,
                                      plot = T)

diagnose(aov.CCSP)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.CCSP)
emmeans(aov.CCSP, 
        pairwise ~ cTreat,
        adjust = "tukey")



hist(birds.trtT$estimate[birds.trtT$Species == "MODO"])
hist(birds.trtT$estimate[birds.trtT$Species == "MODO" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "MODO" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "MODO" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "MODO" & birds.trtT$cTreat == "Heavy"])

aov.MODO <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "MODO"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.MODO,
                                      n = 999,
                                      plot = T)

diagnose(aov.MODO)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.MODO)
emmeans(aov.MODO, 
        pairwise ~ cTreat,
        adjust = "tukey")



hist(birds.trtT$estimate[birds.trtT$Species == "RWBL"])
hist(birds.trtT$estimate[birds.trtT$Species == "RWBL" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "RWBL" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "RWBL" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "RWBL" & birds.trtT$cTreat == "Heavy"])

aov.RWBL <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "RWBL"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.RWBL,
                                      n = 999,
                                      plot = T)

diagnose(aov.RWBL)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.RWBL)
emmeans(aov.RWBL, 
        pairwise ~ cTreat,
        adjust = "tukey")



# hist(birds.trtT$estimate[birds.trtT$Species == "YHBL"])
# hist(birds.trtT$estimate[birds.trtT$Species == "YHBL" & birds.trtT$cTreat == "Rest"])
# hist(birds.trtT$estimate[birds.trtT$Species == "YHBL" & birds.trtT$cTreat == "Moderate"])
# hist(birds.trtT$estimate[birds.trtT$Species == "YHBL" & birds.trtT$cTreat == "Full"])
# hist(birds.trtT$estimate[birds.trtT$Species == "YHBL" & birds.trtT$cTreat == "Heavy"])
# 
# aov.YHBL <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
#                     data = filter(birds.trtT, Species == "YHBL"),
#                     ziformula = ~cTreat,
#                     family = ziGamma(link = "log"))
# 
# simulationOutput <- simulateResiduals(aov.YHBL,
#                                       n = 999,
#                                       plot = T)
# 
# diagnose(aov.YHBL)
# testResiduals(simulationOutput)
# testZeroInflation(simulationOutput)
# testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)
# 
# summary(aov.YHBL)
# emmeans(aov.YHBL, 
#         pairwise ~ cTreat,
#         adjust = "tukey")



hist(birds.trtT$estimate[birds.trtT$Species == "NOPI"])
hist(birds.trtT$estimate[birds.trtT$Species == "NOPI" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "NOPI" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "NOPI" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "NOPI" & birds.trtT$cTreat == "Heavy"])

aov.NOPI <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "NOPI"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))
 
simulationOutput <- simulateResiduals(aov.NOPI,
                                      n = 999,
                                      plot = T)

diagnose(aov.NOPI)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.NOPI)
emmeans(aov.NOPI, 
        pairwise ~ cTreat,
        adjust = "tukey")



hist(birds.trtT$estimate[birds.trtT$Species == "GADW"])
hist(birds.trtT$estimate[birds.trtT$Species == "GADW" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "GADW" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "GADW" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "GADW" & birds.trtT$cTreat == "Heavy"])

aov.GADW <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "GADW"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.GADW,
                                      n = 999,
                                      plot = T)

diagnose(aov.GADW)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.GADW)
emmeans(aov.GADW, 
        pairwise ~ cTreat)



hist(birds.trtT$estimate[birds.trtT$Species == "BWTE"])
hist(birds.trtT$estimate[birds.trtT$Species == "BWTE" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "BWTE" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "BWTE" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "BWTE" & birds.trtT$cTreat == "Heavy"])

aov.BWTE <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "BWTE"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.BWTE,
                                      n = 999,
                                      plot = T)

diagnose(aov.BWTE)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.BWTE)
emmeans(aov.BWTE, 
        pairwise ~ cTreat)


hist(birds.trtT$estimate[birds.trtT$Species == "MALL"])
hist(birds.trtT$estimate[birds.trtT$Species == "MALL" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "MALL" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "MALL" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "MALL" & birds.trtT$cTreat == "Heavy"])

aov.MALL <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "MALL"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.MALL,
                                      n = 999,
                                      plot = T)

diagnose(aov.MALL)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.MALL)
emmeans(aov.MALL, 
        pairwise ~ cTreat)


hist(birds.trtT$estimate[birds.trtT$Species == "NSHO"])
hist(birds.trtT$estimate[birds.trtT$Species == "NSHO" & birds.trtT$cTreat == "Rest"])
hist(birds.trtT$estimate[birds.trtT$Species == "NSHO" & birds.trtT$cTreat == "Moderate"])
hist(birds.trtT$estimate[birds.trtT$Species == "NSHO" & birds.trtT$cTreat == "Full"])
hist(birds.trtT$estimate[birds.trtT$Species == "NSHO" & birds.trtT$cTreat == "Heavy"])

aov.NSHO <- glmmTMB(estimate ~ cTreat + (1|Year/Replicate),
                    data = filter(birds.trtT, Species == "NSHO"),
                    ziformula = ~cTreat,
                    family = ziGamma(link = "log"))

simulationOutput <- simulateResiduals(aov.NSHO,
                                      n = 999,
                                      plot = T)

diagnose(aov.NSHO)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.NSHO)
emmeans(aov.NSHO, 
        pairwise ~ cTreat)
# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute "rm(list = ls(all = TRUE))" - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all = TRUE))
# Then, execute "cleanup(ask = FALSE)" to delete orphaned output
#  files from MARK. Execute "?cleanup" to learn more
cleanup(ask = FALSE)
