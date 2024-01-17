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


# Calculate densities -----------------------------------------------------


# This loop calculates the densities and confidence intervals for all
# of the models listed above

birds.trt <- NULL  # Initialize the combined dataframe

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
          startsWith(Group, "2023") ~ "2023"))
    
    real_dataframe <- complete(real_dataframe,
                               Replicate,
                               Year,
                               fill = list(estimate = 0,
                                           se = 0,
                                           lcl = 0,
                                           ucl = 0))
    
    real_dataframe$Year <- factor(real_dataframe$Year,
                                  level = c("2021", "2022", "2023"))
    
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
          Year == 2021 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2022 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2023 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2021 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2022 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2023 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2021 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2022 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2023 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"]) * 64),
          estimate == 0 ~ 0),
        lcl = case_when(
          Year == 2021 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2022 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2023 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2021 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2022 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2023 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2021 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2022 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2023 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2021 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2022 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2023 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"])) / (lcl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"]) * 64),
          estimate == 0 ~ 0),
        ucl = case_when(
          Year == 2021 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2022 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2023 & Replicate == "NE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE"]) * 64),
          Year == 2021 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2022 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2023 & Replicate == "SE" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE"]) * 64),
          Year == 2021 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2022 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2023 & Replicate == "SW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW"]) * 64),
          Year == 2021 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2022 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW"]) * 64),
          Year == 2023 & Replicate == "NW" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"])) / (ucl^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW"]) * 64),
          estimate == 0 ~ 0),)
    
    # Combine with the existing birds.trt dataframe
    birds.trt <- rbind(birds.trt, group |> mutate(Species = species))
  } else {
    warning(paste0("Objects ", surv_obj, " and/or ", trt_obj, " not found in the environment. Skipping species ", species))
  }
}

(birdY.density <- birds.trt |> 
  group_by(Species) |> 
  summarize(birds_ha = mean(estimate),
            se =  mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl)))


# Plot densities ----------------------------------------------------------


(density.plotY <- ggplot(birds.trt, 
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
   labs(title = "Avian Densities",
        x = NULL, 
        y = "Nests Per Ha"))

ggsave(density.plotY,
       filename = "outputs/figs/AvianDensity_Year.png",
       dpi = "print",
       bg = NULL,
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


# Calculate densities -----------------------------------------------------


# This loop calculates the densities and confidence intervals for all
# of the models listed above

birds.trt <- NULL  # Initialize the combined dataframe

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
          grepl("2023", Group) ~ "2023"),
        Replicate = case_when(
          endsWith(Group, "NE") ~ "NE",
          endsWith(Group, "SE") ~ "SE",
          endsWith(Group, "SW") ~ "SW",
          endsWith(Group, "NW") ~ "NW"))
    
    real_dataframe$cTreat <- factor(real_dataframe$cTreat,
                                    levels = c("Rest", "Moderate", "Full", "Heavy"))
    
    real_dataframe$Year <- factor(real_dataframe$Year,
                                  level = c("2021", "2022", "2023"))
    
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
          Year == 2021 & Replicate == "NE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "Rest"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "Rest"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "Moderate"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "Moderate"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj$)cTreat == "Full"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "Full"]) * 64),
          Year == 2021 & Replicate == "NE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "Heavy"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NE" & get(surv_obj)$cTreat == "Heavy"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "Rest"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "Rest"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "Moderate"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "Moderate"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj$)cTreat == "Full"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "Full"]) * 64),
          Year == 2021 & Replicate == "SE" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "Heavy"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SE" & get(surv_obj)$cTreat == "Heavy"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "Rest"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "Rest"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "Moderate"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "Moderate"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj$)cTreat == "Full"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "Full"]) * 64),
          Year == 2021 & Replicate == "SW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "Heavy"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "SW" & get(surv_obj)$cTreat == "Heavy"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Rest" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "Rest"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "Rest"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Moderate" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "Moderate"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2022" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "Moderate"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Full" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj$)cTreat == "Full"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2023" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "Full"]) * 64),
          Year == 2021 & Replicate == "NW" & cTreat == "Heavy" & estimate != 0 ~ length(unique(get(surv_obj)$id[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "Heavy"])) / (estimate^mean(get(surv_obj)$AgeFound[get(surv_obj)$Year == "2021" & get(surv_obj)$Replicate == "NW" & get(surv_obj)$cTreat == "Heavy"]) * 64),
          estimate == 0 ~ 0),
        lcl = case_when(
          estimate == 0 ~ 0),
        ucl = case_when(
          estimate == 0 ~ 0),))
    
    # Combine with the existing birds.trt dataframe
    birds.trt <- rbind(birds.trt, group |> mutate(Species = species))
    } else {
      warning(paste0("Objects ", surv_obj, " and/or ", trt_obj, " not found in the environment. Skipping species ", species))
    }
}

birdT.density <- birds.trt |> 
  group_by(cTreat,
           Species) |> 
  summarize(birds_ha = mean(estimate),
            se =  mean(se),
            lcl = mean(lcl),
            ucl = mean(ucl))


# Plot densities ----------------------------------------------------------


(density.plotT <- ggplot(birdT.density, 
                         aes(x = cTreat, 
                             y = birds_ha)) +
   geom_point(aes(x = cTreat,
                  y = birds_ha),
              size = 5) +
   geom_errorbar(aes(x = cTreat,
                     ymin = lcl,
                     ymax = ucl),
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
                             colour = "black")) +                                  # change the color of the axis titles
   facet_wrap(~Species,
              nrow = 2,
              ncol = 4) +
   labs(title = "Avian Densities",
        x = NULL, 
        y = "Nests Per Ha"))


# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all = TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
#  files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)
