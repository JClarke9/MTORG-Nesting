# Load libraries ----------------------------------------------------------


library(tidyverse)
library(RMark)


# Read in Data ------------------------------------------------------------


nest <- read.csv("working/RMarknesting.csv")


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


# Run constant survival models --------------------------------------------


WEME.trt <- mark(WEME.surv, 
                 nocc = max(WEME.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1)))

BRBL.trt <- mark(BRBL.surv, 
                 nocc = max(BRBL.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1)))

CCSP.trt <- mark(CCSP.surv, 
                 nocc = max(CCSP.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1)))

MODO.trt <- mark(MODO.surv, 
                 nocc = max(MODO.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1)))

RWBL.trt <- mark(RWBL.surv, 
                 nocc = max(RWBL.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1)))

GADW.trt <- mark(GADW.surv, 
                 nocc = max(GADW.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1)))

BWTE.trt <- mark(BWTE.surv, 
                 nocc = max(BWTE.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1)))

NOPI.trt <- mark(NOPI.surv, 
                 nocc = max(NOPI.surv$LastChecked), 
                 model = "Nest",
                 adjust = FALSE,
                 delete = TRUE,
                 model.parameters = list(S = list(formula =  ~1)))


# Calculate densities -----------------------------------------------------


# This loop calculates the densities and confidence intervals for all
# of the models listed above

# List all objects in the environment
all_objects <- ls()

# Initialize an empty list to store result dataframes
result_dataframes <- list()

# Loop through all objects
for (object in all_objects) {
  
  # Check if the object is a dataframe ending with ".surv"
  if (grepl("\\.surv$", object) && is.data.frame(get(object))) {
    
    # Extract species name
    species <- gsub("\\.surv$", "", object)
    
    # Create a dataframe with ".real" data
    trt_object <- get(paste0(species, ".trt"))
    
    # Check if the trt_object and results$real exist
    if (!is.null(trt_object) && "results" %in% names(trt_object) && "real" %in% names(trt_object$results)) {
      real_dataframe <- as.data.frame(trt_object$results$real)
      
      # Summarize the ".real" dataframe
      summarized_real <- summarize(real_dataframe,
                                   estimate = mean(estimate),
                                   se = mean(se),
                                   lcl = mean(lcl),
                                   ucl = mean(ucl))
      
      # Create model dataframe
      model <- data.frame(
        birds_ha = length(unique(get(object)$id)) / (summarized_real$estimate^mean(get(object)$AgeFound) * 64),
        lcl = length(unique(get(object)$id)) / (summarized_real$lcl^mean(get(object)$AgeFound) * 64),
        ucl = length(unique(get(object)$id)) / (summarized_real$ucl^mean(get(object)$AgeFound) * 64),
        Species = species
      )
      
      # Add the result dataframe to the list
      result_dataframes[[paste0(species, "_result")]] <- model
    }
  }
}

# Combine all result dataframes into a single dataframe
(bird.density <- do.call(rbind, result_dataframes))


# Plot densities ----------------------------------------------------------


(density.plot <- ggplot(bird.density, 
                        aes(x = Species, 
                            y = birds_ha)) +
    geom_point(aes(x = Species,
                   y = birds_ha),
               size = 5) +
    geom_errorbar(aes(x = Species,
                      ymin = lcl,
                      ymax = ucl),
                  linewidth = 1) +
    ggtitle("Avian Densities") +
    theme(panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill = "white",                           # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill = "white",                            # make the outer background transparent
                                         colour = NA),                              # remove any other colors
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text = element_text(size = 20, colour = "black"),                    # color the axis text
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size = 20,                                              # change the size of the axis titles
                            colour = "black")) +                                  # change the color of the axis titles
    labs(x = NULL, y = "Nests Per Ha"))

ggsave(density.plot,
       filename = "outputs/figs/AvianDensity.png",
       dpi = "print",
       bg = "white",
       height = 6,
       width = 15)


# If you want to clean up the mark*.inp, .vcv, .res and .out
#  and .tmp files created by RMark in the working directory,
#  execute 'rm(list = ls(all = TRUE))' - see 2 lines below.
# NOTE: this will delete all objects in the R session.
rm(list = ls(all = TRUE))
# Then, execute 'cleanup(ask = FALSE)' to delete orphaned output
#  files from MARK. Execute '?cleanup' to learn more
cleanup(ask = FALSE)