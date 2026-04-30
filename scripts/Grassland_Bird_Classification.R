

# Library and Data Import -------------------------------------------------


library(readxl)
library(tidyverse)

raw <- read_xlsx("raw/ACAD Global 2024.05.23.xlsx") |> 
  dplyr::select(`Common Name`, `Scientific Name`, USA, Canada, Mexico, `C America`, `Primary Breeding Habitat`, `Secondary Breeding Habitat`) 
  

# Data Wrangling ----------------------------------------------------------


grassland <- raw |> 
  filter(`Primary Breeding Habitat` %in% c("Grasslands:  Chihuahuan", "Grasslands:  Temperate", "Grasslands:  Tropical" ) |
           `Secondary Breeding Habitat` %in% c("Grasslands:  Chihuahuan", "Grasslands:  Temperate", "Grasslands:  Tropical", "Grasslands:  Pampas and Campos") &
           (USA == 1 | Canada == 1 | `C America` == 1))

mtorg <- raw |> 
  filter(`Common Name` %in% c("Brewer's Blackbird", "Red-winged Blackbird", "Yellow-headed Blackbird", 
                              "Western Meadowlark", "Mourning Dove", "Clay-colored Sparrow",
                              "Gadwall", "Northern Pintail", "Blue-winged Teal", "Northern Shoveler",
                              "Mallard"))
