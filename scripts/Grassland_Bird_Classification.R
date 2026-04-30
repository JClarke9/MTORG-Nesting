

# Library and Data Import -------------------------------------------------


library(readxl)
library(tidyverse)

raw <- read_xlsx("raw/ACAD Global 2024.05.23.xlsx")


# Data Wrangling ----------------------------------------------------------


grassland <- raw |> 
  filter(str_detect(`Primary Breeding Habitat`, "Grasslands:") & str_detect(`Secondary Breeding Habitat`, "Grasslands:"))

grassland <- raw |> 
  dplyr::select(`Common Name`, `Scientific Name`, USA, Canada, Mexico, `C America`, `Primary Breeding Habitat`, `Secondary Breeding Habitat`) |> 
  filter(`Primary Breeding Habitat` %in% c("Grasslands:  Chihuahuan", "Grasslands:  Temperate", "Grasslands:  Tropical" ) |
           `Secondary Breeding Habitat` %in% c("Grasslands:  Chihuahuan", "Grasslands:  Temperate", "Grasslands:  Tropical", "Grasslands:  Pampas and Campos") &
           (USA == 1 | Canada == 1 | `C America` == 1))
