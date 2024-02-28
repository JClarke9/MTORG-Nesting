# Reading in Libraries ----------------------------------------------------


library(tidyverse)


# Reading in Data ---------------------------------------------------------


jus_raw <- read_csv("raw/Nesting.csv")
cam_raw <- read_csv("raw/Cameron_Nests.csv")


# Data Wrangling ---------------------------------------------------------------------------------------------


just <- jus_raw |> 
  group_by(id, Spec) |> 
  filter(Treatment == "MTORG" &
           Spec != "UNK" &
           Spec != "UNKN" &
           Spec != "DUCK" &
           id != "EAKI04649545176423") |>
  summarise()

just <- just |> 
  group_by(Spec) |> 
  summarise(Count = n(),
            Research = "Justin") |> 
  rename("Alpha" = "Spec")

cam <- cam_raw |> 
  group_by(Species) |> 
  summarise(Count = sum(across(PBG20:SLG)),
            Research = "Cameron") |> 
  rename("Alpha" = "Species")

cgrec <- bind_rows(just, cam)

unique(cgrec$Alpha)

cgrec <- cgrec |> 
  mutate(Species = case_when(Alpha == "AMBI" ~ "American Bittern",
                             Alpha == "AMCO" ~ "American Coot",
                             Alpha == "AMWI" ~ "American Widgeon",
                             Alpha == "BOBO" ~ "Bobolink",
                             Alpha == "BRBL" ~ "Brewer's Blackbird",
                             Alpha == "BWTE" ~ "Blue-winged Teal",
                             Alpha == "CCLO" ~ "Chestnut-collared Longspur",
                             Alpha == "CCSP" ~ "Clay-colored Sparrow",
                             Alpha == "COGR" ~ "Common Grackle",
                             Alpha == "DICK" ~ "Dickcissel",
                             Alpha == "EAKI" ~ "Eastern Kingbird",
                             Alpha == "GADW" ~ "Gadwall",
                             Alpha == "GRSP" ~ "Grasshopper Sparrow",
                             Alpha == "GWTE" ~ "Green-winged Teal",
                             Alpha == "KILL" ~ "Killdeer",
                             Alpha == "LESC" ~ "Lesser Scaup",
                             Alpha == "MALL" ~ "Mallard",
                             Alpha == "MAWR" ~ "Marsh Wren",
                             Alpha == "MODO" ~ "Mourning Dove",
                             Alpha == "NOPI" ~ "Northern Pintail",
                             Alpha == "NSHO" ~ "Northern Shoveler",
                             Alpha == "PBGR" ~ "Pied-billed Grebe",
                             Alpha == "RWBL" ~ "Red-winged Blackbird",
                             Alpha == "SAVS" ~ "Savannah Sparrow",
                             Alpha == "SORA" ~ "Sora",
                             Alpha == "STGR" ~ "Sharp-tailed Grouse",
                             Alpha == "UPSA" ~ "Upland Sandpiper",
                             Alpha == "WEME" ~ "Western Meadowlark",
                             Alpha == "WILL" ~ "Willet",
                             Alpha == "WIPH" ~ "Wilson's Phalarope",
                             Alpha == "WISN" ~ "Wilson's Snipe",
                             Alpha == "YHBL" ~ "Yellow-headed Blackbird",
                             Alpha == "CAGO" ~ "Canada Goose",
                             Alpha == "CONI" ~ "Common Nighthawk",
                             Alpha == "HOLA" ~ "Horned Lark",
                             Alpha == "MAGO" ~ "Marbled Godwit"),
         .before = Alpha)

cgrec <- cgrec |>
    mutate(Group = case_when(
      Alpha == "AMBI" ~ "FAC",
      Alpha == "AMCO" ~ "WET",
      Alpha == "AMWI" ~ "FAC",
      Alpha == "BOBO" ~ "OBL",
      Alpha == "BRBL" ~ "FAC",
      Alpha == "BWTE" ~ "FAC",
      Alpha == "CCLO" ~ "OBL",
      Alpha == "CCSP" ~ "FAC",
      Alpha == "COGR" ~ "FAC",
      Alpha == "DICK" ~ "OBL",
      Alpha == "EAKI" ~ "FAC",
      Alpha == "GADW" ~ "FAC",
      Alpha == "GRSP" ~ "OBL",
      Alpha == "GWTE" ~ "FAC",
      Alpha == "KILL" ~ "FAC",
      Alpha == "LESC" ~ "FAC",
      Alpha == "MALL" ~ "FAC",
      Alpha == "MAWR" ~ "WET",
      Alpha == "MODO" ~ "FAC",
      Alpha == "NOPI" ~ "FAC",
      Alpha == "NSHO" ~ "FAC",
      Alpha == "PBGR" ~ "WET",
      Alpha == "RWBL" ~ "FAC",
      Alpha == "SAVS" ~ "OBL",
      Alpha == "SORA" ~ "WET",
      Alpha == "STGR" ~ "OBL",
      Alpha == "UPSA" ~ "OBL",
      Alpha == "WEME" ~ "OBL",
      Alpha == "WILL" ~ "FAC",
      Alpha == "WIPH" ~ "WET",
      Alpha == "WISN" ~ "WET",
      Alpha == "YHBL" ~ "WET",
      Alpha == "CAGO" ~ "FAC",
      Alpha == "CONI" ~ "FAC",
      Alpha == "HOLA" ~ "OBL",
      Alpha == "MAGO" ~ "OBL"),
    .after = Alpha)

cgrec_total <- cgrec |> 
  group_by(Species, Alpha, Group) |> 
  summarise(Abundance = sum(Count))

write_csv(cgrec_total,
          "CGREC_Bird_Communities/outputs/CGREC_totals.csv")
