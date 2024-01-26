# Loading Libraries -------------------------------------------------------


library(tidyverse)


# Importing Data ----------------------------------------------------------


totals21 <- read.csv("working/totals21.csv")                             # read in the data set

totals22 <- read.csv("working/totals22.csv")                             # read in the data set

totals23 <- read.csv("working/totals23.csv")


# Calculating species totals ----------------------------------------------


totals <- bind_rows(totals21, totals22, totals23)

groups <- totals |> 
  group_by(Group, Spec) |> 
  summarise()

groups <- groups |> 
  group_by(Group) |> 
  summarise(count = n())

totals <- totals |> 
  group_by(Spec, cTreat, Group) |> 
  summarise(Abundance = sum(Abundance))

totalsW <- pivot_wider(totals,
                       names_from = cTreat,
                       values_from = Abundance,
                       values_fill = 0) |> 
  relocate(Rest, .before = Moderate)

totalsW <- totalsW |> 
  mutate(Totals = rowSums(across(Rest:Heavy)))

totalsW <- bind_rows(totalsW,
                     colSums(totalsW[c("Rest", "Moderate",
                                       "Full", "Heavy", "Totals")]))

totalsW <- rename(totalsW,
                  "Species" = "Spec")

totalsW$Species <- factor(totalsW$Species,
                          levels = c("PBGR", "AMBI", "MALL", "GADW", "NOPI", 
                                     "AMWI", "NSHO", "BWTE", "GWTE", "LESC", 
                                     "STGR", "AMCO", "SORA", "KILL", "WILL", 
                                     "UPSA", "WISN", "WIPH", "MODO", "EAKI", 
                                     "MAWR", "DICK", "CCSP", "GRSP", "SAVS", 
                                     "CCLO", "WEME", "BOBO", "YHBL", "RWBL", 
                                     "BRBL", "COGR"))

totalsW$Species <- case_when(totalsW$Species == "AMBI" ~ "American Bittern",
                             totalsW$Species == "AMCO" ~ "American Coot",
                             totalsW$Species == "AMWI" ~ "American Wigeon",
                             totalsW$Species == "BOBO" ~ "Bobolink",
                             totalsW$Species == "BRBL" ~ "Brewer's Blackbird",
                             totalsW$Species == "BWTE" ~ "Blue-winged Teal",
                             totalsW$Species == "CCLO" ~ "Chestnut-collared Longspur",
                             totalsW$Species == "CCSP" ~ "Clay-colored Sparrow",
                             totalsW$Species == "COGR" ~ "Common Grackle",
                             totalsW$Species == "DICK" ~ "Dickcissel",
                             totalsW$Species == "EAKI" ~ "Eastern Kingbird",
                             totalsW$Species == "GADW" ~ "Gadwall",
                             totalsW$Species == "GRSP" ~ "Grasshopper Sparrow",
                             totalsW$Species == "GWTE" ~ "Green-winged Teal",
                             totalsW$Species == "KILL" ~ "Killdeer",
                             totalsW$Species == "LESC" ~ "Lesser Scaup",
                             totalsW$Species == "MALL" ~ "Mallard",
                             totalsW$Species == "MAWR" ~ "Marsh Wren",
                             totalsW$Species == "MODO" ~ "Mourning Dove",
                             totalsW$Species == "NOPI" ~ "Northern Pintail",
                             totalsW$Species == "NSHO" ~ "Northern Shoveler",
                             totalsW$Species == "PBGR" ~ "Pied-billed Grebe",
                             totalsW$Species == "RWBL" ~ "Red-winged Blackbird",
                             totalsW$Species == "SAVS" ~ "Savannah Sparrow",
                             totalsW$Species == "SORA" ~ "Sora",
                             totalsW$Species == "STGR" ~ "Sharp-tailed Grouse",
                             totalsW$Species == "UPSA" ~ "Upland Sandpiper",
                             totalsW$Species == "WEME" ~ "Western Meadowlark",
                             totalsW$Species == "WILL" ~ "Willet",
                             totalsW$Species == "WIPH" ~ "Wilson's Phalarope",
                             totalsW$Species == "WISN" ~ "Wilson's Snipe",
                             totalsW$Species == "YHBL" ~ "Yellow-headed Blackbird",
                             .default = "Totals")

totalsW$Species <- factor(totalsW$Species,
                          levels = c("Pied-billed Grebe", 
                                     "American Bittern", 
                                     "Mallard", 
                                     "Gadwall", 
                                     "Northern Pintail", 
                                     "American Wigeon", 
                                     "Northern Shoveler", 
                                     "Blue-winged Teal", 
                                     "Green-winged Teal", 
                                     "Lesser Scaup", 
                                     "Sharp-tailed Grouse", 
                                     "American Coot", 
                                     "Sora", 
                                     "Killdeer", 
                                     "Willet", 
                                     "Upland Sandpiper", 
                                     "Wilson's Snipe", 
                                     "Wilson's Phalarope", 
                                     "Mourning Dove", 
                                     "Eastern Kingbird", 
                                     "Marsh Wren", 
                                     "Dickcissel", 
                                     "Clay-colored Sparrow", 
                                     "Grasshopper Sparrow", 
                                     "Savannah Sparrow", 
                                     "Chestnut-collared Longspur", 
                                     "Western Meadowlark", 
                                     "Bobolink", 
                                     "Yellow-headed Blackbird", 
                                     "Red-winged Blackbird", 
                                     "Brewer's Blackbird", 
                                     "Common Grackle",
                                     "Totals"))

rest <- totalsW |> 
  filter(Rest > 0) |> 
  select(Species, Rest)

moderate <- totalsW |> 
  filter(Moderate > 0) |> 
  select(Species, Moderate)

full <- totalsW |> 
  filter(Full > 0) |> 
  select(Species, Full)

heavy <- totalsW |> 
  filter(Heavy > 0) |> 
  select(Species, Heavy)


write_csv(totalsW,
          "outputs/files/Totals.csv")
