
# Reading in Libraries ----------------------------------------------------


library(tidyverse)

SpeciesGrouping <- function(x) {
  x <- x |>                                                                     # select the data frame
    mutate(Group = case_when(                                                   # group together each species and
      endsWith(Spec, "STGR") ~ "OBL",                                           # list each species as FAC or OBL based on Vickery et al. 1999
      endsWith(Spec, "UPSA") ~ "OBL",
      endsWith(Spec, "SAVS") ~ "OBL",
      endsWith(Spec, "GRSP") ~ "OBL",
      endsWith(Spec, "CCLO") ~ "OBL",
      endsWith(Spec, "DICK") ~ "OBL",
      endsWith(Spec, "BOBO") ~ "OBL",
      endsWith(Spec, "WEME") ~ "OBL",
      endsWith(Spec, "GADW") ~ "FAC",
      endsWith(Spec, "AMWI") ~ "FAC",
      endsWith(Spec, "MALL") ~ "FAC",
      endsWith(Spec, "BWTE") ~ "FAC",
      endsWith(Spec, "NSHO") ~ "FAC",
      endsWith(Spec, "NOPI") ~ "FAC",
      endsWith(Spec, "GWTE") ~ "FAC",
      endsWith(Spec, "WILL") ~ "FAC",
      endsWith(Spec, "MODO") ~ "FAC",
      endsWith(Spec, "EAKI") ~ "FAC",
      endsWith(Spec, "CCSP") ~ "FAC",
      endsWith(Spec, "RWBL") ~ "FAC",
      endsWith(Spec, "BRBL") ~ "FAC",
      endsWith(Spec, "KILL") ~ "FAC",
      endsWith(Spec, "COGR") ~ "FAC",
      endsWith(Spec, "WIPH") ~ "WET", 
      endsWith(Spec, "WISN") ~ "WET", 
      endsWith(Spec, "YHBL") ~ "WET",
      endsWith(Spec, "PBGR") ~ "WET",
      endsWith(Spec, "SORA") ~ "WET",
      endsWith(Spec, "AMCO") ~ "WET",
      endsWith(Spec, "MAWR") ~ "WET",
      endsWith(Spec, "AMBI") ~ "FAC",
      endsWith(Spec, "LESC") ~ "FAC"
    ))
}


# Reading in Data ---------------------------------------------------------

raw <- read_csv("raw/Nesting.csv")

y23 <- filter(raw,
              Year == 2023 & 
                Treatment == "MTORG" &
                Spec != "UNK" &
                Spec != "UNKN" &
                Spec != "DUCK")

unique(y23$Spec)

NDreport23 <- y23 |> 
  group_by(id) |> 
  mutate(Exposure = sum(Expos)) |> 
  filter(Visit.Interval == max(Visit.Interval)) |> 
  select(Year, Spec, X, Y, FirstFound, AgeFound, Exposure, Fate, InitBHCO, 
         Clutch, InitClutch) |> 
  ungroup()

NDreport23 <- select(NDreport23, -id)

NDreport23$Fate <- ifelse(NDreport23$Fate != "Survive", "Failed", "Survive")

write.csv(NDreport23, "doc/NestDragging23.csv",
          row.names = F)

Fate_survive <- filter(NDreport23,
                       Fate == "Survive")

Fate_failed <- filter(NDreport23,
                      Fate != "Survive")

length(unique(Fate_survive$id))
length(unique(Fate_failed$id))

species <- y23 |> 
  distinct(id, Spec) |> 
  group_by(Spec) |> 
  summarise(count = n())

species <- SpeciesGrouping(species)

(group <- species |> 
    group_by(Group) |> 
    summarise(count = n()))
