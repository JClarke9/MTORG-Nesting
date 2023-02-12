# Loading Libraries -------------------------------------------------------


library(tidyverse)


# Importing Data ----------------------------------------------------------


totals21 <- read.csv("working/totals21.csv", 
                     row.names = 1)                             # read in the data set

totals22 <- read.csv("working/totals22.csv", 
                     row.names = 1)                             # read in the data set


# Calculating species totals ----------------------------------------------


totals <- bind_rows(totals21, totals22)

totals <- totals |> 
  group_by(Spec, cTreat) |> 
  summarise(Abundance = sum(Abundance))

totalsW <- pivot_wider(totals,
                       names_from = cTreat,
                       values_from = Abundance,
                       values_fill = 0) |> 
  relocate(Heavy, .before = Full)


rest <- totalsW |> 
  filter(Rest > 0) |> 
  select(Spec, Rest)

moderate <- totalsW |> 
  filter(Moderate > 0) |> 
  select(Spec, Moderate)

full <- totalsW |> 
  filter(Full > 0) |> 
  select(Spec, Full)

heavy <- totalsW |> 
  filter(Heavy > 0) |> 
  select(Spec, Heavy)


write_csv(totalsW,
          "outputs/files/Totals.csv")
