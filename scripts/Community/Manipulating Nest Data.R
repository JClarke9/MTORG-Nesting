library(lubridate)
library(vegan)
library(tidyverse)
source("scripts/Functions/Species_Grouping_Code.R")

raw <- read_csv("working/community.csv")                                        # read in the data set


# Data Wrangling ----------------------------------------------------------


raw <- SpeciesGrouping(raw) |> 
  drop_na(Group)

# I need to shorten the data to just show each species rather than each 
# observation (multiple of each species). I also need to assign each of those 
# individuals an abundance so that I can total them and determine how many 
# of each species is in each pasture/treatment


unique <- raw |>                                                                # select the data frame
  group_by(id,                                                                  # group data by Nest.ID in order to remove duplicate visits
           Spec,                                                                # species
           Paddock,                                                             # pasture
           cTreat,                                                              # current treatment
           Replicate,
           Year,
           Group) |>                                                            # treatment
  summarize()                                                                   # summarize the data to only show the grouped variables

unique <- subset(unique,                                                        # this selects only unique observations (i.e. each individual).
                 cTreat!= "")                                                   # it also removes any columns with no id and/or no treatme

# this is the totals for the 2021 data grouped by pasture, species, 
# and current treatment

totals21 <- unique |>                                                           # select the dataframe
  filter(Year=="2021") |>                                                       # filter the 2021 data out
  group_by(Paddock,                                                             # group by pasture
           Spec,                                                                # species
           cTreat,                                                              # current treatment
           Group) |>                                                            # species group (i.e. OBL, FAC, Gen)
  summarize(Abundance=n()) |>                                                   # count all abundances and create a new abundance column
  as.data.frame()                                                               # convert to a data frame

# this is the totals for the 2022 data grouped by pasture, species, 
# and current treatment

totals22 <- unique |>                                                           # select the dataframe
  filter(Year=="2022") |>                                                       # filter the 2021 data out
  group_by(Paddock,                                                             # group by pasture
           Spec,                                                                # species
           cTreat,                                                              # current treatment
           Group) |>                                                            # species group (i.e. OBL, FAC, Gen)
  summarize(Abundance=n()) |>                                                   # count all abundances and create a new abundance column
  as.data.frame()                                                               # convert to a data frame


# this is the totals for the 2023 data grouped by pasture, species, 
# and current treatment

totals23 <- unique |>                                                           # select the dataframe
  filter(Year=="2023") |>                                                       # filter the 2021 data out
  group_by(Paddock,                                                             # group by pasture
           Spec,                                                                # species
           cTreat,                                                              # current treatment
           Group) |>                                                            # species group (i.e. OBL, FAC, Gen)
  summarize(Abundance=n()) |>                                                   # count all abundances and create a new abundance column
  as.data.frame()                                                               # convert to a data frame


# This is to create a totals table for each grazing intensity -------------


paddock <- unique |> 
  group_by(Paddock,
           Spec,
           cTreat,
           Group) |> 
  summarize(Abundance=n()) |> 
  as.data.frame()

totals <- paddock |> 
  group_by(Spec,
           cTreat,
           Group) |> 
  summarize(Abundance=sum(Abundance)) |> 
  as.data.frame()

totals.mat <- pivot_wider(totals[c(1:2,4)],                                     # select the data frame to turn into a matrix
                          names_from = cTreat,                                  # select the column names
                          values_from = Abundance,                              # select the abundance values for each species
                          values_fill = list(Abundance = 0)) |>                 # fill all NA with 0
  column_to_rownames("Spec") |>
  relocate(Rest,
           .before = Moderate)

totals.mat$SpecTotals <- rowSums(totals.mat)

totals.mat <- rbind(totals.mat, "Totals"= colSums(totals.mat))

totals.mat <- rownames_to_column(totals.mat, var = "Spec")

totals.mat <- SpeciesGrouping(totals.mat)

totals.mat <- column_to_rownames(totals.mat, var = "Spec")


# Creating site by species matrices ---------------------------------------



birds21 <- pivot_wider(totals21[c(1:3,5)],                                      # select the data frame to turn into a matrix
                       names_from = Spec,                                       # select the column names
                       values_from = Abundance,                                 # select the abundance values for each species
                       values_fill = list(Abundance = 0)) |>                    # fill all NA with 0
  column_to_rownames("Paddock")                                                 # set column names as the pasture ID

birds22 <- pivot_wider(totals22[c(1:3,5)],                                      # select the data frame to turn into a matrix
                       names_from = Spec,                                       # select the column names
                       values_from = Abundance,                                 # select the abundance values for each species
                       values_fill = list(Abundance = 0)) |>                    # fill all NA with 0
  column_to_rownames("Paddock")                                                 # set column names as the pasture ID

birds23 <- pivot_wider(totals23[c(1:3,5)],                                      # select the data frame to turn into a matrix
                       names_from = Spec,                                       # select the column names
                       values_from = Abundance,                                 # select the abundance values for each species
                       values_fill = list(Abundance = 0)) |>                    # fill all NA with 0
  column_to_rownames("Paddock")                                                 # set column names as the pasture ID


# Saving all of the files -------------------------------------------------

write_csv(totals21, "working/totals21.csv")
write_csv(totals22, "working/totals22.csv")
write_csv(totals23, "working/totals23.csv")
write_csv(totals.mat, "working/totals.csv")
write_csv(birds21, file = "working/birds21.csv")
write_csv(birds22, file = "working/birds22.csv")
write_csv(birds23, file = "working/birds23.csv")