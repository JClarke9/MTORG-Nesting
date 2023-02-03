library(lubridate)
library(vegan)
library(tidyverse)

raw <- read.csv("raw/Nesting.csv")                                              # read in the data set

raw$Date <- as.Date(raw$Date,                                                   # read the Date column as a date
                    "%m/%d/%y") |>                                              # show the format of the existing dates
  yday()                                                                        # calculate the Julian day

raw$aRobel <- rowMeans(raw[,28:31])                                             # create a new column with the robel readings averaged

# create a column for each data frame with the associated group information

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
      endsWith(Spec, "MAWR") ~ "WET"
    ))
}
raw <- SpeciesGrouping(raw) |> 
  drop_na(Group)


# I need to shorten the data to just show each species rather than each 
# observation (multiple of each species). I also need to assign each of those 
# individuals an abundance so that I can total them and determine how many 
# of each species is in each pasture/treatment


unique <- raw |>                                                                # select the data frame
  group_by(id,                                                                  # group data by Nest.ID in order to remove duplicate visits
           Spec,                                                                # species
           Pasture,                                                             # pasture
           cTreat,                                                              # current treatment
           Patch,
           Year,
           Group) |>                                                            # treatment
  summarize()                                                                   # summarize the data to only show the grouped variables

unique <- subset(unique,                                                        # this selects only unique observations (i.e. each individual).
                 Spec!="DUCK" & cTreat!= "")                                    # it also removes any columns with no id and/or no treatment info

# this is the totals for the 2021 data grouped by pasture, species, 
# and current treatment

totals21 <- unique |>                                                           # select the dataframe
  filter(Year=="2021") |>                                                       # filter the 2021 data out
  group_by(Pasture,                                                             # group by pasture
           Spec,                                                                # species
           cTreat,                                                              # current treatment
           Group) |>                                                            # species group (i.e. OBL, FAC, Gen)
  summarize(Abundance=n()) |>                                                   # count all abundances and create a new abundance column
  as.data.frame()                                                               # convert to a data frame

# this is the totals for the 2021 data grouped by pasture, species, 
# and current treatment

totals22 <- unique |>                                                           # select the dataframe
  filter(Year=="2022") |>                                                       # filter the 2021 data out
  group_by(Pasture,                                                             # group by pasture
           Spec,                                                                # species
           cTreat,                                                              # current treatment
           Group) |>                                                            # species group (i.e. OBL, FAC, Gen)
  summarize(Abundance=n()) |>                                                   # count all abundances and create a new abundance column
  as.data.frame()                                                               # convert to a data frame

patch <- unique |> 
  group_by(Patch,
           Spec,
           cTreat,
           Group) |> 
  summarize(Abundance=n()) |> 
  as.data.frame()

totals <- patch |> 
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
  relocate(Heavy,
           .before = Full)

totals.mat$SpecTotals <- rowSums(totals.mat)

totals.mat <- rbind(totals.mat, "Totals"= colSums(totals.mat))

totals.mat <- rownames_to_column(totals.mat, var = "Spec")

totals.mat <- SpeciesGrouping(totals.mat)

totals.mat <- column_to_rownames(totals.mat, var = "Spec")

write.csv(raw, "working/RAWadjusted.csv")
write.csv(totals.mat, "working/totals.csv")
write.csv(totals21, "working/totals21.csv")
write.csv(totals22, "working/totals22.csv")