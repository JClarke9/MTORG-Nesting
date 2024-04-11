
# Load Libraries ----------------------------------------------------------

library(lubridate)
library(vegan)
library(tidyverse)


# Data Import -------------------------------------------------------------


raw <- read_csv("raw/Nesting.csv")

birds <- select(raw, c(Spec, id, Year, Visit.Interval, Expos, Date, FirstFound, 
                       AgeFound, Fate, Fate2, InitClutch, Clutch, Treatment))


# Data Wrangling ----------------------------------------------------------


birds$AgeFound <- as.integer(birds$AgeFound)


birds <- birds |> 
  filter(Spec != "DUCK" & Spec != "UNKN" & Expos != 0 & Fate2 != "Unknown" & !is.na(AgeFound)) |> 
  mutate(site_abbreviation = paste0("CGREC - ", Treatment),
         date_format = "ordinal (1 = January 1)",
         Date = as.POSIXct(Date, format = "%m/%d/%Y") |> yday(),
         FirstFound = as.POSIXct(FirstFound, format = "%m/%d/%Y") |> yday(),
         LastChecked = Date[which.max(Visit.Interval)],
         nid_estimated = FirstFound - AgeFound,
         nid_how = ifelse(Spec %in% c("KILL", "STGR", "WIPH", "UPSA", "WISN",
                                      "WILL", "AMCO", "PBGR", "SORA", "AMBI"),
                          "candling/back dating",
                          "candling"),
         fate = Fate2[which.max(Visit.Interval)],
         clutch_size = ifelse(max(Clutch) >= InitClutch, max(Clutch),
                              ifelse(max(Clutch) < InitClutch, InitClutch,
                                         NA)),
         .by = id)

birds <- birds |> 
  rename("species_abbreviation" = "Spec",
         "year" = "Year",
         "date_found" = "FirstFound",
         "age_found" = "AgeFound")

birds <- birds |> 
  select(id, species_abbreviation, site_abbreviation, year,
         date_format, nid_estimated, nid_how, date_found,
         age_found, fate, clutch_size) |> 
  distinct()

usgs <- birds |> 
  filter(clutch_size > 0) |> 
  summarise(n_nest = n(),
            .by = c(species_abbreviation, site_abbreviation, year,
                    date_format, nid_estimated, nid_how, date_found,
                    age_found, fate, clutch_size))

# usgs$comments <- ifelse(usgs$age_found < 0 & usgs$nid_how == "candling", "Negative ages represent nests in either the nest building or laying stage with an age 0 being the first day of incubation. For nests found in the incubating stage, we candled 2 eggs from the nest to determine age. For nests found in the nestling stage, we determined age based on nestling/feather development. Nest initiation date was considered the day a nest contained a full clutch of eggs and incubation began.",
#                         ifelse(usgs$age_found >= 0 & usgs$nid_how == "candling", "For nests found in the incubating stage, we candled 2 eggs from the nest to determine age. For nests found in the nestling stage, we determined age based on nestling/feather development. Nest initiation date was considered the day a nest contained a full clutch of eggs and incubation began.",
#                                ifelse(usgs$age_found < 0 & usgs$nid_how == "candling/back dating", "Negative ages represent nests in either the nest building or laying stage with an age 0 being the first day of incubation. Since eggs were often hard to see through, we would estimate the age by candling 2 eggs and back date based on hatch date to ensure accurate estimation. Nest initiation date was considered the day a nest contained a full clutch of eggs and incubation began.",
#                                       ifelse(usgs$age_found >= 0 & usgs$nid_how == "candling/back dating", "Since eggs were often hard to see through, we would estimate the age by candling 2 eggs and back date based on hatch date to ensure accurate estimation. Nest initiation date was considered the day a nest contained a full clutch of eggs and incubation began.",
#                                              NA))))

usgs$comments <- ifelse(usgs$nid_how == "candling", "For nests found in the incubating stage, we candled 2 eggs from the nest to determine age. For nests found in the nestling stage, we determined age based on nestling/feather development. Nest initiation date was considered the day a nest contained a full clutch of eggs and incubation began.",
                        ifelse(usgs$nid_how == "candling/back dating", "Since eggs were often hard to see through, we would estimate the age by candling 2 eggs and back date based on hatch date to ensure accurate estimation. Nest initiation date was considered the day a nest contained a full clutch of eggs and incubation began.",
                               NA))

usgs <- usgs |> 
  select(species_abbreviation, site_abbreviation, year, date_format, nid_estimated, nid_how, clutch_size, n_nest, comments)

write_csv(usgs, "USGS_datashare/working/usgs.csv")
