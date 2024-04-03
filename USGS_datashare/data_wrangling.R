
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
                              ifelse(max(Clutch < InitClutch, InitClutch,
                                         NA))),
           max(Clutch),
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

birds <- birds |> 
  summarise(n_nest = n(),
            .by = c(species_abbreviation, site_abbreviation, year,
                    date_format, nid_estimated, nid_how, date_found,
                    age_found, fate, clutch_size))

usgs <- birds |> 
  select(species_abbreviation, site_abbreviation, year, date_format, nid_estimated, nid_how, clutch_size, n_nest)

usgs$comments <- ifelse(usgs$nid_how == "candling",
                        "Candled 2 eggs from the nest for nests found while incubating. For nests found at the nestling stage, we determined age based on feather/nestling development.",
                        ifelse(usgs$nid_how == "candling/back dating",
                               "Since eggs were often hard to see through, we would estimate the age by candling 2 eggs and back date based on hatch date to ensure accurate estimation.",
                               NA))

write_csv(usgs, "USGS_datashare/working/usgs.csv")
