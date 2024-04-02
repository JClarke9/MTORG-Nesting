
# Load Libraries ----------------------------------------------------------

library(lubridate)
library(vegan)
library(tidyverse)

# Data Import -------------------------------------------------------------

raw <- read_csv("raw/Nesting.csv")

totals <- raw |> 
  group_by(Year, id) |> 
  summarise()

totals <- totals |> 
  group_by(Year) |> 
  summarise(count = n())

# Data Wrangling ----------------------------------------------------------

birds <- select(raw, c(Spec, id, Year, Visit.Interval, 
                       Expos, Date, FirstFound, AgeFound, Fate, Fate2, Clutch))

# converting dates to julian day
birds$Date <- as.POSIXct(birds$Date,                       # read the Date column as a dat
                         format = "%m/%d/%Y") |>         # format of the existing dates
  yday()                                               # convert dates to julian days

birds$FirstFound <- as.POSIXct(birds$FirstFound,           # read the date found column as a date
                               format = "%m/%d/%Y") |>   # show the format of the existing dates
  yday()                                               # convert dates a nest was found to julian days

# removing unknown species and birds we found when they hatched because we couldn't estimate a start date.
birds <- filter(birds, Spec != "DUCK" & Spec != "UNKN" & Expos != 0 & Fate2 != "Unknown")

#RMark needs the date that a nest was last checked

LastChecked <- birds |>                                  # select the data frame
  group_by(id) |>                                      # group data by Nest.ID
  filter(Visit.Interval == max(Visit.Interval)) |>     # select the most recent visit interval
  reframe(Date)                                        # select only the two relevant columns

LastChecked$id[duplicated(LastChecked$id)]

birds <- left_join(birds,                                  # select data frame 1
                 LastChecked,                          # select data frame 2
                 by= "id",                             # join the two data frames by the shared "Nest.ID" column
                 keep=FALSE)                           # remove duplicate "Nest.ID" column

birds <- rename(birds, "DateChecked" = "Date.x")
birds <- rename(birds, "LastChecked" = "Date.y")

# This is just used to calculate the date a nest last was confirmed occupied
# its important only for failed nests because successful nests it's the same as 
# the last checked

LastPresent <- birds |>                                  # select the data frame
  group_by(id) |>                                      # group by the Nest.ID column
  filter(Visit.Interval == max(Visit.Interval)) |>     # filter out the rows to show ony the most recent checks
  summarize(id,                                        # keep Nest.ID
            Expos,                                     # Exposure days
            DateChecked)                               # and date 

# Calculate the actual date that a nest was last active
# by subtracting the last check from the last Exposure day
# I needed to join the two lists together so that 
# the new list includes the new LastPresent column

LastPresent$LastPresent <- LastPresent$DateChecked - LastPresent$Expos

birds <- left_join(birds,                                  # select data frame 1
                 LastPresent,                          # select data frame 2
                 by = "id",
                 keep=FALSE)                           # remove duplicate "Nest.ID" column

# I needed to rename the columns to match the input values for RMark

birds <- rename(birds, "DateChecked" = "DateChecked.x")    # rename the Date a nest was checked column
birds <- rename(birds, "Expos" = "Expos.x")

# I needed to create a table to use in order to select only
# successful nests because it was counting the last Exposure day
# as a separate nest for the failed nests since it went from 1 to 0

# I first filtered out the failed nest and used an anti-join because otherwise
# it counted failed nests as successful since they had visits where they were
# successful

f.sort <- birds |>                                       # select the data frame
  filter(Fate == "1") |>                               # filter out only failed nests
  group_by(id,                                         # group by Nest.ID
           Fate) |>                                    # group by Fate2
  summarise() |>                                       # remove unrelevant columns
  as.data.frame()                                      # coerce table into a data frame

f.sort$id[duplicated(f.sort$id)]

# this is to create a list of just the successful nests and 
# all their associated data

success <- anti_join(birds,                              # select data frame 1 to remove values from
                     f.sort,                           # select data frame 2 that values will be removed with
                     by = "id") |>                     # remove rows based on the shared Nest.ID column
  filter(Fate == 0) |>                                 # filter out successful (0) nests
  group_by(id,                                         # group by Nest.ID
           Spec,                                       # species
           LastPresent,                                # date last present
           LastChecked,                                # date last checked
           Fate,                                       # Fate
           FirstFound,                                 # date first found
           AgeFound) |>                                # age first found
  reframe(Expos = sum(Expos),                          # create a column with Expos added together for Nest IDs
          AgeFound) |>
  distinct(id, .keep_all = TRUE) |>                    # select only unique nest ID and keep all other columns
  as.data.frame()                                      # coerce table into a data frame

success$id[duplicated(success$id)]

# this is to create a list of just the failed nests and 
# all their associated data. I couldn't use the previous list
# because it only selects the last check so I can't sum
# the Exposure days. This method lets me remove the successul nests
# and sum the Exposure days that way.

failed <- anti_join(birds,                               # select data frame 1 to remove values from
                    success,                           # select data frame 2 that values will be removed with
                    by = "id") |>                      # remove rows based on the shared Nest.ID column
  filter(Fate == 1) |>                                 # filter out failed (1) nests
  group_by(id,                                         # group by Nest.ID
           Spec,                                       # species
           LastPresent,                                # date last present
           LastChecked,                                # date last checked
           Fate,                                       # Fate
           FirstFound,                                 # date first found
           AgeFound) |>                                # age first found
  reframe(Expos = sum(Expos),                          # create a column with Expos added together for Nest IDs
          AgeFound) |>
  as.data.frame()                                      # coerce table into a data frame

failed$id[duplicated(failed$id)]

nest <- bind_rows(success, failed)                     # combine the list of unique successful and failed nests

nest$id[duplicated(nest$id)]

# this is just to verify that I have the correct number
# of nests and there aren't any duplicates. Don't forget
# to include the number of unknown nests removed when
# checking the total number of nests.

check <- nest |>                                       # select the data frame
  group_by(id) |>                                      # group by Nest.ID to only show unique nests
  summarise()                                          # remove unrelevant columns

# RMark needs successful nest dates to have the date last present
# and date last checked to be equal

nest$LastChecked <- ifelse(nest$Fate == 0,             # if Fate2=0 (i.e. successful)
                           nest$LastPresent,           # set the day last checked equal to last day present
                           nest$LastChecked)           # if Fate2 < 1 don't change the day last checked

data <- select(birds, Year:id, Visit.Interval, DateChecked,
               InitBHCO:InitClutch, Open:pDoD) |>    # select the other data columns from the altered raw data 
  distinct(id, .keep_all = TRUE)                       # select only unique nest ID and keep all other columns

data$id[duplicated(data$id)]

nest <- left_join(nest,                                # select data frame 1
                  data,                                # select data frame 2
                  by = "id",                           # combine the data frames by the shared Nest.ID column
                  keep = FALSE)                        # remove duplicate "Nest.ID" column

nest <- full_join(nest, sched, by = c("Year", "cTreat"="Intensity"))

nest$Stage <- ifelse(nest$Nestling == 1, "Nestling",
                     ifelse(nest$Incubating == 1, "Incubating",
                            "Laying"))

nest <- relocate(nest,
                 Stage,
                 .after = DateChecked)

nest$AgeFound <- as.numeric(nest$AgeFound)

nest$Initiation <- ifelse(is.na(nest$AgeFound), nest$FirstFound,
                          nest$FirstFound - nest$AgeFound)

nest <- relocate(nest,
                 Initiation,
                 .before = FirstFound)

nest$BHCOPres <- ifelse(nest$BHCONum > 0, 1, 0)

# I added one to each of these otherwise a nest with the same day as the start date would have 0 grazing days
nest$grazed <- ifelse(nest$cTreat == "Rest", 0,
                      ifelse(nest$Initiation < nest$start & nest$LastChecked < nest$start, 0,
                             ifelse(nest$Initiation > nest$end & nest$LastChecked > nest$end, 0,
                                    ifelse(nest$Initiation < nest$start & nest$LastChecked <= nest$end, nest$LastChecked - nest$start + 1,
                                           ifelse(nest$Initiation >= nest$start & nest$LastChecked <= nest$end, abs(nest$LastChecked - nest$Initiation) + 1,
                                                  ifelse(nest$Initiation > nest$start & nest$LastChecked > nest$end, abs(nest$end - nest$Initiation) + 1,
                                                         ifelse(nest$Initiation <= nest$start & nest$LastChecked >= nest$end, nest$end - nest$start + 1,
                                                                NA)))))))

nest$grazep <- ifelse(nest$grazed > 0, 1, 0)

# remove columns with NA values in environmental covariates

MISSING <- is.na(nest$KBG) |
  is.na(nest$Smooth.Brome) |
  is.na(nest$Litter) |
  is.na(nest$Bare) |
  is.na(nest$Forb) |
  is.na(nest$Grasslike) |
  is.na(nest$Woody) |
  is.na(nest$Litter.Depth) |
  is.na(nest$Veg.Height) |
  is.na(nest$VOR) |
  is.na(nest$cTreat)

sum(MISSING)

nest <- subset(nest, 
               subset = !MISSING)

nest$cTreat <- as.factor(nest$cTreat)                  # coerce the treatment column into a factor
nest$pTreat <- as.factor(nest$pTreat)
nest$Replicate <- as.factor(nest$Replicate)
nest$Open <- as.factor(nest$Open)
nest$Stage <- as.factor(nest$Stage)
nest$Fate <- as.factor(nest$Fate)
nest$id <- as.factor(nest$id)
nest$Spec <- as.factor(nest$Spec)
nest$Bare <- as.numeric(nest$Bare)
nest$Litter.Depth <- as.numeric(nest$Litter.Depth)
nest$AgeFound <- as.numeric(nest$AgeFound)

str(nest)                                              # check the structure of the data

rm(list = ls()[!ls() %in%  "nest"])

# this loop with pull out each species and calculate relative vegetation cover, 
# frequencies for each encounter history, and 
# recreate the data frame.

spec.nest <- data.frame()

for (i in unique(nest$Spec)) {
  spec.surv <- filter(nest, Spec==i)
  
  spec.surv$TotalVegCover <- spec.surv |> 
    select(KBG:Woody) |> 
    rowSums(na.rm=TRUE)
  
  spec.surv <- spec.surv |> 
    mutate(across(KBG:Woody, ~ .x/TotalVegCover * 100))
  
  spec.surv$Litter.Depth <- as.integer(spec.surv$Litter.Depth)
  
  spec.surv$Veg.Height <- as.integer(spec.surv$Veg.Height)
  spec.surv$AgeFound <- as.numeric(spec.surv$AgeFound)
  
  spec.surv <- select(spec.surv, id:grazep)
  
  # I wanted to group nests into unique encounter histories
  # to create a frequency for each of those histories
  spec.unq <- spec.surv |> 
    group_by(FirstFound,
             LastPresent,
             LastChecked) |> 
    summarise(Freq=n()) |> 
    ungroup()
  
  # joining the two data frames to create a frequency column for the data frame
  spec.surv <- full_join(spec.surv,
                         spec.unq,
                         by = c("FirstFound",
                                "LastPresent",
                                "LastChecked"))
  
  spec.surv <- relocate(spec.surv,
                        Freq,
                        .after = Fate)
  
  #standardizing dates so each column starts at the first day a nest was found
  spec.surv$LastPresent <- spec.surv$LastPresent - min(spec.surv$FirstFound) + 1
  spec.surv$LastChecked <- spec.surv$LastChecked - min(spec.surv$FirstFound) + 1
  spec.surv$FirstFound <- spec.surv$FirstFound - min(spec.surv$FirstFound) + 1
  spec.surv$AgeDay1 <- spec.surv$AgeFound - spec.surv$FirstFound + 1
  
  spec.surv <- rename(spec.surv, LitterD=Litter.Depth)
  spec.surv <- rename(spec.surv, SmoothB=Smooth.Brome)
  
  
  spec.surv$cTreat <- dplyr::recode(spec.surv$cTreat,
                                    "Heavy" = "68",
                                    "Full" = "49",
                                    "Moderate" = "39",
                                    "Rest" = "0") |> 
    as.factor()
  
  spec.surv$pTreat <- dplyr::recode(spec.surv$pTreat,
                                    "Heavy" = "68",
                                    "Full" = "49",
                                    "Moderate" = "39",
                                    "Rest" = "0") |> 
    as.factor()
  
  spec.surv <- select(spec.surv, id:AgeDay1)
  
  spec.nest <- bind_rows(spec.nest, spec.surv)
}

spec.nest <- select(spec.nest, Year, id, Spec, FirstFound, LastPresent, 
                    LastChecked, Fate, Freq, Expos, AgeFound, AgeDay1, Stage, 
                    Laying, Incubating, Nestling, InitBHCO, BHCONum, BHCOPres, 
                    InitClutch, Clutch, Open, KBG:grazep, start, end)

#remove everything but the final dataframe
rm(list = ls()[!ls() %in%  c("spec.nest", "nest")])

write_csv(spec.nest, "working/RMarknesting.csv")
