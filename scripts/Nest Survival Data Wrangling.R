
# Load Libraries ----------------------------------------------------------

library(lubridate)
library(vegan)
library(tidyverse)

# Data Import -------------------------------------------------------------

raw <- read_csv("raw/Nesting.csv")
sched <- read_csv("raw/MTORGschedule.csv")

totals <- raw |> 
  group_by(Year, id) |> 
  summarise()

totals <- totals |> 
  group_by(Year) |> 
  summarise(count = n())

# Data Wrangling ----------------------------------------------------------

check_failed <- filter(raw, Fate == 1 & Fate2 == "Survive")
check_survive <- filter(raw, Fate == 0 & Fate2 != "Survive")

CCSP <- filter(raw, Spec == "CCSP")

raw <- filter(raw, Treatment == "MTORG")
raw <- filter(raw, Spec != "DUCK" & Spec != "UNKN")

sched$start <- as.POSIXct(sched$start,
                    format = "%m/%d/%Y") |> 
  yday()

sched$end <- as.POSIXct(sched$end,
                        format = "%m/%d/%Y") |>
  yday()

raw$Date <- as.POSIXct(raw$Date,                       # read the Date column as a dat
                       format = "%m/%d/%Y") |>         # format of the existing dates
  yday()                                               # convert dates to julian days

raw$FirstFound <- as.POSIXct(raw$FirstFound,           # read the date found column as a date
                             format = "%m/%d/%Y") |>   # show the format of the existing dates
  yday()                                               # convert dates a nest was found to julian days

# double check dates are entered correctly
filter(raw, FirstFound > Date)

head(raw$Date)                                         # show the first few data frames
tail(raw$Date)                                         # show the last few data frames

raw$cTreat <- ifelse(raw$Paddock == 1 & raw$Year == 2021 | raw$Paddock == 5 & raw$Year == 2021 | raw$Paddock == 12 & raw$Year == 2021 | raw$Paddock == 16 & raw$Year == 2021, "Heavy",
                     ifelse(raw$Paddock == 2 & raw$Year == 2021 | raw$Paddock == 6 & raw$Year == 2021 | raw$Paddock == 9 & raw$Year == 2021 | raw$Paddock == 13 & raw$Year == 2021, "Full",
                            ifelse(raw$Paddock == 3 & raw$Year == 2021 | raw$Paddock == 7 & raw$Year == 2021 | raw$Paddock == 10 & raw$Year == 2021 | raw$Paddock == 14 & raw$Year == 2021, "Moderate",
                                   ifelse(raw$Paddock == 4 & raw$Year == 2021 | raw$Paddock == 8 & raw$Year == 2021 | raw$Paddock == 11 & raw$Year == 2021 | raw$Paddock == 15 & raw$Year == 2021, "Rest",
                                          ifelse(raw$Paddock == 1 & raw$Year == 2022 | raw$Paddock == 5 & raw$Year == 2022 | raw$Paddock == 12 & raw$Year == 2022 | raw$Paddock == 16 & raw$Year == 2022, "Rest",
                                                 ifelse(raw$Paddock == 2 & raw$Year == 2022 | raw$Paddock == 6 & raw$Year == 2022 | raw$Paddock == 9 & raw$Year == 2022 | raw$Paddock == 13 & raw$Year == 2022, "Heavy",
                                                        ifelse(raw$Paddock == 3 & raw$Year == 2022 | raw$Paddock == 7 & raw$Year == 2022 | raw$Paddock == 10 & raw$Year == 2022 | raw$Paddock == 14 & raw$Year == 2022, "Full",
                                                               ifelse(raw$Paddock == 4 & raw$Year == 2022 | raw$Paddock == 8 & raw$Year == 2022 | raw$Paddock == 11 & raw$Year == 2022 | raw$Paddock == 15 & raw$Year == 2022, "Moderate",
                                                                      ifelse(raw$Paddock == 1 & raw$Year == 2023 | raw$Paddock == 5 & raw$Year == 2023 | raw$Paddock == 12 & raw$Year == 2023 | raw$Paddock == 16 & raw$Year == 2023, "Moderate",
                                                                             ifelse(raw$Paddock == 2 & raw$Year == 2023 | raw$Paddock == 6 & raw$Year == 2023 | raw$Paddock == 9 & raw$Year == 2023 | raw$Paddock == 13 & raw$Year == 2023, "Rest",
                                                                                    ifelse(raw$Paddock == 3 & raw$Year == 2023 | raw$Paddock == 7 & raw$Year == 2023 | raw$Paddock == 10 & raw$Year == 2023 | raw$Paddock == 14 & raw$Year == 2023, "Heavy",
                                                                                           ifelse(raw$Paddock == 4 & raw$Year == 2023 | raw$Paddock == 8 & raw$Year == 2023 | raw$Paddock == 11 & raw$Year == 2023 | raw$Paddock == 15 & raw$Year == 2023, "Full",
                                                                                                  NA))))))))))))

raw$pTreat <- ifelse(raw$cTreat == "Rest", "Moderate",
                     ifelse(raw$cTreat == "Moderate", "Full",
                            ifelse(raw$cTreat == "Full", "Heavy",
                                   ifelse(raw$cTreat == "Heavy", "Rest",
                                          NA))))

raw$Fate <- as.factor(raw$Fate)                        # coerce survival (0-success, 1-fail) to a factor

unique(raw$R1)
unique(raw$R2)
unique(raw$R3)
unique(raw$R4)

test <- filter(raw, R1 == "13+" | R1 == "18+" |
                 R2 == "18+" | R2 == "15+" | R2 == "13+" |
                 R3 == "18+" | R3 == "13+" | R4 == "18+" |
                 R4 == "13+")

raw$R1 <- recode(raw$R1,
                 "18+" = "20",
                 "13+" = "NA")

raw$R2 <- recode(raw$R2,
                 "18+" = "20",
                 "13+" = "NA",
                 "15+" = "NA")

raw$R3 <- recode(raw$R3,
                 "18+" = "20",
                 "13+" = "NA")

raw$R4 <- recode(raw$R4,
                 "18+" = "20",
                 "13+" = "NA")

test <- filter(raw, R1 == "NA" | R2 == "NA" | R3 == "NA" | R4 == "NA")

raw$R1 <- as.numeric(raw$R1)
raw$R2 <- as.numeric(raw$R2)
raw$R3 <- as.numeric(raw$R3)
raw$R4 <- as.numeric(raw$R4)

raw$VOR <- raw |> 
  select(R1:R4) |> 
  rowMeans()                                           # create a new column with the robel readings averaged

write_csv(raw, "working/community.csv")

unknown <- filter(raw, Fate2 == "Unknown") |>          # sort out nests with an unknown Fate
  group_by(id, Spec) |>                                # group by Nest.ID and Species
  summarise()                                          # remove all unrelevant columns

zero.expo <- filter(raw, Expos == "0") |>              # sort out nests without Exposure days
  group_by(id, Spec) |>                                # group by Nest.ID and Species
  summarise()                                          # remove all unrelevant columns

raw <- anti_join(raw,                                  # select data frame 1 to remove rows from
                 unknown,                              # select data frame 2 to remove rows with
                 by="id")                              # remove data based on the shared Nest.ID column

raw <- anti_join(raw,                                  # select data frame 1 to remove rows from
                 zero.expo,                            # select data frame 2 to remove rows with
                 by="id")                              # remove data based on the shared Nest.ID column

#RMark needs the date that a nest was last checked

LastChecked <- raw |>                                  # select the data frame
  group_by(id) |>                                      # group data by Nest.ID
  filter(Visit.Interval == max(Visit.Interval)) |>     # select the most recent visit interval
  reframe(Date)                                        # select only the two relevant columns

LastChecked$id[duplicated(LastChecked$id)]

raw <- left_join(raw,                                  # select data frame 1
                 LastChecked,                          # select data frame 2
                 by= "id",                             # join the two data frames by the shared "Nest.ID" column
                 keep=FALSE)                           # remove duplicate "Nest.ID" column

raw <- rename(raw, "DateChecked" = "Date.x")
raw <- rename(raw, "LastChecked" = "Date.y")

raw <- select(raw, Year:DateChecked, LastChecked,
              Visit.Interval:Veg.Height, VOR, Paddock:Replicate, 
              cTreat:pTreat)

# This is just used to calculate the date a nest last was confirmed occupied
# its important only for failed nests because successful nests it's the same as 
# the last checked

LastPresent <- raw |>                                  # select the data frame
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

raw <- left_join(raw,                                  # select data frame 1
                 LastPresent,                          # select data frame 2
                 by = "id",
                 keep=FALSE)                           # remove duplicate "Nest.ID" column

# I needed to rename the columns to match the input values for RMark

raw <- rename(raw, "DateChecked" = "DateChecked.x")    # rename the Date a nest was checked column
raw <- rename(raw, "Expos" = "Expos.x")

raw <- select(raw, Year:LastChecked, LastPresent,
              Visit.Interval:Veg.Height, VOR, Paddock:Replicate, 
              cTreat:pTreat)

raw$Nestling <- ifelse(raw$Stage == "Nestling", 
                       1, 0)                           # 1 is nestling, 0 is not

raw$Incubating <- ifelse(raw$Stage == "Nestling", 1,
                         ifelse(raw$Stage == "Incubating", 1,
                                0))                    # 1 is incubating, 0 is not

raw$Laying <- ifelse(raw$Stage == "Nestling", 1, 
                     ifelse(raw$Stage == "Incubating", 1,
                            ifelse(raw$Stage == "Laying", 1,
                                   0)))                # 1 is laying, 0 is not

check_lay <- filter(raw, Stage == "Laying")
check_inc <- filter(raw, Stage == "Incubating")
check_nst <- filter(raw, Stage == "Nestling")

str(raw)                                               # show the structure of the raw data

# I needed to create a table to use in order to select only
# successful nests because it was counting the last Exposure day
# as a separate nest for the failed nests since it went from 1 to 0

# I first filtered out the failed nest and used an anti-join because otherwise
# it counted failed nests as successful since they had visits where they were
# successful

f.sort <- raw |>                                       # select the data frame
  filter(Fate == "1") |>                               # filter out only failed nests
  group_by(id,                                         # group by Nest.ID
           Fate) |>                                    # group by Fate2
  summarise() |>                                       # remove unrelevant columns
  as.data.frame()                                      # coerce table into a data frame

f.sort$id[duplicated(f.sort$id)]

# this is to create a list of just the successful nests and 
# all their associated data

success <- anti_join(raw,                              # select data frame 1 to remove values from
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
          AgeFound,                                    # keep the Age Found column
          BHCONum = max(BHCONum),
          Laying = max(Laying),
          Incubating = max(Incubating),
          Nestling = max(Nestling)) |>
  distinct(id, .keep_all = TRUE) |>                    # select only unique nest ID and keep all other columns
  as.data.frame()                                      # coerce table into a data frame

success$id[duplicated(success$id)]

# this is to create a list of just the failed nests and 
# all their associated data. I couldn't use the previous list
# because it only selects the last check so I can't sum
# the Exposure days. This method lets me remove the successul nests
# and sum the Exposure days that way.

failed <- anti_join(raw,                               # select data frame 1 to remove values from
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
          AgeFound,                                    # keep the Age Found column
          BHCONum = max(BHCONum),
          Laying = max(Laying),
          Incubating = max(Incubating),
          Nestling = max(Nestling)) |>
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

data <- select(raw, Year:id, Visit.Interval, DateChecked,
               InitBHCO:InitClutch, Open:pTreat) |>    # select the other data columns from the altered raw data 
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
