
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


# Checking Missing Data/Mistakes ------------------------------------------


check_failed <- filter(raw, Fate == 1 & Fate2 == "Survive")
check_survive <- filter(raw, Fate == 0 & Fate2 != "Survive")


unknown <- filter(raw, Fate2 == "Unknown") |>          # sort out nests with an unknown Fate
  group_by(id, Spec) |>                                # group by Nest.ID and Species
  summarise()                                          # remove all unrelevant columns

zero.expo <- filter(raw, Expos == "0") |>              # sort out nests without Exposure days
  group_by(id, Spec) |>                                # group by Nest.ID and Species
  summarise()                                          # remove all unrelevant columns


# Data Wrangling ----------------------------------------------------------


# removing unknown species and birds we found when they hatched because we couldn't estimate a start date.

#RMark needs the date that a nest was last checked

# Calculate the actual date that a nest was last active
# by subtracting the last check from the last Exposure day
# I needed to join the two lists together so that 
# the new list includes the new LastPresent column

raw$AgeFound <- as.numeric(raw$AgeFound)

nest <- raw |> 
  filter(Treatment == "MTORG", Spec != "DUCK" & Spec != "UNKN" & Expos != 0 & Fate2 != "Unknown" & id != "YHBL04647295176591") |> 
  mutate(DateChecked = as.POSIXct(Date, format = "%m/%d/%Y") |> yday(),
         FirstFound = as.POSIXct(FirstFound, format = "%m/%d/%Y") |> yday(),
         LastChecked = DateChecked[which.max(Visit.Interval)],
         LastPresent = ifelse(Fate[which.max(Visit.Interval)] == 0, DateChecked[which.max(Visit.Interval)],
                              ifelse(Fate[which.max(Visit.Interval)] == 1, DateChecked[which.max(Visit.Interval)] - Expos[which.max(Visit.Interval)],
                                     NA)),
         Initiation = ifelse(is.na(AgeFound), FirstFound, FirstFound - AgeFound),
         BHCOPres = ifelse(BHCONum > 0, 1, 0),
         Fate = Fate,
         Nestling = ifelse(Stage == "Nestling", 1, 0),
         .by = id)

sched <- sched |> 
  mutate(start = as.POSIXct(sched$start, format = "%m/%d/%Y") |> yday(),
         end = as.POSIXct(sched$end, format = "%m/%d/%Y") |> yday())

head(nest$DateChecked)                                         # show the first few data frames
tail(nest$DateChecked)                                         # show the last few data frames


# Adding Treatment Data ---------------------------------------------------

nest$cTreat <- ifelse(nest$Paddock == 1 & nest$Year == 2021 | nest$Paddock == 5 & nest$Year == 2021 | nest$Paddock == 12 & nest$Year == 2021 | nest$Paddock == 16 & nest$Year == 2021, "Heavy",
                     ifelse(nest$Paddock == 2 & nest$Year == 2021 | nest$Paddock == 6 & nest$Year == 2021 | nest$Paddock == 9 & nest$Year == 2021 | nest$Paddock == 13 & nest$Year == 2021, "Full",
                            ifelse(nest$Paddock == 3 & nest$Year == 2021 | nest$Paddock == 7 & nest$Year == 2021 | nest$Paddock == 10 & nest$Year == 2021 | nest$Paddock == 14 & nest$Year == 2021, "Moderate",
                                   ifelse(nest$Paddock == 4 & nest$Year == 2021 | nest$Paddock == 8 & nest$Year == 2021 | nest$Paddock == 11 & nest$Year == 2021 | nest$Paddock == 15 & nest$Year == 2021, "Rest",
                                          ifelse(nest$Paddock == 1 & nest$Year == 2022 | nest$Paddock == 5 & nest$Year == 2022 | nest$Paddock == 12 & nest$Year == 2022 | nest$Paddock == 16 & nest$Year == 2022, "Rest",
                                                 ifelse(nest$Paddock == 2 & nest$Year == 2022 | nest$Paddock == 6 & nest$Year == 2022 | nest$Paddock == 9 & nest$Year == 2022 | nest$Paddock == 13 & nest$Year == 2022, "Heavy",
                                                        ifelse(nest$Paddock == 3 & nest$Year == 2022 | nest$Paddock == 7 & nest$Year == 2022 | nest$Paddock == 10 & nest$Year == 2022 | nest$Paddock == 14 & nest$Year == 2022, "Full",
                                                               ifelse(nest$Paddock == 4 & nest$Year == 2022 | nest$Paddock == 8 & nest$Year == 2022 | nest$Paddock == 11 & nest$Year == 2022 | nest$Paddock == 15 & nest$Year == 2022, "Moderate",
                                                                      ifelse(nest$Paddock == 1 & nest$Year == 2023 | nest$Paddock == 5 & nest$Year == 2023 | nest$Paddock == 12 & nest$Year == 2023 | nest$Paddock == 16 & nest$Year == 2023, "Moderate",
                                                                             ifelse(nest$Paddock == 2 & nest$Year == 2023 | nest$Paddock == 6 & nest$Year == 2023 | nest$Paddock == 9 & nest$Year == 2023 | nest$Paddock == 13 & nest$Year == 2023, "Rest",
                                                                                    ifelse(nest$Paddock == 3 & nest$Year == 2023 | nest$Paddock == 7 & nest$Year == 2023 | nest$Paddock == 10 & nest$Year == 2023 | nest$Paddock == 14 & nest$Year == 2023, "Heavy",
                                                                                           ifelse(nest$Paddock == 4 & nest$Year == 2023 | nest$Paddock == 8 & nest$Year == 2023 | nest$Paddock == 11 & nest$Year == 2023 | nest$Paddock == 15 & nest$Year == 2023, "Full",
                                                                                                  NA))))))))))))

nest$cDoD <- ifelse(nest$Year == 2021 & nest$cTreat == "Heavy", 76.3,
                   ifelse(nest$Year == 2021 & nest$cTreat == "Full", 52.3,
                          ifelse(nest$Year == 2021 & nest$cTreat == "Moderate", 51.7,
                                 ifelse(nest$Year == 2021 & nest$cTreat == "Rest", 0,
                                        ifelse(nest$Year == 2022 & nest$cTreat == "Heavy", 65.3,
                                               ifelse(nest$Year == 2022 & nest$cTreat == "Full", 45.3,
                                                      ifelse(nest$Year == 2022 & nest$cTreat == "Moderate", 35.2,
                                                             ifelse(nest$Year == 2022 & nest$cTreat == "Rest", 0,
                                                                    ifelse(nest$Year == 2023 & nest$cTreat == "Heavy", NA,
                                                                           ifelse(nest$Year == 2023 & nest$cTreat == "Full", NA,
                                                                                  ifelse(nest$Year == 2023 & nest$cTreat == "Moderate", NA,
                                                                                         ifelse(nest$Year == 2023 & nest$cTreat == "Rest", 0, 
                                                                                                NA))))))))))))

nest$pTreat <- ifelse(nest$cTreat == "Rest", "Moderate",
                     ifelse(nest$cTreat == "Moderate", "Full",
                            ifelse(nest$cTreat == "Full", "Heavy",
                                   ifelse(nest$cTreat == "Heavy", "Rest",
                                          NA))))

nest$pDoD <- ifelse(nest$Year == 2021 & nest$pTreat == "Heavy", 60.6,
                   ifelse(nest$Year == 2021 & nest$pTreat == "Full", 60.2,
                          ifelse(nest$Year == 2021 & nest$pTreat == "Moderate", 33.7,
                                 ifelse(nest$Year == 2021 & nest$pTreat == "Rest", 0,
                                        ifelse(nest$Year == 2022 & nest$pTreat == "Heavy", 76.3,
                                               ifelse(nest$Year == 2022 & nest$pTreat == "Full", 52.3,
                                                      ifelse(nest$Year == 2022 & nest$pTreat == "Moderate", 51.7,
                                                             ifelse(nest$Year == 2022 & nest$pTreat == "Rest", 0,
                                                                    ifelse(nest$Year == 2023 & nest$pTreat == "Heavy", 65.3,
                                                                           ifelse(nest$Year == 2023 & nest$pTreat == "Full", 45.3,
                                                                                  ifelse(nest$Year == 2023 & nest$pTreat == "Moderate", 35.2,
                                                                                         ifelse(nest$Year == 2023 & nest$pTreat == "Rest", 0, 
                                                                                                NA))))))))))))

nest <- full_join(nest, sched, by = c("Year", "cTreat"="Intensity"))

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


# Cleaning Veg Data -------------------------------------------------------


unique(nest$R1)
unique(nest$R2)
unique(nest$R3)
unique(nest$R4)

test <- filter(nest, R1 == "13+" | R1 == "18+" |
                 R2 == "18+" | R2 == "15+" | R2 == "13+" |
                 R3 == "18+" | R3 == "13+" | R4 == "18+" |
                 R4 == "13+")

nest$R1 <- recode(nest$R1,
                 "18+" = "20",
                 "13+" = "NA")

nest$R2 <- recode(nest$R2,
                 "18+" = "20",
                 "13+" = "NA",
                 "15+" = "NA")

nest$R3 <- recode(nest$R3,
                 "18+" = "20",
                 "13+" = "NA")

nest$R4 <- recode(nest$R4,
                 "18+" = "20",
                 "13+" = "NA")

test <- filter(nest, R1 == "NA" | R2 == "NA" | R3 == "NA" | R4 == "NA")

nest$R1 <- as.numeric(nest$R1)
nest$R2 <- as.numeric(nest$R2)
nest$R3 <- as.numeric(nest$R3)
nest$R4 <- as.numeric(nest$R4)

nest$VOR <- nest |> 
  select(R1:R4) |> 
  rowMeans()                                           # create a new column with the robel readings averaged

write_csv(nest, "working/community.csv")

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
nest$id <- as.factor(nest$id)
nest$Spec <- as.factor(nest$Spec)
nest$Bare <- as.numeric(nest$Bare)
nest$Litter.Depth <- as.numeric(nest$Litter.Depth)
nest$AgeFound <- as.numeric(nest$AgeFound)

str(nest)                                              # check the structure of the data


# Creating Stage Columns --------------------------------------------------

# Summarizing the data ----------------------------------------------------


nest <- nest |> 
  group_by(Year, id, Spec, FirstFound, LastPresent, 
           LastChecked, AgeFound, InitBHCO, InitClutch, Open, KBG,
           Smooth.Brome, Litter, Bare, Forb, Grasslike, Woody,
           Litter.Depth, Veg.Height, VOR, Paddock, Replicate, cTreat, 
           pTreat, cDoD, pDoD, grazed, grazep, start, end) |> 
  summarize(Expos = sum(Expos),
            Fate = max(Fate),
            Nestling = max(Nestling),
            BHCONum = max(BHCONum),
            BHCOPres = max(BHCOPres),
            Clutch = max(Clutch),
            Stage = Stage[which.max(Visit.Interval)]) |>
  ungroup()

raw$Fate <- as.factor(raw$Fate)                        # coerce survival (0-success, 1-fail) to a factor

str(nest)                                               # show the structure of the nest data


# Loop to create new columns ----------------------------------------------

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
  
  spec.nest <- bind_rows(spec.nest, spec.surv)
}

spec.nest <- spec.nest |> 
  select(Year, id, Spec, FirstFound, LastPresent, 
         LastChecked, Fate, Freq, Expos, AgeFound, AgeDay1, Stage, 
         Nestling, InitBHCO, BHCONum, BHCOPres, 
         InitClutch, Clutch, Open, KBG:grazep, start, end)

#remove everything but the final dataframe
rm(list = ls()[!ls() %in%  c("spec.nest", "nest")])

write_csv(spec.nest, "working/RMarknesting.csv")

