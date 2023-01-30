######################################## WEME Veg Wrangling #########################################################

library(lubridate)
library(vegan)
library(tidyverse)

######################################### Data import ###############################################################

raw <- read.csv("~/Git/NDSU/RMARK/Raw Data/Nesting.csv")

########################################### Data Wrangling ##########################################################

#This will separate species out into separate dataframes and remove the nest dataframe from the environment
#for(i in nest$Spec) {
#  assign(i, nest |> 
#           filter(Spec == i),
#         envir = .GlobalEnv)
#}

#rm(nest)

raw$Date <- as.Date(raw$Date,                                                   # read the Date column as a date
                    "%m/%d/%y") |>                                              # show the format of the existing dates
  yday()                                                                        # convert dates checked to julian days

raw$FirstFound <- as.Date(raw$FirstFound,                                       # read the date found column as a date
                          "%m/%d/%y") |>                                            # show the format of the existing dates
  yday()                                                                        # convert dates a nest was found to julian days

raw$Survive <- as.factor(raw$Survive)                                           # coerce survival (0-success, 1-fail) to a factor

head(raw$Date)                                                                  # show the first few data frames
tail(raw$Date)                                                                  # show the last few data frames

raw$VOR <- rowMeans(raw[,28:31])                                                # create a new column with the robel readings averaged

unknown <- filter(raw, Survive=="Unknown") |>                                   # sort out nests with an unknown Survive
  group_by(id, Spec) |>                                                         # group by Nest.ID and Species
  summarize()                                                                   # remove all unrelevant columns

zero.expo <- filter(raw, Expos=="0") |>                                         # sort out nests without Exposure days
  group_by(id, Spec) |>                                                         # group by Nest.ID and Species
  summarize()                                                                   # remove all unrelevant columns

raw <- anti_join(raw,                                                           # select data frame 1 to remove rows from
                 unknown,                                                       # select data frame 2 to remove rows with
                 by="id")                                                       # remove data based on the shared Nest.ID column

raw <- anti_join(raw,                                                           # select data frame 1 to remove rows from
                 zero.expo,                                                     # select data frame 2 to remove rows with
                 by="id")                                                       # remove data based on the shared Nest.ID column

#RMark needs the date that a nest was last checked

LastChecked <- raw |>                                                           # select the data frame
  group_by(id) |>                                                               # group data by Nest.ID
  filter(Visit.Interval == max(Visit.Interval)) |>                              # select the most recent visit interval
  summarize(id, Date)                                                           # select only the two relevant columns

raw <- left_join(raw,                                                           # select data frame 1
                 LastChecked,                                                   # select data frame 2
                 by="id",                                                       # join the two data frames by the shared "Nest.ID" column
                 keep=FALSE)                                                    # remove duplicate "Nest.ID" column

raw <- raw[c(1:37,39:40)]                                                       # getting ride of the notes columns

names(raw)[6] <- "DateChecked"                                                  # rename the Date a nest was checked column
names(raw)[39] <- "LastChecked"                                                 # rename the Last Checked column

# This is just used to calculate the date a nest was last occupied

LastPresent <- raw |>                                                           # select the data frame
  group_by(id) |>                                                               # group by the Nest.ID column
  filter(Visit.Interval == max(Visit.Interval)) |>                              # filter out the rows to show ony the most recent checks
  summarize(id,                                                                 # keep Nest.ID
            Expos,                                                              # Exposure days
            DateChecked)                                                        # and date 

# Calculate the actual date that a nest was last active
# by subtracting the last check from the last Exposure day
# I needed to join the two lists together so that 
# the new list includes the new LastPresent column

LastPresent$LastPresent <- LastPresent$DateChecked - LastPresent$Expos          

raw <- left_join(raw,                                                           # select data frame 1
                 LastPresent,                                                   # select data frame 2
                 by="id",                                                       # join by the shared Nest.ID column
                 keep=FALSE)                                                    # remove duplicate "Nest.ID" column

names(raw)[6] <- "DateChecked"                                                  # rename the Date a nest was checked column
names(raw)[8] <- "Expos"                                                        # rename the Exposure column

# I needed to rename the columns to match the input values for RMark

raw <- raw[c(1:39, 42)]                                                         # remove extra columns

str(raw)                                                                        # show the structure of the raw data


# I needed to create a table to use in order to select only
# successful nests because it was counting the last Exposure day
# as a separate nest for the failed nests since it went from 1 to 0

f.sort <- raw |>                                                                # select the data frame
  filter(Fate=="1") |>                                                          # filter out only failed nests
  group_by(id,                                                                  # group by Nest.ID
           Fate) |>                                                             # group by fate
  summarize() |>                                                                # remove unrelevant columns
  as.data.frame()                                                               # coerce table into a data frame

# this is to create a list of just the successful nests and 
# all their associated data

success <- anti_join(raw,                                                       # select data frame 1 to remove values from
                     f.sort,                                                    # select data frame 2 that values will be removed with
                     by="id") |>                                                # remove rows based on the shared Nest.ID column
  filter(Fate==0) |>                                                            # filter out successful (1) nests
  group_by(id,                                                                  # group by Nest.ID
           Spec,                                                                # species
           LastPresent,                                                         # date last present
           LastChecked,                                                         # date last checked
           Fate,                                                                # fate
           FirstFound,                                                          # date first found
           AgeFound) |>                                                         # age first found
  summarize(Expos=sum(Expos),                                                   # create a column with Expos added together for Nest IDs
            AgeFound,                                                           # keep the Age Found column
            BHCONum=max(BHCONum)) |>                                            # keep the Age Found column
  distinct(id, .keep_all=TRUE) |>                                               # select only unique nest ID and keep all other columns
  as.data.frame()                                                               # coerce table into a data frame

# this is to create a list of just the failed nests and 
# all their associated data. I couldn't use the previous list
# because it only selects the last check so I can't sum
# the Exposure days. This method lets me remove the successul nests
# and sum the Exposure days that way.

failed <- anti_join(raw,                                                        # select data frame 1 to remove values from
                    success,                                                    # select data frame 2 that values will be removed with
                    by="id") |>                                                 # remove rows based on the shared Nest.ID column
  filter(Fate==1) |>                                                            # filter out failed (0) nests
  group_by(id,                                                                  # group by Nest.ID
           Spec,                                                                # species
           LastPresent,                                                         # date last present
           LastChecked,                                                         # date last checked
           Fate,                                                                # fate
           FirstFound,                                                          # date first found
           AgeFound) |>                                                         # age first found
  summarize(Expos=sum(Expos),                                                   # create a column with Expos added together for Nest IDs
            AgeFound,
            BHCONum = max(BHCONum)) |>
  as.data.frame()                                                               # coerce table into a data frame

nest <- rbind(success, failed)                                                  # combine the list of unique successful and failed nests

# this is just to verify that I have the correct number
# of nests and there aren't any duplicates. Don't forget
# to include the number of unknown nests removed when
# checking the total number of nests.

check <- nest |>                                                                # select the data frame
  group_by(id) |>                                                               # group by Nest.ID to only show unique nests
  summarize()                                                                   # remove unrelevant columns

# RMark needs successful nest dates to have the date last present
# and date last checked to be equal

nest$LastChecked <- ifelse(nest$Fate==0,                                        # if fate=0 (i.e. successful)                       
                           nest$LastPresent,                                    # set the day last checked equal to last day present
                           nest$LastChecked)                                    # if fate < 1 don't change the day last checked

data <- raw[c(1:2, 12, 14:17, 19:27, 33:34, 36:38)] |>                                     # select the veg data columns from the altered raw data
  distinct(id, .keep_all=TRUE)                                                  # select only unique nest ID and keep all other columns

nest <- left_join(nest,                                                         # select data frame 1
                  data,                                                         # select data frame 2
                  by="id",                                                      # combine the data frames by the shared Nest.ID column
                  keep=FALSE)                                                   # remove duplicate "Nest.ID" column

nest <- relocate(nest,                                                          # select the data frame
                 FirstFound,                                                    # select the column to move
                 .before = LastPresent)                                         # select where to move the column to

nest <- nest |> 
  na.omit(nest$cTreat) |> 
  ungroup()

nest$cTreat <- as.factor(nest$cTreat)                                           # coerce the treatment column into a factor
nest$pTreat <- as.factor(nest$pTreat)

nest$Spec <- as.factor(nest$Spec)

str(nest)                                                                       # check the structure of the data

rm(list = ls()[!ls() %in%  "nest"])

# this loop with pull out each species and calculate relative vegetation cover, 
# frequencies for each encounter history, and 
# recreate the data frame.

spec.nest <- data.frame()

for (i in unique(nest$Spec)) {
  spec.surv <- filter(nest, Spec==i)
  
  spec.surv$TotalVegCover <- rowSums(spec.surv[16:22], na.rm=TRUE)
  spec.surv[16:22] <- spec.surv[16:22]/spec.surv$TotalVegCover *100
  
  spec.surv$Litter.Depth <- as.integer(spec.surv$Litter.Depth)
  
  spec.surv$Veg.Height <- as.integer(spec.surv$Veg.Height)
  spec.surv$AgeFound <- as.numeric(spec.surv$AgeFound)
  
  spec.surv <- spec.surv[c(1:29)]
  # I wanted to group nests into unique encounter histories
  # to create a frequency for each of those histories
  spec.unq <- spec.surv |> 
    group_by(FirstFound,
             LastPresent,
             LastChecked) |> 
    summarize(Freq=n()) |> 
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
  
  spec.surv$Found <- spec.surv$FirstFound
  
  #standardizing dates so each column starts at the first day a nest was found
  spec.surv$LastPresent <- spec.surv$LastPresent - min(spec.surv$FirstFound) + 1
  spec.surv$LastChecked <- spec.surv$LastChecked - min(spec.surv$FirstFound) + 1
  spec.surv$FirstFound <- spec.surv$FirstFound - min(spec.surv$FirstFound) + 1
  spec.surv$AgeDay1 <- spec.surv$AgeFound - spec.surv$FirstFound + 1
  
  
  spec.surv <- rename(spec.surv, LitterD=Litter.Depth)
  spec.surv <- rename(spec.surv, SmoothB=Smooth.Brome)
  
  spec.surv <- na.omit(spec.surv) |> 
    ungroup()
  
  spec.surv$cTreat <- dplyr::recode(spec.surv$cTreat,
                                    "Heavy" = "68",
                                    "Full" = "49",
                                    "Moderate" = "39",
                                    "Rest" = "0") |> 
    as.factor()
  
  spec.surv$start <- ifelse(spec.surv$Year =="2021",
                            "5/1/2021",
                            "5/1/2022") |> 
    as.Date("%m/%d/%y") |> 
    yday()
  
  #standardizing the julian dates
  spec.surv$Julian <- spec.surv$Found - spec.surv$AgeFound - spec.surv$start
  spec.surv$Julian <- spec.surv$Julian - min(spec.surv$Julian) + 1
  
  spec.surv <- spec.surv[c(1:30,32,34)]
  
  #filter21 <- filter(spec.surv, Year == "2021")
  #filter22 <- filter(spec.surv, Year == "2022")
  
  #filter21$Julian <- filter21$Julian - min(filter21$Julian) + 1
  #filter22$Julian <- filter22$Julian - min(filter22$Julian) + 1
  
  #spec.nest <- rbind(spec.nest, filter21)
  spec.nest <- rbind(spec.nest, spec.surv)
}

#remove everything but the final dataframe
rm(list = ls()[!ls() %in%  c("spec.nest", "nest")])

write.csv(spec.nest, "~/Git/NDSU/RMARK/Working Data/RMarknesting.csv")