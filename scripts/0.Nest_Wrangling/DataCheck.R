library(tidyverse)

new <- read_csv("raw/Nesting.csv")
old <- read_csv("raw/Nesting_OLD.csv")

new_f <- filter(new, Year != 2023 & Year != 2024)

new_f <- new_f |> 
  group_by(id) |> 
  summarise(Paddock = mean(Paddock),
            Year = mean(Year),
            Clutch = max(Clutch))

old_filtered <- old_f |> 
  filter(InitClutch > 0 & Fate != "Unknown")

old_f <- old |> 
  group_by(id) |> 
  summarise(Paddock = mean(Pasture),
            Year = mean(Year),
            Clutch = max(Clutch))

missing1 <- anti_join(new_f, old_f,
                      by = c("id", "Year", "Paddock"))
missing2 <- anti_join(old_f, new_f,
                      by = c("id", "Year", "Paddock"))
missing3 <- rbind(missing1, missing2)

missing_f <- filter(missing3, Clutch > 0)

missing_f21 <- filter(missing_f, Year == 2021)
missing_f22 <- filter(missing_f, Year == 2022)
