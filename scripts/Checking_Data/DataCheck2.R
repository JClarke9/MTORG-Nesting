library(tidyverse)

new <- read_csv("raw/Nesting.csv")
old <- read_csv("raw/Nesting_OLD.csv")

bwte.new <- filter(new, Spec == "BWTE" & Year != 2023)
bwte.old <- filter(old, Spec == "BWTE")

bwte.old <- rename(bwte.old,
                   Fate2 = Survive,
                   Paddock = Pasture,
                   Replicate = Patch)

bwte.old$Type <- "old"
bwte.new$Type <- "new"

bwte.old <- select(bwte.old,
                   c(Type, Year, id, X, Y, Spec, Date, Visit.Interval, Expos, Fate, Fate2,
                     FirstFound, Stage, BHCONum, InitBHCO, Clutch, InitClutch, AgeFound,
                     KBG:R4, Paddock, Replicate, Treatment))

bwte.new <- select(bwte.new,
                   c(Type, Year, id, X, Y, Spec, Date, Visit.Interval, Expos, Fate, Fate2,
                     FirstFound, Stage, BHCONum, InitBHCO, Clutch, InitClutch, AgeFound,
                     KBG:R4, Paddock, Replicate, Treatment))

bwte.old$Clutch <- as.numeric(bwte.old$Clutch)
bwte.old$InitClutch <- as.numeric(bwte.old$InitClutch)
bwte.new$AgeFound <- as.numeric(bwte.new$AgeFound)
bwte.new$Bare <- as.numeric(bwte.new$Bare)

bwte.new$R1 <- recode(bwte.new$R1,
                      "18+" = "20",
                      "13+" = "NA")

bwte.new$R2 <- recode(bwte.new$R2,
                      "18+" = "20",
                      "13+" = "NA",
                      "15+" = "NA")

bwte.new$R3 <- recode(bwte.new$R3,
                      "18+" = "20",
                      "13+" = "NA")

bwte.new$R4 <- recode(bwte.new$R4,
                      "18+" = "20",
                      "13+" = "NA")

test <- filter(bwte.new, R1 == "NA" | R2 == "NA" | R3 == "NA" | R4 == "NA")

bwte.new$R1 <- as.numeric(bwte.new$R1)
bwte.new$R2 <- as.numeric(bwte.new$R2)
bwte.new$R3 <- as.numeric(bwte.new$R3)
bwte.new$R4 <- as.numeric(bwte.new$R4)

bwte.diff1 <- anti_join(bwte.old, bwte.new)
bwte.diff2 <- anti_join(bwte.new, bwte.old)

bwte.diff <- bind_rows(bwte.diff1, bwte.diff2)

write_csv(bwte.diff, "outputs/checking.csv")
