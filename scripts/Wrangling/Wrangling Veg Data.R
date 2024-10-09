################################################## Read in libraries and dataframe ##################################################

library(vegan)
library(tidyverse)
library(lubridate)

raw <- read.csv("raw/Transect Vegetation.csv")

########################################## Adding in Treatment, Pasture, and Patch Columns ##########################################

raw$Treatment <- ifelse(raw$Transect >= 97, 
                        "MTORG",
                        ifelse(raw$Transect >= 73, 
                               "PBG",
                               ifelse(raw$Transect >= 25, 
                                      "SLG",
                                      "PBG")))

raw$Pasture <- ifelse(raw$Transect >= 133,
                      "NW",
                      ifelse(raw$Transect >= 121,
                             "SW",
                             ifelse(raw$Transect >= 109,
                                    "SE",
                                    ifelse(raw$Transect >= 97,
                                           "NE",
                                           ifelse(raw$Transect >= 85,
                                                  "Dottie S",
                                                  ifelse(raw$Transect >= 73,
                                                         "Dottie N",
                                                         ifelse(raw$Transect >= 61,
                                                                "Refuge S",
                                                                ifelse(raw$Transect >= 49,
                                                                       "Refuge M",
                                                                       ifelse(raw$Transect >= 37,
                                                                              "Refuge NE",
                                                                              ifelse(raw$Transect >= 27,
                                                                                     "Refuge NW",
                                                                                     ifelse(raw$Transect >= 13,
                                                                                            "Engledorf S",
                                                                                            "Engledorf N")))))))))))

raw$Patch <- ifelse(raw$Transect %in% c(1,2,3,13,14,15,28,29,30,37,38,39,49,50,51,61,64,67,73,74,75,85,86,87),
                    "NW",
                    ifelse(raw$Transect %in% c(4,5,6,16,17,18,31,32,33,40,41,42,55,56,57,62,65,68,76,77,78,88,89,90),
                           "NE",
                           ifelse(raw$Transect %in% c(7,8,9,19,20,21,25,26,27,43,46,47,52,53,54,63,66,69,79,80,81,91,92,93),
                                  "SW",
                                  "SE")))

raw$Patch <- ifelse(raw$Transect >= 142,
                    "16",
                    ifelse(raw$Transect >= 139,
                           "15",
                           ifelse(raw$Transect >= 136,
                                  "14",
                                  ifelse(raw$Transect >= 133,
                                         "13",
                                         ifelse(raw$Transect >= 130,
                                                "12",
                                                ifelse(raw$Transect >= 127,
                                                       "11",
                                                       ifelse(raw$Transect >= 124,
                                                              "10",
                                                              ifelse(raw$Transect >= 121,
                                                                     "9",
                                                                     ifelse(raw$Transect >= 118,
                                                                            "8",
                                                                            ifelse(raw$Transect >= 115,
                                                                                   "7",
                                                                                   ifelse(raw$Transect >= 112,
                                                                                          "6",
                                                                                          ifelse(raw$Transect >= 109,
                                                                                                 "5",
                                                                                                 ifelse(raw$Transect >= 106,
                                                                                                        "4",
                                                                                                        ifelse(raw$Transect >= 103,
                                                                                                               "3",
                                                                                                               ifelse(raw$Transect >= 100,
                                                                                                                      "2",
                                                                                                                      ifelse(raw$Transect >= 97,
                                                                                                                             "1",
                                                                                                                             raw$Patch))))))))))))))))

raw$SubPatch <- ifelse(raw$Treatment == "MTORG" & raw$Year == 2021 & raw$Patch %in% c(1, 5, 12, 16), "Heavy",
                       ifelse(raw$Treatment == "MTORG" & raw$Year == 2021 & raw$Patch %in% c(2, 6, 9, 13), "Full",
                              ifelse(raw$Treatment == "MTORG" & raw$Year == 2021 & raw$Patch %in% c(3, 7, 10, 14), "Moderate",
                                     ifelse(raw$Treatment == "MTORG" & raw$Year == 2021 & raw$Patch %in% c(4, 8, 11, 15), "Rest",
                                            ifelse(raw$Treatment == "MTORG" & raw$Year == 2022 & raw$Patch %in% c(2, 6, 9, 13), "Heavy",
                                                   ifelse(raw$Treatment == "MTORG" & raw$Year == 2022 & raw$Patch %in% c(3, 7, 10, 14), "Full",
                                                          ifelse(raw$Treatment == "MTORG" & raw$Year == 2022 & raw$Patch %in% c(4, 8, 11, 15), "Moderate",
                                                                 ifelse(raw$Treatment == "MTORG" & raw$Year == 2022 & raw$Patch %in% c(1, 5, 12, 16), "Rest",
                                                                        ifelse(raw$Treatment == "MTORG" & raw$Year == 2023 & raw$Patch %in% c(3, 7, 10, 14), "Heavy",
                                                                               ifelse(raw$Treatment == "MTORG" & raw$Year == 2023 & raw$Patch %in% c(4, 8, 11, 15), "Full",
                                                                                      ifelse(raw$Treatment == "MTORG" & raw$Year == 2023 & raw$Patch %in% c(1, 5, 12, 16), "Moderate",
                                                                                             ifelse(raw$Treatment == "MTORG" & raw$Year == 2023 & raw$Patch %in% c(2, 6, 9, 13), "Rest",
                                                                                                    ""))))))))))))

write_csv(raw, "working/VegAdj.csv")