library(scales)
library(cowplot)
library(tidyverse)

raw <- read_csv("raw/Nesting.csv")

raw$ref <- ifelse(raw$Year == 2021, 
                  "2021-01-01",
                  "2022-01-01")

raw$start <- as.POSIXct(raw$FirstFound,
                          format = "%m/%d/%Y") |> 
  yday()

raw$initiation <- raw$start - raw$AgeFound - 1

raw$initiation <- as.Date(raw$initiation,
                     origin = raw$ref)

init <- raw |> 
  filter(Visit.Interval == 1) |> 
  select(Year, id, Spec, FirstFound, initiation, AgeFound)

# filtering out duck species. I'm not including AMWI and GWTE because
# we have a very small sample size for both.
duck <- init |> 
  filter(Spec == "BWTE" | Spec == "GADW" | Spec == "MALL" | 
           Spec == "NOPI" | Spec == "NSHO")

duck <- duck |> 
  group_by(Year, Spec, initiation) |> 
  summarize(count = n()) |> 
  na.omit() |> 
  ungroup()

duck21 <- filter(duck, Year == 2021)
duck22 <- filter(duck, Year == 2022)

(bar21 <- ggplot(duck21, 
                 aes(x = initiation, y = count)) +
    geom_col() + 
    scale_x_date(labels = date_format("%Y-%m-%d"), 
                 breaks = date_breaks(width = "2 weeks")) +
    theme(axis.text.x = element_text(angle = 50,
                                     vjust = 1,
                                     hjust = 1)) +
    facet_wrap(~Spec))

(bar22 <- ggplot(duck22, 
                 aes(x = initiation, y = count)) +
    geom_col() + 
    scale_x_date(labels = date_format("%Y-%m-%d"), 
                 breaks = date_breaks(width = "2 weeks")) +
    theme(axis.text.x = element_text(angle = 50,
                                     vjust = 1,
                                     hjust = 1)) +
    facet_wrap(~Spec))


(bar <- plot_grid(bar21, bar22,
                  scale = .95,
                  labels = c("2021", "2022"),
                  hjust = -9.1))
