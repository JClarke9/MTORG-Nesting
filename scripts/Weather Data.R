# Load libraries ----------------------------------------------------------


library(tidyverse)
library(ggplot2)
library(cowplot)

windowsFonts(my_font = windowsFont("Tahoma"))

# Data import -------------------------------------------------------------

# Daily weather data from ndawn including: air temp (max), air temp (min), air temp (avg), and rainfall (total)
# I also included departures from normal (30 year average)
# dates are January 1, 1990 to present (September 10, 2024)
# all measurements are metric, precip is in mm
raw.ndawn <- read_csv("raw/DailyWeather.csv", skip = 4)


# Data import -------------------------------------------------------------


raw.ndawn$Day <- as.numeric(raw.ndawn$Day)
raw.ndawn$Month <- as.numeric(raw.ndawn$Month)
raw.ndawn$Year <- as.factor(raw.ndawn$Year)


# Data Wrangling ----------------------------------------------------------


ndawn <- raw.ndawn |> 
  dplyr::select(Year, Month, Day, `Max Temp`, `Normal Max Temp`, `Min Temp`, 
         `Normal Min Temp`,`Avg Temp`, `Normal Avg Temp`, Rainfall, 
         `Normal Daily Total Rainfall`) |> 
  filter(Month %in% c(4, 5, 6, 7, 8, 9, 10) & Year %in% c(2021, 2022, 2023, 2024))

ndawn <- ndawn |> 
  rename(max.temp = "Max Temp",
         norm.max.temp = "Normal Max Temp",
         min.temp = "Min Temp",
         norm.min.temp = "Normal Min Temp",
         avg.temp = "Avg Temp",
         norm.avg.temp = "Normal Avg Temp",
         precip = "Rainfall",
         norm.precip = "Normal Daily Total Rainfall") |> 
  mutate(norm.precip = norm.precip/10,
         precip = precip/10)


# Summarizing the Data ----------------------------------------------------


(growing.year <- ndawn %>%                                                   # select the data frame
   group_by(Year) %>%                                                            # group by year
   summarize(min.temp = mean(min.temp),
             max.temp = mean(max.temp),
             temp = mean(avg.temp),                                                # create a new column with temps averaged
             norm.temp = mean(norm.avg.temp),                                    # create a new column with normal temps averaged
             precip = sum(precip),                                               # create a new column with rainfall summed for each year
             norm.precip = sum(norm.precip)))


# average across the whole study
(averages.d <- ndawn %>%                                                  # select the data frame
    summarize(min.temp = mean(min.temp),
              max.temp = mean(max.temp),
              temp = mean(avg.temp),                                                # create a new column with temps averaged
              norm.temp = mean(norm.avg.temp),                                    # create a new column with normal temps averaged
              precip = sum(precip),                                               # create a new column with rainfall summed for each year
              norm.precip = sum(norm.precip)))


# Defining Theme ------------------------------------------------------------------------------


my_theme <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank(),
                  axis.line = element_line(color = "black"),
                  axis.ticks = element_line(color = "black"),
                  plot.title = element_text(family = "my_font",
                                            hjust = 0.5,
                                            size = 14,
                                            vjust = 1,
                                            face = "bold"),
                  axis.title.x = element_text(family = "my_font",
                                              color = "black",
                                              size = 10), 
                  axis.title.y = element_text(family = "my_font",
                                              color = "black",
                                              size = 10),
                  axis.text.x = element_text(family ="my_font",
                                             color = "black",
                                             size = 10,
                                             angle = 45, 
                                             hjust = 1),
                  axis.text.y = element_text(family = "my_font",
                                             color = "black",
                                             size = 10),
                  strip.text.x = element_text(family = "my_font",
                                              color = "black",
                                              size = 12),
                  strip.background = element_blank(),
                  legend.position = "none")


# Plotting the data -------------------------------------------------------


(precip <- ggplot(ndawn[ndawn$Year %in% c("2021", "2022", "2023", "2024"),],
                  aes(x = Month,
                      y = precip)) +
    stat_summary(aes(y = norm.precip),
                 fun = sum,
                 na.rm = T,
                 geom = "bar",
                 fill = "gray90") +
    stat_summary(fun = sum,
                 na.rm = T,
                 geom = "point",
                 colour = "blue",
                 size = 1.5) +
    stat_summary(fun = sum,
                 na.rm = T,
                 geom = "line",
                 colour = "blue",
                 linewidth = 1) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(from = 4,to = 10, by = 1),
                       labels = c("4"  = "Apr",
                                  "5"  = "May",
                                  "6"  = "Jun",
                                  "7"  = "Jul",
                                  "8"  = "Aug",
                                  "9"  = "Sep",
                                  "10" = "Oct")) +
    coord_cartesian(ylim = c(0, 20)) +
    my_theme +
    labs(x = NULL, y = "Cumulative Precipitation (cm)", title = NULL) +
    facet_grid(~Year))

ggsave("outputs/figs/precip.png",
       precip,
       bg = "white",
       dpi = 600,
       height = 3,
       width = 6.45)


