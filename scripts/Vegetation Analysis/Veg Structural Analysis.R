# Load libraries ----------------------------------------------------------


library(tidyverse)
library(ggplot2)
library(cowplot)
library(glmmTMB)
library(DHARMa)
library(emmeans)


# Data import -------------------------------------------------------------


raw <- read.csv("working/VegAdj.csv") |> 
  filter(Transect.ID %in% c(97:144)) |> 
  mutate(Transect.ID = factor(Transect.ID, levels = c(97:144)),
         Year = factor(Year, levels = c("2021", "2022", 
                                        "2023", "2024")),
         Replicate = factor(Replicate, levels = c("NE", "SE",
                                                  "SW", "NW")),
         Patch = factor(Patch, levels = c(1:16)),
         Intensity = factor(Intensity, levels = c("Rest", "Moderate", 
                                                  "Full", "Heavy")),
         Date = as.Date(Date, "%m/%d/%y"))


# Data wrangling ----------------------------------------------------------


MTORG.veg <- raw |> 
  separate(col = Date,
           into = c("Day", "Month", "Year2"),
           sep = "-",
           remove = F) |> 
  filter(Month == "07" | Month == "08") |> 
  select(-c(Day, Month, Year2))

MTORG.veg$a.Robel <- MTORG.veg |> 
  select(R1:R4) |> 
  rowMeans()

MTORG.veg$TotalVegCover <- MTORG.veg |> 
  select(KBG:Woody) |> 
  rowSums(na.rm = TRUE)

MTORG.veg <- MTORG.veg |> 
  mutate(across(KBG:Woody, ~ .x/TotalVegCover * 100))

MTORG.str <- MTORG.veg |> 
  group_by(Year,
           Replicate,
           Patch,
           Intensity) |> 
  summarize(Litter.Depth = mean(Litter.Depth),
            a.Robel = mean(a.Robel)) |> 
  ungroup()

write_csv(MTORG.str, "working/MTORG_str.csv")


# Setting theme -----------------------------------------------------------


windowsFonts(my_font = windowsFont("Gandhi Sans"))

robel_theme <- theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "transparent",
                                                     color = NA), 
                     plot.background = element_rect(fill = "transparent",
                                                    color = NA),
                     axis.line = element_line(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     plot.title = element_text(family = "my_font",
                                               hjust = .5,
                                               vjust = 1,
                                               size = 40,
                                               color = "black"),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(family = "my_font",
                                                 size = 28,
                                                 color = "black",
                                                 angle = 90,
                                                 vjust = 1),
                     axis.text.x = element_blank(),
                     axis.text.y = element_text(family = "my_font",
                                                size = 28,
                                                color = "black"),
                     strip.text.x = element_text(family = 'my_font',
                                                 size = 36,
                                                 vjust = 3),
                     strip.background = element_rect(fill = "transparent",
                                                     color = NA),
                     legend.position = "none")

litter_theme <- theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill = "transparent",
                                                      color = NA), 
                      plot.background = element_rect(fill = "transparent",
                                                     color = NA),
                      axis.line = element_line(color = "black"),
                      axis.ticks = element_line(color = "black"),
                      plot.title = element_text(family = "my_font",
                                                hjust = .5,
                                                vjust = 1,
                                                size = 40,
                                                color = "black"),
                      axis.title.x = element_text(family = "my_font",
                                                  size = 28,
                                                  color = "black",
                                                  vjust = 1),
                      axis.title.y = element_text(family = "my_font",
                                                  size = 28,
                                                  color = "black",
                                                  angle = 90,
                                                  vjust = 1),
                      axis.text.x = element_text(family = "my_font",
                                                 size = 28,
                                                 vjust = 1,
                                                 hjust = 1,
                                                 angle = 45,
                                                 color = "black"),
                      axis.text.y = element_text(family = "my_font",
                                                 size = 28,
                                                 color = "black"),
                      strip.text.x = element_blank(),
                      legend.position = "none")


# Creating Plots ----------------------------------------------------------


(robel <- ggplot(MTORG.str,
                 aes(fill = Intensity,
                     x = Intensity,
                     y = a.Robel)) +
   geom_violin(trim = FALSE,
               draw_quantiles = c(.25, .75),
               color = "black",
               size = 0.5) +
   scale_fill_manual(values = c("#A2A4A2", "lightgoldenrod2", 
                                "#D4A634", "#717F5B")) +
   robel_theme +
   labs(title = NULL, 
        x = "Grazing Intensity", 
        y = "Vegetation Density (dm)") +
   facet_wrap(~Year, nrow = 1, ncol = 4, axes = "all", axis.labels = "margins"))

(litter <- ggplot(MTORG.str,
                  aes(fill = Intensity,
                      x = Intensity,
                      y = Litter.Depth)) +
    geom_violin(trim = FALSE,
                draw_quantiles = c(.25, .75),
                color = "black",
                size = 0.5) +
    scale_fill_manual(values = c("#A2A4A2", "lightgoldenrod2", 
                                 "#D4A634", "#717F5B")) +
    litter_theme +
    labs(title = NULL, 
         x = NULL, 
         y = "Litter Depth (mm)") +
    facet_wrap(~Year, nrow = 1, ncol = 4, axes = "all", axis.labels = "margins"))

(veg_violin <- plot_grid(robel,
                         litter,
                         align = "v",
                         nrow = 2,
                         ncol = 1))

ggsave(veg_violin,
       filename = "outputs/figs/VegStr_violin.png",
       bg = "transparent",
       dpi = 600,
       height = 10.78,
       width = 22.11)


# Robel Models ------------------------------------------------------------


hist(MTORG.veg$a.Robel)
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2021 & MTORG.veg$Intensity == "Rest"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2021 & MTORG.veg$Intensity == "Moderate"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2021 & MTORG.veg$Intensity == "Full"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2021 & MTORG.veg$Intensity == "Heavy"])

aov.robel21 <- glmmTMB(a.Robel ~ Intensity + (1|Replicate/Transect.ID),
                       data = MTORG.veg[MTORG.veg$Year == 2021,],
                       family = Gamma(link = "identity"))

simulationOutput <- simulateResiduals(aov.robel21,
                                      n = 999,
                                      plot = T)

diagnose(aov.robel21)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.robel21)
emmeans(aov.robel21, 
        pairwise ~ Intensity,
        adjust = "tukey")

hist(MTORG.veg$a.Robel)
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2022 & MTORG.veg$Intensity == "Rest"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2022 & MTORG.veg$Intensity == "Moderate"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2022 & MTORG.veg$Intensity == "Full"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2022 & MTORG.veg$Intensity == "Heavy"])


aov.robel22 <- glmmTMB(a.Robel ~ Intensity + (1|Replicate/Transect.ID),
                       data = MTORG.veg[MTORG.veg$Year == 2022,],
                       family = gaussian(link = "log"))

simulationOutput <- simulateResiduals(aov.robel22,
                                      n = 999,
                                      plot = T)

diagnose(aov.robel22)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.robel22)
emmeans(aov.robel22, 
        pairwise ~ Intensity,
        adjust = "tukey")


hist(MTORG.veg$a.Robel)
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2023 & MTORG.veg$Intensity == "Rest"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2023 & MTORG.veg$Intensity == "Moderate"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2023 & MTORG.veg$Intensity == "Full"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2023 & MTORG.veg$Intensity == "Heavy"])

aov.robel23 <- glmmTMB(a.Robel ~ Intensity + (1|Replicate/Transect.ID),
                       data = MTORG.veg[MTORG.veg$Year == 2023,],
                       dispformula = ~Intensity,
                       family = gaussian(link = "identity"))

simulationOutput <- simulateResiduals(aov.robel23,
                                      n = 999,
                                      plot = T)

diagnose(aov.robel23)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.robel23)
emmeans(aov.robel23, 
        pairwise ~ Intensity,
        adjust = "tukey")


hist(MTORG.veg$a.Robel)
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2024 & MTORG.veg$Intensity == "Rest"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2024 & MTORG.veg$Intensity == "Moderate"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2024 & MTORG.veg$Intensity == "Full"])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2024 & MTORG.veg$Intensity == "Heavy"])

aov.robel24 <- glmmTMB(a.Robel ~ Intensity + (1|Replicate/Transect.ID),
                       data = MTORG.veg[MTORG.veg$Year == 2024,],
                       family = gaussian(link = "log"))

simulationOutput <- simulateResiduals(aov.robel24,
                                      n = 999,
                                      plot = T)

diagnose(aov.robel24)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.robel24)
emmeans(aov.robel24, 
        pairwise ~ Intensity,
        adjust = "tukey")


# Litter Depth Models -----------------------------------------------------


hist(MTORG.veg$Litter.Depth)
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2021 & MTORG.veg$Intensity == "Rest"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2021 & MTORG.veg$Intensity == "Moderate"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2021 & MTORG.veg$Intensity == "Full"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2021 & MTORG.veg$Intensity == "Heavy"])

aov.litdep21 <- glmmTMB(Litter.Depth ~ Intensity + (1|Replicate/Transect.ID),
                        data = MTORG.veg[MTORG.veg$Year == 2021,],
                        family = Gamma(link = "log"))

simulationOutput <- simulateResiduals(aov.litdep21,
                                      n = 999,
                                      plot = T)

diagnose(aov.litdep21)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.litdep21)
emmeans(aov.litdep21, 
        pairwise ~ Intensity,
        adjust = "tukey")

hist(MTORG.veg$Litter.Depth)
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2022 & MTORG.veg$Intensity == "Rest"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2022 & MTORG.veg$Intensity == "Moderate"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2022 & MTORG.veg$Intensity == "Full"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2022 & MTORG.veg$Intensity == "Heavy"])


aov.litdep22 <- glmmTMB(Litter.Depth ~ Intensity + (1|Replicate/Transect.ID),
                        data = MTORG.veg[MTORG.veg$Year == 2022,],
                        family = gaussian(link = "identity"))

simulationOutput <- simulateResiduals(aov.litdep22,
                                      n = 999,
                                      plot = T)

diagnose(aov.litdep22)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.litdep22)
emmeans(aov.litdep22, 
        pairwise ~ Intensity,
        adjust = "tukey")


hist(MTORG.veg$Litter.Depth)
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2023 & MTORG.veg$Intensity == "Rest"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2023 & MTORG.veg$Intensity == "Moderate"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2023 & MTORG.veg$Intensity == "Full"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2023 & MTORG.veg$Intensity == "Heavy"])

aov.litdep23 <- glmmTMB(Litter.Depth ~ Intensity + (1|Replicate/Transect.ID),
                        data = MTORG.veg[MTORG.veg$Year == 2023,],
                        family = gaussian(link = "identity"))

simulationOutput <- simulateResiduals(aov.litdep23,
                                      n = 999,
                                      plot = T)

diagnose(aov.litdep23)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.litdep23)
emmeans(aov.litdep23, 
        pairwise ~ Intensity,
        adjust = "tukey")


hist(MTORG.veg$Litter.Depth)
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2024 & MTORG.veg$Intensity == "Rest"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2024 & MTORG.veg$Intensity == "Moderate"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2024 & MTORG.veg$Intensity == "Full"])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2024 & MTORG.veg$Intensity == "Heavy"])

aov.litdep24 <- glmmTMB(Litter.Depth ~ Intensity + (1|Replicate/Transect.ID),
                        data = MTORG.veg[MTORG.veg$Year == 2024,],
                        family = gaussian(link = "identity"))

simulationOutput <- simulateResiduals(aov.litdep24,
                                      n = 999,
                                      plot = T)

diagnose(aov.litdep24)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.litdep24)
emmeans(aov.litdep24, 
        pairwise ~ Intensity,
        adjust = "tukey")
