# Load libraries ----------------------------------------------------------


library(tidyverse)
library(ggplot2)
library(cowplot)
library(glmmTMB)
library(DHARMa)
library(car)
library(emmeans)
library(multcomp)


# Data import -------------------------------------------------------------


raw <- read.csv("working/VegAdj.csv") |> 
  filter(Transect.ID %in% c(97:144)) |> 
  mutate(Transect.ID = factor(Transect.ID, levels = c(97:144)),
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
  dplyr::select(-c(Day, Month, Year2))

MTORG.veg$a.Robel <- MTORG.veg |> 
  dplyr::select(R1:R4) |> 
  rowMeans()

MTORG.veg$TotalVegCover <- MTORG.veg |> 
  dplyr::select(KBG:Woody) |> 
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

summarize(MTORG.str, 
          avg.Litter.Depth = mean(Litter.Depth),
          se.litdep = sd(Litter.Depth) / sqrt(n()),
          avg.Robel = mean(a.Robel), 
          se.robel = sd(a.Robel) / sqrt(n()),
          .by = Year)


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

Anova(aov.robel21, type = 2)
(em.robel21 <- emmeans(aov.robel21, 
        pairwise ~ Intensity,
        adjust = "tukey"))

cld.robel21 <- cld(em.robel21, Letters = LETTERS, alpha = 0.05)
cld.robel21 <- cld.robel21 |> 
  mutate(Year = 2021,
         .group = gsub(" ", "", .group))

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

Anova(aov.robel22, type = 2)
(em.robel22 <- emmeans(aov.robel22, 
                       pairwise ~ Intensity,
                       adjust = "tukey"))

cld.robel22 <- cld(em.robel22, Letters = LETTERS, alpha = 0.05)
cld.robel22 <- cld.robel22 |> 
  mutate(Year = 2022,
         .group = gsub(" ", "", .group))


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

Anova(aov.robel23, type = 2)
(em.robel23 <- emmeans(aov.robel23, 
                       pairwise ~ Intensity,
                       adjust = "tukey"))

cld.robel23 <- cld(em.robel23, Letters = LETTERS, alpha = 0.05)
cld.robel23 <- cld.robel23 |> 
  mutate(Year = 2023,
         .group = gsub(" ", "", .group))


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

Anova(aov.robel24, type = 2)
(em.robel24 <- emmeans(aov.robel24, 
                       pairwise ~ Intensity,
                       adjust = "tukey"))

cld.robel24 <- cld(em.robel24, Letters = LETTERS, alpha = 0.05)
cld.robel24 <- cld.robel24 |> 
  mutate(Year = 2024,
         .group = gsub(" ", "", .group))


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

Anova(aov.litdep21, type = 2)
(em.litdep21 <- emmeans(aov.litdep21, 
                       pairwise ~ Intensity,
                       adjust = "tukey"))

cld.litdep21 <- cld(em.litdep21, Letters = LETTERS, alpha = 0.05)
cld.litdep21 <- cld.litdep21 |> 
  mutate(Year = 2021,
         .group = gsub(" ", "", .group))

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

Anova(aov.litdep22, type = 2)
(em.litdep22 <- emmeans(aov.litdep22, 
                        pairwise ~ Intensity,
                        adjust = "tukey"))

cld.litdep22 <- cld(em.litdep22, Letters = LETTERS, alpha = 0.05)
cld.litdep22 <- cld.litdep22 |> 
  mutate(Year = 2022,
         .group = gsub(" ", "", .group))


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

Anova(aov.litdep23, type = 2)
(em.litdep23 <- emmeans(aov.litdep23, 
                        pairwise ~ Intensity,
                        adjust = "tukey"))

cld.litdep23 <- cld(em.litdep23, Letters = LETTERS, alpha = 0.05)
cld.litdep23 <- cld.litdep23 |> 
  mutate(Year = 2023,
         .group = gsub(" ", "", .group))


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

Anova(aov.litdep24, type = 2)
(em.litdep24 <- emmeans(aov.litdep24, 
                        pairwise ~ Intensity,
                        adjust = "tukey"))

cld.litdep24 <- cld(em.litdep24, Letters = LETTERS, alpha = 0.05)
cld.litdep24 <- cld.litdep24 |> 
  mutate(Year = 2024,
         .group = gsub(" ", "", .group))


# Merging CLD ---------------------------------------------------------------------------------


cld.robel <- bind_rows(cld.robel21,
                       cld.robel22,
                       cld.robel23,
                       cld.robel24)

cld.litdep <- bind_rows(cld.litdep21,
                        cld.litdep22,
                        cld.litdep23,
                        cld.litdep24)

cld.robel <- cld.robel |>
  mutate(.group = case_when(Year == 2021 & Intensity == "Rest" ~ "A",
                            Year == 2021 & Intensity == "Moderate" ~ "AB",
                            Year == 2021 & Intensity == "Full" ~ "AB",
                            Year == 2021 & Intensity == "Heavy" ~ "B",

                            Year == 2022 & Intensity == "Rest" ~ "A",
                            Year == 2022 & Intensity == "Moderate" ~ "A",
                            Year == 2022 & Intensity == "Full" ~ "B",
                            Year == 2022 & Intensity == "Heavy" ~ "B",

                            Year == 2023 & Intensity == "Rest" ~ "A",
                            Year == 2023 & Intensity == "Moderate" ~ "A",
                            Year == 2023 & Intensity == "Full" ~ "A",
                            Year == 2023 & Intensity == "Heavy" ~ "A",

                            Year == 2024 & Intensity == "Rest" ~ "A",
                            Year == 2024 & Intensity == "Moderate" ~ "B",
                            Year == 2024 & Intensity == "Full" ~ "B",
                            Year == 2024 & Intensity == "Heavy" ~ "B",

                            .default = NA))

cld.litdep <- cld.litdep |>
  mutate(.group = case_when(Year == 2021 & Intensity == "Rest" ~ "BC",
                            Year == 2021 & Intensity == "Moderate" ~ "A",
                            Year == 2021 & Intensity == "Full" ~ "AB",
                            Year == 2021 & Intensity == "Heavy" ~ "C",

                            Year == 2022 & Intensity == "Rest" ~ "B",
                            Year == 2022 & Intensity == "Moderate" ~ "A",
                            Year == 2022 & Intensity == "Full" ~ "AB",
                            Year == 2022 & Intensity == "Heavy" ~ "B",

                            Year == 2023 & Intensity == "Rest" ~ "AB",
                            Year == 2023 & Intensity == "Moderate" ~ "A",
                            Year == 2023 & Intensity == "Full" ~ "AB",
                            Year == 2023 & Intensity == "Heavy" ~ "B",

                            Year == 2024 & Intensity == "Rest" ~ "AB",
                            Year == 2024 & Intensity == "Moderate" ~ "A",
                            Year == 2024 & Intensity == "Full" ~ "BC",
                            Year == 2024 & Intensity == "Heavy" ~ "C",

                            .default = NA))


# Setting theme -----------------------------------------------------------


windowsFonts(my_font = windowsFont("Gandhi Sans"))

my_theme <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  plot.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  plot.title = element_text(family = "my_font",
                                            hjust = .5,
                                            vjust = 1,
                                            size = 16,
                                            color = "black"),
                  axis.title.x = element_text(family = "my_font",
                                              size = 10,
                                              color = "black",
                                              vjust = 1),
                  axis.title.y = element_text(family = "my_font",
                                              size = 10,
                                              color = "black",
                                              angle = 90,
                                              vjust = 1,
                                              margin = margin(r = 8, unit = "pt")),
                  axis.text.x = element_text(family = "my_font",
                                             size = 10,
                                             vjust = 1,
                                             hjust = 1,
                                             angle = 45,
                                             color = "black"),
                  axis.text.y = element_text(family = "my_font",
                                             size = 10,
                                             color = "black"),
                  strip.text.x = element_text(family = "my_font", 
                                              size = 12, 
                                              colour = "black"),
                  strip.background = element_blank(),
                  legend.background = element_blank(),
                  legend.position = "none")


# Creating Plots ----------------------------------------------------------


(robel <- ggplot(MTORG.str,
                 aes(fill = Intensity,
                     x = Intensity,
                     y = a.Robel)) +
   geom_violin(trim = FALSE,
               draw_quantiles = c(.25, .75),
               color = "black",
               size = 0.3) +
   geom_text(data = cld.robel, 
             aes(x = Intensity, 
                 y = 13, 
                 label = .group),
             family = "my_font",
             size = 10,
             size.unit = "pt",
             colour = "black") +
   scale_fill_manual(values = c("#A2A4A2", "lightgoldenrod2", 
                                "#D4A634", "#717F5B")) +
   my_theme +
   labs(title = NULL, 
        x = NULL, 
        y = "Vegetation Density (dm)") +
   facet_wrap(~Year, nrow = 1, ncol = 4, axes = "all", axis.labels = "margins"))

(litter <- ggplot(MTORG.str,
                  aes(fill = Intensity,
                      x = Intensity,
                      y = Litter.Depth)) +
    geom_violin(trim = FALSE,
                draw_quantiles = c(.25, .75),
                color = "black",
                size = 0.3) +
    geom_text(data = cld.litdep, 
              aes(x = Intensity, 
                  y = 110, 
                  label = .group),
              family = "my_font",
              size = 10,
              size.unit = "pt",
              colour = "black") +
    scale_fill_manual(values = c("#A2A4A2", "lightgoldenrod2", 
                                 "#D4A634", "#717F5B")) +
    my_theme +
    labs(title = NULL, 
         x = NULL, 
         y = "Litter Depth (mm)") +
    facet_wrap(~Year, nrow = 1, ncol = 4, axes = "all", axis.labels = "margins"))

(veg_violin <- plot_grid(robel + theme(axis.text.x = element_blank()),
                         litter + theme(strip.text.x = element_blank()),
                         align = "v",
                         nrow = 2,
                         ncol = 1))

ggsave(veg_violin,
       filename = "outputs/figs/VegStr_violin.png",
       bg = "white",
       dpi = 600,
       height = 4.5,
       width = 6)

