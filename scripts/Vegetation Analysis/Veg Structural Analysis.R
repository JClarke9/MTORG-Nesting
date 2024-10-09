# Load libraries ----------------------------------------------------------


library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(emmeans)


# Data import -------------------------------------------------------------


raw <- read.csv('working/VegAdj.csv') |> 
  filter(Transect.ID %in% c(97:144)) |> 
  mutate(Transect.ID = factor(Transect.ID, levels = c(97:144)),
         Year = factor(Year, levels = c('2021', '2022', 
                                        '2023', '2024')),
         Pasture = factor(Pasture, levels = c('NE', 'SE',
                                              'SW', 'NW')),
         Patch = factor(Patch, levels = c(1:16)),
         SubPatch = factor(SubPatch, levels = c('Rest', 'Moderate', 
                                                'Full', 'Heavy')),
         Date = as.Date(Date, '%m/%d/%y'))


# Data wrangling ----------------------------------------------------------


MTORG.veg <- raw |> 
  separate(col = Date,
           into = c('Day', 'Month', 'Year2'),
           sep = '-',
           remove = F) |> 
  filter(Month == '07' | Month == '08') |> 
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
           Pasture,
           Patch,
           SubPatch) |> 
  summarize(Litter.Depth = mean(Litter.Depth),
            a.Robel = mean(a.Robel)) |> 
  ungroup()

MTORG.str <- rename(MTORG.str,
                    'Replicate' = 'Pasture',
                    'cTreat' = 'SubPatch')

write.csv(MTORG.str, 'working/MTORG_str.csv')


# Setting theme -----------------------------------------------------------


windowsFonts(my_font = windowsFont('Gandhi Sans'))

my_theme <- theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = 'transparent'), 
                  plot.background = element_rect(fill = 'transparent'),
                  axis.line = element_line(colour = 'black'),
                  axis.ticks = element_line(colour = 'black'),
                  plot.title = element_text(family = 'my_font',
                                            hjust = .5,
                                            vjust = 1,
                                            size = 40,
                                            colour = 'black'),
                  axis.title.x = element_text(family = 'my_font',
                                              size = 20,
                                              colour = 'black',
                                              vjust = 1),
                  axis.title.y = element_text(family = 'my_font',
                                              size = 20,
                                              colour = 'black',
                                              angle = 90,
                                              vjust = 1),
                  axis.text.x = element_text(family = 'my_font',
                                             size = 20, 
                                             colour = 'black'),
                  axis.text.y = element_text(family = 'my_font',
                                             size = 20,
                                             colour = 'black'),
                  legend.position = 'none')


# Creating Plots ----------------------------------------------------------


(robel <- ggplot(MTORG.str,
                 aes(fill = cTreat,
                     x = cTreat,
                     y = a.Robel)) +
   geom_violin(trim = FALSE,
               draw_quantiles = c(.25, .75),
               colour = "black",
               size = 0.5) +
   scale_fill_manual(values = c('#A2A4A2', 'lightgoldenrod2', 
                                '#D4A634', '#717F5B')) +
   my_theme +
   labs(title = 'Vegetation Density', 
        x = 'Grazing Intensity', 
        y = 'Decimeter') +
   facet_wrap(~Year))

(litter <- ggplot(MTORG.str,
                  aes(fill = cTreat,
                      x = cTreat,
                      y = Litter.Depth)) +
    geom_violin(trim = FALSE,
                draw_quantiles = c(.25, .75),
                colour = "black",
                size = 0.5) +
    scale_fill_manual(values = c('#A2A4A2', 'lightgoldenrod2', 
                                 '#D4A634', '#717F5B')) +
    my_theme +
    labs(title = 'Vegetation Density', 
         x = 'Grazing Intensity', 
         y = 'Decimeter') +
    facet_wrap(~Year))


# Robel Models ------------------------------------------------------------


hist(MTORG.veg$a.Robel)
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2021 & MTORG.veg$SubPatch == 'Rest'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2021 & MTORG.veg$SubPatch == 'Moderate'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2021 & MTORG.veg$SubPatch == 'Full'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2021 & MTORG.veg$SubPatch == 'Heavy'])

aov.robel21 <- glmmTMB(a.Robel ~ SubPatch + (1|Pasture/Transect.ID),
                       data = MTORG.veg[MTORG.veg$Year == 2021,],
                       family = Gamma(link = 'identity'))

simulationOutput <- simulateResiduals(aov.robel21,
                                      n = 999,
                                      plot = T)

diagnose(aov.robel21)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.robel21)
emmeans(aov.robel21, 
        pairwise ~ SubPatch,
        adjust = 'tukey')

hist(MTORG.veg$a.Robel)
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2022 & MTORG.veg$SubPatch == 'Rest'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2022 & MTORG.veg$SubPatch == 'Moderate'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2022 & MTORG.veg$SubPatch == 'Full'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2022 & MTORG.veg$SubPatch == 'Heavy'])


aov.robel22 <- glmmTMB(a.Robel ~ SubPatch + (1|Pasture/Transect.ID),
                       data = MTORG.veg[MTORG.veg$Year == 2022,],
                       family = gaussian(link = 'log'))

simulationOutput <- simulateResiduals(aov.robel22,
                                      n = 999,
                                      plot = T)

diagnose(aov.robel22)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.robel22)
emmeans(aov.robel22, 
        pairwise ~ SubPatch,
        adjust = 'tukey')


hist(MTORG.veg$a.Robel)
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2023 & MTORG.veg$SubPatch == 'Rest'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2023 & MTORG.veg$SubPatch == 'Moderate'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2023 & MTORG.veg$SubPatch == 'Full'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2023 & MTORG.veg$SubPatch == 'Heavy'])

aov.robel23 <- glmmTMB(a.Robel ~ SubPatch + (1|Pasture/Transect.ID),
                       data = MTORG.veg[MTORG.veg$Year == 2023,],
                       dispformula = ~SubPatch,
                       family = gaussian(link = 'identity'))

simulationOutput <- simulateResiduals(aov.robel23,
                                      n = 999,
                                      plot = T)

diagnose(aov.robel23)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.robel23)
emmeans(aov.robel23, 
        pairwise ~ SubPatch,
        adjust = 'tukey')


hist(MTORG.veg$a.Robel)
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2024 & MTORG.veg$SubPatch == 'Rest'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2024 & MTORG.veg$SubPatch == 'Moderate'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2024 & MTORG.veg$SubPatch == 'Full'])
hist(MTORG.veg$a.Robel[MTORG.veg$Year == 2024 & MTORG.veg$SubPatch == 'Heavy'])

aov.robel24 <- glmmTMB(a.Robel ~ SubPatch + (1|Pasture/Transect.ID),
                       data = MTORG.veg[MTORG.veg$Year == 2024,],
                       family = gaussian(link = 'log'))

simulationOutput <- simulateResiduals(aov.robel24,
                                      n = 999,
                                      plot = T)

diagnose(aov.robel24)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.robel24)
emmeans(aov.robel24, 
        pairwise ~ SubPatch,
        adjust = 'tukey')


# Litter Depth Models -----------------------------------------------------


hist(MTORG.veg$Litter.Depth)
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2021 & MTORG.veg$SubPatch == 'Rest'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2021 & MTORG.veg$SubPatch == 'Moderate'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2021 & MTORG.veg$SubPatch == 'Full'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2021 & MTORG.veg$SubPatch == 'Heavy'])

aov.litdep21 <- glmmTMB(Litter.Depth ~ SubPatch + (1|Pasture/Transect.ID),
                        data = MTORG.veg[MTORG.veg$Year == 2021,],
                        family = Gamma(link = 'log'))

simulationOutput <- simulateResiduals(aov.litdep21,
                                      n = 999,
                                      plot = T)

diagnose(aov.litdep21)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.litdep21)
emmeans(aov.litdep21, 
        pairwise ~ SubPatch,
        adjust = 'tukey')

hist(MTORG.veg$Litter.Depth)
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2022 & MTORG.veg$SubPatch == 'Rest'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2022 & MTORG.veg$SubPatch == 'Moderate'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2022 & MTORG.veg$SubPatch == 'Full'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2022 & MTORG.veg$SubPatch == 'Heavy'])


aov.litdep22 <- glmmTMB(Litter.Depth ~ SubPatch + (1|Pasture/Transect.ID),
                        data = MTORG.veg[MTORG.veg$Year == 2022,],
                        family = gaussian(link = 'identity'))

simulationOutput <- simulateResiduals(aov.litdep22,
                                      n = 999,
                                      plot = T)

diagnose(aov.litdep22)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.litdep22)
emmeans(aov.litdep22, 
        pairwise ~ SubPatch,
        adjust = 'tukey')


hist(MTORG.veg$Litter.Depth)
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2023 & MTORG.veg$SubPatch == 'Rest'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2023 & MTORG.veg$SubPatch == 'Moderate'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2023 & MTORG.veg$SubPatch == 'Full'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2023 & MTORG.veg$SubPatch == 'Heavy'])

aov.litdep23 <- glmmTMB(Litter.Depth ~ SubPatch + (1|Pasture/Transect.ID),
                        data = MTORG.veg[MTORG.veg$Year == 2023,],
                        family = gaussian(link = 'identity'))

simulationOutput <- simulateResiduals(aov.litdep23,
                                      n = 999,
                                      plot = T)

diagnose(aov.litdep23)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.litdep23)
emmeans(aov.litdep23, 
        pairwise ~ SubPatch,
        adjust = 'tukey')


hist(MTORG.veg$Litter.Depth)
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2024 & MTORG.veg$SubPatch == 'Rest'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2024 & MTORG.veg$SubPatch == 'Moderate'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2024 & MTORG.veg$SubPatch == 'Full'])
hist(MTORG.veg$Litter.Depth[MTORG.veg$Year == 2024 & MTORG.veg$SubPatch == 'Heavy'])

aov.litdep24 <- glmmTMB(Litter.Depth ~ SubPatch + (1|Pasture/Transect.ID),
                        data = MTORG.veg[MTORG.veg$Year == 2024,],
                        family = gaussian(link = 'identity'))

simulationOutput <- simulateResiduals(aov.litdep24,
                                      n = 999,
                                      plot = T)

diagnose(aov.litdep24)
testResiduals(simulationOutput)
testZeroInflation(simulationOutput)
testQuantiles(simulationOutput, quantiles = c(0.25, 0.5, 0.75), plot = T)

summary(aov.litdep24)
emmeans(aov.litdep24, 
        pairwise ~ SubPatch,
        adjust = 'tukey')
