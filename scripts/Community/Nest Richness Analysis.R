# Loading libraries -------------------------------------------------------

library(vegan)
library(tidyverse)
library(RColorBrewer)
library(emmeans)
library(cowplot)
library(wesanderson)

# Data import -------------------------------------------------------------

totals21 <- read.csv("working/totals21.csv", 
                     row.names = 1)                             # read in the data set

totals22 <- read.csv("working/totals22.csv", 
                     row.names = 1)                             # read in the data set

windowsFonts(my_font = windowsFont("Gandhi Sans"))

# 2021 site by species matrix ---------------------------------------------

birds21 <- pivot_wider(totals21[c(1:3,5)],                                      # select the data frame to turn into a matrix
                       names_from = Spec,                                       # select the column names
                       values_from = Abundance,                                 # select the abundance values for each species
                       values_fill = list(Abundance = 0)) |>                    # fill all NA with 0
  column_to_rownames("Pasture")                                                 # set column names as the pasture ID

write.csv(birds21, 
          file = "working/birds21.csv")

obl.birds21 <- totals21 |> 
  filter(Group == "OBL") |> 
  pivot_wider(names_from = Spec,
              values_from = Abundance,
              values_fill = list(Abundance = 0)) |> 
  column_to_rownames("Pasture")

fac.birds21 <- totals21 |> 
  filter(Group == "FAC") |> 
  pivot_wider(names_from = Spec,
              values_from = Abundance,
              values_fill = list(Abundance = 0)) |> 
  column_to_rownames("Pasture")

wet.birds21 <- totals21 |> 
  filter(Group == "WET") |> 
  pivot_wider(names_from = Spec,
              values_from = Abundance,
              values_fill = list(Abundance = 0)) |> 
  column_to_rownames("Pasture")

#calculate species richness for all birds

tot.rich21 <- specnumber(birds21[2:25],
                         birds21$cTreat,
                         MARGIN=1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(tot.rich21) <- c("treat", 
                          "tot.rich21")

tot.rich21

#calculate species richness for obligate birds

obl.rich21 <- specnumber(obl.birds21[3:10],
                         obl.birds21$cTreat,
                         MARGIN=1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(obl.rich21) <- c("treat", 
                          "obl.rich21")

obl.rich21

#calculate species richness for facultative birds

fac.rich21 <- specnumber(fac.birds21[3:15],
                         fac.birds21$cTreat,
                         MARGIN=1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(fac.rich21) <- c("treat", 
                          "fac.rich21")

fac.rich21

#calculate species richness for wetland birds

wet.rich21 <- specnumber(wet.birds21[3:5],
                         wet.birds21$cTreat,
                         MARGIN=1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(wet.rich21) <- c("treat", 
                          "wet.rich21")

wet.rich21

rich21 <- full_join(tot.rich21, 
                    obl.rich21, 
                    by="treat") |> 
  full_join(fac.rich21, 
            by="treat") |> 
  full_join(wet.rich21, 
            by="treat")

rich21

# 2022 site by species matrix ---------------------------------------------

birds22 <- pivot_wider(totals22[c(1:3,5)],                                      # select the data frame to turn into a matrix
                       names_from = Spec,                                      # select the column names
                       values_from = Abundance,                                # select the abundance values for each species
                       values_fill = list(Abundance = 0)) |>                   # fill all NA with 0
  column_to_rownames("Pasture")                                                 # set column names as the pasture ID

write.csv(birds22, file = "working/birds22.csv")

obl.birds22 <- totals22 |> 
  filter(Group == "OBL") |> 
  pivot_wider(names_from = Spec,
              values_from = Abundance,
              values_fill = list(Abundance = 0)) |> 
  column_to_rownames("Pasture")

fac.birds22 <- totals22 |> 
  filter(Group == "FAC") |> 
  pivot_wider(names_from = Spec,
              values_from = Abundance,
              values_fill = list(Abundance = 0)) |> 
  column_to_rownames("Pasture")

wet.birds22 <- totals22 |> 
  filter(Group == "WET") |> 
  pivot_wider(names_from = Spec,
              values_from = Abundance,
              values_fill = list(Abundance = 0)) |> 
  column_to_rownames("Pasture")

# 2022 Richness calculations ----------------------------------------------

#calculate species richness for all birds

tot.rich22 <- specnumber(birds22[2:28],
                         birds22$cTreat,
                         MARGIN=1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(tot.rich22) <- c("treat", 
                          "tot.rich22")

tot.rich22

#calculate species richness for obligate birds

obl.rich22 <- specnumber(obl.birds22[3:10],
                         obl.birds22$cTreat,
                         MARGIN=1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(obl.rich22) <- c("treat", 
                          "obl.rich22")

obl.rich22

#calculate species richness for facultative birds

fac.rich22 <- specnumber(fac.birds22[3:10],
                         fac.birds22$cTreat,
                         MARGIN=1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(fac.rich22) <- c("treat", 
                          "fac.rich22")

fac.rich22

#calculate species richness for generalist birds

wet.rich22 <- specnumber(wet.birds22[3:8],
                         wet.birds22$cTreat,
                         MARGIN=1) |> 
  data.frame() |> 
  rownames_to_column()

colnames(wet.rich22) <- c("treat", 
                          "wet.rich22")

wet.rich22

rich22 <- full_join(tot.rich22, 
                    obl.rich22, 
                    by="treat") |> 
  full_join(fac.rich22, 
            by="treat") |> 
  full_join(wet.rich22, 
            by="treat")

rich22

rich <- full_join(tot.rich21, 
                  tot.rich22, 
                  by = "treat") |> 
  full_join(obl.rich21, 
            by="treat") |> 
  full_join(obl.rich22, 
            by="treat") |> 
  full_join(fac.rich21, 
            by="treat") |> 
  full_join(fac.rich22, 
            by="treat") |> 
  full_join(wet.rich21, 
            by="treat") |> 
  full_join(wet.rich22, 
            by="treat")

rich <- rich |> 
  rename("Grazing Intensity" = "treat",
         "Total Richness 2021" = "tot.rich21",
         "Total Richness 2022" = "tot.rich22",
         "OBL Richness 2021" = "obl.rich21",
         "OBL Richness 2022" = "obl.rich22",
         "FAC Richness 2021" = "fac.rich21",
         "FAC Richness 2022" = "fac.rich22",
         "WET Richness 2021" = "wet.rich21",
         "WET Richness 2022" = "wet.rich22")

write.csv(rich, file = "working/richness.csv")

# Manipulating data -------------------------------------------------------

past.obl21 <- obl.birds21 |> 
  rownames_to_column(var = "Pasture") |> 
  select(-Group)

past.fac21 <- fac.birds21 |> 
  rownames_to_column(var = "Pasture") |> 
  select(-Group)

past.rich21 <- full_join(past.obl21, past.fac21)

past.rich21 <- replace(past.rich21, is.na(past.rich21), 0)

past.obl22 <- obl.birds22 |> 
  rownames_to_column(var = "Pasture") |> 
  select(-Group)

past.fac22 <- fac.birds22 |> 
  rownames_to_column(var = "Pasture") |> 
  select(-Group)

past.rich22 <- full_join(past.obl22, past.fac22)

past.rich22 <- replace(past.rich22, is.na(past.rich22), 0)

write.csv(past.rich21, file = "working/richness21.csv")
write.csv(past.rich22, file = "working/richness22.csv")


obl.rich21 <- data.frame("Richness"=specnumber(past.rich21[3:10],
                                                past.rich21$Pasture,
                                                MARGIN=1),
                          Treat = past.rich21$cTreat,
                          Year=2021,
                          Pasture = past.rich21$Pasture)

obl.aov21 <- glm(Richness ~ Treat, 
                 obl.rich21,
                 family = poisson(link = "log"))
summary(obl.aov21)

emmeans(obl.aov21,
        pairwise~Treat)

fac.rich21 <- data.frame("Richness"=specnumber(past.rich21[10:23],
                                                past.rich21$Pasture,
                                                MARGIN=1),
                          Treat = past.rich21$cTreat,
                          Year=2021,
                          Pasture = past.rich21$Pasture)

fac.aov21 <- glm(Richness ~ Treat, 
               fac.rich21,
               family = poisson(link = "log"))
summary(fac.aov21)

emmeans(fac.aov21,
        pairwise~Treat)

obl.rich22 <- data.frame("Richness"=specnumber(past.rich22[3:10],
                                               past.rich22$Pasture,
                                               MARGIN=1),
                         Treat = past.rich22$cTreat,
                         Year=2022,
                         Pasture = past.rich22$Pasture)

obl.aov22 <- glm(Richness ~ Treat, 
                obl.rich22,
                family = poisson(link = "log"))
summary(obl.aov22)

emmeans(obl.aov22,
        pairwise~Treat)

fac.rich22 <- data.frame("Richness"=specnumber(past.rich22[10:23],
                                               past.rich22$Pasture,
                                               MARGIN=1),
                         Treat = past.rich22$cTreat,
                         Year=2022,
                         Pasture = past.rich22$Pasture)

fac.aov22 <- glm(Richness ~ Treat, 
                 fac.rich22,
                 family = poisson(link = "log"))
summary(fac.aov22)

emmeans(fac.aov22,
        pairwise~Treat)

obl.rich <- rbind(obl.rich21,
                  obl.rich22) |> 
  ungroup()

fac.rich <- rbind(fac.rich21,
                  fac.rich22) |> 
  ungroup()

obl.aov <- glm(Richness ~ Treat,
               obl.rich,
               family = poisson(link = "log"))
summary(obl.aov)

emmeans(obl.aov,
        pairwise~Year)

fac.aov <- glm(Richness ~ Treat,
               fac.rich,
               family = poisson(link = "log"))
summary(fac.aov)

emmeans(fac.aov,
        pairwise~Treat)


past.rich21 <- data.frame("Richness"=specnumber(past.rich21[3:23],
                                                past.rich21$Pasture,
                                                MARGIN=1),
                          Treat = past.rich21$cTreat,
                          Year=2021,
                          Pasture = past.rich21$Pasture)

past.rich22 <- data.frame("Richness"=specnumber(past.rich22[3:23],
                                                past.rich22$Pasture,
                                                MARGIN=1),
                          Treat=past.rich22$cTreat,
                          Year=2022,
                          Pasture=past.rich22$Pasture)

tot.rich <- rbind(past.rich21,
                  past.rich22) |> 
  ungroup()

tot.rich$Year <- as.factor(tot.rich$Year)

tot.aov <- glm(Richness ~ Treat, 
               tot.rich,
               family = poisson(link = "log"))
summary(tot.aov)

emmeans(tot.aov,
        pairwise~Year)

write.csv(tot.rich, "working/TotalRichness.csv")

# Creating species richness plots -----------------------------------------

rich.box <- ggplot(tot.rich,                                                  # select the data to graph
                   aes(fill=Year,
                       x=Treat,                                                # define the x axis
                       y=Richness)) +
  geom_boxplot(colour= "#273046",                                                 # create black outlines around the bar plot
               size=1,
               notch = FALSE, 
               outlier.color = NULL,
               outlier.size = NA,
               outlier.shape = NA,
               fatten=1) +
  #stat_summary(geom = 'pointrange',
  #             fun.data = "mean_sdl",
  #             fun.args = list(mult = 1),
  #             color = 'black',
  #             size = 1.15,
  #             position = position_dodge(0.75)) +
  scale_fill_manual(values=c("#798E87", "#CCC591")) +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill=NA,                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill=NA,                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.y = element_text(size=20, colour = "black"),                    # color the axis text
        axis.text.x = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=24,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  labs(title = "Breeding Bird Richness", 
       x = NULL, 
       y = "Species Richness") +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", 
                            "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, 
                  ylim=c(0,14))

rich.box

ggsave(rich.box, 
       filename = "outputs/figs/RichnessBox.png",  
       dpi = "retina", 
       bg = "white",
       width = 6,
       height = 6.5)

birds.div <- plot_grid(rich.box, simp.box)

birds.div <- birds.div + plot_annotation(theme=theme(panel.grid.major = element_blank(), # remove the vertical grid lines
                                                     panel.grid.minor = element_blank(), # remove the horizontal grid lines
                                                     panel.background = element_rect(fill=NA, # make the interior background transparent
                                                                                     colour = NA), 
                                                     plot.background = element_rect(fill =NA, # make the outer background transparent
                                                                                    colour = NA), 
                                                     axis.line = element_line(colour = "black"), # color the x and y axis
                                                     axis.text = element_text(colour = "black"), # color the axis text
                                                     axis.ticks = element_line(colour = "black"),
                                                     text=element_text(colour = "black"),# change the color of the axis titles
                                                     legend.position = "none"))

birds.div

ggsave(birds.div, 
       filename = "outputs/figs/AvianDiversity.png",  
       dpi = "print", 
       bg = "transparent",
       width = 12.36,
       height = 11)

# Manipulating data -------------------------------------------------------

totals21.bar <- totals21 |> 
  filter(Group != "GEN") |> 
  group_by(Pasture,
           cTreat,
           Group) |> 
  summarize(Abundance = sum(Abundance)) |> 
  ungroup()

totals21.bar <- complete(totals21.bar,
                         nesting(Pasture,
                                 cTreat),
                         Group,
                         fill=list(Abundance=0))

totals.bar <- rbind(totals21.bar,
                    totals22.bar) |> 
  filter(Group != "WET") |> 
  group_by(Pasture,
           cTreat,
           Group) |> 
  summarize(Abundance = sum(Abundance)) |> 
  ungroup()

totals22.bar <- totals22 |> 
  filter(Group != "WET") |> 
  group_by(Pasture,
           cTreat,
           Group) |> 
  summarize(Abundance = sum(Abundance)) |> 
  ungroup()

totals22.bar <- complete(totals22.bar,
                         nesting(Pasture,
                                 cTreat),
                         Group,
                         fill=list(Abundance=0))

totals.bar <- rbind(totals21.bar,
                    totals22.bar) |> 
  filter(Group != "WET") |> 
  group_by(Pasture,
           cTreat,
           Group) |> 
  summarize(Abundance = sum(Abundance)) |> 
  ungroup()

FACtotals.bar <- filter(totals.bar, 
                        Group == "FAC")
OBLtotals.bar <- filter(totals.bar, 
                        Group == "OBL")
FAC21.bar <- filter(totals21.bar, 
                    Group == "FAC")
OBL21.bar <- filter(totals21.bar, 
                    Group == "OBL")
FAC22.bar <- filter(totals22.bar, 
                    Group == "FAC")
OBL22.bar <- filter(totals22.bar, 
                    Group == "OBL")

# Testing for differences between TRT -------------------------------------

obl.model <- glm(Abundance ~ cTreat,                                          # test whether average robel measurements were different in each grazing intensity
                 data=OBLtotals.bar,                                              # select the data
                 family=poisson(link="log"))                                            # use family Poisson because it is count data
summary(obl.model) 

emmeans(obl.model,                                                            # select the model
        pairwise ~ cTreat)    

# This runs the obligate with years broken out

obl.model21 <- glm(Abundance ~ cTreat,                                          # test whether average robel measurements were different in each grazing intensity
                   data=OBL21.bar,                                              # select the data
                   family=poisson(link="log"))                                  # use family Poisson because it is count data
summary(obl.model21) 

emmeans(obl.model21,                                                            # select the model
        pairwise ~ cTreat)                                                      # compare each treatment to each other

obl.model22 <- glm(Abundance ~ cTreat,                                          # test whether average robel measurements were different in each grazing intensity
                   data=OBL22.bar,                                              # select the data
                   family=poisson(link="log"))                                              # use family Gamma because it is count data
summary(obl.model22) 

emmeans(obl.model22,                                                            # select the model
        pairwise ~ cTreat)                                                      # compare each treatment to each other

# This runs the facultative models with are broken down by year

fac.model <- glm(Abundance ~ cTreat,                                            # test whether average robel measurements were different in each grazing intensity
                 data=FACtotals.bar,                                          # select the data
                 family=poisson(link="log"))                                  # use family Gamma because it is count data
summary(fac.model) 

emmeans(fac.model,                                                              # select the model
        pairwise ~ cTreat)    

fac.model21 <- glm(Abundance ~ cTreat,                                          # test whether average robel measurements were different in each grazing intensity
                   data=FAC21.bar,                                              # select the data
                   family=poisson(link="log"))                                              # use family Gamma because it is count data
summary(fac.model21) 

emmeans(fac.model21,                                                            # select the model
        pairwise ~ cTreat)                                                      # compare each treatment to each other


fac.model22 <- glm(Abundance ~ cTreat,                                          # test whether average robel measurements were different in each grazing intensity
                   data=FAC22.bar,                                              # select the data
                   family=poisson(link="log"))                                              # use family Gamma because it is count data
summary(fac.model22) 

emmeans(fac.model22,                                                            # select the model
        pairwise ~ cTreat)                                                      # compare each treatment to each other

# Creating OBG and FAC bar plots ------------------------------------------

obl.plot21 <- ggplot(OBL21.bar,                                                  # select the data to graph
                     aes(x=cTreat,                                                # define the x axis
                         y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod4") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="black",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="black",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.y = element_text(size=20, colour = "black"),                    # color the axis text
        axis.text.x = element_blank(),
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = "2021 Average Avian Abundance", 
       x = NULL, 
       y = "Mean OBL Abundance") +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", 
                            "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0,10))

obl.plot21

fac.plot21 <- ggplot(FAC21.bar,                                                  # select the data to graph
                     aes(x=cTreat,                                                # define the x axis
                         y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod3") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="black",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="black",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = NULL, x=NULL, y = "Mean FAC Abundance") +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0,50))

fac.plot21

birds21.bar <- (obl.plot21/fac.plot21)

birds21.bar <- birds21.bar + plot_annotation(theme=theme(panel.grid.major = element_blank(), # remove the vertical grid lines
                                                         panel.grid.minor = element_blank(), # remove the horizontal grid lines
                                                         panel.background = element_rect(fill="black", # make the interior background transparent
                                                                                         colour = NA), 
                                                         plot.background = element_rect(fill ="black", # make the outer background transparent
                                                                                        colour = NA), 
                                                         axis.line = element_line(colour = "black"), # color the x and y axis
                                                         axis.text = element_text(colour = "black"), # color the axis text
                                                         axis.ticks = element_line(colour = "black"),
                                                         text=element_text(colour = "black"),# change the color of the axis titles
                                                         legend.position = "none"))

birds21.bar

obl.plot22 <- ggplot(OBL22.bar,                                                  # select the data to graph
                     aes(x=cTreat,                                                # define the x axis
                         y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod4") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="black",                           # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="black",                            # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.y = element_blank(),                                          # color the axis text
        axis.text.x = element_blank(),
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_blank(),                                                   # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = "2022 Average Avian Abundance", 
       x = NULL, 
       y = NULL) +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", 
                            "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0,10))

obl.plot22

fac.plot22 <- ggplot(FAC22.bar,                                                  # select the data to graph
                     aes(x=cTreat,                                                # define the x axis
                         y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod3") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="black",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="black",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = NULL, x=NULL, y = NULL) +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0,50))

fac.plot22

# Combing OBL and FAC bar plots -------------------------------------------

birds22.bar <- (obl.plot22/fac.plot22)

birds22.bar <- birds22.bar + plot_annotation(theme=theme(panel.grid.major = element_blank(), # remove the vertical grid lines
                                                         panel.grid.minor = element_blank(), # remove the horizontal grid lines
                                                         panel.background = element_rect(fill="black", # make the interior background transparent
                                                                                         colour = NA), 
                                                         plot.background = element_rect(fill ="black", # make the outer background transparent
                                                                                        colour = NA), 
                                                         axis.line = element_line(colour = "black"), # color the x and y axis
                                                         axis.text = element_text(colour = "black"), # color the axis text
                                                         axis.ticks = element_line(colour = "black"),
                                                         text=element_text(colour = "black"),# change the color of the axis titles
                                                         legend.position = "none"))

birds22.bar

birds.bar <- (birds21.bar|birds22.bar)

birds.bar <- birds.bar + plot_annotation(theme=theme(panel.grid.major = element_blank(), # remove the vertical grid lines
                                                     panel.grid.minor = element_blank(), # remove the horizontal grid lines
                                                     panel.background = element_rect(fill="black", # make the interior background transparent
                                                                                     colour = NA), 
                                                     plot.background = element_rect(fill ="black", # make the outer background transparent
                                                                                    colour = NA), 
                                                     axis.line = element_line(colour = "black"), # color the x and y axis
                                                     axis.text = element_text(colour = "black"), # color the axis text
                                                     axis.ticks = element_line(colour = "black"),
                                                     text=element_text(colour = "black"),# change the color of the axis titles
                                                     legend.position = "none"))

birds.bar

# Creating totals bar plots -----------------------------------------------

obl.plot <- ggplot(OBLtotals.bar,                                                  # select the data to graph
                   aes(x=cTreat,                                                # define the x axis
                       y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod4") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="black",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="black",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = "Mean Avian Abundance", x=NULL, y = "Mean OBL Abundance") +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0,10))

obl.plot

fac.plot <- ggplot(FACtotals.bar,                                                  # select the data to graph
                   aes(x=cTreat,                                                # define the x axis
                       y=Abundance)) +
  geom_boxplot(aes(fill=Group),
               colour= "black",                                                       # create black outlines around the bar plot
               size=.5,
               notch = FALSE, 
               outlier.color = "black",
               outlier.size = 3,
               outlier.shape = NA) +                                                      # not sure what this does
  scale_fill_manual(values="goldenrod3") +                                   # select the color for group 2
  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                  size=30,
                                  hjust=.5),
        panel.grid.major = element_blank(),                                     # remove the vertical grid lines
        panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
        panel.background = element_rect(fill="black",                     # make the interior background transparent
                                        colour = NA),                           # remove any other colors
        plot.background = element_rect(fill="black",                      # make the outer background transparent
                                       colour=NA),                              # remove any other colors
        axis.line = element_line(colour = "black"),                             # color the x and y axis
        axis.text = element_text(size=20, colour = "black"),                    # color the axis text
        axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
        text=element_text(size=20,                                              # change the size of the axis titles
                          colour = "black"),                                    # change the color of the axis titles
        legend.position = "none") +                                             # remove the legend
  geom_point(aes(fill=Group),
             position = position_jitterdodge(jitter.width = 0, 
                                             dodge.width = .75)) +
  labs(title = NULL, x=NULL, y = "Mean FAC Abundance") +                                               # changing axis titles
  scale_x_discrete(breaks=c("Rest", "Moderate", "Full", "Heavy"), 
                   labels=c("Rest"="Rest", "Moderate"="Moderate", "Full"="Full", "Heavy"="Heavy"),
                   limits=c("Heavy", "Full", "Moderate", "Rest")) +
  coord_cartesian(xlim=NULL, ylim=c(0,50))

fac.plot

totalbirds.bar <- (obl.plot/fac.plot)

totalbirds.bar <- totalbirds.bar + plot_annotation(theme=theme(panel.grid.major = element_blank(), # remove the vertical grid lines
                                                               panel.grid.minor = element_blank(), # remove the horizontal grid lines
                                                               panel.background = element_rect(fill="black", # make the interior background transparent
                                                                                               colour = NA), 
                                                               plot.background = element_rect(fill ="black", # make the outer background transparent
                                                                                              colour = NA), 
                                                               axis.line = element_line(colour = "black"), # color the x and y axis
                                                               axis.text = element_text(colour = "black"), # color the axis text
                                                               axis.ticks = element_line(colour = "black"),
                                                               text=element_text(colour = "black"),# change the color of the axis titles
                                                               legend.position = "none"))

totalbirds.bar

ggsave(birds.bar, 
       filename = "outputs/figs/BirdsBarYear.png",  
       dpi = "print", 
       bg = "black",
       width = 15,
       height = 10)

ggsave(totalbirds.bar, 
       filename = "outputs/figs/BirdsBar.png",  
       dpi = "print", 
       bg = "black",
       width = 15,
       height = 10)
