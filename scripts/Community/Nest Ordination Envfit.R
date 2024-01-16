# Load libraries ----------------------------------------------------------


library(ggord)
library(tidyverse)
library(vegan)
library(cowplot)

windowsFonts(my_font = windowsFont("Gandhi Sans"))



# Data Import -------------------------------------------------------------


raw <- read.csv("working/Community.csv")

totals21 <- read.csv("working/totals21.csv")
totals22 <- read.csv("working/totals22.csv")
totals23 <- read.csv("working/totals23.csv")

birds21 <- read.csv("working/birds21.csv")
birds22 <- read.csv("working/birds22.csv")
birds23 <- read.csv("working/birds23.csv")


# Data Wrangling ----------------------------------------------------------


raw$Litter.Depth <- as.numeric(raw$Litter.Depth)
raw$Veg.Height <- as.numeric(raw$Veg.Height)
raw$Bare <- as.integer(raw$Bare)


raw21 <- filter(raw, Year == "2021")
raw22 <- filter(raw, Year == "2022")
raw23 <- filter(raw, Year == "2023")

birds21$cTreat <- factor(birds21$cTreat,
                         levels = c("Rest", "Moderate", "Full", "Heavy"))
birds23$cTreat <- factor(birds23$cTreat,
                         levels = c("Rest", "Moderate", "Full", "Heavy"))
birds23$cTreat <- factor(birds23$cTreat,
                         levels = c("Rest", "Moderate", "Full", "Heavy"))
totals21$cTreat <- factor(totals21$cTreat,
                          levels = c("Rest", "Moderate", "Full", "Heavy"))
totals23$cTreat <- factor(totals23$cTreat,
                          levels = c("Rest", "Moderate", "Full", "Heavy"))
totals23$cTreat <- factor(totals23$cTreat,
                          levels = c("Rest", "Moderate", "Full", "Heavy"))


# Creating 2021 Ordination ------------------------------------------------


birds21.b <- as.data.frame(colSums(birds21[2:25]>0))                            # create a new data frame that counts the number of occurrences (i.e. columns where it was observed)

birds21.b <- birds21.b |>                                                       # select the data frame
  rownames_to_column(var = "Spec") |>                                             # convert the built in row names to be a column called Spec
  filter(colSums(birds21[2:25]>0)>1)                                            # retain only rows where the number of occurrence is greater than 1 in the original data frame

# create a new data frame and filter the original data set 
# using the list of species in the second dataset so that 
# only rows with non-singleton species are left in the new data set

birds21.c <- semi_join(totals21,                                                # select data frame 1
                       birds21.b,                                                 # select data frame 2
                       by = "Spec")                                               # join the two data frames by the shared species column

birds21.perm <- pivot_wider(birds21.c[c(1:3,5)],                                            # select the data frame
                            names_from = Spec,                                    # select the column names
                            values_from = Abundance,                              # select the values to fill the table with
                            values_fill = 0) |>                                     # fill all NA with 0
  column_to_rownames(var = "Paddock")                                             # convert the column names to row names

birds21.ord <- metaMDS(birds21.perm[2:18],                                      # run the ordination
                       distance = "bray",                                         # using bray-curtis dissimilarity matrix
                       k = 3,                                                     # use 3 dimensions to reduce stress
                       trymax = 9999)                                             # run test max of 1000 times to reach a solution

stressplot(birds21.ord)                                                         # plot the ordination to determine fit


# Create 2022 ordination object -------------------------------------------


birds22.b <- as.data.frame(colSums(birds22[2:28]>0))                            # create a new data frame that counts the number of occurrences (i.e. columns where it was observed)

birds22.b <- birds22.b |>                                                       # select the data frame
  rownames_to_column(var = "Spec") |>                                             # convert the built in row names to be a column called Spec
  filter(colSums(birds22[2:28]>0)>1)                                            # retain only rows where the number of occurrence is greater than 1 in the original data frame

# create a new data frame and filter the original data set 
# using the list of species in the second dataset so that 
# only rows with non-singleton species are left in the new data set

birds22.c <- semi_join(totals22,                                                # select data frame 1
                       birds22.b,                                                 # select data frame 2
                       by = "Spec")                                               # join the two data frames by the shared species column

birds22.perm <- pivot_wider(birds22.c[c(1:3, 5)],                                            # select the data frame
                            names_from = Spec,                                    # select the column names
                            values_from = Abundance,                              # select the values to fill the table with
                            values_fill=0) |>                                     # fill all NA with 0
  column_to_rownames(var = "Paddock")                                             # convert the column names to row names

birds22.ord <- metaMDS(birds22.perm[2:21],                                          # run the ordination
                       distance = "bray",                                         # using bray-curtis dissimilarity matrix
                       k = 3,                                                     # use 3 dimensions to reduce stress
                       trymax = 9999)                                             # run test max of 1000 times to reach a solution

stressplot(birds22.ord)                                                           # plot the ordination to determine fit


# Create 2023 ordination object -------------------------------------------


birds23.b <- as.data.frame(colSums(birds23[2:29]>0))                            # create a new data frame that counts the number of occurrences (i.e. columns where it was observed)

birds23.b <- birds23.b |>                                                       # select the data frame
  rownames_to_column(var = "Spec") |>                                             # convert the built in row names to be a column called Spec
  filter(colSums(birds23[2:29]>0)>1)                                            # retain only rows where the number of occurrence is greater than 1 in the original data frame

# create a new data frame and filter the original data set 
# using the list of species in the second dataset so that 
# only rows with non-singleton species are left in the new data set

birds23.c <- semi_join(totals23,                                                # select data frame 1
                       birds23.b,                                                 # select data frame 2
                       by = "Spec")                                               # join the two data frames by the shared species column

birds23.perm <- pivot_wider(birds23.c[c(1:3, 5)],                                            # select the data frame
                            names_from = Spec,                                    # select the column names
                            values_from = Abundance,                              # select the values to fill the table with
                            values_fill=0) |>                                     # fill all NA with 0
  column_to_rownames(var = "Paddock")                                             # convert the column names to row names

birds23.ord <- metaMDS(birds23.perm[2:22],                                          # run the ordination
                       distance = "bray",                                         # using bray-curtis dissimilarity matrix
                       k = 3,                                                     # use 3 dimensions to reduce stress
                       trymax = 9999)                                             # run test max of 1000 times to reach a solution

stressplot(birds23.ord)                                                           # plot the ordination to determine fit


# Creating 2021 Envfit Overlay --------------------------------------------


envfit21 <- raw21 |> 
  select(KBG:Veg.Height, VOR, Paddock:pTreat) |>   # Select relevant columns from the data frame, including numeric columns and grouping columns
  group_by(Paddock, cTreat, pTreat) |>                               # Group the data by specified columns
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) |> # Calculate the mean for each numeric column within each group
  na.omit() |> 
  as.data.frame()                                                    # Convert the result to a data frame

birds.vegfit21 <-envfit(birds21.ord ~                                           # select the ordination to attach it to
                          KBG +
                          Litter +
                          Woody +
                          Litter.Depth +
                          Veg.Height +
                          VOR,
                        envfit21)                                                 # the file the information is being pulled from

birdfit.df21 <- as.data.frame(birds.vegfit21$vectors$arrows*sqrt(birds.vegfit21$vectors$r)) # turning the envfit information into a dataframe

birdfit.df21$parameters <- c("KBG",
                             "Litter",
                             "Woody",
                             "Litter.Depth",
                             "VOR")                                                  # parameter2 
#creating species names to accompany each variable


# Creating 2022 Envfit Overlay --------------------------------------------


envfit22 <- raw22 |> 
  select(KBG:Veg.Height, VOR, Paddock:pTreat) |>   # Select relevant columns from the data frame, including numeric columns and grouping columns
  group_by(Paddock, cTreat, pTreat) |>                               # Group the data by specified columns
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) |> # Calculate the mean for each numeric column within each group
  na.omit() |> 
  as.data.frame()                                                     # Convert the result to a data frame

birds.vegfit22 <-envfit(birds22.ord ~                                           # select the ordination to attach it to
                          KBG +
                          Litter +
                          Woody +
                          Litter.Depth +
                          VOR,
                        envfit22)                                                 # the file the information is being pulled from

birdfit.df22 <- as.data.frame(birds.vegfit22$vectors$arrows*sqrt(birds.vegfit22$vectors$r)) # turning the envfit information into a dataframe

birdfit.df22$parameters <- c("KBG",
                             "Litter",
                             "Woody",
                             "Litter.Depth",
                             "VOR")                                                  # parameter2 
#creating species names to accompany each variable


# Creating 2023 Envfit Overlay --------------------------------------------


envfit23 <- raw23 |> 
  select(KBG:Veg.Height, VOR, Paddock:pTreat) |>   # Select relevant columns from the data frame, including numeric columns and grouping columns
  group_by(Paddock, cTreat, pTreat) |>                               # Group the data by specified columns
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) |> # Calculate the mean for each numeric column within each group
  na.omit() |> 
  as.data.frame()                                                    # Convert the result to a data frame

birds.vegfit23 <-envfit(birds23.ord ~                                           # select the ordination to attach it to
                          KBG +
                          Litter +
                          Woody +
                          Litter.Depth +
                          VOR,
                        envfit23)                                                 # the file the information is being pulled from

birdfit.df23 <- as.data.frame(birds.vegfit23$vectors$arrows*sqrt(birds.vegfit23$vectors$r)) # turning the envfit information into a dataframe

birdfit.df23$parameters <- c("KBG",
                             "Litter",
                             "Woody",
                             "Litter.Depth",
                             "VOR")                                                  # parameter2 
#creating species names to accompany each variable


# Defining Theme ----------------------------------------------------------


my_theme <- theme(plot.title = element_text(family="my_font",                            # select the font for the title
                                            hjust=.5),                                   # center the title
                  panel.grid.major = element_blank(),                                    # remove the vertical grid lines
                  panel.grid.minor = element_blank(),                                    # remove the horizontal grid lines
                  panel.background = element_rect(fill="white",                     # make the interior background transparent
                                                  colour = NA),                           # remove other colors
                  plot.background = element_rect(fill="white",                      # make the outer background transparent
                                                 colour=NA),                              # remove other colors
                  panel.border = element_rect(colour = NA),                               # remove other colors
                  axis.line = element_line(colour = "black"),                             # color the x and y axis
                  axis.text = element_text(family="my_font",                              # set the font for the axis
                                           size =15,                                      # change the size of the axis
                                           colour = "black"),                             # color the axis text
                  axis.ticks = element_line(colour = "black"),                            # color the axis ticks
                  text=element_text(family="my_font",                                     # set the font for the axis labels
                                    size=20,                                              # change the size of the axis labels
                                    colour = "black"),                                    # change the color of the axis titles
                  aspect.ratio = 0.78,                                                    # this helps when setting the size (divide height by width)
                  legend.position="none")                                                 # remove the legend


# Creating ordinations ----------------------------------------------------

(ord.veg21 <- ggord(birds21.ord,                                              # select the ordination file to plot
                    cols = c('#A2A4A2',                                                       # choose colors for each of the groups
                             'lightgoldenrod2',                                               # choose colors for each of the groups
                             '#D4A634',                                                       # choose colors for each of the groups
                             '#717F5B'),                                                      # choose colors for each of the groups
                    grp_in = birds21$cTreat,                                                # group the communities based on grazing intensity
                    ellipse=TRUE,                                         # show confidence ellipses for each group
                    ellipse_pro=0.95,                                     # set confidence intervals to 95%
                    arrow = NULL,                                         # remove the arrow
                    size=10,                                              # change the size of the points
                    alpha_el = 0.4,                                       # change the transparency of the ellipses
                    txt=NULL,                                             # change the size of the text
                    repel=TRUE,                                           # prevent overlap of labels
                    xlim= c(-2, 2),                                       # set the limits of the x axis
                    ylim= c(-2, 2.1)) +                                   # set the limits of the y axis
   my_theme +
   geom_segment(data=birdfit.df21,                                                 # add the environmental overlay
                aes(x=0,                                                         # set the beginning x coordinate of the overlay
                    xend=NMDS1,                                                  # set the end x coordinate of the overlay
                    y=0,                                                         # set the beginning y coordinate of the overlay
                    yend=NMDS2),                                                 # set the end y coordinate of the overlay
                size=1,                                                          # changing the size of the arrows
                arrow=arrow(length=unit(.7,                                      # changing the length of the wings of the arrow
                                        "cm")),                                  # setting the unit
                colour="goldenrod4") +                                           # changing the color of the arrow
   geom_text(data=birdfit.df21,                                                    # adding text to the arrows
             aes(x=NMDS1,                                         # defining the x location of the labels
                 y=NMDS2),                                        # defining the y location of the labels
             family="my_font",                                                   # changing the font of the labels
             label=birdfit.df21$parameters,                                        # setting the labels
             size=5,                                                            # changing the size of the labels
             color="black",                                                      # changing the colors of the labels
             nudge_y = c(0,                                                    # changing the vertical location of average robel
                         -.1),                                                   # changing the vertical location of litter depth
             nudge_x = c(.4,                                                    # changing the horizontal location of average robel 
                         0)) +                                                 # changing the horizontal location of litter depth 
   coord_fixed())                                                                 # forces a specified ratio between physical representation of data units on axis


(ord.veg22 <- ggord(birds22.ord,                                              # select the ordination file to plot
                    cols = c('#A2A4A2',                                                       # choose colors for each of the groups
                             'lightgoldenrod2',                                               # choose colors for each of the groups
                             '#D4A634',                                                       # choose colors for each of the groups
                             '#717F5B'),                                                      # choose colors for each of the groups
                    grp_in = birds22$cTreat,                                                # group the communities based on grazing intensity
                    ellipse=TRUE,                                                             # show confidence ellipses for each group
                    ellipse_pro=0.95,                                                         # set confidence intervals to 95%
                    arrow = NULL,                                                             # remove the arrow
                    size=10,                                                                  # change the size of the points
                    alpha_el = 0.4,                                                           # change the transparency of the ellipses
                    txt=NULL,                                                                 # change the size of the text
                    repel=TRUE,                                                               # prevent overlap of labels
                    xlim= c(-2, 2),                                                           # set the limits of the x axis
                    ylim= c(-2, 2.1)) +                                                       # set the limits of the y axis
    my_theme +
    geom_segment(data=birdfit.df22,                                               # add the environmental overlay
                 aes(x=0,                                                         # set the beginning x coordinate of the overlay
                     xend=NMDS1,                                                  # set the end x coordinate of the overlay
                     y=0,                                                         # set the beginning y coordinate of the overlay
                     yend=NMDS2),                                                 # set the end y coordinate of the overlay
                 size=1,                                                          # changing the size of the arrows
                 arrow=arrow(length=unit(.7,                                      # changing the length of the wings of the arrow
                                         "cm")),                                  # setting the unit
                 colour="goldenrod4") +                                           # changing the color of the arrow
    geom_text(data=birdfit.df22,                                                  # adding text to the arrows
              aes(x=NMDS1,                                         # defining the x location of the labels
                  y=NMDS2),                                        # defining the y location of the labels
              family="my_font",                                                   # changing the font of the labels
              label=birdfit.df22$parameters,                                      # setting the labels
              size=5,                                                             # changing the size of the labels
              color="black",                                                      # changing the colors of the labels
              nudge_y = c(0,                                                      # changing the vertical location of average robel
                          0),                                                     # changing the vertical location of litter depth
              nudge_x = c(0,                                                      # changing the horizontal location of average robel 
                          0)) +                                                   # changing the horizontal location of litter depth 
    coord_fixed())                                                                 # forces a specified ratio between physical representation of data units on axis


(ord.veg23 <- ggord(birds23.ord,                                              # select the ordination file to plot
                    cols = c('#A2A4A2',                                                       # choose colors for each of the groups
                             'lightgoldenrod2',                                               # choose colors for each of the groups
                             '#D4A634',                                                       # choose colors for each of the groups
                             '#717F5B'),                                                      # choose colors for each of the groups
                    grp_in = birds23$cTreat,                                                # group the communities based on grazing intensity
                    ellipse=TRUE,                                                             # show confidence ellipses for each group
                    ellipse_pro=0.95,                                                         # set confidence intervals to 95%
                    arrow = NULL,                                                             # remove the arrow
                    size=10,                                                                  # change the size of the points
                    alpha_el = 0.4,                                                           # change the transparency of the ellipses
                    txt=NULL,                                                                 # change the size of the text
                    repel=TRUE,                                                               # prevent overlap of labels
                    xlim= c(-2, 2),                                                           # set the limits of the x axis
                    ylim= c(-2, 2.1)) +                                                       # set the limits of the y axis
    my_theme +
    geom_segment(data=birdfit.df23,                                               # add the environmental overlay
                 aes(x=0,                                                         # set the beginning x coordinate of the overlay
                     xend=NMDS1,                                                  # set the end x coordinate of the overlay
                     y=0,                                                         # set the beginning y coordinate of the overlay
                     yend=NMDS2),                                                 # set the end y coordinate of the overlay
                 size=1,                                                          # changing the size of the arrows
                 arrow=arrow(length=unit(.7,                                      # changing the length of the wings of the arrow
                                         "cm")),                                  # setting the unit
                 colour="goldenrod4") +                                           # changing the color of the arrow
    geom_text(data=birdfit.df23,                                                  # adding text to the arrows
              aes(x=NMDS1,                                         # defining the x location of the labels
                  y=NMDS2),                                        # defining the y location of the labels
              family="my_font",                                                   # changing the font of the labels
              label=birdfit.df23$parameters,                                      # setting the labels
              size=5,                                                             # changing the size of the labels
              color="black",                                                      # changing the colors of the labels
              nudge_y = c(0,                                                      # changing the vertical location of average robel
                          0),                                                     # changing the vertical location of litter depth
              nudge_x = c(0,                                                      # changing the horizontal location of average robel 
                          0)) +                                                   # changing the horizontal location of litter depth 
    coord_fixed())                                                                 # forces a specified ratio between physical representation of data units on axis


# Combining Plots -----------------------------------------------------------

ord.veg21
ord.veg22
ord.veg23

(envfit.plot <- plot_grid(ord.veg21, ord.veg22, ord.veg23,
                          nrow = 1,
                          ncol = 3))
