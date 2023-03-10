# Load libraries ----------------------------------------------------------

library(ggord)
library(tidyverse)
library(vegan)
library(cowplot)
library(indicspecies)

windowsFonts(my_font = windowsFont("Gandhi Sans"))

#This set of code downloads the ggord package
#
#options(repos = c(
#  fawda123 = 'https://fawda123.r-universe.dev',
#  CRAN = 'https://cloud.r-project.org'))
#
#install.packages('ggord')

# Data import -------------------------------------------------------------

totals21 <- read.csv("working/totals21.csv", 
                     row.names = 1)
totals22 <- read.csv("working/totals22.csv", 
                     row.names = 1)
birds21 <- read.csv("working/birds21.csv", 
                    header = TRUE, 
                    row.names = 1)
birds22 <- read.csv("working/birds22.csv", 
                    header = TRUE, 
                    row.names = 1)


# Defining Theme ----------------------------------------------------------

my_theme <-  theme(plot.title = element_text(family="my_font",                             # select the font for the title
                                             size=30,
                                             hjust=.5),
                   panel.grid = element_blank(),                                     # remove the horizontal grid lines
                   panel.background = element_blank(),
                   plot.background = element_blank(),
                   panel.border = element_rect(colour = NA),                               # remove colors
                   axis.line = element_line(colour = "black"),                             # color the x and y axis
                   axis.text.x = element_text(family = "my_font",
                                              size = 20),                                          # color the axis text
                   axis.text.y = element_text(family = "my_font",
                                              size = 20),
                   axis.ticks = element_line(colour = "black"),                            # set axis tick color
                   text=element_text(family = "my_font",
                                     size=24,                                              # change the size of the axis titles
                                     colour = "black"),                                    # change the color of the
                   legend.position="none") 



# Creating 2021 Ordination ------------------------------------------------


birds.b21 <- as.data.frame(colSums(birds21[2:25]>0))                            # create a new data frame that counts the number of occurrences (i.e. columns where it was observed)

birds.b21 <- birds.b21 |>                                                       # select the data frame
  rownames_to_column(var="Spec") |>                                             # convert the built in row names to be a column called Spec
  filter(colSums(birds21[2:25]>0)>1)                                            # retain only rows where the number of occurrence is greater than 1 in the original data frame

# create a new data frame and filter the original data set 
# using the list of species in the second dataset so that 
# only rows with non-singleton species are left in the new data set

birds.c21 <- semi_join(totals21,                                                # select data frame 1
                       birds.b21,                                                 # select data frame 2
                       by = "Spec")                                               # join the two data frames by the shared species column

birds.perm21 <- pivot_wider(birds.c21[c(1:3,5)],                                            # select the data frame
                            names_from = Spec,                                    # select the column names
                            values_from = Abundance,                              # select the values to fill the table with
                            values_fill=0) |>                                     # fill all NA with 0
  column_to_rownames(var="Pasture")                                             # convert the column names to row names

ord.birds21 <- metaMDS(birds.perm21[2:18],                                      # run the ordination
                       distance="bray",                                         # using bray-curtis dissimilarity matrix
                       k=3,                                                     # use 3 dimensions to reduce stress
                       trymax=9999)                                             # run test max of 1000 times to reach a solution

stressplot(ord.birds21)                                                         # plot the ordination to determine fit

# old colors c('#00CCCC', '#F2E30D', '#cc0029', '#000000')

avian.ordination21 <- ggord(ord.birds21,                                        # plot the ordination
                            cols = c('#A2A4A2',                                                       # color for group 1
                                              'lightgoldenrod2',                                               # color for group 2
                                              '#D4A634',                                                       # color for group 3
                                              '#717F5B'),                                                      # color for group 4
                                              grp_in = birds21$cTreat,                                                  # choose the groupings the communities will be placed in
                            ellipse_pro=0.95,                                                         # set confidence intervals to 95%
                            arrow = NULL,                                                             # remove the arrow
                            size=3,                                                                  # change the size of the points
                            alpha_el = .4,                                                            # change the transparency of the ellipses
                            txt=NULL,                                                                 # change the size of the text
                            repel=TRUE,                                                               # prevent overlap of labels
                            xlim= c(-2, 2.5),                                                           # set the limits of the x axis
                            ylim= c(-2, 2),                                                         # set the limits of the y axis
                            grp_title = "Grazing Intensity") +
  geom_text(aes(x = -1.2,
                y = 1.85,
                label = c("Stress = 0.12")),
            family = "my_font",
            size = 5) +
  geom_text(aes(x = -1.2,
                y = 1.7,
                label = c("p = 0.06")),
            family = "my_font",
            size = 5) +
  my_theme +
  labs(title = "2021")

avian.ordination21                                                                # show the ordination output

ggsave(avian.ordination21,                                                        # select the file
       filename = "outputs/figs/Ord21.png",                                       # name the ordination
       dpi = "retina",                                                           # select the print quality (300 dpi)
       bg = "white",                                                      # set the background color to white
       width = 6.5,
       height = 6.5)

birds.PMr21 <- vegdist(birds.perm21[2:18],                                          # create a resemblance matrix from the filtered data matrix
                       method='bray')                                             # create a new resemblance matrix

birds.dSe21 <- betadisper(birds.PMr21,                                              # create an object that tests the dispersion of factors
                          group = birds.perm21$cTreat,                          # test the dispersion of the treatments
                          type = 'centroid',                                      # use centroids of the spatial median
                          bias.adjust = T)                                        # adjust for small sample bias

disp21 <- permutest(birds.dSe21)                                                   # run the test of dispersion

disp21                                                                           # show the results from the dispersion test

perm21 <- adonis2(birds.PMr21 ~ cTreat,                                        # create a permanova object using the resemblance 
                  # and all factors and their interactions -- the * crosses
                  # factors - testing all individual factors and their interactions
                  # if a factor is not fully crossed use a +
                  data = birds.perm21,                                             # pull the factor data from this data set
                  method = "bray",                                               # use bray-curtis distance
                  by = "terms")                                                  # effects ('marginal') or overall significance ('NULL')

perm21                                                                           # show the results of the permanova

# Create 2022 ordination object -------------------------------------------

birds.b22 <- as.data.frame(colSums(birds22[2:28]>0))                            # create a new data frame that counts the number of occurrences (i.e. columns where it was observed)

birds.b22 <- birds.b22 |>                                                       # select the data frame
  rownames_to_column(var="Spec") |>                                             # convert the built in row names to be a column called Spec
  filter(colSums(birds22[2:28]>0)>1)                                            # retain only rows where the number of occurrence is greater than 1 in the original data frame

# create a new data frame and filter the original data set 
# using the list of species in the second dataset so that 
# only rows with non-singleton species are left in the new data set

birds.c22 <- semi_join(totals22,                                                # select data frame 1
                       birds.b22,                                                 # select data frame 2
                       by = "Spec")                                               # join the two data frames by the shared species column

birds.perm22 <- pivot_wider(birds.c22[c(1:3, 5)],                                            # select the data frame
                            names_from = Spec,                                    # select the column names
                            values_from = Abundance,                              # select the values to fill the table with
                            values_fill=0) |>                                     # fill all NA with 0
  column_to_rownames(var="Pasture")                                             # convert the column names to row names

ord.birds22 <- metaMDS(birds.perm22[2:21],                                          # run the ordination
                       distance="bray",                                         # using bray-curtis dissimilarity matrix
                       k=3,                                                     # use 3 dimensions to reduce stress
                       trymax=9999)                                             # run test max of 1000 times to reach a solution

stressplot(ord.birds22)                                                           # plot the ordination to determine fit

# old colors c('#00CCCC', '#F2E30D', '#cc0029', '#000000')

avian.ordination22 <- ggord(ord.birds22,                                        # plot the ordination
                            cols = c('#A2A4A2',                                                       # color for group 1
                                              'lightgoldenrod2',                                               # color for group 2
                                              '#D4A634',                                                       # color for group 3
                                              '#717F5B'),                                                      # color for group 4
                                              grp_in = birds22$cTreat,                                                  # choose the groupings the communities will be placed in
                            ellipse_pro=0.95,                                                         # set confidence intervals to 95%
                            arrow = NULL,                                                             # remove the arrow
                            size=3,                                                                  # change the size of the points
                            alpha_el = .4,                                                            # change the transparency of the ellipses
                            txt=NULL,                                                                 # change the size of the text
                            repel=TRUE,                                                               # prevent overlap of labels
                            xlim= c(-2, 2.5),                                                           # set the limits of the x axis
                            ylim= c(-2, 2),                                                         # set the limits of the y axis
                            grp_title = "Grazing Intensity") +
  geom_text(aes(x = -1.2,
                y = 1.85,
                label = c("Stress = 0.10")),
            family = "my_font",
            size = 5) +
  geom_text(aes(x = -1.2,
                y = 1.7,
                label = c("p = 0.10")),
            family = "my_font",
            size = 5) +
  my_theme +
  labs(title = "2022")

avian.ordination22                                                                # show the ordination output

ggsave(avian.ordination22,
       filename = "outputs/figs/Ord22.png",
       dpi = "retina",
       bg = "white",
       width = 6.5,
       height = 6.5)

birds.PMr22 <- vegdist(birds.perm22[2:21],                                          # create a resemblance matrix from the filtered data matrix
                       method='bray')                                             # create a new resemblance matrix

birds.dSe22 <- betadisper(birds.PMr22,                                              # create an object that tests the dispersion of factors
                          group = birds.perm22$cTreat,                          # test the dispersion of the treatments
                          type = 'centroid',                                      # use centroids of the spatial median
                          bias.adjust = T)                                        # adjust for small sample bias

disp22 <- permutest(birds.dSe22)                                                   # run the test of dispersion

disp22                                                                           # show the results from the dispersion test

perm22 <- adonis2(birds.PMr22 ~ cTreat,                                        # create a permanova object using the resemblance 
                  # and all factors and their interactions -- the * crosses
                  # factors - testing all individual factors and their interactions
                  # if a factor is not fully crossed use a +
                  data = birds.perm22,                                             # pull the factor data from this data set
                  method = "bray",                                               # use bray-curtis distance
                  by = "terms")                                                  # effects ('marginal') or overall significance ('NULL')

perm22                                                                           # show the results of the permanova


# Combining ordinations ---------------------------------------------------


legend <- get_legend(avian.ordination22)

avianord.year <- plot_grid(avian.ordination21 + theme(legend.position = "none"), 
                           avian.ordination22 + theme(legend.position = "none"))

avianord.year

ggsave(avianord.year, 
       filename = "outputs/figs/AvianOrd.png",  
       dpi = "print", 
       bg = "white",
       width = 13,
       height = 6.5)

# Indicator species analysis ----------------------------------------------

ind.spec21 <- multipatt(birds.perm21[2:18],                                         # select the data frame for the analysis (I excluded the row with treatments)
                        birds.perm21$cTreat,                                     # set the groups that will be compared
                        control = how(nperm=1000),                                 # set the number of permutations the test will run through
                        func = "IndVal.g",                                         # set to '.g' to control for unequal sample sizes
                        duleg = T)                                                 # True avoids considering group combos

summary(ind.spec21,                                                               # show summary table
        alpha = 0.05)                                                           # only report species with a significant p-value

ind.spec22 <- multipatt(birds.perm22[2:21],                                         # select the data frame for the analysis (I excluded the row with treatments)
                        birds.perm22$cTreat,                                     # set the groups that will be compared
                        control = how(nperm=1000),                                 # set the number of permutations the test will run through
                        func = "IndVal.g",                                         # set to '.g' to control for unequal sample sizes
                        duleg = T)                                                 # True avoids considering group combos

summary(ind.spec22,                                                               # show summary table
        alpha = 0.05)                                                           # only report species with a significant p-value

# Creating environmental overlay ------------------------------------------

raw <- read.csv("working/RAWadjusted.csv", 
                row.names = 1)

#creating a table with environmental factors
raw21 <- filter(raw, Year == "2021")
raw21$Litter.Depth <- as.numeric(raw21$Litter.Depth)
raw21$Veg.Height <- as.numeric(raw21$Veg.Height)

envfit21 <- group_by(raw21[,c(19:27,33,36:37,39)],                              # select the environmental factors from the raw data
                     Pasture,                                                        # grouped into pastures
                     cTreat,
                     pTreat) |>                                                      # and 2021 grazing intensity
  na.omit() |>                                                                  # omiting any columns without values
  summarise_each(funs(mean)) |>                                                 # and calculating the mean of each variable as it's grouped
  as.data.frame()

birds.vegfit21 <-envfit(ord.birds21 ~                                           # select the ordination to attach it to
                          Litter.Depth +                                          # parameter 1 to include
                          aRobel,                                                 # parameter 2 to include
                        envfit21)                                                 # the file the information is being pulled from

birdfit.df21 <- as.data.frame(birds.vegfit21$vectors$arrows*sqrt(birds.vegfit21$vectors$r)) # turning the envfit information into a dataframe

birdfit.df21$parameters <- c("Litter Depth",                                    # creating a list of parameter 1 and
                             "VOR")                                                  # parameter2 
#creating species names to accompany each variable

# Creating ordinations ----------------------------------------------------

ordination.veg21 <- ggord(ord.birds21,                                              # select the ordination file to plot
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
  theme( plot.title = element_text(family="my_font",                            # select the font for the title
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
         legend.position="none") +                                                # remove the legend
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
  coord_fixed()                                                                 # forces a specified ratio between physical representation of data units on axis

ordination.veg21                                                                  # show the ordination output

ggsave(ordination.veg21,                                                          # select the file
       filename = "outputs/figs/VEGordination21.png",                                          # name the file
       dpi = "print",                                                           # select the print quality (300 dpi)
       bg = "white",                                                      # set the background to transparent
       height = 12.32, 
       width = 15.5)



#creating a table with environmental factors
raw22 <- filter(raw, Year == "2022")
raw22$Litter.Depth <- as.numeric(raw22$Litter.Depth)
raw22$Veg.Height <- as.numeric(raw22$Veg.Height)

envfit22 <- group_by(raw22[,c(19:27,33,36:37,39)],                              # select the environmental factors from the raw data
                     Pasture,                                                        # grouped into pastures
                     cTreat,
                     pTreat) |>                                                      # and 2022 grazing intensity
  na.omit() |>                                                                  # omiting any columns without values
  summarise_each(funs(mean)) |>                                                 # and calculating the mean of each variable as it's grouped
  as.data.frame()

birds.vegfit22 <-envfit(ord.birds22 ~                                           # select the ordination to attach it to
                          Litter.Depth +                                          # parameter 1 to include
                          aRobel,                                                 # parameter 2 to include
                        envfit22)                                                 # the file the information is being pulled from

birdfit.df22 <- as.data.frame(birds.vegfit22$vectors$arrows*sqrt(birds.vegfit22$vectors$r)) # turning the envfit information into a dataframe

birdfit.df22$parameters <- c("Litter Depth",                                      # creating a list of parameter 1 and
                             "VOR")                                                  # parameter2 
#creating species names to accompany each variable

ordination.veg22 <- ggord(ord.birds22,                                              # select the ordination file to plot
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
  theme( plot.title = element_text(family="my_font",                             # select the font for the title
                                   hjust=.5),                                    # center the title
         panel.grid.major = element_blank(),                                     # remove the vertical grid lines
         panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
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
                           size=15,                                              # change the size of the axis labels
                           colour = "black"),                                    # change the color of the axis titles
         aspect.ratio = 0.78,                                                    # this helps when setting the size (divide height by width)
         legend.position="top") +                                                 # remove the legend
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
  coord_fixed()                                                                 # forces a specified ratio between physical representation of data units on axis

ordination.veg22                                                                # show the ordination output

ggsave(ordination.veg22,                                                        # select the file
       filename = "outputs/figs/VEGordination22.png",                                        # name the file
       dpi = "print",                                                           # select the print quality (300 dpi)
       bg = "white",                                                            # set the background to transparent
       height = 12.32, 
       width = 15.5)
```