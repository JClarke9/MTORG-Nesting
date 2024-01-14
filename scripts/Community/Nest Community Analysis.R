# Load libraries ----------------------------------------------------------


library(ggord)
library(tidyverse)
library(vegan)
library(cowplot)
library(indicspecies)

windowsFonts(my_font = windowsFont("Gandhi Sans"))

# # This set of code downloads the ggord package
# 
# options(repos = c(
#  fawda123 = 'https://fawda123.r-universe.dev',
#  CRAN = 'https://cloud.r-project.org'))
# 
# install.packages('ggord')


# Data import -------------------------------------------------------------


totals21 <- read.csv("working/totals21.csv")
totals22 <- read.csv("working/totals22.csv")
totals23 <- read.csv("working/totals23.csv")

birds21 <- read.csv("working/birds21.csv")                             # read in the data set
birds22 <- read.csv("working/birds22.csv")                             # read in the data set
birds23 <- read.csv("working/birds23.csv")                             # read in the data set

windowsFonts(my_font = windowsFont("Gandhi Sans"))


# Defining Theme ----------------------------------------------------------


my_theme <-  theme(plot.title = element_text(family = "my_font",                             # select the font for the title
                                             size = 30,
                                             hjust = .5),
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
                                     size = 24,                                              # change the size of the axis titles
                                     colour = "black"),                                    # change the color of the
                   legend.position="none") 


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

# old colors c('#00CCCC', '#F2E30D', '#cc0029', '#000000')

(avian21.ordination <- ggord(birds21.ord,                                        # plot the ordination
                             cols = c('#A2A4A2',                                                       # color for group 1
                                      'lightgoldenrod2',                                               # color for group 2
                                      '#D4A634',                                                       # color for group 3
                                      '#717F5B'),                                                      # color for group 4
                             grp_in = birds21$cTreat,                                                  # choose the groupings the communities will be placed in
                             ellipse_pro = 0.95,                                                         # set confidence intervals to 95%
                             arrow = NULL,                                                             # remove the arrow
                             size = 3,                                                                  # change the size of the points
                             alpha_el = .4,                                                            # change the transparency of the ellipses
                             txt = NULL,                                                                 # change the size of the text
                             repel = TRUE,                                                               # prevent overlap of labels
                             xlim = c(-2.5, 2.5),                                                           # set the limits of the x axis
                             ylim = c(-2.5, 2.5),                                                         # set the limits of the y axis
                             grp_title = "Grazing Intensity") +
    geom_text(aes(x = -1.2,
                  y = 1.85,
                  label = c("Stress = 0.11")),
              family = "my_font",
              size = 5) +
    geom_text(aes(x = -1.2,
                  y = 1.7,
                  label = c("p = 0.38")),
              family = "my_font",
              size = 5) +
    my_theme +
    labs(title = "2021"))

ggsave(avian21.ordination,                                                        # select the file
       filename = "outputs/figs/Ord21.png",                                       # name the ordination
       dpi = "retina",                                                           # select the print quality (300 dpi)
       bg = "white",                                                      # set the background color to white
       width = 6.5,
       height = 6.5)

birds21.PMr <- vegdist(birds21.perm[2:18],                                          # create a resemblance matrix from the filtered data matrix
                       method = 'bray')                                             # create a new resemblance matrix

birds21.dSe <- betadisper(birds21.PMr,                                              # create an object that tests the dispersion of factors
                          group = birds21.perm$cTreat,                          # test the dispersion of the treatments
                          type = 'centroid',                                      # use centroids of the spatial median
                          bias.adjust = T)                                        # adjust for small sample bias

(disp21 <- permutest(birds21.dSe))                                                   # run the test of dispersion

(perm21 <- adonis2(birds21.PMr ~ cTreat,                                        # create a permanova object using the resemblance 
                   # and all factors and their interactions -- the * crosses
                   # factors - testing all individual factors and their interactions
                   # if a factor is not fully crossed use a +
                   data = birds21.perm,                                             # pull the factor data from this data set
                   method = "bray",                                               # use bray-curtis distance
                   by = "terms"))                                                  # effects ('marginal') or overall significance ('NULL')


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
                       k = 2,                                                     # use 3 dimensions to reduce stress
                       trymax = 9999)                                             # run test max of 1000 times to reach a solution

stressplot(birds22.ord)                                                           # plot the ordination to determine fit

# old colors c('#00CCCC', '#F2E30D', '#cc0029', '#000000')

(avian22.ordination <- ggord(birds22.ord,                                        # plot the ordination
                            cols = c('#A2A4A2',                                                       # color for group 1
                                              'lightgoldenrod2',                                               # color for group 2
                                              '#D4A634',                                                       # color for group 3
                                              '#717F5B'),                                                      # color for group 4
                                              grp_in = birds22$cTreat,                                                  # choose the groupings the communities will be placed in
                            ellipse_pro = 0.95,                                                         # set confidence intervals to 95%
                            arrow = NULL,                                                             # remove the arrow
                            size = 3,                                                                  # change the size of the points
                            alpha_el = .4,                                                            # change the transparency of the ellipses
                            txt = NULL,                                                                 # change the size of the text
                            repel = TRUE,                                                               # prevent overlap of labels
                            xlim = c(-2, 2.5),                                                           # set the limits of the x axis
                            ylim = c(-2, 2),                                                         # set the limits of the y axis
                            grp_title = "Grazing Intensity") +
  geom_text(aes(x = -1.2,
                y = 1.85,
                label = c("Stress = 0.10")),
            family = "my_font",
            size = 5) +
  geom_text(aes(x = -1.2,
                y = 1.7,
                label = c("p = 0.737")),
            family = "my_font",
            size = 5) +
  my_theme +
  labs(title = "2022"))

ggsave(avian22.ordination,
       filename = "outputs/figs/Ord22.png",
       dpi = "retina",
       bg = "white",
       width = 6.5,
       height = 6.5)

birds22.PMr <- vegdist(birds22.perm[2:21],                                          # create a resemblance matrix from the filtered data matrix
                       method = 'bray')                                             # create a new resemblance matrix

birds22.dSe <- betadisper(birds22.PMr,                                              # create an object that tests the dispersion of factors
                          group = birds22.perm$cTreat,                          # test the dispersion of the treatments
                          type = 'centroid',                                      # use centroids of the spatial median
                          bias.adjust = T)                                        # adjust for small sample bias

(disp22 <- permutest(birds22.dSe))                                                   # run the test of dispersion

(perm22 <- adonis2(birds22.PMr ~ cTreat,                                        # create a permanova object using the resemblance 
                   # and all factors and their interactions -- the * crosses
                   # factors - testing all individual factors and their interactions
                   # if a factor is not fully crossed use a +
                   data = birds22.perm,                                             # pull the factor data from this data set
                   method = "bray",                                               # use bray-curtis distance
                   by = "terms"))                                                  # effects ('marginal') or overall significance ('NULL')


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
                       k = 2,                                                     # use 3 dimensions to reduce stress
                       trymax = 9999)                                             # run test max of 1000 times to reach a solution

stressplot(birds23.ord)                                                           # plot the ordination to determine fit

# old colors c('#00CCCC', '#F2E30D', '#cc0029', '#000000')

(avian23.ordination <- ggord(birds23.ord,                                        # plot the ordination
                             cols = c('#A2A4A2',                                                       # color for group 1
                                      'lightgoldenrod2',                                               # color for group 2
                                      '#D4A634',                                                       # color for group 3
                                      '#717F5B'),                                                      # color for group 4
                             grp_in = birds23$cTreat,                                                  # choose the groupings the communities will be placed in
                             ellipse_pro = 0.95,                                                         # set confidence intervals to 95%
                             arrow = NULL,                                                             # remove the arrow
                             size = 3,                                                                  # change the size of the points
                             alpha_el = .4,                                                            # change the transparency of the ellipses
                             txt = NULL,                                                                 # change the size of the text
                             repel = TRUE,                                                               # prevent overlap of labels
                             xlim = c(-2, 2.5),                                                           # set the limits of the x axis
                             ylim = c(-2, 2),                                                         # set the limits of the y axis
                             grp_title = "Grazing Intensity") +
    geom_text(aes(x = -1.2,
                  y = 1.85,
                  label = c("Stress = 0.10")),
              family = "my_font",
              size = 5) +
    geom_text(aes(x = -1.2,
                  y = 1.7,
                  label = c("p = 0.737")),
              family = "my_font",
              size = 5) +
    my_theme +
    labs(title = "2023"))

ggsave(avian23.ordination,
       filename = "outputs/figs/Ord23.png",
       dpi = "retina",
       bg = "white",
       width = 6.5,
       height = 6.5)

birds23.PMr <- vegdist(birds23.perm[2:22],                                          # create a resemblance matrix from the filtered data matrix
                       method = 'bray')                                             # create a new resemblance matrix

birds23.dSe <- betadisper(birds23.PMr,                                              # create an object that tests the dispersion of factors
                          group = birds23.perm$cTreat,                          # test the dispersion of the treatments
                          type = 'centroid',                                      # use centroids of the spatial median
                          bias.adjust = T)                                        # adjust for small sample bias

(disp23 <- permutest(birds23.dSe))                                                   # run the test of dispersion

(perm23 <- adonis2(birds23.PMr ~ cTreat,                                        # create a permanova object using the resemblance 
                   # and all factors and their interactions -- the * crosses
                   # factors - testing all individual factors and their interactions
                   # if a factor is not fully crossed use a +
                   data = birds23.perm,                                             # pull the factor data from this data set
                   method = "bray",                                               # use bray-curtis distance
                   by = "terms"))                                                  # effects ('marginal') or overall significance ('NULL')


# Combining ordinations ---------------------------------------------------


legend <- get_legend(avian23.ordination)

(avianord.year <- plot_grid(avian21.ordination + theme(legend.position = "none") | 
                           avian22.ordination + theme(legend.position = "none") |
                           avian23.ordination + theme(legend.position = "none")))

ggsave(avianord.year, 
       filename = "outputs/figs/AvianOrd.png",  
       dpi = "print", 
       bg = "white",
       width = 13,
       height = 6.5)

# Indicator species analysis ----------------------------------------------

ind.spec21 <- multipatt(birds21.perm[2:18],                                         # select the data frame for the analysis (I excluded the row with treatments)
                        birds21.perm$cTreat,                                     # set the groups that will be compared
                        control = how(nperm=1000),                                 # set the number of permutations the test will run through
                        func = "IndVal.g",                                         # set to '.g' to control for unequal sample sizes
                        duleg = T)                                                 # True avoids considering group combos

summary(ind.spec21,                                                               # show summary table
        alpha = 0.05)                                                           # only report species with a significant p-value


ind.spec22 <- multipatt(birds22.perm[2:21],                                         # select the data frame for the analysis (I excluded the row with treatments)
                        birds22.perm$cTreat,                                     # set the groups that will be compared
                        control = how(nperm=1000),                                 # set the number of permutations the test will run through
                        func = "IndVal.g",                                         # set to '.g' to control for unequal sample sizes
                        duleg = T)                                                 # True avoids considering group combos

summary(ind.spec22,                                                               # show summary table
        alpha = 0.05)                                                           # only report species with a significant p-value


ind.spec23 <- multipatt(birds23.perm[2:22],                                         # select the data frame for the analysis (I excluded the row with treatments)
                        birds23.perm$cTreat,                                     # set the groups that will be compared
                        control = how(nperm=1000),                                 # set the number of permutations the test will run through
                        func = "IndVal.g",                                         # set to '.g' to control for unequal sample sizes
                        duleg = T)                                                 # True avoids considering group combos

summary(ind.spec23,                                                               # show summary table
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
