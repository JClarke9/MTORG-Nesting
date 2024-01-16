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

birds21 <- read.csv("working/birds21.csv")
birds22 <- read.csv("working/birds22.csv")
birds23 <- read.csv("working/birds23.csv")


# Data Wrangling ----------------------------------------------------------


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


# Defining Theme ----------------------------------------------------------


my_theme <-  theme(plot.title = element_text(family = "my_font",                             # select the font for the title
                                             size = 40,
                                             hjust = .5),
                   panel.grid = element_blank(),                                     # remove the horizontal grid lines
                   panel.background = element_blank(),
                   plot.background = element_blank(),
                   panel.border = element_rect(colour = NA),                               # remove colors
                   axis.line = element_line(colour = "black"),                             # color the x and y axis
                   axis.text.x = element_text(family = "my_font",
                                              size = 30),                                          # color the axis text
                   axis.text.y = element_text(family = "my_font",
                                              size = 30),
                   axis.ticks = element_line(colour = "black"),                            # set axis tick color
                   text=element_text(family = "my_font",
                                     size = 30,                                              # change the size of the axis titles
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
                  y = 1.90,
                  label = c("Stress = 0.12")),
              family = "my_font",
              size = 9) +
    geom_text(aes(x = -1.2,
                  y = 1.65,
                  label = c("p = 0.03")),
              family = "my_font",
              size = 9) +
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

(disp21 <- permutest(birds21.dSe,
                     permutations = 9999))                                                   # run the test of dispersion

(perm21 <- adonis2(birds21.PMr ~ cTreat,                                        # create a permanova object using the resemblance 
                   # and all factors and their interactions -- the * crosses
                   # factors - testing all individual factors and their interactions
                   # if a factor is not fully crossed use a +
                   data = birds21.perm,                                             # pull the factor data from this data set
                   permutations = 9999,
                   method = "bray",                                               # use bray-curtis distance
                   by = "terms"))                                                  # effects ('marginal') or overall significance ('NULL')

library(emmeans)


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
                  y = 1.90,
                  label = c("Stress = 0.10")),
              family = "my_font",
              size = 9) +
    geom_text(aes(x = -1.2,
                  y = 1.65,
                  label = c("p = 0.11")),
              family = "my_font",
              size = 9) +
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

(disp22 <- permutest(birds22.dSe,
                     permutations = 9999))                                                   # run the test of dispersion

(perm22 <- adonis2(birds22.PMr ~ cTreat,                                        # create a permanova object using the resemblance 
                   # and all factors and their interactions -- the * crosses
                   # factors - testing all individual factors and their interactions
                   # if a factor is not fully crossed use a +
                   data = birds22.perm,                                             # pull the factor data from this data set
                   permutations = 9999,
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
                       k = 3,                                                     # use 3 dimensions to reduce stress
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
                  y = 1.90,
                  label = c("Stress = 0.15")),
              family = "my_font",
              size = 9) +
    geom_text(aes(x = -1.2,
                  y = 1.65,
                  label = c("p = 0.304")),
              family = "my_font",
              size = 9) +
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

(disp23 <- permutest(birds23.dSe,
                     permutations = 9999))                                                   # run the test of dispersion

(perm23 <- adonis2(birds23.PMr ~ cTreat,                                        # create a permanova object using the resemblance 
                   # and all factors and their interactions -- the * crosses
                   # factors - testing all individual factors and their interactions
                   # if a factor is not fully crossed use a +
                   data = birds23.perm,                                             # pull the factor data from this data set
                   permutations = 9999,
                   method = "bray",                                               # use bray-curtis distance
                   by = "terms"))                                                  # effects ('marginal') or overall significance ('NULL')


# Combining ordinations ---------------------------------------------------


(avianord.year <- plot_grid(avian21.ordination + theme(legend.position = "none"), 
                            avian22.ordination + theme(legend.position = "none"),
                            avian23.ordination + theme(legend.position = "none"),
                            nrow = 1,
                            ncol = 3))

ggsave(avianord.year, 
       filename = "outputs/figs/AvianOrd.png",  
       dpi = "print", 
       bg = NULL,
       height = 9,
       width = 22.5)

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
