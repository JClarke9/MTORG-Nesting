# Loading libraries -------------------------------------------------------

library(ggcorrplot)
library(GGally)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(tidyverse)

# Data import -------------------------------------------------------------

nest <- read.csv("working/RMarknesting.csv")

# Subsetting data ---------------------------------------------------------

WEME.surv <- filter(nest, Spec == "WEME")

# Creating a PCA ----------------------------------------------------------

veg <- select(WEME.surv, KBG:VOR)

largeHeight <- WEME.surv |> 
  filter(Veg.Height > quantile(WEME.surv$Veg.Height, 0.90, 
                               na.rm = T)) |> 
  select(id, KBG:VOR)

largeLitter <- WEME.surv |> 
  filter(LitterD > quantile(WEME.surv$LitterD, 0.90, 
                            na.rm = T)) |> 
  select(id, KBG:VOR)

largeVOR <- WEME.surv |> 
  filter(VOR > quantile(WEME.surv$VOR, 0.90, 
                        na.rm = T)) |> 
  select(id, KBG:VOR)

# pcaveg.scaled <- scale(veg, center=TRUE, scale=TRUE)
# pcaveg <- PCA(pcaveg.scaled)
# 
# fviz_screeplot(pcaveg, addlabels=TRUE)
# var.pca <- get_pca_var(pcaveg)
# corrplot(var.pca$cos2)
# corrplot(var.pca$cor)
# 
# var.pca$cos2
# 
# fviz_cos2(pcaveg, choice = "var", axes = 1)
# fviz_cos2(pcaveg, choice = "var", axes = 2)
# 
# ind.pca <- get_pca_ind(pcaveg)
# 
# tall.veg <- ind.pca$coord[,1]
# sparse.grass <- ind.pca$coord[,2]
# 
# WEME.veg <- data.frame(tall.veg, sparse.grass)
# WEME.surv <- cbind(WEME.surv, WEME.veg)
# 
# write.csv(WEME.surv, "~/Git/NDSU/RMARK/Working Data/WEMEsurvival.csv")

# Creating correlation function -------------------------------------------

my_custom_cor <- function(                               # begins creation of a function/graph for the upper half
  # which will contain the coefficient of correlation and an
  # indicator of significane
  data,                                                  # variable 1, analagous to ggplot data argument
  mapping,                                               # variable 2, analagous to ggplot mapping argument
  color = I("grey50"),                                   # variable 3, analagous to ggplot color argument
  sizeRange = c(1, 5),                                   # variable that limits maximum and minimum font size
  ...) {                                                 # ... allows for arguments to be passed to ggplot commands
  # { opens the code that will define the function
  
  # get the x and y data to use the other code
  x <- GGally::eval_data_col(data,                       # pulls the values for x-axis from the data
                             mapping$x)
  y <- GGally::eval_data_col(data,                       # pulls the values for y-axis from the data
                             mapping$y)
  
  ct <- cor.test(x,y)                                    # calculates a vector of p-values and correlations
  sig <- symnum(                                         # creates a dataframe where numbers are associated with signs
    ct$p.value,                                          # column to pull numbers from
    corr = FALSE,                                        # no correlation values
    na = FALSE,                                          # no NAs
    cutpoints = c(0,                                     # sets the intervals of numbers to associated with a sign
                  0.001,                                 # there are 6 cutpoints, resulting in 5 total ranges of
                  0.01,                                  # possible values
                  0.05,
                  0.1,
                  1),
    symbols = c("***",                                   # vector of symbols to be associated with the ranges generated
                "**",                                    # by the cutpoints argument, each symbol is associated with each
                "*",                                     # each range by index number, or position in line
                ".",
                " ")
  )
  
  r <- unname(ct$estimate)                               # create a vector of the coefficient of correlation
  rt <- format(r, digits=2)[1]                           # make the vector coercable
  
  cex <- max(sizeRange)                                  # creates a value to be used as font size later
  
  # plot the cor value
  ggally_text(                                           # begin making the graph
    label = as.character(rt),                            # labels the graph
    mapping = aes(),                                     # allows for mapping to be inherited from other commands
    xP = 0.5, yP = 0.5,                                  # Not sure what these do, most likely affect spacing
    size = 4,                                            # sets the font size of the coefficient of correlation
    color = "black",                                     # sets the text color
    ...                                                  # allows for other ggplot2 and GGally arguments to be used
  ) + 
    geom_text(                                           # adds in the significance indicators
      aes_string(                                        # affects the position of the indicators
        x = 0.8,
        y = 0.8
      ),
      label = sig,                                       # use the text/characters that were created earlier
      size = I(cex),                                     # set the size of the indicators
      color = "darkred",                                 # set the color of the indicators
      ...
    ) + 
    theme_classic() +                                    # sets the theme to classic
    theme(                                               # adjusts the theme
      panel.background = element_rect(                   # changes the background of the graph we are creating
        color = color,                                   # sets the color of the background
        linetype = "longdash"                            # surrounds the plot in a dashed line - can be overwritten later
      ), 
      axis.line = element_blank(),                       # removes unnecessary axis elements
      axis.ticks = element_blank(),                      # removes unnecessary axis elements
      axis.text.y = element_blank(),                     # removes unnecessary axis elements
      axis.text.x = element_blank()                      # removes unnecessary axis elements
    )
}


my_custom_cor(iris, aes(Sepal.Length, Sepal.Width))      # test plot using data built into R

my_custom_smooth <- function(                            # create a function/graph for the bottom half of the matrix
  data,                                                  # argument 1, analagous to data from ggplot2
  mapping,                                               # argument 2, analagous to mapping from ggplot2
  ...) {                                                 # ... allows for other arguments to be passed through/inherited
  ggplot(data = data,                                    # create the example graph
         mapping = mapping) +
    geom_point(color = I("blue"),                        # create the points on the scatterplot - colored blue
               shape = 1,                                # set the shape of the dots, 1 for hollow circles, 21 for filled
               position = "jitter")+                     # jitter the points to eliminate point pooling
    geom_smooth(                                         # add a smoother line/trend line
      method = "lm",                                     # model used to generate line
      color = I("black"), ...)                           # set line color and allow for other arguments
}

my_custom_smooth(iris, aes(Sepal.Length, Sepal.Width))   # test plot using data built into R

# Running WEME correlations -----------------------------------------------

a <- ggpairs(as.data.frame(veg),     # create the matrix and put it in an object -- allows it to
             # be modified by ggplot2 commands later
             upper = list(continuous = my_custom_cor),        # set the upper half to be the figure we created
             lower = list(continuous = my_custom_smooth),     # set the lower halff to be the scatterplot we created
             axisLabels = "none",                              # remove the axis labels
             diag = list(continuous = "barDiag"))              # adds a histogram down the diagonal

b <- a + theme_bw() +                                    # change the theme of the graph to be black and white
  theme(strip.text.y = element_text(size = 5.6,          # change the facet font size on the y-axis
                                    angle = -90,         # rotate it -90 degrees
                                    colour = "black"))+  # set the text color to black
  theme(strip.text.x = element_text(size = 6.5,          # change the facet font size on the x-axis
                                    angle = 0,           # ensure that the text is not rotated
                                    colour = "black"))   # set the text color to blacks

b

ggsave("outputs/figs/WEMEcorrelations.jpg", 
       plot = b, 
       width = 10, 
       height = 10)
