# Loading libraries -------------------------------------------------------

library(vegan)
library(tidyverse)
library(lubridate)
library(plotrix)
library(RMark)

# Data import -------------------------------------------------------------

nest <- read.csv("working/RMarknesting.csv",
                 row.names = 1)

# Creating a density analysis function ------------------------------------

# currently I'm using max number of nests for each treatment, 
#mean age, constant DSR for each treatment, and average patch size.

density <- function(df, spec, model, groups, age) {
  
  beta <- data.frame(model$data$group.covariates,                               # create a dataframe with grouping categories (i.e. treatments)
                     model$results$real)                                        # and DSR estimates
  
  temp.count <- df |>                                                           # create a dataframe
    group_by(Pasture = as.factor(Pasture),                                      # grouped by patch (i.e. each quarter)
             cTreat = as.factor({{groups}}),                                    # and the treatment groups ~ specified in function
             Patch = as.factor(Patch),
             Year = as.factor(Year)) |> 
    summarise(count=n()) |>                                                     # count the number of each species and create a new column
    ungroup()                                                                   # ungroup the list
  
  temp.comb <- full_join(temp.count,                                            # join the count dataframe
                         beta,                                                       # with the DSR dataframe
                         by = c("cTreat", "Year"))                                           # by treatment
  
  temp.comb <- complete(temp.comb[1:9],                                         # this takes the dataframe
                        nesting(Year, cTreat),
                        Patch,                                                # as well as any missing treatments
                        fill = list(count = 0,
                                    estimate = 0,
                                    se = 0,
                                    lcl = 0,
                                    ucl = 0))                                   # and fills all count values as 0
  
  temp.comb <- temp.comb |> 
    mutate(Pasture = case_when(
      Year == 2021 & cTreat == 0 & Patch == "NE" ~ 4,
      Year == 2021 & cTreat == 0 & Patch == "SE" ~ 8,
      Year == 2021 & cTreat == 0 & Patch == "SW" ~ 11,
      Year == 2021 & cTreat == 0 & Patch == "NW" ~ 15,
      
      Year == 2021 & cTreat == 39 & Patch == "NE" ~ 3,
      Year == 2021 & cTreat == 39 & Patch == "SE" ~ 7,
      Year == 2021 & cTreat == 39 & Patch == "SW" ~ 10,
      Year == 2021 & cTreat == 39 & Patch == "NW" ~ 14,
      
      Year == 2021 & cTreat == 49 & Patch == "NE" ~ 2,
      Year == 2021 & cTreat == 49 & Patch == "SE" ~ 6,
      Year == 2021 & cTreat == 49 & Patch == "SW" ~ 9,
      Year == 2021 & cTreat == 49 & Patch == "NW" ~ 13,
      
      Year == 2021 & cTreat == 68 & Patch == "NE" ~ 1,
      Year == 2021 & cTreat == 68 & Patch == "SE" ~ 5,
      Year == 2021 & cTreat == 68 & Patch == "SW" ~ 12,
      Year == 2021 & cTreat == 68 & Patch == "NW" ~ 16,
      
      Year == 2022 & cTreat == 0 & Patch == "NE" ~ 1,
      Year == 2022 & cTreat == 0 & Patch == "SE" ~ 5,
      Year == 2022 & cTreat == 0 & Patch == "SW" ~ 12,
      Year == 2022 & cTreat == 0 & Patch == "NW" ~ 16,
      
      Year == 2022 & cTreat == 39 & Patch == "NE" ~ 4,
      Year == 2022 & cTreat == 39 & Patch == "SE" ~ 8,
      Year == 2022 & cTreat == 39 & Patch == "SW" ~ 11,
      Year == 2022 & cTreat == 39 & Patch == "NW" ~ 15,
      
      Year == 2022 & cTreat == 49 & Patch == "NE" ~ 3,
      Year == 2022 & cTreat == 49 & Patch == "SE" ~ 7,
      Year == 2022 & cTreat == 49 & Patch == "SW" ~ 10,
      Year == 2022 & cTreat == 49 & Patch == "NW" ~ 14,
      
      Year == 2022 & cTreat == 68 & Patch == "NE" ~ 2,
      Year == 2022 & cTreat == 68 & Patch == "SE" ~ 6,
      Year == 2022 & cTreat == 68 & Patch == "SW" ~ 9,
      Year == 2022 & cTreat == 68 & Patch == "NW" ~ 13,
    ))
  
  temp.age <- df |>                                                             # create a data frame
    group_by(Patch = as.factor(Patch),                                          # grouped by patch (i.e. each quarter)
             cTreat = as.factor({{groups}}),                                    # and the treatment groups ~ specified in function
             Year = Year) |>                                 
    summarise(meanAge=mean({{age}})) |>                                         # calculate the mean age for each patch
    ungroup()
  
  temp.age <- complete(temp.age,                                                # this takes the dataframe 
                       nesting(Year, cTreat),
                       Patch,                                                   # as well as any missing treatments
                       fill = list(count = 0,
                                   meanAge = 0))                                # and fills all count values as 0
  
  temp.den <- cbind(temp.comb,                                                  # combine the count/DSR dataframe with the mean age information
                    meanAge=temp.age$meanAge)
  
  Spec.Density <- data.frame(Year = temp.den$Year,
                             Patch = temp.den$Patch,                            # create a species dataframe with patch
                             cTreat = temp.den$cTreat,                          # treatment
                             density = ((temp.den$count/(temp.den$estimate^temp.den$meanAge))/16)) # and the density equation
  
  # this calculates the exact confidence interval ~ better for smaller sample sizes
  #exactPoiCI <- function(x, conf.level) {
  #alpha = 1 - conf.level
  #upper <- 0.5 * qchisq((1-(alpha/2)), (2*(x+1)))
  #lower <- 0.5 * qchisq(alpha/2, (2*x))
  #CI <- data.frame(lower, upper)
  #}
  #CI <- exactPoiCI(Spec.Density$density, 0.95)
  #print(CI)
  
  # this is good for large sample sizes (mine is pretty small)
  spec.graph <- Spec.Density |> 
    group_by(cTreat) |> 
    summarize(density=mean(density),
              se = sqrt(mean(density)/n()),
              ucl = mean(density) + (1.96) * sqrt(mean(density)/n()),
              lcl = mean(density) - (1.96) * sqrt(mean(density)/n()))
  
  density.plot <- ggplot(spec.graph, 
                         aes(x = cTreat, 
                             y = density)) +
    geom_point(aes(x = cTreat,
                   y = density),
               size = 5) +
    geom_errorbar(aes(x = cTreat,
                      ymin = lcl,
                      ymax = ucl),
                  width = .5,
                  size = 1) +
    ggtitle({{spec}}) +
    theme(panel.grid.major = element_blank(),                                     # remove the vertical grid lines
          panel.grid.minor = element_blank(),                                     # remove the horizontal grid lines
          panel.background = element_rect(fill="white",                           # make the interior background transparent
                                          colour = NA),                           # remove any other colors
          plot.background = element_rect(fill="white",                            # make the outer background transparent
                                         colour=NA),                              # remove any other colors
          plot.title = element_text(hjust = 0.5),
          axis.line = element_line(colour = "black"),                             # color the x and y axis
          axis.text = element_text(size=20, colour = "black"),                    # color the axis text
          axis.ticks = element_line(colour = "black"),                            # change the colors of the axis tick marks
          text=element_text(size=20,                                              # change the size of the axis titles
                            colour = "black")) +                                  # change the color of the axis titles
    labs(x = NULL, y = "Nests Per Ha") +
    scale_x_discrete(breaks=c("0", "39", "49", "68"),
                     labels=c("0"="Rest", "39"="Moderate", "49"="Full", "68"="Heavy"),
                     limits=c("68", "49", "39", "0"))  
  print(density.plot)
  print(spec.graph)
  print(Spec.Density)
}

# WEME Density Analysis ---------------------------------------------------


WEME.surv <- filter(nest, Spec=="WEME")

WEME.trt <- mark(WEME.surv, nocc=max(WEME.surv$LastChecked), model = "Nest", groups = c("cTreat", "Year"), model.parameters = list(S = list(formula =  ~1 + Bare)))

WEME.nest <- density(WEME.surv, "WEME", WEME.trt, cTreat, AgeFound)

kruskal.test(density ~ cTreat, data = WEME.nest)
#pairwise.wilcox.test(WEME.nest$density, WEME.nest$cTreat, p.adjust.method = "bonf")

rm(list = ls()[!ls() %in%  c("nest", "density", "WEME.nest")])
cleanup(ask = FALSE)

# BRBL Density Analysis ---------------------------------------------------


BRBL.surv <- filter(nest, Spec=="BRBL")

BRBL.trt <- mark(BRBL.surv, nocc=max(BRBL.surv$LastChecked), model = "Nest", groups = c("cTreat", "Year"), model.parameters = list(S = list(formula =  ~1 + Year + Time + VOR)))

BRBL.nest <- density(BRBL.surv, "BRBL", BRBL.trt, cTreat, AgeFound)

kruskal.test(density ~ cTreat, data = BRBL.nest)
#pairwise.wilcox.test(BRBL.nest$density, BRBL.nest$cTreat, p.adjust.method = "bonf")

rm(list = ls()[!ls() %in%  c("nest", "density", "BRBL.nest")])
cleanup(ask = FALSE)

# RWBL Density Analysis ---------------------------------------------------


RWBL.surv <- filter(nest, Spec=="RWBL")

RWBL.trt <- mark(RWBL.surv, nocc=max(RWBL.surv$LastChecked), model = "Nest",groups = c("cTreat", "Year"), model.parameters = list(S = list(formula =  ~-1 + cTreat + Year + Stage + BHCOpres + KBG + LitterD)))

RWBL.nest <- density(RWBL.surv, "RWBL", RWBL.trt, cTreat, AgeFound)

kruskal.test(density ~ cTreat, data = RWBL.nest)
#pairwise.wilcox.test(RWBL.nest$density, RWBL.nest$cTreat, p.adjust.method = "bonf")

rm(list = ls()[!ls() %in%  c("nest", "density", "RWBL.nest")])
cleanup(ask = FALSE)

# CCSP Density Analysis ---------------------------------------------------


CCSP.surv <- filter(nest, Spec=="CCSP")

CCSP.trt <- mark(CCSP.surv, nocc=max(CCSP.surv$LastChecked), model = "Nest", groups = c("cTreat", "Year"), model.parameters = list(S = list(formula =  ~Year + BHCOpres + AgeFound + Veg.Height + cTreat)))

CCSP.nest <- density(CCSP.surv, 
                     "CCSP", 
                     CCSP.trt, 
                     cTreat, 
                     AgeFound)

kruskal.test(density ~ cTreat, data = CCSP.nest)
#pairwise.wilcox.test(CCSP.nest$density, CCSP.nest$cTreat, p.adjust.method = "bonf")

rm(list = ls()[!ls() %in%  c("nest", "density", "CCSP.nest")])
cleanup(ask = FALSE)

#The following nests did not include year in the top model 
#so I included it in order to group them to create nest densities 
#for each treatment and year.

# NOPI Density Analysis ---------------------------------------------------


NOPI.surv <- filter(nest, Spec=="NOPI")

NOPI.trt <- mark(NOPI.surv, nocc=max(NOPI.surv$LastChecked), model = "Nest", groups = c("cTreat", "Year"), model.parameters = list(S = list(formula =  ~Year + Julian + Julian2 + Litter + cTreat)))

NOPI.nest <- density(NOPI.surv, "NOPI", NOPI.trt, cTreat, AgeFound)

kruskal.test(density ~ cTreat, data = NOPI.nest)
#pairwise.wilcox.test(NOPI.nest$density, NOPI.nest$cTreat, p.adjust.method = "bonf")

rm(list = ls()[!ls() %in%  c("nest", "density", "NOPI.nest")])
cleanup(ask = FALSE)

# BWTE Density Analysis ---------------------------------------------------


BWTE.surv <- filter(nest, Spec=="BWTE")

BWTE.trt <- mark(BWTE.surv, nocc=max(BWTE.surv$LastChecked), model = "Nest", groups = c("cTreat", "Year"), model.parameters = list(S = list(formula =  ~-1 + Year + cTreat + Julian + Julian2 + LitterD)))

BWTE.nest <- density(BWTE.surv, "BWTE", BWTE.trt, cTreat, AgeFound)

kruskal.test(density ~ cTreat, data = BWTE.nest)
#pairwise.wilcox.test(BWTE.nest$density, BWTE.nest$cTreat, p.adjust.method = "bonf")

rm(list = ls()[!ls() %in%  c("nest", "density", "BWTE.nest")])
cleanup(ask = FALSE)

# GADW Density Analysis ---------------------------------------------------


GADW.surv <- filter(nest, Spec=="GADW")

GADW.trt <- mark(GADW.surv, nocc=max(GADW.surv$LastChecked), model = "Nest", groups = c("cTreat", "Year"), model.parameters = list(S = list(formula =  ~ Year + Julian + Julian2 + Forb + Veg.Height + cTreat)))

GADW.nest <- density(GADW.surv, "GADW", GADW.trt, cTreat, AgeFound)

kruskal.test(density ~ cTreat, data = GADW.nest)
#pairwise.wilcox.test(GADW.nest$density, GADW.nest$cTreat, p.adjust.method = "bonf")

rm(list = ls()[!ls() %in%  c("nest", "density", "GADW.nest")])
cleanup(ask = FALSE)

# MODO Density Analysis ---------------------------------------------------


MODO.surv <- filter(nest, Spec=="MODO")

MODO.trt <- mark(MODO.surv, nocc=max(MODO.surv$LastChecked), model = "Nest", groups = "cTreat", model.parameters = list(S = list(formula =  ~Julian + Julian2 + Stage + KBG + LitterD + cTreat)))

MODO.nest <- density(MODO.surv, "MODO", MODO.trt, cTreat, AgeFound)

kruskal.test(density ~ cTreat, data = MODO.nest)
#pairwise.wilcox.test(MODO.nest$density, MODO.nest$cTreat, p.adjust.method = "bonf")

rm(list = ls()[!ls() %in%  c("nest", "density", "MODO.nest")])
cleanup(ask = FALSE)