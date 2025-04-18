# Data import -------------------------------------------------------------

raw <- read.csv("working/Community.csv")

cowbird <- raw |>                                                               # select the data frame
  group_by(id,                                                                  # group the data by Nest.ID
           Spec,                                                                # species
           cTreat) |>                                                           # and the 2021 grazing intensity
  summarize(BHCOpres=max(BHCOpres)) |>                                          # only include the highest value (i.e. 1 or presence) for each species
  as.data.frame()                                                               # coerce the table to a data frame

cowbird <- subset(cowbird, 
                  Spec!="DUCK" & cTreat!= "")

cowbird <- cowbird |>                                                           # select the data frame
  group_by(cTreat,                                                              # group the data by the 2021 grazing intensity
           BHCOpres) |>                                                         # and by whether BHCO eggs were present (1) or absent (0)
  summarize(Abundance=n()) |>                                                   # count all nests for presence/absence
  as.data.frame()                                                               # coerce the table to a data frame

library(ggplot2)

model <- glm(BHCOpres ~ cTreat, data=cowbird, family=binomial())                # test cowbird presence between grazing intensities with a glm
summary(model)                                                                  # show the output of the GLM

library(emmeans)

emmeans(model,                                                                  # select the glm
        pairwise ~ cTreat)                                                      # break the results down to show differences between the 2021 grazing intensities

c.matrix <- cowbird |>                                                          # select the data frame
  pivot_wider(names_from = cTreat,                                              # select the column names
              values_from = Abundance,                                          # select the values to fill the table
              values_fill = 0) |>                                               # fill all NA with 0
  column_to_rownames(var="BHCOpres") |>                                         # convert the BHCOpres column to the row names
  as.matrix()                                                                   # convert the data to a matrix so I can convert it to a table

c.table <- as.table(c.matrix)                                                   # convert the matrix to a table

library(gplots)

balloonplot(t(c.table),                                                         # select the table to graph
            main = "Cowbird Parasitism",                                        # title the graph
            xlab="",                                                            # change the label of the x axis (no label)
            ylab="",                                                            # change the label of the y axis (no label)
            label = TRUE,                                                       # show the numbers in each treatment/group
            show.margins = FALSE)                                               # show the total numbers in each column and treatment

library(graphics)

mosaicplot(c.table,                                                             # select the data frame
           shade=TRUE,                                                          # color the bars based on residuals
           las=1,                                                               # convert the table to vertical (1) or horizontal (2)
           main="Cowbird Parasitism")                                           # title the graph

library(vcd)

assoc(head(c.table,                                                             # select the data frame
           5),                                                                  # this does something with the residuals but I'm not sure what
      shade =TRUE,                                                              # color the bars based on residuals
      las=3)                                                                    # not positive what this does

library(car)

cowbird$BHCOpres <- replace(cowbird$BHCOpres,                                   # select the column to replace values for
                            cowbird$BHCOpres == "1",                            # replace a 1 (presence) with
                            "parasitized")                                      # parasitized

cowbird$BHCOpres <- replace(cowbird$BHCOpres,                                   # select the column to replace values for
                            cowbird$BHCOpres == "0",                            # replace a 0 (absent) with
                            "unparasitized")                                    # unparasitized

write.csv(cowbird, "working/BHCOdata.csv")                                      # write a CSV

c.chi <- chisq.test(c.matrix)                                                   # run a chi squared test on the matrix
c.chi                                                                           # display the results of the chi squared test
