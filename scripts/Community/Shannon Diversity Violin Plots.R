# Loading libraries -------------------------------------------------------

library(wesanderson)
library(cowplot)
library(tidyverse)
library(extrafont)
library(ggpattern)

windowsFonts(my_font = windowsFont("Gandhi Sans"))                          # downloading a font to be used for my ordination

# Data import -------------------------------------------------------------

tot.rich <- read.csv("working/TotalRichness.csv")
rich21 <- read.csv("working/TotalRichness21.csv")
rich22 <- read.csv("working/TotalRichness22.csv")
rich23 <- read.csv("working/TotalRichness23.csv")

tot.shan <- read.csv("working/ShannonTotal.csv")
shan21 <- read.csv("working/ShannonTotal21.csv")
shan22 <- read.csv("working/ShannonTotal22.csv")
shan23 <- read.csv("working/ShannonTotal23.csv")

tot.rich$Year <- factor(tot.rich$Year,
                        levels = c("2021", "2022", "2023"))
tot.rich$Treat <- factor(tot.rich$Treat,
                         levels = c("Rest", "Moderate", "Full", "Heavy"))

tot.shan$Year <- factor(tot.shan$Year,
                        levels = c("2021", "2022", "2023"))
tot.shan$Treat <- factor(tot.shan$Treat,
                         levels = c("Rest", "Moderate", "Full", "Heavy"))


# Defining themes ---------------------------------------------------------

rich_theme <- theme(plot.title = element_text(family = "my_font",
                                              hjust = .5,
                                              size = 40,
                                              margin = margin(0, 0, 30, 0)), 
                    panel.grid.major = element_blank(),                           # remove the vertical grid lines
                    panel.grid.minor = element_blank(),                           # remove the horizontal grid lines
                    panel.background = element_rect(fill = "transparent",           # make the interior background transparent
                                                    colour = NA),
                    plot.background = element_rect(fill = "transparent",            # make the outer background transparent
                                                   colour = NA),
                    axis.line = element_line(colour = "black"),                   # color the x and y axis
                    axis.text.x = element_blank(),                 # color the axis text
                    axis.title.y = element_text(family = "my_font",
                                                margin = margin(0, 15, 0, 0)),
                    axis.text.y = element_text(colour = "black",
                                               size = 30),
                    axis.ticks = element_line(colour = "black"),
                    text = element_text(family = "my_font",
                                        size = 30,                                    # change the size of the axis titles
                                        colour = "black"),                          # change the color of the axis titles
                    legend.position = "none")

shan_theme <- theme(plot.title = element_text(family = "my_font",
                                              hjust = .5,
                                              size = 40,
                                              margin = margin(0, 0, 30, 0)), 
                    panel.grid.major = element_blank(),                           # remove the vertical grid lines
                    panel.grid.minor = element_blank(),                           # remove the horizontal grid lines
                    panel.background = element_rect(fill = "transparent",           # make the interior background transparent
                                                    colour = NA),
                    plot.background = element_rect(fill = "transparent",            # make the outer background transparent
                                                   colour = NA),
                    axis.line = element_line(colour = "black"),                   # color the x and y axis
                    axis.text.x = element_text(family = "my_font",
                                               size = 30, 
                                               colour = "black"),                 # color the axis text
                    axis.title.y = element_text(family = "my_font",
                                                margin = margin(0, 15, 0, 0)),
                    axis.text.y = element_text(colour = "black",
                                               size = 30),
                    axis.ticks = element_line(colour = "black"),
                    text = element_text(family = "my_font",
                                        size = 30,                                    # change the size of the axis titles
                                        colour = "black"),                          # change the color of the axis titles
                    legend.position = "bottom",
                    legend.background = element_blank(),
                    legend.title = element_text(family = "my_font",
                                                size = 40),
                    legend.text = element_text(family = "my_font",
                                               size = 30),
                    legend.key.width = unit(2, "cm"))


# Creating violin plots ---------------------------------------------------


(rich.viol21 <- ggplot(rich21,
                       aes(fill = Treat,
                           x=Treat, 
                           y=Richness)) +
   geom_violin(trim=FALSE,
               draw_quantiles=c(.25, 
                                .75),
               colour= "black",
               linewidth=.5) +
   scale_fill_manual(values=c('#A2A4A2', 'lightgoldenrod2', 
                              '#D4A634', '#717F5B')) +
   rich_theme +
   labs(title = "2021 Avian Nesting Richness", 
        x = NULL, 
        y = "Species Richness"))

(rich.viol22 <- ggplot(rich22,
                       aes(fill = Treat,
                           x=Treat, 
                           y=Richness)) +
    geom_violin(trim=FALSE,
                draw_quantiles=c(.25, 
                                 .75),
                colour= "black",
                linewidth=.5) +
    scale_fill_manual(values=c('#A2A4A2', 'lightgoldenrod2', 
                               '#D4A634', '#717F5B')) +
    rich_theme +
    labs(title = "2022 Avian Nesting Richness", 
         x = NULL, 
         y = "Species Richness"))

(rich.viol23 <- ggplot(rich23,
                       aes(fill = Treat,
                           x=Treat, 
                           y=Richness)) +
    geom_violin(trim=FALSE,
                draw_quantiles=c(.25, 
                                 .75),
                colour= "black",
                linewidth=.5) +
    scale_fill_manual(values=c('#A2A4A2', 'lightgoldenrod2', 
                               '#D4A634', '#717F5B')) +
    rich_theme +
    labs(title = "2023 Avian Nesting Richness", 
         x = NULL, 
         y = "Species Richness"))

(rich.viol <- ggplot(tot.rich,
                     aes(fill = Treat,
                         x=Year, 
                         y=Richness)) +
    geom_violin(trim=FALSE,
                draw_quantiles=c(.25, 
                                 .75),
                colour= "black",
                linewidth=.5) +
    scale_fill_manual(values=c('#A2A4A2', 'lightgoldenrod2', 
                               '#D4A634', '#717F5B')) +
    rich_theme +
    labs(title = "Avian Nesting Richness", 
         x = NULL, 
         y = "Species Richness"))


# old colors - "#798E87", "#CCC591"

(shan.viol21 <- ggplot(shan21, 
                       aes(fill = Treat,
                           x=Treat, 
                           y=shannon)) +
    geom_violin(trim=FALSE,
                draw_quantiles=c(.25, 
                                 .75),
                colour= "black",
                linewidth=.5) +
    scale_fill_manual(values=c('#A2A4A2', 'lightgoldenrod2', 
                               '#D4A634', '#717F5B')) +
    guides(fill = guide_legend(byrow = TRUE)) +
    shan_theme +
    labs(title = "2021 Avian Nesting Diversity", 
         x = NULL, 
         y = "Shannon's Diversity"))

(shan.viol22 <- ggplot(shan22, 
                       aes(fill = Treat,
                           x=Treat, 
                           y=shannon)) +
    geom_violin(trim=FALSE,
                draw_quantiles=c(.25, 
                                 .75),
                colour= "black",
                linewidth=.5) +
    scale_fill_manual(values=c('#A2A4A2', 'lightgoldenrod2', 
                               '#D4A634', '#717F5B')) +
    guides(fill = guide_legend(byrow = TRUE)) +
    shan_theme +
    labs(title = "2022 Avian Nesting Diversity", 
         x = NULL, 
         y = "Shannon's Diversity"))

(shan.viol23 <- ggplot(shan23, 
                       aes(fill = Treat,
                           x=Treat, 
                           y=shannon)) +
    geom_violin(trim=FALSE,
                draw_quantiles=c(.25, 
                                 .75),
                colour= "black",
                linewidth=.5) +
    scale_fill_manual(values=c('#A2A4A2', 'lightgoldenrod2', 
                               '#D4A634', '#717F5B')) +
    guides(fill = guide_legend(byrow = TRUE)) +
    shan_theme +
    labs(title = "2023 Avian Nesting Diversity", 
         x = NULL, 
         y = "Shannon's Diversity"))


(shan.viol <- ggplot(tot.shan, 
                     aes(fill = Treat,
                         x=Year, 
                         y=shannon)) +
    geom_violin(trim=FALSE,
                draw_quantiles=c(.25, 
                                 .75),
                colour= "black",
                linewidth=.5) +
    scale_fill_manual(values=c('#A2A4A2', 'lightgoldenrod2', 
                               '#D4A634', '#717F5B')) +
    guides(fill = guide_legend(byrow = TRUE)) +
    shan_theme +
    labs(title = "Avian Nesting Diversity", 
         x = NULL, 
         y = "Shannon's Diversity"))


# Combining violin plots --------------------------------------------------

# legend <- get_legend(shan.viol)
#
# (birds.viol <- plot_grid(rich.viol,
#                          legend,
#                          shan.viol + theme(legend.position = "none"),
#                          nrow = 2, 
#                          ncol = 2,
#                          align = "v",
#                          axis = "l",
#                          rel_widths = c(1, .2, 1, .1)))

(birdsY.viol <- plot_grid(rich.viol21 + theme(plot.title = element_blank()), 
                          rich.viol22 + theme(plot.title = element_blank()), 
                          rich.viol23 + theme(plot.title = element_blank()),
                          shan.viol21 + theme(plot.title = element_blank(),
                                              legend.position = "none"), 
                          shan.viol22 + theme(plot.title = element_blank(),
                                              legend.position = "none"), 
                          shan.viol23 + theme(plot.title = element_blank(),
                                              legend.position = "none"),
                          nrow = 2,
                          ncol = 3,
                          alight = "v"))

(birds.viol <- plot_grid(rich.viol,
                         shan.viol,
                         nrow = 2,
                         ncol = 1,
                         align = "v",
                         axis = "l"))


ggsave(birdsY.viol, 
       filename = "outputs/figs/Shannon_Diversity_Violin_Year.png",  
       dpi = 300,
       bg = NULL,
       height = 15,
       width = 12.23,
       units = "in")

ggsave(birds.viol, 
       filename = "outputs/figs/Shannon_Diversity_Violin.png",  
       dpi = 300,
       bg = NULL,
       height = 15,
       width = 12.23,
       units = "in")
