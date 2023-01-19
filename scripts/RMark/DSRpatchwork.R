
library(ggplot2)
library(cowplot)

windowsFonts(my_font = windowsFont("Gandhi Sans"))

litter.plot <- NOPIlitter.plot
litter.plot

vor.plot <- BRBLvor.plot
vor.plot

kbg.plot <- MODOkbg.plot
kbg.plot

height.plot <- CCSPheight.plot
height.plot

bhco.plot <- CCSPbhco.plot
bhco.plot

bare.plot <- plot_grid(WEMEbare.plot, 
                       RWBLbare.plot, 
                       labels = c("A", "B"))
bare.plot

litdep.plot <- plot_grid(MODOlitdep.plot, 
                         BWTElitdep.plot, 
                         labels = c("A", "B"))
litdep.plot

julian.plot <- plot_grid(BRBLjulian.plot,
                         MODOjulian.plot, 
                         CCSPjulian.plot, 
                         BWTEjulian.plot, 
                         RWBLjulian.plot, 
                         labels = c("AUTO"))
julian.plot

ggsave(litter.plot,
       filename = "~/Git/NDSU/RMARK/Figures/litterplot.png",
       dpi = "print",
       bg = "white")

ggsave(vor.plot,
       filename = "~/Git/NDSU/RMARK/Figures/vorplot.png",
       dpi = "print",
       bg = "white")

ggsave(kbg.plot,
       filename = "~/Git/NDSU/RMARK/Figures/kbgplot.png",
       dpi = "print",
       bg = "white")

ggsave(height.plot,
       filename = "~/Git/NDSU/RMARK/Figures/heightplot.png",
       dpi = "print",
       bg = "white")

ggsave(bhco.plot,
       filename = "~/Git/NDSU/RMARK/Figures/bhcoplot.png",
       dpi = "print",
       bg = "white")

ggsave(bare.plot,
       filename = "~/Git/NDSU/RMARK/Figures/bareplot.png",
       dpi = "print",
       bg = "white")

ggsave(litdep.plot,
       filename = "~/Git/NDSU/RMARK/Figures/litdepplot.png",
       dpi = "print",
       bg = "white")

ggsave(julian.plot,
       filename = "~/Git/NDSU/RMARK/Figures/julianplot.png",
       dpi = "print",
       bg = "white",
       height = 12,
       width=20)