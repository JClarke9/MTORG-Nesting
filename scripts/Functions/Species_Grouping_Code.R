# Species Grouping Function

# create a column for each data frame with the associated group information

SpeciesGrouping <- function(x) {
  x <- x |>                                                                     # select the data frame
    mutate(Group = case_when(                                                   # group together each species and
      endsWith(Spec, "STGR") ~ "OBL",                                           # list each species as FAC or OBL based on Vickery et al. 1999
      endsWith(Spec, "UPSA") ~ "OBL",
      endsWith(Spec, "SAVS") ~ "OBL",
      endsWith(Spec, "GRSP") ~ "OBL",
      endsWith(Spec, "CCLO") ~ "OBL",
      endsWith(Spec, "DICK") ~ "OBL",
      endsWith(Spec, "BOBO") ~ "OBL",
      endsWith(Spec, "WEME") ~ "OBL",
      endsWith(Spec, "GADW") ~ "FAC",
      endsWith(Spec, "AMWI") ~ "FAC",
      endsWith(Spec, "MALL") ~ "FAC",
      endsWith(Spec, "BWTE") ~ "FAC",
      endsWith(Spec, "NSHO") ~ "FAC",
      endsWith(Spec, "NOPI") ~ "FAC",
      endsWith(Spec, "GWTE") ~ "FAC",
      endsWith(Spec, "WILL") ~ "FAC",
      endsWith(Spec, "MODO") ~ "FAC",
      endsWith(Spec, "EAKI") ~ "FAC",
      endsWith(Spec, "CCSP") ~ "FAC",
      endsWith(Spec, "RWBL") ~ "FAC",
      endsWith(Spec, "BRBL") ~ "FAC",
      endsWith(Spec, "KILL") ~ "FAC",
      endsWith(Spec, "COGR") ~ "FAC",
      endsWith(Spec, "WIPH") ~ "WET", 
      endsWith(Spec, "WISN") ~ "WET", 
      endsWith(Spec, "YHBL") ~ "WET",
      endsWith(Spec, "PBGR") ~ "WET",
      endsWith(Spec, "SORA") ~ "WET",
      endsWith(Spec, "AMCO") ~ "WET",
      endsWith(Spec, "MAWR") ~ "WET",
      endsWith(Spec, "AMBI") ~ "FAC",
      endsWith(Spec, "LESC") ~ "FAC"
    ))
}