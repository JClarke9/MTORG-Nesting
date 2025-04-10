# Species Grouping Function

# create a column for each data frame with the associated group information
SpeciesGrouping <- function(x) {
  x <- x |>                                                                     
    mutate(Group = case_when(                                                   
      # Grebes
      endsWith(Spec, "PBGR") ~ "WET",  # Pied-billed Grebe
      
      # Herons & Bitterns
      endsWith(Spec, "AMBI") ~ "FAC",  # American Bittern
      
      # Waterfowl
      endsWith(Spec, "MALL") ~ "FAC",  # Mallard
      endsWith(Spec, "GADW") ~ "FAC",  # Gadwall
      endsWith(Spec, "NOPI") ~ "FAC",  # Northern Pintail
      endsWith(Spec, "AMWI") ~ "FAC",  # American Wigeon
      endsWith(Spec, "NSHO") ~ "FAC",  # Northern Shoveler
      endsWith(Spec, "BWTE") ~ "FAC",  # Blue-winged Teal
      endsWith(Spec, "GWTE") ~ "FAC",  # Green-winged Teal
      endsWith(Spec, "RNDU") ~ "WET",  # Ring-necked Duck
      endsWith(Spec, "LESC") ~ "FAC",  # Lesser Scaup
      
      # Grouse
      endsWith(Spec, "STGR") ~ "OBL",  # Sharp-tailed Grouse
      
      # Rails & Coots
      endsWith(Spec, "AMCO") ~ "WET",  # American Coot
      endsWith(Spec, "SORA") ~ "WET",  # Sora
      
      # Shorebirds
      endsWith(Spec, "KILL") ~ "FAC",  # Killdeer
      endsWith(Spec, "WILL") ~ "FAC",  # Willet
      endsWith(Spec, "UPSA") ~ "OBL",  # Upland Sandpiper
      endsWith(Spec, "WISN") ~ "WET",  # Wilson's Snipe
      endsWith(Spec, "WIPH") ~ "WET",  # Wilson's Phalarope
      
      # Doves
      endsWith(Spec, "MODO") ~ "FAC",  # Mourning Dove
      
      # Kingbirds & Wrens
      endsWith(Spec, "EAKI") ~ "FAC",  # Eastern Kingbird
      endsWith(Spec, "MAWR") ~ "WET",  # Marsh Wren
      
      # Sparrows & Longspurs
      endsWith(Spec, "DICK") ~ "OBL",  # Dickcissel
      endsWith(Spec, "CCSP") ~ "FAC",  # Clay-colored Sparrow
      endsWith(Spec, "GRSP") ~ "OBL",  # Grasshopper Sparrow
      endsWith(Spec, "SAVS") ~ "OBL",  # Savannah Sparrow
      endsWith(Spec, "CCLO") ~ "OBL",  # Chestnut-collared Longspur
      
      # Blackbirds
      endsWith(Spec, "WEME") ~ "OBL",  # Western Meadowlark
      endsWith(Spec, "BOBO") ~ "OBL",  # Bobolink
      endsWith(Spec, "YHBL") ~ "WET",  # Yellow-headed Blackbird
      endsWith(Spec, "RWBL") ~ "FAC",  # Red-winged Blackbird
      endsWith(Spec, "BRBL") ~ "FAC",  # Brewer's Blackbird
      endsWith(Spec, "COGR") ~ "FAC"   # Common Grackle
    ))
}