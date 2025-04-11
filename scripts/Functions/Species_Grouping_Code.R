# Species Grouping Function

# create a column for each data frame with the associated group information
SpeciesGrouping <- function(x, colname) {
  x <- x |>
    mutate(Group = case_when(
      # Grebes
      endsWith(as.character(x[[colname]]), "PBGR") ~ "WET",  # Pied-billed Grebe
      
      # Herons & Bitterns
      endsWith(as.character(x[[colname]]), "AMBI") ~ "FAC",  # American Bittern
      
      # Waterfowl
      endsWith(as.character(x[[colname]]), "MALL") ~ "FAC",  # Mallard
      endsWith(as.character(x[[colname]]), "GADW") ~ "FAC",  # Gadwall
      endsWith(as.character(x[[colname]]), "NOPI") ~ "FAC",  # Northern Pintail
      endsWith(as.character(x[[colname]]), "AMWI") ~ "FAC",  # American Wigeon
      endsWith(as.character(x[[colname]]), "NSHO") ~ "FAC",  # Northern Shoveler
      endsWith(as.character(x[[colname]]), "BWTE") ~ "FAC",  # Blue-winged Teal
      endsWith(as.character(x[[colname]]), "GWTE") ~ "FAC",  # Green-winged Teal
      endsWith(as.character(x[[colname]]), "RNDU") ~ "WET",  # Ring-necked Duck
      endsWith(as.character(x[[colname]]), "LESC") ~ "FAC",  # Lesser Scaup
      
      # Grouse
      endsWith(as.character(x[[colname]]), "STGR") ~ "OBL",  # Sharp-tailed Grouse
      
      # Rails & Coots
      endsWith(as.character(x[[colname]]), "AMCO") ~ "WET",  # American Coot
      endsWith(as.character(x[[colname]]), "SORA") ~ "WET",  # Sora
      
      # Shorebirds
      endsWith(as.character(x[[colname]]), "KILL") ~ "FAC",  # Killdeer
      endsWith(as.character(x[[colname]]), "WILL") ~ "FAC",  # Willet
      endsWith(as.character(x[[colname]]), "UPSA") ~ "OBL",  # Upland Sandpiper
      endsWith(as.character(x[[colname]]), "WISN") ~ "WET",  # Wilson's Snipe
      endsWith(as.character(x[[colname]]), "WIPH") ~ "WET",  # Wilson's Phalarope
      
      # Doves
      endsWith(as.character(x[[colname]]), "MODO") ~ "FAC",  # Mourning Dove
      
      # Kingbirds & Wrens
      endsWith(as.character(x[[colname]]), "EAKI") ~ "FAC",  # Eastern Kingbird
      endsWith(as.character(x[[colname]]), "MAWR") ~ "WET",  # Marsh Wren
      
      # Sparrows & Longspurs
      endsWith(as.character(x[[colname]]), "DICK") ~ "OBL",  # Dickcissel
      endsWith(as.character(x[[colname]]), "CCSP") ~ "FAC",  # Clay-colored Sparrow
      endsWith(as.character(x[[colname]]), "GRSP") ~ "OBL",  # Grasshopper Sparrow
      endsWith(as.character(x[[colname]]), "SAVS") ~ "OBL",  # Savannah Sparrow
      endsWith(as.character(x[[colname]]), "CCLO") ~ "OBL",  # Chestnut-collared Longspur
      
      # Blackbirds
      endsWith(as.character(x[[colname]]), "WEME") ~ "OBL",  # Western Meadowlark
      endsWith(as.character(x[[colname]]), "BOBO") ~ "OBL",  # Bobolink
      endsWith(as.character(x[[colname]]), "YHBL") ~ "WET",  # Yellow-headed Blackbird
      endsWith(as.character(x[[colname]]), "RWBL") ~ "FAC",  # Red-winged Blackbird
      endsWith(as.character(x[[colname]]), "BRBL") ~ "FAC",  # Brewer's Blackbird
      endsWith(as.character(x[[colname]]), "COGR") ~ "FAC"   # Common Grackle
    ))
}