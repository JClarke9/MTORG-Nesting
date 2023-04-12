README text document for Journal of Applied Ecology article titled:
Finer-scale habitat predicts nest survival in grassland birds more than management and landscape: a multi-scale perspective

Authors:
Justin J. Shew, Clayton K. Nielsen, and Donald W. Sparling

Corresponding author: Justin J. Shew (justin.shew@gmail.com)

Note: Contact Justin regarding question with the analyses. In the provided data files I have not attached coordinates to these study fields or nests, because this nesting study was conducted exclusively on private lands in northwestern, Illinois, USA.  If you are interested in collaborating on research that requires such spatially explicit information, feel free to contact me.

Data files included (“Data” folder):
RWBL.csv = centered red-winged blackbird data file (csv)
DICK.csv = centered dickcissel data file (csv)
ABOVE.csv = centered above-ground nesting data file (csv)
GRND.csv = centered ground-nesting data file (csv)

Analysis code files included (“Code” folder):

Variable definitions from data files and text script:

Data files called in script:
RWBL = centered red-winged blackbird data file (csv)
DICK = centered dickcissel data file (csv)
ABOVE = centered above-ground nesting data file (csv)
GRND = centered ground-nesting data file (csv)

Generic variables of importance:
Date = data of nest check
Ref = reference date used to calculate Julian data in excel
expos = calculated exposure period in days
survive = (1 = nest still active or successful during exposure period; 0 = suspected nest depredation during exposure period)
trials = 1 for all data (survive/trials was used as the response; thus, 0 or 1 for the exposure period) (note: this way approaching the response was unnecessary)

Note: below follows the analysis steps taken in the manuscript. Please see article text and Tables 1 and 2 for further clarification and details.

Random effect variables called in script (Step 1)
Nest.ID = individual nest identifier
Field = study field identifier
Landscape = field cluster identifier
Year = study year
Spec = species American Ornithologist Union alpha code

Temporal data variables called in script (Step 2)
Julian = standardized Julian date of nest visit
Egg = dominant stage of nest for the exposure period (EGG or NTLG; see article text for description)

Nest-site data variables called in script (Step 3a)
Nht = Nest.Ht. in article (nest height)
Forb = Forb% in article (% forb cover at nest)
Robel = Robel in article (vegetation density surrounding nest)
DeadVeg = Dead.Veg% in article (% dead vegetation cover at nest)
Dist.Edge = Dist.Edge in article (nest distance to nearest non-grassland edge)

Microhabitat-scale variables called in script (Step 3b)
MeanPC1 = MeanPC1 in article
MeanPC2 = MeanPC2 in article
SDPC1 = SDPC1 in article

Patch-scale variable called in script (Step 3c)
Area = AREA in article (area of field)
For.Edge = FOR.EDGE in article (% of field edge with forest)
Ag.Edge = AG.EDGE in article (% of field edge with agriculture)
Grass.100 = GRASS.100 in article (% of grassland within a 100-m buffer of field outline)
PeriAreaR = PAR in article (perimeter-to-area ratio)

Landscape-scale variable called in script (Step 3d)
For.1600 = FOR.1600 in article (% of forest cover within 1600-m radius circle centered on field)
Grass.1600 = GRASS.1600 in article (% of grassland cover within 1600-m radius circle centered on field)
TE.1600 = EDGE.1600 in article (Total amount of edge within 1600-m radius circle centered on field)

Multi-scale models (Step 4)
Note: same variables as above

Management variables called in script (Step 5)
Dist = MGT.T in article (management combination factor)
Y.Man = %TREAT in article (yearly % amount of management on field)
C.Man = %TOTAL.TR in article (Cumulative % of field managed during study)
GRASS = GRASS in article (a native- or brome-grass dominated field)
MCM = MCM in article (field assigned mid-contract management treatment)

Steps 6 and 7 use variables used previously




