##%######################################################%##
#                                                          #
####           Predicting current and future            ####
####     global distribution of invasive Ligustrum      ####
####          lucidum (W.T. Aiton): assessing           ####
####      emerging risks to biodiversity hotspots       ####
#                                                          #
##%######################################################%##

# Written by: Santiago J. E. Velazco


# Packages
if (!"devtools"%in%installed.packages()){install.packages("devtools")}  
if (!"ENMTML"%in%installed.packages()){devtools::install_github("andrefaa/ENMTML")}  
require(ENMTML)
require(raster)
require(readr)



# We will create M with WWF ecoregion selected by biased uncorrected occ
occ <-
  readr::read_tsv('Occ_Ligustrum_cleanded_2.txt') # Note that this database contain all cleaned records without environmental filtering!
R <-
  raster('./Model/PCA_Current_ClimaSoil_20km/PC1.tif') #Load a raster to be used a base to create another raster with accessible area
Eco <-
  raster::shapefile(
    "C:/Users/santi/OneDrive/Documentos/FORESTAL/NicheModeling/Biomes/WWF/wwf_terr_ecos.shp"
  ) #Olson Terrestrial ecoregions (https://c402277.ssl.cf1.rackcdn.com/publications/15/files/original/official_teow.zip?1349272619)
filt <- extract(Eco, occ[c('x', 'y')])
filt2 <- which(Eco$OBJECTID %in% filt$OBJECTID)
Eco2 <-
  Eco[filt2, ] #Filter from shapefiel those polygon were species have occurrences
plot(Eco2)
points(occ[c('x', 'y')], col = 'red', pch = 19)

Mr <- rasterize(Eco2, R) #Rasterize polygons
Mr[!is.na(Mr[])] <-
  1 #Set value of 1 to cell were will compose calibration area (a.k.a. species accessible area)
dir.create('./Models/CalibrationArea')
writeRaster(Mr,
            './Models/CalibrationArea/Ligustrum_lucidum.tif',
            format = 'GTiff')


# Modeling procedure with ENMTML package using Climate and edaphic variables
?ENMTML::ENMTML

ENMTML::ENMTML(
  pred_dir = './Models/PCA_current_ClimateSoil_20km', # Directory path with predictors
  proj_dir = './Models/PCA_future_ClimateSoil20km/Projection_PCA', # Directory path containing FOLDERS with predictors for different climate scenarios
  result_dir = "Result_30bin_ClimateEdaphic_20km",# Directory path with the folder in which model results will be recorded.
  occ_file = "./Models/Occ_Ligustrum_cleanded_filtered.txt",
  sp = 'sp', #Columns names with species names
  x = 'x', #Columns names with longitude data
  y = 'y', #Columns names with latitude data
  min_occ = 10,
  thin_occ = NULL,
  eval_occ = NULL,
  colin_var = NULL,
  imp_var = FALSE,
  sp_accessible_area=c(method='USER-DEFINED', filepath='./Models/CalibrationArea'), # Method to define calibratio area
  pseudoabs_method =  c(method='GEO_CONST', width=50), #Method to allocate pseudo-absences 
  pres_abs_ratio = 1, 
  part=c(method= 'BLOCK'), #Partition method
  save_part = FALSE,
  save_final = TRUE, 
  algorithm = c('SVM', 'RDF', 'MXD', 'GAU'), #Algorithm used for create models
  thr = c(type='MAX_TSS'), # Threshold use to binarise models
  msdm = NULL,
  ensemble=c(method=c('PCA', 'PCA_SUP','W_MEAN'), metric='TSS'), #Ensemble methods
  extrapolation = TRUE, # Evaluate extrapolation
  cores = 1
)
