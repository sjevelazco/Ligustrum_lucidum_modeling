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

# Modeling procedur with ENMTML package using Climate and edaphic variables
?ENMTML::ENMTML

ENMTML::ENMTML(
  pred_dir = './Models/PCA_current_ClimateSoil_20km', # Directory path with predictors
  proj_dir = './Models/PCA_future_ClimateSoil20km/Projection_PCA', # Directory path containing FOLDERS with predictors for different climate scenarios
  result_dir = "Result_30bin_ClimateEdaphic_20km",# Directory path with the folder in which model results will be recorded.
  occ_file = "./Models/occ_cleaned_2019_FiltClim_30bin.txt",
  sp = 'sp',
  x = 'x',
  y = 'y',
  min_occ = 10,
  thin_occ = NULL,
  eval_occ = NULL,
  colin_var = NULL,
  imp_var = FALSE,
  sp_accessible_area=c(method='USER-DEFINED', filepath='./Models/CalibrationArea'),
  pseudoabs_method =  c(method='GEO_CONST', width=50),
  pres_abs_ratio = 1,
  part=c(method= 'BLOCK'),
  save_part = FALSE,
  save_final = TRUE,
  algorithm = c('SVM', 'RDF', 'MXD', 'GAU'),
  thr = c(type='MAX_TSS'),
  msdm = NULL,
  ensemble=c(method=c('PCA', 'PCA_SUP','W_MEAN'), metric='TSS'),
  extrapolation = TRUE,
  cores = 3
)
