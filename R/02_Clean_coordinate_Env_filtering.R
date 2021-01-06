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
require(ape)
require(CoordinateCleaner)
require(maps)
require(dplyr)
require(raster)
require(readr)

##%######################################################%##
#                                                          #
####                Coordinate cleaning                 ####
#                                                          #
##%######################################################%##
occ <- readr::read_tsv("Occ_Ligustrum.txt")

occ2 <- clean_coordinates(
  occ,
  lon = "x",
  lat = "y",
  species = "sp",
  countries = 'country',
  tests = c("capitals", "centroids", "institutions"),
  capitals_rad = 1000,
  centroids_rad = 10000,
  centroids_detail = "both",
  inst_rad = 1000,
  range_rad = 0,
  zeros_rad = 0.5,
  capitals_ref = NULL,
  centroids_ref = NULL,
  country_ref = NULL,
  inst_ref = NULL,
  range_ref = NULL,
  seas_ref = 10,
  seas_scale = 10,
  urban_ref = NULL,
  value = "spatialvalid",
  verbose = TRUE,
  report = FALSE
)

occ2 <- dplyr::as_tibble(occ2)
table(occ2$.summary)

# Detect duplicated coordinates
occ2$.summary[duplicated(occ[,c('x', 'y')])] <- FALSE
occ2[occ2$.summary, c('x', 'y')]
plot(occ2[!occ2$.summary, c('x', 'y')], col='red', pch=19)
points(occ2[occ2$.summary, c('x', 'y')], col='blue', pch=19)
maps::map(add = T)

# Filter database by those records not flagged
occ3 <- occ2[occ2$.summary, c('x', 'y')]

# Save occurrences records in shapefile format to check all location on GIS
occ3shp <- sp::SpatialPointsDataFrame(coords = occ3[c('x', 'y')], data = occ3)
shapefile(x = occ3shp, "Occ_Ligustrum_cleanded.shp")


##%######################################################%##
#                                                          #
####              Sampling bias correction              ####
#                                                          #
##%######################################################%##
# Load function to perform occurrences environmental filtering
source('https://raw.githubusercontent.com/sjevelazco/Ligustrum_lucidum_modeling/main/R/env_filtering.R') 
args(env_filtering)
# coord: A data.frame with longitude and latitude in the first and second columns.
# variables: A data.frame with environmental conditions. It is possible use two or three variables.
# nbins: A number of classes used to split each environmental condition.
# plot: Plot filtering procedure.

# List Principal Component layer performed with climate and edaphic variables
list_var <- list.files('./PCA_Current_ClimaSoil_20km/',
                       full.names = TRUE,
                       pattern = '.tif')

# Will be used the first three variables (PC1, PC2, and PC3) for occurrence environmental filtering 
var <- raster::stack(list_var[1:3])

# Read final occurrence records
occ <- readr::read_tsv("./GitHub/Occ_Ligustrum_cleanded_2.txt") # This is the final database after the manual revision of the database in a GIS
names(occ)
# [1] "sp"      "country" "year"    "y"       "x"       "ID"      "comment"


# Eliminate records that fall in grid cells with NA values 
filtr <- data.frame(raster::extract(var, occ[,c('x','y')])) %>% 
  complete.cases()
filtr %>% table
occ <- occ %>% dplyr::filter(filtr)


# FILTERING PROCEDURE
# Filtering will be performed using the three first principal components (PC) for current conditions (which explain more than 70% of the original environmental data variability)

# Create a new database with coordinates and PC data
coord_env <- occ  %>% dplyr::select(x, y)  
data <- coord_env %>% extract(var, .) %>% data.frame()
coord_env <- coord_env %>% dplyr::mutate(data) %>% data.frame()
rm(data)
coord_env %>% head

bins <- c(50, 40, 30, 20, 10) #Set number of bins to be tested
filtered <- list() # Create a list to allocate outputs 

for (i in 1:length(bins)) {
  message('Filtering ', i)
  filtered[[i]] <- env_filtering(
    coord = coord_env %>% dplyr::select(x, y), # select columns with coordinates
    variables = coord_env %>% dplyr::select(PC1, PC2, PC3), # select columns with PC
    nbins = bins[i],
    do.plot = TRUE
  )
}

# How many records were left at the end of each filtering?
sapply(filtered, nrow)

names(filtered) <- bins

filtered2 <- filtered

for(i in 1:length(filtered)){
  filtered2[[i]] <- occ %>% 
    dplyr::filter(x%in%filtered[[i]]$x & y%in%filtered[[i]]$y)
}

filtered2


# Autocorrelation analysis and number of records for each bin
imoran <- list()
for(i in 1:length(filtered2)){
  message("Spatial autocorrelation dataset ", i)
  coord <- filtered2[[i]] %>% dplyr::select(x, y)
  data <- data.frame(extract(var, coord))
  distm <- dist(coord)
  distm <- as.matrix(distm)
  distm <- 1/distm
  diag(distm) <- 0
  try(imoran[[i]] <-
        apply(data, 2, function(x)
          ape::Moran.I(x, distm, na.rm = T)[c(1, 4)] %>% unlist) %>% data.frame() %>% as_tibble()
  )
  try(imoran[[i]] <- imoran[[i]][1,])
}

names(imoran) <- bins
imorandf <- dplyr::bind_rows(imoran, .id='bins')
imorandf <- imorandf %>% dplyr::mutate(mean_bin=apply(imorandf[, c('PC1', 'PC2', 'PC3')], 1, mean))
imorandf$nrecords <- sapply(filtered2, nrow)
imorandf

# Final filter
finalfilter <- imorandf %>% 
  dplyr::filter(mean_bin<=quantile(mean_bin)[2]) %>% # Bins set with lower spatial autocorrelation
  dplyr::filter(nrecords==max(nrecords)) %>% # Bins set with higher number of records
  pull(bins)
finalfilter

# Best filtered occurrence set 
finaldatabase <- filtered2[[finalfilter]]
readr::write_tsv(finaldatabase, "Occ_Ligustrum_cleanded_filtered.txt")
readr::write_tsv(imorandf, "Filtering_info.txt")