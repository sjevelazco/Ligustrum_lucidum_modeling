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
require(maps)
require(raster)
require(readr)


##%######################################################%##
#                                                          #
####          PCA with climate & edaphic data           ####
#            for current and future conditions             #
#                                                          #
##%######################################################%##

e <- extent(-180.0001, 179.9999, -55.9000, 83.99986) # Set world extent used to homogenize extent of climate and edaphic layers


# Create direcotries were will be saver PCA outputs
saveDir <-
  './PCA_Current_ClimaSoil_20km' # Folder to save PCA output for current conditions
saveDirF <-
  "./PCA_Future_ClimaSoil20km" # Folder to save PCA output for future conditions
dir.create(saveDir)
dir.create(saveDirF)


# Current conditions
variables1 <-
  raster::stack(list.files(
    './Variables/Chelsa/Current/20km',
    full.names = T,
    pattern = '.tif'
  )) # Climate data

variables2 <-
  raster::stack(
    list.files(
      './Variables/SoilGrids/SoilGrids20km_World_mean',
      full.names = T,
      pattern = '.tif'
    )
  ) # Edaphi data

variables1 <- raster::crop(variables1, e) # crop climate layers
variables2 <- raster::crop(variables2, e) # crop edaphic layers

extent(variables1) # check extent
extent(variables2) # check extent

variables <- raster::stack(variables1, variables2) # stack both databases

# Future
DirP <- "./Variables/Chelsa/Future_climate" # in this folder
DirF <- file.path(DirP, list.files(DirP))
# DirF object is a vector with the paths to folders with different future climate scenarios
# > DirF
# [1] "./Variables/Future_climate/2050_RCP8.5_CCSM4"         
# [2] "./Variables/Future_climate/2050_RCP8.5_HadGEM2-AO"    
# [3] "./Variables/Future_climate/2050_RCP8.5_IPSL-CM5A-LR"  
# [4] "./Variables/Future_climate/2050_RCP8.5_MIROC-ESM-CHEM"
# [5] "./Variables/Future_climate/2080_RCP8.5_CCSM4"         
# [6] "./Variables/Future_climate/2080_RCP8.5_HadGEM2-AO"    
# [7] "./Variables/Future_climate/2080_RCP8.5_IPSL-CM5A-LR"  
# [8] "./Variables/Future_climate/2080_RCP8.5_MIROC-ESM-CHEM"

DirP_PCA <- file.path(saveDirF,'Projection_PCA')
dir.create(DirP_PCA)
FoldersProj <- as.list(file.path(DirP_PCA,basename(DirF)))
lapply(FoldersProj, function(x) dir.create(x)) # Create a sub-folder for each climate scenario


####   Performing PCA   ###
if(class(variables)=="RasterBrick") {
  variables <- stack(variables)
}

df<-rasterToPoints(variables)
dim(df)
df<-na.omit(df)
pca.raw<-df[, -c(1:2)]

# mean and standard deviation (sd) of 
# current conditions used to standardize current and future values (i.e. climate scenarios)
means<-colMeans(pca.raw)
stds<-apply(pca.raw, 2, sd)


# Scaling variables: This procedure will generate a PCA with a CORRELATION MATRIX in order to avoid 
# PCA distortion because different variables scales 
data.scaled <- data.frame(apply(pca.raw, 2, scale))

# Perform the PCA
data.pca <- prcomp(data.scaled, retx = TRUE)
Coef <- data.pca$rotation 

#Variance explained for each PC and their coefficient
varexplained <- round(t(summary(data.pca)$importance), 3)
varexplained <- data.frame(varexplained, round(data.pca$rotation, 3))
readr::write_tsv(
  varexplained,
  file.path(saveDir, 'pca_varexplained_coef.txt')
)

#axis with cumulative variance explanation until 95%
var.96 <- varexplained$Cumulative.Proportion <= 0.95

# Creation of new raster with principal components  
axis<-as.data.frame(data.pca$x)
axis <- round(axis[,var.96], 5)
axis <- data.frame(df[,1:2], axis)
variables3 <- axis
gridded(variables3) <- ~x+y
variables3 <- stack(variables3)

raster::writeRaster(
  variables3,
  filename = file.path(saveDir, names(variables3)),
  format = 'GTiff',
  overwrite = TRUE,
  bylayer = TRUE
)

# Future PCs: here we used the eigenvectors originated of each principal component 
# to calculate the scores for each grid cell for future conditions 
DirP <- as.list(DirF)
names(DirP) <- basename(DirF)

ProjT <-
  lapply(DirP, function(x)
    raster::crop(raster::brick(raster::stack(
      list.files(x, pattern = '.tif', full.names = T)
    )), e))

ProjT <- lapply(ProjT, function(x)
  raster::stack(x, variables2)) # Stack with soils

ProjE <- lapply(ProjT, function(x)
  raster::rasterToPoints(x))

rm(ProjT)

ProjE <- lapply(ProjE, function(x)
  na.omit(x))

ProjER <-
  lapply(ProjE, function(z)
    z[, !(colnames(z) %in% c("x", "y"))])

scale <- lapply(ProjER, function(x) sweep(x, 2, means))
scale <- lapply(scale, function(x) x %*% diag(1 / stds))
PCAFut <- lapply(scale, function(x) x %*% Coef)

PCAFut <-
  lapply(PCAFut, function(x) {
    data.frame(cbind(ProjE[[1]][, (1:2)], x))
  })

PCAFut.95 <- list()
for (j in 1:length(PCAFut)) {
  PCAFut.95[[j]] <- PCAFut[[j]][, c(1, 2, which(var.96) + 2)]
  sp::gridded(PCAFut.95[[j]]) <- ~ x + y
  PCAFut.95[[j]] <- raster::stack(PCAFut.95[[j]])
}

# Save principal components for each future scenario
for (j in 1:length(PCAFut)) {
  raster::writeRaster(
    PCAFut.95[[j]],
    file.path(FoldersProj[[j]], names(PCAFut.95[[j]])),
    bylayer = TRUE,
    format = "GTiff",
    overwrite = TRUE
  )
}
