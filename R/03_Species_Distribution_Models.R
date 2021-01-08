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
require(dplyr)
require(readr)
require(ggplot2)
require(ggrepel)
require(reshape2)


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






##%######################################################%##
#                                                          #
####        MODELING POST PROCESSING AND FIGURES        ####
#                                                          #
##%######################################################%##
source('https://raw.githubusercontent.com/sjevelazco/Ligustrum_lucidum_modeling/main/R/additional_functions.R') 
library(scales)
library(cowplot)

dir.create('./Figures')

# Read models performace txt
perf <-
  './Result_30bin_ClimateEdaphic_20km/Evaluation_Table.txt' %>% read_tsv()
perf <-
  perf %>% dplyr::select(-ends_with('_SD'), -Partition, -Sp) %>% dplyr::filter(Algorithm ==
                                                                                 'PCA')

# Read threshold txt
thr <-
  './Result_30bin_ClimateEdaphic_20km/Thresholds_Ensemble.txt' %>% read_tsv()
thr <- thr %>% dplyr::filter(Algorithm == 'PCA')
thr


# Read rasters for current conditions
ClimaSoil <-
  './Result_30bin_ClimateEdaphic_20km/Ensemble' %>% list.dirs()
ClimaSoil <-
  list.files(ClimaSoil, full.names = T) %>% grep("PCA/", ., value = T)
ens1 <- raster(ClimaSoil[1])
ens2 <-
  ens1 >= thr$THR_VALUE[5] #threshold that maximize Jaccard metric

ext_current <-
  raster('./Result_30bin_ClimateEdaphic_20km/Extrapolation/Ligustrum_lucidum_MOP.tif')

ens2[which((ext_current < 0.91)[])] <-
  0 # threshold selected for extrapolation
plot(ens2)
points(occ[, c('x', 'y')], pch = 19, cex = 0.5)
writeRaster(ens2, "./Figures/CurrentDistributionBin", format = 'GTiff')


# Plot for current conditions - continuous model
p <- plot_cont(ens1, '') + scale_fill_viridis()
ggsave(
  plot = p,
  './Figures/Suitability_Current_Continuous_2020.png',
  dpi = 1000,
  units = "cm",
  width = 15,
  height = 6.5,
  scale = 1.8
)

# Plot for current conditions - binary model
show_col(viridis_pal(option = 'B')(15))
viridis_pal(option = 'B')(15)

pb <-
  plot_cont(ens2, '', cont = F) + scale_fill_manual(values = c("gray20", "#FB9E07FF")) +
  theme(legend.position = 'none')
pb

ggsave(
  plot = ggarrange(p, pb, ncol = 1),
  './Figures/Suitability_Current_ConBin_2020.png',
  dpi = 1000,
  units = "cm",
  width = 15,
  height = 6.5 * 2,
  scale = 1.8
)


# Plot current extrapolation maps
p <- plot_cont(ext_current, '') + scale_fill_viridis(option = 'B')
ggsave(
  plot = p,
  './Figures/Extrapolation_Current_Binary_2020.png',
  dpi = 1000,
  units = "cm",
  width = 15,
  height = 6.5,
  scale = 1.8
)



# Read raster of models for futures Projections
ClimaSoilF <-
  './Result_30bin_ClimateEdaphic_20km/FutureLayers' %>% list.files(full.names = T)
ClimaSoilF <-
  ClimaSoilF %>% grep("PCA_2", ., value = T) %>% grep("mean_con", ., value = T)
ens1f <- stack(ClimaSoilF)
plot(ens1f)
plot(ens1f >= thr$THR_VALUE[5])
points(occ[, c('x', 'y')], pch = 19, cex = 0.5)
ens2f <- ens1f >= thr$THR_VALUE[5] #threshold that maximise Jaccard
plot(ens2f + ens2)


# Process extrapolation raster
extr <-
  './Result_30bin_ClimateEdaphic_20km/FutureLayers' %>% list.files(full.names = T) %>% grep('Extrapolation', ., value = T) %>% stack
extrc <- extr
extr$Extrapolation_2050 <- extr$Extrapolation_2050 < 0.92
extr$Extrapolation_2080 <- extr$Extrapolation_2080 < 0.94
plot(extr)
for (i in 1:2) {
  ens2f[[i]][extr[[i]] == 1] <- 0
}
plot(ens2f %>% sum)
plot(ens2f + ens2 * 2, col = sort(pals::tol(4)))
plot(ens2f)
# Save raster of extrapolation
writeRaster(ens2f$PCA_2050_mean_con,
            './Figures/FutureDistributionBin2080',
            format = 'GTiff')
writeRaster(ens2f$PCA_2080_mean_con,
            './Figures/FutureDistributionBin2050',
            format = 'GTiff')

# Calculate future uncertainty metrics based on standard deviation
FutureUnc <-
  './Result_30bin_ClimateEdaphic_20km/FutureLayers' %>% list.files(full.names = T)
unc <-
  FutureUnc %>% grep("sd", ., value = T) %>% grep("PCA_2", ., value = T) %>% stack
plot(unc, col = pals::viridis(12))

ens1f # Continuous predictions
ens2f # Binary predictions
ens3f <-
  ens2f + (ens2 * 2) # Changes between current and future conditions
unc # Uncertainty


# Figures
# Future continuous
p1 <- plot_cont(ens1f[[1]], '') + scale_fill_viridis()
p2 <- plot_cont(ens1f[[2]], '') + scale_fill_viridis()

ggsave(
  plot = ggarrange(p2, p1, ncol = 1),
  './Figures/Suitability_Future_Continuous_2020.png',
  dpi = 1000,
  units = "cm",
  width = 15,
  height = (6.5 * 2),
  scale = 1.8
)


# Future Binary
fcolor <- c("gray20", '#2E6E8EFF', '#D44842FF', "#FB9E07FF")
plot(ens3f$layer.1, col = fcolor)
p1 <-
  plot_cont(ens3f$layer.1, '', cont = F) + scale_fill_manual(values = fcolor) +
  theme(legend.position = 'none')
p2 <-
  plot_cont(ens3f$layer.2, '', cont = F) + scale_fill_manual(values = fcolor) +
  theme(legend.position = 'none')

ggsave(
  plot = ggarrange(p2, p1, ncol = 1),
  './Figures/Suitability_Future_Change_2020.png',
  dpi = 1000,
  units = "cm",
  width = 15,
  height = 6.5 * 2,
  scale = 1.8
)

# Future uncertainty
unc %>% plot
p1 <-
  plot_cont(unc[[1]], '') + scale_fill_viridis(option = 'B', direction = -1)
p2 <-
  plot_cont(unc[[2]], '') + scale_fill_viridis(option = 'B', direction = -1)

ggsave(
  plot = ggarrange(p2, p1, ncol = 1),
  './Figures/Uncerainity_Future_2020.png',
  dpi = 1000,
  units = "cm",
  width = 15,
  height = 6.5 * 2,
  scale = 1.8
)

# Future ucertainty
unc %>% plot
p1 <-
  plot_cont(extrc[[1]], '') + scale_fill_viridis(option = 'B', direction = 1)
p2 <-
  plot_cont(extrc[[2]], '') + scale_fill_viridis(option = 'B', direction = 1)
ggsave(
  plot = ggarrange(p2, p1, ncol = 1),
  './Figures/Extrapolation_Future_2020.png',
  dpi = 1000,
  units = "cm",
  width = 15,
  height = 6.5 * 2,
  scale = 1.8
)

# Calibration area used to construct species distribuion models
ClimaSoil <-
  './Result_30bin_ClimateEdaphic_20km/Ensemble' %>% list.dirs()
ClimaSoil <-
  list.files(ClimaSoil, full.names = T) %>% grep("PCA/", ., value = T)
ens1 <- raster(ClimaSoil[1])
ens1 <- ens1 >= 0
M <-
  raster('./Model/Result_30bin_ClimateEdaphic_20km/Extent_Masks/Ligustrum_lucidum.tif')
ens1[!is.na(M[])] <- 2
plot(ens1)
p1 <-
  plot_cont(ens1, '', cont = F) + scale_fill_manual(values = fcolor[c(1, 3)])
ggsave(
  plot = p1,
  './Figures/M.png',
  dpi = 400,
  units = "cm",
  width = 15,
  height = 6.5,
  scale = 1.8
)

##%######################################################%##
#                                                          #
####              Potential invasive area               ####
#                                                          #
##%######################################################%##
setwd("./Model")

world <- shapefile('./country_world.shp')
world@data$name %>% unique %>% sort
r <-
  list.files('./Figures', pattern = '.tif', full.names = T) %>% stack
r[[1]][r[[1]][] != 1] <- NA
r[[2]][r[[2]][] != 1] <- NA
r[[3]][r[[3]][] != 1] <- NA

plot(r)

a <- area(r, na.rm = T)
names(a) <- names(r)
plot(a[[1]])

cntr <-
  data.frame(
    Country = world$geounit %>% sort,
    Current = NA,
    'Fut_2050' = NA,
    'Fut_2080' = NA,
    "Fut_2050_Gain" = NA,
    "Fut_2080_Gain" = NA,
    "Fut_2050_Loss" = NA,
    "Fut_2080_Loss" = NA
  ) #area in km2

for(i in 1:nrow(cntr)) {
  pol <- world[which(world$geounit == cntr$Country[i]), ]
  try(cntr$Current[i] <-
        raster::crop(a$CurrentDistributionBin, pol) %>% raster::mask(., pol) %>% cellStats(., sum))
  try(cntr$Fut_2050[i] <-
        raster::crop(a$FutureDistributionBin2050, pol) %>% raster::mask(., pol) %>% cellStats(., sum))
  try(cntr$Fut_2080[i] <-
        raster::crop(a$FutureDistributionBin2080, pol) %>% raster::mask(., pol) %>% cellStats(., sum))
}

# Gain and lost
r2 <-
  list.files('./Figures', pattern = '.tif', full.names = T) %>% stack
plot(r2)
r3 <-
  stack(list(
    r2$FutureDistributionBin2050 + (r2$CurrentDistributionBin * 2),
    r2$FutureDistributionBin2080 + (r2$CurrentDistributionBin *
                                      2)
  ))
r3[[1]][r3[[1]][] == 0] <- NA
r3[[2]][r3[[2]][] == 0] <- NA
crs(r3$layer.1) <- crs(r)
crs(r3$layer.2) <- crs(r)

a2050 <- area(r3$layer.1, na.rm = T)
a2080 <- area(r3$layer.2, na.rm = T)

plot(r3$layer.1, col = pals::parula(3))
mask(a2050, r3$layer.1 == 1, maskvalue = FALSE) %>% plot #gain
mask(a2050, r3$layer.1 == 2, maskvalue = FALSE) %>% plot #loss
mask(a2050, r3$layer.1 == 3, maskvalue = FALSE) %>% plot #stable


# Calc of area for each country for each futur range (i.e. expantion and loss)
i <- 33 #test with Argentina
pol <- world[which(world$geounit == cntr$Country[i]), ]
par(mfrow = c(2, 2))
asd <- mask(a2050, r3$layer.1 == 1, maskvalue = FALSE) %>%
  crop(., pol) %>% mask(., pol)
plot(asd, col = 'red', main = cellStats(asd, sum) %>% round)
plot(pol, add = T)
asd <- mask(a2050, r3$layer.1 == 2, maskvalue = FALSE) %>%
  crop(., pol) %>% mask(., pol)
plot(asd, col = 'red', main = cellStats(asd, sum) %>% round)
plot(pol, add = T)
asd <- mask(a2080, r3$layer.2 == 1, maskvalue = FALSE) %>%
  crop(., pol) %>% mask(., pol)
plot(asd, col = 'red', main = cellStats(asd, sum) %>% round)
plot(pol, add = T)
asd <- mask(a2080, r3$layer.2 == 2, maskvalue = FALSE) %>%
  crop(., pol) %>% mask(., pol)
plot(asd, col = 'red', main = cellStats(asd, sum) %>% round)
plot(pol, add = T)
par(mfrow = c(1, 1))

for(i in 1:nrow(cntr)) {
  pol <- world[which(world$geounit == cntr$Country[i]), ]
  try(cntr$Fut_2050_Gain[i] <-
        mask(a2050, r3$layer.1 == 1, maskvalue = FALSE) %>%
        crop(., pol) %>% mask(., pol) %>% cellStats(., sum))
  try(cntr$Fut_2050_Loss[i] <-
        mask(a2050, r3$layer.1 == 2, maskvalue = FALSE) %>%
        crop(., pol) %>% mask(., pol)  %>% cellStats(., sum))
  try(cntr$Fut_2080_Gain[i] <-
        mask(a2080, r3$layer.2 == 1, maskvalue = FALSE) %>%
        crop(., pol) %>% mask(., pol) %>% cellStats(., sum))
  try(cntr$Fut_2080_Loss[i] <-
        mask(a2080, r3$layer.2 == 2, maskvalue = FALSE) %>%
        crop(., pol) %>% mask(., pol) %>% cellStats(., sum))
}

dim(cntr2)
cntr2 <- cntr[which(rowSums(cntr[, -1] == 0) < 7), ]
cntr2 <- cntr2 %>% arrange((Current))
cntr_order <- cntr2$Country

# world
dim(cntr2)
world <- matrix(NA, nrow = 1, ncol = 8) %>% data.frame()
colnames(world) <- colnames(cntr2)
world[, 2:8] <- c(
  cellStats(a, sum),
  Fut_2050_Gain = cellStats(mask(a2050, r3$layer.1 == 1, maskvalue = FALSE), sum),
  Fut_2050_Loss = cellStats(mask(a2050, r3$layer.1 == 2, maskvalue = FALSE), sum),
  Fut_2080_Gain = cellStats(mask(a2080, r3$layer.1 == 1, maskvalue = FALSE), sum),
  Fut_2080_Loss = cellStats(mask(a2080, r3$layer.1 == 2, maskvalue = FALSE), sum)
)
world$Country <- "World"
cntr2 <- rbind(cntr2, world)
tail(cntr2)
readr::write_tsv(cntr2, './Figures/area_country.txt')


# Figures of invadable area by country
cntr2 <- read_tsv('./Figures/area_country.txt')
cntr3 <- reshape2::melt(cntr2)

cntr_invaded <-
  c(
    'Japan',
    'Argentina',
    'Brazil',
    'Uruguay',
    'Mexico',
    'Spain',
    'United States of America',
    'South Africa',
    'New Zealand',
    'Australia'
  )
cntr_risk <-
  c(
    'Italy',
    'Portugal',
    "United Kingdom",
    'France',
    'Ecuador',
    'Peru',
    'Bolivia',
    'Colombia',
    'Venezuela',
    'Turkey',
    'Sudan',
    'Republic of the Congo',
    'Kenya',
    'Madagascar'
  )
cntr_risk <- sort(cntr_risk)
cntr_risk[!cntr_risk %in% cntr3$Country]
cntr3$Country %>% sort %>% unique
invaded <- cntr3 %>% dplyr::filter(Country %in% cntr_invaded)
invrisk <- cntr3 %>% dplyr::filter(Country %in% cntr_risk)

cl <- c("gray20",  "#7bb3d1", "#016eae", "#cce88b", "#008837", "#f35926", "#b30000") 

A <- ggplot(invaded, aes(x = variable, y = (value), fill = variable)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_wrap(. ~ Country) +
  scale_fill_manual(values = cl) +
  labs(x = element_blank(), y = "Area (Km2)") +
  theme(
    panel.background  = element_rect(fill = 'gray90'),
    legend.position = 'bottom',
    axis.text.x = element_blank(),
    strip.background = element_blank()
  ) + 
  guides(fill = guide_legend(nrow = 1, title = element_blank()))
A

B <- ggplot(invrisk, aes(x = variable, y = (value), fill = variable)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_wrap(. ~ Country) +
  scale_fill_manual(values = cl) +
  labs(x = element_blank(), y = "Area (Km2)") +
  theme(
    panel.background  = element_rect(fill = 'gray90'),
    legend.position = 'bottom',
    axis.text.x = element_blank(),
    strip.background = element_blank()
  ) + 
  guides(fill = guide_legend(nrow = 1, title = element_blank()))

ggsave(
  plot = A,
  filename = './Figures/AreaByCountry_invaded.png',
  width = 14,
  dpi = 300,
  height = 10,
  units = 'cm',
  scale = 1.5
)
ggsave(
  plot = B,
  filename = './Figures/AreaByCountry_risk.png',
  width = 14,
  dpi = 300,
  height = 12,
  units = 'cm',
  scale = 1.5
)


# Area calculation for hotsopts
eco <-
  shapefile("./Ecorregions/Hotspots_2011/hotspots_2011_polygons.shp") # Shapefile with hotspots (source: https://www.arcgis.com/home/item.html?id=e5a7d024c4674cc185f78a99071feb07)
eco <- eco[which(eco$TYPE == 'hotspot_area'),]
eco %>% plot



ecodf <-
  data.frame(
    Hotspot = eco$NAME %>% sort,
    Current = NA,
    'Fut_2050' = NA,
    'Fut_2080' = NA,
    "Fut_2050_Gain" = NA,
    "Fut_2080_Gain" = NA,
    "Fut_2050_Loss" = NA,
    "Fut_2080_Loss" = NA
  ) #area in km2

for(i in 1:nrow(ecodf)) {
  pol <- eco[which(eco$NAME == ecodf$Hotspot[i]), ]
  try(ecodf$Current[i] <-
        crop(a$CurrentDistributionBin, pol) %>% mask(., pol) %>% cellStats(., sum))
  try(ecodf$Fut_2050[i] <-
        crop(a$FutureDistributionBin2050, pol) %>% mask(., pol) %>% cellStats(., sum))
  try(ecodf$Fut_2080[i] <-
        crop(a$FutureDistributionBin2080, pol) %>% mask(., pol) %>% cellStats(., sum))
  try(ecodf$Fut_2050_Gain[i] <-
        mask(a2050, r3$layer.1 == 1, maskvalue = FALSE) %>%
        crop(., pol) %>% mask(., pol) %>% cellStats(., sum))
  try(ecodf$Fut_2050_Loss[i] <-
        mask(a2050, r3$layer.1 == 2, maskvalue = FALSE) %>%
        crop(., pol) %>% mask(., pol)  %>% cellStats(., sum))
  try(ecodf$Fut_2080_Gain[i] <-
        mask(a2080, r3$layer.2 == 1, maskvalue = FALSE) %>%
        crop(., pol) %>% mask(., pol) %>% cellStats(., sum))
  try(ecodf$Fut_2080_Loss[i] <-
        mask(a2080, r3$layer.2 == 2, maskvalue = FALSE) %>%
        crop(., pol) %>% mask(., pol) %>% cellStats(., sum))
}

readr::write_tsv(ecodf, './Figures/area_hotspots.txt')


ecodf <- readr::read_tsv('./Figures/area_hotspots.txt')

{ecodf$Hotspot[7] <- "Chilean Winter Rainfall\nand Valdivian Forests"
  ecodf$Hotspot[8] <- "Coastal Forests\nof Eastern Africa"
  ecodf$Hotspot[9] <-"East Melanesian\nIslands"
  ecodf$Hotspot[18] <-"Madagascar and the\nIndian Ocean Islands" 
  ecodf$Hotspot[12] <-"Guinean Forests\nof West Africa"
  ecodf$Hotspot[19] <-"Madrean\nPine-Oak Woodlands"
  ecodf$Hotspot[20] <-"Maputaland\nPondoland-Albany"
  ecodf$Hotspot[24] <-"Mountains of\nSouthwest China"
  ecodf$Hotspot[23] <- "Mountains of\nCentral Asia"
  ecodf$Hotspot[33] <- "Tumbes\nChoco-Magdalena"  
  ecodf$Hotspot[35] <- "Western Ghats\nand Sri Lanka"   
}

ecodf2 <- reshape2::melt(ecodf)
A <- ggplot(ecodf2, aes(x = variable, y = (value), fill = variable)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  facet_wrap(. ~ Hotspot) +
  scale_fill_manual(values = cl) +
  labs(x = element_blank(), y = "Area (Km2)") +
  theme(
    panel.background  = element_rect(fill = 'gray90'),
    legend.position = 'bottom',
    axis.text.x = element_blank(),
    strip.background = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, title = element_blank()))
A

ggsave(
  plot = A,
  filename = './Figures/AreaByHotspot.png',
  width = 16,
  dpi = 300,
  height = 16,
  units = 'cm',
  scale=1.5
)


# Figure of PCA ------
db <-
  "./Model/PCA_Current_ClimaSoil_20km/Table/pca_varexplained_coef.txt" %>% readr::read_tsv()

db <-
  db %>% mutate(Variable = names(variables)) %>% dplyr::relocate(1:3, Variable)

db2 <- db %>%
  dplyr::select(Variable:PC3) %>%
  mutate(x0 = 0, y0 = 0)

db2 <-
  db2 %>% dplyr::mutate(Variable = gsub("_mean", "", gsub("CHELSA_", "BC", db2$Variable)))

db2$type <- "Climatic"

db2$type[20:nrow(db2)] <- "Edaphic"

readr::write_tsv(db2,
                 './Model/PCA_Current_ClimaSoil_20km/Table/pca_plot.txt')

db2 <-
  vroom::vroom('./Model/PCA_Current_ClimaSoil_20km/Table/pca_plot.txt')



# Distribution suitability
thr <-
  './Model/Result_30bin_ClimateEdaphic_20km/Thresholds_Ensemble.txt' %>% read_tsv()
thr <- thr %>% dplyr::filter(Algorithm == 'PCA')
thr

ClimaSoil <-
  './Model/Result_30bin_ClimateEdaphic_20km/Ensemble' %>% list.dirs()
ClimaSoil <-
  list.files(ClimaSoil, full.names = T) %>% grep("PCA/", ., value = T)
ens1 <- raster(ClimaSoil[1])
ens2 <- ens1 >= thr$THR_VALUE[5] #threshold that maximize Jaccard

ext_current <-
  raster(
    './Model/Result_30bin_ClimateEdaphic_20km/Extrapolation/Ligustrum_lucidum_MOP.tif'
  )

ens2[which((ext_current < 0.91)[])] <- 0 #threshold

pc <-
  stack((
    "./Model/PCA_Current_ClimaSoil_20km/" %>% list.files(., full.names = T)
  )[1:2])
ens1 <- stack(pc, ens1)

m <-
  raster("./Model/Result_30bin_ClimateEdaphic_20km/Extent_Masks/Ligustrum_lucidum.tif")

m <- crop(m, ens1)
ens1 <- mask(ens1, m)
df_suit <- rtodf(ens1) %>% as_tibble()
colnames(df_suit) <- c('x', 'y', 'PC1', 'PC2', 'suit')
nrow(df_suit)
df_suit <- df_suit[df_suit$PC1 < 6,]

A <- ggplot(db2, aes(PC1, PC2, label = Variable)) +
  geom_hline(yintercept = 0, col = 'gray50') +
  geom_vline(xintercept = 0, col = 'gray50') +
  geom_segment(aes(
    x = x0,
    y = y0,
    xend = PC1,
    yend = PC2,
    col = type
  ),
  arrow = arrow(type = "closed", angle = 15, unit(0.30, "cm"))) +
  scale_color_manual(values = pals::jet(10)[c(3, 8)]) +
  ggrepel::geom_text_repel() +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = 'bottom')


B <- ggplot(df_suit[seq(1, nrow(df_suit), by = 2),],
            aes(PC1, PC2, col = suit)) +
  geom_point(alpha = 0.9) +
  scale_color_gradientn(colours = pals::jet(10)) +
  geom_hline(yintercept = 0, col = 'gray50') +
  geom_vline(xintercept = 0, col = 'gray50') +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = 'bottom')
B

ggsave(
  filename = '"./Model/Figures/pca_plot.png',
  plot = A,
  units = 'cm',
  width = 15,
  height = 8
  # # dpi = 400,
  # scale = 1.2
)
