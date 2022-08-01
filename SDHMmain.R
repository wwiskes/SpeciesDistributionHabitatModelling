# ---
# This script requires a functions script titled "SDHMfunctions"
# This script is designed to create ensemble species distribution models from presence point data and predictor layers
# Written by William Wiskes
# Last update 6/22/2021
# ---

#at 32gb ram & 16cpu this takes around 1hr to run
library(tidyverse)
library(raster)
library(sp)
library(sf)
library(dplyr) # data manipulation
#options(scipen = 999)#disabling scientific notation
#load in functions
source("SDHMfunctions.R")
#The data must be in EPSG 5070. The proj4string below is for 5070, if you are running this script in the cloud you must use proj4
proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "

#Make your list of rasters from the document found here (lists can be any length greater than 1):
#Shown below are the 'standard' stacks for all stacked variables in the format extent/resolution
# https://storage.googleapis.com/predictors_structure/structure.csv
#Utah/1km
rasterList <- c("terrestrial/gradientmetrics/topo/1km/gm_allvars_topo_ut.tif",
                "terrestrial/landfire/dist/1km/lf_allvars_dist_ut.tif",
                "terrestrial/landfire/topo/1km/lf_allvars_topo_ut.tif",
                "terrestrial/landfire/veg/1km/lf_allvars_veg_ut.tif",
                "terrestrial/nlcd/1km/nlcd_allvars_veg_ut.tif",
                "terrestrial/polaris/soil/1km/polaris_allvars_ut.tif",
                "terrestrial/usfws/km_1/prism/ppt_30yr_normal_stack_ut.tif",
                "terrestrial/usfws/km_1/prism/tmax_30yr_normal_stack_ut.tif",
                "terrestrial/usfws/km_1/prism/tmin_30yr_normal_stack_ut.tif",
                "terrestrial/usfws/km_1/usfws_allvars_clim_ut.tif", 
                "terrestrial/usfws/km_1/usfws_allvars_geo_ut.tif",
                "terrestrial/usfws/km_1/usfws_allvars_hydro_ut.tif",
                "terrestrial/usfws/km_1/usfws_allvars_soil_ut.tif",
                "terrestrial/usfws/km_1/usfws_allvars_topo_ut.tif",
                "terrestrial/usfws/km_1/usfws_allvars_veg_ut.tif")
# #Utah/100m
# rasterList <- c("terrestrial/gradientmetrics/topo/100m/gm_allvars_topo_ut.tif",
#                 "terrestrial/landfire/topo/100m/lf_allvars_topo_ut.tif",
#                 "terrestrial/landfire/veg/100m/lf_allvars_veg_ut.tif",
#                 "terrestrial/nlcd/100m/nlcd_allvars_ut.tif",
#                 "terrestrial/polaris/soil/100m/polaris_allvars_ut.tif")
# #Utah/100m fws only
# rasterList <- c("terrestrial/usfws/m_100/ppt_stack_100m_ut.tif",
#                 'terrestrial/usfws/m_100/tmax_stack_100m_ut.tif',
#                 "terrestrial/usfws/m_100/tmin_stack_100m_ut.tif",
#                 "terrestrial/usfws/m_100/usfws_climate_100m_stack_ut.tif",
#                 "terrestrial/usfws/m_100/usfws_geo_100m_stack_ut.tif",
#                 "terrestrial/usfws/m_100/usfws_hydro_100m_stack_ut.tif",
#                 "terrestrial/usfws/m_100/usfws_topo_100m_stack_ut.tif")
# # #Utah/30m
# rasterList <- c("terrestrial/gradientmetrics/topo/30m/gm_allvars_topo_ut.tif",
#                 "terrestrial/landfire/topo/30m/lf_allvars_topo_ut.tif",
#                 "terrestrial/landfire/veg/30m/lf_allvars_veg_ut.tif",
#                 "terrestrial/nlcd/30m/nlcd_allvars_ut.tif")
# #WNA/1km
# rasterList <- c("terrestrial/gradientmetrics/topo/1km/gm_allvars_topo_wna.tif",
#                 "terrestrial/landfire/dist/1km/lf_allvars_dist_wna.tif",
#                 "terrestrial/landfire/topo/1km/lf_allvars_topo_wna.tif",
#                 "terrestrial/landfire/veg/1km/lf_allvars_veg_wna.tif",
#                 "terrestrial/usfws/km_1/usfws_allvars_clim_wna.tif",
#                 "terrestrial/usfws/km_1/usfws_allvars_geo_wna.tif",
#                 "terrestrial/usfws/km_1/usfws_allvars_hydro_wna.tif",
#                 "terrestrial/usfws/km_1/usfws_allvars_soil_wna.tif",
#                 "terrestrial/usfws/km_1/usfws_allvars_topo_wna.tif",
#                 "terrestrial/usfws/km_1/usfws_allvars_veg_wna.tif",
#                 "terrestrial/usfws/km_1/prism/tmax_30yr_normal_stack_wna.tif",
#                 "terrestrial/usfws/km_1/prism/ppt_30yr_normal_stack_wna.tif",
#                 "terrestrial/usfws/km_1/prism/t,tmin_30yr_normal_stack_wna.tif",
#                 "terrestrial/polaris/soil/1km/polaris_allvars_wna.tif")
# #WNA/100m
# rasterList <- c("terrestrial/usfws/m_100/usfws_climate_100m_stack_wna.tif",
#                 "terrestrial/usfws/m_100/tmin_stack_100m_wna.tif",
#                 "terrestrial/usfws/m_100/usfws_geo_100m_stack_wna.tif",
#                 "terrestrial/usfws/m_100/tmax_stack_100m_wna.tif",
#                 "terrestrial/usfws/m_100/usfws_hydro_100m_stack_wna.tif",
#                 "terrestrial/usfws/m_100/ppt_stack_100m_wna.tif",
#                 "terrestrial/usfws/m_100/usfws_topo_100m_stack_wna.tif")
# #WNA/30m - we do not yet have a template for this resolution at WNA scale
# rasterList <- c("terrestrial/gradientmetrics/topo/gm_allvars_topo_wna.tif",
#                 "terrestrial/landfire/dist/lf_allvars_dist_wna.tif",
#                 "terrestrial/landfire/topo/lf_allvars_topo_wna.tif",
#                 "terrestrial/landfire/veg/lf_allvars_veg_wna.tif")


#Load in Utah extent simple feature
#Replace this with your own study area if needed. MUST BE IN EPSG 5070
blob <- st_read("https://storage.googleapis.com/predictors_public/bounds/huc_aea.geojson") 
st_crs(blob) <- proj 

#1000 random points within the state of utah have been provided as a placeholder
pointData <- read.csv("randompoints.csv")
#They are in the projection 5070, and so must be anything you replace them with.
#The dataset must have 3 columns. 
#The first column is required but not used in these functions. But could be a 'keep' column to filter off of
#The second column must be X (longitude)
#The third column must be Y (latitude)
########################
# it is also an option to pull data from a postgres database
# pointData <- queryPostgres("dkm") 
# queryPostgres code will not run without a env script.
# it is also an option to pull data from a bigquery database
# pointData <- queryBiobase("randompoints") 
# queryBiobase is intended to be used in a cloud environment

# +
# for testing, only use a small subset of data
# pointData <- head(pointData,100) 
# pointData
# -

#They are in the projection 5070, please set x&y coords accordingly
coords <- colnames(head(pointData)[,2:3])
pointSF <- st_as_sf(pointData, coords = coords, crs = proj)
head(pointSF)

#Set Extent
#The extent is set from the study area named "Blob", if you do not have a study area you 
#can generate one from your point data below
ext <- extent(pointSF)
blob <- st_as_sf(as(ext, "SpatialPolygons"))
st_crs(blob) <- proj 
#set buffer, buffer distance is in meters WARNING, make certain buffer does not extend outside raster extent
blob <- st_buffer(blob, 5000)
#from study area:
ext <- extent(blob)

#Load in modelling extent. Currently do not have 30m for WNA
#1km works for WNA and Utah
temp <- raster("/vsicurl/https://storage.googleapis.com/predictors_public/bounds/1km_template.tif")
#100m works for WNA and Utah
#temp <- raster("/vsicurl/https://storage.googleapis.com/predictors_public/bounds/100m_template_wna.tif")
#30m works ONLY for Utah
#temp <- raster("/vsicurl/https://storage.googleapis.com/predictors_public/bounds/30m_template.tif")
#crop modelling extent template
temp <- crop(temp, ext)

#Generate pseudo absences from the extent of the study area
#Needs the points simple feature, the blob geojson, and the template raster
pointPseudo <- pseudoFunction(pointSF, blob, temp)

head(pointPseudo) #x&y columns always lat/long

#Run your point data against the rasters to extract the prediction values
data <- extractStack(pointPseudo, rasterList)
column_names <- colnames(head(data)[,9:ncol(data)])
## Remove columns and rows with more than 50% NA
data<- data[which(rowMeans(!is.na(data[,9:ncol(data)])) > 0.5), which(colMeans(!is.na(data)[,9:ncol(data)]) > 0.5)]
head(data)

#set the threshold at which a column is no longer statistically relevant
cut <- 0.8
# declare list of columns to preserve or remove regardless of spearmans. To leave blank use: <- c("")
preserve <- c("")
remove <- c("")
#remove warnings
options (warn = - 1)
#remove columns without statistical correlation, with exception to preserve/remove
cutF <- cutFunction(data, cut, preserve, remove)
#for checking the output of the spearmans the table is provided here, this is not used in subsequent steps
spearmans <- cutF[[2]]
cutData <- cutF[[1]]
head(cutData)
boxFunction(cutData)

#this function not only makes the raster stack, but renames each layer to the correct corresponding column
rasters <- rasterStack(cutData, rasterList, column_names)
#optional, crop rasters to the modelling extent.
rasters <- crop(rasters, ext)

##########################################START TESTING
source("SDHMfunctions.R")
############################################END TESTING

#Generalized linear model
glm <- glmFunction(cutData, rasters)

plot(glm)
plot(blob$geometry, add =T)

#Generalized additive model
gam <- gamFunction(cutData, rasters)

plot(gam)
plot(blob$geometry, add =T)

#MaxEnt model
max <- maxFunction(cutData, rasters)

plot(max)
plot(blob$geometry, add =T)

#Random Forest model
raf <- rafFunction(cutData, rasters)

plot(raf)
plot(blob$geometry, add =T)

#Boosted Regression Tree model
brt <- brtFunction(cutData, rasters)

plot(brt)
plot(blob$geometry, add =T)

#Ensemble model
ens <- stack(glm,gam,raf,max,brt)
ensemble <- sum(ens)
plot(ensemble, axes = T, main = "concordance Class map: Ramp")
plot(blob$geometry, add =T)


writeRaster(ensemble, filename = "ensemble.img", overwrite = T)

