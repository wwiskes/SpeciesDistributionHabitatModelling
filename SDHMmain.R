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
options(scipen = 999)#disabling scientific notation - important for table readability 
#load in functions
source("SDHMfunctions.R")
#The data must be in EPSG 5070. The proj4string below is for 5070, if you are running this script in the cloud you must use proj4
proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "
#Load in Utah 1km bounds - modelling extent
temp <- raster("/vsicurl/https://storage.googleapis.com/predictors_public/bounds/1km_template.tif")
#Load in Utah Fnet
fnetSF <-st_read("/vsicurl/https://storage.googleapis.com/predictors_public/bounds/km1_fnet.geojson", crs = proj)

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
# pointData <- queryPostgres("ybcu") 
# queryPostgres code will not run without a env script.

# +
# for testing, only use a small subset of data
pointData <- head(pointData,100) 
pointData
# -

#They are in the projection 5070, please set x&y coords accordingly
coords <- colnames(head(pointData)[,2:3])
pointSF <- st_as_sf(pointData, coords = coords, crs = proj)
head(pointSF)

# +
#Set Extent
#The extent is set from the study area named "Blob", if you do not have a study area you 
#can generate one from your point data below
ext <- extent(pointSF)
blob <- st_as_sf(as(ext, "SpatialPolygons"))
st_crs(blob) <- proj 
#set buffer, buffer distance is in meters
blob <- st_buffer(blob, 5000)
#from study area:
ext <- extent(blob)
#crop modelling extent template
temp <- crop(temp, ext)
# -

#Generate pseudo absences from the extent of the study area
#Needs the points simple feature, the blob geojson, and the template raster
pointPseudo <- pseudoFunction(pointSF, blob, temp)

head(pointPseudo) #x&y columns always lat/long

#Make your list of rasters from the document found here (lists can be any length greater than 1):
#Shown below is the 'standard' stack for all 1km resolution variables in the state of Utah
# https://storage.googleapis.com/predictors_structure/structure.csv
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

#Run your point data against the rasters to extract the prediction values
data <- extractStack(pointPseudo, rasterList)
column_names <- colnames(head(data)[,9:ncol(data)])
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




