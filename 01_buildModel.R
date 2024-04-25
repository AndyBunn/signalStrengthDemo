################################################################################
#
# This reads in a text file of 100 samples with x,y coordinates containing the 
# signal strengths from three stations (e.g., station01) and three predictor 
# variables (e.g., predictor01). These are used in a regression kriging model
# to predict surfaces of signal strength using arrays (rasters stored as tifs,
# e.g., `predictor01.tif`) of the predictor variables over the study extent.
# Those surfaces are then saved.
#
# Inputs: signalStrengthSamples.csv
#         predictor01.tif, predictor02.tif, predictor03.tif
# 
# Outputs: station01_surface.tif, station02_surface.tif, station03_surface.tif
#
################################################################################
library(tidyverse)
library(terra)
library(tidyterra)
library(gstat)

# read in sample data (e.g., the data that comes from the GPS-enabled tags)
samples <- read_csv("data/signalStrengthSamples.csv")
head(samples)

# read in the predictor variables
predictor01 <- rast("data/predictor01.tif")
predictor02 <- rast("data/predictor02.tif")
predictor03 <- rast("data/predictor03.tif")
allPredictors <- c(predictor01,predictor02,predictor03)
names(allPredictors) <- c("predictor01","predictor02","predictor03")
allPredictors_df <- as.data.frame(allPredictors,xy=TRUE)

##########################
### Process station01

# Step 1. Create a gstat object with the regression model, include coordinate info
# this is a naive prediction without using the spatial structure of the
# model residuals
station01_gstat <- gstat(id="signalStrength",
  formula = station01~predictor01+predictor02+predictor03,
  locations = ~x+y,
  data=samples)
# initial skill
station01_gstat_skill <- gstat.cv(object=station01_gstat,nfold=10)
# plot
plot(station01_gstat_skill$observed, station01_gstat_skill$station01.pred)
plot(station01_gstat_skill$residual)
# r2
cor(station01_gstat_skill$observed, station01_gstat_skill$signalStrength.pred)^2 

# Step 2. Fit an emiprical variogram to the model residuals and plot it.
station01_gstat_obsVariogram <- variogram(station01_gstat)
plot(station01_gstat_obsVariogram)
# in this case a spherical model is probably best.
# NB. The variogram model will need to be fit to data and spherical might not be the most
# appropriate.

# Step 3. Fit the variogram model and plot it
initial_variogramModel <- vgm(model = "Sph")
station01_gstat_fittedVariogram <- fit.variogram(object = station01_gstat_obsVariogram, 
                                                 model = initial_variogramModel)
plot(station01_gstat_obsVariogram, model = station01_gstat_fittedVariogram) 

# Step 4. Update the gstat object with the fitted variogram model and assess skill
station01_vgm_gstat <- gstat(g = station01_gstat, 
                         id="signalStrength", 
                         model=station01_gstat_fittedVariogram)
# Assess skill using k-fold cross validation
station01_gstat_skill <- gstat.cv(object=station01_vgm_gstat,nfold=10)
#r2
cor(station01_gstat_skill$observed, station01_gstat_skill$signalStrength.pred)^2 

# Step 5. Predict the model using the predictor rasters
# station 01 predicted surface
station01_PredictedSurface <- predict(station01_vgm_gstat,newdata=allPredictors_df)
# convert to a raster for plotting and exporting
station01_PredictedSurface <- rast(station01_PredictedSurface)
writeRaster(x = station01_PredictedSurface, filename = "data/station01_surface.tif",overwrite=TRUE)

# example plot
ggplot() +
  geom_spatraster(data = station01_PredictedSurface$signalStrength.pred) +
  geom_spatraster_contour(data = station01_PredictedSurface$signalStrength.pred,
                          bins=10,color="white") +
  coord_equal() +
  scale_fill_continuous(type = "viridis") +
  labs(title = "Station01 Predicted Signal Strength") +
  theme_void()

ggplot() +
  geom_spatraster(data = station01_PredictedSurface$signalStrength.var) +
  geom_spatraster_contour(data = station01_PredictedSurface$signalStrength.var,
                          bins=10,color="white") +
  coord_equal() +
  scale_fill_continuous(type = "viridis") +
  labs(title = "Station01 Variance of Predicted Signal Strength") +
  theme_void()


##########################
### Process station02

# Step 1. Create a gstat object with the regression model, include coordinate info
# this is a naive prediction without using the spatial structure of the
# model residuals
station02_gstat <- gstat(id="signalStrength",
                         formula = station02~predictor01+predictor02+predictor03,
                         locations = ~x+y,
                         data=samples)
# initial skill
station02_gstat_skill <- gstat.cv(object=station02_gstat,nfold=10)
# plot
plot(station02_gstat_skill$observed, station02_gstat_skill$station02.pred)
plot(station02_gstat_skill$residual)
# r2
cor(station02_gstat_skill$observed, station02_gstat_skill$signalStrength.pred)^2 

# Step 2. Fit an emiprical variogram to the model residuals and plot it.
station02_gstat_obsVariogram <- variogram(station02_gstat)
plot(station02_gstat_obsVariogram)
# in this case a spherical model is probably best.
# NB. The variogram model will need to be fit to data and spherical might not be the most
# appropriate.

# Step 3. Fit the variogram model and plot it
initial_variogramModel <- vgm(model = "Sph")
station02_gstat_fittedVariogram <- fit.variogram(object = station02_gstat_obsVariogram, 
                                                 model = initial_variogramModel)
plot(station02_gstat_obsVariogram, model = station02_gstat_fittedVariogram) 

# Step 4. Update the gstat object with the fitted variogram model and assess skill
station02_vgm_gstat <- gstat(g = station02_gstat, 
                             id="signalStrength", 
                             model=station02_gstat_fittedVariogram)
# Assess skill using k-fold cross validation
station02_gstat_skill <- gstat.cv(object=station02_vgm_gstat,nfold=10)
#r2
cor(station02_gstat_skill$observed, station02_gstat_skill$signalStrength.pred)^2 

# Step 5. Predict the model using the predictor rasters
# station 01 predicted surface
station02_PredictedSurface <- predict(station02_vgm_gstat,newdata=allPredictors_df)
# convert to a raster for plotting and exporting
station02_PredictedSurface <- rast(station02_PredictedSurface)
writeRaster(x = station02_PredictedSurface, filename = "data/station02_surface.tif",overwrite=TRUE)

# example plot
ggplot() +
  geom_spatraster(data = station02_PredictedSurface$signalStrength.pred) +
  geom_spatraster_contour(data = station02_PredictedSurface$signalStrength.pred,
                          bins=10,color="white") +
  coord_equal() +
  scale_fill_continuous(type = "viridis") +
  labs(title = "Station02 Predicted Signal Strength") +
  theme_void()

ggplot() +
  geom_spatraster(data = station02_PredictedSurface$signalStrength.var) +
  geom_spatraster_contour(data = station02_PredictedSurface$signalStrength.var,
                          bins=10,color="white") +
  coord_equal() +
  scale_fill_continuous(type = "viridis") +
  labs(title = "Station02 Variance of Predicted Signal Strength") +
  theme_void()


##########################
### Process station03

# Step 1. Create a gstat object with the regression model, include coordinate info
# this is a naive prediction without using the spatial structure of the
# model residuals
station03_gstat <- gstat(id="signalStrength",
                         formula = station03~predictor01+predictor02+predictor03,
                         locations = ~x+y,
                         data=samples)
# initial skill
station03_gstat_skill <- gstat.cv(object=station03_gstat,nfold=10)
# plot
plot(station03_gstat_skill$observed, station03_gstat_skill$station03.pred)
plot(station03_gstat_skill$residual)
# r2
cor(station03_gstat_skill$observed, station03_gstat_skill$signalStrength.pred)^2 

# Step 2. Fit an emiprical variogram to the model residuals and plot it.
station03_gstat_obsVariogram <- variogram(station03_gstat)
plot(station03_gstat_obsVariogram)
# in this case a spherical model is probably best.
# NB. The variogram model will need to be fit to data and spherical might not be the most
# appropriate.

# Step 3. Fit the variogram model and plot it
initial_variogramModel <- vgm(model = "Sph")
station03_gstat_fittedVariogram <- fit.variogram(object = station03_gstat_obsVariogram, 
                                                 model = initial_variogramModel)
plot(station03_gstat_obsVariogram, model = station03_gstat_fittedVariogram) 

# Step 4. Update the gstat object with the fitted variogram model and assess skill
station03_vgm_gstat <- gstat(g = station03_gstat, 
                             id="signalStrength", 
                             model=station03_gstat_fittedVariogram)
# Assess skill using k-fold cross validation
station03_gstat_skill <- gstat.cv(object=station03_vgm_gstat,nfold=10)
#r2
cor(station03_gstat_skill$observed, station03_gstat_skill$signalStrength.pred)^2 

# Step 5. Predict the model using the predictor rasters
# station 01 predicted surface
station03_PredictedSurface <- predict(station03_vgm_gstat,newdata=allPredictors_df)
# convert to a raster for plotting and exporting
station03_PredictedSurface <- rast(station03_PredictedSurface)
writeRaster(x = station03_PredictedSurface, filename = "data/station03_surface.tif",overwrite=TRUE)

# example plot
ggplot() +
  geom_spatraster(data = station03_PredictedSurface$signalStrength.pred) +
  geom_spatraster_contour(data = station03_PredictedSurface$signalStrength.pred,
                          bins=10,color="white") +
  coord_equal() +
  scale_fill_continuous(type = "viridis") +
  labs(title = "Station03 Predicted Signal Strength") +
  theme_void()

ggplot() +
  geom_spatraster(data = station03_PredictedSurface$signalStrength.var) +
  geom_spatraster_contour(data = station03_PredictedSurface$signalStrength.var,
                          bins=10,color="white") +
  coord_equal() +
  scale_fill_continuous(type = "viridis") +
  labs(title = "Station03 Variance of Predicted Signal Strength") +
  theme_void()

