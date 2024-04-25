# this makes three surfaces of predictors and three surfaces of signal strength stations.
# then we sample the station data to be used as example data.
# we save both the predictors and the samples.
# this is the base data for the demo
rm(list=ls())
library(tidyverse)
library(terra)
library(tidyterra)

set.seed(3625)
# set grid dimensions
n <- 100
nCells <- n*n

# predictor01
predictor01 <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure01 <- cbind(x = runif(50, 2, 99),
                   y=runif(50, 2, 99),
                   z = rnorm(50))
predictor01 <- interpIDW(predictor01, structure01, radius=30, power=1, smooth=1)

# predictor02
predictor02 <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure02 <- cbind(x = runif(100, 2, 99),
                   y=runif(100, 2, 99),
                   z = rnorm(100))
predictor02 <- interpIDW(predictor02, structure02, radius=20, power=1, smooth=1)

# predictor03
predictor03 <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure03 <- cbind(x = runif(1000, 2, 99),
                     y=runif(1000, 2, 99),
                     z = rnorm(1000))
predictor03 <- interpIDW(predictor03, structure03, radius=10, power=1, smooth=1)

# noise
epsilon <- rast(ncol=n, nrow=n, crs="", 
                xmin=1, xmax=100, ymin=1, ymax=100,
                vals=rnorm(n^2))
plot(predictor01)
plot(predictor02)
plot(predictor03)
plot(epsilon)


# now make station surfaces
# empty raster
r  <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)

# station01
structure01 <- cbind(x = runif(100, 75, 85),
                     y=runif(100, 75, 85),
                     z = rnorm(100,mean=5,sd=5))
stationNoise <- cbind(x = runif(100, 0, 100),
                      y=runif(100, 0, 100),
                      z = rnorm(100,mean=3,sd=1))

structure01 <- rbind(structure01,
                     stationNoise)

station01 <- interpIDW(r, structure01, radius=80, power=5,fill=0,smooth=10)

station01 <- station01 + 0.25 * predictor01 + 0.2 * predictor02 + 0.13 * predictor03 + 0.05 * epsilon

# station02
r  <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure02 <- cbind(x = runif(100, 25, 35),
                     y=runif(100, 55, 65),
                     z = rnorm(100,mean=5,sd=5))
stationNoise <- cbind(x = runif(100, 0, 100),
                      y=runif(100, 0, 100),
                      z = rnorm(100,mean=3,sd=1))

structure02 <- rbind(structure02,
                     stationNoise)

station02 <- interpIDW(r, structure02, radius=80, power=5,fill=0,smooth=10)

station02 <- station02 + -0.1 * predictor01 + 0.25 * predictor02 + -0.1 * predictor03 + 0.06 * epsilon

# station03
r  <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure03 <- cbind(x = runif(100, 48, 52),
                     y=runif(100, 55, 65),
                     z = rnorm(100,mean=5,sd=5))
stationNoise <- cbind(x = runif(100, 0, 100),
                      y=runif(100, 0, 100),
                      z = rnorm(100,mean=3,sd=1))

structure03 <- rbind(structure03,
                     stationNoise)

station03 <- interpIDW(r, structure03, radius=80, power=5,fill=0,smooth=10)

station03 <- station03 + -0.01 * predictor01 + -0.1 * predictor02 + 0.3 * predictor03 + 0.04 * epsilon


ggplot() +
  geom_spatraster(data = station01) +
  geom_spatraster_contour(data = station01,bins=10,color="white") +
  coord_equal() +
  scale_fill_continuous(type = "viridis") +
  theme_void()

## Now create samples

dataStack <- c(station01,station02,station03,predictor01,predictor02,predictor03)
names(dataStack) <- c("station01","station02","station03","predictor01","predictor02","predictor03")
samples <- spatSample(dataStack, size=100, xy=TRUE)
head(samples)

# look at the variogram
samples_sf <- st_as_sf(samples,coords = c("x","y"))

station01_gstat <- gstat(id="station01",
  formula = station01~predictor01+predictor02+predictor03,
  locations = ~x+y,
  data=samples)

# Step 2. Fit an emiprical variogram to the model residuals and plot it.
station01_gstat_obsVariogram <- variogram(station01_gstat)
plot(station01_gstat_obsVariogram)


write_csv(samples,file = "data/signalStrengthSamples.csv")
 
## And save pred layers
writeRaster(predictor01,filename = "data/predictor01.tif",overwrite=TRUE)
writeRaster(predictor02,filename = "data/predictor02.tif",overwrite=TRUE)
writeRaster(predictor03,filename = "data/predictor03.tif",overwrite=TRUE)
