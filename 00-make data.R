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

# pred01
pred01 <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure01 <- cbind(x = runif(50, 2, 99),
                   y=runif(50, 2, 99),
                   z = rnorm(50))
pred01 <- interpIDW(pred01, structure01, radius=30, power=1, smooth=1)

# pred02
pred02 <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure02 <- cbind(x = runif(100, 2, 99),
                   y=runif(100, 2, 99),
                   z = rnorm(100))
pred02 <- interpIDW(pred02, structure02, radius=20, power=1, smooth=1)

# pred03
pred03 <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure03 <- cbind(x = runif(1000, 2, 99),
                     y=runif(1000, 2, 99),
                     z = rnorm(1000))
pred03 <- interpIDW(pred03, structure03, radius=10, power=1, smooth=1)

# noise
epsilon <- rast(ncol=n, nrow=n, crs="", 
                xmin=1, xmax=100, ymin=1, ymax=100,
                vals=rnorm(n^2))
plot(pred01)
plot(pred02)
plot(pred03)
plot(epsilon)


# now make station surfaces

# station01
r  <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure01 <- cbind(x = runif(100, 75, 85),
                     y=runif(100, 75, 85),
                     z = rnorm(100,mean=5,sd=5))
structure01 <- rbind(structure01,
                     cbind(x = runif(10, 50, 100),
                     y=runif(10, 50, 100),
                     z = rnorm(10,mean=3,sd=1)))

station01 <- interpIDW(r, structure01, radius=80, power=10,fill=0,smooth=10)

station01 <- station01 + 0.1 * pred01 + 0.2 * pred02 + 0.1 * pred03 + 0.01 * epsilon


# station02
r  <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure02 <- cbind(x = runif(100, 25, 35),
                     y=runif(100, 55, 65),
                     z = rnorm(100,mean=5,sd=5))
structure02 <- rbind(structure02,
                     cbind(x = runif(10, 50, 100),
                           y=runif(10, 50, 100),
                           z = rnorm(10,mean=3,sd=1)))

station02 <- interpIDW(r, structure02, radius=80, power=10,fill=0,smooth=10)

station02 <- station02 + -0.1 * pred01 + 0.25 * pred02 + -0.1 * pred03 + 0.01 * epsilon

# station03
r  <- rast(ncol=n, nrow=n, crs="", xmin=1, xmax=100, ymin=1, ymax=100)
structure03 <- cbind(x = runif(100, 48, 52),
                     y=runif(100, 55, 65),
                     z = rnorm(100,mean=5,sd=5))
structure03 <- rbind(structure03,
                     cbind(x = runif(10, 0, 100),
                           y=runif(10, 0, 100),
                           z = rnorm(10,mean=3,sd=1)))

station03 <- interpIDW(r, structure03, radius=80, power=10,fill=0,smooth=10)

station03 <- station03 + -0.01 * pred01 + -0.1 * pred02 + 0.3 * pred03 + 0.01 * epsilon


ggplot() +
  geom_spatraster(data = station01) +
  geom_spatraster_contour(data = station01,bins=10,color="white") +
  coord_equal() +
  scale_fill_continuous(type = "viridis") +
  theme_void()



## Now create samples

stationStack <- c(station01,station02,station03)
samps <- spatSample(stationStack, size=500, xy=TRUE)
names(samps) <- c("x","y","station01","station02","station03")
write_csv(samps,file = "data/signalStrengthSamples.csv")

## And save pred layers
writeRaster(pred01,filename = "data/predictor01.tif")
writeRaster(pred02,filename = "data/predictor02.tif")
writeRaster(pred03,filename = "data/predictor03.tif")


