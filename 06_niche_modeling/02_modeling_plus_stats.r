
library(rJava)
library(rgeos)
library(sp)
library(maptools)
library(raster)
library(dismo)
library(rgdal)

countries <- readOGR("./ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp")
Ethiopia <- countries[countries$NAME == "Ethiopia", ]

# raster directories at 2.5 arc minutes
worldclim_directory <- "04_worldclim_layers"
ccsm_directory <- "06_ccsm_layers"
miroc_directory <- "05_miroc_layers"

# read in layers
worldclim_layers <- stack(list.files(worldclim_directory, full.names=T))
names(worldclim_layers) <- gsub("wc2.0_", "", names(worldclim_layers))
names(worldclim_layers) <- gsub("_2.5m_", "", names(worldclim_layers))
ccsm_layers <- stack(list.files(ccsm_directory, full.names=T))
names(ccsm_layers) <- gsub("cclgmbi", "bio", names(ccsm_layers))
miroc_layers <- stack(list.files(miroc_directory, full.names=T))
names(miroc_layers) <- gsub("mrlgmbi", "bio", names(miroc_layers))

# clip layers to eastern Africa
# max lat = 20, min lat = -20, max long = 55, min long = 10
clip_east_africa <- cbind(c(10, 55, 55, 10), c(20, 20, -20, -20))
clip_east_africa <- SpatialPolygons(list(Polygons(list(Polygon(clip_east_africa)), 1)))
worldclim_layers <- crop(worldclim_layers, clip_east_africa)
ccsm_layers <- crop(ccsm_layers, clip_east_africa)
miroc_layers <- crop(miroc_layers, clip_east_africa)

# find correlation of layers in worldclim dataset
layerStats(worldclim_layers, 'pearson', na.rm=T)$'pearson correlation coefficient'
# manually inspect output and choose the layers to keep
# here, we will keep all layers with less than 0.65 correlation
# bioclim layers = 1, 2, 3, 12, 14, 15
bioclim_keep <- c("bio1", "bio2", "bio3", "bio12", "bio14", "bio15")
worldclim_layers <- dropLayer(worldclim_layers, seq(1:length(names(worldclim_layers)))[is.na(match(names(worldclim_layers),bioclim_keep))])
ccsm_layers <- dropLayer(ccsm_layers, seq(1:length(names(ccsm_layers)))[is.na(match(names(ccsm_layers),bioclim_keep))])
miroc_layers <- dropLayer(miroc_layers, seq(1:length(names(miroc_layers)))[is.na(match(names(miroc_layers),bioclim_keep))])
# need to multiply current worldclim bio 1 and bio 2 * 10
worldclim_layers[[1]] <- worldclim_layers[[1]] * 10
worldclim_layers[[5]] <- worldclim_layers[[5]] * 10

# set up all lists for the entire process
input_files <- list()
training_points <- list()
testing_points <- list()
species_names <- list()
clipping_masks <- list()
raster_stacks <- list()
bioclim_values <- list()
niche_models_R <- list()
thresholds <- list()
projections_contemporary <- list()
projections_miroc <- list()
projections_ccsm <- list()
reclass_contemporary <- list()
reclass_miroc <- list()
reclass_ccsm <- list()
tests <- list()

# files of occurrence data
input_folder <- "03_final_points"
input_file_list <- list.files(input_folder, pattern="*rarefy.csv", full.name=T)
for(a in 1:length(input_file_list)) {
	input_files[[a]] <- read.csv(input_file_list[a])
	species_names[[a]] <- sapply(strsplit(sapply(strsplit(input_file_list[a], "/"), "[[", 2), "_"), "[[", 1)
}

# set up training and testing points
# set up percent of training and testing points
test_proportion <- 0.2
for(a in 1:length(input_files)) {
	a_rep <- input_files[[a]]
	a_rep <- data.frame(Long=a_rep[,2], Lat=a_rep[,1])
	a_test_points <- sample(1:nrow(a_rep), floor(nrow(a_rep) * test_proportion))
	training_points[[a]] <- a_rep[-a_test_points, ]
	testing_points[[a]] <- a_rep[a_test_points, ]
	
	# convert points to spatial points and add CRS WGS84
	coordinates(training_points[[a]]) <- c("Long", "Lat")
	proj4string(training_points[[a]]) <- CRS("+init=epsg:4326")
	coordinates(testing_points[[a]]) <- c("Long", "Lat")
	proj4string(testing_points[[a]]) <- CRS("+init=epsg:4326")
}

# set up training region for each species 
# here we will use 250 km (250 *1000 m)
model_buffer <- 250 * 1000
for(a in 1:length(input_files)) {
	a_rep <- training_points[[a]]
	
	# create a projection of the points in CRS 
	a_rep <- spTransform(a_rep, CRS=CRS("+init=epsg:3857"))
	
	# buffer the points by the model buffer and transform back to WGS84
	clipping_masks[[a]] <- spTransform(gBuffer(a_rep, width=model_buffer), CRS=CRS("+init=epsg:4326"))

	# create a training raster stack for each species
	raster_stacks[[a]] <- crop(worldclim_layers, clipping_masks[[a]])
	raster_stacks[[a]] <- mask(raster_stacks[[a]], clipping_masks[[a]])
}

# extract bioclim values for observed points and 500 background points 
for(a in 1:length(input_files)) {
	bioclim_values[[a]] <- extract(raster_stacks[[a]], training_points[[a]])
	a_rep <- sampleRandom(raster_stacks[[a]], 500)
	a_rep <- data.frame(rbind(cbind(presence=1, bioclim_values[[a]]), cbind(presence=0, a_rep)))
	bioclim_values[[a]] <- a_rep
}

# perform maxent models for each species
for(a in 1:length(input_files)) {
	niche_models_R[[a]] <- maxent(x=bioclim_values[[a]][ ,2:ncol(bioclim_values[[a]])], p=bioclim_values[[a]][ ,1], path="./07_maxent_output")
	projections_contemporary[[a]] <- predict(niche_models_R[[a]], worldclim_layers)
	projections_ccsm[[a]] <- predict(niche_models_R[[a]], ccsm_layers)
	projections_miroc[[a]] <- predict(niche_models_R[[a]], miroc_layers)
}

# model reclassify function
# model is already read in as a raster
# threshold is numeric between 0 and 1
reclass_ENM <- function(model, threshold) {
  require(raster)
  y <- model
  y[y >= threshold] <- 1
  y[y < threshold] <- NA
  return(y)
}

# test the models accuracy with testing points and reclassify simultaneously
for(a in 1:length(input_files)) {
  train_values <- extract(projections_contemporary[[a]], training_points[[a]])
  thresholds[[a]] <- sort(train_values)[floor(length(train_values) - length(train_values) * 0.9)]
  reclass_contemporary[[a]] <- reclass_ENM(projections_contemporary[[a]], thresholds[[a]])
  reclass_ccsm[[a]] <- reclass_ENM(projections_ccsm[[a]], thresholds[[a]])
  reclass_miroc[[a]] <- reclass_ENM(projections_miroc[[a]], thresholds[[a]])
  
  # pull out whether test points are predicted or not
  # and also null hypothesis of "How many points are predicted out of total points?"
  perc_predicted <- as.vector(table(reclass_contemporary[[a]]@data@values)[1] / (reclass_contemporary[[a]]@ncols * reclass_contemporary[[a]]@nrows))
  test_values <- extract(reclass_contemporary[[a]], testing_points[[a]])
  test_values[is.na(test_values)] <- 0
  tests[[a]] <- binom.test(length(test_values[test_values == 1]), length(test_values), p = perc_predicted, alternative="g")
}

# crop the thresholded models to just ethiopia
worldclim_ethiopia <- crop(worldclim_layers[[1]], Ethiopia)
worldclim_ethiopia <- mask(worldclim_ethiopia, Ethiopia)
for(a in 1:length(input_files)) {
  reclass_contemporary[[a]] <- crop(reclass_contemporary[[a]], Ethiopia)
  reclass_contemporary[[a]][is.na(reclass_contemporary[[a]])] <- 0
  reclass_contemporary[[a]] <- mask(reclass_contemporary[[a]], Ethiopia)
  reclass_ccsm[[a]] <- crop(reclass_ccsm[[a]], Ethiopia)
  reclass_ccsm[[a]][is.na(reclass_ccsm[[a]])] <- 0
  reclass_ccsm[[a]] <- mask(reclass_ccsm[[a]], Ethiopia)
  reclass_miroc[[a]] <- crop(reclass_miroc[[a]], Ethiopia)
  reclass_miroc[[a]][is.na(reclass_miroc[[a]])] <- 0
  reclass_miroc[[a]] <- mask(reclass_miroc[[a]], Ethiopia)
}



par(mfrow=c(2,3))
for(a in 1:6) {
  plot(reclass_contemporary[[a]])
  lines(Ethiopia)
}
for(a in 1:6) {
  plot(reclass_ccsm[[a]])
  lines(Ethiopia)
}
for(a in 1:6) {
  plot(reclass_miroc[[a]])
  lines(Ethiopia)
}

# make a polygon that connects the populations east and west of the GRV 
# with a 100 km buffer
grv_points <- rbind(c(7.095, 39.789), c(8.967, 38.572))
grv_points <- data.frame(Long=grv_points[,2], Lat=grv_points[,1])
coordinates(grv_points) <- c("Long", "Lat")
proj4string(grv_points) <- CRS("+init=epsg:4326")
grv_lines <- spLines(grv_points)
proj4string(grv_lines) <- CRS("+init=epsg:4326")
grv_lines2 <- spTransform(grv_lines, CRS=CRS("+init=epsg:3857"))
grv_points2 <- spTransform(grv_points, CRS=CRS("+init=epsg:3857"))
grv_buffer_size <- 100 * 1000
grv_buffer <- spTransform(gBuffer(grv_lines2, width=grv_buffer_size), CRS=CRS("+init=epsg:4326"))
east_buffer <- spTransform(gBuffer(grv_points2[1], width=grv_buffer_size), CRS=CRS("+init=epsg:4326"))
west_buffer <- spTransform(gBuffer(grv_points2[2], width=grv_buffer_size), CRS=CRS("+init=epsg:4326"))


plot(reclass_contemporary[[a]])
points(grv_points, pch=19, cex=0.8)
lines(grv_buffer)
lines(east_buffer)
lines(west_buffer)
lines(Ethiopia)

# calculate statistics through time of the models
model_stats <- c()
# columns of output = 
# 1. genus name 
# 2. area contemporary in Ethiopia
# 3. area contemporary in close region
# 4. area contemporary east
# 5. area contemporary west
# 6. mean area past in Ethiopia
# 7. mean area past in close region
# 8. mean area past east
# 9. mean area past west
for(a in 1:length(input_files)) {
  # 1 and 2
  a_rep <- c(species_names[[a]], 
              cellStats(area(reclass_contemporary[[a]]) * reclass_contemporary[[a]], 
              stat='sum'))
  # 3
  a_temp <- mask(crop(reclass_contemporary[[a]], grv_buffer), grv_buffer)
  a_rep <- c(a_rep, cellStats(area(a_temp) * a_temp, stat='sum'))
  # 4
  a_temp <- mask(crop(reclass_contemporary[[a]], east_buffer), east_buffer)
  a_rep <- c(a_rep, cellStats(area(a_temp) * a_temp, stat='sum'))
  # 5
  a_temp <- mask(crop(reclass_contemporary[[a]], west_buffer), west_buffer)
  a_rep <- c(a_rep, cellStats(area(a_temp) * a_temp, stat='sum'))
  # 6
  a_temp1 <- cellStats(area(reclass_ccsm[[a]]) * reclass_ccsm[[a]], stat='sum')
  a_temp2 <- cellStats(area(reclass_miroc[[a]]) * reclass_miroc[[a]], stat='sum')
  a_rep <- c(a_rep, mean(c(a_temp1, a_temp2)))
  # 7
  a_temp1 <- mask(crop(reclass_ccsm[[a]], grv_buffer), grv_buffer)
  a_temp1 <- cellStats(area(a_temp1) * a_temp1, stat='sum')
  a_temp2 <- mask(crop(reclass_miroc[[a]], grv_buffer), grv_buffer)
  a_temp2 <- cellStats(area(a_temp2) * a_temp2, stat='sum')
  a_rep <- c(a_rep, mean(c(a_temp1, a_temp2)))
  # 8
  a_temp1 <- mask(crop(reclass_ccsm[[a]], east_buffer), east_buffer)
  a_temp1 <- cellStats(area(a_temp1) * a_temp1, stat='sum')
  a_temp2 <- mask(crop(reclass_miroc[[a]], east_buffer), east_buffer)
  a_temp2 <- cellStats(area(a_temp2) * a_temp2, stat='sum')
  a_rep <- c(a_rep, mean(c(a_temp1, a_temp2)))
  # 9
  a_temp1 <- mask(crop(reclass_ccsm[[a]], west_buffer), west_buffer)
  a_temp1 <- cellStats(area(a_temp1) * a_temp1, stat='sum')
  a_temp2 <- mask(crop(reclass_miroc[[a]], west_buffer), west_buffer)
  a_temp2 <- cellStats(area(a_temp2) * a_temp2, stat='sum')
  a_rep <- c(a_rep, mean(c(a_temp1, a_temp2)))
  # add to total
  model_stats <- rbind(model_stats, a_rep)
}
model_stats <- data.frame(Species=as.character(model_stats[,1]),
                          Current_Ethiopia=as.numeric(model_stats[,2]),
                          Current_close=as.numeric(model_stats[,3]),
                          Current_east=as.numeric(model_stats[,4]),
                          Current_west=as.numeric(model_stats[,5]),
                          LGM_Ethiopia=as.numeric(model_stats[,6]),
                          LGM_close=as.numeric(model_stats[,7]),
                          LGM_east=as.numeric(model_stats[,8]),
                          LGM_west=as.numeric(model_stats[,9]))
write.table(model_stats, file="model_output.txt", quote=F, row.names=F, sep="\t")

# save
save.image("modeling.RData")
