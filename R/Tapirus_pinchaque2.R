
#code from: http://blog.remote-sensing-conservation.org/species-distribution-modelling-to-show-potential-distribution-of-fagus-sylvatica-2013-and-prediction-for-2050-and-2080/

# install and load packages 

#install and load packages 

#install.packages
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(rms)
library(ellipse)
library(gam)
library(biomod2)
#setwd("C:/Users/Moses/Desktop/Project")

######################################
### Creating map of Europe

# Map of the World
#install.packages("maptools")
# library(maptools)
# worldmap<-readOGR("/worldmap/world_adm0","world_adm0") 
# 
# summary(worldmap)
# plot(worldmap, axes=T, col="lightgrey")
# 
# # Map of Europe
# Europe<-worldmap[worldmap$REGION=="Europe",]
# 
# # Europe without Russia
# Europe<-Europe[!Europe$NAME=="Russia",]
# 
# #plot:
# plot(Europe,axes=T,col="lightgrey")

window<- readShapePoly("data/Dissolved_Cliped.shp") #get the limit window
summary(window) # notice no projection, but see Dissolved_Cliped.shp
# lets define projection
geo<-CRS("+proj=longlat +datum=WGS84")
proj4string(window)<- geo
###############################################
### DOWNLOADING Fagus Sylvatica occurence data 

# if (file.exists("./GIS/Fagussylvatica/Fagussylvatica.mif")) {
#   fagus <- readOGR("./GIS/Fagussylvatica/Fagussylvatica.mif", layer = "Fagussylvatica")
# } else {
  # Download species location data from gbif
  Tpinchaque <- gbif("Tapirus", "pinchaque", sp=TRUE, geo=FALSE)
  
#   # Add projection information
#   proj4string(fagus) <- CRS("+proj=longlat +datum=WGS84")
#   
#   # Save species records in mif-format (preserves full column names)
#   writeOGR(sylvatica, "./Fagus present occurence/Fagussylvatica", 
#            "Fagussylvatica", driver="MapInfo File", dataset_options="FORMAT=MIF")
# }

proj4string(Tpinchaque) <- CRS("+proj=longlat +datum=WGS84")
plot(window)
plot(Tpinchaque,add=T)
head(Tpinchaque)
summary(Tpinchaque)
class(Tpinchaque)

fieldpoints_PerEcuaCol<-read.csv("C:/Mapas/Tapirus_pinchaque/PerEcuaCol.csv")

#filter field points
fieldpoints<-fieldpoints_PerEcuaCol[, 1:2]
colnames(fieldpoints)<-c("lat","lon")
gbiftapir<-cbind(Tpinchaque$lat, Tpinchaque$lon)
colnames(gbiftapir)<-c("lat","lon")
pinchaque<-rbind(fieldpoints,gbiftapir)
dup<-duplicated(pinchaque[,c('lon','lat')])
pinchaque.clean<-pinchaque[!dup,] #elimina duplicados
pinchaque.point<-pinchaque.clean
# change from data frame to Object of class SpatialPoints
coordinates(pinchaque.point)<-c("lat","lon") # 
proj4string(pinchaque.point)<- latlon

# #select the records that have longitude and latitude data
# colnames(Tpinchaque)
# CORANT <- subset(Tpinchaque, institution=="Corantioquia") # elimina Corantioquia 
# dim(Tpinchaque) # to get the number columns of the raster object sylvatica with long and lat information
# ####

##########################################
### Make a map of the occurrence localities of Fagus Sylvatica
## --> it is important to make such maps to assure that the points are, at least
## roughly, in the right location.
#' ###Distribution of Fagus Sylvatica From GBIF###
### plot points
#x11()
plot(window,axes=T,col="lightgrey")
points(pinchaque.clean$lon, pinchaque.clean$lat, col="red",pch=".", cex=3)

#The wrld_simpl dataset contains rough country outlines. You can use other
#datasets of polygons (or lines or points) as well. For example, you can download
#higher resolution data country and subnational administrative boundaries data
#with the getData function of the raster package. You can also read your own
#shapefile data into R using the shapefile function in the raster package


#########################################
### Download climate element data 
### Here we make a list of worldclim climate and then create a rasterStack from 
### these, show the names of each layer, and finally plot them all.   
### http://www.worldclim.org/formats

library(raster)
#tmin <- raster::getData("worldclim", var="tmin", res=10)
tmin.tile1<-getData("worldclim",var="tmin",res=10,lon=-81, lat=6)
tmin.tile2<-getData("worldclim",var="tmin",res=10,lon=-81, lat=-6)
tmin<-merge(tmin.tile1, tmin.tile2)
names(tmin)<-names(tmin.tile1)
# res= 10 -->this will download global data on minimum temperature at 10 min resolution;
# Temperature data is in units of °C * 10 because that allows worldclim to store the data as Integer
# values which is more efficient than storing the data as Real values

# plot(tmin.tile1, ext=extent(window))
# names(tmin.tile1)

#tmax <- raster::getData("worldclim", var="tmax", res=10)
tmax.tile1<-getData("worldclim",var="tmax",res=10,lon=-81, lat=6)
tmax.tile2<-getData("worldclim",var="tmax",res=10,lon=-81, lat=-6)
tmax<-merge(tmax.tile1, tmax.tile2)
names(tmax)<-names(tmax.tile1)

#prec <- raster::getData("worldclim", var="prec", res=10)
prec.tile1<-getData("worldclim",var="prec",res=10,lon=-81, lat=6)
prec.tile2<-getData("worldclim",var="prec",res=10,lon=-81, lat=-6)
prec<-merge(prec.tile1, prec.tile2)
names(prec)<-names(prec.tile1)

plot(prec,ext=extent(window))
#head(prec)


###################
# get vegetation map
# donwload a copy of the Global Land Cover 2000 map:
library(GSIF)
# see http://worldgrids.org/doku.php?id=wiki:functions
URI = "http://wps.worldgrids.org/pywps.cgi"
server <- list(URI=URI, request="execute", version="version=1.0.0", service.name="service=wps",
                 +   identifier="identifier=sampler_local1pt_nogml")
biocl15.wps <- new("WPS", server=server, inRastername="biocl15")
str(biocl15.wps)
##############################
# predictors

predictors <- stack(tmin, tmax, prec)
#
plot(predictors,1 , ext=extent(window))
points(pinchaque.clean$lon,pinchaque.clean$lat,col="red",pch=".", cex=3)
#
#plot(predictors,ext=extent(Europe))
#
#plot(raster(predictors, 1))
#plot(sylvatica, add = TRUE)
#

names(predictors)  # lists the different predictor layers that are selected
#' Select species records for which environmental information is available
#' -------------------------------
T_pinchaque <- pinchaque.point[complete.cases(extract(predictors, pinchaque.point)), ]

#' Collinearity
#' -----------------------------
#' ### Visual inspection of collinearity ###
cm <- cor(getValues(predictors), use = "complete.obs")
plotcorr(cm, col=ifelse(abs(cm) > 0.7, "red", "grey"))
#
#' ### Select an uncorrelated subset of environmental variables ###
varnames <- c("tmin1", "tmax7", "prec1","prec7")
env <- subset(predictors, varnames)
#
#' Sampling of (pseudo-)absence points
#' ====================================================
#' The function randomPoints in package dismo allows to 
#' randomly select a certain number of random points,
#' and to adjust the probability of selecting a cell
#' according to its size, which is relevant in lat-lon-grids,
#' where cells are of differing size

#' Selecting 2000 random background points, excluding cells where
#' the species is present
set.seed(2)
background <- randomPoints(env, 2000, pinchaque.point)
#' Select only one presence record in each cell of the environmental layer
presence <- gridSample(pinchaque.point, env, n = 1)
#' 
#' Now we combine the presence and background points, adding a 
#' column "species" that contains the information about presence (1)
#' and background (0)
fulldata <- SpatialPointsDataFrame(rbind(presence, background),
                                   data = data.frame("pinchaque.point" = rep(c(1,0), 
                                   c(nrow(presence), nrow(background)))),
                                   match.ID = FALSE,
                                   proj4string = CRS(projection(env)))
#' Add information of environmental conditions at point locations
fulldata@data <- cbind(fulldata@data, extract(env, fulldata))
#' 
# Split data set into a training and test data set
set.seed(2)
fold <- kfold(fulldata, k = 5)
traindata <- fulldata[fold != 1, ]
testdata <- fulldata[fold == 1, ]
#' We can now use a range of statistical methods to estimate the
#' probability of species occurrence.
#' Unfortunately, there are often subtle differences in how the models
#' are specified and in which data formats are useable
library(gam)
## Generalized additive models
gammodel <- gam(pinchaque.point ~ s(tmin1) + s(tmax7) + s(prec1) + s(prec7),
                family="binomial", data=traindata)

summary(gammodel)
#
# Now we should do model selection: bio14 does not contribute to the fit

# Evaluate model on test data
# a) Predict to test data
gamtest <- predict(gammodel, newdata = testdata, type = "response")
# b) Calculate performance indices #does not work?
#  val.prob(gamtest, testdata[["pinchaque.point"]])
#
# Variable importance work using random forest???????????????????????????//

# gamimp <- varImpBiomod(gammodel, varnames,
#                        traindata)
# barplot(100 * gamimp/sum(gamimp), ylab = "Variable importance (%)")

# Response functions
plot(gammodel, pages = 1)
# Prediction map
gammap <- predict(env, gammodel, type = "response")
plot(gammap,ext=extent(window)) # ?????????
######
## Random forest
library(randomForest)
# randomForest requires the dependent variable to be a factor
# if we want to do classification
rftraindata <- as(traindata, "data.frame")
rftraindata$pinchaque.point <- factor(rftraindata$pinchaque.point)
#rm(rftraindata$sylvatica)
# TODO: check proper settings of random forest algorithm
rftraindata <- rftraindata[complete.cases(rftraindata), ]
rfmodel <- randomForest(pinchaque.point ~ tmin1 + tmax7 + prec1 + prec7, data = rftraindata)

# Evaluate model on test data
# a) Predict to test data
rftest <- predict(rfmodel, newdata = testdata, type = "prob")[,2]
# b) Calculate performance indices
val.prob(rftest, testdata[["pinchaque.point"]])

# Variable importance
rfImp <- importance(rfmodel)
varImpPlot(rfmodel)

# Response functions
par(mfrow=c(3,2))
for (i in seq_along(varnames)) {
  partialPlot(rfmodel, rftraindata, varnames[i], xlab = varnames[i], main="")  
}
# Prediction map
rfmap <- predict(env, rfmodel, type = "prob", index = 2)
par(mfrow=c(1, 1))
plot(rfmap,ext=extent(Europe))

## Maxent
# The following code assumes that the column with the species information
# is in the first position
library(dismo)
maxentmodel <- maxent(traindata@data[, -1], traindata[["pinchaque.point"]], 
                      args = c("nothreshold", 
                               "nohinge"))

# Model evaluation on test data
maxenttest <- predict(maxentmodel, testdata)
#val.prob(maxenttest, testdata[["sylvatica"]])

# Alternatively, we can use the evaluate function
maxente <- evaluate(p = maxenttest[testdata[["pinchaque.point"]] == 1],
                    a = maxenttest[testdata[["pinchaque.point"]] == 0])

# Show variable importance
plot(maxentmodel)

# Plot response functions
response(maxentmodel)

# Prediction map
maxentmap <- predict(maxentmodel, env)
plot(maxentmap,ext=extent(Europe))


#' ###Prediction for 2050 #####
setwd("./wc2050")

varnames2050 <- c("tmin_1", "tmax_7", "prec_1","prec_7")
predictors2050 <- stack(paste(varnames2050, ".asc", sep=""))
plot(predictors2050)

env2050 <- subset(predictors2050, varnames2050)
# to rename predictors
names(env2050) <- gsub("_", "", names(env2050))
env2050 <- crop(env2050, extent(Europe))
plot(env2050)
#'### gam prediction 2050 ######
gampred2050 <- predict(env2050, gammodel, type = "response")
plot(gampred2050,main="GAM Prediction 2050")

#'###RandomForest Prediction 2050 ####
rfpred2050 <- 1 - predict(env2050, rfmodel, type = "prob")
plot(rfpred2050,main="Random Forest Prediction 2050")
#' ### maxent prediction 2050 ###
#' 
maxent2050<-predict(env2050, maxentmodel,type="prob")
plot(maxent2050,main="Maxent Prediction 2050")


#-----------------2080-----------------------#
#' ###### Predictions for 2080#####
setwd("./wc2080")
varnames2080 <- c("tmin_1", "tmax_7", "prec_1","prec_7")
predictors2080 <- stack(paste(varnames2080, ".asc", sep=""))
plot(predictors2080)

env2080 <- subset(predictors2080, varnames2080)
# to rename predictors
names(env2080) <- gsub("_", "", names(env2080))
env2080 <- crop(env2080, extent(Europe))
plot(env2080)
#'### gam prediction 2080 ####
gampred2080 <- predict(env2080, gammodel, type = "response")
plot(gampred2080,main="GAM Prediction 2080")
#' ##### Random forest Prediction 2080####
rfpred2080 <-1- predict(env2080, rfmodel, type = "prob")
plot(rfpred2080,main="Random Forest Prediction 2080")
#'###maxent prediction 2080##
maxent2080<-predict(env2080, maxentmodel,type="prob")
plot(maxent2080,main="Maxent Prediction 2080")
# plots comparing same type of models#
par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
brks <- seq(0, 1, by = 0.1)
arg <- list(at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))
col <- rev(terrain.colors(length(brks) - 1))
plot(gampred2050, breaks = brks, col = col, axis.args = arg,main="GAM prediction 2050")
plot(gampred2080, breaks = brks, col = col, axis.args = arg,main="GAM prediction 2080")
north.arrow(xb=80, yb=65, len=3, lab="N")
plot(rfpred2050, breaks = brks, col = col, axis.args = arg,main="Random Forest Prediction 2050")
plot(rfpred2080, breaks = brks, col = col, axis.args = arg,main="Random Forest Prediction 2080")
plot(maxent2050, breaks = brks, col = col, axis.args = arg,main="Maxent Prediction 2050")
plot(maxent2080, breaks = brks, col = col, axis.args = arg,main="Maxent Prediction 2080")

##################################################
## comparing  models###
#2050#
par(mfrow = c(3, 1), mar = c(3, 3, 1, 1))
brks <- seq(0, 1, by = 0.1)
arg <- list(at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))
col <- rev(terrain.colors(length(brks) - 1))
plot(gampred2080, breaks = brks, col = col, axis.args = arg,main="GAM prediction 2050")
north.arrow(xb=100, yb=65, len=3, lab="N")
plot(rfpred2080, breaks = brks, col = col, axis.args = arg,main="Random Forest Prediction 2050")
plot(maxent2080, breaks = brks, col = col, axis.args = arg,main="Maxent Prediction 2050")

# 2080#
par(mfrow = c(3, 1), mar = c(3, 3, 1, 1))
brks <- seq(0, 1, by = 0.1)
arg <- list(at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))
col <- rev(terrain.colors(length(brks) - 1))
plot(gampred2080, breaks = brks, col = col, axis.args = arg,main="GAM prediction 2080")
north.arrow(xb=100, yb=65, len=3, lab="N")
plot(rfpred2080, breaks = brks, col = col, axis.args = arg,main="Random Forest Prediction 2080")
plot(maxent2080, breaks = brks, col = col, axis.args = arg,main="Maxent Prediction 2080")