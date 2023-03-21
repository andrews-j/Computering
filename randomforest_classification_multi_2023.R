#############################################################################
# The script reads an ESRI Shapefile (defined by the "shapefile" variable) with 
# training polygons and then either select all pixels or randomly select a 
# user-determined number of samples (defied by the "numsamps" variable) from 
# each land cover type. A multilayer image that contains spectral, other 
# continuous data or categorical data is also input (defined by the inImage 
# variable). For each randomly selected sample the data values for that pixel 
# are determined and these data are used to run the Random Forest model. 
#
# After building the model the multilayer image is read, and up to three output 
# layers (classImage, probImage, threshImage) can be selected for the output image. 
#      "classImage" classifies all of the pixels.
#
#     "probImage" outputs the class probability of the class that got the most votes 
#      (i.e., the class that was selected for the classImage layer). 
#
#     "threshImage" is the same as "classImage" except all pixels with a class probability 
#      of the class that got the most votes below the "probThreshold" parameter are set to 0. 
#      This is useful to identify pixels with inter-class confusion.
#
# The image is written out (name and location is defined by the "outImage variable) 
# using the GeoTIFF format. A variable importance plot is displayed to provide information 
# about the influence of each variable. An error rate estimate and confusion matrix are also 
# printed to provide information about classification accuracy.
#
# There is an option to assess the quality of the training data. The metric for this 
# is the margin. The margin of a training point is the proportion of votes for the correct 
# class minus maximum proportion of votes for the other classes for that segment. Positive margin 
# values represent correct classification, and vice versa. The margin data are written to a 
# point ESRI Shapefile so they can be overlaid on the image and training polygons to assess which 
# points need to be removed and relabeled in the training data and it can help determine which 
# classes needs additional training segments. If this output is not needed you can enter two 
# double or single-quotes ("" or '') for the variable outPointsFile.
#
# There is also an option to output a feature space plot using two bands of your choice.
# If a feature space plot is not needed then enter "0" for the variables xBand and/or yBand.
# 
# Set the variables below in the "SET VARIABLES HERE" section of the script. 
#
# The original script was written by Ned Horning [horning@amnh.org]
# Support for writing and maintaining this script comes from The John D. and 
# Catherine T. MacArthur Foundation and Google.org.
#
# The script has been altered by Sean Griffin
#
# This script is free software; you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software Foundation
# either version 2 of the Licenase, or ( at your option ) any later version.                               *
# ******************************************************************************************
# In August 2014, changes were made to the original script allowing these improvements:
# 1.inclusion of the rgdal library to read in the spatial reference from the Shapefile
#   allowing use of different projections
# 2.creating an ENVI Classification directly from the outputs. This change requires a small
#   CSV file to be referenced which contains the color table and classes of the outputs. An 
#   example of its contents is shown below:     
#     Unclassified,0,0,0
#     Soybean,255,255,0
#     Corn,154,0,0
#     Fallow,132,121,66
#     Urban,132,124,67
#     Water,0,0,159
#     Rice,0,103,0
#     Cloud,255,255,255
#     Forest,0,255,0
#     Forest slope,0,255,0
#     Riparian,90,150,150
#     Edge,0,0,0
#   To enable this output, provide the full path to the variable, "legend", otherwise,
#   keep it as a blank string, (i.e. "").
#
#
#############################################################################
#Load libraries

install.packages("maptools")
install.packages("sp")
install.packages("randomForest")
install.packages("raster")
install.packages("rgdal")

require(maptools)
require(sp)
require(randomForest)
require(raster)
require(rgdal)
#
cat("Set variables and start processing\n")
#
#############################   SET VARIABLES HERE  ###################################
# Name and path for the input satellite image for Training 
filepath='C:/Users/jgand/OneDrive - Clark University/Advanced RS/Lab 5 Classification/Data/InputData'
inImageTrain <-'C:/Users/jgand/OneDrive - Clark University/Advanced RS/Lab 5 Classification/Data/InputData/Spectral_LT7_19990707.tif'
# Name and path for the Shapefile
shapefile    <-'C:/Users/jgand/OneDrive - Clark University/Advanced RS/Lab 5 Classification/Data/CalibrationValidation/calibrationPolys_Subset.shp'
# Name of the attribute that holds the integer land cover type identifyer
attName <- 'GRIDCODE'
# No-data value for the input image
nd <- -9999
# Approximate number of training samples to be randomly selected for each land cover class
# If numsamps is set to "0" then all pixels in all of the polygons will be used as training samples
numsamps <- 50
# Name and path for the input satellite image for classifying
inImageClass <-'C:/Users/jgand/OneDrive - Clark University/Advanced RS/Lab 5 Classification/Data/InputData/Spectral_LT7_19990707.tif'
# Name and path for the names of each band in the input image composite for replacing the original column names
inImageBNames <-'C:/Users/jgand/OneDrive - Clark University/Advanced RS/Lab 5 Classification/Data/InputData/inputtable.csv'

# Name and path of the output image - WITH NO EXTENSION!
outImageName <-'C:/Users/jgand/OneDrive - Clark University/Advanced RS/Lab 5 Classification/Data/Output/Spectral_out_image'
outFormat <- 'GTiff'
# CSV file containing class names and color map to create ENVI Classification output header, 
# if NULL or '', then header of output file is left alone
legend <- ''
########## You should not have to change any of the parameters below #################
# Name and location of the output Shapefile point file that will be created. If this output 
# is not needed you can enter two double or single-quotes ("" or '')
# Note that if this file exists the write will fail with the message "Creation of output file failed"  
outMarginFile <- 'C:/Users/jgand/OneDrive - Clark University/Advanced RS/Lab 5 Classification/Data/Output/Spectral_out_margin'
# Output classification layer without applying threshold (enter TRUE or FALSE)
classImage <- TRUE
# Output probability image layer (enter TRUE or FALSE)
probImage <- TRUE
# Output classification layer and set pixels with probability less than "probThreshold" to 0 (enter TRUE or FALSE)
threshImage <- TRUE
# Enter threshold probability in percent (values must be between 0 and 100) only used if threshImage=TRUE
probThreshold <- 75
# Layer number (band number) for the X and Y axis of the feature space plot. 
# If you do not want to calculate a feature plot enter 0 as the layer number
xBand <- 0
yBand <- 0
###########################################################
#############################
# Start processing
startTime <- Sys.time()
cat("Start time", format(startTime),"\n")

# Read the Shapefile using rgdal (better than maptools for projections)
tmp <- basename(shapefile)
tmp <- unlist(strsplit(tmp,".shp"))
vec <- readOGR(dsn=dirname(shapefile),layer=tmp)

# import table with the variable names
bnames <- read.csv(file = inImageBNames)
bnames_spectral<- bnames[7:12,]
#
# Load the training image then flag all no-data values (nd) so they are not processed
satImageTrainall <- stack(inImageTrain)
#satImageTrain <-satImageTrainall[[7:12]]  # only use the spectral bands as input, try other combinations as well
satImageTrain <-satImageTrainall
# replace variable names
names(satImageTrain) <- c(bnames_spectral$Name)
for (b in 1:nlayers(satImageTrain)) { NAvalue(satImageTrain@layers[[b]]) <- nd }

# Load the classifying image then flag all no-data values (nd) so they are not processed
satImageClassall <- stack(inImageClass)
#satImageClass <- satImageClassall[[7:12]] # only use the spectral bands as input, try other combinations as well 
satImageClass <- satImageClassall
# replace variable names
names(satImageClass) <- c(bnames_spectral$Name)
for (b in 1:nlayers(satImageClass)) { NAvalue(satImageClass@layers[[b]]) <- nd }
#


# Create vector of unique land cover attribute values
allAtt <- slot(vec, "data")
tabAtt <-table(allAtt[[attName]])
uniqueAtt <-as.numeric(names(tabAtt))

# Check if there are data in uniqueAtt
if (is.na(uniqueAtt[1])) {
  cat("\n*************************No attributes were found**************************** \n")
  stop("Check the attName variable in the variable settings\n", call.=FALSE)
}
# If all pixels in a polygon are to be used process this block
if (numsamps == 0) {
  # Create input data from a Shapefile using all training data 
  cat("Create training data using all pixels in training polygons\n")
  predictors <- data.frame()
  response <- numeric()
  for (x in 1:length(uniqueAtt)) {
    # Get the metadata for all polygons for a particular class (based on the uniqueAtt variable)
    class_data<- vec[vec[[attName]]==uniqueAtt[x],]
    # Extract and combine predictor and response variables for each polygon within a class
    for (i in 1:dim(class_data)[1]) {
      satValues <- extract(satImageTrain, class_data[i,])
      satValues <- as.data.frame(do.call(rbind,satValues))
      attributeVector <- rep.int(uniqueAtt[x],nrow(satValues))
      
      predictors <- rbind(predictors, satValues)
      response <- c(response, attributeVector)
    }
  }
  trainvals <- cbind(response, predictors)
} else {
  # Create input data from a Shapefile by sampling training data 
  cat("Create training data by sampling", numsamps, "pixels for each class\n")
  for (x in 1:length(uniqueAtt)) {
    # Get the metadata for all polygons for a particular class (based on the uniqueAtt variable)
    class_data<- vec[vec[[attName]]==uniqueAtt[x],]
    # Get the area of each polygon for a particular class
    areas <- sapply(slot(class_data, "polygons"), slot, "area")
    # Calculate the number of samples for each polygon based on the area in proportion to total area for a class
    nsamps <- ceiling(numsamps*(areas/sum(areas)))
    # Use random sampling to select training points (proportial based on area) from each polygon for a given class 
    for (i in 1:dim(class_data)[1]) {
      xy_class <- spsample(class_data[i,], type="random", n=nsamps[i], iter=10)
      # Add coordinates to create a list of random points for all polygons
      if (i == 1) cpts <- xy_class
      else cpts <- rbind(cpts, xy_class)
    }
    # The number of points might not match numsamps exactly.
    classpts <- cpts
    if (x == 1) {
      xy_allClasses<- classpts
    } else {
      xy_allClasses<- rbind(xy_allClasses, classpts)
    }
  } 
  # Get class number for each sample point for responce variable
  response <- over(xy_allClasses, vec)[[attName]]
  # Get pixel DNs from the image for each sample point
  trainvals <- cbind(response, as.data.frame(extract(satImageTrain, xy_allClasses))) # SY modified
}
# Test if feature space plot is needed
if (xBand != 0 & yBand != 0) {
  #Plot feature space and samples
  #
  continue <- "c"
  while (continue == "c") {
    plotImage <- stack(satImageTrain[[xBand]], satImageTrain[[yBand]])
    # Get pixel values from the image under each sample point and create a table with 
    # observed and predicted values
    cat("Getting pixel values to create feature space plot\n\n")
    featurePlotPoints <- sampleRegular(plotImage,100000 )
  
    # Remove NA values from trainvals table created above
    featurePlotPoints <- na.omit(featurePlotPoints)
  
    minBand1 <- min(featurePlotPoints[,1])
    maxBand1 <- max(featurePlotPoints[,1])
    minBand2 <- min(featurePlotPoints[,2])
    maxBand2 <- max(featurePlotPoints[,2])
    rangeBand1 <- maxBand1 - minBand1 + 1
    rangeBand2 <- maxBand2 - minBand2 + 1
  
    xAxisLabel <- paste("Layer", xBand, sep=" ")
    yAxisLabel <- paste("Layer", yBand, sep=" ")
  
    plot(featurePlotPoints[,1], featurePlotPoints[,2], col="lightgrey", xlab=xAxisLabel, ylab=yAxisLabel)
  
    uniqueValues <- unique(trainvals[,1])
    for (v in 1:length(uniqueValues)) {
      points(trainvals[which(trainvals[,1]==uniqueValues[v]), xBand+1], trainvals[which(trainvals[,1]==uniqueValues[v]), yBand+1], col=v, pch=20)
    }
  
    legend(minBand1, maxBand2, col=1:v, pch=20, title="Classes", legend=as.character(uniqueValues))
  
    continue <- readline(prompt="Type n to stop, c to change feature space bands or any other key to continue with randome forests model creation and prediciton: \n\n")
  
    if (substr(continue, 1,1) == "n") {
      stop("Processing stopped at users request \n\n", call.=FALSE)
    }
    if (substr(continue, 1,1) == "c") {
      xBand <- as.numeric(readline(prompt="Enter the band number for the x axis: \n"))
      yBand <- as.numeric(readline(prompt="Enter the band number for the y axis: \n"))
    }
  }
}

# Remove NA values 
trainvals <- na.omit(trainvals)


# Run Random Forest
cat("Calculating random forest object\n")
randfor <- randomForest(as.factor(response) ~., data=trainvals, importance=TRUE, na.action=na.omit)
# Start predictions
cat("Starting predictions\n")




### Partial Dependence Plot ###
## Specify which class you want to focus on 
class <- 1
imp <- importance(randfor)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)][1:6]
op <- par(mfrow=c(2,3))

dev.off()

for (i in seq_along(impvar)){
  partialPlot(randfor, trainvals,impvar[i],class, xlab=impvar[i],
              main=paste("Partial Dependence on",impvar[i])
              )
}
mtext(paste("Partial Dependence Plots for Class ", class),side=3,line=-18,outer=TRUE)




# Calculate how many bands the output image should have
numBands <- classImage + probImage + threshImage
#numBands <- classImage # SY modified
# Calculate the image block size for processing
bs <- blockSize(satImageClass)

# Create the output raster block
outImage <- brick(satImageClass, values=FALSE, nl=numBands)
outImage <- writeStart(outImage, filename=outImageName, progress='text', format=outFormat, datatype='INT1U', overwrite=TRUE)

# "Trick" the random forest object to think its terms are the bands from
# image to be classified
names(satImageClass) = names(satImageTrain)

# Loop though each of the image blocks to calculate the output layers selected in the variables section
for (i in 1:bs$n) {
  cat("processing block", i, "of", bs$n, "\r")
  imageBlock <-  getValuesBlock(satImageClass, row=bs$row[i], nrows=bs$nrows[i])
  predValues <- predict(randfor, imageBlock, type='response')
  classValues <- as.numeric(levels(predValues))[predValues]
  outMatrix <- matrix(nrow=nrow(imageBlock), ncol=0)
  if (classImage) {
    outMatrix <- cbind(outMatrix, classValues)
  }
  if (probImage || threshImage) { 
    predProbs <- as.data.frame(predict(randfor, imageBlock, type='prob'))
    maxProb <- round(apply(predProbs, 1, max) * 100)
    if (probImage) { 
      outMatrix <- cbind(outMatrix, maxProb)
    }
    if (threshImage) {
      threshValues <- classValues
      threshValues[which(maxProb <= probThreshold)] <- 0
      outMatrix <- cbind(outMatrix,threshValues)
    }
  }
  writeValues(outImage, outMatrix, bs$row[i])
}

# Stop writing and close the file
outImage <- writeStop(outImage)

# Now that the image has been written, take out the first band of the output - the classification itself -
# and write to its own ENVI format file
# dimensions of output image
cols <- ncol(outImage)
rows <- nrow(outImage)
# first, read in the full image that was just produced
gdalInput <- paste(outImageName,'.tif',sep='')
result = readGDAL(gdalInput,band=1)
# swap out the filename extension with *class.dat
className <- gsub('.envi','class.dat',gdalInput)
classHeader <- gsub('.dat','.hdr',className)
# write out an ENVI format file as 8-bit unsigned byte and using 255 for flags, which
# seems to be compulsory
writeGDAL(result,className,drivername=outFormat,type='Byte', mvFlag=255)

## convert this ENVI format file into an ENVI classification
if (file.exists(legend) == T){
  # read in new header
  headertxt <- readLines(classHeader, n = -1L)
  numLines <- length(headertxt)
  
  # change the ENVI file type from Standard to Classification
  for (i in 1:length(headertxt)){ 
    headertxt[i] = gsub('Standard','Classification',headertxt[i])
  }  
  
  # read in legend color table as CSV
  legendtxt <- read.csv(legend, header=F, strip.white = T)  
  #create list of class names
  classnames = paste('class names = {',legendtxt[1,1], sep='')
  for (i in 2:nrow(legendtxt)){
    classnames = paste(classnames,legendtxt[i,1],sep = ', ')
  }
  classnames = paste(classnames,'}',sep = '')  
  
  # create color map as vector
  cmap = as.vector(aperm(as.matrix(legendtxt[,2:4])))
  # make permuted transpose 
  classlookup = paste('class lookup = {',as.character(cmap[1]),sep = '')
  for (i in 2:length(cmap)){
    classlookup = paste(classlookup,as.character(cmap[i]),sep = ', ')
  }
  classlookup = paste(classlookup,'}',sep = '')
  
  # assign number of classes to next line of header
  classesnum = paste('classes = ',as.character(nrow(legendtxt)))
  headertxt[numLines+1] = classesnum
  # assign class names to next line of header
  headertxt[numLines+2] = classnames
  # assign class colors to next line of header
  headertxt[numLines+3] = classlookup
  # create header file
  writeLines(headertxt,classHeader)
} # end of changing ENVI header file

# Plotting variable importance plot
varImpPlot(randfor)

# Print error rate and confusion matrix for this classification
confMatrix <- randfor$confusion
cat("#################################################################################\n")
cat("OOB error rate estimate\n", 1 - (sum(diag(confMatrix)) / sum(confMatrix[,1:ncol(confMatrix)-1])), "%\n\n", sep="")
cat("Confusion matrix\n")
print(randfor$confusion)
cat("\n")

if (outMarginFile != "") {
  # Calculate margin (proportion of votes for correct class minus maximum proportion of votes for other classes)
  marginData <- margin(randfor)
  trainingAccuracy <- cbind(marginData[order(marginData)], trainvals[order(marginData),1])
  
  # Add column names to attributes table
  colnames(trainingAccuracy) <- c("margin", "classNum")  
  # Calculate X and Y coordinates for training data points
  xyCoords <- xy_allClasses@coords
  xyCoords <- xyCoords[order(marginData),]
  
  # Create and write point Shapefile with margin information to help improve training data
  pointVector <- SpatialPointsDataFrame(xyCoords, as.data.frame(trainingAccuracy), coords.nrs = numeric(0), proj4string = satImageClass@crs)
  writeOGR(pointVector, outMarginFile, "layer", driver="ESRI Shapefile", check_exists=TRUE)
}

# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("\nProcessing time", format(timeDiff), "\n")


