### Preparing input matrices ###



#### Neighbor matrix ####

### Open Atlantic shapefile
library(raster)

setwd("C:\\Users\\Andre\\Google Drive\\SampSim\\Data\\Basic\\Shapes")
Atlantic <- shapefile("Atlantic.shp")


### Creat grid
r <- raster(extent(Atlantic))
res(r) <- 5
rAtlantic <- rasterize(Atlantic, r, 0, getCover=T) #5>25:392 / 1>0:9718 / 1>25:9360

plot(rAtlantic > 25)
plot(Atlantic, add=T, border="red")

cls <- which(rAtlantic[] <= 25) 
rAtlantic[][cls] <- NA

AtlGrid <- as(rAtlantic, "SpatialPolygonsDataFrame")
plot(AtlGrid)


### Create neighbor matrix
library(rgeos)
vizinhos <- gTouches(AtlGrid, byid=TRUE)
vizinhos[1:10,1:5]

vizinhos <- vizinhos*1
vizinhos[1:10,1:5]
#text(coordinates(AtlGrid)[,1], coordinates(AtlGrid)[,2], vizinhos[,1])


### Export neighbor matrix
setwd("C:\\Users\\Andre\\Google Drive\\RangeSimu\\Dados")
write.table(vizinhos, file="NeigCellMat.txt",row.names = F, col.names = F)



#### Latitude information ####

### Load Ophiuroids data and select occurrences records within the Atlantic
load("C:\\Users\\Andre\\Google Drive\\RangeSimu\\Dados\\Basic\\Ophiuroids.RData")

ophiPoints <- Ophiuroidea
coordinates(ophiPoints) <- ~decimalLongitude+decimalLatitude

ii <- which(!is.na(over(x = ophiPoints, Atlantic)))
ophiAtlant <- Ophiuroidea[ii,]
nrow(ophiAtlant) # 69165


### Calculate number of unique samplings events
unic <- UnicSamp(matriz = ophiAtlant)
unic <- unic[-which(unic==0)]


### Identify cells ID to each latitude
cellLat <- matrix(nrow =length(AtlGrid), ncol = 2)
for(i in 1:length(AtlGrid))
{
  cellLat[i,1] <- i-1
  cellLat[i,2] <- AtlGrid@polygons[[i]]@Polygons[[1]]@coords[3,2]
}
#plot(AtlGrid)
#text(coordinates(AtlGrid)[,1], coordinates(AtlGrid)[,2], cellLat[,2], cex=.5, col="blue")


# This matrix have 3 columns: latitude, ID of the 1° cell, total of cells
CellLatMat <- matrix(nrow = length(unique(cellLat[,2])), ncol = 2)
CellLatMat <- cbind(unique(cellLat[,2]), CellLatMat)
for(i in 1:length(unique(cellLat[,2])))
{
  x <- which(cellLat[,2]==unique(cellLat[,2])[i])
  CellLatMat[i,2] <- cellLat[,1][x][1]
  CellLatMat[i,3] <- length(x)
}


### Combine with sampling information
CellLatMat <- cbind(CellLatMat,rev(unic))


### Export data
setwd("C:\\Users\\Andre\\Google Drive\\RangeSimu\\Dados")
write.table(CellLatMat, file="CellLatMat.txt",row.names = F, col.names = F)

#### End













#### Records by sampling event ####

### Run the 'SampUnic' function below

### Total Ophiuroidea
x <- SampUnic(matriz = Ophiuroidea)
mean(x$mRec)


### Atlantic
y <- SampUnic(matriz = ophiAtlant)
mean(y$mRec, na.rm = T)


### Function to calculate average number of records by sampling event
SampUnic <- function(matriz)
{
  require('reshape')
  
  # Add the latitude intervals
  matriz[,7] <- cut(matriz[,3], breaks = seq(-80, 85, 5))
  lev <- levels(matriz[,7])
  
  # Start the progress bar
  pb <- txtProgressBar(min = 0, max = length(lev), style = 3)
  
  # Create the vectors to save the results
  UnicRec <- numeric()
  totRec <- numeric()
  mRec <- numeric()
  
  # For each latitude...
  for(i in 1:length(lev))
  {
    # Create a submatrix with the data of this latitudinal band
    z <- which(matriz[,7]==lev[i])
    if(length(z)==0)
    {
      UnicRec[i] <- NA
      totRec[i] <- NA
      mRec[i] <- NA
    }
    else
    {
      tempMat <- matriz[z,]
      tempMat2 <- cbind(tempMat[,1:4], count=rep(1, nrow(tempMat)))
      
      # Transform the submatrix in a community matrix
      CommMat <- cast(tempMat2, decimalLatitude + decimalLongitude + eventDate ~ species, value = "count")
      CommMat2 <- as.data.frame(CommMat)
      CommMat3 <- CommMat2[,-c(1,2,3)]
      CommMat3 <- as.matrix(CommMat3)
      CommMat3[which(is.na(CommMat3))] <- 0
      
      # Count how many sampling events record only one species
      x <- table(rowSums(CommMat3))
      pos <- which(names(x)=="1")
      
      if(length(pos)>0)
      {
        UnicRec[i] <- table(rowSums(CommMat3))[pos]
      }
      else
      {
        UnicRec[i] <- 0
      }
      
      # Calculate the mean and total number of records by sampling event at each latitude
      mRec[i] <- mean(rowSums(CommMat3))
      totRec[i] <- sum(rowSums(CommMat3))
    }
    setTxtProgressBar(pb, i)
  }  
  return(list(UnicRec=UnicRec,totRec=totRec,mRec=mRec))
}
