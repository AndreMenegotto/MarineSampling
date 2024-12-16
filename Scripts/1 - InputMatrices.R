
#### Preparing input matrices ####
library(terra)

### Open Atlantic shapefile
load('./Data/Empirical/Atlantic.RData')
Atlantic <- unwrap(w_Atlantic)
AtlGrid <- unwrap(w_AtlGrid)



#### Neighbor matrix -> NeigCellMat ####

### Reduce Atlantic size a bit to avoid creating unrealistic spatial connections
b <- buffer(Atlantic, width=-7000)
x <- intersect(AtlGrid, b)
length(x)==length(AtlGrid)


### Get neighbor matrix ('adjacent' does not work good for this grid)
vizinhos <- nearby(x, distance=.01, centroids=F, symmetrical=F)
vizinhos <- data.frame(vizinhos, Value=1)
vizinhos <- with(vizinhos, tapply(Value, list(from, to), FUN = identity))
vizinhos <- ifelse(is.na(vizinhos),0,1)
#plot(AtlGrid, col='grey90')
#text(crds(centroids(AtlGrid))[,1], crds(centroids(AtlGrid))[,2], vizinhos[,1], cex=.5)


### Export neighbor matrix
#write.table(vizinhos, file="./Data/SimInput/NeigCellMat.txt",row.names = F, col.names = F)



#### Latitudinal information -> CellLatMat ####

### Load Ophiuroids data and select occurrence records within the Atlantic
# This requires the raw data set of occurrence records
#load("./Data/Empirical/Raw/Ophiuroids.RData")
load("./Data/Empirical/Functions.RData")

### Calculate the number of unique sampling events
unic <- UnicSamp(matriz = ophiAtlant)
unic <- unic[-which(unic==0)]


### Associate cells ID to each latitude
cellLat <- matrix(nrow = length(AtlGrid), ncol = 2)
for(i in 1:length(AtlGrid))
{
  # ID starts in 0 (zero) to be compatible with Object Pascal language
  cellLat[i,1] <- i-1
  cellLat[i,2] <- crds(centroids(AtlGrid[i,]))[,2] #AtlGrid@polygons[[i]]@Polygons[[1]]@coords[3,2]
}
cellLat[,2] <- round(cellLat[,2], 1)
#plot(AtlGrid)
#text(crds(centroids(AtlGrid))[,1], crds(centroids(AtlGrid))[,2], cellLat[,2], cex=.5, col="blue")


# This matrix has 3 columns: latitude, ID of the first cell, total number of cells
CellLatMat <- matrix(nrow = length(unique(cellLat[,2])), ncol = 2)
CellLatMat <- cbind(unique(cellLat[,2]), CellLatMat)
for(i in 1:length(unique(cellLat[,2])))
{
  x <- which(cellLat[,2]==CellLatMat[i,1])
  CellLatMat[i,2] <- cellLat[,1][x][1]
  CellLatMat[i,3] <- length(x)
}


### Combine with sampling information
CellLatMat <- cbind(CellLatMat,unic)


### Stop here and export data if the number of specimens per sampling event will be fixed
#write.table(CellLatMat, file="./Data/SimInput/CellLatMat.txt",row.names = F, col.names = F)



#### Records per sampling event ####

### Find the parameters of the beta distribution that better replicate the empirical number of records per sampling event
matriz <- ophiAtlant
matriz[,9] <- cut(matriz[,3], breaks = seq(-65, 65, 5))
lev <- levels(matriz[,9])

# Define the latitude
j=1
z <- which(matriz[,9]==lev[j])
tempMat <- matriz[z,]
nEvents <- sum(duplicated(tempMat[,2:4])==F)
nRecs <- nrow(tempMat)

tempMat2 <- cbind(tempMat[,1:4], count=rep(1, nrow(tempMat)))
spBySamp <- as.vector(table(paste(tempMat2[,2], tempMat2[,3], tempMat2[,4], sep = '_')))
m <- max(spBySamp)

# Find probable parameters, if possible, and adjust manually
#qtl <- quantile(spBySamp/max(spBySamp), probs = c(.025,.5,.975))
#ab <- get.beta.par(q = qtl, plot = F)
#round(ab/m, 2)

{
  frame()
  plot.window(xlim=c(0,40), ylim=c(1,2))
  a <- boxplot(spBySamp, horizontal = T, add=T, col="deepskyblue", outpch=21, outbg="deepskyblue", at=1.25)
  
  recs <- numeric()
  alfa <- 1/100
  for(i in 1:100)
  {
    x <- rbeta(nEvents, .63,  4)*m
    x <- ceiling(x)
    recs[i] <- sum(x)
    
    b <- boxplot(x, horizontal = F, at=1.75, add=T, frame.plot=F, axes=F, col="deepskyblue", outbg="deepskyblue", plot=F)
    boxplot(x, horizontal = T, at=1.75, add=T, frame.plot=F, axes=F, outline=F, col=rgb(red = .5, green = .5, blue = .5, alpha = alfa), border=rgb(red = 0, green = 0, blue = 0, alpha = alfa), medcol=rgb(red = 1, green = 0, blue = 0, alpha = alfa))
    points(x=unique(b$out), y=rep(1.75, length(unique(b$out))), pch=21, bg=rgb(red = .5, green = .5, blue = .5, alpha = alfa), col=rgb(red = 0, green = 0, blue = 0, alpha = alfa))
  }
  print(paste('Empiric: ',nRecs,sep=''))
  print(paste('Simulated: ',mean(recs),sep=''),)
}


### Reference matrix to replicate the number of species sampled at each sampling event
# alpha, beta, multiplier (max records per sampling event)
RefMat <- rbind(
  c(.63, 4,  18), #1
  c(.36, 2,  22), #2
  c(.35, 1.5, 16),#3
  c(.7,  3.5, 6), #4
  c(.6,  7,  13), #5
  c(.85, 4,  5),  #6
  c(.18, 1,  4),  #7
  c(.17, 1,  8),  #8
  c(.69, 5,  12), #9
  c(.63, 2,  9),  #10
  c(.5,  3,  7),  #11
  c(.7,  3,  5),  #12
  c(.75, 5,  8),  #13
  c(.1, .6,  4),  #14
  c(.75, 5,  13), #15
  c(.35, 4,  22), #16
  c(.2,  4.5, 39),#17
  c(.21, 3,  25), #18
  c(.3,  3.3, 10),#19
  c(.35, 3.3, 10),#20
  c(.46, 5.5, 10),#21
  c(.2,  5.8, 18),#22
  c(.74, 15,  19),#23
  c(.5,  4.55, 8),#24
  c(.7,  7,   13),#25
  c(.6,  3.4, 11) #26
)


# . Fig. S1 ----
### Plot the empirical distribution of records per sampling event
par(mar=c(5.1,2.1,4.1,2.1))
frame()
plot.window(xlim=c(0,52), ylim=c(0,40))
axis(side=1, at=c(0,26+.39,52+.5), labels=NA, tck=-.015, pos=0)
text(x=c(0,26+.39,52+.5), y=par("usr")[3], labels=c('65S','0','65N'), cex=.7, xpd=T, adj=c(.5,1))
axis(side=2, at=seq(0,40,10), labels=NA, las=1, tck=-.01, pos=0)
text(x=par("usr")[1]+1.4, y=seq(0,40,10), labels=seq(0,40,10), cex=.7, xpd=T, adj=c(1,.5))

matriz <- ophiAtlant
matriz[,9] <- cut(matriz[,3], breaks = seq(-65, 65, 5))
lev <- levels(matriz[,9])
pos <- -1

# For each latitude...
for(i in 1:length(lev))
{
  pos <- pos+2
  
  # Create a submatrix with the data of this latitudinal band
  z <- which(matriz[,9]==lev[i])
  tempMat <- matriz[z,]
  tempMat2 <- cbind(tempMat[,1:4], count=rep(1, nrow(tempMat)))
  
  # Transform the submatrix in a community matrix
  spBySamp <- as.vector(table(paste(tempMat2[,2], tempMat2[,3], tempMat2[,4], sep = '_')))
  print(paste(lev[i],': max = ', max(spBySamp), sep = ''))
  
  # I plot only one outlier per position (i.e. when they overlap) because otherwise the figure becomes too heavy
  #boxplot(spBySamp, horizontal = F, add=T, at=pos, frame.plot=F, axes=F, col="deepskyblue", outpch=21, outbg="deepskyblue", outcex=.5, boxlwd = .25, outlwd=.25, lwd=.35, width=.6)
  a <- boxplot(spBySamp, horizontal = F, add=T, at=pos, frame.plot=F, axes=F, col="deepskyblue", outpch=21, outbg="deepskyblue", outcex=.5, boxlwd = .25, outlwd=.25, lwd=.35, width=.6, plot=F)
  p1 <- recordPlot(boxplot(spBySamp, horizontal = F, add=T, at=pos, frame.plot=F, axes=F, col="deepskyblue", boxlwd = .25, lwd=.35, width=.6, outline=F))
  p1 <- recordPlot(points(y=unique(a$out), x=rep(pos, length(unique(a$out))), pch=21, bg="deepskyblue", cex=.5, lwd=.25))
}
# Run again the last row to save the result of the last command
p1 <- recordPlot(points(y=unique(a$out), x=rep(pos, length(unique(a$out))), pch=21, bg="deepskyblue", cex=.5, lwd=.25))


### Plot the simulated distribution of records per sampling event
pos <- -0.2
for(i in 1:length(lev))
{
  pos <- pos+2
  z <- which(matriz[,9]==lev[i])
  tempMat <- matriz[z,]
  
  tempMat2 <- cbind(tempMat[,1:4], count=rep(1, nrow(tempMat)))
  spBySamp <- as.vector(table(paste(tempMat2[,2], tempMat2[,3], tempMat2[,4], sep = '_')))
  
  # We can plot the distribution based on a single simulation
  if(TRUE)
  {
    x <- rbeta(length(spBySamp), RefMat[i,1], RefMat[i,2])*RefMat[i,3]
    x <- ceiling(x)
    
    a <- boxplot(x, horizontal = F, add=T, at=pos, frame.plot=F, axes=F, outpch=21, outcex=.5, boxlwd = .25, outlwd=.25, lwd=.35, width=.6, plot=F)
    p2 <- recordPlot(boxplot(x, horizontal = F, add=T, at=pos, frame.plot=F, axes=F, boxlwd = .25, lwd=.35, width=.6, col='grey', medcol='red', outline=F))
    p2 <- recordPlot(points(y=unique(a$out), x=rep(pos, length(unique(a$out))), pch=21, cex=.5, lwd=.25, bg='grey'))
    #p2 <- recordPlot(points(y=unique(a$out), x=rep(pos, length(unique(a$out))), pch=21, cex=.5, lwd=.25, bg="grey"))
  }
  else
  {
    # Or we can plot the average distribution based on one hundred simulations
    alfa <- 1/100
    for(j in 1:100)
    {
      x <- rbeta(length(spBySamp), RefMat[i,1], RefMat[i,2])*RefMat[i,3]
      x <- ceiling(x)
      
      a <- boxplot(x, horizontal = F, add=T, at=pos, frame.plot=F, axes=F, outpch=21, outbg="grey", outcex=.5, boxlwd = .25, outlwd=.25, lwd=.35, width=.6, plot=F)
      p2 <- recordPlot(boxplot(x, horizontal = F, add=T, at=pos, frame.plot=F, axes=F, boxlwd = .25, lwd=.35, width=.6, outline=F, col=rgb(red = .5, green = .5, blue = .5, alpha = alfa), border=rgb(red = 0, green = 0, blue = 0, alpha = alfa), medcol=rgb(red = 1, green = 0, blue = 0, alpha = alfa)))
      p2 <- recordPlot(points(y=unique(a$out), x=rep(pos, length(unique(a$out))), pch=21, cex=.5, lwd=.25, bg=rgb(red = .5, green = .5, blue = .5, alpha = alfa), col=rgb(red = 0, green = 0, blue = 0, alpha = alfa)))
      #p2 <- recordPlot(points(y=unique(a$out), x=rep(pos, length(unique(a$out))), pch=21, cex=.5, lwd=.25, bg=rgb(red = .5, green = .5, blue = .5, alpha = alfa), col=rgb(red = 0, green = 0, blue = 0, alpha = alfa)))
    }
  }
}
# Run again the last code row to save the result of the last command


### Final plot (B)
p2


### Plot total number of records per latitude

# Empirical records
vec <- numeric()
for(i in 1:length(lev))
{
  z <- which(matriz[,9]==lev[i])
  vec <- c(vec, length(z))
}

# Simulated records (100 simulations)
matX <- matrix(nrow=100, ncol = length(lev))
for(j in 1:length(lev))
{
  z <- which(matriz[,9]==lev[j])
  tempMat <- matriz[z,]
  
  for(i in 1:100)
  {
    # How many specimens I expect to record at each latitude applying the same number of sampling events
    x <- rbeta(sum(duplicated(tempMat[,2:4])==F), RefMat[j,1], RefMat[j,2])*RefMat[j,3]
    x <- ceiling(x)
    matX[i,j] <- sum(x)
  }
  print(j)
}


### final plot (A)
par(lwd=.3)
a <- barplot(rbind(vec,colMeans(matX)), col=c('deepskyblue','grey'), ylim=c(0,50000), space=c(0,1), xaxs="i", xaxt="n", axes=F, beside=T)
axis(side=2, at=seq(0,50000,5000), labels=NA, tck=-.01, hadj=1)
text(x=-.06, y=seq(0,50000,5000), labels=seq(0,50,5), cex=.7, xpd=T, adj=c(1,.5))
axis(side=1, at=c(1,39.5,78), labels=NA, tck=-.015, hadj=1)
text(x=c(1,39.5,78), y=par("usr")[3]-1800, labels=c('65S','0','65N'), cex=.7, xpd=T, adj=c(.5,1))
arrows(a[2,],colMeans(matX)-apply(matX, 2, sd),a[2,],colMeans(matX)+apply(matX, 2, sd),lwd = .3, angle = 90, code = 3, length = 0.015)


### Combine the full reference matrix (Latitude, first cell's ID, number of cells, sampling events, parameters to replicate records / sample [alpha, beta, multipler])
CellLatMat <- cbind(CellLatMat,RefMat)


### Stop here and export data if grid cells at each latitude will be sampled at random
#write.table(CellLatMat, file="./Data/SimInput/CellLatMat.txt",row.names = F, col.names = F)



#### Sampled cells per latitude ####

### Create a spatial point for each sampling event
colnames(ophiAtlant)
ophiPoints <- unique(ophiAtlant[,c(2:4)])
ophiPoints <- vect(ophiPoints, geom=c("decimalLongitude","decimalLatitude"), crs='epsg:4326')

### Get the number of sampling events by grid cell (concentration)
sampByCell <- relate(AtlGrid, ophiPoints, relation='intersects', pairs=F)
sampByCell <- rowSums(sampByCell)

### Get the grid cells with sampling events (distribution)
sampCells <- ifelse(sampByCell > 0, 1, 0)

### Calculate the total distribution and spatial concentration of sampling events by latitude
S_dist <- numeric()
S_dist90 <- numeric()
for(i in 1:nrow(CellLatMat))
{
  x <- (CellLatMat[i,2]+1):(CellLatMat[i,2]+CellLatMat[i,3])
  S_dist[i] <- sum(sampCells[x])/length(x)
  
  xSum <- cumsum(sort(sampByCell[x]/sum(sampByCell[x]), decreasing = T))
  S_dist90[i] <- which(xSum>.9)[1]/length(x)
}


### Find the parameters that best describe the sampling distribution by latitude
pCellLat <- numeric()
for(i in 1:nrow(CellLatMat))
{
  # For each latitude
  cont <- 0
  prm <- 1
  dif <- 1
  
  while(cont < 5)
  {
    # Calculate the percentage of grid cells sampled using the parameter 'prm'
    pCells <- numeric()
    for(j in 1:2000)
    {
      x <- floor((rbeta(n = CellLatMat[i,4], shape1 = prm, shape2 = prm)*CellLatMat[i,3])+CellLatMat[i,2])
      pCells[j] <- length(unique(x))/CellLatMat[i,3]
    }
    
    # If this parameters produces, on average, the closest percentage to the real distribution, save the information
    tDif <- mean(abs(pCells-S_dist[i]))
    if(tDif <= dif)
    {
      dif <- tDif
      pCellLat[i] <- prm
    }
    else
    {
      # Else, inform that the new prm is not increasing the similarity
      cont <- cont+1
    }
    
    prm <- prm+.5
  }
  
  print(i)
}


# . Fig. S8a ----
### Simulate the spatial distribution and concentration of sampling events by latitude
meanSD <- numeric()
meanSD90 <- numeric()
for(i in 1:nrow(CellLatMat))
{
  distTemp <- numeric()
  distTemp90 <- numeric()
  for(j in 1:100)
  {
    x <- floor((rbeta(n = CellLatMat[i,4], shape1 = pCellLat[i], shape2 = pCellLat[i])*CellLatMat[i,3])+CellLatMat[i,2])
    distTemp[j] <- length(table(x))/CellLatMat[i,3]
    
    xSum <- cumsum(sort(table(x)/sum(table(x)), decreasing = T))
    distTemp90[j] <- which(xSum>.9)[1]/CellLatMat[i,3]
  }
  meanSD <- rbind(meanSD, c(mean(distTemp), sd(distTemp)))
  meanSD90 <- rbind(meanSD90, c(mean(distTemp90), sd(distTemp90)))
}
meanSD <- meanSD[nrow(meanSD):1,]
meanSD[meanSD[,2]==0,2] <- NA

meanSD90 <- meanSD90[nrow(meanSD90):1,]
meanSD90[meanSD90[,2]==0,2] <- NA


### Plot
# Real-world distribution and concentration
x <- barplot(rev(S_dist), space=0, axes=F, col=rgb(0,.749,1,.2), ylim=c(0,1.05)) #deepskyblue
barplot(rev(S_dist90), add=T, space=0, axes=F, col=rgb(.008,.318,.588,.9)) ##025196
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65N','0','65S'), cex.axis=1.4)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=3, pos=-.5, cex.axis=1.2)

# Simulated distribution
lines(x, meanSD[,1], lty=2)
arrows(x, meanSD[,1]-meanSD[,2], x, meanSD[,1]+meanSD[,2], lwd=1.5, angle=90, code=3, length=0.02)
points(x=x, y=meanSD[,1], pch=21, bg='lightgrey', col="black", cex=1.5)

# Simulated concentration
lines(x, meanSD90[,1], lty=2)
arrows(x, meanSD90[,1]-meanSD90[,2], x, meanSD90[,1]+meanSD90[,2], lwd=1.5, angle=90, code=3, length=0.02)
points(x=x, y=meanSD90[,1], pch=21, bg='#5d5d5d', col="black", cex=1.5)


### Compare with a random distribution (use the code above to plot)
meanSD <- numeric()
meanSD90 <- numeric()
for(i in 1:nrow(CellLatMat))
{
  distTemp <- numeric()
  distTemp90 <- numeric()
  for(j in 1:100)
  {
    x <- floor((rbeta(n = CellLatMat[i,4], shape1 = 1, shape2 = 1)*CellLatMat[i,3])+CellLatMat[i,2])
    distTemp[j] <- length(table(x))/CellLatMat[i,3]
    
    xSum <- cumsum(sort(table(x)/sum(table(x)), decreasing = T))
    distTemp90[j] <- which(xSum>.9)[1]/CellLatMat[i,3]
  }
  meanSD <- rbind(meanSD, c(mean(distTemp), sd(distTemp)))
  meanSD90 <- rbind(meanSD90, c(mean(distTemp90), sd(distTemp90)))
}
meanSD <- meanSD[nrow(meanSD):1,]
meanSD[meanSD[,2]==0,2] <- NA

meanSD90 <- meanSD90[nrow(meanSD90):1,]
meanSD90[meanSD90[,2]==0,2] <- NA


### Combine the full reference matrix (Latitude, first cell's ID, number of cells, sampling events, parameters to replicate records / sample [alpha, beta, multipler], parameters to replicate sampled cells / latitude)
CellLatMat <- cbind(CellLatMat,pCellLat)


### Export CellLatMat matrix
#write.table(CellLatMat, file="./Data/SimInput/CellLatMat.txt",row.names = F, col.names = F)



#### END ####