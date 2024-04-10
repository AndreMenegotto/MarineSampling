
#### Analyses ####

### Load packages
library(vioplot)
library(viridis)
library(raster)
library(iNEXT)
library(rgeos)


### Load basic RData
load("./Data/Empirical/Atlantic.RData")
load("./Data/Empirical/Functions.RData")
CellLatMat <- read.table("./Data/SimInput/CellLatMat.txt", h=F, na.strings = -1)

# Original data from https://www.naturalearthdata.com/
land <- shapefile("./Data/Empirical/ne_110m_land/ne_110m_land.shp")
landCut <- suppressWarnings(gIntersection(land, as(extent(AtlGrid), "SpatialPolygons"), byid = TRUE))

# Lat values
OphiLat <- read.table("./Data/Empirical/OphiLat.txt", h=T)
SPrich <- OphiLat$SPrich
SPgaps <- OphiLat$SPgaps
USamp <- OphiLat$USamp
OphiSC <- OphiLat[,4:6]

# Grid cell values
OphiCell <- read.table("./Data/Empirical/OphiCell.txt", h=T)

# Average grid cell values by latitude
levs <- rev(levels(cut(0, breaks = seq(-65, 65, 5))))

USampCell <- numeric()
SPrichCell <- numeric()
OphiSCCell <- numeric()
for(i in 1:length(levs))
{
  pos <- which(OphiCell$Lat==rev(levs)[i])
  
  USampCell <- rbind(USampCell, c(mean(OphiCell$USamp[pos], na.rm = T), sd(OphiCell$USamp[pos], na.rm=T)))
  SPrichCell <- rbind(SPrichCell, c(mean(OphiCell$SPrich[pos], na.rm = T), sd(OphiCell$SPrich[pos], na.rm=T)))
  OphiSCCell <- rbind(OphiSCCell, c(mean(OphiCell$SC[pos], na.rm = T), sd(OphiCell$SC[pos], na.rm=T)))
}


### Function to average grid cell values by latitude
grid2lat <- function(estMat, CellLatMat=CellLatMat)
{
  meanSD <- numeric()
  for(i in 1:nrow(CellLatMat))
  {
    # Find the grid cells within this specific latitude
    rowPos <- (CellLatMat[i,2]+1):(CellLatMat[i,2]+CellLatMat[i,3])
    
    # Latitudinal average of the grid cells averaged across simulations
    avVec <- mean(colMeans(estMat[rowPos,], na.rm = T), na.rm = T)
    
    # Latitudinal SD averaged across simulations: margin argument = 2
    n <- apply(estMat[rowPos,], 2, function(x) sum(!is.na(x)))
    n <- n[n>1]
    sdVec <- na.omit(apply(estMat[rowPos,], 2, sd, na.rm=T))
    sdVec <- sqrt(sum((n-1)*sdVec^2)/(sum(n)-length(sdVec)))
    
    # Combine values
    meanSD <- rbind(meanSD, c(avVec, sdVec))
  }
  return(meanSD)
}

### Notice
# .Real-world data sets are structured S -> N
# ..Simulated data sets are structured N -> S


# . Figure 1 ----

## Plot the sampling map
plot(AtlGrid, col="grey95", border="transparent")

# This requires the raw data set of occurrence records
#load("./Data/Empirical/Raw/Ophiuroids.RData")
#SampSite <- unique(ophiAtlant[,2:3])
#coordinates(SampSite) <- ~decimalLongitude+decimalLatitude
#points(SampSite, pch=20, col="#187bcd", cex=.2)

plot(landCut, col="grey80", border="transparent", add=T)
plot(as(extent(AtlGrid), "SpatialPolygons"), add=T)
segments(x0 = extent(AtlGrid)[1], x1 = extent(AtlGrid)[2], y0 = 0, y1 = 0, lty = 2)


## Plot the richness maps (choose the variable to plot)
x <- OphiCell$SPrich

MatRichCelly <- read.table("./Data/SimOutput/RSFD39_gcY/MatRichCell.txt", h=F, na.strings = -1)
x <- rowMeans(MatRichCelly, na.rm = T)
MatObsCelly <- read.table("./Data/SimOutput/RSFD39_gcY/MatObsCell.txt", h=F, na.strings = c(-1,0))
x <- rowMeans(MatObsCelly, na.rm = T)

# For the area model I used only one simulation instead of the simulations' average to better visualize the resulting random pattern
MatRichCelln <- read.table("./Data/SimOutput/RSFD39_gcN/MatRichCell.txt", h=F, na.strings = -1)
x <- MatRichCelln[,3]
MatObsCelln <- read.table("./Data/SimOutput/RSFD39_gcN/MatObsCell.txt", h=F, na.strings = c(-1,0))
x <- MatObsCelln[,3]

# Note: zero richness values are treated as NA in the Observed Matrix because grid cells lacking values are typically omitted from analyses with real-world data sets.
# To include zero values in the analyses, we must ensure that they signify the absence of species, not the absence of sampling.


# Create the color scale
escala <- cut(x, seq(floor(min(x, na.rm=T)/5)*5, ceiling(max(x, na.rm=T)/5)*5, 5))

colCode <- numeric(length(AtlGrid@data[,1]))
for(i in 1:length(levels(escala)))
{
  pos <- which(escala==levels(escala)[i])
  colCode[pos] <- viridis(length(levels(escala)))[i]
}

# Plot the map
plot(AtlGrid, col=colCode, border="transparent")
plot(landCut, col="grey80", border="transparent", add=T)
plot(as(extent(AtlGrid), "SpatialPolygons"), add=T)
segments(x0 = extent(AtlGrid)[1], x1 = extent(AtlGrid)[2], y0 = 0, y1 = 0, lty = 2)


## Gamma diversity
x <- barplot(rev(SPrich/max(SPrich)), col=rgb(.6,.6,.6), space=0, ylim=c(0,1), axes=F, axisnames=F)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65N','0','65S'), cex.axis=1.4)
axis(side=4, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=0, pos=26+1, lwd=1.4, cex.axis=1.4)

MatRichLaty <- read.table("./Data/SimOutput/RSFD39_gcY/MatRichLat.txt", h=F, na.strings = -1)
MatObsLaty <- read.table("./Data/SimOutput/RSFD39_gcY/MatObsLat.txt", h=F, na.strings = c(-1,0))

MatRichLatn <- read.table("./Data/SimOutput/RSFD39_gcN/MatRichLat.txt", h=F, na.strings = -1)
MatObsLatn <- read.table("./Data/SimOutput/RSFD39_gcN/MatObsLat.txt", h=F, na.strings = c(-1,0))

# These plots are reversed because they are used in the vertical, next to the maps
x <- barplot(rowMeans(MatRichLaty)/max(rowMeans(MatRichLaty)), col=rgb(.95,.95,.95), space=0, ylim=c(0,1), axes=F)
barplot(rowMeans(MatObsLaty)/max(rowMeans(MatRichLaty)), col=rgb(.6,.6,.6), space=0, axes=F, add=T)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65N','0','65S'), cex.axis=1.4)
axis(side=4, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=0, pos=26+1, lwd=1.4, cex.axis=1.4)

barplot(rowMeans(MatRichLatn)/max(rowMeans(MatRichLatn)), col=rgb(.95,.95,.95), space=0, ylim=c(0,1), axes=F)
barplot(rowMeans(MatObsLatn)/max(rowMeans(MatRichLatn)), col=rgb(.6,.6,.6), space=0, axes=F, add=T)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65N','0','65S'), cex.axis=1.4)
axis(side=4, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=0, pos=26+1, lwd=1.4, cex.axis=1.4)


# Correlations
plot(rowMeans(MatRichLaty), rowMeans(MatRichLatn))
round(cor(rowMeans(MatRichLaty), rowMeans(MatRichLatn), method = 'p'), 2) #-0.54
plot(rowMeans(MatObsLaty), rowMeans(MatObsLatn))
round(cor(rowMeans(MatObsLaty), rowMeans(MatObsLatn), method = 'p'), 2) #0.93

plot(rowMeans(MatObsLaty), rowMeans(MatRichLaty))
round(cor(rowMeans(MatObsLaty), rowMeans(MatRichLaty), method = 'p'), 2) #-0.03
plot(rowMeans(MatObsLatn), rowMeans(MatRichLatn))
round(cor(rowMeans(MatObsLatn), rowMeans(MatRichLatn), method = 'p'), 2) #0.43

plot(rev(rowMeans(MatObsLaty)), SPrich)
round(cor(rev(rowMeans(MatObsLaty)), SPrich, method = 'p'), 2) #0.84
plot(rev(rowMeans(MatObsLatn)), SPrich)
round(cor(rev(rowMeans(MatObsLatn)), SPrich, method = 'p'), 2) #0.68


## Alpha diversity
x <- barplot(rev(SPrichCell[,1]/max(SPrichCell[,1])), col=rgb(.6,.6,.6), space=0, ylim=c(0,1), axes=F)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65N','0','65S'), cex.axis=1.4)
axis(side=4, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=0, pos=26+1, lwd=1.4, cex.axis=1.4)

avRichy <- grid2lat(estMat = MatRichCelly, CellLatMat = CellLatMat)[,1]
avObsy <- grid2lat(estMat = MatObsCelly, CellLatMat = CellLatMat)[,1]

avRichn <- grid2lat(estMat = MatRichCelln, CellLatMat = CellLatMat)[,1]
avObsn <- grid2lat(estMat = MatObsCelln, CellLatMat = CellLatMat)[,1]

x <- barplot(avRichy/max(avRichy), col=rgb(.95,.95,.95), space=0, ylim=c(0,1), axes=F)
barplot(avObsy/max(avRichy), col=rgb(.6,.6,.6), space=0, axes=F, add=T)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65N','0','65S'), cex.axis=1.4)
axis(side=4, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=0, pos=26+1, lwd=1.4, cex.axis=1.4)

x <- barplot(avRichn/max(avRichn), col=rgb(.95,.95,.95), space=0, ylim=c(0,1), axes=F)
barplot(avObsn/max(avRichn), col=rgb(.6,.6,.6), space=0, axes=F, add=T)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65N','0','65S'), cex.axis=1.4)
axis(side=4, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=0, pos=26+1, lwd=1.4, cex.axis=1.4)


# Correlations
plot(avRichy, avRichn)
round(cor(avRichy, avRichn, method = 'p'), 2) #-0.21
plot(avObsy, avObsn)
round(cor(avObsy, avObsn, method = 'p'), 2) #0.98

plot(avObsy, avRichy)
round(cor(avObsy, avRichy, method = 'p'), 2) #-0.27
plot(avObsn, avRichn)
round(cor(avObsn, avRichn, method = 'p'), 2) #-0.09

plot(rev(avObsy), SPrichCell[,1])
round(cor(rev(avObsy), SPrichCell[,1], method = 'p'), 2) #0.54
plot(rev(avObsn), SPrichCell[,1])
round(cor(rev(avObsn), SPrichCell[,1], method = 'p'), 2) #0.48




# . Figure 2 ----

## Sampling effort
SampData <- numeric()
for(i in 1:length(levs))
{
  x <- which(OphiCell$Lat==rev(levs)[i])
  SampData <- rbind(SampData, cbind(OphiCell$USamp[x], i))
}
SampData <- as.data.frame(SampData)
colnames(SampData) <- c('samp','lat')

frame()
plot.window(xlim=c(1,26), ylim=c(0,5))
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=1.2, lwd=2, padj = .5)
axis(side=2, at=seq(0,5,1), labels=parse(text = paste(10, "^", seq(0,5,1), sep="")), las=1, cex.axis=1.2, lwd=2)

vioplot(log10(SampData$samp)~SampData$lat, frame.plot=F, col='grey90', colMed = rgb(253,179,56, maxColorValue = 255), cex=1.3, add=T)
points(aggregate(log10(samp)~lat, data=SampData, median)[,2], cex=1.5)
lines(log10(USamp), x = 1:26, lty=2)
points(log10(USamp), pch=21, bg=rgb(2,75,122, maxColorValue = 255), col="black", cex=2.5, lwd=1.5)


## Empirical gaps and sample coverage
x <- barplot(SPgaps/max(SPrich), col=rgb(.8,.8,.8), space=0, ylim=c(0,1.025), axes=F, axisname=F)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.2)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.2)

polygon(c(.5:25.5, rev(.5:25.5)), c(OphiSC[,3], rev(OphiSC[,2])), col=rgb(2,81,150,50, maxColorValue = 255), border=NA)
polygon(c(.5:9.5, rev(.5:9.5)), c(OphiSCCell[1:10,1]+OphiSCCell[1:10,2], rev(OphiSCCell[1:10,1]-OphiSCCell[1:10,2])), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
polygon(c((11.5:25.5), rev(11.5:25.5)), c(OphiSCCell[12:26,1]+OphiSCCell[12:26,2], rev(OphiSCCell[12:26,1]-OphiSCCell[12:26,2])), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
lines(OphiSC[,1], x = x, lty=2)
lines(OphiSCCell[,1], x = x, lty=2)
points(OphiSC[,1], x = x, pch=21, bg=rgb(2,81,150, maxColorValue = 255), col="black", cex=2)
points(OphiSCCell[,1], x = x, pch=21, bg=rgb(253,179,56, maxColorValue = 255), col="black", cex=2)


## Geometric constraints model
MatObsLaty <- read.table("./Data/SimOutput/RSFD39_gcY/MatObsLat.txt", h=F, na.strings = c(-1,0))
MatGapsLaty <- read.table("./Data/SimOutput/RSFD39_gcY/MatGapsLat.txt", h=F, na.strings = -1)
SCLaty <- read.table("./Data/SimOutput/RSFD39_gcY/MatCoverageLat.txt", h=F, na.strings = -1)
SCCelly <- read.table("./Data/SimOutput/RSFD39_gcY/MatCoverageCell.txt", h=F, na.strings = -1)

MatGapsLatST <- rowMeans(MatGapsLaty)/max(rowMeans(MatObsLaty))
meanSD1 <- cbind(rowMeans(SCLaty, na.rm=T), apply(SCLaty, 1, sd, na.rm=T))
meanSD2 <- grid2lat(estMat = SCCelly, CellLatMat = CellLatMat)

x <- barplot(rev(MatGapsLatST), col=rgb(.8,.8,.8), space=0, ylim=c(0,1.025), axes=F)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.2)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.2)
polygon(c(.5:25.5, rev(.5:25.5)), c(rev(meanSD1[,1]+meanSD1[,2]), meanSD1[,1]-meanSD1[,2]), col=rgb(2,81,150,50, maxColorValue = 255), border=NA)
polygon(c(.5:25.5, rev(.5:25.5)), c(rev(meanSD2[,1]+meanSD2[,2]), meanSD2[,1]-meanSD2[,2]), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
lines(rev(meanSD1[,1]), x = x, lty=2)
lines(rev(meanSD2[,1]), x = x, lty=2)
points(rev(meanSD1[,1]), x = x, pch=21, bg=rgb(2,81,150, maxColorValue = 255), col="black", cex=2)
points(rev(meanSD2[,1]), x = x, pch=21, bg=rgb(253,179,56, maxColorValue = 255), col="black", cex=2)


## Area model
MatObsLatn <- read.table("./Data/SimOutput/RSFD39_gcN/MatObsLat.txt", h=F, na.strings = c(-1,0))
MatGapsLatn <- read.table("./Data/SimOutput/RSFD39_gcN/MatGapsLat.txt", h=F, na.strings = -1)
SCLatn <- read.table("./Data/SimOutput/RSFD39_gcN/MatCoverageLat.txt", h=F, na.strings = -1)
SCCelln <- read.table("./Data/SimOutput/RSFD39_gcN/MatCoverageCell.txt", h=F, na.strings = -1)

MatGapsLatST <- rowMeans(MatGapsLatn)/max(rowMeans(MatObsLatn))
meanSD1 <- cbind(rowMeans(SCLatn, na.rm=T), apply(SCLatn, 1, sd, na.rm=T))
meanSD2 <- grid2lat(estMat = SCCelln, CellLatMat = CellLatMat)

x <- barplot(rev(MatGapsLatST), col=rgb(.8,.8,.8), space=0, ylim=c(0,1.025), axes=F)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.2)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.2)
polygon(c(.5:25.5, rev(.5:25.5)), c(rev(meanSD1[,1]+meanSD1[,2]), meanSD1[,1]-meanSD1[,2]), col=rgb(2,81,150,50, maxColorValue = 255), border=NA)
polygon(c(.5:12.5, rev(.5:12.5)), c(rev(meanSD2[,1]+meanSD2[,2])[1:13], (meanSD2[,1]-meanSD2[,2])[14:26]), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
polygon(c(14.5:25.5, rev(14.5:25.5)), c(rev(meanSD2[,1]+meanSD2[,2])[15:26], (meanSD2[,1]-meanSD2[,2])[1:12]), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
lines(rev(meanSD1[,1]), x = x, lty=2)
lines(rev(meanSD2[,1]), x = x, lty=2)
points(rev(meanSD1[,1]), x = x, pch=21, bg=rgb(2,81,150, maxColorValue = 255), col="black", cex=2)
points(rev(meanSD2[,1]), x = x, pch=21, bg=rgb(253,179,56, maxColorValue = 255), col="black", cex=2)


# Correlations (completeness between models)
round(cor(rowMeans(SCLaty, na.rm=T), rowMeans(SCLatn, na.rm=T), method = 'p'), 2) #0.99

meanSDy <- grid2lat(estMat = SCCelly, CellLatMat = CellLatMat)[,1]
meanSDn <- grid2lat(estMat = SCCelln, CellLatMat = CellLatMat)[,1]
round(cor(meanSDy, meanSDn, method = 'p'), 2) #0.99


# Correlations (completeness vs. sampling)
plot(USamp, rev(rowMeans(SCLaty, na.rm=T)))
round(cor(rev(rowMeans(SCLaty, na.rm=T)), USamp, method = 'spearman'), 2)  #0.98
round(cor(rev(rowMeans(SCLatn, na.rm=T)), USamp, method = 'spearman'), 2)  #0.99

plot(aggregate(SampData$samp, by=list(SampData$lat), mean)[,2], rev(meanSDy))
round(cor(rev(meanSDy), aggregate(SampData$samp, by=list(SampData$lat), mean)[,2], method = 'spearman'), 2) #0.85
round(cor(rev(meanSDn), aggregate(SampData$samp, by=list(SampData$lat), mean)[,2], method = 'spearman'), 2)  #0.85


# Correlations (completeness|sampling vs. gaps)
plot(rowMeans(MatGapsLaty), rowMeans(SCLaty, na.rm=T))
round(cor(rowMeans(SCLaty, na.rm=T), rowMeans(MatGapsLaty), method = 'p'), 2)  #-0.99
round(cor(rowMeans(SCLatn, na.rm=T), rowMeans(MatGapsLatn), method = 'p'), 2)  #-0.99

plot(USamp, rev(rowMeans(MatGapsLaty)))
round(cor(rev(rowMeans(MatGapsLaty)), USamp, method = 'spearman'), 2) #-0.92
round(cor(rev(rowMeans(MatGapsLatn)), USamp, method = 'spearman'), 2) #-0.93


# Correlations (completeness: simulated vs. real world)
plot(rev(rowMeans(SCLaty, na.rm=T)), OphiSC[,1])
round(cor(OphiSC[,1], rev(rowMeans(SCLaty, na.rm=T)), method = 'spearman'), 2) #0.96
round(cor(OphiSC[,1], rev(rowMeans(SCLatn, na.rm=T)), method = 'spearman'), 2) #0.94

rmNA <- which(is.na(OphiSCCell[,1]))
plot(OphiSCCell[-rmNA,1], rev(meanSDy[-rmNA]))
round(cor(rev(meanSDy[-rmNA]), OphiSCCell[-rmNA,1], method = 'spearman'), 2) #0.35
round(cor(rev(meanSDn[-rmNA]), OphiSCCell[-rmNA,1], method = 'spearman'), 2) #0.33

plot(rev(rowMeans(MatGapsLaty)), SPgaps)
round(cor(SPgaps,rev(rowMeans(MatGapsLaty)), method = 'p'), 2) #0.36
round(cor(SPgaps,rev(rowMeans(MatGapsLatn)), method = 'p'), 2) #0.40


# Correlations (sampling vs. completeness vs. gaps - real world)
plot(USamp, OphiSC[,1])
round(cor(OphiSC[,1], USamp, method = 'spearman'), 2) #0.94

plot(SPgaps, OphiSC[,1])
round(cor(OphiSC[,1], SPgaps, method = 'spearman'), 2) #-0.45

plot(SPgaps, USamp)
round(cor(USamp, SPgaps, method = 'spearman'), 2) #-0.28




# . Figure 3 ----
myRd <- colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee090'))
myBu <- colorRampPalette(c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
color <- c(myRd(7)[-7],'#ffdb58', myBu(7)[-1]) #e5e5ab #ffdb58
color <- c(rev(color), color)

## Latitude (choose the variable to plot)
MatRichLat <- read.table("./Data/SimOutput/RSFD39_gcY/MatRichLat.txt", h=F, na.strings = -1)
MatChao2Lat <- read.table("./Data/SimOutput/RSFD39_gcY/MatChao2Lat.txt", h=F, na.strings = -1)
MatInExtLat <- read.table("./Data/SimOutput/RSFD39_gcY/MatInExtLat.txt", h=F, na.strings = -1)
MatES50Lat <- read.table("./Data/SimOutput/RSFD39_gcY/MatES50Lat.txt", h=F, na.strings = -1)

estMat <- MatRichLat
estMat <- MatChao2Lat
estMat <- MatInExtLat
estMat <- MatES50Lat

estSTD <- rowMeans(estMat, na.rm=T)
estSTD <- apply(cbind(estSTD/max(estSTD, na.rm=T), apply(estMat, 1, sd, na.rm=T)/max(estSTD, na.rm=T)), 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

arrows(1:26, estSTD[,1], 1:26, estSTD[,1]-estSTD[,2],lwd = 1.5, angle = 90,code = 3, length = 0.025)
lines(estSTD[,1], lwd=2, lty=2)
points(estSTD[,1], pch=21, bg=rgb(.45,.45,.45), col="black", cex=2.5, lwd=1.5)


## Grid cell (average)
MatRichCell <- read.table("./Data/SimOutput/RSFD39_gcY/MatRichCell.txt", h=F, na.strings = -1)
MatChao2Cell <- read.table("./Data/SimOutput/RSFD39_gcY/MatChao2Cell.txt", h=F, na.strings = -1)
MatInExtCell <- read.table("./Data/SimOutput/RSFD39_gcY/MatInExtCell.txt", h=F, na.strings = -1)
MatES50Cell <- read.table("./Data/SimOutput/RSFD39_gcY/MatES50Cell.txt", h=F, na.strings = -1)

estMat <- MatRichCell
estMat <- MatChao2Cell
estMat <- MatInExtCell
estMat <- MatES50Cell

estSTD <- grid2lat(estMat = estMat, CellLatMat = CellLatMat)
estSTD <- apply(estSTD/max(estSTD[,1], na.rm = T), 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

arrows(1:26, estSTD[,1], 1:26, estSTD[,1]-estSTD[,2], lwd=2.5, angle=90, code=3, length=0.03)
lines(estSTD[,1], lwd=2, lty=2)
points(estSTD[,1], pch=21, bg=rgb(.45,.45,.45), col="black", cex=2.5, lwd=1.5)


## Grid cell (same variable used before)
estSTD <- rowMeans(estMat, na.rm=T)
estSTD <- estSTD/max(estSTD, na.rm=T)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

for(i in nrow(CellLatMat):1)
{
  rowPos <- (CellLatMat[i,2]+1):(CellLatMat[i,2]+CellLatMat[i,3])
  points(x=rep(27-i, length(rowPos)), y=estSTD[rowPos], pch=21, bg=color[i], col="black", cex=2.5, lwd=.5)
}




# . Figure 4 ----
# Open files
MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat/MatRichCell.txt", h=F, na.strings = -1)
MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop99/MatRichCell.txt", h=F, na.strings = -1)
MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatRichCell.txt", h=F, na.strings = -1)
x <- rowMeans(MatRichCelly, na.rm = T)

# Create color scale
escala <- cut(x, seq(floor(min(x, na.rm=T)), ceiling(max(x, na.rm=T)), 1))

colCode <- numeric(length(AtlGrid@data[,1]))
for(i in 1:length(levels(escala)))
{
  pos <- which(escala==levels(escala)[i])
  colCode[pos] <- viridis(length(levels(escala)))[i]
}

# Plot the map
plot(AtlGrid, col=colCode, border="transparent")
plot(landCut, col="grey80", border="transparent", add=T)
plot(as(extent(AtlGrid), "SpatialPolygons"), add=T)
segments(x0 = extent(AtlGrid)[1], x1 = extent(AtlGrid)[2], y0 = 0, y1 = 0, lty = 2)


# Histograms
MatRichLaty <- read.table("./Data/SimOutput/RSFD_nat/MatRichLat.txt", h=F, na.strings = -1)
MatObsLaty <- read.table("./Data/SimOutput/RSFD_nat/MatObsLat.txt", h=F, na.strings = c(-1,0))
MatGapsLaty <- read.table("./Data/SimOutput/RSFD_nat/MatGapsLat.txt", h=F, na.strings = -1)

MatRichLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop99/MatRichLat.txt", h=F, na.strings = -1)
MatObsLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop99/MatObsLat.txt", h=F, na.strings = c(-1,0))
MatGapsLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop99/MatGapsLat.txt", h=F, na.strings = -1)

MatRichLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatRichLat.txt", h=F, na.strings = -1)
MatObsLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatObsLat.txt", h=F, na.strings = c(-1,0))
MatGapsLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatGapsLat.txt", h=F, na.strings = -1)


MatRichLatyST <- rowMeans(MatRichLaty)/max(rowMeans(MatObsLaty))
MatObsLatyST <- rowMeans(MatObsLaty)/max(rowMeans(MatObsLaty))
ExpectedST <- rowMeans(MatObsLaty+MatGapsLaty)/max(rowMeans(MatObsLaty))
maxV <- max(MatRichLatyST)

x <- barplot(MatRichLatyST/maxV, col=rgb(.95,.95,.95), space=0, ylim=c(0,1.05), axes=F)
barplot(ExpectedST/maxV, col=rgb(.8,.8,.8), space=0, axes=F, add=T)
barplot(MatObsLatyST/maxV, col=rgb(.6,.6,.6), space=0, axes=F, add=T)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65N','0','65S'))
axis(side=4, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=0, pos=26+1, lwd=1.5, cex.axis=2, padj=.2)

lines(x = (1:26)-.5, (rev(SPrich/max(SPrich)))/maxV, lty=2)
points(x = (1:26)-.5, (rev(SPrich/max(SPrich)))/maxV, pch=21, bg=rgb(.35,.35,.35, maxColorValue = 1), col="black", cex=1.5)




# . Figure 5 ----
myRd <- colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee090'))
myBu <- colorRampPalette(c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
color <- c(myRd(7)[-7],'#ffdb58', myBu(7)[-1]) #e5e5ab #ffdb58
color <- c(rev(color), color)

int_f <- function(x, mu1, mu2, sd1, sd2)
{
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}

## Latitude (choose the variable to plot)
MatRichLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatRichLat.txt", h=F, na.strings = -1)
MatChao2Lat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatChao2Lat.txt", h=F, na.strings = -1)
MatInExtLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatInExtLat.txt", h=F, na.strings = -1)
MatES50Lat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatES50Lat.txt", h=F, na.strings = -1)

estMat <- MatChao2Lat
estMat <- MatInExtLat
estMat <- MatES50Lat

# Plot estimate from first sampling scenario
estSTD <- estMat[,1]
estSTD <- apply(cbind(estSTD/max(estSTD, na.rm=T), estMat[,2]/max(estSTD, na.rm=T)), 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.5, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=1.5, lwd=2, padj = .5)

refMeanSD <- cbind(rowMeans(MatRichLat), apply(MatRichLat, 1, sd))
refMeanSD <- refMeanSD/max(refMeanSD[,1])
polygon(c(.5:25.5, rev(.5:25.5)), c(rev(refMeanSD[,1]+refMeanSD[,2]), refMeanSD[,1]-refMeanSD[,2]), col=rgb(.35,.35,.35,.2, maxColorValue = 1), border=NA)

arrows(1:26, estSTD[,1], 1:26, estSTD[,1]-estSTD[,2],lwd = 1.5, angle = 90,code = 3, length = 0.025)
lines(estSTD[,1], lwd=2, lty=2)
points(estSTD[,1], pch=21, bg=rgb(.45,.45,.45), col="black", cex=2.5, lwd=1.5)


# Standardize variable and optimize the plot order based on changes between sampling scenarios
refMean <- rowMeans(MatRichLat)
refSD <- apply(MatRichLat, 1, sd)

mu1 <- refMean/max(refMean)
sd1 <- refSD/max(refMean)

matMean <- estMat[,seq(1, ncol(estMat), 3)]
matSD <- estMat[,seq(2, ncol(estMat), 3)]

xn <- matMean[,ncol(matMean)]/apply(matMean, 2, max)[ncol(matMean)]
x1 <- matMean[,1]/apply(matMean, 2, max)[1]
n <- order(xn-x1)

# Plot richness estimation at each latitude
frame()
plot.window(xlim=log10(c(25,6400)), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.5, lwd=1.5)
axis(side=1, at=log10(50*(2^(-1:7))), labels=c(0,50*(2^(0:7)))/100, cex.axis=1.5, lwd=1.5)
segments(x0 = log10(50*(2^(-1:7))), x1 = log10(50*(2^(-1:7))), y0 = 0, y1 = 1, lty=2, lwd=1.5, col='grey')
segments(x0 = log10(50*(2^-1)), x1 = log10(50*(2^7)), y0 = seq(0,1,.2), y1 = seq(0,1,.2), lty=3, lwd=1.5, col='grey90')

for(j in 1:nrow(CellLatMat))
{
  overlVec <- numeric()
  for(i in 1:ncol(matMean))
  {
    mu2 <- matMean[,i]/max(matMean[,i], na.rm=T)
    sd2 <- matSD[,i]/max(matMean[,i], na.rm=T)
    
    if(!is.na(sd2[n[j]]))
    {
      overlVec[i] <- integrate(int_f, -Inf, Inf, mu1=mu1[n[j]], mu2=mu2[n[j]], sd1=sd1[n[j]], sd2=sd2[n[j]])[[1]]
    }
  }
  lines(x=log10(50*(2^(-1:7))), y=overlVec, lwd=3, col=color[n[j]])
}


## Grid cell (average)
MatRichCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatRichCell.txt", h=F, na.strings = -1)
MatChao2Cell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatChao2Cell2Lat.txt", h=F, na.strings = -1)
MatInExtCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatInExtCell2Lat.txt", h=F, na.strings = -1)
MatES50Cell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatES50Cell2Lat.txt", h=F, na.strings = -1)

# Plot estimate from first sampling scenario
estMat <- MatChao2Cell
estMat <- MatInExtCell
estMat <- MatES50Cell

estSTD <- estMat[,1:2]
estSTD <- apply(estSTD/max(estSTD[,1], na.rm = T), 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.5, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=1.5, lwd=2, padj = .5)

refData <- grid2lat(estMat = MatRichCell, CellLatMat = CellLatMat)
refData <- refData/max(refData[,1])
polygon(c(.5:25.5, rev(.5:25.5)), c(rev(refData[,1]+refData[,2]), refData[,1]-refData[,2]), col=rgb(.35,.35,.35,.2, maxColorValue = 1), border=NA)

arrows(1:26, estSTD[,1], 1:26, estSTD[,1]-estSTD[,2], lwd=2.5, angle=90, code=3, length=0.03)
lines(estSTD[,1], lwd=2, lty=2)
points(estSTD[,1], pch=21, bg=rgb(.45,.45,.45), col="black", cex=2.5, lwd=1.5)


# Standardize variable and optimize the plot order based on changes between sampling scenarios
refMeanSD <- grid2lat(estMat = MatRichCell, CellLatMat = CellLatMat)

mu1 <- refMeanSD[,1]/max(refMeanSD[,1])
sd1 <- refMeanSD[,2]/max(refMeanSD[,1])

meanLat <- estMat[,seq(1, ncol(estMat), 3)]
sdLat <- estMat[,seq(2, ncol(estMat), 3)]

xn <- meanLat[,ncol(meanLat)]/apply(meanLat, 2, max)[ncol(meanLat)]
x1 <- meanLat[,1]/apply(meanLat, 2, max, na.rm=T)[1]
n <- order(xn-x1)

# Plot richness estimation at each latitude
frame()
plot.window(xlim=log10(c(25,6400)), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.5, lwd=1.5)
axis(side=1, at=log10(50*(2^(-1:7))), labels=c(0,50*(2^(0:7)))/100, cex.axis=1.5, lwd=1.5)
segments(x0 = log10(50*(2^(-1:7))), x1 = log10(50*(2^(-1:7))), y0 = 0, y1 = 1, lty=2, lwd=1.5, col='grey')
segments(x0 = log10(50*(2^-1)), x1 = log10(50*(2^7)), y0 = seq(0,1,.2), y1 = seq(0,1,.2), lty=3, lwd=1.5, col='grey90')

for(lat in 1:nrow(CellLatMat))
{
  overlVec <- numeric()
  for(i in 1:ncol(meanLat))
  {
    mu2 <- meanLat[,i]/max(meanLat[,i], na.rm=T)
    sd2 <- sdLat[,i]/max(meanLat[,i], na.rm=T)
    
    if(!is.na(sd2[n[lat]]))
    {
      overlVec[i] <- integrate(int_f, -Inf, Inf, mu1=mu1[n[lat]], mu2=mu2[n[lat]], sd1=sd1[n[lat]], sd2=sd2[n[lat]])[[1]]
    }
    else
    {
      overlVec[i] <- NA
    }
  }
  lines(x=log10(50*(2^(-1:7))), y=overlVec, lwd=3, col=color[n[lat]])
}




#### Supplementar ####


# . Figure S1 ----
# On script '1 - InputMatrices.R'




# . Figure S2 ----

### subGrid method

## Plot lat
chaoInput <- function(matriz)
{
  # Data is already round to the nearest 0.01º
  unic <- unique(matriz)
  lastC <- ncol(unic)+1
  
  # Create latitudinal intervals
  unic[,lastC] <- cut(unic[,3], breaks = seq(-65, 65, 5))
  lev <- levels(unic[,lastC])
  
  biota <- list()
  for(i in 1:length(lev))
  {
    # Create a submatrix for this latitudinal band
    Recs <- which(unic[,lastC]==lev[i])
    
    tempMat <- unic[Recs,]
    SampUn <- nrow(unique(tempMat[,c(2:(lastC-1))]))
    regist <- table(as.character(tempMat[,1]))
    
    # List the output
    x <- c(SampUn, as.numeric(regist))
    biota[[i]] <- x
  }  
  return(biota)
}

# This requires the raw data set of occurrence records
#load("./Data/Empirical/Raw/Ophiuroids.RData")
OphiINEXT <- chaoInput(matriz = ophiAtlant[,1:3])

OphiSC_sg <- numeric()
for(i in 1:length(OphiINEXT))
{ 
  # Only for those latitudes with at least 6 samples and number of incidences > number of unique species
  if(OphiINEXT[[i]][1]>5 & any(OphiINEXT[[i]][-1] > 1))
  {
    out <- iNEXT(OphiINEXT[[i]], q=0, datatype="incidence_freq", size = OphiINEXT[[i]][1], se = T, nboot = 50)
    OphiSC_sg <- rbind(OphiSC_sg, out$iNextEst$size_based[1,8:10])
  }
  else
  {
    out <- rep(NA, 3)
    names(out) <- c('SC','SC.LCL','SC.UCL')
    OphiSC_sg <- rbind(OphiSC_sg, out)
  }
  print(i)
}


frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=1.4, lwd=2, padj = .5)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.4, lwd=2)

polygon(c(.5:25.5, rev(.5:25.5)), c(OphiSC[,3], rev(OphiSC[,2])), col=rgb(205,221,234, maxColorValue = 255), border=NA)
lines(OphiSC[,1], x = .5:25.5, lty=2)
points(OphiSC[,1], x = .5:25.5, pch=21, bg=rgb(2,81,150, maxColorValue = 255), col="black", cex=2)

lines(OphiSC_sg[,1], x = .5:25.5, lty=2, col='red')
arrows(.5:25.5,OphiSC_sg[,3],.5:25.5,OphiSC_sg[,2],lwd = 1.5, angle = 90,code = 3, length = 0.05, col='grey30')
points(OphiSC_sg[,1], x = .5:25.5, pch=21, bg="red", cex=1)
#cor(OphiSC[,1], OphiSC_sg[,1]) #0.9926979


## Plot cell
ophiAtlantSG <- unique(ophiAtlant[,1:3])
ophiPoints <- ophiAtlantSG
coordinates(ophiPoints) <- ~decimalLongitude+decimalLatitude

SampGrid <- AtlGrid
SampCell <- data.frame(SC=rep(NA,length(SampGrid)), Lat=rep(NA,length(SampGrid)))
for(i in 1:length(SampGrid))
{
  # For each grid cell...
  singleCell <- spPolygons(SampGrid@polygons[[i]]@Polygons[[1]]@coords)
  
  # Find all of its sampling events 
  ii <- which(!is.na(over(x = ophiPoints, singleCell)))
  if(length(ii)>0)
  {
    tempMat <- unique(ophiAtlantSG[ii,])
    
    # Used to avoid duplicating points in the latitude border (this does not occur with the longitude because grid cell's longitude has five digits)
    cats <- cut(tempMat[,3], breaks = seq(-65, 65, 5))
    realCat <- names(which(table(cats)==max(table(cats))))
    tempMat <- tempMat[which(cats==realCat),]
    
    # Use this sub-matrix to calculate sample coverage
    SampUn <- nrow(unique(tempMat[,c(2,3)]))
    regist <- table(as.character(tempMat[,1]))
    OphiINEXT <- c(SampUn, as.numeric(regist))
    
    # Only for those grid cells with at least 6 samples and number of incidences > number of unique species
    if(OphiINEXT[1]>5 & any(OphiINEXT[-1] > 1))
    {
      out <- DataInfo(OphiINEXT, datatype="incidence_freq")
      SampCell$SC[i] <- out$SC
    }
    SampCell$Lat[i] <- realCat
  }
}

# Average values by latitude
levs <- rev(levels(cut(ophiAtlant$decimalLatitude, breaks = seq(-65, 65, 5))))

OphiSCCell_sg <- numeric()
for(i in 1:length(levs))
{
  pos <- which(SampCell$Lat==rev(levs)[i])
  OphiSCCell_sg <- rbind(OphiSCCell_sg, c(mean(SampCell$SC[pos], na.rm = T), sd(SampCell$SC[pos], na.rm = T)))
}


frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=1.4, lwd=2, padj = .5)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.4, lwd=2)

polygon(c(.5:9.5, rev(.5:9.5)), c(OphiSCCell[1:10,1]+OphiSCCell[1:10,2], rev(OphiSCCell[1:10,1]-OphiSCCell[1:10,2])), col=rgb(255,240,216, maxColorValue = 255), border=NA)
polygon(c((11.5:25.5), rev(11.5:25.5)), c(OphiSCCell[12:26,1]+OphiSCCell[12:26,2], rev(OphiSCCell[12:26,1]-OphiSCCell[12:26,2])), col=rgb(255,240,216, maxColorValue = 255), border=NA)
lines(OphiSCCell[,1], x = .5:25.5, lty=2)
points(OphiSCCell[,1], x = .5:25.5, pch=21, bg=rgb(253,179,56, maxColorValue = 255), col="black", cex=2)

lines(OphiSCCell_sg[,1], x = .5:25.5, lty=2, col='red')
arrows(.5:25.5,OphiSCCell_sg[,1]+OphiSCCell_sg[,2],.5:25.5,OphiSCCell_sg[,1]-OphiSCCell_sg[,2],lwd = 1.5, angle = 90,code = 3, length = 0.05, col='grey30')
points(OphiSCCell_sg[,1], x = .5:25.5, pch=21, bg="red", cex=1)
#cor(OphiSCCell[,1], OphiSCCell_sg[,1], use='complete.obs') #0.9366629




# . Figure S3 ----
MatObsLat <- read.table("./Data/SimOutput/MatObsLat_10000.txt", h=F)

resMat <- numeric()
for(i in 3:ncol(MatObsLat))
{
  vec <- apply(MatObsLat[,1:i], 1, sd)
  resMat <- cbind(resMat, vec)
  
  if((i%%100)==0)
  {
    print(i)
  }
}
range(resMat)


# Plot
frame()
plot.window(xlim=c(0,10000), ylim=c(0,25))
axis(side=1, at=seq(0,10000,2000), labels=seq(0,10,2))
axis(side=2, at=seq(0,25,5), labels=seq(0,25,5), las=2)

color <- c(plasma(nrow(CellLatMat)/2), rev(plasma(nrow(CellLatMat)/2)))
for(i in 1:nrow(CellLatMat))
{
  lines(1:ncol(resMat), resMat[i,], col=color[i], lwd=2)
}
abline(v=2000, lty=2)




# . Figure S4 ----
myRd <- colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee090'))
myBu <- colorRampPalette(c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
color <- c(myRd(7)[-7],'#ffdb58', myBu(7)[-1]) #e5e5ab #ffdb58
color <- c(rev(color), color)

## Latitude (choose the variable to plot)
MatRichLat <- read.table("./Data/SimOutput/RSFD39_gcN/MatRichLat.txt", h=F, na.strings = -1)
MatChao2Lat <- read.table("./Data/SimOutput/RSFD39_gcN/MatChao2Lat.txt", h=F, na.strings = -1)
MatInExtLat <- read.table("./Data/SimOutput/RSFD39_gcN/MatInExtLat.txt", h=F, na.strings = -1)
MatES50Lat <- read.table("./Data/SimOutput/RSFD39_gcN/MatES50Lat.txt", h=F, na.strings = -1)

estMat <- MatRichLat
estMat <- MatChao2Lat
estMat <- MatInExtLat
estMat <- MatES50Lat

estSTD <- rowMeans(estMat, na.rm=T)
estSTD <- apply(cbind(estSTD/max(estSTD, na.rm=T), apply(estMat, 1, sd, na.rm=T)/max(estSTD, na.rm=T)), 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

arrows(1:26, estSTD[,1], 1:26, estSTD[,1]-estSTD[,2],lwd = 1.5, angle = 90,code = 3, length = 0.025)
lines(estSTD[,1], lwd=2, lty=2)
points(estSTD[,1], pch=21, bg=rgb(.45,.45,.45), col="black", cex=2.5, lwd=1.5)


## Grid cell (average)
MatRichCell <- read.table("./Data/SimOutput/RSFD39_gcN/MatRichCell.txt", h=F, na.strings = -1)
MatChao2Cell <- read.table("./Data/SimOutput/RSFD39_gcN/MatChao2Cell.txt", h=F, na.strings = -1)
MatInExtCell <- read.table("./Data/SimOutput/RSFD39_gcN/MatInExtCell.txt", h=F, na.strings = -1)
MatES50Cell <- read.table("./Data/SimOutput/RSFD39_gcN/MatES50Cell.txt", h=F, na.strings = -1)

estMat <- MatRichCell
estMat <- MatChao2Cell
estMat <- MatInExtCell
estMat <- MatES50Cell

estSTD <- grid2lat(estMat = estMat, CellLatMat = CellLatMat)
estSTD <- apply(estSTD/max(estSTD[,1], na.rm = T), 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

arrows(1:26, estSTD[,1], 1:26, estSTD[,1]-estSTD[,2], lwd=2.5, angle=90, code=3, length=0.03)
lines(estSTD[,1], lwd=2, lty=2)
points(estSTD[,1], pch=21, bg=rgb(.45,.45,.45), col="black", cex=2.5, lwd=1.5)


## Grid cell (same variable used before)
estSTD <- rowMeans(estMat, na.rm=T)
estSTD <- estSTD/max(estSTD, na.rm=T)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

for(i in nrow(CellLatMat):1)
{
  rowPos <- (CellLatMat[i,2]+1):(CellLatMat[i,2]+CellLatMat[i,3])
  points(x=rep(27-i, length(rowPos)), y=estSTD[rowPos], pch=21, bg=color[i], col="black", cex=2.5, lwd=.5)
}




# . Figure S5 ----
## Calculate the observed range size frequency distribution (RSFD) of the species

# This requires the raw data set of occurrence records
#load("./Data/Empirical/Raw/Ophiuroids.RData")
species <- unique(ophiAtlant$species)
rangeSize <- numeric()

pb.MX <- txtProgressBar(min = 0, max = length(species), style = 3)
for(i in 1:length(species))
{
  # Define the species
  x <- which(ophiAtlant$species==species[i])
  sp1 <-ophiAtlant[x,]
  
  # If there are more than on point of occurrence...
  if(nrow(sp1)>1)
  {
    # Create a convex hull polygon
    ch <- chull(sp1[,2:3])
    coords <- sp1[c(ch, ch[1]),2:3]
    poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)),1)))
    crs(poly) <- crs(AtlGrid)
    
    # Find the grid cells within the area potentially occupied by the species
    ii <- over(AtlGrid, poly)
    ii <- which(!is.na(ii))
    
    # Count the number of cells (1 if occurrence outside the grid)
    rangeSize[i] <- ifelse(length(ii)==0, 1, length(ii))
  }
  else
  {
    rangeSize[i] <- 1
  }
  setTxtProgressBar(pb.MX, i)
}


### Compare the observed and simulated RSFD
frame()
plot.window(xlim=c(0,4), ylim=c(0,350))
axis(side=1, at=c(1:4), labels=NA, tck=-.015, pos=0)
text(x=seq(1.5, 3.5, 1), y=par("usr")[3], labels=c('Obs',expression(paste(alpha,'=0.3,',beta,'=5')),expression(paste(alpha,'=3,',beta,'=9'))), cex=.7, xpd=T, adj=c(.5,1))
axis(side=2, at=seq(0,350,50), labels=NA, las=1, tck=-.01, pos=1)
text(x=par("usr")[1]+1.1, y=seq(0,350,50), labels=seq(0,350,50), cex=.7, xpd=T, adj=c(1,.5))

a <- boxplot(rangeSize, horizontal = F, add=T, col="deepskyblue", outpch=21, outbg="deepskyblue", frame.plot=F, axes=F, at=1.5)

alfa <- 1/100
for(i in 1:100)
{
  x <- ceiling(rbeta(675,.3,5)*396)
  
  b <- boxplot(x, horizontal = F, at=2, add=T, frame.plot=F, axes=F, col="deepskyblue", outbg="deepskyblue", plot=F)
  boxplot(x, horizontal = F, at=2.5, add=T, frame.plot=F, axes=F, outline=F, col=rgb(red = .5, green = .5, blue = .5, alpha = alfa), border=rgb(red = 0, green = 0, blue = 0, alpha = alfa), medcol=rgb(red = 1, green = 0, blue = 0, alpha = alfa))
  points(x=rep(2.5, length(unique(b$out))), y=unique(b$out), pch=21, bg=rgb(red = .5, green = .5, blue = .5, alpha = alfa), col=rgb(red = 0, green = 0, blue = 0, alpha = alfa))
}

for(i in 1:100)
{
  x <- ceiling(rbeta(675,3,9)*396)
  
  b <- boxplot(x, horizontal = F, at=3.5, add=T, frame.plot=F, axes=F, col="deepskyblue", outbg="deepskyblue", plot=F)
  boxplot(x, horizontal = F, at=3.5, add=T, frame.plot=F, axes=F, outline=F, col=rgb(red = .5, green = .5, blue = .5, alpha = alfa), border=rgb(red = 0, green = 0, blue = 0, alpha = alfa), medcol=rgb(red = 1, green = 0, blue = 0, alpha = alfa))
  points(x=rep(3.5, length(unique(b$out))), y=unique(b$out), pch=21, bg=rgb(red = .5, green = .5, blue = .5, alpha = alfa), col=rgb(red = 0, green = 0, blue = 0, alpha = alfa))
}




# . Figure S6 ----

## barplot Alpha
MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat/MatRichCell.txt", h=F, na.strings = -1)
MatObsCelly <- read.table("./Data/SimOutput/RSFD_nat/MatObsCell.txt", h=F, na.strings = c(-1,0))

MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop99/MatRichCell.txt", h=F, na.strings = -1)
MatObsCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop99/MatObsCell.txt", h=F, na.strings = c(-1,0))

MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatRichCell.txt", h=F, na.strings = -1)
MatObsCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatObsCell.txt", h=F, na.strings = c(-1,0))

avRichy <- grid2lat(estMat = MatRichCelly, CellLatMat = CellLatMat)[,1]
avObsy <- grid2lat(estMat = MatObsCelly, CellLatMat = CellLatMat)[,1]
maxV <- max(avRichy/max(avObsy))

x <- barplot((rev(avRichy)/max(avObsy))/maxV, col=rgb(.95,.95,.95), space=0, ylim=c(0,1), axes=F)
barplot((rev(avObsy)/max(avObsy))/maxV, col=rgb(.6,.6,.6), space=0, axes=F, add=T)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.4)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.4)

lines(x = (1:26)-.5, (SPrichCell[,1]/max(SPrichCell[,1]))/maxV, lty=2)
points(x = (1:26)-.5, (SPrichCell[,1]/max(SPrichCell[,1]))/maxV, pch=21, bg=rgb(.35,.35,.35, maxColorValue = 1), col="black", cex=1.5)


## Completeness Gamma
SCLaty <- read.table("./Data/SimOutput/RSFD_nat/MatCoverageLat.txt", h=F, na.strings = -1)
SCCelly <- read.table("./Data/SimOutput/RSFD_nat/MatCoverageCell.txt", h=F, na.strings = -1)

SCLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop99/MatCoverageLat.txt", h=F, na.strings = -1)
SCCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop99/MatCoverageCell.txt", h=F, na.strings = -1)

SCLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatCoverageLat.txt", h=F, na.strings = -1)
SCCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatCoverageCell.txt", h=F, na.strings = -1)


meanSD1 <- cbind(rowMeans(SCLaty, na.rm=T), apply(SCLaty, 1, sd, na.rm=T))
meanSD1 <- apply(meanSD1, 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.4)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.4)

polygon(c(.5:25.5, rev(.5:25.5)), c(meanSD1[,1]+meanSD1[,2], rev(meanSD1[,1]-meanSD1[,2])), col=rgb(.35,.35,.35,.2, maxColorValue = 1), border=NA)
lines(meanSD1[,1], x = x, lty=2)
points(meanSD1[,1], x = x, pch=22, bg=rgb(.35,.35,.35, maxColorValue = 1), col="black", cex=2)

polygon(c(.5:25.5, rev(.5:25.5)), c(OphiSC[,3], rev(OphiSC[,2])), col=rgb(2,81,150,50, maxColorValue = 255), border=NA)
lines(OphiSC[,1], x = x, lty=2)
points(OphiSC[,1], x = x, pch=21, bg=rgb(2,81,150, maxColorValue = 255), col="black", cex=1.5)


## Completeness Alpha
meanSD2 <- grid2lat(estMat = SCCelly, CellLatMat = CellLatMat)
meanSD2 <- apply(meanSD2, 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.4)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.4)

polygon(c(.5:25.5, rev(.5:25.5)), c(meanSD2[,1]+meanSD2[,2], rev(meanSD2[,1]-meanSD2[,2])), col=rgb(.35,.35,.35,.2, maxColorValue = 1), border=NA)
lines(meanSD2[,1], x = x, lty=2)
points(meanSD2[,1], x = x, pch=22, bg=rgb(.35,.35,.35, maxColorValue = 1), col="black", cex=2)

polygon(c(.5:9.5, rev(.5:9.5)), c(OphiSCCell[1:10,1]+OphiSCCell[1:10,2], rev(OphiSCCell[1:10,1]-OphiSCCell[1:10,2])), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
polygon(c(11.5:25.5, rev(11.5:25.5)), c(OphiSCCell[12:26,1]+OphiSCCell[12:26,2], rev(OphiSCCell[12:26,1]-OphiSCCell[12:26,2])), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
lines(OphiSCCell[,1], x = x, lty=2)
points(OphiSCCell[,1], x = x, pch=21, bg=rgb(253,179,56, maxColorValue = 255), col="black", cex=1.5)




# . Figure S7 ----
myRd <- colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee090'))
myBu <- colorRampPalette(c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
color <- c(myRd(7)[-7],'#ffdb58', myBu(7)[-1]) #e5e5ab #ffdb58
color <- c(rev(color), color)

## Latitude (choose the variable to plot)
MatCoverageLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatCoverageLat.txt", h=F, na.strings = -1)
estMat <- MatCoverageLat

MatCoverageCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatCoverageCell2Lat.txt", h=F, na.strings = -1)
estMat <- MatCoverageCell

MatGapsLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatGapsLat.txt", h=F, na.strings = -1)
estMat <- MatGapsLat/max(MatGapsLat[,1], na.rm = T)

# Standardize variable and optimize the plot order based on changes between sampling scenarios
meanSTD <- estMat[,seq(1, ncol(estMat), 3)]
xn <- meanSTD[,ncol(meanSTD)]
x1 <- meanSTD[,1]
n <- order(xn-x1)

# Plot richness estimation at each latitude
frame()
plot.window(xlim=log10(c(25,6400)), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.5, lwd=1.5)
axis(side=1, at=log10(50*(2^(-1:7))), labels=c(0,50*(2^(0:7)))/100, cex.axis=1.5, lwd=1.5)
segments(x0 = log10(50*(2^(0:7))), x1 = log10(50*(2^(0:7))), y0 = -.05, y1 = 1, lty=2, lwd=1.5, col='grey')
segments(x0 = log10(50*(2^-1)), x1 = log10(50*(2^7)), y0 = seq(0,1,.2), y1 = seq(0,1,.2), lty=3, lwd=1.5, col='grey90')

for(i in 1:nrow(CellLatMat))
{
  lines(x=log10(50*(2^(-1:7))), y=meanSTD[n[i],], lwd=3, col=color[n[i]])
}
format(apply(meanSTD, 2, min, na.rm=T)*100, scientific = F, digits = 3)
format(apply(meanSTD, 2, max, na.rm=T)*100, scientific = F, digits = 3)


## Grid cell
MatCoverageCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatCoverageCell.txt", h=F, na.strings = -1)
estMat <- MatCoverageCell
meanSTD <- estMat[,seq(1, ncol(estMat), 3)]

# Plot richness estimation at each latitude (for ES50 set F)
frame()
plot.window(xlim=log10(c(25,6400)), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.5, lwd=1.5)
axis(side=1, at=log10(50*(2^(-1:7))), labels=c(0,50*(2^(0:7)))/100, cex.axis=1.5, lwd=1.5)
segments(x0 = log10(50*(2^(0:7))), x1 = log10(50*(2^(0:7))), y0 = -.05, y1 = 1, lty=2, lwd=1.5, col='grey')
segments(x0 = log10(50*(2^-1)), x1 = log10(50*(2^7)), y0 = seq(0,1,.2), y1 = seq(0,1,.2), lty=3, lwd=1.5, col='grey90')

for(i in 1:nrow(CellLatMat))
{
  pos <- (CellLatMat[i,2]+1):(CellLatMat[i,2]+CellLatMat[i,3])
  
  for(j in 1:length(pos))
  {
    lines(x=log10(50*(2^(-1:7))), y=meanSTD[pos[j],], lwd=1.5, col=color[i])
  }
}
format(apply(meanSTD, 2, min, na.rm=T)*100, scientific = F, digits = 3)




# . Figure S8 ----

### PartA: here and script '2 - InputMatrices.R'
MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatRichCell.txt", h=F, na.strings = -1)
MatObsCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatObsCell.txt", h=F, na.strings = c(-1,0))
x <- rowMeans(MatObsCelly, na.rm = T)

MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatRichCell.txt", h=F, na.strings = -1)
MatObsCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatObsCell.txt", h=F, na.strings = c(-1,0))
x <- rowMeans(MatObsCelly, na.rm = T)

# Create color scale
escala <- cut(x, seq(floor(min(x, na.rm=T)/5)*5, ceiling(max(x, na.rm=T)/5)*5, 5))
levels(escala)

colCode <- numeric(length(AtlGrid@data[,1]))
for(i in 1:length(levels(escala)))
{
  pos <- which(escala==levels(escala)[i])
  colCode[pos] <- viridis(length(levels(escala)))[i]
}

# Plot the map
plot(AtlGrid, col=colCode, border="transparent")
plot(AtlGrid[colCode=='0',], border="black", lwd=.2, add=T)
plot(landCut, col="grey80", border="transparent", add=T)
plot(as(extent(AtlGrid), "SpatialPolygons"), add=T)
segments(x0 = extent(AtlGrid)[1], x1 = extent(AtlGrid)[2], y0 = 0, y1 = 0, lty = 2)


### PartB:
## Gamma diversity
MatRichLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatRichLat.txt", h=F, na.strings = -1)
MatObsLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatObsLat.txt", h=F, na.strings = c(-1,0))
MatGapsLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatGapsLat.txt", h=F, na.strings = -1)

MatRichLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatRichLat.txt", h=F, na.strings = -1)
MatObsLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatObsLat.txt", h=F, na.strings = c(-1,0))
MatGapsLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatGapsLat.txt", h=F, na.strings = -1)


MatRichLatyST <- rowMeans(MatRichLaty)/max(rowMeans(MatObsLaty))
MatObsLatyST <- rowMeans(MatObsLaty)/max(rowMeans(MatObsLaty))
ExpectedST <- rowMeans(MatObsLaty+MatGapsLaty)/max(rowMeans(MatObsLaty))
maxV <- max(MatRichLatyST)

x <- barplot(rev(MatRichLatyST)/maxV, col=rgb(.95,.95,.95), space=0, ylim=c(0,1.05), axes=F)
barplot(rev(ExpectedST)/maxV, col=rgb(.8,.8,.8), space=0, axes=F, add=T)
barplot(rev(MatObsLatyST)/maxV, col=rgb(.6,.6,.6), space=0, axes=F, add=T)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.4)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.4)

lines(x = (1:26)-.5, (SPrich/max(SPrich))/maxV, lty=2)
points(x = (1:26)-.5, (SPrich/max(SPrich))/maxV, pch=21, bg=rgb(.35,.35,.35, maxColorValue = 1), col="black", cex=1.5)


## Alpha diversity
MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatRichCell.txt", h=F, na.strings = -1)
MatObsCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatObsCell.txt", h=F, na.strings = c(-1,0))

MatRichCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatRichCell.txt", h=F, na.strings = -1)
MatObsCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatObsCell.txt", h=F, na.strings = c(-1,0))


avRichy <- grid2lat(estMat = MatRichCelly, CellLatMat = CellLatMat)[,1]
avObsy <- grid2lat(estMat = MatObsCelly, CellLatMat = CellLatMat)[,1]
maxV <- max(avRichy/max(avObsy))

x <- barplot((rev(avRichy)/max(avObsy))/maxV, col=rgb(.95,.95,.95), space=0, ylim=c(0,1), axes=F)
barplot((rev(avObsy)/max(avObsy))/maxV, col=rgb(.6,.6,.6), space=0, axes=F, add=T)
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.4)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.4)

lines(x = (1:26)-.5, (SPrichCell[,1]/max(SPrichCell[,1]))/maxV, lty=2)
points(x = (1:26)-.5, (SPrichCell[,1]/max(SPrichCell[,1]))/maxV, pch=21, bg=rgb(.35,.35,.35, maxColorValue = 1), col="black", cex=1.5)


## Completeness
SCLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatCoverageLat.txt", h=F, na.strings = -1)
SCCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79/MatCoverageCell.txt", h=F, na.strings = -1)

SCLaty <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatCoverageLat.txt", h=F, na.strings = -1)
SCCelly <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc/MatCoverageCell.txt", h=F, na.strings = -1)

# Gamma
meanSD1 <- cbind(rowMeans(SCLaty, na.rm=T), apply(SCLaty, 1, sd, na.rm=T))
frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.4)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.4)

polygon(c(.5:25.5, rev(.5:25.5)), c(rev(meanSD1[,1]+meanSD1[,2]), meanSD1[,1]-meanSD1[,2]), col=rgb(.35,.35,.35,.2, maxColorValue = 1), border=NA)
lines(rev(meanSD1[,1]), x = x, lty=2)
points(rev(meanSD1[,1]), x = x, pch=22, bg=rgb(.35,.35,.35, maxColorValue = 1), col="black", cex=2)

polygon(c(.5:25.5, rev(.5:25.5)), c(OphiSC[,3], rev(OphiSC[,2])), col=rgb(2,81,150,50, maxColorValue = 255), border=NA)
lines(OphiSC[,1], x = x, lty=2)
points(OphiSC[,1], x = x, pch=21, bg=rgb(2,81,150, maxColorValue = 255), col="black", cex=1.5)

# Alpha
meanSD2 <- grid2lat(estMat = SCCelly, CellLatMat = CellLatMat)
meanSD2 <- apply(meanSD2, 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=1, at=c(x[1]-0.5,mean(x[c(13,14)]),x[26]+0.5), labels=c('65S','0','65N'), cex.axis=1.4)
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, pos=-.5, lwd=1.5, cex.axis=1.4)

polygon(c(.5:25.5, rev(.5:25.5)), c(meanSD2[,1]+meanSD2[,2], rev(meanSD2[,1]-meanSD2[,2])), col=rgb(.35,.35,.35,.2, maxColorValue = 1), border=NA)
lines(meanSD2[,1], x = x, lty=2)
points(meanSD2[,1], x = x, pch=22, bg=rgb(.35,.35,.35, maxColorValue = 1), col="black", cex=2)

which(is.na(OphiSCCell[,2]))
polygon(c(.5:9.5, rev(.5:9.5)), c(OphiSCCell[1:10,1]+OphiSCCell[1:10,2], rev(OphiSCCell[1:10,1]-OphiSCCell[1:10,2])), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
polygon(c(11.5:25.5, rev(11.5:25.5)), c(OphiSCCell[12:26,1]+OphiSCCell[12:26,2], rev(OphiSCCell[12:26,1]-OphiSCCell[12:26,2])), col=rgb(253,179,56,50, maxColorValue = 255), border=NA)
lines(OphiSCCell[,1], x = x, lty=2)
points(OphiSCCell[,1], x = x, pch=21, bg=rgb(253,179,56, maxColorValue = 255), col="black", cex=1.5)


### PartC:
myRd <- colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee090'))
myBu <- colorRampPalette(c('#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
color <- c(myRd(7)[-7],'#ffdb58', myBu(7)[-1]) #e5e5ab #ffdb58
color <- c(rev(color), color)

int_f <- function(x, mu1, mu2, sd1, sd2)
{
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}

## Latitude (choose the variable to plot)
MatRichLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatRichLat.txt", h=F, na.strings = -1)
MatChao2Lat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatChao2Lat.txt", h=F, na.strings = -1)
MatInExtLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatInExtLat.txt", h=F, na.strings = -1)
MatES50Lat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatES50Lat.txt", h=F, na.strings = -1)

estMat <- MatChao2Lat
estMat <- MatInExtLat
estMat <- MatES50Lat

# Plot estimate from first sampling scenario
estSTD <- estMat[,1]
estSTD <- apply(cbind(estSTD/max(estSTD, na.rm=T), estMat[,2]/max(estSTD, na.rm=T)), 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

refMeanSD <- cbind(rowMeans(MatRichLat), apply(MatRichLat, 1, sd))
refMeanSD <- refMeanSD/max(refMeanSD[,1])
polygon(c(.5:25.5, rev(.5:25.5)), c(rev(refMeanSD[,1]+refMeanSD[,2]), refMeanSD[,1]-refMeanSD[,2]), col=rgb(222,222,222, maxColorValue = 255), border=NA)

arrows(1:26, estSTD[,1], 1:26, estSTD[,1]-estSTD[,2],lwd = 1.5, angle = 90,code = 3, length = 0.025)
lines(estSTD[,1], lwd=2, lty=2)
points(estSTD[,1], pch=21, bg=rgb(.45,.45,.45), col="black", cex=2.5, lwd=1.5)


# Standardize variable and optimize the plot order based on changes between sampling scenarios
refMean <- rowMeans(MatRichLat)
refSD <- apply(MatRichLat, 1, sd)

mu1 <- refMean/max(refMean)
sd1 <- refSD/max(refMean)

matMean <- estMat[,seq(1, ncol(estMat), 3)]
matSD <- estMat[,seq(2, ncol(estMat), 3)]

xn <- matMean[,ncol(matMean)]/apply(matMean, 2, max)[ncol(matMean)]
x1 <- matMean[,1]/apply(matMean, 2, max)[1]
n <- order(xn-x1)

# Plot richness estimation at each latitude
frame()
plot.window(xlim=log10(c(25,6400)), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.5, lwd=1.5)
axis(side=1, at=log10(50*(2^(-1:7))), labels=c(0,50*(2^(0:7)))/100, cex.axis=1.5, lwd=1.5)
segments(x0 = log10(50*(2^(-1:7))), x1 = log10(50*(2^(-1:7))), y0 = 0, y1 = 1, lty=2, lwd=1.5, col='grey')
segments(x0 = log10(50*(2^-1)), x1 = log10(50*(2^7)), y0 = seq(0,1,.2), y1 = seq(0,1,.2), lty=3, lwd=1.5, col='grey90')

for(j in 1:nrow(CellLatMat))
{
  overlVec <- numeric()
  for(i in 1:ncol(matMean))
  {
    mu2 <- matMean[,i]/max(matMean[,i], na.rm=T)
    sd2 <- matSD[,i]/max(matMean[,i], na.rm=T)
    
    if(!is.na(sd2[n[j]]))
    {
      overlVec[i] <- integrate(int_f, -Inf, Inf, mu1=mu1[n[j]], mu2=mu2[n[j]], sd1=sd1[n[j]], sd2=sd2[n[j]])[[1]]
    }
  }
  lines(x=log10(50*(2^(-1:7))), y=overlVec, lwd=3, col=color[n[j]])
}


## Grid cell (average)
MatRichCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatRichCell.txt", h=F, na.strings = -1)
MatChao2Cell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatChao2Cell2Lat.txt", h=F, na.strings = -1)
MatInExtCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatInExtCell2Lat.txt", h=F, na.strings = -1)
MatES50Cell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatES50Cell2Lat.txt", h=F, na.strings = -1)

# Plot estimate from first sampling scenario
estMat <- MatChao2Cell
estMat <- MatInExtCell
estMat <- MatES50Cell

estSTD <- estMat[,1:2]
estSTD[estSTD==0] <- NA
estSTD <- apply(estSTD/max(estSTD[,1], na.rm = T), 2, rev)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

refData <- grid2lat(estMat = MatRichCell, CellLatMat = CellLatMat)
refData <- refData/max(refData[,1])
polygon(c(.5:25.5, rev(.5:25.5)), c(rev(refData[,1]+refData[,2]), refData[,1]-refData[,2]), col=rgb(222,222,222, maxColorValue = 255), border=NA)

arrows(1:26, estSTD[,1], 1:26, estSTD[,1]-estSTD[,2], lwd=2.5, angle=90, code=3, length=0.03)
lines(estSTD[,1], lwd=2, lty=2)
points(estSTD[,1], pch=21, bg=rgb(.45,.45,.45), col="black", cex=2.5, lwd=1.5)


# Standardize variable and optimize the plot order based on changes between sampling scenarios
refMeanSD <- grid2lat(estMat = MatRichCell, CellLatMat = CellLatMat)

mu1 <- refMeanSD[,1]/max(refMeanSD[,1])
sd1 <- refMeanSD[,2]/max(refMeanSD[,1])

meanLat <- estMat[,seq(1, ncol(estMat), 3)]
sdLat <- estMat[,seq(2, ncol(estMat), 3)]

xn <- meanLat[,ncol(meanLat)]/apply(meanLat, 2, max)[ncol(meanLat)]
x1 <- meanLat[,1]/apply(meanLat, 2, max, na.rm=T)[1]
n <- order(xn-x1)

# Plot richness estimation at each latitude
frame()
plot.window(xlim=log10(c(25,6400)), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=1.5, lwd=1.5)
axis(side=1, at=log10(50*(2^(-1:7))), labels=c(0,50*(2^(0:7)))/100, cex.axis=1.5, lwd=1.5)
segments(x0 = log10(50*(2^(-1:7))), x1 = log10(50*(2^(-1:7))), y0 = 0, y1 = 1, lty=2, lwd=1.5, col='grey')
segments(x0 = log10(50*(2^-1)), x1 = log10(50*(2^7)), y0 = seq(0,1,.2), y1 = seq(0,1,.2), lty=3, lwd=1.5, col='grey90')

for(lat in 1:nrow(CellLatMat))
{
  overlVec <- numeric()
  for(i in 1:ncol(meanLat))
  {
    mu2 <- meanLat[,i]/max(meanLat[,i], na.rm=T)
    sd2 <- sdLat[,i]/max(meanLat[,i], na.rm=T)
    
    if(!is.na(sd2[n[lat]]))
    {
      overlVec[i] <- integrate(int_f, -Inf, Inf, mu1=mu1[n[lat]], mu2=mu2[n[lat]], sd1=sd1[n[lat]], sd2=sd2[n[lat]])[[1]]
    }
    else
    {
      overlVec[i] <- NA
    }
  }
  lines(x=log10(50*(2^(-1:7))), y=overlVec, lwd=3, col=color[n[lat]])
}




# . Figure S9 ----

# Color pallete
myPalette <- c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')


### Latitude
# Open files
MatRichLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatRichLat.txt", h=F, na.strings = -1)
MatChao2Lat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatChao2Lat.txt", h=F, na.strings = -1)
MatInExtLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatInExtLat.txt", h=F, na.strings = -1)
MatES50Lat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatES50Lat.txt", h=F, na.strings = -1)

MatRichLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatRichLat.txt", h=F, na.strings = -1)
MatChao2Lat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatChao2Lat.txt", h=F, na.strings = -1)
MatInExtLat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatInExtLat.txt", h=F, na.strings = -1)
MatES50Lat <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatES50Lat.txt", h=F, na.strings = -1)


## Chao2 estimate
# This requires the raw data set of occurrence records
#load("./Data/Empirical/Raw/Ophiuroids.RData")
OphiINEXT <- chaoInput(ophiAtlant)

chao2 <- numeric()
for(i in 1:length(OphiINEXT))
{
  chao2[i] <- ChaoRichness(x = OphiINEXT[[i]], datatype = 'incidence_freq')[1,2]
}

# Define the reliability of the estimates
refMean <- rowMeans(MatRichLat)
refSD <- apply(MatRichLat, 1, sd)
mu1 <- refMean/max(refMean)
sd1 <- refSD/max(refMean)

matMean <- MatChao2Lat[,1]
matSD <- MatChao2Lat[,2]
mu2 <- matMean/max(matMean)
sd2 <- matSD/max(matMean)

overlVec <- numeric()
for(j in 1:nrow(CellLatMat))
{
  overlVec[j] <- integrate(int_f, -Inf, Inf, mu1=mu1[j], mu2=mu2[j], sd1=sd1[j], sd2=sd2[j])[[1]]
}
overlVec <- rev(overlVec)
mycol <- myPalette[cut(overlVec, breaks = seq(0,1,.1))]

# Plot
frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

lines(chao2/max(chao2), lty=2)
points(chao2/max(chao2), pch=21, bg=mycol, cex=2.5)


## iNext estimate
Csize <- numeric()
for(i in 1:length(OphiINEXT))
{
  x <- estimateD(x = OphiINEXT[[i]], datatype = 'incidence_freq', base='size', level = OphiINEXT[[i]][1]*2, q = 0, nboot=0)
  Csize[i] <- x$SC
  print(i)
}

inext <- numeric()
for(i in 1:length(OphiINEXT))
{
  x <- estimateD(x = OphiINEXT[[i]], datatype = 'incidence_freq', base='coverage', level = min(Csize), q = 0, nboot=0)
  inext[i] <- x$qD
  print(x$Method)
}

# Define the reliability of the estimates
matMean <- MatInExtLat[,1]
matSD <- MatInExtLat[,2]
mu2 <- matMean/max(matMean)
sd2 <- matSD/max(matMean)

overlVec <- numeric()
for(j in 1:nrow(CellLatMat))
{
  overlVec[j] <- integrate(int_f, -Inf, Inf, mu1=mu1[j], mu2=mu2[j], sd1=sd1[j], sd2=sd2[j])[[1]]
}
overlVec <- rev(overlVec)
mycol <- myPalette[cut(overlVec, breaks = seq(0,1,.1))]

# Plot
frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

lines(inext/max(inext), lty=2)
points(inext/max(inext), pch=21, bg=mycol, cex=2.5)


## ES50 estimate
es50SAM <- numeric()
for(i in 1:length(OphiINEXT))
{
  x=OphiINEXT[[i]][-1]
  sample=50
  x <- x[x > 0]
  J <- sum(x)
  J <- OphiINEXT[[i]][1]
  ldiv <- lchoose(J, sample)
  p1 <- ifelse(J - x < sample, 0, exp(lchoose(J - x, sample) - ldiv))
  
  if(J>=sample){es50SAM[i] <- sum(1 - p1)}
}

# Define the reliability of the estimates
matMean <- MatES50Lat[,1]
matSD <- MatES50Lat[,2]
mu2 <- matMean/max(matMean, na.rm = T)
sd2 <- matSD/max(matMean, na.rm = T)

overlVec <- numeric()
for(j in 1:nrow(CellLatMat))
{
  if(!is.na(mu2[j]))
  {
    overlVec[j] <- integrate(int_f, -Inf, Inf, mu1=mu1[j], mu2=mu2[j], sd1=sd1[j], sd2=sd2[j])[[1]]
  }
}
overlVec <- rev(overlVec)
mycol <- myPalette[cut(overlVec, breaks = seq(0,1,.1))]

# Plot
frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

lines(es50SAM/max(es50SAM, na.rm = T), lty=2)
points(es50SAM/max(es50SAM, na.rm = T), pch=21, bg=mycol, cex=2.5)



### Grid cell

# Estimates
ophiPoints <- ophiAtlant
coordinates(ophiPoints) <- ~decimalLongitude+decimalLatitude

SampGrid <- AtlGrid
estimCell <- data.frame(chao2=rep(NA,length(SampGrid)), SCdup=rep(NA,length(SampGrid)), inExt=rep(NA,length(SampGrid)), ES50=rep(NA,length(SampGrid)), Lat=rep(NA,length(SampGrid)))
for(i in 1:length(SampGrid))
{
  # For each grid cell...
  singleCell <- spPolygons(SampGrid@polygons[[i]]@Polygons[[1]]@coords)
  
  # Find all of its sampling events 
  ii <- which(!is.na(over(x = ophiPoints, singleCell)))
  if(length(ii)>0)
  {
    tempMat <- unique(ophiAtlant[ii,1:4])
    
    # Used to avoid duplicating points in the latitude border (this does not occur with the longitude because it have five digits)
    cats <- cut(tempMat[,3], breaks = seq(-65, 65, 5))
    realCat <- names(which(table(cats)==max(table(cats))))
    tempMat <- tempMat[which(cats==realCat),]
    
    # Use this sub-matrix to calculate sample coverage
    SampUn <- nrow(unique(tempMat[,c(2,3,4)]))
    regist <- table(as.character(tempMat[,1]))
    OphiINEXT <- c(SampUn, as.numeric(regist))
    
    # Only for those grid cells with at least 6 samples and number of incidences > number of unique species
    if(OphiINEXT[1]>5 & any(OphiINEXT[-1] > 1))
    {
      estimCell$chao2[i] <- ChaoRichness(x = OphiINEXT, datatype = 'incidence_freq')[1,2]
      
      x <- estimateD(x = OphiINEXT, datatype = 'incidence_freq', base='size', level = OphiINEXT[1]*2, q = 0, nboot=0)
      estimCell$SCdup[i] <- x$SC
    }
    
    if(OphiINEXT[1] >= 50)
    {
      x=OphiINEXT[-1]
      sample=50
      x <- x[x > 0]
      J <- OphiINEXT[1]
      ldiv <- lchoose(J, sample)
      p1 <- ifelse(J - x < sample, 0, exp(lchoose(J - x, sample) - ldiv))
      
      estimCell$ES50[i] <- sum(1 - p1)
    }
    
    estimCell$Lat[i] <- realCat
  }
  print(i)
}

for(i in 1:length(SampGrid))
{
  # For each grid cell...
  singleCell <- spPolygons(SampGrid@polygons[[i]]@Polygons[[1]]@coords)
  
  # Find all of its sampling events 
  ii <- which(!is.na(over(x = ophiPoints, singleCell)))
  if(length(ii)>0)
  {
    tempMat <- unique(ophiAtlant[ii,1:4])
    
    # Used to avoid duplicating points in the latitude border (this does not occur with the longitude because it have five digits)
    cats <- cut(tempMat[,3], breaks = seq(-65, 65, 5))
    realCat <- names(which(table(cats)==max(table(cats))))
    tempMat <- tempMat[which(cats==realCat),]
    
    # Use this sub-matrix to calculate sample coverage
    SampUn <- nrow(unique(tempMat[,c(2,3,4)]))
    regist <- table(as.character(tempMat[,1]))
    OphiINEXT <- c(SampUn, as.numeric(regist))
    
    # Only for those grid cells with at least 6 samples and number of incidences > number of unique species
    if(OphiINEXT[1]>5 & any(OphiINEXT[-1] > 1))
    {
      x <- estimateD(x = OphiINEXT, datatype = 'incidence_freq', base='coverage', level = min(estimCell$SCdup, na.rm = T), q = 0, nboot = 0)
      estimCell$inExt[i] <- x$qD
    }
  }
  print(i)
}

# Average values by latitude
levs <- rev(levels(cut(ophiAtlant$decimalLatitude, breaks = seq(-65, 65, 5))))

chao2 <- numeric()
inExt <- numeric()
es50 <- numeric()
for(i in 1:length(levs))
{
  pos <- which(estimCell$Lat==rev(levs)[i])
  
  chao2 <- rbind(chao2, c(mean(estimCell$chao2[pos], na.rm = T), sd(estimCell$chao2[pos], na.rm=T)))
  inExt <- rbind(inExt, c(mean(estimCell$inExt[pos], na.rm = T), sd(estimCell$inExt[pos], na.rm=T)))
  es50 <- rbind(es50, c(mean(estimCell$ES50[pos], na.rm = T), sd(estimCell$ES50[pos], na.rm=T)))
}

# Open files
MatRichCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatRichCell.txt", h=F, na.strings = -1)
MatChao2Cell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatChao2Cell2Lat.txt", h=F, na.strings = -1)
MatInExtCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatInExtCell2Lat.txt", h=F, na.strings = -1)
MatES50Cell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seINC/MatES50Cell2Lat.txt", h=F, na.strings = -1)

MatRichCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatRichCell.txt", h=F, na.strings = -1)
MatChao2Cell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatChao2Cell2Lat.txt", h=F, na.strings = -1)
MatInExtCell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatInExtCell2Lat.txt", h=F, na.strings = -1)
MatES50Cell <- read.table("./Data/SimOutput/RSFD_nat_Trop79_seConc_seINC/MatES50Cell2Lat.txt", h=F, na.strings = -1)


# Define the reliability of the estimates
refMeanSD <- grid2lat(estMat = MatRichCell, CellLatMat = CellLatMat)
mu1 <- refMeanSD[,1]/max(refMeanSD[,1])
sd1 <- refMeanSD[,2]/max(refMeanSD[,1])

estMat <- MatChao2Cell
estReal <- chao2

matMean <- estMat[,1]
matSD <- estMat[,2]
mu2 <- matMean/max(matMean, na.rm = T)
sd2 <- matSD/max(matMean, na.rm = T)
sd2[is.na(sd2)] <- 0

overlVec <- numeric()
for(j in 1:nrow(CellLatMat))
{
  if(!is.na(mu2[j]))
  {
    overlVec[j] <- integrate(int_f, -Inf, Inf, mu1=mu1[j], mu2=mu2[j], sd1=sd1[j], sd2=sd2[j])[[1]]
  }
}
overlVec <- rev(overlVec)
mycol <- myPalette[cut(overlVec, breaks = seq(0,1,.1))]
mycol[is.na(mycol)] <- 'white'
  
# Plot
estReal <- estReal/max(estReal[,1], na.rm = T)

frame()
plot.window(xlim=c(1,26), ylim=c(0,1))
axis(side=2, at=seq(0,1,.2), labels=sprintf("%1.1f", seq(0,1,.2)), las=1, cex.axis=2, lwd=2)
axis(side=1, at=c(1-0.5,14-0.5,26+0.5), labels=c('65S','0','65N'), cex.axis=2, lwd=2, padj = .5)

lines(estReal[,1], lty=2)
points(estReal[,1], pch=21, bg=mycol, cex=2.5)




#### END ####