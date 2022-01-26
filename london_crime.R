#London Crime

library (sp)
library (raster)
library (spatstat)
library (maptools)
library (plotrix)
library(sf)
library(viridisLite)

setwd("D:/Study/Copernicus Master/Semester 2/Spatial Statistics/lesson 5/2018-06")
data <- read.csv("2018-06-metropolitan-street.csv")

head(data)
str(data)
data <- data[!is.na(data$Longitude)&!is.na(data$Latitude),]
row(data)

#Point Pattern analysis
#stochastic process for which we observe its results, or events, only in a specific region, part of spatstat
#has it own format (ppp)

#remove duplicate locations
coordinates(data)=~Longitude+Latitude
zero <- zerodist(data)
length(unique(zero[,1]))

#download shapefile area of study
download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip",destfile="ne_10m_admin_1_states_provinces.zip")
unzip("ne_10m_admin_1_states_provinces.zip",exdir="NaturalEarth")
setwd("D:/Study/Copernicus Master/Semester 2/Spatial Statistics/lesson 5")
border <- shapefile("ne_10m_admin_1_states_provinces.shp")

head(border)

#extract only Greater London

GreaterLondon <- border[paste(border$region)=="Greater London",]

typeof(GreaterLondon)

#overlay with crime data and eliminate point outside Greater London
projection(data)=projection(border)
overlay <- over(data,GreaterLondon)

data$over <- overlay$adm1_code
data.London <- data[!is.na(data$over),]

data.London$Crime.type



jpeg("PP_plot.jpg",2500,2000,res=300)
plot(data.London,pch="+",cex=0.5,main="Greater London Crimes per Type May 2018",col=terrain.colors(length(unique(data.London$Crime.type))))
plot(GreaterLondon,add=T)
legend(x=-0.53,y=51.41,pch="+",col=terrain.colors(length(unique(data.London$Crime.type))),legend=unique(data.London$Crime.type),cex=0.4)
dev.off()

jpeg("PP_plot_try.jpg",2500,2000,res=300)
plot(data.London,pch="+",cex=0.5,main="Greater London Crimes per Type June 2018",col= rep(1:length(unique(data.London$Crime.type))))
plot(GreaterLondon,add=T)
legend(x=-0.53,y=51.41,pch="+",col=turbo(length(unique(data.London$Crime.type))),legend=unique(data.London$Crime.type),cex=0.4)
dev.off()


#Descriptive Statistics 2d data

mean_centerX <- mean(data.London@coords[,1])
mean_centerY <- mean(data.London@coords[,2])

sdX <- sd(data.London@coords[,1])
sdY <- sd(data.London@coords[,2])

standard_distance <- sqrt(sum(((data.London@coords[,1]-mean_centerX)^2+(data.London@coords[,2]-mean_centerY)^2))/(nrow(data.London)))

jpeg("PP_Circle.jpeg",2500,2000,res=300)
plot(data.London,pch="+",cex=0.5,main="")
plot(GreaterLondon,add=T)
points(mean_centerX,mean_centerY,col="red",pch=16)
draw.circle(mean_centerX,mean_centerY,radius=standard_distance,border="red",lwd=2)
dev.off()

jpeg("PP_Ellipse_june.jpeg",2500,2000,res=300)
plot(data.London,pch="+",cex=0.5,main="PP Ellipse Greater London June 2018")
plot(GreaterLondon,add=T)
points(mean_centerX,mean_centerY,col="red",pch=16)
draw.ellipse(mean_centerX,mean_centerY,a=sdX,b=sdY,border="red",lwd=2)
dev.off()

jpeg("PP_circle_june.jpeg",2500,2000,res=300)
plot(data.London,pch="+",cex=0.5,main="PP Circlee Greater London June 2018")
plot(GreaterLondon,add=T)
points(mean_centerX,mean_centerY,col="red",pch=16)
draw.circle(mean_centerX,mean_centerY,radius=standard_distance,border="red",lwd=2)
dev.off()


#working with spatstat

#Cleaning up Public order values 
Public_order <- data.London[data.London$Crime.type==unique(data.London$Crime.type)[5],]
Public_order <- remove.duplicates(Public_order)
#change the CRS format into UTM format
Public_order <- spTransform(Public_order,CRS("+init=epsg:32630"))
Public_order

#change area of study CRS format into UTM
GreaterLondonUTM <- spTransform(GreaterLondon,CRS("+init=epsg:32630"))
window <- as.owin(GreaterLondonUTM)
window
Public_order.ppp <- ppp(x=Public_order@coords[,1],y=Public_order@coords[,2],window=window)
plot(Public_order.ppp)
Public_order.ppp$n/sum(sapply(slot(GreaterLondonUTM, "polygons"), slot, "area"))

jpeg("PP_QuadratCounting.jpeg",2500,2000,res=300)
plot(Public_order.ppp,pch="+",cex=0.5,main="Public Order")
plot(quadratcount(Public_order.ppp, nx = 4, ny = 4),add=T,col="blue")
dev.off()

#since divided in quadrat counting give meaningless area bonderies, we might change into borough bondary
Local.Intensity <- data.frame(Borough=factor(),Number=numeric())
for(i in unique(GreaterLondonUTM$name)){
  sub.pol <- GreaterLondonUTM[GreaterLondonUTM$name==i,]
  
  sub.ppp <- ppp(x=Public_order.ppp$x,y=Public_order.ppp$y,window=as.owin(sub.pol))
  Local.Intensity <- rbind(Local.Intensity,data.frame(Borough=factor(i,levels=GreaterLondonUTM$name),Number=sub.ppp$n))
}

colorScale <- color.scale(Local.Intensity[order(Local.Intensity[,2]),2],color.spec="rgb",extremes=c("green","red"),alpha=0.8)

jpeg("PP_BoroughCounting.jpeg",2000,2000,res=300)
par(mar=c(5,13,4,2)) 
barplot(Local.Intensity[order(Local.Intensity[,2]),2],names.arg=Local.Intensity[order(Local.Intensity[,2]),1],horiz=T,las=2,space=1,col=colorScale)
dev.off()

#bandwidth test

jpeg("Kernel_Density.jpeg",2500,2000,res=300)
par(mfrow=c(2,2))
plot(density.ppp(Public_order.ppp, sigma = bw.diggle(Public_order.ppp),edge=T),main=paste("h =",round(bw.diggle(Public_order.ppp),2)))
plot(density.ppp(Public_order.ppp, sigma = bw.ppl(Public_order.ppp),edge=T),main=paste("h =",round(bw.ppl(Public_order.ppp),2)))
plot(density.ppp(Public_order.ppp, sigma = bw.scott(Public_order.ppp)[2],edge=T),main=paste("h =",round(bw.scott(Public_order.ppp)[2],2)))
plot(density.ppp(Public_order.ppp, sigma = bw.scott(Public_order.ppp)[1],edge=T),main=paste("h =",round(bw.scott(Public_order.ppp)[1],2)))
dev.off()

#G function
jpeg("GFunction.jpeg",2500,2000, res=300)
plot(Gest(Public_order.ppp), main="Public order Related Crimes")
dev.off()

#Drugs data

#Cleaning up Drugs values 
data.London$Crime.type[3]
Drugs <- data.London[data.London$Crime.type==unique(data.London$Crime.type)[3],]
Drugs
Drugs <- remove.duplicates(Drugs)
#change the CRS format into UTM format
Drugs <- spTransform(Drugs,CRS("+init=epsg:32630"))
Drugs

#change area of study CRS format into UTM
GreaterLondonUTM <- spTransform(GreaterLondon,CRS("+init=epsg:32630"))
window <- as.owin(GreaterLondonUTM)
window
Drugs.ppp <- ppp(x=Drugs@coords[,1],y=Drugs@coords[,2],window=window)
plot(Drugs.ppp)
Drugs.ppp$n/sum(sapply(slot(GreaterLondonUTM, "polygons"), slot, "area"))

jpeg("PP_QuadratCounting_Drugs.jpeg",2500,2000,res=300)
plot(Drugs.ppp,pch="+",cex=0.5,main="Drugs June 2018")
plot(quadratcount(Drugs.ppp, nx = 4, ny = 4),add=T,col="blue")
dev.off()

#since divided in quadrat counting give meaningless area bonderies, we might change into borough bondary
Local.Intensity <- data.frame(Borough=factor(),Number=numeric())
for(i in unique(GreaterLondonUTM$name)){
  sub.pol <- GreaterLondonUTM[GreaterLondonUTM$name==i,]
  
  sub.ppp <- ppp(x=Drugs.ppp$x,y=Drugs.ppp$y,window=as.owin(sub.pol))
  Local.Intensity <- rbind(Local.Intensity,data.frame(Borough=factor(i,levels=GreaterLondonUTM$name),Number=sub.ppp$n))
}

colorScale <- color.scale(Local.Intensity[order(Local.Intensity[,2]),2],color.spec="rgb",extremes=c("green","red"),alpha=0.8)

jpeg("PP_BoroughCounting_Drugs.jpeg",2000,2000,res=300)
par(mar=c(5,13,4,2)) 
barplot(Local.Intensity[order(Local.Intensity[,2]),2],names.arg=Local.Intensity[order(Local.Intensity[,2]),1],horiz=T,las=2,space=1,col=colorScale)
dev.off()

#bandwidth test

jpeg("Kernel_Density_drugs.jpeg",2500,2000,res=300)
par(mfrow=c(2,2))
plot(density.ppp(Drugs.ppp, sigma = bw.diggle(Drugs.ppp),edge=T),main=paste("h =",round(bw.diggle(Drugs.ppp),2)))
plot(density.ppp(Drugs.ppp, sigma = bw.ppl(Drugs.ppp),edge=T),main=paste("h =",round(bw.ppl(Drugs.ppp),2)))
plot(density.ppp(Drugs.ppp, sigma = bw.scott(Drugs.ppp)[2],edge=T),main=paste("h =",round(bw.scott(Drugs.ppp)[2],2)))
plot(density.ppp(Drugs.ppp, sigma = bw.scott(Drugs.ppp)[1],edge=T),main=paste("h =",round(bw.scott(Drugs.ppp)[1],2)))
dev.off()

#G function
jpeg("GFunction_drugs.jpeg",2500,2000, res=300)
plot(Gest(Drugs.ppp), main="Drugs Related Crimes")
dev.off()
