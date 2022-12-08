## 12/8/2022
## figure 2 and 3
## download data from data folder first

library(sf)
library(tcltk)
library(raster)
library(rasterVis)
library(ggplot2)
library(rgdal)

setwd()

chn.shp = readOGR(dsn=yourdatafolder, layer="wc2.1_10m_prec.1970-2000")
forestarea.df = read.csv("forest_area_by_GridID2.csv")
forestarea.df$ForestArea[which(is.na(forestarea.df$ForestArea))] = 0

chn.shp@data = merge(chn.shp@data, forestarea.df, by.x = "ID", by.y = "GridID")

forestarea.grd = data.frame(chn.shp@data$x, chn.shp@data$y,chn.shp@data$ForestArea)
forestarea.grd = rasterFromXYZ(forestarea.grd)

chn.shp2 = readOGR(dsn=yourdatafolder, layer="CHN_adm1_wTaiwan")

# figure 2b - 2c
yrs = 1986:2020

mean.stack = stack()
median.stack = stack()
freq.stack = stack()

for (i in 1:length(yrs)){

size1 = data.frame(read.csv(paste0("yourdatafolder/patch_size_", yrs[i],".csv")))[,c(2:8)]
size2 = merge(chn.shp@data, size1, by.x = "ID", by.y = "GridID", all.x = TRUE)

tmp = data.frame(long = size2$x, lat = size2$y, mean1 = size2$mean/10000)
tmp2 = data.frame(long = size2$x, lat = size2$y, mean1 = size2$median/10000)
tmp4 = data.frame(long = size2$x, lat = size2$y, mean1 = size2$count/size2$ForestArea)

tmp.grd = rasterFromXYZ(tmp)
tmp2.grd = rasterFromXYZ(tmp2)
tmp4.grd = rasterFromXYZ(tmp4)

tmp.grd[tmp.grd == 0] = NA
tmp2.grd[tmp2.grd == 0] = NA
tmp4.grd[tmp4.grd == 0] = NA

mean.stack = addLayer(mean.stack, tmp.grd)
median.stack = addLayer(median.stack, tmp2.grd)
freq.stack = addLayer(freq.stack, tmp4.grd)


}

list.stack = list(mean.stack, median.stack, 1000*freq.stack)
list.title = list("Mean Patch Size (ha)", "Median Patch Size (ha)",  "Disturbance Frequence (# 1000 km-2)")

png.title = list("MeanPatchSize", "MedianPatchSize",  "DisturbanceFrequence")

# produce figures using a for loop
for (i in 1:3){

r1.grd = calc(list.stack[[i]], mean, na.rm = TRUE)
r1.grd[forestarea.grd < 5] = NA

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.01,0.99), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(min_, max_) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

png(paste0(yourdatafolder,png.title[[i]],".png"),height = 1500, width = 2000, res = 300, units = "px")

p.strip <- list(cex=1, lines=2, fontface='bold')
rasterVis::levelplot(r1.grd, 
          zlim = r1.range,
		  margin = FALSE,
		  main = list.title[[i]],
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
		  col.regions = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(255),
          # xlab=list(label = "Longtitude", cex=1),ylab=list(label = "Latitude", cex=1),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  # layout=c(5, 1),
		  colorkey=list(space="right", 
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) +  
		  latticeExtra::layer(sp.polygons(chn.shp2, col = "black", lwd = 2))  


dev.off()

# create histgram plot 
myData <- as.vector(as.matrix(r1.grd))
binwidth <- as.numeric((r1.range[2]-r1.range[1])/20)
   
library(ggplot2)   # CRAN version 2.2.1 used
n_bins <- length(ggplot2:::bin_breaks_width(range(myData, na.rm = TRUE), width = binwidth)$breaks) - 1L
ggplot() + geom_histogram(aes(x = myData), binwidth = binwidth, fill = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(n_bins)) + 
aes(y=100*stat(count)/sum(stat(count)))+
geom_vline(aes(xintercept = mean(myData, na.rm = TRUE)),col='black',size=2, lty = 2)+
theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("% of Forests Area") + 
  xlab(list.title[[i]]) +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") 

ggsave(paste0(yourdatafolder,png.title[[i]],"_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")

}

########### at trend ##################################
# install.packages("spatialEco")
library(spatialEco)

#length of NAs
lenNA = function(x){return(length(which(is.na(x))))}

# Write a function to calcuate the trend
Trend1 = function(stkgrd){
require(spatialEco)
stkgrd_lenNA=calc(stkgrd,lenNA)
stkgrd[stkgrd_lenNA > 10] = -9999 #replace NAs with -9999 to use the raster.kendall function
stkgrd_trend= raster.kendall(stkgrd,p.value = TRUE, na.rm = TRUE)
stkgrd_trend[stkgrd_lenNA > 10] = NA
return(stkgrd_trend)
}

# calculate trends: Kendall trend test 
mean_trend = Trend1(mean.stack)
median_trend = Trend1(median.stack)
freq_trend = Trend1(freq.stack)

# plot

list.stack2 = list(mean_trend, median_trend, 1000*freq_trend)
list.title2 = list("Trend for Mean Patch Size (ha yr-1)", "Trend for Median Patch Size (ha yr-1)",  "Trend for Disturbance Frequence (# 1000 km-2 yr-1)")

png.title2 = list("MeanPatchSize_trend", "MedianPatchSize_trend",  "DisturbanceFrequence_trend")

for (i in 1:3){
r1.grd = list.stack2[[i]][[1]]
r1.grd[forestarea.grd < 5] = NA

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.01,0.99), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(-max(abs(min_),abs(max_)),max(abs(min_),abs(max_))) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

# plot(r1.grd)
## extract p value < 0.05, and add to plot
# Ex.pts was from https://github.com/liuzh811/ForestDegradationWestAfrica/blob/master/function.r
pts.sp.sig1 = Ex.pts(list.stack2[[i]][[2]], sig.level = 0.05) #extract significant relation points

png(paste0(yourdatafolder,png.title2[[i]],".png"),height = 1500, width = 2000, res = 300, units = "px")

p.strip <- list(cex=1, lines=2, fontface='bold')
rasterVis::levelplot(r1.grd, 
          zlim = r1.range,
		  margin = FALSE,
		  main = list.title2[[i]],
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
          col.regions = colorRampPalette(c("blue", "white", "red"))(255),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  colorkey=list(space="bottom", 
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) +  
		  latticeExtra::layer(sp.polygons(chn.shp2, col = "black", lwd = 1)) +
		  latticeExtra::layer(sp.points(pts.sp.sig1, col = "black", pch = 20, cex = 0.5))  


dev.off()


myData <- as.vector(as.matrix(r1.grd))
binwidth <- as.numeric((r1.range[2]-r1.range[1])/20)

# create plot    

n_bins <- length(ggplot2:::bin_breaks_width(range(myData, na.rm = TRUE), width = binwidth)$breaks) - 1L
ggplot() + geom_histogram(aes(x = myData), binwidth = binwidth, fill = colorRampPalette(c("blue", "white", "red"))(n_bins)) + 
	aes(y=100*stat(count)/sum(stat(count)))+
	geom_vline(aes(xintercept = mean(myData, na.rm = TRUE)),col='black',size=2, lty = 2)+
	theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("% of Forests Area") + 
  xlab(list.title2[[i]]) +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") 

ggsave(paste0(yourdatafolder,png.title2[[i]],"_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")

}

##########  severity #################
yrs = 1987:2020

meanNDVI.stack = stack()
medianNDVI.stack = stack()

for (i in 1:length(yrs)){

size1 = data.frame(read.csv(paste0("F:/ChinaWork/DisturbMap/GEE_output/patch_severity_ndvi_", yrs[i],".csv")))[,c(2:8)]
size2 = merge(chn.shp@data, size1, by.x = "ID", by.y = "GridID", all.x = TRUE)

tmp = data.frame(long = size2$x, lat = size2$y, mean1 = size2$mean)
tmp2 = data.frame(long = size2$x, lat = size2$y, mean1 = size2$median)


tmp.grd = rasterFromXYZ(tmp)
tmp2.grd = rasterFromXYZ(tmp2)

tmp.grd[tmp.grd == 0] = NA
tmp2.grd[tmp2.grd == 0] = NA

meanNDVI.stack = addLayer(meanNDVI.stack, tmp.grd)
medianNDVI.stack = addLayer(medianNDVI.stack, tmp2.grd)

}

list.stack3 = list(meanNDVI.stack/1000, medianNDVI.stack/1000)
list.title3 = list("Mean NDVI", "Median NDVI")

png.title3 = list("MeanNDVI", "MedianNDVI")

## produce figures using for loop  
for (i in 1:2){
r1.grd = calc(list.stack3[[i]], mean, na.rm = TRUE)
r1.grd[forestarea.grd < 5] = NA

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.05,0.95), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(min_, max_) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

png(paste0(yourdatafolder,png.title3[[i]],".png"),height = 1500, width = 2000, res = 300, units = "px")

p.strip <- list(cex=1, lines=2, fontface='bold')
rasterVis::levelplot(r1.grd, 
          zlim = r1.range,
		  margin = FALSE,
		  main = list.title3[[i]],
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
		  col.regions = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(255),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  colorkey=list(space="right", 
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) +  
		  latticeExtra::layer(sp.polygons(chn.shp2, col = "black", lwd = 2))  


dev.off()

r1.grd = calc(list.stack3[[i]], mean, na.rm = TRUE)
r1.grd[forestarea.grd < 5] = NA

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.05,0.95), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(min_, max_) 
# r1.grd[r1.grd > r1.range[2]] = r1.range[2]
# r1.grd[r1.grd < r1.range[1]] = r1.range[1]
r1.grd[r1.grd > r1.range[2]] = NA
r1.grd[r1.grd < r1.range[1]] = NA

myData <- as.vector(as.matrix(r1.grd))
binwidth <- as.numeric((r1.range[2]-r1.range[1])/20)

# create plot    
library(ggplot2)   # CRAN version 2.2.1 used
n_bins <- length(ggplot2:::bin_breaks_width(range(myData, na.rm = TRUE), width = binwidth)$breaks) - 1L
ggplot() + geom_histogram(aes(x = myData), binwidth = binwidth, fill = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(n_bins)) + 
aes(y=100*stat(count)/sum(stat(count)))+
geom_vline(aes(xintercept = mean(myData, na.rm = TRUE)),col='black',size=2, lty = 2)+
theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("% of Forests Area") + 
  xlab(list.title3[[i]]) +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") 

ggsave(paste0(yourdatafolder,png.title3[[i]],"_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")


}

################ calculate trends: Kendall trend test 
meanNDVI_trend = Trend1(meanNDVI.stack)
medianNDVI_trend = Trend1(medianNDVI.stack)

list.stack4 = list(meanNDVI_trend/1000, medianNDVI_trend/1000)
list.title4 = list("Mean NDVI trend", "Median NDVI trend")

png.title4 = list("Mean_NDVI_trend", "Median_NDVI_trend")

for (i in 1:2){
r1.grd = list.stack4[[i]][[1]]
r1.grd[forestarea.grd < 5] = NA

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.01,0.99), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(-max(abs(min_),abs(max_)),max(abs(min_),abs(max_))) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

pts.sp.sig1 = Ex.pts(list.stack4[[i]][[2]], sig.level = 0.05) #extract significant relation points

png(paste0(yourdatafolder,png.title4[[i]],".png"),height = 1500, width = 2000, res = 300, units = "px")

p.strip <- list(cex=1, lines=2, fontface='bold')
levelplot(r1.grd, 
          zlim = r1.range,
		  margin = FALSE,
		  main = list.title4[[i]],
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
          col.regions = colorRampPalette(c("blue", "white", "red"))(255),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  colorkey=list( space="bottom", 
						# x= 0.5, y = 0.5,
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) +  
		  latticeExtra::layer(sp.polygons(chn.shp2, col = "black", lwd = 1)) +
		  latticeExtra::layer(sp.points(pts.sp.sig1, col = "black", pch = 20, cex = 0.5))  


dev.off()

myData <- as.vector(as.matrix(r1.grd))
binwidth <- as.numeric((r1.range[2]-r1.range[1])/20)

# create plot    

n_bins <- length(ggplot2:::bin_breaks_width(range(myData, na.rm = TRUE), width = binwidth)$breaks) - 1L
ggplot() + geom_histogram(aes(x = myData), binwidth = binwidth, fill = colorRampPalette(c("blue", "white", "red"))(n_bins)) + 
	aes(y=100*stat(count)/sum(stat(count)))+
	geom_vline(aes(xintercept = mean(myData, na.rm = TRUE)),col='black',size=2, lty = 2)+
	theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("% of Forests Area") + 
  xlab(list.title4[[i]]) +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") 

ggsave(paste0(yourdatafolder,png.title4[[i]],"_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")

}

#####################################################
############## plot total disturbed area #################
chn.shp = readOGR(dsn=yourdatafolder, layer="wc2.1_10m_prec.1970-2000")
forestarea.df = read.csv("forest_area_by_GridID2.csv")
forestarea.df$ForestArea[which(is.na(forestarea.df$ForestArea))] = 0

chn.shp@data = merge(chn.shp@data, forestarea.df, by.x = "ID", by.y = "GridID")
forestarea.grd = data.frame(chn.shp@data$x, chn.shp@data$y,chn.shp@data$ForestArea)
forestarea.grd = rasterFromXYZ(forestarea.grd)

chn.shp2 = readOGR(dsn=yourdatafolder, layer="CHN_adm1_wTaiwan")

## lt
lt =  data.frame(read.csv(paste0("istubrance_area_by_GridID.csv")))
lt = merge(chn.shp@data, lt, by.x = "ID", by.y = "GridID", all.x = TRUE)

yrs = 1986:2020

sum.stack = stack()
burnrate.stack = stack()

for (i in 8:41){

tmp = data.frame(long = lt$x, lat = lt$y, mean1 = (lt[,i]/1000000)/lt$ForestArea)
tmp.grd = rasterFromXYZ(tmp)
tmp.grd[tmp.grd == 0] = NA

tmp2 = data.frame(long = lt$x, lat = lt$y, mean1 = lt[,i]/1000000)
tmp2.grd = rasterFromXYZ(tmp2)
tmp2.grd[tmp2.grd == 0] = NA

sum.stack = addLayer(sum.stack, tmp2.grd)
burnrate.stack = addLayer(burnrate.stack, tmp.grd)

}


list.stack8 = list(sum.stack, 100*burnrate.stack)
list.title8 = list( "Annually Disturbed Area (km2 yr-1)", "Percent Forest Disturbed Annually (% yr-1)")

png.title8 = list("TotalDisturbedSize_final", "PercentDisturbed_final")

for (i in 1:2){

r1.grd = calc(list.stack8[[i]], mean, na.rm = TRUE)

r1.grd[forestarea.grd < 5] = NA

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.001,0.999), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(min_, max_) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

png(paste0(yourdatafolder,png.title8[[i]],".png"),height = 1500, width = 2000, res = 300, units = "px")

p.strip <- list(cex=1, lines=2, fontface='bold')
rasterVis::levelplot(r1.grd, 
          zlim = r1.range,
		  margin = FALSE,
		  main = list.title8[[i]],
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
		  col.regions = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(255),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  colorkey=list(space="right", 
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) +  
		  latticeExtra::layer(sp.polygons(chn.shp2, col = "black", lwd = 2))  


dev.off()


myData <- as.vector(as.matrix(r1.grd))
binwidth <- as.numeric((r1.range[2]-r1.range[1])/20)

# create plot    
library(ggplot2)   # CRAN version 2.2.1 used
n_bins <- length(ggplot2:::bin_breaks_width(range(myData, na.rm = TRUE), width = binwidth)$breaks) - 1L
ggplot() + geom_histogram(aes(x = myData), binwidth = binwidth, fill = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(n_bins)) + 
aes(y=100*stat(count)/sum(stat(count)))+
geom_vline(aes(xintercept = mean(myData, na.rm = TRUE)),col='black',size=2, lty = 2)+
theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("% of Forests Area") + 
  xlab(list.title8[[i]]) +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") 

ggsave(paste0(yourdatafolder,png.title8[[i]],"_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")

}

# #############################look at trend ##################################
# install.packages("spatialEco")
library(spatialEco)
#length of NAs
lenNA = function(x){return(length(which(is.na(x))))}

# Write a function to calcuate the trend
Trend1 = function(stkgrd){
require(spatialEco)
stkgrd_lenNA=calc(stkgrd,lenNA)
stkgrd[stkgrd_lenNA > 30] = -9999 #replace NAs with -9999 to use the raster.kendall function
stkgrd_trend= raster.kendall(stkgrd,p.value = TRUE, na.rm = TRUE)
stkgrd_trend[stkgrd_lenNA > 10] = NA
return(stkgrd_trend)
}

# calculate trends: Kendall trend test 
sum_trend = Trend1(sum.stack)
burnrate_trend = Trend1(burnrate.stack)

# plot

list.stack9 = list(sum_trend, 100*burnrate_trend)
list.title9 = list("Trend for Annually Disturbed Area (km2 yr-2)", "Trend for Percent Forest Disturbed Annually (% yr-2)")

png.title9 = list("TotalDisturbedSize_trend_final", "PercentDisturbed_trend_final")

for (i in 1:2){

r1.grd = list.stack9[[i]][[1]]
r1.grd[forestarea.grd < 5] = NA

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.001,0.999), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(-max(abs(min_),abs(max_)),max(abs(min_),abs(max_))) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

plot(r1.grd)

## extract p value < 0.05, and add to plot
pts.sp.sig1 = Ex.pts(list.stack9[[i]][[2]], sig.level = 0.05) #extract significant relation points


png(paste0(yourdatafolder,png.title9[[i]],".png"),height = 1500, width = 2000, res = 300, units = "px")

p.strip <- list(cex=1, lines=2, fontface='bold')
rasterVis::levelplot(r1.grd, 
          zlim = r1.range,
		  margin = FALSE,
		  main = list.title9[[i]],
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
          col.regions = colorRampPalette(c("blue", "white", "red"))(255),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  colorkey=list( space="bottom", 
						# x= 0.5, y = 0.5,
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) +  
		  latticeExtra::layer(sp.polygons(chn.shp2, col = "black", lwd = 1)) +
		  latticeExtra::layer(sp.points(pts.sp.sig1, col = "black", pch = 20, cex = 0.5))  


dev.off()



myData <- as.vector(as.matrix(r1.grd))
binwidth <- as.numeric((r1.range[2]-r1.range[1])/20)

# create plot    

n_bins <- length(ggplot2:::bin_breaks_width(range(myData, na.rm = TRUE), width = binwidth)$breaks) - 1L
ggplot() + geom_histogram(aes(x = myData), binwidth = binwidth, fill = colorRampPalette(c("blue", "white", "red"))(n_bins)) + 
	aes(y=100*stat(count)/sum(stat(count)))+
	geom_vline(aes(xintercept = mean(myData, na.rm = TRUE)),col='black',size=2, lty = 2)+
	theme_bw() + 
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("% of Forests Area") + 
  xlab(list.title9[[i]]) +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.position = "none") 

ggsave(paste0(yourdatafolder,png.title9[[i]],"_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")

}

