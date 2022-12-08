
## Fig 5

library(sf)
library(tcltk)
library(raster)
library(rasterVis)
library(ggplot2)
library(rgdal)

#####################################################################################################
## compare disturbance before 2000 and after 2000 [natural forest projection project]
###################################################################################################
nfpp.shp = readOGR(dsn="F:/ChinaWork/DisturbMap/NFPP", layer="nfpp_bound")
chn0.shp = readOGR(dsn="F:/ChinaWork/Jiajia_Su090420/00_boundary", layer="CHN_adm0_wTaiwan")

## lt
lt =  data.frame(read.csv(paste0("F:/ChinaWork/DisturbMap/GEE_output/distubrance_area_by_GridID.csv")))
lt = merge(chn.shp@data, lt, by.x = "ID", by.y = "GridID", all.x = TRUE)

lt.stack = stack()
for (i in 8:ncol(lt)){

tmp = data.frame(long = lt$x, lat = lt$y, mean1 = (lt[,i]/1000000)/lt$ForestArea)
tmp.grd = rasterFromXYZ(tmp)

tmp.grd[tmp.grd == 0] = NA

Q1 =  quantile(as.matrix(tmp.grd),probs = c(0.01,0.99), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

tmp.grd[tmp.grd > max_ | tmp.grd < min_] = NA

lt.stack = addLayer(lt.stack, tmp.grd)
}

lt.stack[is.na(lt.stack)] = 0
names(lt.stack) <- names(lt)[c(8:ncol(lt))]
# 
r1 = 100*calc(lt.stack[[c(1:15)]], mean, na.rm = TRUE) # 1986-2000
r2 =  100*calc(lt.stack[[c(16:34)]], mean, na.rm = TRUE)

r1.grd = stack(r1, r2)
names(r1.grd) <- c("Before NFPP", "After NFPP")
# smooth here
for(i in 1:2){
# r0 = focal(r1.grd[[i]], w=matrix(1/9,nrow=3,ncol=3), na.rm = TRUE)
# r0.df = data.frame(long = chn.shp@data$x, lat = chn.shp@data$y, dat = extract(r0, data.frame(long = chn.shp@data$x, lat = chn.shp@data$y)))
r0.df = data.frame(long = chn.shp@data$x, lat = chn.shp@data$y, dat = extract(r1.grd[[i]], data.frame(long = chn.shp@data$x, lat = chn.shp@data$y)))
r0 = rasterFromXYZ(r0.df)
r1.grd[[i]] <- r0

}
# smooth end here

r1.grd[forestarea.grd < 5] = NA

Q1 =  quantile(as.matrix(r1.grd),probs = c(0.001,0.999), na.rm = TRUE)
min_ = Q1[1]
max_ = Q1[2]

r1.range = c(min_, max_) 
r1.grd[r1.grd > r1.range[2]] = r1.range[2]
r1.grd[r1.grd < r1.range[1]] = r1.range[1]

names(r1.grd) <- c("Before NFPP", "After NFPP")

plot(r1.grd)
# plot(chn.shp2, add = TRUE)

png(paste0("F:/ChinaWork/DisturbMap/Results/","beforeandafter2000",".png"),height = 2500, width = 2000, res = 300, units = "px")

p.strip <- list(cex=1, lines=2, fontface='bold')
rasterVis::levelplot(r1.grd, 
          zlim = r1.range,
		  # main = list.title[[i]],
		  # main = expression("Net C Uptake Trend"~ "("~gC ~ m^{-2} ~ yr ^{-1}~ yr ^{-1}~ "1980-2017" ~")"),		  
          maxpixels = nrow(r1.grd)*ncol(r1.grd),
          # col.regions = rev(colorRampPalette(c("blue", "white", "red"))(255)),
		  col.regions = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(255),
          # scales=list(x=list(cex=1),y=list(cex=1)),
          # xlab=list(label = "Longtitude", cex=1),ylab=list(label = "Latitude", cex=1),
		  xlab=NULL, ylab=NULL,
		  scales=list(draw=FALSE),
		  layout=c(1, 2),
		  colorkey=list(space="right", 
		                at=seq(r1.range[1], r1.range[2], length.out=100), 
						height=1, 
						labels = list(cex = 1.5)),
          par.strip.text=p.strip) +  
		  latticeExtra::layer(sp.polygons(chn0.shp, col = "black", lwd = 2)) +
		  latticeExtra::layer(sp.polygons(nfpp.shp, col = "blue", lwd = 2)) 


dev.off()


######## plot histgram before and after NFPP
# http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
library(plyr)

d1 = as.vector(as.matrix(r1.grd[[1]]))
d1 = d1[which(!is.na(d1))]
d2 = as.vector(as.matrix(r1.grd[[2]]))
d2 = d2[which(!is.na(d2))]

df <- data.frame(
  Period = factor(c(rep("BeforeNFPP", length(d1)), rep("AfterNFPP", length(d2)))),
  disturb = c(d1,d2)
  )

mu <- ddply(df, "Period", summarise, grp.mean=median(disturb))
sd <- ddply(df, "Period", summarise, grp.mean=sd(disturb))

p <- ggplot(df, aes(x=disturb, fill=Period, color=Period)) +
  geom_histogram(position="identity", alpha=0.5) + 
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Period),
             linetype="dashed", size = 2)


 p + theme_bw() + 
  theme(legend.position = c(0.75, 0.75))+
  # scale_fill_discrete(breaks=shortnames,
                      # name="",
                      # labels=shortnames)+
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("") + 
  xlab("% Forest Disturbed") +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.title = element_text(colour="black", size=18, face="bold"))
 
# ggsave(paste0("F:/ChinaWork/DisturbMap/Results/","NFPP_national","_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")
 
 
### only look at the NFPP regions
nfpp.grd = rasterize(nfpp.shp,  r1.grd) 
nfpp.grd = nfpp.grd >= 1
 
r0.grd = r1.grd*nfpp.grd
 
d1 = as.vector(as.matrix(r0.grd[[1]]))
d1 = d1[which(!is.na(d1))]
d2 = as.vector(as.matrix(r0.grd[[2]]))
d2 = d2[which(!is.na(d2))]

df <- data.frame(
  Period = factor(c(rep("BeforeNFPP", length(d1)), rep("AfterNFPP", length(d2)))),
  disturb = c(d1,d2)
  )

mu <- ddply(df, "Period", summarise, grp.mean=median(disturb))
sd <- ddply(df, "Period", summarise, grp.mean=sd(disturb))

p <- ggplot(df, aes(x=disturb, fill=Period, color=Period)) +
  geom_histogram(position="identity", alpha=0.5) + 
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Period),
             linetype="dashed", size = 2)


 p + theme_bw() + 
  theme(legend.position = c(0.75, 0.75))+
  # scale_fill_discrete(breaks=shortnames,
                      # name="",
                      # labels=shortnames)+
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("") + 
  xlab("% Forest Disturbed") +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.title = element_text(colour="black", size=18, face="bold"))
  
# ggsave(paste0("F:/ChinaWork/DisturbMap/Results/","NFPP_region","_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")
 
### outside NFPP regions
r0.grd = r1.grd
r0.grd[nfpp.grd == 1] = NA

d1 = as.vector(as.matrix(r0.grd[[1]]))
d1 = d1[which(!is.na(d1))]
d2 = as.vector(as.matrix(r0.grd[[2]]))
d2 = d2[which(!is.na(d2))]

df <- data.frame(
  Period = factor(c(rep("BeforeNFPP", length(d1)), rep("AfterNFPP", length(d2)))),
  disturb = c(d1,d2)
  )

mu <- ddply(df, "Period", summarise, grp.mean=median(disturb))
sd <- ddply(df, "Period", summarise, grp.mean=sd(disturb))


p <- ggplot(df, aes(x=disturb, fill=Period, color=Period)) +
  geom_histogram(position="identity", alpha=0.5) + 
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Period),
             linetype="dashed", size = 2)


 p + theme_bw() + 
  theme(legend.position = c(0.75, 0.75))+
  # scale_fill_discrete(breaks=shortnames,
                      # name="",
                      # labels=shortnames)+
  theme(legend.text = element_text(size = 18))+
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),axis.text.x  = element_text(colour="black",size=18))+
  theme(axis.title.y = element_text(face="bold", colour="black", size=18),axis.text.y  = element_text(colour="black",size=18))+
  ylab("") + 
  xlab("% Forest Disturbed") +
  theme(strip.text.x = element_text(size=18), strip.text.y = element_text(size=18)) + 
  theme(legend.title = element_text(colour="black", size=18, face="bold"))
  
 # ggsave(paste0("F:/ChinaWork/DisturbMap/Results/","NFPP_outside","_histgram.png"), width = 6, height = 4, units = "in", bg = "transparent")
