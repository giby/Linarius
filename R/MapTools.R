#############################################################################################################################################################################
####                                                            			         Anamorphisis	                                                                     ####
#############################################################################################################################################################################

#Basic
#' @title  Adapt distence of GPS coordinates according to distences to a single points
#'
#' @description 
#' This function aims to prepare anamorphosis maps or cartograms.
#' This function adapts GPS coordinates according to a focal point. The final Coordinates reflect the distence to the focal point but the distence between two image points is meaningless. 
#' 
#' @param xx Coordinates of the points to transform 
#' @param focal Focal points of the anamorphosis
#' @param dist A vector 
#' @param Dist.ref Murphy's Constant
#'
#' 
#' @examples
#' data(Leymus)
#' dist.pop(Leymus,Leypop, ploidy=8)->dist.ley
#' expend.coord(locations[2:14,],locations[1,],dist.ley[1:13],dist.coord(locations[2:14,]))
#' 
#' 

expend.coord<-function(xx,focal=c(64.13806,-21.92861),dist,Dist.Ref=0.8) {
	if(class(xx) != "data.frame") xx<-as.data.frame(xx)
	cbind((xx[,1]-focal[1]),(xx[,2]-focal[2]))->centre
(centre*dist/Dist.Ref)->redu
cbind(rownames(xx),xx,(redu[,1]+focal[1]),(redu[,2]+focal[2]))->NeoCoord
colnames(NeoCoord)<-c("ID","lon","lat","long.ana","lat.ana")
return(NeoCoord)

}

#multiple
reshape.coord<-function(xx,dist) {
	cbind(rownames(xx),xx,xx)->yy
	for( i in 1:nrow(xx)){
	#expend.coord(yy[,c(4,5)],focal=as.numeric(yy[i,c(4,5)]),dist[i,],Dist.Ref=nrow(xx)*mean(dist[i,])/(nrow(xx)-1))->yy	
	expend.coord(yy[,c(4,5)],focal=as.numeric(yy[i,c(4,5)]),dist[i,],Dist.Ref=mean(dist[i,]))->yy
	}

cbind(rownames(xx),xx,yy[,c(4,5)])->NeoCoord
colnames(NeoCoord)<-c("ID","lat","lon","lat.ana","long.ana")
return(NeoCoord)

}

#Boot 
#' @title  Generate non-focal anamorphosis 
#'
#' @description 
#' This function aims to prepare anamorphosis maps or cartograms.
#' 
#' @param xx A data frame contening the coordinates of the points to transform. First colomn is longitude, second is latitude. 
#' @param dist A dist object with the distence to outline by anamorphosis
#' @param nboot Times you need for your coffee break in seconds 
#' @examples
#' data(Leymus)
#' dist.pop(Leymus,Leypop, ploidy=8)->dist.ley
#' boot.reshape.coord(locations,dist.ley,nboot=100)
#' 
#'
#' 
#' 
#' 


boot.reshape.coord<-function(xx,dist,nboot=10) {
	        if(class(xx) != "data.frame") xx<-as.data.frame(xx)
			if(class(dist) != "matrix") dist <-as.matrix(dist)
			cbind(xx,dist)->yy
			
			NULL->BOOT.lat# <-matrix(data = NA, nrow = nrow(xx), ncol = nboot)
			NULL->BOOT.lon# <-matrix(data = NA, nrow = nrow(xx), ncol = nboot)

	for( i in 1:nboot){
zz<-yy[sample(nrow(yy)),]
(reshape.coord(zz[,c(1,2)],as.matrix(zz[,3:ncol(zz)])))->zzz
#rownames(zz)->rownames(zzz)
BOOT.lat<-(cbind(BOOT.lat,(zzz[with(zzz,order(rownames(zzz))),])[,4]))
BOOT.lon<-(cbind(BOOT.lon,(zzz[with(zzz,order(rownames(zzz))),])[,5]))
	}
#rownames(xx)->rownames(BOOT.lon)
#cbind(rownames(xx),xx,yy[,c(4,5)])->NeoCoord
#colnames(NeoCoord)<-c("ID","lat","lon","lat.ana","long.ana")

cbind(rowMeans(BOOT.lat),rowMeans(BOOT.lon))->lati
rownames(lati)<-rownames(xx)
#return(lati)
cbind(rownames(xx),xx,lati)->NeoCoord
colnames(NeoCoord)<-c("ID","lon","lat","long.ana","lat.ana")
return(NeoCoord)

}

#cbind(uniqueGPS(datgps),DIST)[sample(nrow(cbind(uniqueGPS(datgps),DIST))),]
#' @title  Computes Murphy's constant 
#'
#' @description Computes Murphy's constant 
#' 
#' @param xx A datafram of coordinates 
#'
#' 
#'
#' 
#' 
#' 

dist.coord<-function(xx) {
	focal<-colMeans(xx)
	cbind((xx[,1]-focal[1]),(xx[,2]-focal[2]))->centre
#(centre*dist/Dist.Ref)->redu
#cbind(rownames(xx),xx,(redu[,1]+focal[1]),(redu[,2]+focal[2]))->NeoCoord
#colnames(NeoCoord)<-c("ID","lat","lon","lat.ana","long.ana")
#return(NeoCoord)
return(mean(sqrt((centre[,1]^2)+(centre[,2]^2))))
}

#FullTrap
Full.reshape.coord<-function(xx,dist) {
			cbind(xx,dist)->zz
			index<-permutations(nrow(xx),nrow(xx),rownames(xx))
			NULL->BOOT.lat# <-matrix(data = NA, nrow = nrow(xx), ncol = nboot)
			NULL->BOOT.lon# <-matrix(data = NA, nrow = nrow(xx), ncol = nboot)

	for( i in 1:nrow(xx)){
#zz<-cbind(index[i,],yy)
zz<-(zz[with(zz,order(index[i,])),])
#zz<-(zz[,2:ncol(zz)])
(reshape.coord(zz[,c(1,2)],as.matrix(zz[,3:ncol(zz)])))->zzz
rownames(zz)->rownames(zzz)
BOOT.lat<-(cbind(BOOT.lat,(zzz[with(zzz,order(rownames(zzz))),])[,4]))
BOOT.lon<-(cbind(BOOT.lon,(zzz[with(zzz,order(rownames(zzz))),])[,5]))
	}

cbind(rowMeans(BOOT.lat),rowMeans(BOOT.lon))->lati
rownames(lati)<-rownames(xx)
return(lati)

}
#' @title  Harmonize anamorphogram with original map 
#'
#' @description 
#' This function aims to prepare anamorphosis maps or cartograms.
#' This function centers and homogenized the distances and positions of image points according to real points. 
#' @param xx Anamorphosis points coordinate
#' @param yy Real points coordinates
#' @note Use this function if the difference of position and/or size between image and source is too big. This function is only a cosmetic feature usefull if your public can misinterprate the cartogram. 
#' 
#'
#' 
#' 
#' 

harmonize.coord<-function(xx,yy) {
	focal<-colMeans(xx)
	cbind((xx[,1]-focal[1]),(xx[,2]-focal[2]))->centre
(centre*dist.coord(yy)/dist.coord(xx))->redu 
focal<-colMeans(yy)
cbind((redu[,1]+focal[1]),(redu[,2]+focal[2]))->NeoCoord
#colnames(NeoCoord)<-c("ID","lat","lon","lat.ana","")
return(NeoCoord)

}


#############################################################################################################################################################################
####                                                            		         Tripolar Heat map	                                                                     ####
#############################################################################################################################################################################

#' @title Transform genotype in RGB colors 
#'
#' @description 
#' This function aims to create a tri-polar heatmap. 
#' A tripolar heat map aims to indicate relatedness between populations. 
#' 
#' @param XX A binary datafram of genotypes, individuals as row and alleles as column
#' @param pop A vector informing of population every sample belongs to 
#' @param red,green,blue A vector of alleles to put in a single color group
#'
#' @examples
#' require(sp)
#' data(Leymus)
#' RGB.pop(Leymus,pop=Leypop,red=c("GreenExtra.160","GreenExtra.105","GreenExtra.176","Green1.188"),green=c("BlueExtra.127","Blue.add.128","BlueExtra.222","Blue.add.129","BlueExtra.182","BlueExtra.131") ,blue=c("Green1.120","GreenExtra.136","Green1.193","Green1.169","Green1.133"))->col.RGB
#' plot(west)
#' points(locations,col=rgb(col.RGB, max = 255) ,pch=16,cex=2)
#' 
#' 
#' 

RGB.pop<- function(XX,pop,red,green,blue) 
{ 
	NULL->RGB
	#unique(pop)->POP
	#return(POP)
for(i in unique(pop)) RGB<-rbind(RGB, c(sum(XX[pop==i,red])/(length(XX[pop==i,red])*nrow(XX[pop==i,red])),sum(XX[pop==i,green])/(length(XX[pop==i,green])*nrow(XX[pop==i,green])),sum(XX[pop==i,blue])/(length(XX[pop==i,blue])*nrow(XX[pop==i,blue])))) 
unique(pop)->rownames(RGB)
#rebase
for(j in 1:3) RGB[,j]<-RGB[,j]-min(RGB[,j])
for(j in 1:3) RGB[,j]<-(RGB[,j]/max(RGB[,j]))*255
#class(RGB) <- "hexmode" 
#as.character(RGB)->RGB
#paste(RGB[,1],RGB[,2],RGB[,3],sep="")->RGB
#unique(pop)->names(RGB)
colnames(RGB)<-c("R","G","B")
return(RGB)
 
}


#Usage: 
#RGB.pop(LeymusFull,pop=Leydata$Pop,red=c("Green1.187","Green1.151"),green=c("BlueExtra.222","BlueExtra.131"),c("Green1.207","Green1.111" ))
#      Al       Dk       Lj       Ml       Rf      Sa1      Sa2      Sa3       Sb      Sc1      Sc2       Sg       Sj       Sr       Th 
#"3096f4" "d39149" "51c8b2" "a80049" "ff001a" "ffe400" "a8ff49" "ffb632" "a0e785" "d94864" "ffda1a" "5cb6ad" "c18dba" "3491d6" "3fbbd1" 
#      Tl       Vn       Vs 
#"0112fd" "0072fe" "0091ff"

#############################################################################################################################################################################
####                                                            		         Heat map	                                                                             ####
#############################################################################################################################################################################
#' @title  Returns values of color of points for heat-map in RGB system 
#'
#' @description Makes heat map 
#' 
#' @param XX something
#' @param MAX I guess it was maximum 
#'
#' 
#'
#' @examples
#' data(Kenya) 
#' heat.values(Kenya$Asconc)
#' 


heat.values<- function(XX,MAX=NULL) #,pop,red,green,blue)
{
 max(XX)->MAX
    matrix(nrow = length(XX), ncol = 3)->RGB
RGB[,1]<-(XX/MAX)
RGB[,2]<-(MAX-XX)/MAX
RGB[,3]<-0

    colnames(RGB)<-c("R","G","B")
    return(RGB)
    
}

#' @title  Makes heat map 
#'
#' @description Plot color points according to an alternate variable.  
#' 
#' @param XX something
#' @param MAX I guess it was maximum 
#' @param loc Position of scalebare 
#' @param length Length of scalebare
#' @param unit Unit of scalebare
#' 
#'
#' @examples
#' require(sp)
#' data(Kenya) 
#' plot(Geotermal.Area)
#' plot(Hells.NP,add=T,col="grey98") 
#' heat.points(Kenya$Asconc, coord ,loc=c(205000,9896000),offset=0.1,length=3000)
#' 


heat.points<- function(XX, coord, MAX=max(XX), loc=c(0,0),length=1000 , unit="" ,pch=16 , cex=par("cex"), digits=0 , offset=0) #,pop,red,green,blue)
{

    matrix(nrow = length(XX), ncol = 3)->RGB
RGB[,1]<-(XX/MAX)
RGB[,2]<-(MAX-XX)/MAX
RGB[,3]<-0

    colnames(RGB)<-c("R","G","B")
   # return(RGB)
    points(coord[,1],coord[,2],col=rgb(RGB[,1], RGB[,2], RGB[,3], max = 1),pch=pch) 
    ### Bar
x<- (0:255*length/255)+loc[1]
y <- c(0,length/(10*3:1))+loc[2]
 cols<-NULL
    cols <- for (i in 0:255) c(cols,rgb((i/255),(255-i)/255,0))

x <- ((length*0:255)/155)+loc[1] 
y <- c(0,length/(2*3:1))+loc[2] 

  
    	for (i in 1:256) rect(x[i],y[1],x[i+1],y[2],col= c(cols,rgb(((i-1)/255),(255-(i-1))/255,0) ),border = NA)
    	rect(x[1],y[1],x[256],y[2])
	 for (i in c(1,128,256)) segments(x[i],y[2],x[i],y[3]) 
	
labels <-round(c(0,MAX/2,MAX ), digits = digits)
text(x[c(1,128,256)],y[3],labels=labels,cex,pos=3, offset=offset) 	
}


#############################################################################################################################################################################
####                                                            		         Scale				                                                                     ####
#############################################################################################################################################################################
#' @title  Draw a scalebar
#'
#' @description Draw a scalebar on a map, 
#' 
#' @param loc Where to put the bar
#' @param lenght how long should it be
#' @param unit The unit to use on for the scalebare. For autocomputed size, several units are possible: kilometer: "km", mile "mi", nautical mile "nmi", meter "m" or league "li"
#' @param degrees A logical argument TRUE is the shape file is in degrees and FALSE is the shapefile is projected 
#' @param scale Ratio between projected map unit and physical distence unit. To be use if degrees = FALSE 
#' 
#' @examples 
#'
#' require(sp)
#' data(Leymus)
#' plot(west)
#' scalebar(c(-20.33,63.3),8)
#' 

scalebar <- function(loc,length,unit="km" ,degrees = TRUE,cex=par("cex"),scale=1, ...) { 
	if(missing(loc)) stop("loc is missing") 
	if(missing(length)) stop("length is missing") 
	 UNIT <- c("km", "mi", "nmi","m","li")
    Unit <- pmatch(unit, UNIT)
    if ((degrees == TRUE)&&(is.na(Unit))) 
        stop("invalid distance unit")
    if ((degrees == TRUE)&&(Unit == -1)) 
        stop("ambiguous distance unit")
	z <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1] 
	if(degrees == FALSE) length<-length/scale
	if((degrees == TRUE)&&(unit=="m"))  length<-length/1000
	if(degrees == TRUE) length<-180*(acos((cos(length/6371)-sin(loc[2]*pi/180)*sin(loc[2]*pi/180))/(cos(loc[2]*pi/180)*cos(loc[2]*pi/180))))/pi
	if((degrees == TRUE)&&(unit=="mi"))  length<-length*1.609347
	if((degrees == TRUE)&&(unit=="nmi"))  length<-length*1.852
	if((degrees == TRUE)&&(unit=="li"))  length<-length*3.248
	x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1] 
	y <- c(0,length/(10*3:1))+loc[2] 
	cols <- rep(c("black","white"),2) 
	for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
	 for (i in 1:5) segments(x[i],y[2],x[i],y[3]) 
	 labels <- z[c(1,3)]-loc[1]
labels <- append(labels,paste(z[5]-loc[1],unit)) 
text(x[c(1,3,5)],y[4],labels=labels,cex,pos=3,offset=0) }

# Test of a more accuarate scale: No visible change - discarted 
#scalebar2 <- function(loc,length,unit="km" ,degrees = TRUE,cex=par("cex"),scale=1, ...) { 
#	if(missing(loc)) stop("loc is missing") 
#	if(missing(length)) stop("length is missing") 
#	z <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1] 
#	if(degrees == FALSE) length<-length/scale
#	if((degrees == TRUE)&&(unit=="m"))  length<-length/1000
#	Rad<-(abs(loc[2])*6356.752+(90-abs(loc[2]))*6378.137)/90
#	if(degrees == TRUE) length<-180*(acos((cos(length/Rad)-sin(loc[2]*pi/180)*sin(loc[2]*pi/180))/(cos(loc[2]*pi/180)*cos(loc[2]*pi/180))))/pi
#	if((degrees == TRUE)&&(unit=="mi"))  length<-length*1.609347
#	if((degrees == TRUE)&&(unit=="nmi"))  length<-length*1.852
#	if((degrees == TRUE)&&(unit=="li"))  length<-length*3.248
#	x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1] 
#	y <- c(0,length/(10*3:1))+loc[2] 
#	cols <- rep(c("black","white"),2) 
#	for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
#	 for (i in 1:5) segments(x[i],y[2],x[i],y[3]) 
#	 labels <- z[c(1,3)]-loc[1]
#labels <- append(labels,paste(z[5]-loc[1],unit)) 
#text(x[c(1,3,5)],y[4],labels=labels,cex,pos=3,offset=0) }

#############################################################################################################################################################################
####                                                            		         Darcy export				                                                             ####
#############################################################################################################################################################################
#' @title  Export anamorphosis coordinate to Darcy
#'
#' @description Creates files to use in the software Darcy for creating anamorphosis maps. Creates 2 files foo-in.txt to use as source and foo-out.txt to use as image. 
#
#' @param xx An object generated by one of the function expend.coord() or boot.reshape.coord()
#' @param file A text to use in filenames  
#' 
Darcy.export<- function(xx,file) 
{
write.table(cbind(xx$lon,xx$lat), file = paste0(file,"-in.txt", sep = ""), append = FALSE, quote = FALSE, sep = " ",    eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = c("#Created with Linarius",""))
write.table(cbind(xx$long.ana,xx$lat.ana), file = paste0(file,"-out.txt", sep = ""), append = FALSE, quote = FALSE, sep = " ",    eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = c("#Created with Linarius",""))
    
}



