#############################################################################################################################################################################
####                                                            			         Anamorphisis	                                                                     ####
#############################################################################################################################################################################

#Basic
expend.coord<-function(xx,focal=c(64.13806,-21.92861),dist,Dist.Ref=0.8) {
	cbind((xx[,1]-focal[1]),(xx[,2]-focal[2]))->centre
(centre*dist/Dist.Ref)->redu
cbind(rownames(xx),xx,(redu[,1]+focal[1]),(redu[,2]+focal[2]))->NeoCoord
colnames(NeoCoord)<-c("ID","lat","lon","lat.ana","long.ana")
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

boot.reshape.coord<-function(xx,dist,nboot=10) {
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
return(lati)

}

#cbind(uniqueGPS(datgps),DIST)[sample(nrow(cbind(uniqueGPS(datgps),DIST))),]

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

harmonize.coord<-function(xx,yy) {
	focal<-colMeans(xx)
	cbind((xx[,1]-focal[1]),(xx[,2]-focal[2]))->centre
(centre*dist.coord(yy)/dist.coord(xx))->redu 
focal<-colMeans(yy)
cbind((redu[,1]+focal[1]),(redu[,2]+focal[2]))->NeoCoord
#colnames(NeoCoord)<-c("ID","lat","lon","lat.ana","long.ana")
return(NeoCoord)

}


#############################################################################################################################################################################
####                                                            		         Tripolar Heat map	                                                                     ####
#############################################################################################################################################################################


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
