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
