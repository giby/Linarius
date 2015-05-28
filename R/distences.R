#############################################################################################################################################################################
####                                                                   Distence Index for populations                                                                    ####
#############################################################################################################################################################################

#Reynolds Index 
Reynolds.dist<-function(xx,yy,plox,ploy){
	AF1<-allele.frec(xx,plox)
	AF2<-allele.frec(yy,ploy)
	AF1<-(AF1[,ncol(AF1)])/100
	AF2<-(AF2[,ncol(AF2)])/100
	J1<-mean(AF1^2)
	J2<-mean(AF2^2)
	J12<-mean(AF1*AF2)
	#calcaulation= 

sqrt(sum((AF1-AF2)^2)/(2*sum(1-(AF1*AF2))))->Rey.dist 

	return(Rey.dist)
	}
	
	Roger.dist<-function(xx,yy,plox,ploy){
	AF1<-allele.frec(xx,plox)
	AF2<-allele.frec(yy,ploy)
	AF1<-(AF1[,ncol(AF1)])/100
	AF2<-(AF2[,ncol(AF2)])/100 
	J1<-mean(AF1^2)
	J2<-mean(AF2^2)
	J12<-mean(AF1*AF2)
	#calcaulation= 

sqrt(sum((AF1-AF2)^2)/2)/length(AF1)->Rog.dist 

	return(Rog.dist)
	}

Euclide.dist<-function(xx,yy,plox,ploy){
	AF1<-allele.frec(xx,plox)
	AF2<-allele.frec(yy,ploy)
	AF1<-(AF1[,ncol(AF1)])/100
	AF2<-(AF2[,ncol(AF2)])/100
	J1<-mean(AF1^2)
	J2<-mean(AF2^2)
	J12<-mean(AF1*AF2)
	#calcaulation= 

sqrt(sum((AF1-AF2)^2))->Euc.dist 

	return(Euc.dist)
	}

Nei.dist<-function(xx,yy,plox,ploy){
	AF1<-allele.frec(xx,plox)
	AF2<-allele.frec(yy,ploy)
	AF1<-(AF1[,ncol(AF1)])/100
	AF2<-(AF2[,ncol(AF2)])/100
	J1<-mean(AF1^2)
	J2<-mean(AF2^2)
	J12<-mean(AF1*AF2)
	#calcaulation= 
	-log(J12/sqrt(J1*J2),base = exp(1)  )->NDist
	return(NDist)
	
	
}

#############################################################################################################################################################################
####                                                                         Physical Distence                                                                           ####
#############################################################################################################################################################################
### Lambert Dist
	dist.gps.Lambert<-
function(lat1,lon1,lat2,lon2)
{
Req<-6378.1370
Rpo<-6356.7523142

	f<-(Req-Rpo)/Req
	ee<-2*f-f^2
	rpl1<-atan(sqrt(1-ee)*tan(lat1))
	rpl2<-atan(sqrt(1-ee)*tan(lat2))
	P<-mean(c(rpl1, rpl2))
	Q<-(rpl2- rpl1)/2
	CA<-(acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos((lon1-lon2)*pi/180)))
	X<-((CA-sin(CA))*sin(P)*cos(Q)*sin(P)*cos(Q)/(cos(CA/2)*cos(CA/2)))
	Y<-((CA-sin(CA))*sin(Q)*cos(P)*sin(Q)*cos(P)/(sin(CA/2)*sin(CA/2)))
	dphy<-Req*(CA-(f*(X+Y)/2))
return(dphy)
}
###Spherical
dist.gps.cst<-
function(lat1,lon1,lat2,lon2)
{
dphy<-(acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos((lon1-lon2)*pi/180))*6371)
return(dphy)
}
###Potato
dist.gps.var<-
function(lat1,lon1,lat2,lon2)
{
	a<-6378.1370
	b<-6356.7523142
	if (lat1*lat2 >= 0) 
	Rad<-mean(
c((abs(lat1)*b+(90-abs(lat1))*a+(90-abs(lat1))*a)/(180-abs(lat1)),
(abs(lat2)*b+(90-abs(lat2))*a+(90-abs(lat2))*a)/(180-abs(lat2)),
max(c((abs(lat1)*b+(90-abs(lat1))*a+(90-abs(lat1))*a)/(180-abs(lat1)),
(abs(lat2)*b+(90-abs(lat2))*a+(90-abs(lat2))*a)/(180-abs(lat2)))))
)
	if (lat1*lat2 <= 0) 
	Rad<-mean(c((abs(lat1/3)*b+(90-abs(lat1/3))*a)/90,(abs(lat2/3)*b+(90-abs(lat2/3))*a)/90))
dphy<-(acos(sin(lat1*pi/180)*sin(lat2*pi/180)+cos(lat1*pi/180)*cos(lat2*pi/180)*cos((lon1-lon2)*pi/180))*Rad)
return(dphy)
}

