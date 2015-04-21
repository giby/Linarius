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
