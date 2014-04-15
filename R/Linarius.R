#############################################################################################################################################################################
###                                                                            Linarius.R                                                                                ####
#############################################################################################################################################################################

#Alleles frequency & heterozygocy 
allele.count<- function(xx,ploidy=2) {
plolev<-sort(unique(ploidy))
NULL->alco
NULL->noms
for(i in plolev) cbind(alco,apply(as.matrix(xx[ploidy==i,])==0,2,sum),apply(as.matrix(xx[ploidy==i,])==1,2,sum))->alco
for(j in plolev) noms<-c(noms,paste(j,"x-absence", sep='',collapse=''),paste(j,"x-presence",sep='', collapse=''))
colnames(alco)<-noms
return(alco)
}
#############################################################################################################################################################################
#sep ploidy al frec
allele.frec<-function(xx,ploidy=2) {
plolev<-sort(unique(ploidy))
NULL->alco
NULL->noms
NULL->glob
#for(i in plolev) cbind(alco,(1-(apply(as.matrix(xx[ploidy==i,])==0,2,sum)/(apply(as.matrix(xx[ploidy==i,])==0,2,sum)+apply(as.matrix(xx[ploidy==i,])==1,2,sum)))^(1/i))*100)->alco
for(i in plolev) cbind(alco, colMeans(xx[ploidy == i, ]^(1/i)*100) -> alco
for(j in plolev) noms<-c(noms,paste(j,"x-frequence(%)", sep="",collapse=""))
colnames(alco)<-noms
for(k in plolev) glob<-c(glob,nrow(xx[ploidy==k,]))
cbind(alco,apply(alco*glob/sum(glob),1,sum))->alco
colnames(alco)<-c(noms,"Frequency")
return(alco)
}
#############################################################################################################################################################################
#sep ploidy al frec 
allele.hetero<-function(xx,ploidy=2) {
plolev<-sort(unique(ploidy))
NULL->alco
NULL->noms
for(i in plolev) cbind(alco,100*(1-((1-(apply(as.matrix(xx[ploidy==i,])==0,2,sum)/(apply(as.matrix(xx[ploidy==i,])==0,2,sum)+apply(as.matrix(xx[ploidy==i,])==1,2,sum)))^(1/i))^i + (apply(as.matrix(xx[ploidy==i,])==0,2,sum)/(apply(as.matrix(xx[ploidy==i,])==0,2,sum)+apply(as.matrix(xx[ploidy==i,])==1,2,sum))))))->alco
for(j in plolev) noms<-c(noms,paste(j,"x-heterozygocy(%)", sep="",collapse=""))
colnames(alco)<-noms
return(alco)
}
#############################################################################################################################################################################
#generate Fake Data
datagen<-function(frec,ploidy) {
fakedata<-matrix(data = NA, nrow = length(ploidy), ncol = length(frec), byrow = FALSE,  dimnames = NULL)
plolev<-unique(ploidy)
lingen<-function(xx){
	matrix(data = NA, nrow = length(ploidy), ncol = 1, byrow = FALSE,  dimnames = NULL)->allfrec
	 for(i in plolev)
	sample(c(0,1),length(ploidy[ploidy==i]),replace=T,prob=c((1-xx)^i,1-((1-xx)^i)))->allfrec[ploidy==i,]
	return(allfrec)
	}
for(j in 1:length(frec)) lingen(frec[j])-> fakedata[,j]
	return(fakedata)
}

#############################################################################################################################################################################
#Counting 
zerodet<- function(X,factor,PV=5,ID=1) {
	#clean NA
X<-X[!is.na(factor),]
factor<-factor[!is.na(factor)]
    # Count
rbind(colSums(X[factor==ID,]),colSums(X[factor!=ID,]),sum(factor==ID)-colSums(X[factor==ID,]),sum(factor!=ID)-colSums(X[factor!=ID,]))->Contingence
rownames(Contingence)<-c("Pres.in","Pres.out","Abs.in","Abs.out")
apply(Contingence,2, min)->limit
Contingence[,limit<=PV]

}