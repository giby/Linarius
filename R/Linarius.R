#############################################################################################################################################################################
###                                                                            Linarius.R                                                                                ####
#############################################################################################################################################################################

  .onAttach <- function(...) {
packageStartupMessage("This is NOT a free software, not reading the licence is a violation of the licence, by continuing using it, you are considered aware of the terms of Licence.l2k and accepting them") 
}


#Alleles frequency & heterozygocy 
#' @title  Function to count allele presence and absence of dominant markers according to ploidy levels
#'
#' @description This function aims at counting allele presence and absence in a ballance population, this considering the ploidy level.
#' 
#' @param xx A binary datafram of genotypes, individuals as row and alleles as column 
#' @param ploidy Interger or vector, ploidy level of the population or individuals 
#' 
#'@examples 
#' data(Birch)
#' allele.count(Betula,code$Ploidy)
#' 

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
#' @title  Calculation of Allele hetherozygocy with Dominant markers and mixed ploidy levels
#'
#' @description This function aims at evaluating allele hetherozygocy in a ballance population, this considering the ploidy level.  Then It calculated an weighted avarage
#' 
#' @param xx A binary datafram of genotypes, individuals as row and alleles as column 
#' @param ploidy Interger or vector, ploidy level of the population or individuals 
#' 
#' @note There is a little biais on small sized sample.
#' 
#' @examples 
#' data(Birch)
#' allele.frec(Betula,code$Ploidy)
#' 
allele.frec<-function(xx,ploidy=2) {
plolev<-sort(unique(ploidy))
NULL->alco
NULL->noms
NULL->glob
for(i in plolev) cbind(alco,(1-(apply(as.matrix(xx[ploidy==i,])==0,2,sum)/(apply(as.matrix(xx[ploidy==i,])==0,2,sum)+apply(as.matrix(xx[ploidy==i,])==1,2,sum)))^(1/i))*100)->alco
#for(i in plolev) cbind(alco, colMeans(xx[ploidy == i, ]^(1/i)*100) -> alco
for(j in plolev) noms<-c(noms,paste(j,"x-frequence(%)", sep="",collapse=""))
colnames(alco)<-noms
for(k in plolev) glob<-c(glob,nrow(xx[ploidy==k,]))
cbind(alco,apply(alco*glob/sum(glob),1,sum))->alco
colnames(alco)<-c(noms,"Frequency")
return(alco)
}
#############################################################################################################################################################################
#sep ploidy al frec 
#' @title  Calculation of Allele frequency with Dominant markers and mixed ploidy levels
#'
#' @description This function aims at evaluating allele frequency in a ballance population, this considering the ploidy level. Then It calculated an weighted avarage
#' 
#' @param xx A binary datafram of genotypes, individuals as row and alleles as column 
#' @param ploidy Interger or vector, ploidy level of the population or individuals 
#' @note There is a little biais on small sized sample.
#' 
#'@examples 
#' data(Birch)
#' allele.hetero(Betula,code$Ploidy)
#' 
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
#' @title Generate a fake dataset fitting allele frequencies and ploidies. 
#'
#' @description 
#' This function generates a random binary datafram with in row: simulated individuals, and in colomn: alleles. 
#' The aim of this function is to test hypothesis, and shall not be used for publishing imaginary results (See licence for more information). 
#' @param frec Allele frequencies: a vector providing the allele frequency of all alleles, every single frequency have to be between 0 and 1
#' @param ploidy Ploidy level of every individuals: a vector providing ploidy level of every single individual. Normally integers. 
#' 
#' @note 
#' Whenever you make something out of this function you have to precise that there are generated data. 
#' You shall not use this for getting fake experimental data (it constitues a violation of this package licence ).
#' Whenever you use that function for bad purposes I will find you, get you and kick your ass! 

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
#############################################################################################################################################################################
#Identity of population 
BjIC<-function(XX,factor,ploidy=2) { 
	allele.frec(XX,ploidy=ploidy)->frec_all
	frec_all<-(frec_all[,ncol(frec_all)])
	levels<-sort(unique(factor))
	0->BIC_VALUE
	for(i in levels){
		
		allele.frec(XX[factor==i,],ploidy=ploidy[factor==i])->BT
		BT<-BT[,ncol(BT)]
BIC_VALUE<-BIC_VALUE+(BT-frec_all)^2
	} 
	BIC_VALUE<-sum(BIC_VALUE)
return(BIC_VALUE)	
	
	}
	#############################################################################################################################################################################
#boot of BIC
Boot.BjIC<-function(XX,factor,ploidy,nboot=100) {
	boot<-vector(mode = "numeric", length = nboot)
	len<-(1:nboot)
	allele.frec(XX,ploidy=ploidy)->frec_all
	frec_all<-(frec_all[,ncol(frec_all)])
	for(i in len) {
		datagen(frec=frec_all/100,ploidy=ploidy)->fake
		colnames(fake)<-names(XX)
		BIC(fake,factor=factor,ploidy=ploidy)->boot[i]
		#Bet <- as.genclone(df2genind(fake, ind.names = row.names(gelplo), type = "PA"))
#sethierarchy(Bet) <- codeplo
#setpop(Bet) <- ~Ploidy
#res <- poppr.amova(Bet, hier = ~Ploidy)
#res[1:4]$componentsofcovariance[1,2]->boot[i]	
}
return(boot)
}

	#############################################################################################################################################################################
#' @title Coloring terminal branch of trees
#'   
#'
#' @description Generate colors for terminal edge branch of a phylogenic tree
#' 
#' @param phy A tree of class phylo
#' @param groups A vector indicating what should be of the same color
#' @param color A vector indicating the color to choose
#' @param root.color A single color, to use for non terminal edge 
#' @note This function should be used with the package ade4 
#' @seealso ade4
#' @examples 
#' data(Birch)
#' require(vegan)
#' require(ade4)
#' require(ape)
#' dist=vegdist(Betula, method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE)
#' dendro=hclust(dist,"complete")
#' dendro2<-as.phylo(dendro,cex=0.5)
#' plot(dendro2,,type="u",lab4ut="axial",font=1, adj=1,label.offset = 0.01,edge.color = col.edge.phylo(dendro2,code$Ploidy),tip.color = code$Ploidy)




col.edge.phylo<-function (phy,groups,color=c(2,3,4,5),root.color=1)
{
	#col.edge<-rep(root.color,length(phy$tip))
	col.edge<-rep(root.color,length(phy$edge[,1]))
	cbind(phy$edge,col.edge)->Edge
	cbind(Edge,(1:length(phy$edge[,1])))-> Edge
	(Edge[order(as.integer(Edge[,2]), decreasing=F),])->Edge
	groups->Edge[1:length(phy$tip),3]
	(Edge[order(as.integer(Edge[,4]), decreasing=F),])->Edge
	return(Edge[,3])
}
	
	#############################################################################################################################################################################
#' @title Compute genetic distances between populations 
#'   
#'
#' @description Returns a dist object with genetic distences between populations. 4 indexes are availlable: "Euclidean", "Reynolds", "Roger" and "Nei".
#' 
#' @param xx A binary datafram of genotypes, individuals as row and alleles as column 
#' @param pop A vector informing of population every sample belongs to 
#' @param method A method for computing genetic distence between populations, can be any unambigous abreviation of "euclidean", "reynolds", "roger" or "nei"
#' @param ploidy Interger or vector, ploidy level of the population or individuals 
#' @examples 
#' data(Birch)
#' dist.pop(Betula, code$Pop,ploidy=code$Ploidy,method="rey") 

	
	dist.pop<-function(xx,pop, method="reynolds" ,ploidy=2)
{
    if (length(ploidy) == 1)
    ploidy <- rep(ploidy, length(rownames(xx)))
	if (length(rownames(xx)) != length(ploidy)) 
        stop("Mismatch of number of samples and number of ploidy ID")
	if (!is.na(pmatch(method, "euclidian"))) 	
	if (length(rownames(xx)) != length(pop)) 
        stop("Mismatch of number of samples and number of population ID")
	if (!is.na(pmatch(method, "euclidian"))) 
        method <- "euclidean"
    METHODS <- c("euclidean", "reynolds", "roger", "nei")
    method <- pmatch(method, METHODS)
    if (is.na(method)) 
        stop("invalid distance method")
    if (method == -1) 
        stop("ambiguous distance method")
        FUN<-c(Euclide.dist,Reynolds.dist,Roger.dist,Nei.dist)[[method]]
if (method >= 5) 
      stop("Huston, we've got a problem")   
	pop.lev<-unique(pop)
  dMat <- matrix(0,ncol=length(pop.lev),nrow=length(pop.lev))
  colnames(dMat) <- pop.lev -> rownames(dMat)
  
  for(i in 1:(length(pop.lev)-1))
  {
    for(j in (i+1):length(pop.lev))
    {
      a <- FUN(xx[pop==pop.lev[i],],xx[pop==pop.lev[j],], ploidy[pop==pop.lev[i]], ploidy[pop==pop.lev[j]])
      dMat[i,j] <- a
      dMat[j,i] <- a
     # return(ploidy[pop==pop.lev[j]])
    }
  }
  dMat<-as.dist(dMat)
return(dMat)
}

	#############################################################################################################################################################################
#' @title Compute distance between GPS points 
#'   
#'
#' @description Returns a dist object with  distences between location in km 
#' 
#' @param lon A vector with longitudes 
#' @param lat A vector with latitudes
#' @param names A vector with names of position
#' @param model A model for computing distence between places, can be any unambigous abreviation of "lambert", "spherical" or "potato"
#' @examples 
#' data(Birch)
#' distGPS(code$Lat,code$Long,rownames(code))
distGPS <-function(lat,lon,names,model="lambert")
{
		if (length(lat) != length(lon)) 
        stop("Latitudes and Longitudes difer in size")
        if (length(lat) != length(names)) 
        stop("Incorect number of names")
 #Methode
     METHODS <- c("lambert", "spherical", "potato")
    model <- pmatch(model, METHODS)
    if (is.na(model)) 
        stop("invalid distance method")
    if (model == -1) 
        stop("ambiguous distance method")
FUN<-c(dist.gps.Lambert,dist.gps.cst,dist.gps.var)[[model]]
if (model >= 4) 
      stop("Huston, we've got a problem") 
#Computation             
	  dMat <- matrix(0,ncol=length(lon),nrow=length(lon))
  colnames(dMat) <- names
  rownames(dMat) <- names
  for(i in 1:(length(lon)-1))
  {
    for(j in (i+1):length(lon))
    {
      a <- FUN(lat[i],lon[i],lat[j],lon[j])
      if((lat[i] == lat[j])&&(lon[i]==lon[j])) a<-0
      dMat[i,j] <- a
      dMat[j,i] <- a
    }
  }
dMat<-as.dist(dMat)  
return(dMat)
}

