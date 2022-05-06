#=============================================================================#
# Functions adapted from those used in ArrayAnalysis QC modules               #
# modules for quality control and pre-processing of array data                #
#                                                                             #
# Copyright 2010-2016 Dept. of Bioinformatics-BiGCaT, Maastricht University   #
#                                                                             #
# Licensed under the Apache License, Version 2.0 (the "License");             #
# you may not use this file except in compliance with the License.            #
# You may obtain a copy of the License at                                     #
#                                                                             #
# http://www.apache.org/licenses/LICENSE-2.0                                  #
#                                                                             #
# Unless required by applicable law or agreed to in writing, software         #
# distributed under the License is distributed on an "AS IS" BASIS,           #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
# See the License for the specific language governing permissions and         #
# limitations under the License.                                              #
#=============================================================================#


####################
## colorsByFactor ##
####################

#create colors for the plots and the legends
#-------------------------------------------

colorsByFactor <- function(experimentFactor) {
  
  #check whether a factor has been provided
  if(class(experimentFactor)!="factor") stop("Parameter 'experimentFactor' must be of class 'factor'")

  if(length(levels(experimentFactor))==1 | length(levels(experimentFactor))==length(experimentFactor)) {
  #if there is only one group (or no groups are provided) or each sample is in its own group
	#take equally spread colors over the rainbow palette
    plotColors <- rainbow(length(experimentFactor),s=.8,v=.7)
	#set group legend color to white, as there is not a specific group color
	legendColors <- "white"
  } else {
    #compute the number of colors needed for each class
    tab.tmp <- table(experimentFactor)

    #set the two extreme colors for each class
    colors.light <- rainbow(length(levels(experimentFactor)),s=1-sapply(tab.tmp,min,5)*.1)
    colors.dark <- rainbow(length(levels(experimentFactor)),v=1-sapply(tab.tmp,min,5)*.14)

    #create the colors to plot, and colors for the legend (average one per experimental group)
    plotColors <- NULL
    legendColors <- NULL
    for(l in 1:length(levels(experimentFactor))) {
      colorFun <- colorRampPalette(c(colors.light[l],colors.dark[l]))
      tmpColors <- colorFun(tab.tmp[l])
      plotColors[experimentFactor==levels(experimentFactor)[l]] <- tmpColors
      legendColors[l] <- tmpColors[ceiling(length(tmpColors)/2)]
    }
  }
  return(list(plotColors=plotColors,legendColors=legendColors))

}


#####################
## symbolsByFactor ##
#####################

#create symbols for the plots and the legends
#--------------------------------------------

symbolsByFactor <- function(experimentFactor) {
  
  #check whether a factor has been provided
  if(class(experimentFactor)!="factor") stop("Parameter 'experimentFactor' must be of class 'factor'")

  #21 selected symbols, in order of preference to be used
  symbolsSelect <- c(15:18,43,42,0:14)

  #create the symbols to plot, and symbols for the legend
  plotSymbols <- symbolsSelect[as.numeric(experimentFactor)]
  legendSymbols <- symbolsSelect[1:max(as.numeric(experimentFactor))]
	
  return(list(plotSymbols=plotSymbols,legendSymbols=legendSymbols))

}


################
## boxplotFun ##
################

boxplotFun <- function(Data, experimentFactor=NULL, plotColors=NULL, legendColors=NULL, normMeth="", postfix="", 
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {
  
  if(is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
  if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
  

  if(normMeth=="") {
    tmain <- "Boxplot of signal intensities"
    tmtext2 <- "Intensity\n\n\n"
  } else {
    library(gdata) #for trim function
	tmain <- paste("Boxplot of signal intensities after ", trim(normMeth), " normalisation", sep="")
    tmtext2 <- "Normalised intensity\n\n\n"
  }

  png(file = paste("Boxplot",ifelse(postfix!="","_",""),postfix,".png", sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)  
  par(oma=c(17,0,0,0), cex.axis=1) 
  suppressWarnings(boxplot(Data, col=plotColors ,main=tmain, axes=FALSE, pch = 20, cex=0.7))
  if(length(levels(experimentFactor))>1){ 
    legend("topright", levels(experimentFactor),
       col=legendColors,fill=legendColors, cex = 0.7, bg = "white", bty = "o")
  }
	if(length(colnames(Data))<MAXARRAY){
		cexval <- 0.65
	}else{
		cexval <- 0.45
	}  
  axis(1,at=1:length(colnames(Data)),las=2,labels=colnames(Data), cex.axis=cexval)        
  axis(2, cex.axis=0.7)  
  mtext(tmtext2, side=2, cex=0.8)  	   
  #mtext("Distributions should be comparable between samples\n", side=3, font=1, cex=0.7)
  dev.off()
}


################
## densityFun ##
################

densityFun <- function(Data, plotColors=NULL, normMeth="", postfix="",
  WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41){           
  
  if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
  
  png(file = paste("DensityHistogram",ifelse(postfix!="","_",""),postfix,".png", sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
  if(length(names(Data))<MAXARRAY){
	cexval <- 0.65
	par(oma=c(12,0,0,0) )
  }else{
	cexval <- 0.45
	par(oma=c(0.1,0,0,0) )
  }    
  
  #determine ranges
  rangeX <- range(apply(Data,2,function(s) range(density(s)$x)))
  rangeY <- range(apply(Data,2,function(s) range(density(s)$y)))
  
  #create empty plot
  if(normMeth=="") {
    plot(0, main="Density histogram of signal intensities", cex.axis = 0.7, cex.lab=0.8,
		xlim=rangeX, ylim=rangeY, type="n", xlab="Value", ylab="Density")
  } else {
    library(gdata) #for trim function
	plot(0, main=paste("Density histogram of signal intensities\nafter ", trim(normMeth)," normalisation\n", sep=""), 
		cex.axis = 0.7, cex.lab=0.8, xlim=rangeX, ylim=rangeY, type="n", xlab="Value", ylab="Density")
  } 
  
  #add density curve for each sample
  for (s in 1:(dim(Data)[2])) {
    dens <- density(Data[,s])
	lines(dens, lty=s, col=plotColors[s], lwd=3)
  }
  
  legend("topright", substr(colnames(Data),1,20), lwd=3, lt = 1:length(colnames(Data)),
    col = plotColors, cex = cexval, bty = "n")           
  #mtext( "Curves should be comparable between samples\n", side=3, font=1, cex=0.7)
  dev.off()
}


###########
## maFun ##
###########

maFun <- function(Data, experimentFactor=NULL, perGroup=FALSE, normMeth="", postfix="", 
  WIDTH=1000, HEIGHT=1414, MAXARRAY=41){                       
  
  library(affyPLM) #MAplot function
  
  if(normMeth=="") {
    tmain <- paste("MA plots",ifelse(perGroup,", computed for group",""),sep="")
  } else {
    library(gdata) #for trim function
	tmain <- paste("MA plots after ",trim(normMeth)," normalisation",ifelse(perGroup,", computed for group ",""),sep="")	
  }
  
  #check whether MA plots have to be created per experimental group or for the whole data set at once
  if(perGroup){
    if(is.null(experimentFactor)) stop("When selecting MA plots per experimental group, 'experimentFactor' must be provided")
  } else {
    experimentFactor <- factor(rep("",length(colnames(Data))))
  }

  #plot the images
  for(k in (levels(experimentFactor))) {
    x2 <- Data[,experimentFactor==k]   
	
   #amount of plots per page
	if(length(colnames(x2))<=12) { nPerPage<-6; nCol<- 2; cexmain<-1.9 ; cexval<-1.8} # max 3p of 6 plots
	if(length(colnames(x2))>12) {nPerPage<-15; nCol<- 3; cexmain<-1.3 ; cexval<-1.3}	
	if(length(colnames(x2))>45) {nPerPage<-24 ;nCol<- 4; cexmain<-1.4 ; cexval<-1.2}		  
    
	#subselection of data
    if(length(colnames(x2)) < 2) {
	  warning("group ",ifelse(perGroup,k,"")," consists of only one array: no MA plot can be made",
	    ifelse(perGroup,": consider selecting 'dataset' as option for the MA plots",""))
	} else {
      nplots <- ceiling(length(colnames(x2))/nPerPage)
	  for(l in 1:nplots){     
        png(paste("MAplot",ifelse(postfix!="","_",""),postfix, ifelse(k!="","_",""), gsub("[:/|\\?<>*\"]","_",k), 
		  ifelse(nplots>1,"_",""), ifelse(nplots>1,l,""), ".png", sep=""), width=WIDTH, height=HEIGHT)
		par(mfrow = c((nPerPage/nCol),nCol),oma=c(0,0,6,0),adj=0)  
        from <- nPerPage*l - (nPerPage-1)
        to <-   min(nPerPage*l, length(colnames(x2)))
		eSet_tmp <- ExpressionSet(x2) #MAplot only works on specific objects
        MAplot(eSet_tmp, pairs = FALSE, which= from:to, plot.method = "smoothScatter", lwd=3, cex.main=cexmain,cex=cexval)    
		mtext(paste(tmain,ifelse(perGroup,k,""),ifelse(nplots>1,paste(l,"/",nplots),"")),side = 3, outer = TRUE, font = 2, cex = 2)  
		rm(eSet_tmp)
        dev.off()
      }
    }	
  }

}


###############
## correlFun ##
###############

correlFun <- function(Data, clusterOption1="pearson", clusterOption2="ward.D2", normMeth="", postfix="", 
  experimentFactor=NULL, legendColors=NULL, WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41){  
  
  library(bioDist) #distance functions
  library(gplots) #heatmap.2 function

  if(is.null(experimentFactor)) stop("the 'exerimentFactor' parameter is required")	
  if(is.null(legendColors)) stop("the 'legendColors' parameter is required")	
  
  if(normMeth == "") {
    text1 <- "Sample data correlation plot"
  } else {
    library(gdata) #for trim function
	text1 <- paste("Sample data correlation plot\nafter",trim(normMeth),"normalisation")
  }
  
  if(length(colnames(Data))<2) {
      warning("Only one array in dataset, no correlation plot made")
  } else {         
	  png(file = paste("Correlation",ifelse(postfix!="","_",""),postfix,".png",sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
 	  if(length(colnames(Data))<MAXARRAY) {
	    par(oma=c(17,0,0,0),cex.axis=0.7,cex.main=0.8)
	    #subval <- 10
	  } else {
	    par(oma=c(17,0,0,0),srt=90,las=2,cex.axis=0.5,cex.main=0.8)
	    #subval <- 16
	  }        
	  
    #note: for computing array correlation, euclidean would not make sense
    #only use euclidean distance to compute the similarity of the correlation vectors for the arrays
    COpt1 <- "pearson"
    if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
    crp <- cor(Data, use="complete.obs", method=COpt1)
    
    text1 <- paste(text1,"\ncorrelation method:",COpt1,"\ncluster method:",clusterOption2)
    
    switch(tolower(clusterOption1), 
      "pearson" = {
        my.dist <- function(x) cor.dist(x, abs=FALSE)
      },
      "spearman" = {
        my.dist <- function(x) spearman.dist(x, abs=FALSE)
      },
      "euclidean" = {
        my.dist <- function(x) euc(x)
      }
    )

    my.hclust <- function(d) hclust(d, method=clusterOption2)

    #in order to create some space to put colored symbols as well
    #colnames(Data) <- paste(colnames(Data)," ")
    
    sideColors <- legendColors[as.numeric(experimentFactor)]
    
    heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", symm=TRUE, density.info="density",
              main=text1, dendrogram="row", ColSideColors=sideColors)

    #correlationPlot(Data)    
    #axis(1,side=3,at=seq(from=0.5, to=(length(colnames(Data)))-0.5,by=1),
    #    labels=substr(as.character(colnames(Data)),1,subval),las=2)
    #par(srt=0) 
	  #plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE, 
	  #    frame.plot = FALSE, xlim = c(0, 2), ylim = c(0,2))
	  #text(1,1,text1,cex=1)  
  	dev.off()
  }
}


################
## clusterFun ##
################

clusterFun <- function(Data, experimentFactor=NULL, clusterOption1="pearson", clusterOption2="ward.D2", normMeth="", postfix="", 
  plotColors=NULL, legendColors=NULL, plotSymbols=NULL, legendSymbols=NULL, WIDTH=1000, HEIGHT=1414, POINTSIZE=24, MAXARRAY=41) {

  library(bioDist) #distance functions
  
  if(is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
  if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if(is.null(plotSymbols)) stop("the 'plotSymbols' parameter is required")
  if(is.null(legendSymbols)) stop("the 'legendSymbols' parameter is required")
  
  if(normMeth=="") {
    main <- "Cluster dendrogram of samples"
  } else {
    library(gdata) #for trim function
	main <- paste("Cluster dendrogram of samples\nafter",trim(normMeth),"normalisation")
  }  
  if(length(colnames(Data))<3) {
    warning("Only ",length(colnames(Data))," sample(s) in dataset, no clustering plot made")
  } else {  
    switch(tolower(clusterOption1), 
      "pearson" = {
        correl <- cor.dist(t(Data),abs=FALSE)
      },
      "spearman" = {
        correl <- spearman.dist(t(Data),abs=FALSE)
      },
      "euclidean" = {
        correl <- euc(t(Data))
      }
    )
    if(tolower(clusterOption2)!="ward.d2") {
	  clust <- hclust(correl, method = tolower(clusterOption2))
	} else {
	  clust <- hclust(correl, method = "ward.D2")
	}
		
    png(file = paste("Cluster_",clusterOption1,ifelse(postfix!="","_",""),postfix,".png",sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
	  if(length(colnames(Data))<MAXARRAY) {
		  cexval1 <- 0.75
		  cexval2 <- 1.23
		  cexval3 <- 0.55
	  } else {
  	  cexval1 <- 0.55
  	  cexval2 <- 1.6
      cexval3 <- 0.41
	  }   
    par(cex=cexval1,oma=c(14,1,0,0))	
    par(cex.axis=cexval2,cex.lab=cexval2,cex.main=cexval2)	       
    plot(clust, hang=-1, main=main, xlab=paste("distance:",clusterOption1), sub=paste(" cluster method:",clusterOption2))
    points(1:length(clust$order),rep(0,length(clust$order)),pch=15,col="white",cex=1.5)
    points(1:length(clust$order),rep(0,length(clust$order)),pch=plotSymbols[clust$order],col=plotColors[clust$order],cex=1.25)
    if(length(levels(experimentFactor))>1) { 
      legend("topright",levels(experimentFactor),pch=legendSymbols,col=legendColors,cex=1.25)
    }
    par(cex=cexval3)    
    dev.off()
  }
}


############
## pcaFun ##
############

pcaFun <- function(Data, experimentFactor=NULL, normMeth="", postfix="", scaled_pca=TRUE, plotColors=NULL, 
   legendColors=NULL, plotSymbols=NULL, legendSymbols=NULL, namesInPlot=FALSE, WIDTH=1000, HEIGHT=1414, POINTSIZE=24){
  # Scaled PCA by default
  if(is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
  if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if(is.null(plotSymbols)) stop("the 'plotSymbols' parameter is required")
  if(is.null(legendSymbols)) stop("the 'legendSymbols' parameter is required")

  if(length(colnames(Data))<3) {
    warning("Only",length(colnames(Data)),"sample(s) in dataset, no PCA plot made")
  } else { 

    if(normMeth==""){
      #raw data
      tmain <- "PCA analysis of sample data"
    } else{
      library(gdata) #for trim function
	  tmain <- paste("PCA analysis of sample data\nafter", trim(normMeth), "normalisation")
    }
  
    pca1 <- NULL  
    try(pca1 <- prcomp(t(Data[apply(Data,1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=scaled_pca),TRUE)
    if(is.null(pca1) & scaled_pca) {
      try(pca1 <- prcomp(t(Data[apply(Data,1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=FALSE),TRUE)
      if(!is.null(pca1)) warning("pca with scaling unsuccessful, successfully retried without scaling")
    }
    if(!is.null(pca1)) {
      perc_expl1 <- round(((pca1$sdev[1:3]^2)/sum(pca1$sdev^2))*100,2)

      cex.circle <- 1.5
      cex.text <- 0.7
	    cex.legend <- 0.75
      tcol <- "#444444"
	  
	    png(file = paste("PCAanalysis",ifelse(postfix!="","_",""),postfix,".png",sep=""), width=WIDTH+200*(!namesInPlot), height=HEIGHT+283*(!namesInPlot),
	      pointsize=POINTSIZE)

      if(!namesInPlot) {
	     layout(rbind(c(1,1,2,2,5),c(3,3,4,4,5)))
	    } else {
	     layout(rbind(c(1,1,2,2),c(1,1,2,2),c(3,3,4,4),c(3,3,4,4)))
	    }
	    par(oma=c(20,0,5,0))
        plot(pca1$x[,1],pca1$x[,2],cex=cex.circle,pch=plotSymbols,
          col=plotColors,xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
          ylab=paste("PC2 (",perc_expl1[2],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,1],pca1$x[,2], colnames(Data),pos=4,cex=cex.text,col=tcol) 
        plot(pca1$x[,1],pca1$x[,3],cex=cex.circle,pch=plotSymbols,
        col=plotColors,xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
        ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,1],pca1$x[,3], colnames(Data),pos=4,cex=cex.text,col=tcol)
        plot(pca1$x[,2],pca1$x[,3],cex=cex.circle,pch=plotSymbols,
        col=plotColors,xlab=paste("PC2 (",perc_expl1[2],"%)",sep=""),
        ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
        if(namesInPlot) text(pca1$x[,2],pca1$x[,3], colnames(Data),pos=4,cex=cex.text,col=tcol)
        barplot((100*pca1$sdev^2)/sum(pca1$sdev^2),xlab="components",ylab="% of total variance explained")
        
		  if(namesInPlot) {
        if(length(levels(experimentFactor))>1){ 
          legend("topright",levels(experimentFactor),
          pch=legendSymbols,col=legendColors,cex=cex.legend)
        }
      } else {
			  par(mar=c(0,0,0,0))	
			  plot(1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
			  if(length(levels(experimentFactor))>1) {
		      legend("topleft",c(levels(experimentFactor),"",colnames(Data)),
 #             pch=c(rep(20,length(unique(experimentFactor))+1),plotSymbols,
			      pch=c(legendSymbols,20,plotSymbols),
		          col=c(legendColors,"white",plotColors),cex=(cex.legend+0.1)
 #             ,fill=c(legendColors,rep("white",length(experimentFactor)+1)),
 #             border=c(legendColors,rep("white",length(experimentFactor)+1))
			    )
			  } else {
			    legend("topleft",colnames(Data),pch=plotSymbols,
			      col=plotColors,cex=0.7, bty = "n")
			  }
      }
        
		  mtext(tmain, side = 3, outer = TRUE, font = 2, cex = 1.2)
      dev.off()
    } else {
	    warning("PCA on the",Type,"data set unsuccessful, image skipped")
    }
  }
}


###################
## createQCPlots ##
###################

#function to create all (except MA) QC plots in one go, coloured for one or several factors

createQCPlots <- function(Data, factors, Table=NULL, normMeth="", postfix="") {
# Data: the data matrix or data frame (no non-numerical columns)
#   must be provided
# factors: the vector of factor variable names to be taken from the *parent environment* to create QC plots
#   must be provided
# Table: an optional table (data.frame, matrix) to take the factor variables from, by selecting columns
# normMeth: the name of the normalisation method, if applicable
#   default empty string
#   note: the function does not perform normalisation, this is just to indicate what has been applied
# postfix: label to be added to the plot file names
#   default empty string
#   note: the name is automatically extended with the name of each factor

	for (expF in factors) {
		
		#load the factor variable
		#in case a data table is given to get factors from, adapt the value to be evaluated
		if(is.null(Table)) {
			experimentFactor <- eval(parse("",-1,expF))  #get(expF)
		} else {
			experimentFactor <- eval(parse("",-1,paste("Table[,\"",expF,"\"]",sep="")))
		}

		#make numeric variables into factors as well, giving a warning
		if(!is.factor(experimentFactor)) {
			experimentFactor <- as.factor(experimentFactor)
			warning("Variable ",expF," is not factorial, automatically used as a factor")
		}
		
		#create colours for the plots
		colList <- colorsByFactor(experimentFactor)
		plotColors <- colList$plotColors
		legendColors <- colList$legendColors
		rm(colList)

		#create symbols for the plots
		symbList <- symbolsByFactor(experimentFactor)
		plotSymbols <- symbList$plotSymbols
		legendSymbols <- symbList$legendSymbols
		rm(symbList)		
		
		#construct postfix for file names, containing factor name
		postfix_ext <- paste(postfix,sub("_","",expF),sep="_")
		
		#create boxplot
		try(boxplotFun(Data, experimentFactor, plotColors, legendColors, normMeth, postfix_ext),TRUE)
		
		#create density plot
		try(densityFun(Data, plotColors, normMeth, postfix_ext),TRUE)
		
		#create cluster plot
		try(clusterFun(Data, experimentFactor, clusterOption1="pearson", clusterOption2="ward.D2", normMeth, postfix_ext, 
				  plotColors, legendColors, plotSymbols, legendSymbols),TRUE)

		#create heatmap (of array correlation)
		try(correlFun(Data, clusterOption1="pearson", clusterOption2="ward.D2", normMeth, postfix_ext, 
				 experimentFactor, legendColors),TRUE)

		#create PCA plot
		try(pcaFun(Data, experimentFactor, normMeth, postfix_ext, scaled_pca=TRUE, plotColors, legendColors,
			   plotSymbols, legendSymbols, namesInPlot=FALSE),TRUE)
			   
		try(pcaFun(Data, experimentFactor, normMeth, paste(postfix_ext,"2",sep=""), scaled_pca=TRUE, plotColors, legendColors,
			   plotSymbols, legendSymbols, namesInPlot=TRUE),TRUE)

	}
	
}


####################
## saveStatOutput ##
####################

saveStatOutput <- function(design,fit,annotation=NULL,na.rm=TRUE,includeIntercept=FALSE,includeB=FALSE,createHists=TRUE,postfix="") {
  files <- NULL
  for (i in 1:length(colnames(design))) {
	if(includeIntercept | length(grep("intercept",colnames(design)[i],ignore.case=TRUE))==0) {
		cat(paste("--[[ Saving table for coefficient ", colnames(design)[i], " ]]--\n", sep="\t"))
		toptab <- topTable(fit, adjust.method="BH", coef=i, number=dim(fit)[1], resort.by="P")
		#toptab[,1] <- substring(toptab[,1],1,15)
		if(!includeB) {
		  toptab <- toptab[,-match("B",colnames(toptab))]
		}
		if(!is.null(annotation)) {
		  toptab <- cbind(toptab,annotation[match(rownames(toptab),rownames(annotation)),])
		}
		fc <- as.matrix(2^toptab[,"logFC"])
		fc[(toptab[, "logFC"] < 0)&(!is.na(toptab[, "logFC"]))] <- -1/fc[(toptab[, "logFC"] < 0)&(!is.na(toptab[, "logFC"]))]
		colnames(fc) <- "Fold Change"
		m.col <- grep("logFC",colnames(toptab))
		toptab <- cbind(toptab[,1:m.col,drop=FALSE],fc,toptab[,(m.col+1):dim(toptab)[2],drop=FALSE])
		if(na.rm) {
		  cat("----[[ ",sum(is.na(toptab[,"P.Value"]))," probes with NA estimates removed ]]\n")
		  toptab <- toptab[!is.na(toptab[,"P.Value"]),]
		}
		filename <- paste("table_",colnames(design)[i],ifelse(postfix!="","_",""),postfix,".tab",sep="")
		write.table(toptab,file=filename,sep="\t",col.names=NA)
		files <- c(files, filename)
		if (createHists) {
		  #also save p value histogram
		  png(paste("pvalue_hist_",colnames(design)[i],"_",postfix,".png",sep=""),width=1000,height=1000)
		  hist(toptab[,"P.Value"],main=paste("p-value histogram for",colnames(design)[i]),xlab="p-values",col="blue",breaks=120,cex.axis=1.2,cex.lab=1.2)
		  dev.off()
		  #and of adapted FC
		  png(paste("FC_hist_",colnames(design)[i],"_",postfix,".png",sep=""),width=1000,height=1000)
		  hist(toptab[,"Fold Change"],main=paste("adapted fold change histogram for",colnames(design)[i]),xlab="adapted fold changes",col="green3",breaks=120,cex.axis=1.2,cex.lab=1.2)
		  dev.off()
		}
	}
  }
  return(files)
}


##########################
## saveStatOutputDESeq2 ##
##########################

saveStatOutputDESeq2 <- function(design,dds,annotation=NULL,na.rm=TRUE,includeIntercept=FALSE,createHists=TRUE,postfix="") {
  files <- NULL
  for (i in 1:length(colnames(design))) {
	if(includeIntercept | length(grep("intercept",colnames(design)[i],ignore.case=TRUE))==0) {
		cat(paste("--[[ Saving table for coefficient ", colnames(design)[i], " ]]--\n", sep="\t"))
		res <- results(dds, contrast=design[,i], format="DataFrame")
		toptab <- as.data.frame(res)
		if(!is.null(annotation)) {
		  toptab <- cbind(toptab,annotation[match(rownames(toptab),rownames(annotation)),])
		}
		fc <- as.matrix(2^toptab[,"log2FoldChange"])
		fc[(toptab[, "log2FoldChange"] < 0)&(!is.na(toptab[, "log2FoldChange"]))] <- -1/fc[(toptab[, "log2FoldChange"] < 0)&(!is.na(toptab[, "log2FoldChange"]))]
		colnames(fc) <- "FoldChange"
		m.col <- grep("log2FoldChange",colnames(toptab))
		toptab <- cbind(toptab[,1:m.col,drop=FALSE],fc,toptab[,(m.col+1):dim(toptab)[2],drop=FALSE])
		if(na.rm) {
		  cat("----[[ ",sum(is.na(toptab[,"pvalue"]))," probes with NA estimates removed ]]\n")
		  toptab <- toptab[!is.na(toptab[,"pvalue"]),]
		}
		filename <- paste("table_",colnames(design)[i],ifelse(postfix!="","_",""),postfix,".tab",sep="")
		write.table(toptab,file=filename,sep="\t",col.names=NA)
		files <- c(files, filename)
		if (createHists) {
		  #also save p value histogram
		  png(paste("pvalue_hist_",colnames(design)[i],"_",postfix,".png",sep=""),width=1000,height=1000)
		  hist(toptab[,"pvalue"],main=paste("p-value histogram for",colnames(design)[i]),xlab="p-values",col="blue",breaks=120,cex.axis=1.2,cex.lab=1.2)
		  dev.off()
		  #and of adapted FC
		  png(paste("FC_hist_",colnames(design)[i],"_",postfix,".png",sep=""),width=1000,height=1000)
		  hist(toptab[,"FoldChange"],main=paste("adapted fold change histogram for",colnames(design)[i]),xlab="adapted fold changes",col="green3",breaks=120,cex.axis=1.2,cex.lab=1.2)
		  dev.off()
		}
	}
  }
  return(files)
}


#####################################
###Create significant pvalue table###
#####################################

createPvalTab <- function(files,postfix="",pvaluelist=c(0.001,0.01,0.05,0.1),
					adjpvaluelist=0.05,foldchangelist=c(1.1,1.2,1.5,2),
					namePVal="P.Value",nameAdjPVal="adj.P.Val",
					nameFC="Fold.Change",nameLogFC="logFC",html=FALSE) {

  if(is.null(files)) stop("vector of one or more file names has to be provided")

  if(html) library(R2HTML)

  #extract data from the comparison tables
  for(i in 1:length(files)){
	#read file
	tab <- read.delim(files[i],header=TRUE,as.is=TRUE)
	
	#create headers
	if(i==1) {
	  numberOfGenes <- dim(tab)[1]
	  ### P-values
	  #expected
	  expPval <- as.vector(apply(as.data.frame(pvaluelist),1,function(x) { 
				c(y <- ceiling(numberOfGenes*x),ceiling(0.5 * y),ceiling(0.5 * y))}))

	  P.Values <- c("Comparisons",paste("pVal <",rep(pvaluelist,each=3),c("total","up","down")),"NAs")
	  Expected <- c("Number expected",expPval,0)
	  P.Values <- cbind(P.Values,Expected)
	   
	  ### Adjusted p-Values
	  Adj.pValues <- c("Comparisons",paste("adj.pVal <",rep(adjpvaluelist,each=3),c("total","up","down")),"NAs")

	  ### Fold changes
	  Fold.Changes <- c("Comparisons",paste(rep(c("|FC| >=","FC >=","FC <="),3),rep(foldchangelist,each=3),c("total","up","down")),"NAs")
	}
	
	#remove extension from file name for display in table
	nameNoExt <- paste(strsplit(files[i],".",fixed=TRUE)[[1]][-length(strsplit(files[1],".",fixed=TRUE)[[1]])],collapse=".")

	### P-values
	rowsP <- NULL
	for (count in 1:length(pvaluelist)){
		noOfPval <- up <- down <- 0
		#compare p-value with entered p-values 
		noOfPval <- sum(tab[,namePVal] < pvaluelist[count],na.rm=TRUE)
		up <- sum((tab[,namePVal] < pvaluelist[count]) & (tab[,nameLogFC] > 0.58),na.rm=TRUE)
		down <- sum((tab[,namePVal] < pvaluelist[count]) & (tab[,nameLogFC] < (-0.58)),na.rm=TRUE)
		rowsP <- c(rowsP,noOfPval,up,down)
	}
	P.Values<-cbind(P.Values,c(nameNoExt,rowsP,sum(is.na(tab[,namePVal]))))

	### Adjusted p-Values
	rowsAP <- NULL
	for (count in 1:length(adjpvaluelist)){
		adjNoOfPval <- adjUp <- adjDown <- 0
		#compare adjusted p-values with entered adjusted p-values 
		adjNoOfPval <- sum(tab[,nameAdjPVal] < adjpvaluelist[count],na.rm=TRUE)	
		adjUp <- sum((tab[,nameAdjPVal] < adjpvaluelist[count]) & (tab[,nameLogFC] > 0.58),na.rm=TRUE)
		adjDown <- sum((tab[,nameAdjPVal] < adjpvaluelist[count]) & (tab[,nameLogFC] < (-0.58) ),na.rm=TRUE)
		rowsAP <- c(rowsAP,adjNoOfPval,adjUp,adjDown)
	}
	Adj.pValues<-cbind(Adj.pValues,c(nameNoExt,rowsAP,sum(is.na(tab[,nameAdjPVal]))))

	### Fold changes
	rowsFC <- NULL
	for (count in 1:length(foldchangelist)){
		FCTot <- FCUp <- FCDown <- 0
		#compare FC with entered FC 
		FCTot <- sum(abs(tab[,nameFC]) >= foldchangelist[count] & (!is.na(tab[,nameFC])))
		FCUp <- length(tab[(tab[,nameFC] >= foldchangelist[count]) & (!is.na(tab[,nameFC])),1.5])
		FCDown <- length(tab[(tab[,nameFC] <= (-foldchangelist[count])) & (!is.na(tab[,nameFC])),1.5])
		rowsFC <- c(rowsFC,FCTot,FCUp,FCDown)
	}
	Fold.Changes<-cbind(Fold.Changes,c(nameNoExt,rowsFC,sum(is.na(tab[,nameFC]))))
	
  }

  colnames(P.Values) <- NULL
  colnames(Adj.pValues) <- NULL
  colnames(Fold.Changes) <- NULL
  
  #write table to tab delimited text file 
  filename <- paste("Summary_tables",ifelse(postfix!="","_",""),postfix,".tab",sep="")
  cat (paste("--[[ Saving",filename,"]]--\n"))
  write.table(t(cbind(c("Total number of measurements",numberOfGenes),rep("",2))),file=filename,sep="\t",
			row.names=FALSE,col.names=FALSE, quote=FALSE)
  write.table(t(P.Values),file=filename,sep="\t",
			row.names=FALSE,col.names=FALSE,append=TRUE, quote=FALSE)
  write.table(t(cbind(rep("",dim(Adj.pValues)[1]),Adj.pValues)),file=filename,sep="\t",
			row.names=FALSE,col.names=FALSE,append=TRUE, quote=FALSE)
  write.table(t(cbind(rep("",dim(Fold.Changes)[1]),Fold.Changes)),file = filename,sep="\t",
			row.names=FALSE,col.names=FALSE,append=TRUE, quote=FALSE)
  
  #if requested, write table to HTML file
  if(html) {
    filenameHTML <- paste("Summary_tables",ifelse(postfix!="","_",""),postfix,".html",sep="")  
	cat (paste("--[[ Saving",filenameHTML,"]]--\n"))
    HTML(t(c("Total number of measurements",numberOfGenes)),file=filenameHTML,
			row.names=TRUE,append=FALSE)
	HTML(t(P.Values),file=filenameHTML,row.names=TRUE)
	HTML(t(Adj.pValues),file=filenameHTML,row.names=TRUE)
	HTML(t(Fold.Changes),file=filenameHTML,row.names=TRUE)
  }
  
}
