#methods-CuffData.R
#
#Author: Loyal A. Goff
#
#
####################

##################
#Initialize
##################
setMethod("initialize","CuffData",
			function(.Object,
					DB,
					tables=list(mainTable = "",dataTable = "",expDiffTable = "",featureTable = "",countTable = "", replicateTable = ""),
					filters=list(),
					type = c("genes","isoforms","TSS","CDS"),
					idField,
					... ){
				.Object<-callNextMethod(.Object,
						DB = DB,
						tables = tables,
						filters = filters,
						type = type,
						idField = idField,
						...)				
		}
)

setValidity("CuffData",function(object){
		TRUE
		}
)			

################
#Class Methods
################
setMethod("show","CuffData",
		function(object){
			size<-dim(object)
			cat(class(object), "instance with:\n\t",size[1],"features and",size[2],"samples\n")
		}
)

setMethod("dim","CuffData",
		function(x){
			countQuery<-paste("SELECT COUNT(",x@idField,") as n FROM ",x@tables$mainTable)
			nIds<-dbGetQuery(x@DB,countQuery)
			sampleQuery<-("SELECT COUNT(sample_name) as n FROM samples")
			nSamples<-dbGetQuery(x@DB,sampleQuery)
			c(nIds$n,nSamples$n)
		}
)

#####################
#Feature Table
#####################

.addFeatures<-function(object,features,...){
	if(!is.data.frame(features)){
		stop("features must be a data.frame")
	}
	colnames(features)[1]<-object@idField
	colnames(features)<-make.db.names(object@DB,colnames(features),unique=T)
	dbWriteTable(object@DB,object@tables$featureTable,features,row.names=F,overwrite=T)
	indexQuery<-paste("CREATE INDEX ",object@idField," ON ", object@tables$featureTable," (",object@idField,")",sep="")
	res<-dbGetQuery(object@DB,indexQuery)
}

setMethod("addFeatures",signature="CuffData",.addFeatures)

###################
#Accessors
###################
.features<-function(object){
	featureQuery<-paste("SELECT * FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField,sep="")
	dbGetQuery(object@DB, featureQuery)
}

setMethod("features","CuffData",.features)

.featureNames<-function(object){
	featureQuery<-paste("SELECT ",object@idField," FROM ",object@tables$mainTable, sep="")
	res<-dbGetQuery(object@DB,featureQuery)
	res[,1]
}

setMethod("featureNames","CuffData",.featureNames)

.samples<-function(object){
	res<-dbReadTable(object@DB,'samples')
	res<-res$sample_name
	res
}

setMethod("samples","CuffData",.samples)

.replicates<-function(object){
	res<-dbReadTable(object@DB,'replicates')
	res<-res$rep_name
	res
}

setMethod("replicates","CuffData",.replicates)

.fpkm<-function(object,features=FALSE,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	#Sample Search String (SQL)
	sampleString<-'('
	for (i in myLevels){
		sampleString<-paste(sampleString,"'",i,"',",sep="")
	}
	sampleString<-substr(sampleString,1,nchar(sampleString)-1)
	sampleString<-paste(sampleString,")",sep="")
	
	if(!features){
		FPKMQuery<-paste("SELECT * FROM ",object@tables$dataTable," WHERE sample_name IN ",sampleString,sep="")
	}else{
		FPKMQuery<-paste("SELECT xf.*,xm.*,x.sample_name,x.fpkm,x.conf_hi,x.conf_lo FROM ",object@tables$dataTable," x LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField," LEFT JOIN ",object@tables$mainTable," xm ON x.",object@idField,"=xm.",object@idField," WHERE x.sample_name IN ",sampleString,sep="")
		#print(FPKMQuery)
	}
	res<-dbGetQuery(object@DB,FPKMQuery)
	res$sample_name<-factor(res$sample_name,levels=getLevels(object))
	res
}

setMethod("fpkm","CuffData",.fpkm)

.repFpkm<-function(object,features=FALSE,repIdList){
	#Sample subsetting
	if(!missing(repIdList)){
		if(.checkReps(object@DB,repIdList)){
			myLevels<-repIdList
		}else{
			stop("Replicate does not exist!")
		}
	}else{
		myLevels<-getRepLevels(object)
	}
	#Sample Search String (SQL)
	sampleString<-'('
	for (i in myLevels){
		sampleString<-paste(sampleString,"'",i,"',",sep="")
	}
	sampleString<-substr(sampleString,1,nchar(sampleString)-1)
	sampleString<-paste(sampleString,")",sep="")
	
	if(!features){
		FPKMQuery<-paste("SELECT * FROM ",object@tables$replicateTable," WHERE rep_name IN ",sampleString,sep="")
	}else{
		FPKMQuery<-paste("SELECT xf.*,xm.*,x.rep_name,x.raw_frags,x.internal_scaled_frags,x.external_scaled_frags,x.fpkm,x.effective_length,x.status FROM ",object@tables$replicateTable," x LEFT JOIN ",object@tables$featureTable," xf on x.",object@idField,"=xf.",object@idField," LEFT JOIN ",object@tables$mainTable," xm ON x.",object@idField,"=xm.",object@idField," WHERE x.rep_name IN ",sampleString,sep="")
	}
	#print(FPKMQuery)
	res<-dbGetQuery(object@DB,FPKMQuery)
	res$rep_name<-factor(res$rep_name,levels=getRepLevels(object))
	res
}

setMethod("repFpkm","CuffData",.repFpkm)

.count<-function(object,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	#Sample Search String (SQL)
	sampleString<-'('
	for (i in myLevels){
		sampleString<-paste(sampleString,"'",i,"',",sep="")
	}
	sampleString<-substr(sampleString,1,nchar(sampleString)-1)
	sampleString<-paste(sampleString,")",sep="")
	
	CountQuery<-FPKMQuery<-paste("SELECT * FROM ",object@tables$countTable," WHERE sample_name IN ",sampleString,sep="")
	
	res<-dbGetQuery(object@DB,CountQuery)
	res$sample_name<-factor(res$sample_name,levels=getLevels(object))
	res
	
}

setMethod("count","CuffData",.count)

.fpkmMatrix<-function(object,fullnames=FALSE,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	
	samp<-samples(object)
	FPKMMatQuery<-paste("select x.",object@idField,", ",sep="")
	for (i in samp){
		FPKMMatQuery<-paste(FPKMMatQuery,"sum(case when xd.sample_name ='",i,"' then fpkm end) as ",i,",",sep="")
	}
	FPKMMatQuery<-substr(FPKMMatQuery, 1, nchar(FPKMMatQuery)-1)
	FPKMMatQuery<-paste(FPKMMatQuery," from ",object@tables$mainTable," x LEFT JOIN ",object@tables$dataTable," xd on x.",object@idField," = xd.",object@idField," group by x.",object@idField,sep="")
	res<-dbGetQuery(object@DB,FPKMMatQuery)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(sampleIdList)){
		res<-data.frame(res[,sampleIdList],row.names=rownames(res))
		colnames(res)<-sampleIdList
	}
	res
}

setMethod("fpkmMatrix","CuffData",.fpkmMatrix)

.repFpkmMatrix<-function(object,fullnames=FALSE,repIdList){
	#Sample subsetting
	if(!missing(repIdList)){
		if(.checkReps(object@DB,repIdList)){
			myLevels<-repIdList
		}else{
			stop("Replicate does not exist!")
		}
	}else{
		myLevels<-getRepLevels(object)
	}
	
	samp<-replicates(object)
	FPKMMatQuery<-paste("select x.",object@idField,", ",sep="")
	for (i in samp){
		FPKMMatQuery<-paste(FPKMMatQuery,"sum(case when xd.rep_name ='",i,"' then fpkm end) as ",i,",",sep="")
	}
	FPKMMatQuery<-substr(FPKMMatQuery, 1, nchar(FPKMMatQuery)-1)
	FPKMMatQuery<-paste(FPKMMatQuery," from ",object@tables$mainTable," x LEFT JOIN ",object@tables$replicateTable," xd on x.",object@idField," = xd.",object@idField," group by x.",object@idField,sep="")
	res<-dbGetQuery(object@DB,FPKMMatQuery)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(repIdList)){
		res<-data.frame(res[,repIdList],row.names=rownames(res))
		colnames(res)<-repIdList
	}
	res
}

setMethod("repFpkmMatrix","CuffData",.repFpkmMatrix)

.countMatrix<-function(object,fullnames=FALSE,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	
	samp<-samples(object)
	CountMatQuery<-paste("select x.",object@idField,", ",sep="")
	for (i in samp){
		CountMatQuery<-paste(CountMatQuery,"sum(case when xd.sample_name ='",i,"' then count end) as ",i,",",sep="")
	}
	CountMatQuery<-substr(CountMatQuery, 1, nchar(CountMatQuery)-1)
	CountMatQuery<-paste(CountMatQuery," from ",object@tables$mainTable," x LEFT JOIN ",object@tables$countTable," xd on x.",object@idField," = xd.",object@idField," group by x.",object@idField,sep="")
	res<-dbGetQuery(object@DB,CountMatQuery)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(sampleIdList)){
		res<-data.frame(res[,sampleIdList],row.names=rownames(res))
		colnames(res)<-sampleIdList
	}
	res
}

setMethod("countMatrix","CuffData",.countMatrix)

#This needs a lot of work...
#TODO: Change this to remove lnFcCutoff but make sure that functions that rely on diffData have their own FC cutoff so that plotting doesn't suffer
.diffData<-function(object,x,y,features=FALSE){
	if(missing(x) && missing(y)){
		if(!features){
			diffQuery<-paste("SELECT * FROM ",object@tables$expDiffTable,sep="")
		}else{
			diffQuery<-paste("SELECT * FROM ",object@tables$expDiffTable," x LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField,sep="")
		}
	}else if (missing(x) || missing(y)){
		stop("You must supply both x and y or neither.")
	}else{
		if(!features){
			diffQuery<-paste("SELECT x.",object@idField,", xed.* FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$expDiffTable," xed on x.",object@idField," = xed.",object@idField," WHERE ((sample_1 = '",x,"' AND sample_2 = '",y,"') OR (sample_1 = '",y,"' AND sample_2 = '",x,"'))",sep="")
		}else{
			diffQuery<-paste("SELECT x.",object@idField,", xed.*, xf.* FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$expDiffTable," xed on x.",object@idField," = xed.",object@idField," LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField," WHERE ((sample_1 = '",x,"' AND sample_2 = '",y,"') OR (sample_1 = '",y,"' AND sample_2 = '",x,"'))",sep="")
		}
	}
	dat<-dbGetQuery(object@DB,diffQuery)
	#diffQuery
	dat
}

setMethod("diffData",signature(object="CuffData"),.diffData)

.getMA<-function(object,x,y,logMode=T,pseudocount=1){
	if (missing(x) || missing(y)){
		stop("You must supply both x and y.")
	}else{
		sql<-paste("SELECT x.",object@idField,", sum(case when x.sample_name = '",x,"' then x.fpkm end) AS 'x', sum(case when x.sample_name = '",y,"' then x.fpkm end) AS 'y' FROM ",object@tables$dataTable," x GROUP BY x.",object@idField,";",sep="")
		#print(sql)
		dat<-dbGetQuery(object@DB,sql)
		
		if(logMode){
			dat$x<-log10(dat$x+pseudocount)
			dat$y<-log10(dat$y+pseudocount)
		}
		dat$A<-(dat$x+dat$y)/2
		dat$M<-dat$x/dat$y
		res<-dat[,c(1,4:5)]
		res
	}
}

.getCountMA<-function(object,x,y,logMode=T,pseudocount=1){
	if (missing(x) || missing(y)){
		stop("You must supply both x and y.")
	}else{
		sql<-paste("SELECT x.",object@idField,", sum(case when x.sample_name = '",x,"' then x.count end) AS 'x', sum(case when x.sample_name = '",y,"' then x.count end) AS 'y' FROM ",object@tables$countTable," x GROUP BY x.",object@idField,";",sep="")
		dat<-dbGetQuery(object@DB,sql)
		
		if(logMode){
			dat$x<-log10(dat$x+pseudocount)
			dat$y<-log10(dat$y+pseudocount)
		}
		dat$A<-(dat$x+dat$y)/2
		dat$M<-dat$x/dat$y
		res<-dat[,c(1,4:5)]
		res
	}
}

setMethod("DB","CuffData",function(object){
		return(object@DB)
		})

setMethod("tables","CuffData",function(object){
		return(object@tables)
		})

setMethod("filters","CuffData",function(object){
		return(object@filters)
		})

setMethod("type","CuffData",function(object){
		return(object@type)
		})

setMethod("idField","CuffData",function(object){
		return(object@idField)
		})

##################
#Setters
##################


##################
#Subsetting
##################
#Example query
#"SELECT * FROM genes WHERE gene_id in ('XLOC_000005','XLOC_000015','XLOC_000055','XLOC_000595','XLOC_005998','ucscCodingXLOC_018816')"
.getLevels<-function(object){
	levelsQuery<-'SELECT s.sample_name FROM samples s ORDER BY s.sample_index ASC'
	levels<-dbGetQuery(object@DB,levelsQuery)$sample_name
	levels
}

setMethod("getLevels",signature(object="CuffData"),.getLevels)

.getRepLevels<-function(object){
	levelsQuery<-'SELECT r.rep_name FROM replicates r JOIN samples s ON r.sample_name=s.sample_name ORDER BY s.sample_index ASC'
	levels<-dbGetQuery(object@DB,levelsQuery)$rep_name
	levels
}

setMethod("getRepLevels",signature(object="CuffData"),.getRepLevels)

#Useful SQL commands

#SELECT g.gene_id, g.class_code, g.nearest_ref_id, g.gene_short_name, g.locus, g.length, g.coverage, g.status, gd.sample_name, gd.fpkm, gd.conf_hi, gd.conf_lo FROM genes g LEFT JOIN geneData gd ON g.gene_id = gd.gene_id WHERE (g.gene_id = 'XLOC_000001');

#SELECT g.gene_id, ged.* FROM genes g LEFT JOIN geneExpDiffData ged on g.gene_id = ged.gene_id WHERE ((sample_1 = 'H1_hESC' AND sample_2 = 'Fibroblasts') OR (sample_1 = 'Fibroblasts' AND sample_2 = 'H1_hESC')) AND ged.log2_fold_change>-20 AND ged.log2_fold_change<20 ;

#Pivot table SQL for scatterplots
#select g.*, sum(case when gd.sample_name = 'Fibroblasts' then fpkm end) as Fibroblasts, sum(case when gd.sample_name = 'H1_hESC' then fpkm end) as H1_hESC from genes g LEFT JOIN geneData gd on g.gene_id = gd.gene_id group by g.gene_id;


##################
#Plotting
##################

.density<-function(object, logMode = TRUE, pseudocount=1.0, labels, features=FALSE, replicates=FALSE,...){
	if(is(object,'CuffData')) {
		if(replicates){
			dat<-repFpkm(object,features=features)
			colnames(dat)[colnames(dat)=="rep_name"]<-"condition"
		}else{
			dat<-fpkm(object,features=features)
			colnames(dat)[colnames(dat)=="sample_name"]<-"condition"
		}
	} else {
		stop('Un-supported class of object.')
	}
	if(logMode) dat$fpkm<-dat$fpkm+pseudocount
	p<-ggplot(dat)
		if(logMode) {
			p<-p+geom_density(aes(x= log10(fpkm),group=condition,color=condition,fill=condition),alpha=I(1/3))
		}else{
			p<-p+geom_density(aes(x=fpkm,group=condition,color=condition,fill=condition),alpha=I(1/3))
		}
	
	p<-p + opts(title=object@tables$mainTable)
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
	
	#TODO: Add label callout
	p
}

setMethod("csDensity",signature(object="CuffData"),.density)

.scatter<-function(object,x,y,logMode=TRUE,pseudocount=1.0,labels,smooth=FALSE,colorByStatus=FALSE, drawRug=TRUE, ...){
	dat<-fpkmMatrix(object)
	samp<-samples(object)
	
	#check to make sure x and y are in samples
	if (!all(c(x,y) %in% samp)){
		stop("One or more values of 'x' or 'y' are not valid sample names!")
	}
	
	#add pseudocount if necessary
	if(logMode){
		for (i in samp){
			dat[[i]]<-dat[[i]]+pseudocount
		}
	}
	
	#make plot object
	p<-ggplot(dat)
	p<- p + aes_string(x=x,y=y)
	
	#Right now, this does nothing, because 'significant' is not returned from fpkmMatrix object so I don't have this as a feature to draw
	if(colorByStatus){
		p<- p + geom_point(size=1.2,alpha=I(1/3))
	}else{
		p<- p + geom_point(size=1.2,alpha=I(1/3))
	}
	#Add symmetry line
	p<- p + geom_abline(intercept=0,slope=1,linetype=2) 
	
	#Add rug
	if(drawRug){
		p<- p + geom_rug(size=0.5,alpha=0.01)
	}
	
	#add smoother
	if(smooth){
		p <- p + stat_smooth(method="lm",fill="blue",alpha=0.2)
	}
	
	#Add highlights from labels
#	if(!missing(labels)){
#		labelIdx<-fData(object)$gene_short_name %in% labels
#		labelfp<-fp[labelIdx,]
#		labelfp$gene_short_name<-fData(object)$gene_short_name[labelIdx]
#		#print(head(labelfp))
#		p <- p + geom_point(data=labelfp,size=1.2,color="red")
#		p <- p + geom_text(data=labelfp,aes(label=gene_short_name),color="red",hjust=0,vjust=0,angle=45,size=2)
#	}
#	
	#logMode
	if(logMode){
		p <- p + scale_y_log10() + scale_x_log10()
	}
	
	#Add title & Return value
	p<- p + opts(title=object@tables$mainTable)
	p
}

setMethod("csScatter",signature(object="CuffData"), .scatter)

.volcano<-function(object,x,y,alpha=0.05,showSignificant=TRUE,features=FALSE,xlimits=c(-20,20),...){
	samp<-samples(object)
	
	#check to make sure x and y are in samples
	if (!all(c(x,y) %in% samp)){
		stop("One or more values of 'x' or 'y' are not valid sample names!")
	}
	
	dat<-diffData(object=object,features=features)
	
	#subset dat for samples of interest
	mySamples<-c(x,y)
	dat<-dat[(dat$sample_1 %in% mySamples & dat$sample_2 %in% mySamples),]
	dat$significant <- 'no'
	dat$significant[dat$q_value<=alpha]<-'yes'
	s1<-unique(dat$sample_1)
	s2<-unique(dat$sample_2)
	
	p<-ggplot(dat)
	if(showSignificant){
		p<- p + geom_point(aes(x=log2_fold_change,y=-log10(p_value),color=significant),size=0.8)
	}else{
		p<- p + geom_point(aes(x=log2_fold_change,y=-log10(p_value)),size=0.8)
	}
	#Add title and return
	p<- p + opts(title=paste(object@tables$mainTable,": ",s2,"/",s1,sep=""))
	p<- p + scale_colour_manual(values=c("red", "steelblue"))
	
	#Set axis limits
	p<- p + scale_x_continuous(limits=xlimits)
	
	p <- p + xlab(bquote(paste(log[2],"(fold change)",sep=""))) + 
	    ylab(bquote(paste(-log[10],"(p value)",sep="")))
	p
}

setMethod("csVolcano",signature(object="CuffData"), .volcano)

.boxplot<-function(object,logMode=TRUE,pseudocount=0.0001,replicates=FALSE,...){
	if(replicates){
		dat<-repFpkm(object)
		colnames(dat)[colnames(dat)=="rep_name"]<-"condition"
	}else{
		dat<-fpkm(object)
		colnames(dat)[colnames(dat)=="sample_name"]<-"condition"
	}
	if(logMode) {
		dat$fpkm<-dat$fpkm+pseudocount
		p<-ggplot(dat)
		p<-p+geom_boxplot(aes(x=condition,y=log10(fpkm),fill=condition),size=0.3,alpha=I(1/3))
	} else {
		p <- ggplot(dat)
		p<-p+geom_boxplot(aes(x=condition,y=fpkm,fill=condition),alpha=I(1/3),size=0.3)
	}
	p<- p + opts(axis.text.x=theme_text(angle=-90, hjust=0))
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200)
	
	p
}

setMethod("csBoxplot",signature(object="CuffData"),.boxplot)

.dendro<-function(object,logMode=T,pseudocount=1,replicates=FALSE){
	if(replicates){
		fpkmMat<-repFpkmMatrix(object)
	}else{
		fpkmMat<-fpkmMatrix(object)
	}
	if(logMode){
		fpkmMat<-log10(fpkmMat+pseudocount)
	}
	res<-JSdist(makeprobs(fpkmMat))
	#colnames(res)<-colnames(fpkmMat)
	
	#res<-as.dist(res)
	res<-as.dendrogram(hclust(res))
	plot(res,main=paste("All",deparse(substitute(object)),sep=" "))
	res
}

setMethod("csDendro",signature(object="CuffData"),.dendro)

.MAplot<-function(object,x,y,logMode=T,pseudocount=1,smooth=FALSE,useCount=FALSE){
	if(useCount){
		dat<-.getCountMA(object,x,y,logMode=logMode,pseudocount=pseudocount)
	}else{
		dat<-.getMA(object,x,y,logMode=logMode,pseudocount=pseudocount)
	}
	p<-ggplot(dat)
	p<-p+geom_point(aes(x=A,y=log2(M)))
	
	#add smoother
	if(smooth){
		p <- p + stat_smooth(aes(x=A,y=log2(M)),fill="blue")
	}
	#add baseline
	p <- p + geom_hline(yintercept=0)
	
	p
}

setMethod("MAplot",signature(object="CuffData"),.MAplot)

.dispersionPlot<-function(object){
	dat<-count(object)
	p<-ggplot(dat)
	p<-p+geom_point(aes(x=count,y=dispersion,color=sample_name)) + scale_x_log10() + scale_y_log10()
	p
}

setMethod("dispersionPlot",signature(object="CuffData"),.dispersionPlot)

#############
# Other Methods
#############
.specificity<-function(object,logMode=T,pseudocount=1,relative=FALSE,...){
	fpkms<-fpkmMatrix(object,...)
	if(logMode){
		fpkms<-log10(fpkms+pseudocount)
	}
	fpkms<-t(makeprobs(t(fpkms)))
	d<-diag(ncol(fpkms))
	res<-apply(d,MARGIN=1,function(q){
				JSdistFromP(fpkms,q)
			})
	colnames(res)<-paste(colnames(fpkms),"_spec",sep="")
	
	if(relative){
		res<-res/max(res)
	}
	1-res
}

setMethod("csSpecificity",signature(object="CuffData"),.specificity)

.makeRNK<-function(object,x,y,filename){
	#This will make a .rnk file for GSEA pre-ranked analysis using the log2_fold_change value from diffData.  Writes to 'filename'
	
}



#############
#Utility functions
#############
.checkSamples<-function(dbConn,sampleIdList){
	dbSamples<-dbReadTable(dbConn,"samples")
	if (all(sampleIdList %in% dbSamples$sample_name)){
		return(TRUE)
	}else{
		return(FALSE)
	}
}

.checkReps<-function(dbConn,repIdList){
	dbReps<-dbReadTable(dbConn,"replicates")
	if (all(repIdList %in% dbReps$rep_name)){
		return(TRUE)
	}else{
		return(FALSE)
	}
}