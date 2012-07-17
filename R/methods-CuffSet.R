##################
#methods-CuffSet.R
#
#Introduces the CuffSet Class for analysis, manipulation, and plotting of Cufflinks data
#
#Author: Loyal A. Goff
#
##################

#Initialize
setMethod("initialize","CuffSet",
		function(.Object,
				DB,
				runInfo=data.frame(),
				phenoData=data.frame(),
				conditions=data.frame(),
				genes,
				isoforms,
				TSS,
				CDS,
				promoters,
				splicing,
				relCDS,
				...){
			.Object<-callNextMethod(.Object,
					DB = DB,
					runInfo=runInfo,
					phenoData=phenoData,
					conditions = conditions,
					genes = genes,
					isoforms = isoforms,
					TSS = TSS,
					CDS = CDS,
					promoters = promoters,
					splicing = splicing,
					relCDS = relCDS,
					...)				
		}
)

##################
#Class Methods
##################
setMethod("show","CuffSet",
		function(object){
			cat(class(object), "instance with:\n\t",
					dim(object@genes)[2],"samples\n\t",
					dim(object@genes)[1],"genes\n\t",
					dim(object@isoforms)[1],"isoforms\n\t",
					dim(object@TSS)[1],"TSS\n\t",
					dim(object@CDS)[1],"CDS\n\t",
					dim(object@promoters)[1],"promoters\n\t",
					dim(object@splicing)[1],"splicing\n\t",
					dim(object@relCDS)[1],"relCDS\n"
					)
		}
)

#This does not subset appropriately yet
#TODO: Fix for multiple values of i
#
#Solution is to test i to determine if it is of type 'numeric' (index), list (multi-index), or 'character' (gene_ids)
#
#TODO: 	- Add 'j' to select on sampleNames as well
#		- Add ability to search on gene_short_name(s) or featureIDs
		
#setMethod("[",signature(x="CuffSet"),function(x, i, ...){
#			featureIDs<-featureNames(x@genes)[i]
#			if(length(featureIDs)==1){
#				res<-getGene(x,featureID)
#			}else{
#				res<-getGenes(x,featureIDs)
#			}
#			res
#		}
#)


setValidity("CuffSet",
		function(object){
		TRUE	
		}
)

############
#Accessors
############
.samples<-function(object){
	sampleQuery<-"SELECT * FROM samples s LEFT JOIN phenoData p on s.sample_name = p.sample_name"
	dbGetQuery(object@DB,sampleQuery)
}

setMethod("samples",signature(object="CuffSet"),.samples)

.replicates<-function(object){
	replicateQuery<-"SELECT * FROM replicates r"
	dbGetQuery(object@DB,replicateQuery)
}

setMethod("replicates",signature(object="CuffSet"),.replicates)

setMethod("DB","CuffSet",function(object){
		return(object@DB)
		})

.runInfo<-function(object){
	runInfoQuery<-"SELECT * FROM runInfo"
	dbGetQuery(object@DB,runInfoQuery)
}

setMethod("runInfo","CuffSet",.runInfo)

#setMethod("phenoData","CuffSet",function(object){
#			return(object@phenoData)
#		})

setMethod("conditions","CuffSet",function(object){
		return(object@conditions)
		})

setMethod("genes","CuffSet",function(object){
		return(object@genes)
		})

setMethod("isoforms","CuffSet",function(object){
		return(object@isoforms)
		})

setMethod("TSS","CuffSet",function(object){
		return(object@TSS)
		})

setMethod("CDS","CuffSet",function(object){
		return(object@CDS)
		})

setMethod("promoters","CuffSet",function(object){
		return(object@promoters)
		})

setMethod("splicing","CuffSet",function(object){
		return(object@splicing)
		})

setMethod("relCDS","CuffSet",function(object){
		return(object@relCDS)
		})


#make CuffGene objects from a gene_ids
.getGene<-function(object,geneId,sampleIdList=NULL){
	
	#Get gene_id from geneId (can use any identifier now to select gene)
	geneId<-getGeneId(object,geneId)
	
	if(length(geneId)>1){
		stop("More than one gene_id found for given query. Please use getGenes() instead.")
	}
	
	#Sample subsetting
	if(!is.null(sampleIdList)){
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
	sampleString<-paste(sampleString,')',sep="")
	
	whereString = paste("WHERE (x.gene_id ='",geneId,"' OR x.gene_short_name = '",geneId,"')",sep="")
	whereStringFPKM = paste("WHERE (x.gene_id ='",geneId,"' OR x.gene_short_name = '",geneId,"')",' AND (y.sample_name IN ',sampleString,')',sep="")
	whereStringDiff = paste("WHERE (x.gene_id ='",geneId,"' OR x.gene_short_name = '",geneId,"')",' AND (y.sample_1 IN ',sampleString,' AND y.sample_2 IN ',sampleString,')',sep="")
	whereStringRep = paste("JOIN replicates r ON y.rep_name=r.rep_name WHERE (x.gene_id ='",geneId,"' OR x.gene_short_name = '",geneId,"')",' AND (r.sample_name IN ',sampleString,')',sep="")
	
	#dbQueries
	geneAnnotationQuery<-paste("SELECT * from genes x ",whereString,sep="")
	geneFPKMQuery<-paste("SELECT y.* from genes x JOIN geneData y ON x.gene_id=y.gene_id ",whereStringFPKM,sep="")
	#print(geneFPKMQuery)
	geneDiffQuery<-paste("SELECT y.* from genes x JOIN geneExpDiffData y ON x.gene_id=y.gene_id ",whereStringDiff,sep="")
	#print(geneDiffQuery)
	geneRepFPKMQuery<-paste("SELECT y.* from genes x JOIN geneReplicateData y ON x.gene_id=y.gene_id ",whereStringRep,sep="")
	geneCountQuery<-paste("SELECT y.* from genes x JOIN geneCount y ON x.gene_id=y.gene_id ",whereStringFPKM,sep="")
	
	isoformAnnotationQuery<-paste("SELECT * from isoforms i JOIN genes x ON i.gene_id = x.gene_id ",whereString,sep="")
	isoformFPKMQuery<-paste("SELECT y.* from isoforms i JOIN isoformData y ON i.isoform_id = y.isoform_id JOIN genes x ON i.gene_id = x.gene_id ",whereStringFPKM,sep="")
	isoformDiffQuery<-paste("SELECT y.* from isoforms i JOIN isoformExpDiffData y ON i.isoform_id = y.isoform_id JOIN genes x ON i.gene_id = x.gene_id ",whereStringDiff,sep="")
	isoformRepFPKMQuery<-paste("SELECT y.* from isoforms i JOIN isoformReplicateData y ON i.isoform_id = y.isoform_id JOIN genes x ON i.gene_id = x.gene_id ",whereStringRep,sep="")
	isoformCountQuery<-paste("SELECT y.* from isoforms i JOIN isoformCount y ON i.isoform_id = y.isoform_id JOIN genes x ON i.gene_id = x.gene_id ",whereStringFPKM,sep="")
	
	TSSAnnotationQuery<-paste("SELECT * from TSS t JOIN genes x ON t.gene_id = x.gene_id ",whereString,sep="")
	TSSFPKMQuery<-paste("SELECT y.* from TSS t JOIN TSSData y ON t.TSS_group_id=y.TSS_group_id JOIN genes x ON t.gene_id = x.gene_id ",whereStringFPKM,sep="")
	TSSDiffQuery<-paste("SELECT y.* from TSS t JOIN TSSExpDiffData y ON t.TSS_group_id=y.TSS_group_id JOIN genes x ON t.gene_id = x.gene_id ",whereStringDiff,sep="")
	TSSRepFPKMQuery<-paste("SELECT y.* from TSS t JOIN TSSReplicateData y ON t.TSS_group_id=y.TSS_group_id JOIN genes x ON t.gene_id = x.gene_id ",whereStringRep,sep="")
	TSSCountQuery<-paste("SELECT y.* from TSS t JOIN TSSCount y ON t.TSS_group_id=y.TSS_group_id JOIN genes x ON t.gene_id = x.gene_id ",whereStringFPKM,sep="")
	
	
	CDSAnnotationQuery<-paste("SELECT * from CDS c JOIN genes x ON c.gene_id = x.gene_id ",whereString,sep="")
	CDSFPKMQuery<-paste("SELECT y.* from CDS c JOIN CDSData y ON c.CDS_id = y.CDS_id JOIN genes x ON c.gene_id = x.gene_id ",whereStringFPKM,sep="")
	CDSDiffQuery<-paste("SELECT y.* from CDS c JOIN CDSExpDiffData y ON c.CDS_id = y.CDS_id JOIN genes x ON c.gene_id = x.gene_id ",whereStringDiff,sep="")
	CDSRepFPKMQuery<-paste("SELECT y.* from CDS c JOIN CDSReplicateData y ON c.CDS_id = y.CDS_id JOIN genes x ON c.gene_id = x.gene_id ",whereStringRep,sep="")
	CDSCountQuery<-paste("SELECT y.* from CDS c JOIN CDSCount y ON c.CDS_id = y.CDS_id JOIN genes x ON c.gene_id = x.gene_id ",whereStringFPKM,sep="")
	
	begin<-dbSendQuery(object@DB,"BEGIN;")
	
	#fetch records
	#genes
	genes.fpkm<-dbGetQuery(object@DB,geneFPKMQuery)
	genes.fpkm$sample_name<-factor(genes.fpkm$sample_name,levels=myLevels)
	genes.diff<-dbGetQuery(object@DB,geneDiffQuery)
	genes.diff$sample_1<-factor(genes.diff$sample_1,levels=myLevels)
	genes.diff$sample_2<-factor(genes.diff$sample_2,levels=myLevels)
	genes.repFpkm<-dbGetQuery(object@DB,geneRepFPKMQuery)
	genes.count<-dbGetQuery(object@DB,geneCountQuery)
	
	#isoforms
	isoform.fpkm<-dbGetQuery(object@DB,isoformFPKMQuery)
	isoform.fpkm$sample_name<-factor(isoform.fpkm$sample_name,levels=myLevels)
	isoform.diff<-dbGetQuery(object@DB,isoformDiffQuery)
	isoform.diff$sample_1<-factor(isoform.diff$sample_1,levels=myLevels)
	isoform.diff$sample_2<-factor(isoform.diff$sample_2,levels=myLevels)
	isoform.repFpkm<-dbGetQuery(object@DB,isoformRepFPKMQuery)
	isoform.count<-dbGetQuery(object@DB,isoformCountQuery)
	
	#CDS
	CDS.fpkm<-dbGetQuery(object@DB,CDSFPKMQuery)
	CDS.fpkm$sample_name<-factor(CDS.fpkm$sample_name,levels=myLevels)
	CDS.diff<-dbGetQuery(object@DB,CDSDiffQuery)
	CDS.diff$sample_1<-factor(CDS.diff$sample_1,levels=myLevels)
	CDS.diff$sample_2<-factor(CDS.diff$sample_2,levels=myLevels)
	CDS.repFpkm<-dbGetQuery(object@DB,CDSRepFPKMQuery)
	CDS.count<-dbGetQuery(object@DB,CDSCountQuery)
	
	#TSS
	TSS.fpkm<-dbGetQuery(object@DB,TSSFPKMQuery)
	TSS.fpkm$sample_name<-factor(TSS.fpkm$sample_name,levels=myLevels)
	TSS.diff<-dbGetQuery(object@DB,TSSDiffQuery)
	TSS.diff$sample_1<-factor(TSS.diff$sample_1,levels=myLevels)
	TSS.diff$sample_2<-factor(TSS.diff$sample_2,levels=myLevels)
	TSS.repFpkm<-dbGetQuery(object@DB,TSSRepFPKMQuery)
	TSS.count<-dbGetQuery(object@DB,TSSCountQuery)
	
	res<-new("CuffGene",
			id=geneId,
			annotation=dbGetQuery(object@DB,geneAnnotationQuery),
			fpkm=genes.fpkm,
			diff=genes.diff,
			repFpkm=genes.repFpkm,
			count=genes.count,
			isoforms=new("CuffFeature",
					annotation=dbGetQuery(object@DB,isoformAnnotationQuery),
					fpkm=isoform.fpkm,
					diff=isoform.diff,
					repFpkm=isoform.repFpkm,
					count=isoform.count
					),
			TSS=new("CuffFeature",
					annotation=dbGetQuery(object@DB,TSSAnnotationQuery),
					fpkm=TSS.fpkm,
					diff=TSS.diff,
					repFpkm=TSS.repFpkm,
					count=TSS.count
			),
			CDS=new("CuffFeature",
					annotation=dbGetQuery(object@DB,CDSAnnotationQuery),
					fpkm=CDS.fpkm,
					diff=CDS.diff,
					repFpkm=CDS.repFpkm,
					count=CDS.count
			)

			
		)
	end<-dbSendQuery(object@DB,"END;")
		
		res
}

setMethod("getGene",signature(object="CuffSet"),.getGene)
	
.getGenes<-function(object,geneIdList,sampleIdList=NULL){
	
	#Determine gene_id from geneIdList
	#This is useful so that people can pass, for example, isoform_id to geneIdList and getGenes will return full genes
	geneIdList<-getGeneId(object=object,geneIdList)
	
	#Sample subsetting
	if(!is.null(sampleIdList)){
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
	
	#ID Search String (SQL)
	idString<-'('
	for (i in geneIdList){
		idString<-paste(idString,"'",i,"',",sep="")
	}
	idString<-substr(idString,1,nchar(idString)-1)
	idString<-paste(idString,")",sep="")
	
	whereStringGene<-paste('WHERE (x.gene_id IN ',idString,' OR x.gene_short_name IN ',idString,')',sep="")
	whereStringGeneFPKM<-paste('WHERE (x.gene_id IN ',idString,' OR x.gene_short_name IN ',idString,')',sep="")
	whereStringGeneDiff<-paste('WHERE (x.gene_id IN ',idString,' OR x.gene_short_name IN ',idString,')',sep="")
	whereString<-paste('WHERE (x.gene_id IN ',idString,' OR g.gene_short_name IN ',idString,')',sep="")
	whereStringFPKM<-paste('WHERE (x.gene_id IN ',idString,' OR g.gene_short_name IN ',idString,')',sep="")
	whereStringDiff<-paste('WHERE (x.gene_id IN ',idString,' OR g.gene_short_name IN ',idString,')',sep="")
	
	if(!is.null(sampleIdList)){
		whereStringGene<-whereStringGene
		whereStringGeneFPKM<-paste(whereStringGeneFPKM,' AND y.sample_name IN ',sampleString,sep="")
		whereStringGeneDiff<-paste(whereStringGeneDiff,' AND (y.sample_1 IN ',sampleString,' AND y.sample_2 IN ',sampleString,')',sep="")
		whereString<-whereString
		whereStringFPKM<-paste(whereStringFPKM, ' AND y.sample_name IN ',sampleString,sep="")
		whereStringDiff<-paste(whereStringDiff,' AND (y.sample_1 IN ',sampleString,' AND y.sample_2 IN ',sampleString,')',sep="")
		
	}
	
	#dbQueries
	idQuery<-paste("SELECT DISTINCT gene_id from genes x ",whereStringGene,sep="")
	
	geneAnnotationQuery<-paste("SELECT * from genes x ", whereStringGene,sep="")
	geneFPKMQuery<-paste("SELECT y.* from genes x JOIN geneData y ON x.gene_id=y.gene_id ", whereStringGeneFPKM,sep="")
	geneDiffQuery<-paste("SELECT y.* from genes x JOIN geneExpDiffData y ON x.gene_id=y.gene_id ", whereStringGeneDiff,sep="")
	geneRepFPKMQuery<-paste("SELECT y.* from genes x JOIN geneReplicateData y on x.gene_id=y.gene_id ", whereStringGeneFPKM,sep="")
	geneCountQuery<-paste("SELECT y.* from genes x JOIN geneCount y on x.gene_id=y.gene_id ", whereStringGeneFPKM,sep="")
	
	isoformAnnotationQuery<-paste("SELECT x.* from isoforms x LEFT JOIN isoformFeatures xf ON x.isoform_id=xf.isoform_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	isoformFPKMQuery<-paste("SELECT y.* from isoforms x JOIN isoformData y ON x.isoform_id = y.isoform_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	isoformDiffQuery<-paste("SELECT y.* from isoforms x JOIN isoformExpDiffData y ON x.isoform_id = y.isoform_id JOIN genes g on x.gene_id=g.gene_id ", whereStringDiff,sep="")
	isoformRepFPKMQuery<-paste("SELECT y.* from isoforms x JOIN isoformReplicateData y on x.isoform_id=y.isoform_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	isoformCountQuery<-paste("SELECT y.* from isoforms x JOIN isoformCount y on x.isoform_id=y.isoform_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	
	TSSAnnotationQuery<-paste("SELECT x.* from TSS x LEFT JOIN TSSFeatures xf ON x.TSS_group_id=xf.TSS_group_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	TSSFPKMQuery<-paste("SELECT y.* from TSS x JOIN TSSData y ON x.TSS_group_id=y.TSS_group_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	TSSDiffQuery<-paste("SELECT y.* from TSS x JOIN TSSExpDiffData y ON x.TSS_group_id=y.TSS_group_id JOIN genes g on x.gene_id=g.gene_id ", whereStringDiff,sep="")
	TSSRepFPKMQuery<-paste("SELECT y.* from TSS x JOIN TSSReplicateData y on x.TSS_group_id=y.TSS_group_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	TSSCountQuery<-paste("SELECT y.* from TSS x JOIN TSSCount y on x.TSS_group_id=y.TSS_group_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	
	CDSAnnotationQuery<-paste("SELECT x.* from CDS x LEFT JOIN CDSFeatures xf ON x.CDS_id=xf.CDS_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	CDSFPKMQuery<-paste("SELECT y.* from CDS x JOIN CDSData y ON x.CDS_id = y.CDS_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	CDSDiffQuery<-paste("SELECT y.* from CDS x JOIN CDSExpDiffData y ON x.CDS_id = y.CDS_id JOIN genes g on x.gene_id=g.gene_id ", whereStringDiff,sep="")
	CDSRepFPKMQuery<-paste("SELECT y.* from CDS x JOIN CDSReplicateData y on x.CDS_id=y.CDS_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	CDSCountQuery<-paste("SELECT y.* from CDS x JOIN CDSCount y on x.CDS_id=y.CDS_id JOIN genes g on x.gene_id=g.gene_id ", whereStringFPKM,sep="")
	
	promotersDistQuery<-paste("SELECT x.* FROM promoterDiffData x LEFT JOIN genes g ON x.gene_id=g.gene_id ", whereString,sep="")
	splicingDistQuery<-paste("SELECT x.* FROM splicingDiffData x LEFT JOIN genes g ON x.gene_id=g.gene_id ", whereString,sep="")
	CDSDistQuery<-paste("SELECT x.* FROM CDSDiffData x LEFT JOIN genes g ON x.gene_id=g.gene_id ", whereString,sep="")
	
	begin<-dbSendQuery(object@DB,"BEGIN;")
	
	#fetch records
	#genes
	write("Getting gene information:",stderr())
	write("\tFPKM",stderr())
	genes.fpkm<-dbGetQuery(object@DB,geneFPKMQuery)
	genes.fpkm$sample_name<-factor(genes.fpkm$sample_name,levels=myLevels)
	write("\tDifferential Expression Data",stderr())
	genes.diff<-dbGetQuery(object@DB,geneDiffQuery)
	genes.diff$sample_1<-factor(genes.diff$sample_1,levels=myLevels)
	genes.diff$sample_2<-factor(genes.diff$sample_2,levels=myLevels)
	write("\tAnnotation Data",stderr())
	genes.annot<-dbGetQuery(object@DB,geneAnnotationQuery)
	write("\tReplicate FPKMs",stderr())
	genes.repFpkm<-dbGetQuery(object@DB,geneRepFPKMQuery)
	write("\tCounts",stderr())
	genes.count<-dbGetQuery(object@DB,geneCountQuery)
	
	#isoforms
	write("Getting isoforms information:",stderr())
	write("\tFPKM",stderr())
	isoform.fpkm<-dbGetQuery(object@DB,isoformFPKMQuery)
	isoform.fpkm$sample_name<-factor(isoform.fpkm$sample_name,levels=myLevels)
	write("\tDifferential Expression Data",stderr())
	isoform.diff<-dbGetQuery(object@DB,isoformDiffQuery)
	isoform.diff$sample_1<-factor(isoform.diff$sample_1,levels=myLevels)
	isoform.diff$sample_2<-factor(isoform.diff$sample_2,levels=myLevels)
	write("\tAnnotation Data",stderr())
	isoform.annot<-dbGetQuery(object@DB,isoformAnnotationQuery)
	write("\tReplicate FPKMs",stderr())
	isoform.repFpkm<-dbGetQuery(object@DB,isoformRepFPKMQuery)
	write("\tCounts",stderr())
	isoform.count<-dbGetQuery(object@DB,isoformCountQuery)
	
	#CDS
	write("Getting CDS information:",stderr())
	write("\tFPKM",stderr())
	CDS.fpkm<-dbGetQuery(object@DB,CDSFPKMQuery)
	CDS.fpkm$sample_name<-factor(CDS.fpkm$sample_name,levels=myLevels)
	write("\tDifferential Expression Data",stderr())
	CDS.diff<-dbGetQuery(object@DB,CDSDiffQuery)
	CDS.diff$sample_1<-factor(CDS.diff$sample_1,levels=myLevels)
	CDS.diff$sample_2<-factor(CDS.diff$sample_2,levels=myLevels)
	write("\tAnnotation Data",stderr())
	CDS.annot<-dbGetQuery(object@DB,CDSAnnotationQuery)
	write("\tReplicate FPKMs",stderr())
	CDS.repFpkm<-dbGetQuery(object@DB,CDSRepFPKMQuery)
	write("\tCounts",stderr())
	CDS.count<-dbGetQuery(object@DB,CDSCountQuery)
	
	#TSS
	write("Getting TSS information:",stderr())
	write("\tFPKM",stderr())
	TSS.fpkm<-dbGetQuery(object@DB,TSSFPKMQuery)
	TSS.fpkm$sample_name<-factor(TSS.fpkm$sample_name,levels=myLevels)
	write("\tDifferential Expression Data",stderr())
	TSS.diff<-dbGetQuery(object@DB,TSSDiffQuery)
	TSS.diff$sample_1<-factor(TSS.diff$sample_1,levels=myLevels)
	TSS.diff$sample_2<-factor(TSS.diff$sample_2,levels=myLevels)
	write("\tAnnotation Data",stderr())
	TSS.annot<-dbGetQuery(object@DB,TSSAnnotationQuery)
	write("\tReplicate FPKMs",stderr())
	TSS.repFpkm<-dbGetQuery(object@DB,TSSRepFPKMQuery)
	write("\tCounts",stderr())
	TSS.count<-dbGetQuery(object@DB,TSSCountQuery)
	
	#Promoters
	write("Getting promoter information:", stderr())
	write("\tdistData",stderr())
	promoters.distData<-dbGetQuery(object@DB,promotersDistQuery)
	promoters.distData$sample_1<-factor(promoters.distData$sample_1,levels=myLevels)
	promoters.distData$sample_2<-factor(promoters.distData$sample_2,levels=myLevels)
	
	#Splicing
	write("Getting splicing information:", stderr())
	write("\tdistData",stderr())
	splicing.distData<-dbGetQuery(object@DB,splicingDistQuery)
	splicing.distData$sample_1<-factor(splicing.distData$sample_1,levels=myLevels)
	splicing.distData$sample_2<-factor(splicing.distData$sample_2,levels=myLevels)

	#relCDS
	write("Getting relCDS information:", stderr())
	write("\tdistData",stderr())
	CDS.distData<-dbGetQuery(object@DB,CDSDistQuery)
	CDS.distData$sample_1<-factor(CDS.distData$sample_1,levels=myLevels)
	CDS.distData$sample_2<-factor(CDS.distData$sample_2,levels=myLevels)

	
	res<-new("CuffGeneSet",
			#TODO: Fix ids so that it only displays those genes in CuffGeneSet
			ids=as.character(dbGetQuery(object@DB,idQuery)),
			annotation=genes.annot,
			fpkm=genes.fpkm,
			diff=genes.diff,
			repFpkm=genes.repFpkm,
			count=genes.count,
			isoforms=new("CuffFeatureSet",
					annotation=isoform.annot,
					fpkm=isoform.fpkm,
					diff=isoform.diff,
					repFpkm=isoform.repFpkm,
					count=isoform.count
					),
			TSS=new("CuffFeatureSet",
					annotation=TSS.annot,
					fpkm=TSS.fpkm,
					diff=TSS.diff,
					repFpkm=TSS.repFpkm,
					count=TSS.count
					),
			CDS=new("CuffFeatureSet",
					annotation=CDS.annot,
					fpkm=CDS.fpkm,
					diff=CDS.diff,
					repFpkm=CDS.repFpkm,
					count=CDS.count
					),
			promoters=new("CuffFeatureSet",
					annotation=genes.annot,
					fpkm=genes.fpkm,
					diff=promoters.distData
					),
			splicing=new("CuffFeatureSet",
					annotation=TSS.annot,
					fpkm=TSS.fpkm,
					diff=splicing.distData
					),
			relCDS=new("CuffFeatureSet",
					annotation=genes.annot,
					fpkm=genes.fpkm,
					diff=CDS.distData
					)
			)
	end<-dbSendQuery(object@DB,"END;")		
	res
}

setMethod("getGenes",signature(object="CuffSet"),.getGenes)

.getGeneId<-function(object,idList){
	#Query that takes list of any identifier and retrieves gene_id values from db
	searchString<-"("
	for(i in idList){
		searchString<-paste(searchString,"'",i,"',",sep="")
	}
	searchString<-substr(searchString,1,nchar(searchString)-1)
	searchString<-paste(searchString,")",sep="")
	
	geneIdQuery<-paste("SELECT DISTINCT g.gene_id FROM genes g LEFT JOIN isoforms i on g.gene_id=i.gene_id LEFT JOIN TSS t on g.gene_id=t.gene_id LEFT JOIN CDS c ON g.gene_id=c.gene_id WHERE (g.gene_id IN ",searchString," OR g.gene_short_name IN ",searchString," OR i.isoform_id IN ",searchString," OR t.tss_group_id IN ",searchString," OR c.CDS_id IN ",searchString,")",sep="")
	#print(geneIdQuery)
	res<-dbGetQuery(object@DB,geneIdQuery)
	as.vector(res[,1])
}

setMethod("getGeneId",signature(object="CuffSet"),.getGeneId)

.getFeatures<-function(object,featureIdList,sampleIdList=NULL,level='isoforms'){
	#Sample subsetting
	if(!is.null(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	
	sampleString<-'('
	for (i in myLevels){
		sampleString<-paste(sampleString,"'",i,"',",sep="")
	}
	sampleString<-substr(sampleString,1,nchar(sampleString)-1)
	sampleString<-paste(sampleString,")",sep="")
	
	#ID Search String (SQL)
	idString<-'('
	for (i in featureIdList){
		idString<-paste(idString,"'",i,"',",sep="")
	}
	idString<-substr(idString,1,nchar(idString)-1)
	idString<-paste(idString,")",sep="")
	
	whereString<-paste(' WHERE (x.',slot(object,level)@idField,' IN ',idString,')',sep="")
	whereStringFPKM<-paste(' WHERE (x.',slot(object,level)@idField,' IN ',idString,')',sep="")
	whereStringDiff<-paste(' WHERE (x.',slot(object,level)@idField,' IN ',idString,')',sep="")
	
	if(!is.null(sampleIdList)){
		whereString<-whereString
		whereStringFPKM<-paste(whereStringFPKM, ' AND y.sample_name IN ',sampleString,sep="")
		whereStringDiff<-paste(whereStringDiff,' AND (y.sample_1 IN ',sampleString,' AND y.sample_2 IN ',sampleString,')',sep="")
	}
	
	
	AnnotationQuery<-paste("SELECT x.* from ",slot(object,level)@tables$mainTable," x LEFT JOIN ",slot(object,level)@tables$featureTable," xf ON x.",slot(object,level)@idField,"=xf.",slot(object,level)@idField, whereString,sep="")
	FPKMQuery<-paste("SELECT y.* from ",slot(object,level)@tables$mainTable," x JOIN ",slot(object,level)@tables$dataTable," y ON x.",slot(object,level)@idField," = y.",slot(object,level)@idField,whereStringFPKM,sep="")
	DiffQuery<-paste("SELECT y.* from ",slot(object,level)@tables$mainTable," x JOIN ",slot(object,level)@tables$expDiffTable," y ON x.",slot(object,level)@idField," = y.",slot(object,level)@idField,whereStringDiff,sep="")
	repFPKMQuery<-paste("SELECT y.* from ",slot(object,level)@tables$mainTable," x JOIN ",slot(object,level)@tables$replicateTable," y ON x.",slot(object,level)@idField," = y.",slot(object,level)@idField,whereStringFPKM,sep="")
	countQuery<-paste("SELECT y.* from ",slot(object,level)@tables$mainTable," x JOIN ",slot(object,level)@tables$countTable," y ON x.",slot(object,level)@idField," = y.",slot(object,level)@idField,whereStringFPKM,sep="")
	
	#print(AnnotationQuery)
	#print(FPKMQuery)
	#print(DiffQuery)
	#print(repFPKMQuery)
	#print(countQuery)
	
	begin<-dbSendQuery(object@DB,"BEGIN;")	
	res<-new("CuffFeatureSet",
			annotation=dbGetQuery(object@DB,AnnotationQuery),
			fpkm=dbGetQuery(object@DB,FPKMQuery),
			diff=dbGetQuery(object@DB,DiffQuery),
			repFpkm=dbGetQuery(object@DB,repFPKMQuery),
			count=dbGetQuery(object@DB,countQuery)
		)
	end<-dbSendQuery(object@DB,"END;")		
	res
	
}

setMethod("getFeatures",signature(object="CuffSet"),.getFeatures)
	



#getGeneIds from featureIds
#SELECT DISTINCT g.gene_id from genes g LEFT JOIN isoforms i on g.gene_id=i.gene_id LEFT JOIN TSS t on g.gene_id=t.gene_id LEFT JOIN CDS c on g.gene_id=c.gene_id WHERE (g.gene_id IN ('$VAL') OR i.isoform_id IN ('$VAL') OR t.tss_group_id IN ('$VAL') OR c.CDS_id IN ('$VAL') OR g.gene_short_name IN ('$VAL'));


#getSig() returns a list vectors of significant features by pairwise comparisons
#Depricated in favor of .getSig2
#.getSig<-function(object,x,y,level="genes",testTable=FALSE){
#	mySamp<-samples(slot(object,level))
#	sigGenes<-list()
#	if(level %in% c('promoters','splicing','relCDS')){
#		diffTable<-slot(object,level)@table
#	}else{
#		diffTable<-slot(object,level)@tables$expDiffTable
#	}
#	
#	#Restrict samples to those provided as x and y
#	if(!missing(x) && !missing(y)){
#		mySamp<-c(x,y)
#		if(!all(mySamp %in% samples(slot(object,level)))){
#			stop("One or more values of 'x' or 'y' are not valid sample names!")
#		}
#	}
#	
#	for (ihat in c(1:(length(mySamp)-1))){
#		for(jhat in c((ihat+1):length(mySamp))){
#			i<-mySamp[ihat]
#			j<-mySamp[jhat]
#			testName<-paste(i,j,sep="vs")
#			queryString<-paste("('",i,"','",j,"')",sep="")
#			sql<-paste("SELECT ",slot(object,level)@idField," from ", diffTable," WHERE sample_1 IN ",queryString," AND sample_2 IN ",queryString, " AND significant='yes'",sep="")
#			sig<-dbGetQuery(object@DB,sql)
#			sigGenes[[testName]]<-sig[,1]
#		}
#	}
#	#TODO: Add conditional return for if x & y are not null, to just return that test...
#	if(testTable){
#		tmp<-reshape2:::melt.list(sigGenes)
#		return(cast(tmp,value~...,length))
#	}else{
#		return(sigGenes)
#	}
#
#}

#Depricated in favor of .getSig
#.getSig2<-function(object,x,y,level="genes",testTable=FALSE,alpha=0.05){
#	mySamp<-samples(slot(object,level))
#	sigGenes<-list()
#	if(level %in% c('promoters','splicing','relCDS')){
#		diffTable<-slot(object,level)@table
#	}else{
#		diffTable<-slot(object,level)@tables$expDiffTable
#	}
#	
#	#Restrict samples to those provided as x and y
#	if(!missing(x) && !missing(y)){
#		mySamp<-c(x,y)
#		if(!all(mySamp %in% samples(slot(object,level)))){
#			stop("One or more values of 'x' or 'y' are not valid sample names!")
#		}
#	}
#	
#	for (ihat in c(1:(length(mySamp)-1))){
#		for(jhat in c((ihat+1):length(mySamp))){
#			i<-mySamp[ihat]
#			j<-mySamp[jhat]
#			testName<-paste(i,j,sep="vs")
#			queryString<-paste("('",i,"','",j,"')",sep="")
#			sql<-paste("SELECT ",slot(object,level)@idField,",p_value,q_value from ", diffTable," WHERE sample_1 IN ",queryString," AND sample_2 IN ",queryString, " AND STATUS='OK'",sep="")
#			sig<-dbGetQuery(object@DB,sql)
#			
#			#recalculate q-values for all tests in single pairwise comparison
#			if(!missing(x) && !(missing(y))) {
#				sig$q_value<-p.adjust(sig$p_value,method="BH")
#			}
#			#Filter on alpha
#			sig<-sig[sig$q_value<=alpha,]
#			sigGenes[[testName]]<-sig[,1]
#		}
#	}
#
#	if(testTable){
#		tmp<-reshape2:::melt.list(sigGenes)
#		return(cast(tmp,value~...,length))
#	}else{
#		return(sigGenes)
#	}
#	
#}

.getSig<-function(object,x,y,alpha=0.05,level='genes'){
	mySamp<-samples(slot(object,level))
	
	if(level %in% c('promoters','splicing','relCDS')){
		diffTable<-slot(object,level)@table
	}else{
		diffTable<-slot(object,level)@tables$expDiffTable
	}
	
	#Restrict samples to those provided as x and y
	if(!missing(x) && !missing(y)){
		mySamp<-c(x,y)
		if(!all(mySamp %in% samples(slot(object,level)))){
			stop("One or more values of 'x' or 'y' are not valid sample names!")
		}
	}
	
	queryString<-paste("(",paste(mySamp,collapse="','",sep=""),")",sep="'")
	sql<-paste("SELECT ",slot(object,level)@idField,",p_value from ", diffTable," WHERE sample_1 IN ",queryString," AND sample_2 IN ",queryString, " AND STATUS='OK'",sep="")
	#print(sql)
	sig<-dbGetQuery(object@DB,sql)
	sig$q_value<-p.adjust(sig$p_value,method="BH")
	sig<-sig[sig$q_value<=alpha,]
	sigGenes<-unique(sig[[slot(object,level)@idField]])
	sigGenes
}

setMethod("getSig",signature(object="CuffSet"),.getSig)


.getSigTable<-function(object,alpha=0.05,level='genes'){
	
	if(level %in% c('promoters','splicing','relCDS')){
		diffTable<-slot(object,level)@table
	}else{
		diffTable<-slot(object,level)@tables$expDiffTable
	}
	
	sql<-paste("SELECT ",slot(object,level)@idField,", sample_1, sample_2, p_value from ", diffTable," WHERE status='OK';")
	sig<-dbGetQuery(object@DB,sql)
	sig$q_value<-p.adjust(sig$p_value,method='BH')
	sig$testName<-paste(sig$sample_1,"vs",sig$sample_2,sep="")
	
	#filter on alpha and clean table
	sig$testResult<-0
	sig$testResult[sig$q_value<=alpha]<-1
	
	fieldsNeeded<-c('gene_id','testName','testResult')
	sig<-sig[,fieldsNeeded]
	
	#recast
	sig.table<-acast(sig,gene_id~testName,value='testResult')
	
	#remove genes that do not reject null in any test
	sig.table<-sig.table[rowSums(sig.table,na.rm=T)>0,]
	
	sig.table
}

setMethod("getSigTable",signature(object="CuffSet"),.getSigTable)

#Find similar genes
.findSimilar<-function(object,x,n,distThresh,returnGeneSet=TRUE,...){
	#x can be either a gene_id, gene_short_name or a vector of FPKM values (fake gene expression profile)
	#TODO: make findSimilar work with all levels
	if(is.character(x)){
		myGene<-getGene(object,x)
		sig<-makeprobsvec(fpkmMatrix(myGene,...)[1,])
	}else if(is.vector(x)){
		sig<-makeprobsvec(x)
	}
	allGenes<-fpkmMatrix(object@genes,...)
	allGenes<-t(makeprobs(t(allGenes)))
	compare<-function(q){
		JSdistVec(sig,q)
	}
	myDist<-apply(allGenes,MARGIN=1,compare)
	
	if(!missing(distThresh)){
		myDist<-myDist[myDist<=distThresh]
	}
	myDist<-sort(myDist)
	
	if(!missing(n)){
		myDist<-myDist[1:n]
	}
	
	mySimilarIds<-names(myDist)
	
	if(returnGeneSet){
		mySimilarGenes<-getGenes(object,mySimilarIds,...)
		return(mySimilarGenes)
	}else{
		res<-as.data.frame(myDist)
		colnames(res)<-c("distance")
		return(res)
	}
	
}
setMethod("findSimilar",signature(object="CuffSet"),.findSimilar)

############
#SQL access
############


################
#Misc Utilities
################
.getLevels<-function(object){
	levelsQuery<-'SELECT s.sample_name FROM samples s ORDER BY s.sample_index ASC'
	levels<-dbGetQuery(object@DB,levelsQuery)$sample_name
	levels
}

setMethod("getLevels",signature(object="CuffSet"),.getLevels)

.getRepLevels<-function(object){
	levelsQuery<-'SELECT r.rep_name FROM replicates r LEFT JOIN samples s ON r.sample_name=s.sample_name ORDER BY s.sample_index ASC'
	levels<-dbGetQuery(object@DB,levelsQuery)$rep_name
	levels
}

setMethod("getRepLevels",signature(object="CuffSet"),.getRepLevels)

.checkSamples<-function(dbConn,sampleIdList){
	dbSamples<-dbReadTable(dbConn,"samples")
	if (all(sampleIdList %in% dbSamples$sample_name)){
		return(TRUE)
	}else{
		return(FALSE)
	}
}


#####################
#Add FeatureData    #
#####################
.addFeatures<-function(object,features,level="genes",...){
	if(!is.data.frame(features)){
		stop("features must be a data.frame")
	}
	colnames(features)[1]<-slot(object,level)@idField
	colnames(features)<-make.db.names(object@DB,colnames(features),unique=T)
	dbWriteTable(object@DB,slot(object,level)@tables$featureTable,features,row.names=F,overwrite=T)
	indexQuery<-paste("CREATE INDEX ",slot(object,level)@idField," ON ", slot(object,level)@tables$featureTable," (",slot(object,level)@idField,")",sep="")
	res<-dbGetQuery(object@DB,indexQuery)
}

setMethod("addFeatures",signature(object="CuffSet"),.addFeatures)

#TODO: Add method to purge existing feature data table to allow 'refresh' of feature level data

##############
#Reporting
##############
#runReport<-function(){
#	if(!file.exists(".output")){
#		dir.create(".output")
#	}
#	file.copy(system.file("reports/runReport.Rnw", package="cummeRbund"),paste(".output/","runReport.Rnw",sep=""),overwrite=T)
#	myWD<-getwd()
#	setwd(".output")
#	Sweave("runReport.Rnw")
#	tools::texi2dvi("runReport.tex",pdf=TRUE)
#	setwd(myWD)
#}

