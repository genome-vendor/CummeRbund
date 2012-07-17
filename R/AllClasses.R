# AllClasses.R
# 
# Author: lgoff
###############################################################################

#TODO: I get the distinct feeling that these two should be nested environments, but I don't really know what that means.

#CuffRun class describes information about a Cuffdiff run?
#setClass("CuffRun",
#		representation(DB = "SQLiteConnection",
#						sampleInfo = "data.frame",
#						runParams = "data.frame"
#						)
#		)

#CuffData Class is a 'pointer' container to a group of features in a cufflinks dataset
setClass("CuffData",
		representation(DB = "SQLiteConnection",
						tables = "list",
						filters = "list",
						type = "character",
						idField = "character"
						)
		)

#CuffDist Class is a 'pointer' container to one of the distribution test data sets (promoters.diff, splicing.diff, cds.diff)
setClass("CuffDist",
		representation(DB = "SQLiteConnection",
						table = "character",
						type = "character",
						idField = "character"
						)
		)
		
#CuffSet Class is a 'pointer' container to a group of CuffData elements in a cufflinks dataset
setClass("CuffSet",
		representation(DB = "SQLiteConnection",
						runInfo = "data.frame",
						phenoData = "data.frame",
						conditions = "data.frame",
						genes = "CuffData",
						isoforms = "CuffData",
						TSS = "CuffData",
						CDS = "CuffData",
						promoters = "CuffDist",
						splicing = "CuffDist",
						relCDS = "CuffDist"
						)

)

#CuffFeature is a 'data' container for all information linked to a single 'idField' (cufflinks class agnostic)
setClass("CuffFeature",
		representation(annotation="data.frame",
						fpkm="data.frame",
						diff="data.frame",
						repFpkm="data.frame",
						count="data.frame"
				)
		)

#CuffGene is a 'data' container for all information linked to a single 'gene_id'
setClass("CuffGene",
		representation(id = "character",
						isoforms = "CuffFeature",
						TSS = "CuffFeature",
						CDS = "CuffFeature",
						promoters = "CuffFeature",
						splicing = "CuffFeature",
						relCDS = "CuffFeature"),
		contains="CuffFeature"
)


#CuffFeatureSet is a 'data' container for all information from a set of features
#This allows for plotting of gene set information
setClass("CuffFeatureSet",
		representation(annotation="data.frame",
				fpkm="data.frame",
				diff="data.frame",
				repFpkm="data.frame",
				count="data.frame"
			)
)

#CuffGene is a 'data' container for all information from a set of genes
#This allows for plotting of gene set information
setClass("CuffGeneSet",
		representation(ids = "character",
				isoforms = "CuffFeatureSet",
				TSS = "CuffFeatureSet",
				CDS= "CuffFeatureSet",
				promoters= "CuffFeatureSet",
				splicing= "CuffFeatureSet",
				relCDS= "CuffFeatureSet"),
		contains = "CuffFeatureSet"
)
