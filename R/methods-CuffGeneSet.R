########################
#methods-CuffGeneSet.R
#
#Author: Loyal A. Goff
#
#Date created: 5-17-2011
#
#Description: Defines a class of cufflinks data for multiple genes
#########################

#################
#Initialize		#
#################


#################
#Validate		#
#################


#################
#Class Methods	#
#################
setMethod("show","CuffGeneSet",function(object){
			cat(class(object),"instance for ", length(object), " genes\n",
					"\nSlots:\n\t annotation\n\t fpkm\n\t repFpkm\n\t diff\n\t count\n\t",
					"isoforms\t",class(object@isoforms),"instance of size",length(object@isoforms),"\n\t",
					"TSS\t\t",class(object@TSS),"instance of size",length(object@TSS),"\n\t",
					"CDS\t\t",class(object@CDS),"instance of size",length(object@CDS),"\n\t",
					"promoters\t\t",class(object@promoters),"instance of size",length(object@promoters),"\n\t",
					"splicing\t\t",class(object@splicing),"instance of size",length(object@splicing),"\n\t",
					"relCDS\t\t",class(object@relCDS),"instance of size",length(object@relCDS),"\n"

			)			
		}
)

#################
#Accessors
#################
#isoforms
setMethod("isoforms","CuffGeneSet",function(object){
			return(object@isoforms)	
		})
#TSS
setMethod("TSS","CuffGeneSet",function(object){
			return(object@TSS)
		})

#CDS
setMethod("CDS","CuffGeneSet",function(object){
			return(object@CDS)
		})

#promoters
setMethod("promoters","CuffGeneSet",function(object){
			return(object@promoters)
		})
#splicing
setMethod("splicing","CuffGeneSet",function(object){
			return(object@splicing)
		})
#relCDS
setMethod("relCDS","CuffGeneSet",function(object){
			return(object@relCDS)
		})

#################
#Subsetting		#
#################


#################
#Plotting		#
#################


#################
#Misc			#
#################