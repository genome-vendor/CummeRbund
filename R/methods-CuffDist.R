#methods-CuffDist.R
#
#Author: Loyal A. Goff
#
#
####################

##################
#Initialize
##################
setMethod("initialize","CuffDist",
			function(.Object,
					DB,
					table="",
					type = c("promoter","splicing","relCDS"),
					idField = c("gene_id","tss_group_id"),
					... ){
				.Object<-callNextMethod(.Object,
						DB = DB,
						table = table,
						type = type,
						idField = idField,
						...)				
		}
)

setValidity("CuffDist",function(object){
		TRUE
		}
)			

################
#Class Methods
################
setMethod("show","CuffDist",
		function(object){
			size<-dim(object)
			cat(class(object), "instance with:\n\t",size[1]," ",object@type," records\n")
		}
)

setMethod("dim","CuffDist",
		function(x){
			countQuery<-paste("SELECT COUNT(",x@idField,") as n FROM ",x@table)
			nIds<-dbGetQuery(x@DB,countQuery)
			c(nIds$n)
		}
)

###################
#Accessors
###################
.values<-function(object){
	valueQuery<-paste("SELECT * FROM ",object@table,sep="")
	dbGetQuery(object@DB, valueQuery)
}

setMethod("distValues","CuffDist",.values)

setMethod("DB","CuffDist",function(object){
		return(object@DB)
		})

#setMethod("table","CuffDist",function(object){
#		return(object@table)
#		})

setMethod("type","CuffDist",function(object){
		return(object@type)
		})

setMethod("idField","CuffDist",function(object){
		return(object@idField)
		})

.samples<-function(object){
	res<-dbReadTable(object@DB,'samples')
	res<-res$sample_name
	res
}

setMethod("samples","CuffDist",.samples)

##################
#Setters
##################


##################
#Subsetting
##################


##################
#Plotting
##################

