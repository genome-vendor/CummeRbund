# TODO: Add comment
# 
# Author: lgoff
###############################################################################


#####################
#File Archetype parsing
#####################

#RunInfo
loadRunInfo<-function(runInfoFile,
		dbConn,
		path,
		fileArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors=FALSE,
		row.names=NULL,
		...) {
	
	#Setup and reporting
	write(paste("Reading Run Info File ",runInfoFile,sep=""),stderr())
	fileArgs$file = runInfoFile
	
	#Read Run Info file
	runInfo = as.data.frame(do.call(read.table,fileArgs))
	
	#Parsing
	#not needed...
	
	#Load into database (runInfo table)
	write("Writing runInfo Table",stderr())
	insert_SQL<-'INSERT INTO runInfo VALUES(:param, :value)'
	bulk_insert(dbConn,insert_SQL,runInfo)
	
}

#ReplicateTable
loadRepTable<-function(repTableFile,
		dbConn,
		path,
		fileArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors=FALSE,
		row.names=NULL,
		...) {
		
	#Setup and reporting
	write(paste("Reading Read Group Info  ",repTableFile,sep=""),stderr())
	fileArgs$file = repTableFile
	
	#Read Run Info file
	full = as.data.frame(do.call(read.table,fileArgs))
	#print(head(full))
	
	#Fix sample_names
	full$condition<-make.db.names(dbConn,as.character(full$condition),unique=FALSE)
	
	#Parsing
	#For now, I need to concatenate condition and replicate number
	full$rep_name<-paste(full$condition,full$replicate_num,sep="_")
	
	#Load into database (replicates table)
	write("Writing replicates Table",stderr())
	insert_SQL<-'INSERT INTO replicates VALUES(:file, :condition, :replicate_num, :rep_name, :total_mass, :norm_mass, :internal_scale, :external_scale)'
	bulk_insert(dbConn,insert_SQL,full)
}

#Genes
loadGenes<-function(fpkmFile,
		diffFile,
		promoterFile,
		countFile,
		replicateFile,
		dbConn,
		path,
		#Arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		diffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		promoterArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		countArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		replicateArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE,
		row.names=NULL,
		...) {

	#Error Trapping
	if (missing(fpkmFile))
		stop("fpkmFile cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	#TODO test dbConn connection and database structure
	
	idCols = c(1:9)

	#Read primary file
	write(paste("Reading ",fpkmFile,sep=""),stderr())
	fpkmArgs$file = fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	
	########
	#Handle Sample Names
	########

	
	#Check that samples table is populated
	write("Checking samples table...",stderr())
	samples<-getSamplesFromColnames(full)
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	dbSamples<-dbReadTable(dbConn,"samples")
	if (dim(dbSamples)[1]>0) {
		if (all(samples %in% dbSamples$sample_name)){
			write ("OK!",stderr())
		}else{
			stop("Sample mismatch!")
		}
	}else{
		write("Populating samples table...",stderr())
		populateSampleTable(samples,dbConn)
	}
	
	######
	#Populate genes table
	######
	genesTable<-full[,c(1:3,5,7:9)]
	write("Writing genes table",stderr())
	#dbWriteTable(dbConn,'genes',genesTable,row.names=F,append=T)
	insert_SQL<-'INSERT INTO genes VALUES(:tracking_id, :class_code, :nearest_ref_id, :gene_short_name, :locus, :length, :coverage)'
	bulk_insert(dbConn,insert_SQL,genesTable)
	
	######
	#Populate geneData table
	######
	write("Reshaping geneData table",stderr())
	genemelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
	colnames(genemelt)[colnames(genemelt)=='variable']<-'sample_name'
	#Clean up and normalize data
	genemelt$measurement = ""
	
	genemelt$measurement[grepl("_FPKM$",genemelt$sample_name)] = "fpkm"
	genemelt$measurement[grepl("_conf_lo$",genemelt$sample_name)] = "conf_lo"
	genemelt$measurement[grepl("_conf_hi$",genemelt$sample_name)] = "conf_hi"
	genemelt$measurement[grepl("_status$",genemelt$sample_name)] = "status"

	genemelt$sample_name<-gsub("_FPKM$","",genemelt$sample_name)
	genemelt$sample_name<-gsub("_conf_lo$","",genemelt$sample_name)
	genemelt$sample_name<-gsub("_conf_hi$","",genemelt$sample_name)
	genemelt$sample_name<-gsub("_status$","",genemelt$sample_name)
	
	#Adjust sample names with make.db.names
	genemelt$sample_name <- make.db.names(dbConn,as.vector(genemelt$sample_name),unique=FALSE)
	
	#Recast
	write("Recasting",stderr())
	genemelt<-as.data.frame(dcast(genemelt,...~measurement))
	
	#debugging
	#write(colnames(genemelt),stderr())
	
	#Write geneData table
	write("Writing geneData table",stderr())
	#dbWriteTable(dbConn,'geneData',as.data.frame(genemelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)
	insert_SQL<-'INSERT INTO geneData VALUES(:tracking_id,:sample_name,:fpkm,:conf_hi,:conf_lo,:status)'
	bulk_insert(dbConn,insert_SQL,genemelt[,c(1:2,5,3,4,6)])
	
	#######
	#Handle gene_exp.diff
	#######
	
	if(file.exists(diffFile)){
		#Read diff file
		write(paste("Reading ",diffFile,sep=""),stderr())
		diffArgs$file = diffFile
		#Something like this to make sure sample names are treated as character values and not numeric, logical, etc.
		#diffArgs$colClasses<-c(rep('character',7),rep('numeric',6),'character')
		diff<-as.data.frame(do.call(read.table,diffArgs))
		if(dim(diff)[1]>0){
			#Adjust sample names with make.db.names
			diff$sample_1<-make.db.names(dbConn,as.vector(diff$sample_1),unique=FALSE)
			diff$sample_2<-make.db.names(dbConn,as.vector(diff$sample_2),unique=FALSE)
			
			write("Writing geneExpDiffData table",stderr())
			diffCols<-c(1,5:14)
			
			#debugging
			#write(colnames(diff[,diffCols]),stderr())
			
			#dbWriteTable(dbConn,'geneExpDiffData',diff[,diffCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO geneExpDiffData VALUES(:test_id,:sample_1,:sample_2,:status,:value_1,:value_2,?,:test_stat,:p_value,:q_value,:significant)"
			bulk_insert(dbConn,insert_SQL,diff[,diffCols])
		}else{
			write(paste("No records found in", diffFile),stderr())
		}
	
	}
	
	########
	#TODO: Handle promoters.diff
	########
	if(file.exists(promoterFile)){
		#Read promoterFile
		write(paste("Reading ",promoterFile,sep=""),stderr())
		promoterArgs$file = promoterFile
		promoter<-as.data.frame(do.call(read.table,promoterArgs))
		
		write("Writing promoterDiffData table",stderr())
		promoterCols<-c(2,5:14)
		if(dim(promoter)[1]>0){
			#dbWriteTable(dbConn,'promoterDiffData',promoter[,promoterCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO promoterDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,promoter[,promoterCols])
		}else{
			write(paste("No records found in", promoterFile),stderr())
		}
	}
	
	#########
	#Handle Feature Data (this will actually be done on CuffData objects instead...but I may include something here as well)
	#########
	
	###########
	#Handle Counts .count_tracking
	###########
	if(file.exists(countFile)){
		
		idCols = c(1)
		
		#Read countFile
		write(paste("Reading ", countFile,sep=""),stderr())
		countArgs$file = countFile
		counts<-as.data.frame(do.call(read.table,countArgs))
		
		if(dim(counts)[1]>0){
			#Reshape geneCount table
			write("Reshaping geneCount table",stderr())
			countmelt<-melt(counts,id.vars=c("tracking_id"),measure.vars=-idCols)
			colnames(countmelt)[colnames(countmelt)=='variable']<-'sample_name'
			
			countmelt$measurement = ""
			
			countmelt$measurement[grepl("_count$",countmelt$sample_name)] = "count"
			countmelt$measurement[grepl("_count_variance$",countmelt$sample_name)] = "variance"
			countmelt$measurement[grepl("_count_uncertainty_var$",countmelt$sample_name)] = "uncertainty"
			countmelt$measurement[grepl("_count_dispersion_var$",countmelt$sample_name)] = "dispersion"
			countmelt$measurement[grepl("_status$",countmelt$sample_name)] = "status"
			
			countmelt$sample_name<-gsub("_count$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_variance$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_uncertainty_var$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_dispersion_var$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_status$","",countmelt$sample_name)
			
			#Adjust sample names with make.db.names
			countmelt$sample_name <- make.db.names(dbConn,as.vector(countmelt$sample_name),unique=FALSE)
			
			#Recast
			write("Recasting",stderr())
			countmelt<-as.data.frame(dcast(countmelt,...~measurement))
			
			#debugging
			#write(colnames(countmelt),stderr())
			
	
			#Write geneCount table
			write("Writing geneCount table",stderr())
			insert_SQL<-'INSERT INTO geneCount VALUES(:tracking_id,:sample_name,:count,:variance,:uncertainty,:dispersion,:status)'
			bulk_insert(dbConn,insert_SQL,countmelt)
		}else{
			write(paste("No records found in", countFile),stderr())
		}
		
	}
		
		
	###########
	#Handle Replicates .rep_tracking
	###########
	if(file.exists(replicateFile)){

		idCols = 1
		#Read countFile
		write(paste("Reading read group info in ", replicateFile,sep=""),stderr())
		replicateArgs$file = replicateFile
		reps<-as.data.frame(do.call(read.table,replicateArgs))
		#print(head(reps))
		
		if(dim(reps)[1]>0){
		
			#Adjust sample names with make.db.names
			reps$condition <- make.db.names(dbConn,as.character(reps$condition),unique=FALSE)
		
			#Create unique rep name
			reps$rep_name<-paste(reps$condition,reps$replicate,sep="_")
			colnames(reps)[colnames(reps)=="condition"]<-"sample_name"
			
			#Write geneReplicateData table
			write("Writing geneReplicateData table",stderr())
			insert_SQL<-'INSERT INTO geneReplicateData VALUES(:tracking_id,:sample_name,:replicate,:rep_name,:raw_frags,:internal_scaled_frags,:external_scaled_frags,:FPKM,:effective_length,:status)'
			bulk_insert(dbConn,insert_SQL,reps)
		}else{
			write(paste("No records found in", replicateFile),stderr())
		}
		
	}
	
}
	
#Isoforms
loadIsoforms<-function(fpkmFile,
		diffFile,
		countFile,
		replicateFile,
		dbConn,
		path,
		#Arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		diffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		countArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		replicateArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE,
		row.names=NULL,
		...) {
	
	#Error Trapping
	if (missing(fpkmFile))
		stop("fpkmFile cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	#TODO test dbConn connection and database structure
	
	idCols = c(1:9)
	
	#Read primary file
	write(paste("Reading ",fpkmFile,sep=""),stderr())
	fpkmArgs$file = fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	
	########
	#Handle Sample Names
	########
	
	
	#Check that samples table is populated
	write("Checking samples table...",stderr())
	samples<-getSamplesFromColnames(full)
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	dbSamples<-dbReadTable(dbConn,"samples")
	if (dim(dbSamples)[1]>0) {
		if (all(samples %in% dbSamples$sample_name)){
			write ("OK!",stderr())
		}else{
			write(samples,stderr())
			stop("Sample mismatch!")
		}
	}else{
		write("Populating samples table...",stderr())
		populateSampleTable(samples,dbConn)
	}
	
	######
	#Populate isoforms table
	######
	isoformCols<-c(1,4,5,6,2,3,7:9)
	isoformsTable<-full[,isoformCols]
	
	#This is a temporary fix until p_id is added to the 'isoforms.fpkm_tracking' file
	isoformsTable<-cbind(isoformsTable[,1:2],data.frame(CDS_id=rep("NA",dim(isoformsTable)[1])),isoformsTable[,-c(1:2)])
	#print (head(isoformsTable))
	
	write("Writing isoforms table",stderr())
	#dbWriteTable(dbConn,'isoforms',as.data.frame(isoformsTable),row.names=F,append=T)
	insert_SQL<-'INSERT INTO isoforms VALUES(?,?,?,?,?,?,?,?,?,?)'
	bulk_insert(dbConn,insert_SQL,isoformsTable)
	
	######
	#Populate isoformData table
	######
	write("Reshaping isoformData table",stderr())
	isoformmelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
	colnames(isoformmelt)[colnames(isoformmelt)=='variable']<-'sample_name'
	#Clean up and normalize data
	isoformmelt$measurement = ""
	
	isoformmelt$measurement[grepl("_FPKM$",isoformmelt$sample_name)] = "fpkm"
	isoformmelt$measurement[grepl("_conf_lo$",isoformmelt$sample_name)] = "conf_lo"
	isoformmelt$measurement[grepl("_conf_hi$",isoformmelt$sample_name)] = "conf_hi"
	isoformmelt$measurement[grepl("_status$",isoformmelt$sample_name)] = "status"
	
	isoformmelt$sample_name<-gsub("_FPKM$","",isoformmelt$sample_name)
	isoformmelt$sample_name<-gsub("_conf_lo$","",isoformmelt$sample_name)
	isoformmelt$sample_name<-gsub("_conf_hi$","",isoformmelt$sample_name)
	isoformmelt$sample_name<-gsub("_status$","",isoformmelt$sample_name)
	
	#Adjust sample names with make.db.names
	isoformmelt$sample_name <- make.db.names(dbConn,as.vector(isoformmelt$sample_name),unique=FALSE)
	
	#Recast
	write("Recasting",stderr())
	isoformmelt<-as.data.frame(dcast(isoformmelt,...~measurement))
	
	#Write geneData table
	write("Writing isoformData table",stderr())
	#dbWriteTable(dbConn,'isoformData',as.data.frame(isoformmelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)
	insert_SQL<-"INSERT INTO isoformData VALUES(?,?,?,?,?,?)"
	bulk_insert(dbConn,insert_SQL,isoformmelt[,c(1:2,5,3,4,6)])
	
	#######
	#Handle isoform_exp.diff
	#######
	
	if(file.exists(diffFile)){
		#Read diff file
		write(paste("Reading ",diffFile,sep=""),stderr())
		diffArgs$file = diffFile
		diff<-as.data.frame(do.call(read.table,diffArgs))
		if(dim(diff)[1]>0){
			#Adjust sample names with make.db.names
			diff$sample_1<-make.db.names(dbConn,as.vector(diff$sample_1),unique=FALSE)
			diff$sample_2<-make.db.names(dbConn,as.vector(diff$sample_2),unique=FALSE)
		
			write("Writing isoformExpDiffData table",stderr())
			diffCols<-c(1,5:14)
			#dbWriteTable(dbConn,'isoformExpDiffData',diff[,diffCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO isoformExpDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,diff[,diffCols])
		}else{
			write(paste("No records found in",diffFile),stderr())
		}
	}
	
	###########
	#Handle Counts .count_tracking
	###########
	if(file.exists(countFile)){
		
		idCols = c(1)
		
		#Read countFile
		write(paste("Reading ", countFile,sep=""),stderr())
		countArgs$file = countFile
		counts<-as.data.frame(do.call(read.table,countArgs))
		
		if(dim(counts)[1]>0){
		
			#Reshape isoformCount table
			write("Reshaping isoformCount table",stderr())
			countmelt<-melt(counts,id.vars=c("tracking_id"),measure.vars=-idCols)
			colnames(countmelt)[colnames(countmelt)=='variable']<-'sample_name'
			
			countmelt$measurement = ""
			
			countmelt$measurement[grepl("_count$",countmelt$sample_name)] = "count"
			countmelt$measurement[grepl("_count_variance$",countmelt$sample_name)] = "variance"
			countmelt$measurement[grepl("_count_uncertainty_var$",countmelt$sample_name)] = "uncertainty"
			countmelt$measurement[grepl("_count_dispersion_var$",countmelt$sample_name)] = "dispersion"
			countmelt$measurement[grepl("_status$",countmelt$sample_name)] = "status"
			
			countmelt$sample_name<-gsub("_count$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_variance$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_uncertainty_var$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_dispersion_var$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_status$","",countmelt$sample_name)
			
			#Adjust sample names with make.db.names
			countmelt$sample_name <- make.db.names(dbConn,as.vector(countmelt$sample_name),unique=FALSE)
			
			
			#Recast
			write("Recasting",stderr())
			countmelt<-as.data.frame(dcast(countmelt,...~measurement))
			
			#debugging
			#write(colnames(countmelt),stderr())
			
			
			#Write isoformCount table
			write("Writing isoformCount table",stderr())
			insert_SQL<-'INSERT INTO isoformCount VALUES(:tracking_id,:sample_name,:count,:variance,:uncertainty,:dispersion,:status)'
			bulk_insert(dbConn,insert_SQL,countmelt)
		}else{
			write(paste("No records found in",countFile),stderr())
		}
	}
	
	
	###########
	#Handle Replicates .rep_tracking
	###########
	if(file.exists(replicateFile)){
		
		idCols = 1
		#Read countFile
		write(paste("Reading read group info in ", replicateFile,sep=""),stderr())
		replicateArgs$file = replicateFile
		reps<-as.data.frame(do.call(read.table,replicateArgs))
		#print(head(reps))
		
		if(dim(reps)[1]>0){
			
			#Adjust sample names with make.db.names
			reps$condition <- make.db.names(dbConn,as.character(reps$condition),unique=FALSE)
		
			#Create unique rep name
			reps$rep_name<-paste(reps$condition,reps$replicate,sep="_")
			colnames(reps)[colnames(reps)=="condition"]<-"sample_name"
			
			#Write isoformReplicateData table
			write("Writing isoformReplicateData table",stderr())
			insert_SQL<-'INSERT INTO isoformReplicateData VALUES(:tracking_id,:sample_name,:replicate,:rep_name,:raw_frags,:internal_scaled_frags,:external_scaled_frags,:FPKM,:effective_length,:status)'
			bulk_insert(dbConn,insert_SQL,reps)
		}else{
			write(paste("No records found in",replicateFile),stderr())
		}	
	}
	
}

#TSS groups
loadTSS<-function(fpkmFile,
		diffFile,
		splicingFile,
		countFile,
		replicateFile,
		dbConn,
		path,
		#Arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		diffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		splicingArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		countArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		replicateArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE,
		row.names=NULL,
		...) {
	
	#Error Trapping
	if (missing(fpkmFile))
		stop("fpkmFile cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	#TODO test dbConn connection and database structure
	
	idCols = c(1:9)
	
	#Read primary file
	write(paste("Reading ",fpkmFile,sep=""),stderr())
	fpkmArgs$file = fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	
	########
	#Handle Sample Names
	########
	
	
	#Check that samples table is populated
	write("Checking samples table...",stderr())
	samples<-getSamplesFromColnames(full)
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	dbSamples<-dbReadTable(dbConn,"samples")
	if (dim(dbSamples)[1]>0) {
		if (all(samples %in% dbSamples$sample_name)){
			write ("OK!",stderr())
		}else{
			stop("Sample mismatch!")
		}
	}else{
		write("Populating samples table...",stderr())
		populateSampleTable(samples,dbConn)
	}
	
	######
	#Populate TSS table
	######
	tssTable<-full[,c(1:5,7:9)]
	write("Writing TSS table",stderr())
	#dbWriteTable(dbConn,'TSS',tssTable,row.names=F,append=T)
	if (nrow(tssTable)>0){
		insert_SQL<-"INSERT INTO TSS VALUES(?,?,?,?,?,?,?,?)"
		bulk_insert(dbConn,insert_SQL,tssTable)
		
		######
		#Populate geneData table
		######
		write("Reshaping TSSData table",stderr())
		tssmelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
		colnames(tssmelt)[colnames(tssmelt)=='variable']<-'sample_name'
		#Clean up and normalize data
		tssmelt$measurement = ""
		
		tssmelt$measurement[grepl("_FPKM$",tssmelt$sample_name)] = "fpkm"
		tssmelt$measurement[grepl("_conf_lo$",tssmelt$sample_name)] = "conf_lo"
		tssmelt$measurement[grepl("_conf_hi$",tssmelt$sample_name)] = "conf_hi"
		tssmelt$measurement[grepl("_status$",tssmelt$sample_name)] = "status"
		
		tssmelt$sample_name<-gsub("_FPKM$","",tssmelt$sample_name)
		tssmelt$sample_name<-gsub("_conf_lo$","",tssmelt$sample_name)
		tssmelt$sample_name<-gsub("_conf_hi$","",tssmelt$sample_name)
		tssmelt$sample_name<-gsub("_status$","",tssmelt$sample_name)
		
		#Adjust sample names with make.db.names
		tssmelt$sample_name <- make.db.names(dbConn,as.vector(tssmelt$sample_name),unique=FALSE)
		
		#Recast
		write("Recasting",stderr())
		tssmelt<-as.data.frame(dcast(tssmelt,...~measurement))
		
		#Write geneData table
		write("Writing TSSData table",stderr())
		#dbWriteTable(dbConn,'TSSData',as.data.frame(tssmelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)

		insert_SQL<-"INSERT INTO TSSData VALUES(?,?,?,?,?,?)"
		bulk_insert(dbConn,insert_SQL,tssmelt[,c(1:2,5,3,4,6)])
	}else{
		write(paste("No records found in",fpkmFile),stderr())
		write("TSS FPKM tracking file was empty.",stderr())
	}
	#######
	#Handle tss_groups_exp.diff
	#######
	
	if(file.exists(diffFile)){
		#Read diff file
		write(paste("Reading ",diffFile,sep=""),stderr())
		diffArgs$file = diffFile
		diff<-as.data.frame(do.call(read.table,diffArgs))
		
		if(dim(diff)[1]>0){
			#Adjust sample names with make.db.names
			diff$sample_1<-make.db.names(dbConn,as.vector(diff$sample_1),unique=FALSE)
			diff$sample_2<-make.db.names(dbConn,as.vector(diff$sample_2),unique=FALSE)

			write("Writing TSSExpDiffData table",stderr())
			diffCols<-c(1,5:14)
			#dbWriteTable(dbConn,'TSSExpDiffData',diff[,diffCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO TSSExpDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,diff[,diffCols])
		}else{
			write(paste("No records found in",diffFile),stderr())
		}
	}
	
	#########
	#TODO: Handle splicing.diff
	########
	if(file.exists(splicingFile)){
		#Read promoterFile
		write(paste("Reading ",splicingFile,sep=""),stderr())
		splicingArgs$file = splicingFile
		splicing<-as.data.frame(do.call(read.table,splicingArgs))
		
		if(dim(splicing)[1]>0){
			write("Writing splicingDiffData table",stderr())
			splicingCols<-c(1:2,5:14)
			#dbWriteTable(dbConn,'splicingDiffData',splicing[,splicingCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO splicingDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,splicing[,splicingCols])
		}else{
			write(paste("No records found in",splicingFile),stderr())
		}
	}
	
	###########
	#Handle Counts .count_tracking
	###########
	if(file.exists(countFile)){
		
		idCols = c(1)
		
		#Read countFile
		write(paste("Reading ", countFile,sep=""),stderr())
		countArgs$file = countFile
		counts<-as.data.frame(do.call(read.table,countArgs))
		
		if(dim(counts)[1]>0){
		
			#Reshape TSSCount table
			write("Reshaping TSSCount table",stderr())
			countmelt<-melt(counts,id.vars=c("tracking_id"),measure.vars=-idCols)
			colnames(countmelt)[colnames(countmelt)=='variable']<-'sample_name'
			
			countmelt$measurement = ""
			
			countmelt$measurement[grepl("_count$",countmelt$sample_name)] = "count"
			countmelt$measurement[grepl("_count_variance$",countmelt$sample_name)] = "variance"
			countmelt$measurement[grepl("_count_uncertainty_var$",countmelt$sample_name)] = "uncertainty"
			countmelt$measurement[grepl("_count_dispersion_var$",countmelt$sample_name)] = "dispersion"
			countmelt$measurement[grepl("_status$",countmelt$sample_name)] = "status"
			
			countmelt$sample_name<-gsub("_count$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_variance$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_uncertainty_var$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_dispersion_var$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_status$","",countmelt$sample_name)
			
			#Adjust sample names with make.db.names
			countmelt$sample_name <- make.db.names(dbConn,as.vector(countmelt$sample_name),unique=FALSE)
			
			
			#Recast
			write("Recasting",stderr())
			countmelt<-as.data.frame(dcast(countmelt,...~measurement))
			
			#debugging
			#write(colnames(countmelt),stderr())
			
			
			#Write TSSCount table
			write("Writing TSSCount table",stderr())
			insert_SQL<-'INSERT INTO TSSCount VALUES(:tracking_id,:sample_name,:count,:variance,:uncertainty,:dispersion,:status)'
			bulk_insert(dbConn,insert_SQL,countmelt)
		}else{
			write(paste("No records found in",countFile),stderr())
		}
	}
	
	
	###########
	#Handle Replicates .rep_tracking
	###########
	if(file.exists(replicateFile)){
		
		idCols = 1
		#Read countFile
		write(paste("Reading read group info in ", replicateFile,sep=""),stderr())
		replicateArgs$file = replicateFile
		reps<-as.data.frame(do.call(read.table,replicateArgs))
		#print(head(reps))
		
		if(dim(reps)[1]>0){
				
			#Adjust sample names with make.db.names
			reps$condition <- make.db.names(dbConn,as.character(reps$condition),unique=FALSE)
		
			#Create unique rep name
			reps$rep_name<-paste(reps$condition,reps$replicate,sep="_")
			colnames(reps)[colnames(reps)=="condition"]<-"sample_name"
			
			#Write TSSReplicateData table
			write("Writing TSSReplicateData table",stderr())
			insert_SQL<-'INSERT INTO TSSReplicateData VALUES(:tracking_id,:sample_name,:replicate,:rep_name,:raw_frags,:internal_scaled_frags,:external_scaled_frags,:FPKM,:effective_length,:status)'
			bulk_insert(dbConn,insert_SQL,reps)
		}else{
			write(paste("No records found in",replicateFile),stderr())
		}
		
	}
	
}

#CDS
loadCDS<-function(fpkmFile,
		diffFile,
		CDSDiff,
		countFile,
		replicateFile,
		dbConn,
		path,
		#Arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		diffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		CDSDiffArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		countArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		replicateArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		sep="\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE,
		row.names=NULL,
		...) {
	
	#Error Trapping
	if (missing(fpkmFile))
		stop("fpkmFile cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	#TODO test dbConn connection and database structure
	
	idCols = c(1:9)
	
	#Read primary file
	write(paste("Reading ",fpkmFile,sep=""),stderr())
	fpkmArgs$file = fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	
	########
	#Handle Sample Names
	########
	
	
	
	#Check that samples table is populated
	write("Checking samples table...",stderr())
	samples<-getSamplesFromColnames(full)
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	dbSamples<-dbReadTable(dbConn,"samples")
	if (dim(dbSamples)[1]>0) {
		if (all(samples %in% dbSamples$sample_name)){
			write ("OK!",stderr())
		}else{
			stop("Sample mismatch!")
		}
	}else{
		write("Populating samples table...",stderr())
		populateSampleTable(samples,dbConn)
	}
	
	######
	#Populate CDS table
	######
	cdsTable<-full[,c(1:5,6:9)]
	write("Writing CDS table",stderr())
	#dbWriteTable(dbConn,'CDS',cdsTable,row.names=F,append=T)
	if (nrow(cdsTable)>0){
		insert_SQL<-"INSERT INTO CDS VALUES(?,?,?,?,?,?,?,?,?)"
		bulk_insert(dbConn,insert_SQL,cdsTable)
		
		######
		#Populate geneData table
		######
		write("Reshaping CDSData table",stderr())
		cdsmelt<-melt(full,id.vars=c("tracking_id"),measure.vars=-idCols,variable_name="sample_name")
		colnames(cdsmelt)[colnames(cdsmelt)=='variable']<-'sample_name'
		#Clean up and normalize data
		cdsmelt$measurement = ""
		
		cdsmelt$measurement[grepl("_FPKM$",cdsmelt$sample_name)] = "fpkm"
		cdsmelt$measurement[grepl("_conf_lo$",cdsmelt$sample_name)] = "conf_lo"
		cdsmelt$measurement[grepl("_conf_hi$",cdsmelt$sample_name)] = "conf_hi"
		cdsmelt$measurement[grepl("_status$",cdsmelt$sample_name)] = "status"
		
		cdsmelt$sample_name<-gsub("_FPKM$","",cdsmelt$sample_name)
		cdsmelt$sample_name<-gsub("_conf_lo$","",cdsmelt$sample_name)
		cdsmelt$sample_name<-gsub("_conf_hi$","",cdsmelt$sample_name)
		cdsmelt$sample_name<-gsub("_status$","",cdsmelt$sample_name)
		
		#Adjust sample names with make.db.names
		cdsmelt$sample_name <- make.db.names(dbConn,as.vector(cdsmelt$sample_name),unique=FALSE)
		
		#Recast
		write("Recasting",stderr())
		cdsmelt<-as.data.frame(dcast(cdsmelt,...~measurement))
		
		#Write geneData table
		write("Writing CDSData table",stderr())
		#dbWriteTable(dbConn,'CDSData',as.data.frame(cdsmelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)
		insert_SQL<-"INSERT INTO CDSData VALUES(?,?,?,?,?,?)"
		bulk_insert(dbConn,insert_SQL,cdsmelt[,c(1:2,5,3,4,6)])
	
	}else {
		write(paste("No records found in",fpkmFile),stderr())
		write("CDS FPKM tracking file was empty.",stderr())
	}
	
	
	#######
	#Handle cds_groups_exp.diff
	#######
	
	if(file.exists(diffFile)){
		#Read diff file
		write(paste("Reading ",diffFile,sep=""),stderr())
		diffArgs$file = diffFile
		diff<-as.data.frame(do.call(read.table,diffArgs))
		
		if(dim(diff)[1]>0){
			#Adjust sample names with make.db.names
			diff$sample_1<-make.db.names(dbConn,as.vector(diff$sample_1),unique=FALSE)
			diff$sample_2<-make.db.names(dbConn,as.vector(diff$sample_2),unique=FALSE)
			
			write("Writing CDSExpDiffData table",stderr())
			diffCols<-c(1,5:14)
			#dbWriteTable(dbConn,'CDSExpDiffData',diff[,diffCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO CDSExpDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,diff[,diffCols])
		}else{
			write(paste("No records found in",diffFile),stderr())
		}
	}
	
	#########
	#TODO: Handle CDS.diff
	########
	if(file.exists(CDSDiff)){
		#Read promoterFile
		write(paste("Reading ",CDSDiff,sep=""),stderr())
		CDSDiffArgs$file = CDSDiff
		CDS<-as.data.frame(do.call(read.table,CDSDiffArgs))
		if(dim(CDS)[1]>0){
			write("Writing CDSDiffData table",stderr())
			CDSCols<-c(2,5:14)
			#dbWriteTable(dbConn,'CDSDiffData',CDS[,CDSCols],row.names=F,append=T)
			insert_SQL<-"INSERT INTO CDSDiffData VALUES(?,?,?,?,?,?,?,?,?,?,?)"
			bulk_insert(dbConn,insert_SQL,CDS[,CDSCols])
		}else{
			write(paste("No records found in",CDSDiff),stderr())
		}
	}
	
	###########
	#Handle Counts .count_tracking
	###########
	if(file.exists(countFile)){
		
		idCols = c(1)
		
		#Read countFile
		write(paste("Reading ", countFile,sep=""),stderr())
		countArgs$file = countFile
		counts<-as.data.frame(do.call(read.table,countArgs))
		
		if(dim(counts)[1]>0){
		
			#Reshape CDSCount table
			write("Reshaping CDSCount table",stderr())
			countmelt<-melt(counts,id.vars=c("tracking_id"),measure.vars=-idCols)
			colnames(countmelt)[colnames(countmelt)=='variable']<-'sample_name'
			
			countmelt$measurement = ""
			
			countmelt$measurement[grepl("_count$",countmelt$sample_name)] = "count"
			countmelt$measurement[grepl("_count_variance$",countmelt$sample_name)] = "variance"
			countmelt$measurement[grepl("_count_uncertainty_var$",countmelt$sample_name)] = "uncertainty"
			countmelt$measurement[grepl("_count_dispersion_var$",countmelt$sample_name)] = "dispersion"
			countmelt$measurement[grepl("_status$",countmelt$sample_name)] = "status"
			
			countmelt$sample_name<-gsub("_count$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_variance$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_uncertainty_var$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_count_dispersion_var$","",countmelt$sample_name)
			countmelt$sample_name<-gsub("_status$","",countmelt$sample_name)
			
			#Adjust sample names with make.db.names
			countmelt$sample_name <- make.db.names(dbConn,as.vector(countmelt$sample_name),unique=FALSE)
			
			
			#Recast
			write("Recasting",stderr())
			countmelt<-as.data.frame(dcast(countmelt,...~measurement))
			
			#debugging
			#write(colnames(countmelt),stderr())
			
			
			#Write CDSCount table
			write("Writing CDSCount table",stderr())
			insert_SQL<-'INSERT INTO CDSCount VALUES(:tracking_id,:sample_name,:count,:variance,:uncertainty,:dispersion,:status)'
			bulk_insert(dbConn,insert_SQL,countmelt)
		}else{
			write(paste("No records found in",countFile),stderr())
		}
	}
	
	
	###########
	#Handle Replicates .rep_tracking
	###########
	if(file.exists(replicateFile)){
		
		idCols = 1
		#Read countFile
		write(paste("Reading read group info in ", replicateFile,sep=""),stderr())
		replicateArgs$file = replicateFile
		reps<-as.data.frame(do.call(read.table,replicateArgs))
		#print(head(reps))
		
		if(dim(reps)[1]>0){
				
			#Adjust sample names with make.db.names
			reps$condition <- make.db.names(dbConn,as.character(reps$condition),unique=FALSE)
		
			#Create unique rep name
			reps$rep_name<-paste(reps$condition,reps$replicate,sep="_")
			colnames(reps)[colnames(reps)=="condition"]<-"sample_name"
			
			#Write CDSReplicateData table
			write("Writing CDSReplicateData table",stderr())
			insert_SQL<-'INSERT INTO CDSReplicateData VALUES(:tracking_id,:sample_name,:replicate,:rep_name,:raw_frags,:internal_scaled_frags,:external_scaled_frags,:FPKM,:effective_length,:status)'
			bulk_insert(dbConn,insert_SQL,reps)
		}else{
			write(paste("No records found in",replicateFile),stderr())
		}
		
	}
	
}

########################
#Add FeatureData
########################


#####################
#Database Setup Functions
#####################

createDB_noIndex<-function(dbFname="cuffData.db",driver="SQLite") {
	#Builds sqlite db at 'dbFname' and returns a dbConnect object pointing to newly created database.
	#No indexes are present
	
	drv<-dbDriver(driver)
	db <- dbConnect(drv,dbname=dbFname)
	
	schema.text<-'
-- Creator:       MySQL Workbench 5.2.33/ExportSQLite plugin 2009.12.02
-- Author:        Loyal Goff
-- Caption:       New Model
-- Project:       Name of the project
-- Changed:       2012-04-30 22:21
-- Created:       2011-05-02 12:52
PRAGMA foreign_keys = OFF;

-- Schema: cuffData
BEGIN;
DROP TABLE IF EXISTS "genes";
CREATE TABLE "genes"(
  "gene_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT
);
DROP TABLE IF EXISTS "biasData";
CREATE TABLE "biasData"(
  "biasData_id" INTEGER PRIMARY KEY NOT NULL
);
DROP TABLE IF EXISTS "samples";
CREATE TABLE "samples"(
  "sample_index" INTEGER NOT NULL,
  "sample_name" VARCHAR(45) PRIMARY KEY NOT NULL
);
DROP TABLE IF EXISTS "TSS";
CREATE TABLE "TSS"(
  "TSS_group_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45) NOT NULL,
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_TSS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSData";
CREATE TABLE "TSSData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_TSSData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "CDS";
CREATE TABLE "CDS"(
  "CDS_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_CDS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_CDS_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
DROP TABLE IF EXISTS "CDSData";
CREATE TABLE "CDSData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_CDSData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "splicingDiffData";
CREATE TABLE "splicingDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_splicingDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_splicingDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSExpDiffData";
CREATE TABLE "TSSExpDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_TSSExpDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_TSSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "CDSDiffData";
CREATE TABLE "CDSDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "CDSExpDiffData";
CREATE TABLE "CDSExpDiffData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSExpDiffData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "promoterDiffData";
CREATE TABLE "promoterDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_promoterDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_promoterDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_promoterDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "geneFeatures";
CREATE TABLE "geneFeatures"(
  "gene_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_geneFeatures_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "TSSFeatures";
CREATE TABLE "TSSFeatures"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_TSSFeatures_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
DROP TABLE IF EXISTS "CDSFeatures";
CREATE TABLE "CDSFeatures"(
  "CDS_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_CDSFeatures_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id")
);
DROP TABLE IF EXISTS "model_transcripts";
CREATE TABLE "model_transcripts"(
  "model_transcript_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL
);
DROP TABLE IF EXISTS "geneCount";
CREATE TABLE "geneCount"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "uncertainty" FLOAT,
  "dispersion" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_geneCount_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "CDSCount";
CREATE TABLE "CDSCount"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "uncertainty" FLOAT,
  "dispersion" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_CDSCount_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "TSSCount";
CREATE TABLE "TSSCount"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "uncertainty" FLOAT,
  "dispersion" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_TSSCount_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "replicates";
CREATE TABLE "replicates"(
  "file" INTEGER NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" INT NOT NULL,
  "rep_name" VARCHAR(45) PRIMARY KEY NOT NULL,
  "total_mass" FLOAT,
  "norm_mass" FLOAT,
  "internal_scale" FLOAT,
  "external_scale" FLOAT,
  CONSTRAINT "fk_replicates_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "geneReplicateData";
CREATE TABLE "geneReplicateData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" INTEGER,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneData_genes10"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneReplicateData_replicates1"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_geneReplicateData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "CDSReplicateData";
CREATE TABLE "CDSReplicateData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" INTEGER,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates100"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_CDSReplicateData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSReplicateData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "TSSReplicateData";
CREATE TABLE "TSSReplicateData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" VARCHAR(45),
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates10000"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_TSSReplicateData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSReplicateData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "runInfo";
CREATE TABLE "runInfo"(
  "param" VARCHAR(45),
  "value" TEXT
);
DROP TABLE IF EXISTS "geneData";
CREATE TABLE "geneData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_geneData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "phenoData";
CREATE TABLE "phenoData"(
  "sample_name" VARCHAR(45) NOT NULL,
  "parameter" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45),
  CONSTRAINT "fk_phenoData_samples"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "geneExpDiffData";
CREATE TABLE "geneExpDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_geneExpDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_geneExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "isoforms";
CREATE TABLE "isoforms"(
  "isoform_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "gene_id" VARCHAR(45),
  "CDS_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_isoforms_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_isoforms_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_isoforms_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
DROP TABLE IF EXISTS "isoformData";
CREATE TABLE "isoformData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT NOT NULL,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_isoformData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
DROP TABLE IF EXISTS "isoformExpDiffData";
CREATE TABLE "isoformExpDiffData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "log2_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_isoformExpDiffData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "isoformFeatures";
CREATE TABLE "isoformFeatures"(
  "isoform_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_isoformFeatures_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
DROP TABLE IF EXISTS "features";
CREATE TABLE "features"(
--   GTF Features (all lines/records from reference .gtf file)
  "feature_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  "genes_gene_id" VARCHAR(45) NOT NULL,
  "isoforms_isoform_id" VARCHAR(45) NOT NULL,
  "seqname" VARCHAR(45) NOT NULL,
  "source" VARCHAR(45) NOT NULL,
  "type_id" INTEGER,
  "start" INTEGER,
  "end" INTEGER,
  "score" FLOAT,
  "strand" VARCHAR(45),
  "frame" VARCHAR(45),
  CONSTRAINT "fk_features_genes1"
    FOREIGN KEY("genes_gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_features_isoforms1"
    FOREIGN KEY("isoforms_isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
DROP TABLE IF EXISTS "attributes";
CREATE TABLE "attributes"(
  "attribute_lookup_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  "feature_id" INTEGER NOT NULL,
  "attribute" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_attribute_lookup_features1"
    FOREIGN KEY("feature_id")
    REFERENCES "features"("feature_id")
);
DROP TABLE IF EXISTS "isoformCount";
CREATE TABLE "isoformCount"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "uncertainty" FLOAT,
  "dispersion" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_isoformCount_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
DROP TABLE IF EXISTS "isoformReplicateData";
CREATE TABLE "isoformReplicateData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" INTEGER,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates10"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_isoformReplicateData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformReplicateData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
COMMIT;


			'
	create.sql <- strsplit(schema.text, "\n")[[1]]
	create.sql <- paste(collapse="\n", create.sql)
	create.sql <- strsplit(create.sql, ";")[[1]]
	create.sql <- create.sql[-length(create.sql)] #nothing to run here
	
	tmp <- sapply(create.sql,function(x) sqliteQuickSQL(db,x))
	db
}


createIndices<-function(dbFname="cuffData.db",driver="SQLite",verbose=F){
	
	drv<-dbDriver(driver)
	db <- dbConnect(drv,dbname=dbFname)
	
	index.text<-
'CREATE INDEX "genes.gsn_index" ON "genes"("gene_short_name");
CREATE INDEX "genes.cc_index" ON "genes"("class_code");
CREATE INDEX "TSS.fk_TSS_genes1" ON "TSS"("gene_id");
CREATE INDEX "TSSData.fk_TSSData_TSS1" ON "TSSData"("TSS_group_id");
CREATE INDEX "TSSData.fk_TSSData_samples1" ON "TSSData"("sample_name");
CREATE INDEX "CDS.fk_CDS_genes1" ON "CDS"("gene_id");
CREATE INDEX "CDS.fk_CDS_TSS1" ON "CDS"("TSS_group_id");
CREATE INDEX "CDSData.fk_CDSData_CDS1" ON "CDSData"("CDS_id");
CREATE INDEX "CDSData.fk_CDSData_samples1" ON "CDSData"("sample_name");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples1" ON "splicingDiffData"("sample_1");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples2" ON "splicingDiffData"("sample_2");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_TSS1" ON "splicingDiffData"("TSS_group_id");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_genes1" ON "splicingDiffData"("gene_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_TSS1" ON "TSSExpDiffData"("TSS_group_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples1" ON "TSSExpDiffData"("sample_1");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples2" ON "TSSExpDiffData"("sample_2");
CREATE INDEX "TSSExpDiffData.TSSExpDiffData_sig_index" ON "TSSExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples1" ON "CDSDiffData"("sample_1");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples2" ON "CDSDiffData"("sample_2");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_genes1" ON "CDSDiffData"("gene_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_CDS1" ON "CDSExpDiffData"("CDS_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples1" ON "CDSExpDiffData"("sample_1");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples2" ON "CDSExpDiffData"("sample_2");
CREATE INDEX "CDSExpDiffData.CDSExpDiffData_sig_index" ON "CDSExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_genes1" ON "promoterDiffData"("gene_id");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples1" ON "promoterDiffData"("sample_1");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples2" ON "promoterDiffData"("sample_2");
CREATE INDEX "geneFeatures.fk_geneFeatures_genes1" ON "geneFeatures"("gene_id");
CREATE INDEX "TSSFeatures.fk_TSSFeatures_TSS1" ON "TSSFeatures"("TSS_group_id");
CREATE INDEX "CDSFeatures.fk_CDSFeatures_CDS1" ON "CDSFeatures"("CDS_id");
CREATE INDEX "geneCount.fk_geneCount_samples1" ON "geneCount"("sample_name");
CREATE INDEX "geneCount.fk_geneCount_genes1" ON "geneCount"("gene_id");
CREATE INDEX "CDSCount.fk_CDSCount_CDS1" ON "CDSCount"("CDS_id");
CREATE INDEX "CDSCount.fk_CDSCount_samples1" ON "CDSCount"("sample_name");
CREATE INDEX "TSSCount.fk_TSSCount_TSS1" ON "TSSCount"("TSS_group_id");
CREATE INDEX "TSSCount.fk_TSSCount_samples1" ON "TSSCount"("sample_name");
CREATE INDEX "replicates.fk_replicates_samples1" ON "replicates"("sample_name");
CREATE INDEX "geneReplicateData.fk_geneReplicateData_genes1" ON "geneReplicateData"("gene_id");
CREATE INDEX "geneReplicateData.fk_geneReplicateData_replicates1" ON "geneReplicateData"("rep_name");
CREATE INDEX "geneReplicateData.fk_geneReplicateData_samples1" ON "geneReplicateData"("sample_name");
CREATE INDEX "CDSReplicateData.fk_CDSReplicateData_replicates1" ON "CDSReplicateData"("rep_name");
CREATE INDEX "CDSReplicateData.fk_CDSReplicateData_CDS1" ON "CDSReplicateData"("CDS_id");
CREATE INDEX "CDSReplicateData.fk_CDSReplicateData_samples1" ON "CDSReplicateData"("sample_name");
CREATE INDEX "TSSReplicateData.fk_TSSReplicateData_replicates1" ON "TSSReplicateData"("rep_name");
CREATE INDEX "TSSReplicateData.fk_TSSReplicateData_TSS1" ON "TSSReplicateData"("TSS_group_id");
CREATE INDEX "TSSReplicateData.fk_TSSReplicateData_samples1" ON "TSSReplicateData"("sample_name");
CREATE INDEX "geneData.fk_geneData_genes1" ON "geneData"("gene_id");
CREATE INDEX "geneData.fk_geneData_samples1" ON "geneData"("sample_name");
CREATE INDEX "phenoData.fk_phenoData_samples" ON "phenoData"("sample_name");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_genes1" ON "geneExpDiffData"("gene_id");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples1" ON "geneExpDiffData"("sample_1");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples2" ON "geneExpDiffData"("sample_2");
CREATE INDEX "geneExpDiffData.geneExpDiff_status_index" ON "geneExpDiffData"("status");
CREATE INDEX "geneExpDiffData.geneExpDiff_sig_index" ON "geneExpDiffData"("significant","p_value","q_value","test_stat");
CREATE INDEX "isoforms.fk_isoforms_TSS1" ON "isoforms"("TSS_group_id");
CREATE INDEX "isoforms.fk_isoforms_CDS1" ON "isoforms"("CDS_id");
CREATE INDEX "isoforms.fk_isoforms_genes1" ON "isoforms"("gene_id");
CREATE INDEX "isoformData.fk_isoformData_samples1" ON "isoformData"("sample_name");
CREATE INDEX "isoformData.fk_isoformData_isoforms1" ON "isoformData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_isoforms1" ON "isoformExpDiffData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples1" ON "isoformExpDiffData"("sample_1");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples2" ON "isoformExpDiffData"("sample_2");
CREATE INDEX "isoformExpDiffData.isoformExpDiffData_sig_index" ON "isoformExpDiffData"("test_stat","p_value","q_value","significant");
CREATE INDEX "isoformFeatures.fk_isoformFeatures_isoforms1" ON "isoformFeatures"("isoform_id");
CREATE INDEX "features.features_seqname_index" ON "features"("seqname");
CREATE INDEX "features.features_type_index" ON "features"("type_id");
CREATE INDEX "features.features_strand_index" ON "features"("strand");
CREATE INDEX "features.features_start_end_index" ON "features"("start","end");
CREATE INDEX "features.fk_features_genes1" ON "features"("genes_gene_id");
CREATE INDEX "features.fk_features_isoforms1" ON "features"("isoforms_isoform_id");
CREATE INDEX "attributes.fk_attributes_feature_id" ON "attributes"("feature_id");
CREATE INDEX "attributes.attributes_attribute_index" ON "attributes"("attribute");
CREATE INDEX "attributes.attributes_value_index" ON "attributes"("value");
CREATE INDEX "isoformCount.fk_isoformCount_isoforms1" ON "isoformCount"("isoform_id");
CREATE INDEX "isoformCount.fk_isoformCount_samples1" ON "isoformCount"("sample_name");
CREATE INDEX "isoformReplicateData.fk_isoformReplicateData_replicates1" ON "isoformReplicateData"("rep_name");
CREATE INDEX "isoformReplicateData.fk_isoformReplicateData_isoforms1" ON "isoformReplicateData"("isoform_id");
CREATE INDEX "isoformReplicateData.fk_isoformReplicateData_samples1" ON "isoformReplicateData"("sample_name");
'
	create.sql <- strsplit(index.text,"\n")[[1]]
	
	tmp <- sapply(create.sql,function(x){
			if (verbose){
						write(paste(x,sep=""),stderr())
					}
			sqliteQuickSQL(db,x)
	})
}


getSamples<-function(fpkmDF){
	sample_name<-unique(fpkmDF$sample)
	#sample_name<-as.data.frame(sample_name)
}

getSamplesFromColnames<-function(fpkmDF){
	samples<-gsub("_FPKM$","",colnames(fpkmDF)[grepl("_FPKM$",colnames(fpkmDF))])
}

populateSampleTable<-function(samples,dbConn){
	samples<-make.db.names(dbConn,samples,unique=FALSE)
	samples<-data.frame(index=c(1:length(samples)),sample_name=samples)
	dbWriteTable(dbConn,'samples',samples,row.names=F,append=T)
}

bulk_insert <- function(dbConn,sql,bound.data)
{
	dbBeginTransaction(dbConn)
	dbGetPreparedQuery(dbConn, sql, bind.data = bound.data)
	dbCommit(dbConn)
}

#############
#readCufflinks
#############
#TODO: Add count and replicate files
readCufflinks<-function(dir = getwd(),
						dbFile="cuffData.db",
						runInfoFile="run.info",
						repTableFile="read_groups.info",
						geneFPKM="genes.fpkm_tracking",
						geneDiff="gene_exp.diff",
						geneCount="genes.count_tracking",
						geneRep="genes.read_group_tracking",
						isoformFPKM="isoforms.fpkm_tracking",
						isoformDiff="isoform_exp.diff",
						isoformCount="isoforms.count_tracking",
						isoformRep="isoforms.read_group_tracking",
						TSSFPKM="tss_groups.fpkm_tracking",
						TSSDiff="tss_group_exp.diff",
						TSSCount="tss_groups.count_tracking",
						TSSRep="tss_groups.read_group_tracking",
						CDSFPKM="cds.fpkm_tracking",
						CDSExpDiff="cds_exp.diff",
						CDSCount="cds.count_tracking",
						CDSRep="cds.read_group_tracking",
						CDSDiff="cds.diff",
						promoterFile="promoters.diff",
						splicingFile="splicing.diff",
						driver = "SQLite",
						rebuild = FALSE,
						verbose = FALSE,
						...){
	
	#Set file locations with directory
	dbFile=file.path(dir,dbFile)
	runInfoFile=file.path(dir,runInfoFile)
	repTableFile=file.path(dir,repTableFile)
	geneFPKM=file.path(dir,geneFPKM)
	geneDiff=file.path(dir,geneDiff)
	geneCount=file.path(dir,geneCount)
	geneRep=file.path(dir,geneRep)
	isoformFPKM=file.path(dir,isoformFPKM)
	isoformDiff=file.path(dir,isoformDiff)
	isoformCount=file.path(dir,isoformCount)
	isoformRep=file.path(dir,isoformRep)
	TSSFPKM=file.path(dir,TSSFPKM)
	TSSDiff=file.path(dir,TSSDiff)
	TSSCount=file.path(dir,TSSCount)
	TSSRep=file.path(dir,TSSRep)
	CDSFPKM=file.path(dir,CDSFPKM)
	CDSExpDiff=file.path(dir,CDSExpDiff)
	CDSCount=file.path(dir,CDSCount)
	CDSRep=file.path(dir,CDSRep)
	CDSDiff=file.path(dir,CDSDiff)
	promoterFile=file.path(dir,promoterFile)
	splicingFile=file.path(dir,splicingFile)
					
					
	#Check to see whether dbFile exists
	if (!file.exists(dbFile) || rebuild == TRUE){
		#if not, create it
		write(paste("Creating database ",dbFile,sep=""),stderr())
		dbConn<-createDB_noIndex(dbFile)
		
		#populate DB
		if(file.exists(runInfoFile)){
			loadRunInfo(runInfoFile,dbConn)
		}
		
		if(file.exists(repTableFile)){
			loadRepTable(repTableFile,dbConn)
		}
		
		loadGenes(geneFPKM,geneDiff,promoterFile,countFile=geneCount,replicateFile=geneRep,dbConn)
		loadIsoforms(isoformFPKM,isoformDiff,isoformCount,isoformRep,dbConn)
		loadTSS(TSSFPKM,TSSDiff,splicingFile,TSSCount,TSSRep,dbConn)
		loadCDS(CDSFPKM,CDSExpDiff,CDSDiff,CDSCount,CDSRep,dbConn)
		
		#Create Indexes on DB
		write("Indexing Tables...",stderr())
		createIndices(dbFile,verbose=verbose)
		
		#load Distribution Tests
		#loadDistTests(promoterFile,splicingFile,CDSDiff)
		
	}
	dbConn<-dbConnect(dbDriver(driver),dbFile)
	return (
			new("CuffSet",DB = dbConn,
					#TODO: need to add replicate and count tables here and in AllClasses.R
					
					genes = new("CuffData",DB = dbConn, tables = list(mainTable = "genes",dataTable = "geneData",expDiffTable = "geneExpDiffData",featureTable = "geneFeatures",countTable="geneCount",replicateTable="geneReplicateData"), filters = list(),type = "genes",idField = "gene_id"),
					isoforms = new("CuffData", DB = dbConn, tables = list(mainTable = "isoforms",dataTable = "isoformData",expDiffTable = "isoformExpDiffData",featureTable = "isoformFeatures",countTable="isoformCount",replicateTable="isoformReplicateData"), filters = list(),type="isoforms",idField = "isoform_id"),
					TSS = new("CuffData", DB = dbConn, tables = list(mainTable = "TSS",dataTable = "TSSData",expDiffTable = "TSSExpDiffData",featureTable = "TSSFeatures",countTable="TSSCount",replicateTable="TSSReplicateData"), filters = list(),type = "TSS",idField = "TSS_group_id"),
					CDS = new("CuffData", DB = dbConn, tables = list(mainTable = "CDS",dataTable = "CDSData",expDiffTable = "CDSExpDiffData",featureTable = "CDSFeatures",countTable="CDSCount",replicateTable="CDSReplicateData"), filters = list(),type = "CDS",idField = "CDS_id"),
					promoters = new("CuffDist", DB = dbConn, table = "promoterDiffData",type="promoter",idField="gene_id"),
					splicing = new("CuffDist", DB = dbConn, table = "splicingDiffData",type="splicing",idField="TSS_group_id"),
					relCDS = new("CuffDist", DB = dbConn, table = "CDSDiffData",type="relCDS",idField="gene_id")
			)
	)	
							
}

############
# Handle GTF file
############
loadGTF<-function(gtfFile,dbConn) {
	
	#Error Trapping
	if (missing(gtfFile))
		stop("GTF file cannot be missing!")
	
	if (missing(dbConn))
		stop("Must provide a dbConn connection")
	
	gtf<-read.table(gtfFile,sep="\t",header=F)
	
	
	attributes<-melt(strsplit(as.character(gtf$V9),"; "))
	colnames(attributes)<-c("attribute","featureID")
	attributes<-paste(attributes$attribute,attributes$featureID)
	attributes<-strsplit(as.character(attributes)," ")
	attributes<-as.data.frame(do.call("rbind",attributes))
	
	colnames(attributes)<-c("attribute","value","featureID")
	attributes<-attributes[,c(3,1,2)]
	
	#Grab only gene_ID and transcript_ID to add to features table
	id.attributes<-attributes[attributes$attribute %in% c("gene_id","transcript_id"),]
	id.attributes$featureID<-as.numeric(as.character(id.attributes$featureID))
	id.attributes<-dcast(id.attributes,...~attribute)
	
	#Main features table
	features<-gtf[,c(1:8)]
	colnames(features)<-c("seqname","source","type","start","end","score","strand","frame")
	features$featureID<-as.numeric(as.character(rownames(features)))
	
	#Merge features and id.attributes
	features<-merge(features,id.attributes,by.x='featureID',by.y='featureID')
	features<-features[,c(1,10:11,2:9)]
	
	#strip gene_id and transcript_id from attributes
	attributes<-attributes[!(attributes$attribute %in% c("gene_id","transcript_id")),]
	
	#Write features table
	write("Writing features table",stderr())
	#dbWriteTable(dbConn,'geneData',as.data.frame(genemelt[,c(1:2,5,3,4,6)]),row.names=F,append=T)
	dbWriteTable(dbConn,'features',as.data.frame(features),append=F)
	
	#Write features table
	write("Writing feature attributes table",stderr())
	dbWriteTable(dbConn,'attributes',as.data.frame(attributes),append=F)
	
}
	

#######
#Unit Test
#######

#dbConn<-createDB()
#date()
#loadGenes("genes.fpkm_tracking","gene_exp.diff",dbConn)
#loadIsoforms("isoforms.fpkm_tracking","isoform_exp.diff",dbConn)
#loadTSS("tss_groups.fpkm_tracking","tss_group_exp.diff",dbConn)
#loadCDS("cds.fpkm_tracking","cds_exp.diff",dbConn)
#date()