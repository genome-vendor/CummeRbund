\name{readCufflinks}
\alias{readCufflinks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readCufflinks
}
\description{
This initializes the backend SQLite table and provides a DB connection for all downstream data analysis.
}
\usage{
readCufflinks(dir = getwd(), dbFile = "cuffData.db", 
				runInfoFile = "run.info", repTableFile = "read_groups.info",
				geneFPKM = "genes.fpkm_tracking", geneDiff = "gene_exp.diff", geneCount="genes.count_tracking", geneRep="genes.read_group_tracking",
				isoformFPKM = "isoforms.fpkm_tracking", isoformDiff = "isoform_exp.diff", isoformCount="isoforms.count_tracking", isoformRep="isoforms.read_group_tracking",
				TSSFPKM = "tss_groups.fpkm_tracking", TSSDiff = "tss_group_exp.diff", TSSCount="tss_groups.count_tracking", TSSRep="tss_groups.read_group_tracking",
				CDSFPKM = "cds.fpkm_tracking", CDSExpDiff = "cds_exp.diff", CDSCount="cds.count_tracking", CDSRep="cds.read_group_tracking",
				CDSDiff = "cds.diff",
				promoterFile = "promoters.diff", 
				splicingFile = "splicing.diff", 
				driver = "SQLite", rebuild = FALSE,verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
Directory in which all CuffDiff output files can be located. Defaults to current working directory.
}
  \item{dbFile}{
Name of backend database. Default is 'cuffData.db'
}
  \item{runInfoFile}{
run.info file
}
  \item{repTableFile}{
read_groups.info file
}
  \item{geneFPKM}{
genes.fpkm_tracking file
}
  \item{geneDiff}{
gene_exp.diff file
}
  \item{geneCount}{
genes.count_tracking file
}
  \item{geneRep}{
genes.read_group_tracking file
}
  \item{isoformFPKM}{
isoforms.fpkm_tracking file
}
  \item{isoformDiff}{
isoform_exp.diff file
}
  \item{isoformCount}{
isoforms.count_tracking file
}
  \item{isoformRep}{
isoform.read_group_tracking file
}
  \item{TSSFPKM}{
tss_groups.fpkm_tracking file
}
  \item{TSSDiff}{
tss_group_exp.diff file
}
  \item{TSSCount}{
tss_groups.count_tracking file
}
  \item{TSSRep}{
tss_groups.read_group_tracking file
}
  \item{CDSFPKM}{
cds.fpkm_tracking file
}
  \item{CDSExpDiff}{
cds_exp.diff file
}
  \item{CDSCount}{
cds.count_tracking file
}
  \item{CDSRep}{
cds.read_group_tracking file
}
  \item{CDSDiff}{
cds.diff file (distribution tests on CDS)
}
  \item{promoterFile}{
promoters.diff file (distribution tests on promoters)
}
  \item{splicingFile}{
splicing.diff (distribution tests on isoforms)
}
  \item{driver}{
Driver for backend database. (Currently only "SQLite" is supported).
}
  \item{rebuild}{
A logical argument to rebuild database backend.
}
  \item{verbose}{
A logical argument for super verbose reporting (As if it wasn't enough already!)
}
  \item{\dots}{
Additional arguments to readCufflinks
}
}
\details{
This is the initialization function for the cummeRbund package. It creates the SQLite backend database, populates the data, and provides a connection object for all future
interactions with the dataset. Once the initial build is complete, this function will default to using the database for all future sessions. 
IMPORTANT:  - Each R session should begin with a call to readCufflinks to re-open the connection to the database.
			- Should any connectivity issues to the database arise, another call to readCufflinks should create a new connection object and repair any issue.
			- The database can always be rebuild (using rebuild=TRUE) from the original CuffDiff output files.
}		
\value{
A CuffSet object.  A 'pointer' class that allows interaction with cufflinks/cuffdiff data via a SQLite database backend.
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None.
}


\examples{

	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data in sample directory and creates CuffSet object 'a'

}

