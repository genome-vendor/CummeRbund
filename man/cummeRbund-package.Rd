\name{cummeRbund-package}
\alias{cummeRbund-package}
\alias{cummeRbund}
\docType{package}
\title{
cummeRbund: The finishing touch on your Tuxedo workflow. Analysis, manipulation, and visualization of Cufflinks HTS data.
~~ package title ~~
}
\description{
Allows for persistent storage, access, and manipulation of Cufflinks high-throughput sequencing data.  In addition, provides numerous plotting functions for commonly used visualizations.
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab cummeRbund\cr
Version: \tab 0.1.3\cr
Suggests: \tab \cr
Depends: \tab R (>= 2.7.0), RSQLite, reshape2, ggplot2, methods\cr
License: \tab MIT License\cr
Collate: \tab AllGenerics.R
AllClasses.R
database-setup.R
methods-CuffSet.R
methods-CuffData.R
methods-CuffDist.R
methods-CuffGeneSet.R
methods-CuffFeatureSet.R
methods-CuffGene.R
methods-CuffFeature.R
tools.R\cr
LazyLoad: \tab yes\cr
biocViews: \tab HighThroughputSequencing, HighThroughputSequencingData, RNAseq, RNAseqData, GeneExpression, DifferentialExpression, Infrastructure, DataImport, DataRepresentation, Visualization, Bioinformatics, Clustering, MultipleComparisons, QualityControl\cr
Packaged: \tab 2011-08-05 18:03:50 UTC; lgoff\cr
Built: \tab R 2.12.1; ; 2011-08-05 18:03:57 UTC; unix\cr
}

Index:
\preformatted{
CuffData-class          Class "CuffData"
CuffDist-class          Class "CuffDist"
CuffFeature-class       Class "CuffFeature"
CuffFeatureSet-class    Class "CuffFeatureSet"
CuffGene-class          Class "CuffGene"
CuffGeneSet-class       Class "CuffGeneSet"
CuffSet-class           Class "CuffSet"
JSdist                  Jensen-Shannon distance on columns
addFeatures             addFeatures
addFeatures-methods     ~~ Methods for Function addFeatures in Package
                        'cummeRbund' ~~
createDB                createDB
csBoxplot               csBoxplot
csBoxplot-methods       ~~ Methods for Function csBoxplot in Package
                        'cummeRbund' ~~
csDensity               Density plot of CuffData
csDensity-methods       ~~ Methods for Function csDensity in Package
                        'cummeRbund' ~~
csHeatmap               csHeatmap
csHeatmap-methods       ~~ Methods for Function csHeatmap in Package
                        'cummeRbund' ~~
csScatter               Scatter Plot
csScatter-methods       ~~ Methods for Function csScatter in Package
                        'cummeRbund' ~~
csVolcano               Volcano Plot
csVolcano-methods       ~~ Methods for Function csVolcano in Package
                        'cummeRbund' ~~
cummeRbund-package      cummeRbund: The finishing touch on your Tuxedo
                        workflow. Analysis, manipulation, and
                        visualization of Cufflinks HTS data.
diffData                Differential comparison data
diffData-methods        ~~ Methods for Function diffData in Package
                        'cummeRbund' ~~
dim-methods             ~~ Methods for Function dim in Package 'base'
                        ~~
expressionBarplot       Barplot
expressionBarplot-methods
                        ~~ Methods for Function expressionBarplot in
                        Package 'cummeRbund' ~~
expressionPlot          Expression Plot
expressionPlot-methods
                        ~~ Methods for Function expressionPlot in
                        Package 'cummeRbund' ~~
featureNames            Feature names
featureNames-methods    ~~ Methods for Function featureNames in Package
                        'cummeRbund' ~~
features                Features
features-methods        ~~ Methods for Function features in Package
                        'cummeRbund' ~~
fpkm                    Retrieve FPKM values
fpkm-methods            ~~ Methods for Function fpkm in Package
                        'cummeRbund' ~~
fpkmMatrix              Retrieve FPKM values as matrix
fpkmMatrix-methods      ~~ Methods for Function fpkmMatrix in Package
                        'cummeRbund' ~~
getGene                 getGene
getGene-methods         ~~ Methods for Function getGene in Package
                        'cummeRbund' ~~
getGenes                getGenes
getGenes-methods        ~~ Methods for Function getGenes in Package
                        'cummeRbund' ~~
getLevels               getLevels
getLevels-methods       ~~ Methods for Function getLevels in Package
                        'cummeRbund' ~~
length-methods          ~~ Methods for Function length in Package
                        'base' ~~
makeprobs               Transform a matrix into probabilities by
                        columns
readCufflinks           readCufflinks
samples                 Get sample list from CuffData object
samples-methods         ~~ Methods for Function samples in Package
                        'cummeRbund' ~~
shannon.entropy         Shannon entropy
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{cummeRbund-manual} \tab An R package for visualization and analysis of Cufflinks high-throughput sequencing data (source, pdf)\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
L. Goff, C. Trapnell

Maintainer: Loyal A. Goff <lgoff@csail.mit.edu>

}
\references{
~~ Literature or other references for background information ~~
}


\keyword{ package }

