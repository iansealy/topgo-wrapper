library(topGO)

Args           <- commandArgs()
geneListFile   <- Args[4]
mappingFile    <- Args[5]
ontology       <- Args[6]
outputFile     <- Args[7]
graphPdfPrefix <- Args[8]
sigLevel       <- as.numeric( Args[9] )

# Gene selection function
topDiffGenes <- function(allScore) {
    return(allScore < sigLevel)
}

# Get genes along with p values
geneList     <- read.table(geneListFile)
genes        <- geneList[,2]
names(genes) <- geneList[,1]

# Get GO annotation
geneID2GO <- readMappings(file = mappingFile)

# Make topGOdata object
# nodeSize=10 : prune GO hierarchy of terms associated with < 10 genes
GOdata <- new("topGOdata", ontology=ontology, allGenes=genes,
    geneSel=topDiffGenes, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=10)

# Run topGO
resultKS.elim <- runTest(GOdata, algorithm="elim", statistic="ks")
nodecount <- length(score(resultKS.elim))
allRes <- GenTable(GOdata, elimKS=resultKS.elim, topNodes=nodecount)

# Write results
write.table( allRes, file=outputFile, quote=FALSE, row.names=FALSE, sep="\t" )
printGraph(GOdata, resultKS.elim, firstSigNodes=5, fn.prefix=graphPdfPrefix,
    useInfo="all", pdfSW=TRUE)
