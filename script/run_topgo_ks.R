suppressWarnings(library(tcltk))
suppressWarnings(suppressMessages(library(topGO)))
suppressPackageStartupMessages(library(Rgraphviz))
suppressPackageStartupMessages(library(grid))

Args           <- commandArgs()
geneListFile   <- Args[4]
mappingFile    <- Args[5]
ontology       <- Args[6]
algorithm      <- Args[7]
outputPrefix   <- Args[8]
inputSigLevel  <- as.numeric( Args[9] )
outputSigLevel <- as.numeric( Args[10] )

# Gene selection function
topDiffGenes <- function(allScore) {
    return(allScore < inputSigLevel)
}

# Get genes along with p values
geneList     <- read.table(geneListFile)
genes        <- geneList[,2]
names(genes) <- geneList[,1]

# Get GO annotation
geneID2GO <- readMappings(file = mappingFile)

# Make topGOdata object
# nodeSize=10 : prune GO hierarchy of terms associated with < 10 genes
GOdata <- suppressMessages(new("topGOdata", ontology=ontology, allGenes=genes,
    geneSel=topDiffGenes, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=10))

# Run topGO
resultKS <- suppressMessages(runTest(GOdata, algorithm=algorithm,
    statistic="ks"))
nodecount <- length(score(resultKS))
allRes <- GenTable(GOdata, resultKS, topNodes=nodecount)
colNames <- names(allRes)
colNames[6] <- 'pval'
names(allRes) <- colNames
# Horrible way to get all the genes associated with each term
allRes$Genes <- sapply(allRes$GO.ID,
    function(x) gsub('[c()" \n]', '', genesInTerm(GOdata, x)))
sigRes <- allRes[suppressWarnings(as.numeric(allRes$pval)) < outputSigLevel,]

# Write results
write.table( allRes, file=paste0(outputPrefix, ".all.tsv"), quote=FALSE,
    row.names=FALSE, sep="\t" )

# Write PDF
pdf(paste0(outputPrefix, ".pdf"))
try(suppressWarnings(showSigOfNodes(GOdata, score(resultKS),
    firstSigNodes=5, useInfo="all")), silent=TRUE)
try(suppressWarnings(showSigOfNodes(GOdata, score(resultKS),
    firstSigNodes=10, useInfo="all")), silent=TRUE)
try(suppressWarnings(showSigOfNodes(GOdata, score(resultKS),
    firstSigNodes=nrow(sigRes), useInfo="all")), silent=TRUE)
try(suppressWarnings(lapply(sigRes[,1],
    function(x) showGroupDensity(GOdata, x))), silent=TRUE)
dev.off()
