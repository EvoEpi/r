if (!require(topGO)) {
  install.packages("topGO")
  library(topGO)
}

#load 'universe' of go terms; a two column table: (1) gene model name and (2) a common separated list of GO terms
geneID2GO<-readMappings("<UNIVERSE>")
geneNames<-names(geneID2GO)

#load Genes Of Interest (GOI); a two column table: (1) gene model name and (2) a common separated list of GO terms
d<-read.table("<GOI>", header=F, sep="\t")
myInterestingGenes<-d$V1
geneList<-factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList)<-geneNames

#Biological Process (BP), Cellular Component (CC), or Molecular Function (MF) ontology enrichment
GOdata<-new("topGOdata", ontology="<BP or CC or MF>", allGenes=geneList, nodeSize=10, annot=annFUN.gene2GO, gene2GO=geneID2GO)
resultFisher<-runTest(GOdata, algorithm="classic", statistic="fisher")
resultWeight<-runTest(GOdata, algorithm="weight", statistic="fisher")
resultKS.elim<-runTest(GOdata, algorithm = "elim", statistic="ks")
allRes<-GenTable(GOdata, classic=resultFisher, weight=resultWeight, elimKS=resultKS.elim, orderBy="weight", ranksOf="weight", topNodes=100)
write.table(allRes, file="<OUTFILE>", sep="\t")
