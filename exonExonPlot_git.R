#!/usr/bin/Rscript

##Input is bed12 file generated from bam file that was filtered to only contain reads with junctions as identified by STAR with the jI flag/SAM attribute

##THIS IS SPECIFICALLY TAILORED TO APP
library(ggplot2)

#****
bedFile <- "/path/to/file.bed12"
dataBed <- data.frame(read.delim(bedFile, header=F, sep='\t'))
colnames(dataBed) <- c("Chr", "Start", "End", "ReadName", "BEDscore", "strand", "thickStart", "thickEnd", "RGB", "#blocks", "blockSizes", "blockStarts")
#***
exonFile <- "/path/to/exonsAPP.bed"
exonData <- data.frame(read.delim(exonFile, header=F, sep='\t'))
colnames(exonData) <- c("Chr", "Start", "End", "exonNum", "RefSeq", "newStart", "newEnd", "col", "label")
#Filter out reads in bed12 file that fall within the APP gene
appData <- dataBed[which(dataBed$Chr =="chr21" & ((dataBed$Start >= 25880550 & dataBed$Start <= 26170820) | (dataBed$End >= 25880550 & dataBed$End <= 26170820))),]

junctionData <- c()
junctionData3_1 <- c()
junctionData3_2 <- c()
junctionData3_3 <- c()
junctionData3_4 <- c()
junctionData3_5 <- c()

for (i in 1:nrow(appData)){
  numJ = appData[i, c("#blocks")]
  if (numJ=="2"){
    segA_start <- appData[i,c("Start")]
    segA_end <- appData[i,c("Start")] + strtoi(unlist(strsplit(toString(appData[i,c("blockSizes")]), ","))[1])
    segB_start <- appData[i,c("Start")] + strtoi(unlist(strsplit(toString(appData[i,c("blockStarts")]), ","))[2])
    segB_end <- appData[i,c("Start")] + strtoi(unlist(strsplit(toString(appData[i,c("blockStarts")]), ","))[2]) + strtoi(unlist(strsplit(toString(appData[i,c("blockSizes")]), ","))[2])
    exonA=0
    exonB=0
    exonTotal=0
    for (j in 1:nrow(exonData)){
      exonStart <- exonData[j, c("Start")]
      exonEnd <- exonData[j, c("End")]
      exonNum <- exonData[j, c("exonNum")]
      if ((segA_start >= exonStart & segA_start <= exonEnd) | (segA_end >= exonStart & segA_end <= exonEnd)){
        exonA <- exonNum
        exonTotal=exonTotal+1
      }
      if ((segB_start >= exonStart & segB_start <= exonEnd) | (segB_end >= exonStart & segB_end <= exonEnd)){
        exonB <- exonNum
        exonTotal=exonTotal+1
      }
    }
    if (exonTotal > 1 & exonA != exonB){
      dfName <- paste(paste("exon", exonA, sep=""), paste("exon", exonB, sep=""), sep="_") #exon#_exon#
      if (!dfName %in% junctionData){
        junctionData <- c(junctionData, dfName)
        do.call("<-", list(dfName, appData[i,]))
      }
      else {
        do.call("<-", list(as.symbol(dfName), rbind(eval(as.symbol(dfName)), appData[i,])))
      }
    }
  }
  if (numJ=="3"){
    segA_start <- appData[i,c("Start")]
    segA_end <- appData[i,c("Start")] + strtoi(unlist(strsplit(toString(appData[i,c("blockSizes")]), ","))[1])
    segB_start <- appData[i,c("Start")] + strtoi(unlist(strsplit(toString(appData[i,c("blockStarts")]), ","))[2])
    segB_end <- appData[i,c("Start")] + strtoi(unlist(strsplit(toString(appData[i,c("blockStarts")]), ","))[2]) + strtoi(unlist(strsplit(toString(appData[i,c("blockSizes")]), ","))[2])
    segC_start <- appData[i,c("Start")] + strtoi(unlist(strsplit(toString(appData[i,c("blockStarts")]), ","))[3])
    segC_end <- segC_start <- appData[i,c("Start")] + strtoi(unlist(strsplit(toString(appData[i,c("blockStarts")]), ","))[3]) + strtoi(unlist(strsplit(toString(appData[i,c("blockSizes")]), ","))[3])
    exonA=0
    exonB=0
    exonC=0
    exonTotal=0
    for (j in 1:nrow(exonData)){
      exonStart <- exonData[j, c("Start")]
      exonEnd <- exonData[j, c("End")]
      exonNum <- exonData[j, c("exonNum")]
      if ((segA_start >= exonStart & segA_start <= exonEnd) | (segA_end >= exonStart & segA_end <= exonEnd)){
        exonA <- exonNum
        exonTotal=exonTotal+1
      }
      if ((segB_start >= exonStart & segB_start <= exonEnd) | (segB_end >= exonStart & segB_end <= exonEnd)){
        exonB <- exonNum
        exonTotal=exonTotal+1
      }
      if ((segC_start >= exonStart & segC_start <= exonEnd) | (segC_end >= exonStart & segC_end <= exonEnd)){
        exonC <- exonNum
        exonTotal=exonTotal+1
      }
    }
    if (exonTotal == 3 & exonA != exonB & exonA != exonC & exonB != exonC){ #case1
      dfName <- paste(paste(paste("exon", exonA, sep=""), paste("exon", exonB, sep=""), sep="_"), paste("exon", exonC, sep=""), sep="_") #exon#_exon#_exon#
      if (!dfName %in% junctionData3_1){
        junctionData3_1 <- c(junctionData3_1, dfName)
        do.call("<-", list(dfName, appData[i,]))
      }
      else {
        do.call("<-", list(as.symbol(dfName), rbind(eval(as.symbol(dfName)), appData[i,])))
      }
    }
    if (exonTotal==3 & exonA==exonB & exonA!=exonC){ #case2
      dfName <- paste(paste("exon", exonB, sep=""), paste("exon", exonC, sep=""), sep="_") #exon#_exon#
      if (!dfName %in% junctionData3_2){
        junctionData3_2 <- c(junctionData3_2, dfName)
        do.call("<-", list(dfName, appData[i,]))
      }
      else {
        do.call("<-", list(as.symbol(dfName), rbind(eval(as.symbol(dfName)), appData[i,])))
      }
    }
    if (exonTotal==3 & exonB==exonC & exonA!=exonB){ #case3
      dfName <- paste(paste("exon", exonA, sep=""), paste("exon", exonB, sep=""), sep="_") #exon#_exon#
      if (!dfName %in% junctionData3_3){
        junctionData3_3 <- c(junctionData3_3, dfName)
        do.call("<-", list(dfName, appData[i,]))
      }
      else {
        do.call("<-", list(as.symbol(dfName), rbind(eval(as.symbol(dfName)), appData[i,])))
      }
    }
    if (exonTotal==2 & exonA!=exonB & exonA!=0 & exonB!=0){ #case4
      dfName <- paste(paste("exon", exonA, sep=""), paste("exon", exonB, sep=""), sep="_") #exon#_exon#
      if (!dfName %in% junctionData3_4){
        junctionData3_4 <- c(junctionData3_4, dfName)
        do.call("<-", list(dfName, appData[i,]))
      }
      else {
        do.call("<-", list(as.symbol(dfName), rbind(eval(as.symbol(dfName)), appData[i,])))
      }
    }
    if (exonTotal==2 & exonB!=exonC & exonB!=0 & exonC!=0){ #case5
      dfName <- paste(paste("exon", exonB, sep=""), paste("exon", exonC, sep=""), sep="_") #exon#_exon#
      if (!dfName %in% junctionData3_5){
        junctionData3_5 <- c(junctionData3_5, dfName)
        do.call("<-", list(dfName, appData[i,]))
      }
      else {
        do.call("<-", list(as.symbol(dfName), rbind(eval(as.symbol(dfName)), appData[i,])))
      }
    }
  }
}

##Adjustment factors
adjust <- c(list(exon18=25881669, exon17=25891618, exon16=25897321, exon15=25904670, exon14=25911332, exon13=25953958, exon12=25954894, exon11=25974207, exon10=25974931, exon9=25981245, exon7=25998781, exon6=26020437, exon5=26049393, exon4=26051434, exon3=26088027, exon2=26109932, exon1=26168348))

##Adjusting the coordinates to graph
for (k in 1:length(junctionData)){
  newExonA <- unlist(strsplit(junctionData[k], "_"))[1]
  Anum <- strtoi(unlist(strsplit(newExonA, ""))[length(unlist(strsplit(newExonA, "")))])
  newExonB <- unlist(strsplit(junctionData[k], "_"))[2]
  adjustA <- adjust[[newExonA]]
  adjustB <- adjust[[newExonB]]
  newStart <- c(eval(as.symbol(junctionData[k]))$Start - adjustA)
  newEnd <- c(get(junctionData[k])$End - adjustB)
  if ((Anum %% 2) == 0){
    rowNum <- c((seq.int(nrow(get(junctionData[k])))*10))
  }
  else {
    rowNum <- c((seq.int(nrow(get(junctionData[k])))*10)-5)  
  }
  do.call("<-", list(as.symbol(junctionData[k]), cbind(eval(as.symbol(junctionData[k])), newStart, newEnd, rowNum)))
}

for (m in 1:length(junctionData3_1)){
  newExonA <- unlist(strsplit(junctionData3_1[m], "_"))[1]
  Anum <- strtoi(unlist(strsplit(newExonA, ""))[length(unlist(strsplit(newExonA, "")))])
  newExonC <- unlist(strsplit(junctionData3_1[m], "_"))[3]
  adjustA <- adjust[[newExonA]]
  adjustC <- adjust[[newExonC]]
  newStart <- c(eval(as.symbol(junctionData3_1[m]))$Start - adjustA)
  newEnd <- c(get(junctionData3_1[m])$End - adjustC)
  if ((Anum %% 2) == 0){
    rowNum <- c((seq.int(nrow(get(junctionData3_1[m])))*-10)-30)
  }
  else {
    rowNum <- c((seq.int(nrow(get(junctionData3_1[m])))*-10)-25)  
  }
  do.call("<-", list(as.symbol(junctionData3_1[m]), cbind(eval(as.symbol(junctionData3_1[m])), newStart, newEnd, rowNum)))
}

##Combine all dataframes into one
graphData <- data.frame()
for (l in 1:length(junctionData)){
  graphData <- rbind(graphData, eval(as.symbol(junctionData[l])))
}
graph3Data <- data.frame()
for(n in 1:length(junctionData3_1)){
  graph3Data <- rbind(graph3Data, eval(as.symbol(junctionData3_1[n])))
}


##GRAPH
colorValues <- c(A="#8752a2", B="#e1058b", C="#71cfeb", D="#2d5980", E="#00864b", F="#fef200", G="#fdb813", H="#ee3f76", I="#5656ab", J="#2d5980", K="#70bf54", L="#f7941d", M="#da1f3e", N="#6a489d", O="#008dd1", P="#134c69", Q="#505050",Asnv="#0AD328", Csnv="#1B4BE1", Tsnv="#CD1919", Gsnv="#F2C621")
#***
pngFile <- "/path/to/outfile.png"

plot <- ggplot(data=exonData, aes(x=Start, y=End))+ xlim(-2300, 0) + ylim(-200, 300) + geom_segment(data = exonData, aes(x = -1*newStart, y = -10, xend = -1*newEnd, yend = -10, colour=col), size=30) + geom_text(data=exonData, aes(x=(-1*newStart-newEnd)/2, y=-25, label=label, fontface="bold"), size=12) + geom_segment(data=graphData, aes(x=-1*newStart, y=rowNum, xend=-1*newEnd, yend=rowNum), color="#B0B0B0", size=7) + geom_segment(data=graph3Data, aes(x=-1*newStart, y=rowNum, xend=-1*newEnd, yend=rowNum), color="#B0B0B0", size=7)+ geom_segment(data=snvs, aes(x=-1*newPos, y=rowNum, xend=-1*(newPos+1), yend=rowNum, colour=SNV), size=7)+scale_color_manual(values=colorValues)+ theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
png(pngFile, width=4000, height=2000)
print(plot)
dev.off()
