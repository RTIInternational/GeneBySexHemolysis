#!/share/apps/R-3.5.1-a/bin/Rscript


args <- commandArgs(TRUE)

loop = TRUE
while (loop) {
        if (args[1] == "--up_table") {
                up_table = args[2]
        }

        if (args[1] == "--down_table") {
                down_table = args[2]
        }

        if (args[1] == "--out_plot") {
                out_plot = args[2]
        }
        if (length(args)>1) {
                args<-args[2:length(args)]
		} else{
			loop = FALSE
		}
}  

p_up <- read.table(up_table,header=T,stringsAsFactors=F)
p_up$p_log <- -1*log10(p_up$P)

p_down <- read.table(down_table,header=T,stringsAsFactors=F)
p_down$p_log <- log10(p_down$P)

p_all <- rbind(p_up,p_down)

inData=p_all
colChr = "CHR"
colPos = "POSITION"
colSE = "p_log"
manhattanCols = c(colChr,colPos,colSE)
manhattanData = inData[manhattanCols]
attach(manhattanData)
manhattanData = manhattanData[complete.cases(manhattanData),]
manhattanData = manhattanData[manhattanData[,colSE] != 0,]
manhattanData = manhattanData[order(manhattanData[,colChr], manhattanData[,colPos]),]
manhattanData[manhattanData[,colChr] %% 2 == 0, "color"] = "blue"
manhattanData[manhattanData[,colChr] %% 2 != 0, "color"] = "red"
manhattanData[manhattanData[,colChr] %% 2 == 0, "color"] = "grey30"
manhattanData[manhattanData[,colChr] %% 2 != 0, "color"] = "grey"
manhattanData[,"pch"] = 1
snpNum = table(manhattanData[,colChr])
snpStart = cumsum(snpNum)
snpStart = c(1,snpStart[-length(snpStart)] + 1)
snpEnd = snpStart + snpNum - 1
chrLen = (as.numeric(manhattanData[snpEnd,colPos]) - as.numeric(manhattanData[snpStart,colPos]))
chrStart = cumsum(chrLen)
chrStart = c(1, chrStart[-length(chrStart)]+1)
lobs = manhattanData[,colSE]
xpos = chrStart[as.numeric(as.factor(manhattanData[,colChr]))] + as.numeric(manhattanData[,colPos])- as.numeric(manhattanData[snpStart[as.numeric(as.factor(manhattanData[,colChr]))],colPos])

manhattanPlot = function(fileOut, hasHighlightList = FALSE, nonHighlightIndices = NULL, highlightIndices = NULL) {
  png(fileOut, width = 850,height=380, type="cairo")
  par(oma=c(0,0,0,0),
      mar=c(5,3,1,1), 
      ps=20,cex=1,cex.main=1,cex.lab=1,cex.axis=1,cex.sub=1)
  ylim1=max(abs(min(lobs)),abs(max(lobs)))
  plot(c(min(xpos),max(xpos)),c(min(lobs),max(lobs)),type='n',
       ylab="",xlab='Chromosome Position',xaxt='n',
       ylim=c(-1*ylim1-1,ylim1+1))
  points(xpos[nonHighlightIndices],
         lobs[nonHighlightIndices],
         col=manhattanData[nonHighlightIndices,"color"],
         cex=0.5,pch=manhattanData[nonHighlightIndices,"pch"])
  abline(h=0)
  abline(h=7.3, lty=3)
  abline(h=-7.3, lty=3)
  plotNum=c(1:23)
  loc = (chrStart+chrLen/2)[plotNum]
  mtext(c(as.numeric(names(snpNum))[1:22],"X")[plotNum],side=1,at=loc,line=0.5,par(las=3))
  dev.off()
}
nonHighlightIndices = rep(TRUE,nrow(manhattanData))
manhattanPlot(out_plot, nonHighlightIndices=nonHighlightIndices)

detach(manhattanData)


