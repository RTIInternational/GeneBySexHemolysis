rm(list=ls())

#RBC Omics

csvdir <- 'Y:\\REDSIII\\data_files\\csv\ files\\'
linkDat <- read.csv(paste0(csvdir, 'linked-v6.csv'), stringsAsFactors = FALSE)

setwd("Y:/GWA_G_by_sex/Osmotic_SNP_check")

all<-list()
for(chr in c(12,17,2,20,7)){
  cauc <- read.table(paste0("raw/rbc.CAUCASIAN.chr",chr,".raw"),header=T,stringsAsFactors = F)
  afr <- read.table(paste0("raw/rbc.AFRAMRCN.chr",chr,".raw"),header=T,stringsAsFactors = F)
  asian <- read.table(paste0("raw/rbc.ASIAN.chr",chr,".raw"),header=T,stringsAsFactors = F)
  hisp <- read.table(paste0("raw/rbc.HISPANIC.chr",chr,".raw"),header=T,stringsAsFactors = F)
  other <- read.table(paste0("raw/rbc.OTHER.chr",chr,".raw"),header=T,stringsAsFactors = F)
  
  
  sum(colnames(cauc)==colnames(asian))
  sum(colnames(cauc)==colnames(afr))
  sum(colnames(cauc)==colnames(hisp))
  sum(colnames(cauc)==colnames(hisp))
  sum(colnames(afr)==colnames(asian))
  sum(colnames(afr)==colnames(hisp))
  for(i in 7:dim(cauc)[2]){
    if(colnames(afr)[i] != colnames(cauc)[i]){
      afr[,i]=abs(afr[,i]-2)
      colnames(afr)[i]=colnames(cauc)[i]
    }
    if(colnames(asian)[i] != colnames(cauc)[i]){
      asian[,i]=abs(asian[,i]-2)
      colnames(asian)[i]=colnames(cauc)[i]
    }
    if(colnames(hisp)[i] != colnames(cauc)[i]){
      hisp[,i]=abs(hisp[,i]-2)
      colnames(hisp)[i]=colnames(cauc)[i]
    }
    if(colnames(other)[i] != colnames(cauc)[i]){
      other[,i]=abs(other[,i]-2)
      colnames(other)[i]=colnames(cauc)[i]
    }
  }
  all[[paste0("chr",chr)]]<-rbind(cauc,afr,asian,hisp,other)
}

all_merged <- Reduce(merge,all)

out_var <- c("BSI.Subj.ID","HUB","Gender","numDonWithinTwoYear",
             "SCRN.pink_pct_hemol.adj",
             "SCRN.Oxidative_pct_hemolysis.adj",
             "SCRN.storage_pct_hemol.adj",
             "WBC", "RBC","HGB", "HCT", "MCV","RDW","PLT")


all_hemolysis <- merge(all_merged,linkDat[,out_var],by.x="IID",by.y="BSI.Subj.ID",all.x=T)

# genetically defined ethnicity groups
evdir="Y:/data"
easian_ev <- read.table(paste0(evdir,"/East.Asian.ex_related.EVs.txt"),header=T,stringsAsFactors = F)
sasian_ev <- read.table(paste0(evdir,"/South.Asian.ex_related.EVs.txt"),header=T,stringsAsFactors = F)
cauc_ev <- read.table(paste0(evdir,"/Caucasian.ex_related.EVs.txt"),header=T,stringsAsFactors = F)
afr_ev <- read.table(paste0(evdir,"/Afr.Amrcn.ex_related.EVs.txt"),header=T,stringsAsFactors = F)
hs1_ev <- read.table(paste0(evdir,"/Hispanic1.ex_related.EVs.txt"),header=T,stringsAsFactors = F)
hs2_ev <- read.table(paste0(evdir,"/Hispanic2.ex_related.EVs.txt"),header=T,stringsAsFactors = F)

# add race variable
idx_cau <- all_hemolysis$IID %in% cauc_ev$sample.id
idx_afr <- all_hemolysis$IID %in% afr_ev$sample.id
idx_easian <- all_hemolysis$IID %in% easian_ev$sample.id

all_hemolysis$Race<-rep(NA,dim(all_hemolysis)[1])
all_hemolysis$Race[idx_cau]="CAUCASIAN"
all_hemolysis$Race[idx_afr]="AFRAMRCN"
all_hemolysis$Race[idx_easian]="East ASIAN"
table(all_hemolysis$Race)

cauc_hemolysis <- all_hemolysis[idx_cau,]


exp <- read.table("GTEx_blood_tpm.txt",header=T,stringsAsFactors = F)
colnames(exp) <- gsub(".","-",fixed=T,colnames(exp))

pheno <- read.table("GTEx_v7_Annotations_SubjectPhenotypesDS.txt",header=T,sep="\t")

SEX <- unlist(lapply(colnames(exp),function(v) pheno$SEX[pheno$SUBJID %in% paste(unlist(strsplit(v,split="-",fixed=T))[1:2],collapse="-")]))
maleCol <- which(SEX==1)
femaleCol <- which(SEX==2)


library(ggplot2)
library(ggsignif)
library(plyr)
library(gtable)
library(grid)
library(gridExtra)

oriDat=cauc_hemolysis
box_plot_snp <- function(oriDat,snp,var,snpname,measurename,tag){
  #homo ref Vs. het Vs. homo alt
  oriDat_sub=oriDat[!is.na(oriDat[,snp]) & !is.na(oriDat$Gender) & (oriDat[,snp]==0 | oriDat[,snp]==1 | oriDat[,snp]==2),]
  dat=data.frame(
    measure=oriDat_sub[,var],
    gender=oriDat_sub$Gender,
    snp=factor(oriDat_sub[,snp])
  )
  gender_sum <- ddply(dat,.(gender),summarize,
                      N0=sum(snp==0 & !is.na(measure)),
                      N1=sum(snp==1 & !is.na(measure)),
                      N2=sum(snp==2 & !is.na(measure)))
  pvalues <- ddply(dat,.(gender),summarize,
                   sig=formatC(summary(aov(measure~snp))[[1]][["Pr(>F)"]][1],format="e",digits=2))
  p_annotations <- ifelse(as.numeric(pvalues[,2])<0.001,"***",
                          ifelse(as.numeric(pvalues[,2])<0.01,"**",
                                 ifelse(as.numeric(pvalues[,2])<0.05,"*","")))
  minorA=unlist(strsplit(snpname,'_',fixed=T))[2]
  majorA=ifelse(unlist(strsplit(unlist(strsplit(snpname,'_',fixed=T))[1],'.',fixed=T))[3]==minorA,
                unlist(strsplit(unlist(strsplit(snpname,'_',fixed=T))[1],'.',fixed=T))[4],
                unlist(strsplit(unlist(strsplit(snpname,'_',fixed=T))[1],'.',fixed=T))[3])
  rsID=unlist(strsplit(snpname,'.',fixed=T))[1]
  num_snp <- c(paste0("Females\n",
                      majorA,majorA,":",majorA,minorA,":",minorA,minorA,"\n",
                      gender_sum[1,2],":",gender_sum[1,3],":",gender_sum[1,4]),
               paste0("Males\n",
                      majorA,majorA,":",majorA,minorA,":",minorA,minorA,"\n",
                      gender_sum[2,2],":",gender_sum[2,3],":",gender_sum[2,4]))
  ylim1 = boxplot.stats(dat$measure)$stats[c(1, 5)]
  ystep=(ylim1[2]-ylim1[1])/10
  pt_size=16
  p=ggplot(dat,aes(x=gender,y=measure,fill=snp))+
    geom_boxplot(outlier.shape = NA,width=0.5)+
    scale_x_discrete(labels=num_snp)+
    theme(legend.position="none",
          text=element_text(size=pt_size,face="bold"),
          axis.text.x=element_text(size=pt_size,face="bold"),
          axis.text.y=element_text(size=pt_size,face="bold"),
          axis.title.y=element_text(size=pt_size,margin=margin(t=0,r=5,b=0,l=0),face="bold"),
          axis.title.x=element_text(size=pt_size,margin=margin(t=0,r=0,b=0,l=0)),
          plot.margin = unit(c(0.1,0.1,0,0.1),"cm")) +  #t,r,b,l
    #    ylim(ylim1[1],ylim1[2]+1.5*ystep) +
    labs(tag=tag)+
    ylab(measurename) + xlab("") +
    scale_y_continuous(limits=c(ylim1[1],ylim1[2]+1.5*ystep),
                       breaks=c(0,20,40,60),
                       labels=c("0","20%","40%","60%"))+
    geom_signif(xmin=c(0.7,1.7), xmax=c(1.3,2.3), y_position=ylim1[2]+0.8*ystep,
                annotation=p_annotations,tip_length = 0, size=0.6) 
    #    annotate("text",label=paste0("P=",pvalue),x=2,y=ylim1[2]+0.9*ystep,size=mm_size)+
    # ggtitle((paste0(rsID," in Caucasians")))
  # return(ggplotGrob(p))
  return(p)
}

plot_gene_tpm <- function(gene,exp,SEX,tag){
  dat <- data.frame(Sex=ifelse(SEX==1,"Males","Females"),
                    Expression=as.numeric(exp[exp$gene==GENE,1:537]))
  pvalue <- formatC(t.test(Expression~Sex,dat)$p.value,format="e",digits = 2)
  p_annotation <- ifelse(as.numeric(pvalue)<0.001,"***",
                         ifelse(as.numeric(pvalue)<0.01,"**",
                                ifelse(as.numeric(pvalue)<0.05,"*","")))
  ylim1 = boxplot.stats(dat$Expression)$stats[c(1, 5)]
  ystep=(ylim1[2]-ylim1[1])/10
  num_tpm <- c(paste0("Females\n(",
                      sum(dat$Sex=="Females",na.rm=T),")"),
               paste0("Males\n(",
                      sum(dat$Sex=="Males",na.rm=T),")"))
  pt_size=16
  p=ggplot(dat,aes(x=Sex,y=Expression,fill=Sex))+
    geom_boxplot(outlier.shape = NA,width=0.5)+
    scale_x_discrete(labels=num_tpm)+
    theme(legend.position="none",
          text=element_text(size=pt_size,face="bold"),
          title=element_text(size=pt_size,face="bold"),
          axis.text.x=element_text(size=17,face="bold",lineheight=1.2,margin=margin(t=10)),
          axis.text.y=element_text(size=pt_size),
          axis.title.y=element_text(size=pt_size,margin=margin(t=0,r=2,b=0,l=0)),
          plot.margin = unit(c(0.1,0.2,0,0),"cm")) +
    ylim(ylim1[1]-ystep,ylim1[2]+1.5*ystep)+ylab("TPM")+xlab("")+
    labs(tag=tag)+
    geom_signif(xmin=1, xmax=2, y_position=ylim1[2]+0.8*ystep,
                annotation=p_annotation,tip_length = 0, size=0.6) 
    # ggtitle((paste0(GENE," in GTEx")))
  return(p)
}

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

snp="rs13306780.42329004.A.C_A"
snpname=snp
GENE="SLC4A1"
var="SCRN.pink_pct_hemol.adj"
measurename="Osmotic hemolysis"

p_snp <- box_plot_snp(all_hemolysis,snp,"SCRN.pink_pct_hemol.adj",snpname,"Osmotic hemolysis","C")
p_gene <- plot_gene_tpm(GENE,exp,SEX,"D")

png(paste0("RBCOmics_osmotic_",GENE,"_",snpname,".png"),width = 8.5, height =3.5, units="in",res=100,pointsize=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 5)))
print(p_snp,  vp = vplayout(1, 1:3))
print(p_gene,  vp = vplayout(1, 4:5))
dev.off()
