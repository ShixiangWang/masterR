## Run TIMER codes in batch
## Statistical analysis of tumor-infiltrating immune cell inferences
## Wrap up for CancerImmunePipeline.R, iteratively analyze each tumor type and make summary tables and plots
## Created on Jan 31, 2015
library(ppcor)
library(reshape2)
#setwd('/Users/libo/Desktop/DFCI/microEnv/scripts/')
immune.cell.names=c('B cell','CD4 T cell','CD8 T cell','Neutrophil','Macrophage','Dendritic cell')
cancers=c('kich','blca','brca','cesc','gbm','hnsc','kirp','lgg','lihc','luad','lusc','prad','sarc','pcpg','paad','tgct','ucec','ov','skcm','dlbc','kirc','acc','meso','thca','uvm','ucs','thym','esca','stad','read','coad','chol')
cc.remove=c('sarc','pcpg','paad','tgct','meso','uvm','thym','esca','chol')
Checkpoint.genes=c('PDCD1','CD274','CTLA4') ## PD1, PDL1 and CTLA4
CT.root=c('MAGE','BAGE','GAGE','CTAG','SAGE','TAGE','CAGE','XAGE','SPAG') ## root for cancer/testis antigen genes
load('CancerTestis.genes.Rdata')
CancerTestis.genes0=CancerTestis.genes
MHC.genes=c('HLA-A','HLA-B','HLA-C','HLA-E','HLA-F','HLA-G',"HLA-DOA" ,"HLA-DMB" , "HLA-DRB1", "HLA-DRA",  "HLA-DQA1", "HLA-DOB","HLA-DMA" , "HLA-DPB1" ,"HLA-DQB1", "HLA-DPA1", "HLA-DRB6")

ImmuneGeneEnrichment=rep(NA,length(cancers))    ## odds ratio for immune gene enrichment in tumor microenvironment
NeutrophilNecrosis.cor=rep(NA,length(cancers))  ## correlation between neutrophils and necrosis
filteredImmuneGenes=rep(NA,length(cancers)) ## number of selected immune genes for analysis
SurvivalSingle=matrix(NA,ncol=length(cancers),nrow=6)   ## single variable survival analysis
SurvivalFull=matrix(NA,ncol=length(cancers),nrow=6) ## full survival analysis model
SurvivalSingle.p=matrix(NA,ncol=length(cancers),nrow=6)   ## single variable survival analysis, purity corrected
SurvivalFull.p=matrix(NA,ncol=length(cancers),nrow=6) ## full survival analysis model, purity corrected
Age.cor=matrix(NA,ncol=length(cancers),nrow=6)  ## correlation with age
Stage.cor=matrix(NA,ncol=length(cancers),nrow=6)    ## correlation with clinical stage
Grade.cor=matrix(NA,ncol=length(cancers),nrow=6)
TumorStatus.cor=matrix(NA,ncol=length(cancers),nrow=6)  ## fold change for tumor status
Race.cor=matrix(NA,ncol=length(cancers),nrow=6) ## fold change for race
AdjacentTissue=matrix(NA,ncol=length(cancers),nrow=6)   ## fold change for tumor vs adjacent tissue
SurvivalSingleWithout=matrix(NA,ncol=length(cancers),nrow=6)
Gender.cor=Race.cor
Top.genes=c()
CYT.list=c()
CYT.HR=c()
sim.COR.mat=c()
CheckPoint.Blockade.PurityMat=matrix(NA,ncol=length(cancers),nrow=length(Checkpoint.genes))
MHC.PurityMat=matrix(NA,ncol=length(cancers),nrow=length(MHC.genes))
rownames(CheckPoint.Blockade.PurityMat)=Checkpoint.genes
rownames(MHC.PurityMat)=MHC.genes
CheckPoint.list=c()
CancerTestis=c()

getID <- function(sID,num.res=3){
    mm=c()
    for(id in sID){
        tmp=unlist(strsplit(id,'-'))
        if(length(tmp)==1){
            tmp=unlist(strsplit(id,'\\.'))
        }
        ll='TCGA'
        for(j in 2:num.res){
            ll=paste(ll,tmp[j],sep='-')
        }
        mm=c(mm,ll)
    }
    return(mm)
}


names(ImmuneGeneEnrichment)=names(NeutrophilNecrosis.cor)=names(filteredImmuneGenes)=colnames(SurvivalSingle)=colnames(SurvivalFull)=colnames(SurvivalSingle.p)=colnames(SurvivalFull.p)=colnames(Age.cor)=colnames(Stage.cor)=colnames(TumorStatus.cor)=colnames(Race.cor)=colnames(AdjacentTissue)=colnames(SurvivalSingleWithout)=colnames(CheckPoint.Blockade.PurityMat)=colnames(MHC.PurityMat)=colnames(Grade.cor)=colnames(Gender.cor)=cancers
rownames(SurvivalSingle)=rownames(SurvivalFull)=rownames(SurvivalSingle.p)=rownames(SurvivalFull.p)=rownames(Age.cor)=rownames(Stage.cor)=rownames(TumorStatus.cor)=rownames(Race.cor)=rownames(AdjacentTissue)=rownames(SurvivalSingleWithout)=rownames(Grade.cor)=rownames(Gender.cor)=c('B_cell','T_cell.CD4','T_cell.CD8','Neutrophil','Macrophage','DC')

FMAT.all=c()
AdjacentTissue1=AdjacentTissue
#MarkPvalue <- function(x,p){
#    if(p<=0.05)x=paste(x,'*')
#    if(p>0.05&p<=0.1)x=paste(x,'.')
#    return(x)
#}

MarkPvalue <- function(x,p)return(paste(x,p,sep=' '))

SurvivalDataList=c()

for(cc in cancers){
    source('CancerImmunePipeline.R')
    
    sim.COR.mat=c(sim.COR.mat,list(sim.COR))
    Top.genes=c(Top.genes,list(top.gg))
    FMAT.all=c(FMAT.all,list(Fmat))
    SurvivalDataList=c(SurvivalDataList,list(dat.surv))
    
    ## Enrichment analysis
    n1=length(vv)   ## immune and correlated
    n2=length(vv.t)-n1 ## non-immune and correlated
    tmp.ss0=intersect(rownames(curated.ref.genes.agg.br),rownames(dd.br))
    tmp.immune=intersect(marker.list.genes,tmp.ss0)
    n3=length(tmp.immune)-n1 # immune and non-correlated
    n4=length(tmp.ss0)-n1-n2-n3 # non-immune and non-correlated
    tmp=fisher.test(matrix(c(n1,n2,n3,n4),2,2))
    xx=tmp$estimate
    tmp.p=tmp$p.value
    xx=MarkPvalue(xx,tmp.p)
    ImmuneGeneEnrichment[cc]=xx
    filteredImmuneGenes[cc]=n.immune
    
    ## Necrosis
    tmp.ss=intersect(rownames(necrosisData),rownames(Fmat0.p))
    if('Neutrophil' %in% colnames(Fmat0.p))tmp=cor.test(as.numeric(necrosisData[tmp.ss,8]),Fmat0.p[tmp.ss,'Neutrophil'],use='complete.obs') else tmp=list(estimate=NA,p.value=1)
    xx=tmp$estimate
    tmp.p=tmp$p.value
    xx=MarkPvalue(xx,tmp.p)
    NeutrophilNecrosis.cor[cc]=xx
    
    ## Cytolytic Activity
    CYT.list=c(CYT.list,CYT)
    tmp=CYT[[1]]
    tmp[which(!is.finite(tmp))]=NA
    tmp=tmp[grep(cc.type,names(tmp))]
    names(tmp)=getID(names(tmp))
    tmp.ss=intersect(rownames(dd.surv),names(tmp))
    tmp.Surv=Surv(dat.surv[tmp.ss,2],dat.surv[tmp.ss,1])
    if(length(cc.stages)>0){
        if(length(cc.grades)>0){
            tmp.HR=summary(coxph(tmp.Surv~tmp[tmp.ss]+cc.ages[tmp.ss]+cc.stages[tmp.ss]+cc.grades[tmp.ss]))$coefficients
        }else{
            tmp.HR=summary(coxph(tmp.Surv~tmp[tmp.ss]+cc.ages[tmp.ss]+cc.stages[tmp.ss]))$coefficients
        }
    }else{
        if(length(cc.grades)>0){
            tmp.HR=summary(coxph(tmp.Surv~tmp[tmp.ss]+cc.ages[tmp.ss]+cc.grades[tmp.ss]))$coefficients
        }else tmp.HR=summary(coxph(tmp.Surv~tmp[tmp.ss]+cc.ages[tmp.ss]))$coefficients
    }
    CYT.HR=c(CYT.HR,MarkPvalue(tmp.HR[1,2],tmp.HR[1,5]))
    
    ## Checkpoint blockade genes
    tmp.ss0=intersect(rownames(Fmat),rownames(AGP))
    tmp.CB=c()
    for(i in Checkpoint.genes){
        if(!i %in% rownames(dd.br)){
            tmp.CB=c(tmp.CB,NA)
            next
        }
        tmp=cor.test(dd.br[i,tmp.ss0],AGP[tmp.ss0,2])
        tmp.p=tmp$p.value
        tmp.c=tmp$estimate
        tmp.c=MarkPvalue(tmp.c,tmp.p)
        tmp.CB=c(tmp.CB,tmp.c)
    }
    CheckPoint.Blockade.PurityMat[,cc]=tmp.CB
    
    tmp.mat=matrix(NA,nrow=4,ncol=length(Checkpoint.genes))
    rownames(tmp.mat)=c('T_cell.CD4','T_cell.CD8','Macrophage','DC')
    colnames(tmp.mat)=Checkpoint.genes
    for(i in Checkpoint.genes){
        if(!i %in% rownames(dd.br))next
        for(j in c('T_cell.CD4','T_cell.CD8','Macrophage','DC')){
            tmp=pcor.test(dd.br[i,tmp.ss0],Fmat[tmp.ss0,j],AGP[tmp.ss0,2])
            tmp.mat[j,i]=MarkPvalue(tmp$estimate,tmp$p.value)
        }
    }
    CheckPoint.list=c(CheckPoint.list,list(tmp.mat))
    
    ## Cancer Testis Antigen
    #    CancerTestis.genes1=c()
    #for(rr in CT.root){
    #   CancerTestis.genes1=c(CancerTestis.genes1,rownames(dd)[grep(rr,rownames(dd))])
    #}
    #CancerTestis.genes1=unique(c(CancerTestis.genes1,CancerTestis.genes0))
    #CancerTestis.genes=intersect(CancerTestis.genes1,rownames(dd))
    #tmp.ss0=intersect(rownames(AGP),colnames(dd))
    #tmp.corE=c()
    #tmp.corP=c()
    #for(tmp.gg in CancerTestis.genes){
    #    tmp=cor.test(dd[tmp.gg,tmp.ss0],AGP[tmp.ss0,2])
    #    tmp.corE=c(tmp.corE,tmp$estimate)
    #    tmp.corP=c(tmp.corP,tmp$p.value)
    #}
    #tmp.vv=which(tmp.corE>= 0.05 & tmp.corP <=0.05)
#tmp.vv=which(cor(t(dd[CancerTestis.genes,tmp.ss0]),AGP[tmp.ss0,2])>= 0)
#CT.genes=CancerTestis.genes[tmp.vv]

#   tmp.mat=matrix(NA,nrow=6,ncol=length(CT.genes))
#   rownames(tmp.mat)=colnames(Fmat)
#   colnames(tmp.mat)=CT.genes
#   for(i in CT.genes){
#       for(j in colnames(Fmat)){
#           tmp=pcor.test(dd[i,tmp.ss0],Fmat[tmp.ss0,j],AGP[tmp.ss0,2],method='s')
#           tmp.mat[j,i]=MarkPvalue(tmp$estimate,tmp$p.value)
#       }
#   }
#   CancerTestis=c(CancerTestis,list(tmp.mat))
    
    ## MHC genes
    tmp.ss0=intersect(rownames(Fmat),rownames(AGP))
    tmp.CB=c()
    for(i in MHC.genes){
        if(!i %in% rownames(dd.br)){
            tmp.CB=c(tmp.CB,NA)
            next
        }
        tmp=cor.test(dd.br[i,tmp.ss0],AGP[tmp.ss0,2])
        tmp.p=tmp$p.value
        tmp.c=tmp$estimate
        tmp.c=MarkPvalue(tmp.c,tmp.p)
        tmp.CB=c(tmp.CB,tmp.c)
    }
    MHC.PurityMat[,cc]=tmp.CB

    ## Survival Analysis without purity correction
    n.death=length(which(dat.surv[,1]==1))
    if(n.death>= 1){
    tmp.ss=intersect(rownames(dat.surv),rownames(Fmat0.p))
    tmp.Surv=Surv(dat.surv[tmp.ss,2],dat.surv[tmp.ss,1])
    for(i in 1:ncol(Fmat0.p)){
        if(length(cc.stages)>0)tmp=summary(coxph(tmp.Surv~Fmat0.p[tmp.ss,i]+cc.ages[tmp.ss]+cc.stages[tmp.ss]))$coefficients else tmp=summary(coxph(tmp.Surv~Fmat0.p[tmp.ss,i]+cc.ages[tmp.ss]))$coefficients
        xx=tmp[1,2]
        tmp.p=tmp[1,5]
        xx=MarkPvalue(xx,tmp.p)
        SurvivalSingle[colnames(Fmat0.p)[i],cc]=xx
    }
    if(length(cc.stages)>0){
        if(length(cc.grades)>0){
            tmp=summary(coxph(tmp.Surv~Fmat0.p[tmp.ss,]+cc.ages[tmp.ss]+cc.stages[tmp.ss]+cc.grades[tmp.ss]))$coefficients
        }else{
            tmp=summary(coxph(tmp.Surv~Fmat0.p[tmp.ss,]+cc.ages[tmp.ss]+cc.stages[tmp.ss]))$coefficients
        }
    }else{
        if(length(cc.grades)>0){
            tmp=summary(coxph(tmp.Surv~Fmat0.p[tmp.ss,]+cc.ages[tmp.ss]+cc.grades[tmp.ss]))$coefficients
        }else tmp=summary(coxph(tmp.Surv~Fmat0.p[tmp.ss,]+cc.ages[tmp.ss]))$coefficients
    }
    for(i in 1:ncol(Fmat0.p)){
        xx=tmp[i,2]
        tmp.p=tmp[i,5]
        xx=MarkPvalue(xx,tmp.p)
        SurvivalFull[colnames(Fmat0.p)[i],cc]=xx
    }
    
    # Survival Analysis single without any adjustment
    for(i in 1:ncol(Fmat0.p)){
        tmp=summary(coxph(tmp.Surv~Fmat0.p[tmp.ss,i]))$coefficients
        xx=tmp[1,2]
        tmp.p=tmp[1,5]
        xx=MarkPvalue(xx,tmp.p)
        SurvivalSingleWithout[colnames(Fmat0.p)[i],cc]=xx
    }
    }
    
    ## Age correlation
    tmp.ss=intersect(rownames(dat.surv),rownames(Fmat0.p))
    for(i in 1:ncol(Fmat0.p)){
        tmp=cor.test(Fmat0.p[tmp.ss,i],cc.ages[tmp.ss],method='s',use='complete.obs')
        xx=tmp$estimate
        tmp.p=tmp$p.value
        xx=MarkPvalue(xx,tmp.p)
        Age.cor[colnames(Fmat0.p)[i],cc]=xx
    }
    
    ## Stage correlation
    if(length(cc.stages)>0){
        for(i in 1:ncol(Fmat0.p)){
            tmp=cor.test(Fmat0.p[tmp.ss,i],as.numeric(as.factor(cc.stages[tmp.ss])),method='s',use='complete.obs')
            xx=tmp$estimate
            tmp.p=tmp$p.value
            xx=MarkPvalue(xx,tmp.p)
            Stage.cor[colnames(Fmat0.p)[i],cc]=xx
        }
    }
    
    ## Grade correlation
    if(length(cc.grades)>0){
        for(i in 1:ncol(Fmat0.p)){
            tmp=cor.test(Fmat0.p[tmp.ss,i],as.numeric(as.factor(cc.grades[tmp.ss])),method='s',use='complete.obs')
            xx=tmp$estimate
            tmp.p=tmp$p.value
            xx=MarkPvalue(xx,tmp.p)
            Grade.cor[colnames(Fmat0.p)[i],cc]=xx
        }
    }

    ## Tumor Status
    for(i in 1:ncol(Fmat0.p)){
        if(length(cc.stages)>0)tmp=summary(glm(tumor_status.c[tmp.ss]~Fmat0.p[tmp.ss,i]+cc.ages[tmp.ss]+cc.stages[tmp.ss]))$coefficients else tmp=summary(glm(tumor_status.c[tmp.ss]~Fmat0.p[tmp.ss,i]+cc.ages[tmp.ss]))$coefficients
        xx=tmp[2,1]
        tmp.p=tmp[2,4]
        xx=MarkPvalue(xx,tmp.p)
        TumorStatus.cor[colnames(Fmat0.p)[i],cc]=xx
    }
    
    ## Gender
    gender=c()
    tmp.vv=grep('gender',colnames(dd.surv))
    if(length(tmp.vv)>1)tmp.vv=tmp.vv[1]
    gender=dd.surv[,tmp.vv]
    gender[grep('\\[',gender)]=NA
    names(gender)=rownames(dd.surv)
    if(length(unique(gender))>1){
        for(i in 1:ncol(Fmat0.p)){
            tmp=wilcox.test(Fmat0.p[tmp.ss,i]~gender[tmp.ss])
            Gender.cor[i,cc]=MarkPvalue(tmp$estimate,tmp$p.value)
        }
    }
    
    ## Race
    race=c()
    tmp.vv=grep('race',colnames(dd.surv))
    if(length(tmp.vv)>1)tmp.vv=tmp.vv[2]
    race=dd.surv[,tmp.vv]
    race[grep('\\[',race)]=NA
    if(length(race)>0 & length(unique(na.omit(race)))>1){
        for(i in 1:ncol(Fmat0.p)){
            tmp=summary(aov(lm(Fmat0.p[tmp.ss,i]~as.factor(race[tmp.ss]))))
            tmp.p= tmp[[1]][5][1,1]
            Race.cor[colnames(Fmat0.p)[i],cc]=tmp.p
        }
    }
    
    ## Adjacent Tissue
    if(cc!='skcm')tmp.11A=grep('11A',rownames(Fmat)) else tmp.11A=grep('01A',rownames(Fmat))
    tmp.01A=grep(cc.type,rownames(Fmat))
    tmp.tissue=rep(0,length(c(tmp.01A,tmp.11A)))
    names(tmp.tissue)=colnames(dd.br)[c(tmp.01A,tmp.11A)]
    tmp.tissue[colnames(dd.br)[tmp.11A]]=1
    if(length(tmp.11A)>0){
        gg='PECAM1'
    for(i in 1:ncol(Fmat)){
        tmp=wilcox.test(Fmat[tmp.01A,i],Fmat[tmp.11A,i])
        xx=tmp$estimate[1]/tmp$estimate[2]
        tmp.p=tmp$p.value
        #if(! gg %in% rownames(dd.br)){
        #   AdjacentTissue[colnames(Fmat0.p)[i],cc]=NA
        #   tmp=wilcox.test(Fmat[tmp.01A,i],Fmat[tmp.11A,i])
        #   xx=MarkPvalue(median(Fmat[tmp.01A,i],na.rm=T)/median(Fmat[tmp.11A,i],na.rm=T),tmp$p.value)
        #   AdjacentTissue1[colnames(Fmat0.p)[i],cc]=xx
        #   next
        #}
        #tmp=summary(glm(tmp.tissue~Fmat[names(tmp.tissue),i]+dd.br[gg,names(tmp.tissue)],family=binomial(link=logit)))$coefficients
        #xx=exp(-tmp[2,1])
        #tmp.p=tmp[2,4]
        xx=MarkPvalue(xx,tmp.p)
        AdjacentTissue[colnames(Fmat0.p)[i],cc]=xx
        
        tmp=wilcox.test(Fmat[tmp.01A,i],Fmat[tmp.11A,i])
        xx=MarkPvalue(median(Fmat[tmp.01A,i],na.rm=T)/median(Fmat[tmp.11A,i],na.rm=T),tmp$p.value)
        AdjacentTissue1[colnames(Fmat0.p)[i],cc]=xx
    }
    }
    
    ## Survival Analysis with purity correction
    tmp.ss=intersect(rownames(dat.surv),rownames(Fmat.res))
    n.death=length(which(dat.surv[tmp.ss,1]==1))
    if(n.death>=1){
    tmp.Surv=Surv(dat.surv[tmp.ss,2],dat.surv[tmp.ss,1])
    for(i in 1:ncol(Fmat.res)){
        if(length(cc.stages)>0)tmp=summary(coxph(tmp.Surv~Fmat.res[tmp.ss,i]+cc.ages[tmp.ss]+cc.stages[tmp.ss]))$coefficients else tmp=summary(coxph(tmp.Surv~Fmat.res[tmp.ss,i]+cc.ages[tmp.ss]))$coefficients
        xx=tmp[1,2]
        tmp.p=tmp[1,5]
        xx=MarkPvalue(xx,tmp.p)
        SurvivalSingle.p[colnames(Fmat.res)[i],cc]=xx
    }
    if(length(cc.stages)>0)tmp=summary(coxph(tmp.Surv~Fmat.res[tmp.ss,]+cc.ages[tmp.ss]+cc.stages[tmp.ss]))$coefficients else tmp=summary(coxph(tmp.Surv~Fmat.res[tmp.ss,]+cc.ages[tmp.ss]))$coefficients
    for(i in 1:ncol(Fmat.res)){
        xx=tmp[i,2]
        tmp.p=tmp[i,5]
        xx=MarkPvalue(xx,tmp.p)
        SurvivalFull.p[colnames(Fmat.res)[i],cc]=xx
    }
    }
}

names(SurvivalDataList)=names(FMAT.all)=names(Top.genes)=cancers
###############------------------This section is for generating plots using above results----------------------##########

## heatmap of Hazard ratios, 10 color levels, <0.1 or >10 will be aggregated to one color, add HR values to the cell

HeatMap.0 <- function(x,scale=c('HR','Cor','OR','FC'),title.lab='',sort=T,FDR=0.1){

if(length(title.lab)==0)title.lab=deparse(substitute(x))

redblue.dat=read.table('/Users/libo/Desktop/bo/redblue100.txt',header=T)
redblue6=redblue.dat[seq(1,101,length.out=6),]
redblue.colors=rgb(redblue6[,1],redblue6[,2],redblue6[,3])

if(scale=='Cor10'){
    redblue10=redblue.dat[seq(1,101,length.out=10),]
    redblue.colors=rgb(redblue10[,1],redblue10[,2],redblue10[,3])
}
require(qvalue)
vv.NA=which(is.na(x))
x[vv.NA]=0
colnames(x)=toupper(colnames(x))
tmp.x=x[,order(colnames(x))]
x=x[,order(colnames(x))]
if(typeof(x)!='double'){
    tmp.x=gsub(' .+','',tmp.x)
    mode(tmp.x)='numeric'
}
tmp.x[which(!is.finite(tmp.x))]=9999
tmp.x[which(tmp.x>9999)]=9999
if(sort){
    if(scale=='HR')tmp.hh=heatmap(log(tmp.x+0.001),keep.dendro=T) else tmp.hh=heatmap(tmp.x,keep.dendro=T)
    tmp.x=tmp.x[,tmp.hh$colInd]
    x=x[,tmp.hh$colInd]
}
x[which(x==0)]=NA
tmp.x[which(tmp.x==0)]=NA

x[which(is.na(x))]='NA NA'
tmp=strsplit(as.character(melt(x[,])[,3]),' ')
sig.levels=c()
#for(i in tmp){
#    if(length(i)==1)sig.levels=c(sig.levels,'') else sig.levels=c(sig.levels,i[2])
#}
pp=sapply(tmp,function(x)x[[2]])
pp=as.numeric(pp)
vv.pp.Nna=which(!is.na(pp))
sig.levels=rep('',length(pp))
if(FDR){
    #pp=p.adjust(pp,method='BH')
    pp[vv.pp.Nna]=qvalue(pp[vv.pp.Nna])$qvalue
    sig.levels[vv.pp.Nna][which(pp[vv.pp.Nna]<=FDR)]='*'
}else{
    sig.levels[vv.pp.Nna][which(pp[vv.pp.Nna]<=0.05)]='*'
}

if(!is.null(nrow(tmp.x)))df.merge=cbind(melt(tmp.x),sig.levels) else df.merge=cbind(melt(x),sig.levels)
colnames(df.merge)=c('Cell_Type','Cancer_Type','Value','Signif')

if(scale=='HR')Value=cut(log10(as.numeric(as.character(df.merge$Value))+0.0000001),breaks=c(-Inf,seq(-1,1,length.out=5),Inf),label=c('<0.1','0.1-0.32','0.32-1.0','1.0-3.2','3.2-10','>10'))
if(scale=='Cor')Value=cut(as.numeric(as.character(df.merge$Value))+0.0000001,breaks=c(-1,seq(-0.2,0.2,length.out=5),1),label=c('-1 - -0.2','-0.2 - -0.1','-0.1 - 0','0 - 0.1','0.1 - 0.2','0.2 - 1'))
if(scale=='Cor10')Value=cut(as.numeric(as.character(df.merge$Value))+0.0000001,breaks=c(-1,seq(-0.8,0.8,length.out=9),1),label=c('-1 - -0.8','-0.8 - -0.6','-0.6 - -0.4','-0.4 - -0.2','-0.2 - 0','0 - 0.2','0.2 - 0.4','0.4 - 0.6','0.6 - 0.8','0.8 - 1'))
if(scale=='OR')Value=cut(as.numeric(as.character(df.merge$Value))+0.0000001,breaks=c(-Inf,seq(-0.8,0.8,length.out=5),Inf),label=c('<0.45','0.45-0.67','0.67-1.0','1.0-1.5','1.5-2.2','>2.2'))
if(scale=='FC')Value=cut(as.numeric(as.character(df.merge$Value))+0.0000001,breaks=c(0,seq(0.2,1.8,length.out=5),Inf),label=c('0-0.2','0.2-0.6','0.6-1.0','1.0-1.4','1.4-1.8','>1.8'))

df.merge$Value=Value
names(redblue.colors)=levels(Value)
colScale <- scale_fill_manual(name='Value',values=redblue.colors)
if(length(which(is.na(x)))>0)ggplot(df.merge,aes(x=Cancer_Type,y=Cell_Type,label=Signif))+geom_tile(data=subset(df.merge, !is.na(Value)),aes(fill=Value))+geom_tile(data=subset(df.merge,is.na(Value)),aes(colour= 'NA' ),fill='gray',linetype=0,alpha=0.8)+geom_text()+colScale+theme(axis.text.x=element_text(size=12,angle=90,vjust=0.7),axis.text.y=element_text(size=12),axis.title=element_text(size=0))+labs(title=title.lab) else ggplot(df.merge,aes(x=Cancer_Type,y=Cell_Type,label=Signif))+geom_tile(data=subset(df.merge, !is.na(Value)),aes(fill=Value))+geom_text()+colScale+theme(axis.text.x=element_text(size=12,angle=90,vjust=0.7),axis.text.y=element_text(size=12),axis.title=element_text(size=0))+labs(title=title.lab)
}

HeatMap.p <- function(x,scale=c('Expr','Cor'),title.lab='',FDR=0.1){
    require(reshape2)
    require(RColorBrewer)
    require(qvalue)
    tmp.col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10)
    
    mx=melt(x)
    tmp=strsplit(as.character(mx[,3]),' ')
    tmp.v=as.numeric(sapply(tmp,function(x)x[[1]]))
    if(scale=='Expr'){
        tmp.v=-tmp.v
        qqs=quantile(tmp.v,seq(0,1,length.out=11))
        print(qqs)
        Value=cut(tmp.v,breaks=qqs,label=1:10)
    }else{
        Value=cut(tmp.v,breaks=seq(-1,1,length.out=11),label=1:10)
    }
    #tmp.p=p.adjust(as.numeric(sapply(tmp,function(x)x[[2]])),method='BH')
    tmp.p=qvalue(as.numeric(sapply(tmp,function(x)x[[2]])))$qvalue
    mx=cbind(mx[,1:2],Value,tmp.p)
    
    Nx=ncol(x)
    Ny=nrow(x)
    par(mar=c(2,3,2,0.2),mgp=c(3,0.1,0),las=2)
    plot(1,1,cex=0,xlim=c(0.5,Nx+2),ylim=c(0.5,Ny+0.5),axes=F,bty='n',xlab='',ylab='',xaxs="i", yaxs="i",main=title.lab,cex.main=0.8)
    for(i in 1:Nx){
        for(j in 1:Ny){
            tmp.vv=which(mx[,1]==rownames(x)[j] & mx[,2]==colnames(x)[i])
            ccll=tmp.col[as.numeric(mx[tmp.vv,3])]
            if(mx[tmp.vv,4]<=FDR){
                rect(i-0.5,j-0.5,i+0.5,j+0.5,col=ccll)
            }else{
                rect(i-0.5,j-0.5,i+0.5,j+0.5,border=ccll,lwd=1)
                segments(i-0.5,j-0.5,i+0.5,j+0.5,col=ccll,lwd=1)
            }
        }
    }
    for(i in 1:10){
        rect(Nx+1,Ny-i+0.5,Nx+2,Ny-i+1.5,col=tmp.col[i],border=NA)
    }
    rect(Nx+1,Ny-12+0.5,Nx+2,Ny-12+1.5,col='gray10')
    rect(Nx+1,Ny-13+0.5,Nx+2,Ny-13+1.5,border='gray10')
    segments(Nx+1,Ny-13+0.5,Nx+2,Ny-13+1.5,col='gray10')
    axis(2,at=1:Ny,rownames(x),cex.axis=0.8,tick=F)
    axis(1,at=1:Nx,toupper(colnames(x)),cex.axis=0.8,tick=F)
}

HeatMap.COR <- function(x,title.lab='',FDR=0.2){

    if(length(title.lab)==0)title.lab=deparse(substitute(x))
    redblue.dat=read.table('/Users/libo/Desktop/bo/redblue100.txt',header=T)
    redblue2=c('dodgerblue4','firebrick3')
    
    vv.NA=which(is.na(x))
    x[vv.NA]='0 1'
    Value.Mat=matrix(NA,ncol=ncol(x),nrow=nrow(x))
    P.Mat=Value.Mat
    rownames(P.Mat)=rownames(Value.Mat)=rownames(x)
    colnames(P.Mat)=colnames(Value.Mat)=colnames(x)
    for(i in 1:nrow(x)){
        tmp=strsplit(x[i,],' ')
        vv=sapply(tmp,function(x)x[[1]])
        pp=sapply(tmp,function(x)x[[2]])
        Value.Mat[i,]=vv
        P.Mat[i,]=pp
    }
    P.Mat.Adjusted=matrix(p.adjust(P.Mat,method='BH'),ncol=ncol(x),nrow=nrow(x),byrow=F)
    rownames(P.Mat.Adjusted)=rownames(P.Mat)
    colnames(P.Mat.Adjusted)=colnames(P.Mat)
    mode(Value.Mat)='numeric'
    Symbol.Mat=matrix('',ncol=ncol(x),nrow=nrow(x))
    Symbol.Mat[which(P.Mat.Adjusted<=FDR)]='*'
    
    tmp.hist=heatmap(Value.Mat)
    Value.Mat=Value.Mat[,tmp.hist$colInd]
    Symbol.Mat=Symbol.Mat[,tmp.hist$colInd]
    
    layout(matrix(c(rep(1,8),2),nrow=1))
    par(mar=c(0,6,6.5,0),las=2,cex=0.8)
    plot(0,0,xlim=c(0.8,ncol(x)+0.2),ylim=c(0.8,6.2),cex=0,xlab='',ylab='',axes=F,main=title.lab)
    rect(0.7,0.7,ncol(x)+0.3,6.3,col=rgb(30,30,30,alpha=30,max=255),border=NA)
    axis(3,at=1:ncol(x),toupper(colnames(x)),font=2)
    axis(2,at=1:nrow(x),rownames(x),font=2)
    for(i in 1:ncol(x))points(rep(i,6),1:6,cex=1.5,pch=15,col=redblue2[(sign(Value.Mat[,i])+3)/2])
    for(i in 1:ncol(x))text(rep(i,6),1:6,Symbol.Mat[,i],cex=1,col='yellow')
    #par(mar=c(0,0,0,0),las=1)
    #plot(0,0,cex=0,axes=0,xlab='',ylab='',xlim=c(-1,1),ylim=c(-1,1))
    #legend(-1,0.8,legend=c('P>=0.1','0.05<P<=0.1','P<=0.05'),pch=15,pt.cex=c(0.3,0.6,0.9)*1.5,cex=0.7)
    #legend(-1,0.2,legend=levels(Value),pch=15,col=redblue.colors,cex=0.7)
}

HeatMap.1 <- function(x,scale=c('HR'),title.lab='',display.dendro=T, FDR=0.15,exclude=0){
    quartz(width=10,height=4.3)
    if(length(title.lab)==0)title.lab=deparse(substitute(x))
    redblue.dat=read.table('/Users/libo/Desktop/bo/redblue100.txt',header=T)
    redblue6=redblue.dat[seq(1,101,length.out=6),]
    redblue.colors=rgb(redblue6[,1],redblue6[,2],redblue6[,3])
    require(qvalue)
    
    vv.NA=which(is.na(x))
    x[vv.NA]=0
    tmp.x=x[,order(colnames(x))]
    x=x[,order(colnames(x))]
    if(typeof(x)!='double'){
    tmp.x=gsub(' .+','',tmp.x)
    mode(tmp.x)='numeric'
    }
    tmp.x[which(!is.finite(tmp.x))]=9999
    tmp.x[which(tmp.x>9999)]=9999
    if(scale=='HR')tmp.hh=heatmap(log(tmp.x+0.001),keep.dendro=T) else tmp.hh=heatmap(tmp.x,keep.dendro=T)
    tmp.x=tmp.x[,tmp.hh$colInd]
    x=x[,tmp.hh$colInd]
    x[which(x==0)]=NA
    tmp.x[which(tmp.x==0)]=NA
    
    x[which(is.na(x))]='NA NA'
    tmp=strsplit(as.character(melt(x[,])[,3]),' ')
    sig.levels=c()
    #for(i in tmp){
    #    if(length(i)==1)sig.levels=c(sig.levels,'') else sig.levels=c(sig.levels,i[2])
    #}
    pp=sapply(tmp,function(x)x[[2]])
    pp=as.numeric(pp)
    vv.pp.Nna=which(!is.na(pp))
    sig.levels=rep('',length(pp))
    if(FDR){
        #pp=p.adjust(pp,method='BH')
        pp[vv.pp.Nna]=qvalue(pp[vv.pp.Nna])$qvalue
        sig.levels[vv.pp.Nna][which(pp[vv.pp.Nna]<=FDR)]='*'
    }else{
        sig.levels[vv.pp.Nna][which(pp[vv.pp.Nna]<=0.05)]='*'
    }
    sig.cex=c()
    for(i in sig.levels){
        if(i=='')sig.cex=c(sig.cex,0.3)
        if(i=='*')sig.cex=c(sig.cex,0.9)
    }
    sig.cex=matrix(sig.cex,byrow=F,ncol=ncol(x))
    
    if(!is.null(nrow(tmp.x)))df.merge=cbind(melt(tmp.x),sig.levels) else df.merge=cbind(melt(x),sig.levels)
    colnames(df.merge)=c('Cell_Type','Cancer_Type','Value','Signif')
    
    if(scale=='HR')Value=cut(log10(as.numeric(as.character(df.merge$Value))+0.0000001),breaks=c(-Inf,seq(-1,1,length.out=5),Inf),label=c('<0.1','0.1-0.32','0.32-1.0','1.0-3.2','3.2-10','>10'))
    if(scale=='Cor')Value=cut(as.numeric(as.character(df.merge$Value))+0.0000001,breaks=c(-1,seq(-0.2,0.2,length.out=5),1),label=c('-1 - -0.2','-0.2 - -0.1','-0.1 - 0','0 - 0.1','0.1 - 0.2','0.2 - 1'))
    if(scale=='OR')Value=cut(as.numeric(as.character(df.merge$Value))+0.0000001,breaks=c(-Inf,seq(-0.8,0.8,length.out=5),Inf),label=c('<0.45','0.45-0.67','0.67-1.0','1.0-1.5','1.5-2.2','>2.2'))
    if(scale=='FC')Value=cut(as.numeric(as.character(df.merge$Value))+0.0000001,breaks=c(0,seq(0.2,1.8,length.out=5),Inf),label=c('0-0.2','0.2-0.6','0.6-1.0','1.0-1.4','1.4-1.8','>1.8'))

    Value.mat=matrix(as.numeric(as.factor(Value)),byrow=F,ncol=ncol(x))
    if(!display.dendro){
        layout(matrix(c(rep(1,8),2),nrow=1))
        par(mar=c(4,6,2,0),las=2,cex=1.1)
        plot(0,0,xlim=c(0.8,ncol(x)+0.2),ylim=c(0.8,nrow(x)+.2),cex=0,xlab='',ylab='',axes=F,main=title.lab)
        rect(0.7,0.7,ncol(x)+0.3,nrow(x)+.3,col=rgb(30,30,30,alpha=30,max=255),border=NA)
        axis(1,at=1:ncol(x),toupper(colnames(x)),font=2)
        axis(2,at=1:nrow(x),rownames(x),font=2)
        for(i in 1:ncol(x))points(rep(i,nrow(x)),1:nrow(x),cex=sig.cex[,i]*1.5,pch=15,col=redblue.colors[Value.mat[,i]])
        par(mar=c(0,0,0,0),las=1)
        plot(0,0,cex=0,axes=0,xlab='',ylab='',xlim=c(-1,1),ylim=c(-1,1))
        if(FDR)legend(-1,0.8,legend=paste(c('q>','q<='),c(FDR,FDR),sep=''),pch=15,pt.cex=c(0.3,0.9)*1.5,cex=0.7)
        
        if(!FDR)legend(-1,0.8,legend=c('P>0.05','P<=0.05'),pch=15,pt.cex=c(0.3,0.9)*1.5,cex=0.7)
        legend(-1,0.2,legend=levels(Value),pch=15,col=redblue.colors,cex=0.7)
        return()
    }
    layout(matrix(c(1,0,2,3),nrow=2,byrow=T),width=c(8,1),height=c(1,2))
    par(mai=c(0,1.4,0.4,0))
    plot(tmp.hh$Colv,leaflab='none',yaxt='n',main=title.lab,xlim=c(0.8,ncol(x)+0.2))
    par(mai=c(0.9,1.4,0,0),las=2,cex=1.1)
    plot(0,0,xlim=c(0.8,ncol(x)+0.2),ylim=c(0.8,nrow(x)+.2),cex=1,xlab='',ylab='',axes=F)
    rect(0.7,0.7,ncol(x)+0.3,nrow(x)+.3,col=rgb(30,30,30,alpha=30,max=255),border=NA)
    axis(1,at=1:ncol(x),toupper(colnames(x)),font=1)
    axis(2,at=1:nrow(x),rownames(x),font=1)
    for(i in 1:ncol(x))points(rep(i,nrow(x)),1:nrow(x),cex=sig.cex[,i]*1.5,pch=15,col=redblue.colors[Value.mat[,i]])
    par(mai=c(0,0,0,0),las=1)
    plot(0,0,cex=0,axes=0,xlab='',ylab='',xlim=c(-1,1),ylim=c(-1,1))
    if(FDR)legend(-1,0.8,legend=paste(c('q>','q<='),c(FDR,FDR),sep=''),pch=15,pt.cex=c(0.3,0.9)*1.5,cex=0.7)
    if(!FDR)legend(-1,0.8,legend=c('P>0.05','P<=0.05'),pch=15,pt.cex=c(0.3,0.9)*1.5,cex=0.7)
    legend(-1,0.4,legend=levels(Value),pch=15,col=redblue.colors,cex=0.7)
}

BarPlot.0 <- function(x,scale=c('OR','Cor'),title.lab=''){
    if(length(title.lab)==0)title.lab=deparse(substitute(x))
    tmp.x=x
    names(x)=toupper(names(x))
    tmp.x=gsub(' .*','',tmp.x)
    x=x[order(tmp.x,decreasing=T)]
    tmp.x=as.numeric(tmp.x)
    tmp.x=tmp.x[order(tmp.x,decreasing=T)]
    tmp=strsplit(x,' ')
    sig.levels=c()
    for(i in tmp){
        if(length(i)==1)sig.levels=c(sig.levels,'') else {
            pp=as.numeric(i[2])
            if(pp<=0.05)sig.levels=c(sig.levels,'*')
            if(pp>0.05 & pp<=0.1)sig.levels=c(sig.levels,'.')
            if(pp>0.1)sig.levels=c(sig.levels,'')
        }
    }
    library(ggplot2)
    df=data.frame(x=tmp.x,l=names(x),s=sig.levels)
    df$l=factor(df$l,levels=df$l[order(df$x,decreasing=T)])
    
    if(scale=='OR'){
        ggplot(df,aes(l,x,label=s))+geom_bar(stat='identity',fill='darkgreen',colour='pink')+geom_text()+theme(axis.text.x=element_text(size=14,angle=90,vjust=0.7),axis.text.y=element_text(size=14),axis.title=element_text(size=0),title=element_text(size=14))+geom_abline(intercept=1,slope=0,colour='lightblue')+scale_y_continuous(breaks=c(0,1,3,6,9,12),labels=c(0,1,3,6,9,12))+labs(title=title.lab)
    }
 
    if(scale=='Cor'){
        cols=rep('lightblue',nrow(df))
        cols[which(df$x<0)]='orange'
        df$cols=cols
        ggplot(df,aes(l,x,label=s))+geom_bar(stat='identity',fill=cols,colour='pink')+geom_text()+theme(axis.text.x=element_text(size=18,angle=45,vjust=0.7),axis.text.y=element_text(size=18),axis.title=element_text(size=0),title=element_text(size=18))+geom_abline(intercept=1,slope=0,colour='lightblue')+labs(title=title.lab)
    }

}

