if (T) {
 
  dir.create("Files")
  dir.create("Figures")
  dir.create("00_origin_datas/GEO",recursive = T)
  dir.create("00_origin_datas/TCGA")
  
}
library(stringr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)


my_mutiviolin=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                       #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                       bw=T,xlab='',ylab='score',title='',size=3,angle = 45, hjust = 1,
                       legend.position='top',fill='group',notch=F){
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(dat.melt,aes(x=type, y=value,fill=Group)) +
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill='Group')+
    geom_violin(trim = F,position=position_dodge(0.8))+  
    scale_fill_manual(values = group_cols)+
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust))
  return(p)
}
plotMutiBar <-plotMutiBar <-function(dat=Age1_compare,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T,color){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) 
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_npg()+theme(legend.position = "bottom")
  pg=pg+ggsci::scale_fill_npg()+scale_fill_manual(values = color, 
                                                  breaks = paste0('R', 1:nrow(dat)), 
                                                  labels = lbr, 
                                                  name = legTitle) 
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") +xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(CHl-Squared p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

wb_beeswarm_plot <- function(dat = NULL,
                             show_compare = T,
                             xlab = 'Groups',
                             ylab = '',
                             method = c('t.test', 'wilcox.test')[1],
                             col = mycolor,
                             leg.pos = c('top','left','right','bottom','none')[1],
                             title = NULL,
                             group = 'Cluster') {
  library(ggbeeswarm)
  colnames(dat) <- c('Cluster', 'Feature')
  
  
  p1 <- ggplot(dat, aes(Cluster, Feature, color = Cluster)) + geom_quasirandom(method = "frowney") +
    ggtitle(title) + scale_color_manual(values = col[1:length(unique(dat$Cluster))]) +
    xlab(xlab) + ylab(ylab) + guides(color=guide_legend(title = group)) + theme_classic() +
    theme(legend.position=leg.pos)
  
  
  if(show_compare){
    uni.group = as.character(unique(dat$Cluster))
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,
                                     method = method,
                                     label= "p.signif", 
                                     step_increase = 0.0)
  }
  return(p1)
}
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin(trim = F)+  
    scale_fill_manual(values = group_cols)+
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}
my_boxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                    fill= "Group",label=c("p.format",'p.signif')[1],notch=F,
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=ARGs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +
    scale_fill_manual(values = group_cols)+   #
    # if(length(names(table(group)))>2){
    #   test_method=''
    # }
    ggpubr::stat_compare_means(aes(group=Group), label = label, method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 
  return(p)
}

my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    # theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","grey"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(
      legend.title = element_blank(),#
      legend.position = leg.pos,
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#
  return(p)
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}

coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
####TCGA-LUAD################
tcga.cli1<-read.delim('00_origin_datas/TCGA/Merge_LUAD_clinical.txt',sep='\t',header = T)
colnames(tcga.cli1)[1:20]
tcga.cli1$number_pack_years_smoked
table(tcga.cli1$tobacco_smoking_history)
tcga.cli1=data.frame(Samples=tcga.cli1$A0_Samples,
                     Age=tcga.cli1$A17_Age,
                     Gender=tcga.cli1$A18_Sex,
                     Smoking_history=tcga.cli1$tobacco_smoking_history,
                     T.stage=tcga.cli1$A3_T,
                     N.stage=tcga.cli1$A4_N,
                     M.stage=tcga.cli1$A5_M,
                     Stage=tcga.cli1$A6_Stage)
tcga.cli1$Samples=paste0(tcga.cli1$Samples,'-01')
rownames(tcga.cli1)=tcga.cli1$Samples
head(tcga.cli1)
table(tcga.cli1$Smoking_history)

table(tcga.cli1$T.stage)
tcga.cli1$T.stage=gsub('[ab]','',tcga.cli1$T.stage)
tcga.cli1$T.stage[tcga.cli1$T.stage=='TX']<-NA

table(tcga.cli1$N.stage)
tcga.cli1$N.stage[tcga.cli1$N.stage=='NX'|tcga.cli1$N.stage=='']<-NA

table(tcga.cli1$M.stage)
tcga.cli1$M.stage=gsub('[ab]','',tcga.cli1$M.stage)
tcga.cli1$M.stage[tcga.cli1$M.stage=='MX'|tcga.cli1$M.stage=='']<-NA

table(tcga.cli1$Stage)
tcga.cli1$Stage=gsub('[AB]','',tcga.cli1$Stage)
tcga.cli1$Stage[tcga.cli1$Stage=='']<-NA
tcga.cli1$Stage=gsub('Stage ','',tcga.cli1$Stage)


tcga.pancancer.cli=read.xlsx
head(tcga.pancancer.cli)
tcga.cli2=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='LUAD'),]
head(tcga.cli2)
tcga.cli2=data.frame(Samples=paste0(tcga.cli2$bcr_patient_barcode,'-01'),
                     tcga.cli2[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
head(tcga.cli2)
tcga.cli2$OS.time
tcga.cli2=tcga.cli2 %>% drop_na(OS.time)
####
tcga.cli2=tcga.cli2[tcga.cli2$OS.time>0,]
dim(tcga.cli2)

tcga.cli=merge(tcga.cli1,tcga.cli2,by='Samples')
rownames(tcga.cli)=tcga.cli$Samples
tcga.cli=as.data.frame(tcga.cli)
fivenum(as.numeric(tcga.cli$Age))
tcga.cli$Age1=ifelse(as.numeric(tcga.cli$Age)>66,'>66','<=66')
dim(tcga.cli)
# 509  17
head(tcga.cli)


#
tcga_data<-read.delim('00_origin_datas/TCGA/Merge_RNA_seq_FPKM _LUAD.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga_data[1:4,1:4]
table(substr(colnames(tcga_data),14,15))
tcga_data <- exp_ensg2symbol(tcga_data)

sample_T=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==1)]#
sample_N=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==11)]#
length(sample_N)
length(sample_T)
tcga_type=data.frame(Samples=c(sample_T,sample_N),type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$type)
# Normal  Tumor 
# 59    513

genecode=read.delim('data/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]


range(tcga_data)
tcga.exp.all=log2(tcga_data[intersect(rownames(tcga_data),mrna_genecode$SYMBOL),tcga_type$Samples]+1)
range(tcga.exp.all)
tcga.exp=tcga.exp.all[,intersect(tcga.cli$Samples,sample_T)]
dim(tcga.exp)
# 19503   500
tcga.cli=tcga.cli[intersect(tcga.cli$Samples,sample_T),]
dim(tcga.cli)
saveRDS(tcga.exp.all,file = '00_pre_datas/TCGA/LUAD_FPKM_tcga.exp.all.RDS')
saveRDS(tcga.exp,file = '00_pre_datas/TCGA/LUAD_FPKM_tcga.exp.RDS')
########GSE31210#############
GSE31210 <- getGEOExpData('GSE31210')
saveRDS(GSE31210,file = "00_origin_datas/GEO/GSE31210.RDS")
#load('00_origin_datas/GEO/GSE31210.RData')

GSE31210.cli=GSE31210$Sample
GSE31210.cli=data.frame(Samples=GSE31210.cli$Acc,
                        Age=GSE31210.cli$`age (years)`,
                        Gender=GSE31210.cli$gender,
                        Status=GSE31210.cli$death,
                        OS.time=GSE31210.cli$`days before death/censor`)
rownames(GSE31210.cli)=GSE31210.cli$Samples
table(GSE31210.cli$Status)
GSE31210.cli=GSE31210.cli[which(GSE31210.cli$Status!='NULL'),]
GSE31210.cli$OS=ifelse(GSE31210.cli$Status=='alive',0,1)
GSE31210.cli <- GSE31210.cli[GSE31210.cli$OS.time>0,]
range(GSE31210.cli$OS.time)


GSE31210.exp=GSE31210$Exp$GPL570_54675_Data_col1
range(GSE31210.exp)
GSE31210.exp=log2(GSE31210.exp+1)
GSE31210.exp=exp_probe2symbol_v2(GSE31210.exp,GPL ='GPL570' )
range(GSE31210.exp)
dim(GSE31210.exp)
dim(GSE31210.cli)
GSE31210.exp=GSE31210.exp[,GSE31210.cli$Samples]
dim(GSE31210.exp)
#20549   226



####01.##########
dir.create('01_cluster')
######1.1 #######
NAM.gene <- read.table("01_cluster/NAM NAMbolism-related genes_PMID 36969219.txt",sep = "\t")
length(NAM.gene)
#42
NAM.gene <- NAM.gene$V1
tcga.NAM.ssGSEA=ssGSEAScore_by_genes(gene.exp = tcga.exp.all,genes =NAM.gene)
saveRDS(tcga.NAM.ssGSEA,file = "01_cluster/tcga.NAM.ssGSEA.RDS")
tcga.NAM.ssGSEA <- readRDS("01_cluster/tcga.NAM.ssGSEA.RDS")
range(tcga.NAM.ssGSEA)

df=crbind2DataFrame(t(tcga.NAM.ssGSEA))

# fig1a=my_violin(dat=t(tcga.NAM.ssGSEA),group=tcga_type$type)
# ggsave("01_cluster/fig1a.pdf",width = 6,height = 6)
######1.2#########
NAM.cox=cox_batch(tcga.exp[NAM.gene,tcga.cli$Samples],time = tcga.cli$OS.time,event = tcga.cli$OS)
NAM.cox[order(NAM.cox$p.value),]
table(NAM.cox$p.value<0.05)
#FALSE  TRUE 
#36     6 
NAM.cox.fit=NAM.cox[NAM.cox$p.value<0.05,]
dim(NAM.cox.fit)
NAM.cox.fit$Type=ifelse(NAM.cox.fit$HR>1,'Risk','Portect')
sig_fit_gene=rownames(NAM.cox.fit)
length(sig_fit_gene)#6
writeMatrix(NAM.cox.fit,outpath = '01_cluster/sig_cox_res.txt')
writeMatrix(NAM.cox.fit,outpath = 'Files/sig_cox_res.txt')


########
pdf('01_cluster/sig_cox_bioForest.pdf',height = 8,width = 8)
bioForest(rt = NAM.cox.fit,col=c('blue','red'))
dev.off()


#####1.3 #######################
library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[1]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[2]
consen_gene=rownames(NAM.cox.fit)
length(consen_gene)#6
tcga_consen_data=as.matrix(tcga.exp[consen_gene,tcga.cli$Samples])
tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   #11,2  21,2
#tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   #
#tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, mean))#
# tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, median))  #
#tcga_consen_data=as.dist(1-cor(tcga_consen_data,method = 'pearson'))
tcga_consen_data=as.matrix(tcga_consen_data)
dim(tcga_consen_data)
tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "TCGA_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = T
                                           , seed = 123456)

k=2
cluster.color=pal_nejm()(8)[c(5,8)]

tcga.subtype <- data.frame(Samples = names(tcga_clust_subtype[[k]]$consensusClass),
                           Cluster=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Cluster=paste0('C',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
write.csv(tcga.subtype[order(tcga.subtype$Cluster),],'01_cluster/TCGA_subtype.csv',row.names = F)
tcga.subtype.cli=merge(tcga.subtype,tcga.cli,by='Samples')
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples


tcga.subtype.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Cluster,
                                        data = tcga.subtype.cli),
                           data=tcga.subtype.cli,
                           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                           title='TCGA-LUAD',#ggtheme=custom_theme(),
                           linetype = c("solid", "dashed","strata")[1],
                           palette = cluster.color,
                           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                           # legend = c(0.8,0.75), # 
                           legend.title = "")
tcga.subtype.km=mg_merge_plot(tcga.subtype.km$plot,tcga.subtype.km$table,ncol = 1,nrow = 2,heights = c(3,1),align = 'v')
tcga.subtype.km
ggsave("01_cluster/tcga.subtype.km.pdf",width = 6,height = 6)




######1.#######
GSE31210_conse<-GSE31210.exp[intersect(consen_gene,rownames(GSE31210.exp)),]
dim(GSE31210_conse)
GSE31210_conse=t(scale(t(GSE31210_conse),scale = F))

GSE31210_subtype <- ConsensusClusterPlus(GSE31210_conse
                                         , maxK = 10, reps = 500, pItem = 0.8
                                         , pFeature = 1
                                         , title = "GSE31320_subtype"
                                         , clusterAlg = clusterAlg_name
                                         , distance = distance_name
                                         # , innerLinkage = 'ward.D2'
                                         , plot = "pdf"
                                         , writeTable = T
                                         , seed = 123456)
GSE31210.subtype <- data.frame(Samples = names(GSE31210_subtype[[k]]$consensusClass),Cluster=GSE31210_subtype[[k]]$consensusClass)
GSE31210.subtype$Cluster=paste0('C',GSE31210.subtype$Cluster)
rownames(GSE31210.subtype)=GSE31210.subtype$Samples
GSE31210.subtype$Cluster=gsub('C1','IS2',GSE31210.subtype$Cluster)
GSE31210.subtype$Cluster=gsub('C2','IS1',GSE31210.subtype$Cluster)
GSE31210.subtype$Cluster=gsub('IS','C',GSE31210.subtype$Cluster)
table(GSE31210.subtype$Cluster)
writeMatrix(GSE31210.subtype,row=T,header=T,outpath = '04_model/GSE31210.subtype.txt')
#GSE31210.subtype=readMatrix(inpath = 'GSE31210.subtype.txt',row=T,header=T)
#KM
GSE31210.subtype.cli <- cbind(GSE31210.cli, GSE31210.subtype[GSE31210.cli$Samples, ])
GSE31210.subtype.cli=crbind2DataFrame(GSE31210.subtype.cli)
write.table(GSE31210.subtype.cli,'04_model/GSE31210.subtype.cli.txt',quote = F,row.names = F,sep='\t')
GSE31210.subtype.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Cluster,
                                            data = GSE31210.subtype.cli),
                               data=GSE31210.subtype.cli,
                               conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                               title='GSE31210',#ggtheme=custom_theme(),
                               linetype = c("solid", "dashed","strata")[1],
                               palette = cluster.color,
                               legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                               # legend = c(0.8,0.75), # 
                               legend.title = "")
GSE31210.subtype.km=mg_merge_plot(GSE31210.subtype.km$plot,GSE31210.subtype.km$table,ncol = 1,nrow = 2,heights = c(3,1),align = 'v')
GSE31210.subtype.km
ggsave("01_cluster/GSE31210.subtype.km.pdf",width = 6,height = 6)



# 02  ######################
dir.create('02_cluster_Clinical_TMB')
tcga.mut.dat <- getTCGAMAFByCode('LUAD')
tcga.mut.dat <- as.data.frame(tcga.mut.dat@data)
tcga.mut.dat <- tcga.mut.dat[, c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]
table(tcga.mut.dat$Variant_Classification)
tcga.mut.dat$Variant_Classification <- 1

#
tcga.mut.dat <- reshape2::dcast(data = tcga.mut.dat, Hugo_Symbol ~ Tumor_Sample_Barcode)
class(tcga.mut.dat)
rownames(tcga.mut.dat) <- tcga.mut.dat$Hugo_Symbol
tcga.mut.dat <- tcga.mut.dat[, -1]


#
colnames(tcga.mut.dat) <- paste0(colnames(tcga.mut.dat), '-01')
mut.samples <- intersect(colnames(tcga.mut.dat), tcga.subtype.cli$Samples)
rownames(tcga.subtype.cli)
tcga.mut.dat <- tcga.mut.dat[, mut.samples]
tcga_mut_cli <- tcga.subtype.cli[mut.samples, ]
tcga.mut.dat <- ifelse(tcga.mut.dat == 0, 0, 1)


# mut.genes <- openxlsx::read.xlsx('00_origin_datas/cancer.driver.genes.PMC4160307.xlsx', sheet = 1)[, 1]
# mut.genes <- intersect(mut.genes, rownames(tcga.mut.dat))
# tcga.mut.dat<- tcga.mut.dat[mut.genes, ]

gene.freq <- as.data.frame(rowSums(tcga.mut.dat))
head(gene.freq)
colnames(gene.freq) <- 'Freq'
gene.freq$Genes <- rownames(gene.freq)
gene.freq <- arrange(gene.freq, desc(Freq))
write.csv(gene.freq,
          file = '02_cluster_Clinical_TMB/tcga.subtype.mut.gene.csv',
          quote = F)

mut.genes <- gene.freq$Genes[gene.freq$Freq > 2]
length(mut.genes)
#11776
mut.res <- data.frame(`NAM_C1` = NA,
                      `NAM_C2` = NA)


mut.p <- c()
tmp.ge <- c()
for (ge in mut.genes) {
  print(ge)
  tmp <- table(tcga.mut.dat[ge, ], tcga_mut_cli$Cluster)
  if (dim(tmp)[1] > 1) {
    pvalue <- fisher.test(tmp)
    mut.p <- c(mut.p, pvalue$p.value)
    mut.res <- rbind(mut.res, tmp[2, ])
    tmp.ge <- c(tmp.ge, ge)
  }
}
mut.res <- na.omit(mut.res)
rownames(mut.res) <- tmp.ge
class(mut.res)
mut.res$P.value <- mut.p
dim(mut.res)
head(mut.res)
write.csv(mut.res,
          file = '02_cluster_Clinical_TMB/mut.res.csv',
          quote = F)

table(mut.res$P.value < 0.05)
# FALSE  TRUE 
# 11407   369 
mut.res.filtered <- mut.res[mut.res$P.value < 0.05, ]
head(mut.res.filtered)
dim(mut.res.filtered)
# write.csv(mut.res.filtered,
#           file = 'results/tcga.mut.dat.csv')
write.csv(mut.res.filtered,
          file = '02_cluster_Clinical_TMB/tcga.genes.csv')

mut.plot.dat <- tcga.mut.dat[rownames(mut.res.filtered)[1:20], ]
mut.plot.dat <- ifelse(mut.plot.dat == 1, 'Mutant', NA)


#####################
# snv.genes.filter1 <- intersect(gene.freq$Genes, snv.genes.filter)
# 
# mut.plot.dat <- tcga.mut.dat[snv.genes.filter1, ]
# mut.plot.dat[mut.plot.dat == 'WildType'] <- NA

library(scales)
library(ggsci)
library(ComplexHeatmap)
alter_graphic <- function (graphic = c("rect", "point"), width = 1, 
                           height = 1, horiz_margin = unit(1, "pt"),
                           vertical_margin = unit(1, "pt"), 
                           fill = "red", col = NA, pch = 16, 
                           ...) 
{
  graphic = match.arg(graphic)[1]
  if (graphic == "rect") {
    if (!is.numeric(width)) {
      stop_wrap("`width` should be nummeric.")
    }
    if (!is.numeric(height)) {
      stop_wrap("`height` should be nummeric.")
    }
    if (width != 1) {
      if (missing(horiz_margin)) {
        horiz_margin = unit(0, "pt")
      }
    }
    if (height != 1) {
      if (missing(vertical_margin)) {
        vertical_margin = unit(0, "pt")
      }
    }
    fun = function(x, y, w, h) {
      w = w * width
      h = h * height
      grid.rect(x, y, w - horiz_margin * 2, h - vertical_margin * 
                  2, gp = gpar(fill = fill, col = col, ...))
    }
  }
  else if (graphic == "point") {
    fun = function(x, y, w, h) {
      grid.points(x, y, pch = pch, gp = gpar(fill = fill, 
                                             col = col, ...))
    }
  }
  return(fun)
}

#### 
col = c("Mutant" = '#006633')

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  Mutant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["Mutant"], col = NA))
  }
)

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),	
  Mutant = alter_graphic("rect", fill = col["Mutant"])
)
heatmap_legend_param = list(title = "Alterations", 
                            at = c("Mutant"), 
                            labels = c("Mutant"))


tcga_mut_cli <- tcga.subtype.cli[mut.samples, ]
tcga_mut_cli <- cbind(tcga_mut_cli,NAM_SCOER=as.matrix(tcga.NAM.ssGSEA[,tcga_mut_cli$Samples ]))
table(tcga_mut_cli$Cluster)
###
color_cluster = cluster.color
names(color_cluster) = c('C1', 'C2')
mycolor=c4a('brewer.dark2',7)
###
library(dplyr)
colnames(tcga_mut_cli)
# tcga_mut_cli <- arrange(tcga_mut_cli, Cluster, `cancer stem cell`)

Gender_color <- mycolor[1:2]
names(Gender_color) <- na.omit(unique(tcga_mut_cli$Gender))

A3_T_color <- mycolor[1:4]
names(A3_T_color) <- na.omit(unique(tcga_mut_cli$T.stage))

A4_N_color <- mycolor[1:4]
names(A4_N_color) <- na.omit(unique(tcga_mut_cli$N.stage))

A5_M_color <- mycolor[1:2]
names(A5_M_color) <- na.omit(unique(tcga_mut_cli$M.stage))

Stage_color <- mycolor[1:4]
names(Stage_color) <- na.omit(unique(tcga_mut_cli$Stage))


# Grade_color <- mycolor[1:3]
# names(Grade_color) <- na.omit(unique(tcga_mut_cli$Grade))

tcga_mut_cli <- crbind2DataFrame(tcga_mut_cli)
tcga_NAM_pheatmap <- oncoPrint(as.matrix(mut.plot.dat[, tcga_mut_cli$Samples]),
                               # row_order = rownames(c1.mut.plot.dat),
                               # column_order = tcga_mut_cli$Samples,
                               alter_fun = alter_fun, 
                               col = col, 
                               column_title = "TCGA",
                               heatmap_legend_param = heatmap_legend_param,
                               pct_digits = 2,
                               column_split = factor(tcga_mut_cli$Cluster),
                               top_annotation = HeatmapAnnotation(NAM.Score=anno_barplot(tcga_mut_cli$NAM_SCOER,
                                                                                         border = F,gp = gpar(fill = "black"),bar_width = 1, height = unit(3, "cm")),
                                                                  df=tcga_mut_cli[, c("Cluster", "Gender", "Age", "T.stage",  "N.stage","M.stage", "Stage")],
                                                                  col=list("Cluster"=color_cluster,
                                                                           "Gender"=Gender_color,
                                                                           "T.stage"=A3_T_color, 
                                                                           "N.stage"=A4_N_color,
                                                                           "M.stage"=A5_M_color, 
                                                                           "Stage"=Stage_color),
                                                                  show_annotation_name = TRUE,
                                                                  gap=unit(1, "mm"),
                                                                  na_col="grey"),
                               pct_side = "right", row_names_side = "left")
tcga_NAM_pheatmap
dev.off()
pdf('02_cluster_Clinical_TMB/tcga_mRNAsi_mut_pheatmap.pdf', width = 12, height = 8)
tcga_NAM_pheatmap
dev.off()



#######   ##############
# TCGA

Age1_compare <- table(tcga.subtype.cli$Age1, tcga.subtype.cli$Cluster)
Age1_compare
Age1_compare_bar <- plotMutiBar(Age1_compare, legTitle = 'Age', showValue = T)
Age1_compare_bar

Gender_compare <- table(tcga.subtype.cli$Gender, tcga.subtype.cli$Cluster)
Gender_compare
Gender_compare_bar <- plotMutiBar(Gender_compare, legTitle = 'Gender', showValue = T)
Gender_compare_bar

A3_T_compare <- table(tcga.subtype.cli$T.stage, tcga.subtype.cli$Cluster)
A3_T_compare
A3_T_compare_bar <- plotMutiBar(A3_T_compare, legTitle = 'T Stage', showValue = T)
A3_T_compare_bar

A4_N_compare <- table(tcga.subtype.cli$N.stage, tcga.subtype.cli$Cluster)
A4_N_compare
A4_N_compare_bar <- plotMutiBar(A4_N_compare, legTitle = 'N Stage', showValue = T)
A4_N_compare_bar

A5_M_compare <- table(tcga.subtype.cli$M.stage, tcga.subtype.cli$Cluster)
A5_M_compare
A5_M_compare_bar <- plotMutiBar(A5_M_compare, legTitle = 'M Stage', showValue = T)
A5_M_compare_bar

A6_Stage_compare <- table(tcga.subtype.cli$Stage, tcga.subtype.cli$Cluster)
A6_Stage_compare
A6_Stage_compare_bar <- plotMutiBar(A6_Stage_compare, legTitle = 'Stage', showValue = T)
A6_Stage_compare_bar

tcga_cluster_bar <- cowplot::plot_grid(
                                       Age1_compare_bar,
                                       Gender_compare_bar,
                                       A3_T_compare_bar,
                                       A4_N_compare_bar,
                                       A5_M_compare_bar,
                                       A6_Stage_compare_bar,
                                       align = 'hv',
                                       ncol = 4,
                                       labels = LETTERS[1:9])
tcga_cluster_bar
ggsave(plot = tcga_cluster_bar,
       filename = '02_cluster_Clinical_TMB/tcga_cluster_bar.pdf',
       width = 14, height = 10)


tcga_cluster_bar1 <- cowplot::plot_grid(A5_M_compare_bar,
                                        A3_T_compare_bar,
                                        A4_N_compare_bar,
                                        A6_Stage_compare_bar,
                                        align = 'hv',
                                        ncol = 4,
                                        labels = LETTERS[1:4])
tcga_cluster_bar1
ggsave(plot = tcga_cluster_bar1,
       filename = '02_cluster_Clinical_TMB/tcga_cluster_bar1.pdf',
       width = 16, height = 5)



####03###########
dir.create('03_diff_genes')
tcga.limma=mg_limma_DEG(exp =tcga.exp[,tcga.subtype.cli$Samples],
                        group=tcga.subtype.cli$Cluster,ulab='C1',dlab = 'C2')
tcga.limma$Summary

#######
tcga.cluster.degs=tcga.limma$DEG
tcga.cluster.degs=tcga.cluster.degs[abs(tcga.cluster.degs$logFC)>log2(1.5) & tcga.cluster.degs$adj.P.Val<0.05,]
dim(tcga.cluster.degs)
write.csv(tcga.cluster.degs,'03_diff_genes/tcga.cluster.degs.csv')

##########

degs_dat=tcga.limma$DEG
degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                            ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
p_cutoff = 0.05
fc_cutoff =log2(1.5)
col=c("red","blue","grey")
ylab='-log10 (adj.PVal)'
xlab='log2 (FoldChange)'
leg.pos='right'

library(ggbreak)
library(ggplot2)
library(ggprism)
p_cutoff = 0.05
fc_cutoff =log2(1.5)
col=c("red","blue","grey")
ylab='-log10 (adj.PVal)'
xlab='log2 (FoldChange)'
leg.pos='right'
plot1 <- ggplot(degs_dat, aes(x=logFC, y=-log10(adj.P.Val), color=type)) +
  geom_point() +
  scale_color_manual(values=col) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab(ylab) +
  xlab(xlab) +
  geom_vline(xintercept=c(-fc_cutoff,fc_cutoff), lty=3, col="black", lwd=0.5) +
  geom_hline(yintercept = -log10(p_cutoff), lty=3, col="black", lwd=0.5) +
  #coord_cartesian(ylim=c(0, 25)) # 
  scale_y_break(c(25,100),#
                space = 0.3,#
                scales = 0.5)+#
  theme_prism(palette = "black_and_white",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16,
              base_line_size = 0.8,
              axis_text_angle = 0)
ggsave('03_diff_genes/fig3a.pdf',plot1,height = 6,width = 6)
######3.2 ########
tcga.geneList=getGeneFC(gene.exp=tcga.exp[,tcga.subtype.cli$Samples],group=tcga.subtype.cli$Cluster
                        ,ulab='C1',dlab = 'C2')

h.all.gmt<-read.gmt("data/h.all.v2023.1.Hs.entrez.gmt")
tcga.hallmark.gsea<-GSEA(tcga.geneList,TERM2GENE = h.all.gmt,seed=T)
library(enrichplot)
library(ggplot2)
gsea.p=clusterProfiler::dotplot(tcga.hallmark.gsea,split=".sign",showCategory=nrow(tcga.hallmark.gsea@result),
                                title='C1 vs C2',font.size=10)+facet_grid(~.sign)+
  theme(text = element_text(family = 'Times'))+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
write.csv(tcga.hallmark.gsea@result,file = '03_diff_genes/cluster_GSEA.csv',quote = F)


ggsave('03_diff_genes/fig3b.pdf',gsea.p,height = 6,width = 8)
#######04.########
dir.create('04_model')

sig.gene.cox=cox_batch(dat = tcga.exp[intersect(rownames(tcga.cluster.degs),rownames(tcga.exp)),tcga.cli$Samples],
                        time = tcga.cli$OS.time,event = tcga.cli$OS)
sig.gene.cox
write.csv(sig.gene.cox,'04_model/sig.cox.csv')


table(sig.gene.cox$p.value<0.05)
# FALSE  TRUE 
# 151    88 
pre.genes=rownames(sig.gene.cox[sig.gene.cox$p.value<0.05,])
length(pre.genes)#88
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

######LASSO
tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[,-c(1:2)],
                               os = tcga_model_data$OS,
                               os.time = tcga_model_data$OS.time)
length(tcga.lasso$lasso.gene)#4
tcga.lasso$plot

#####
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(tcga.lasso$lasso.gene,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
#cox=step(cox)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)
####
######geneCoefficients#####
lan.dataframe <- as.data.frame(lan)

lan.dataframe$gene <- rownames(lan.dataframe) 
lan.dataframe$gene <- factor(lan.dataframe$gene,levels = rownames(lan.dataframe)[order(lan.dataframe$lan)])
# 
lan.dataframe$color_group <- ifelse(lan.dataframe$lan > 0, "Positive", "Negative")
library(ggplot2)

# 


# 
p <- ggplot(lan.dataframe, aes(x=gene, y=lan,fill=color_group)) +
  geom_bar(stat="identity") +
  xlab("Gene Name") +
  ylab("Coefficient") +
  ggtitle("Gene Coefficients") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#99CC66", "Negative" = "#CC9999")) +
  theme_bw()+
  guides(fill=FALSE)
p1 <- p+geom_text(aes(label=sprintf("%.3f", lan)), hjust=-0.2, size=3, color="black")
ggsave('04_module/gene_ Coefficients.pdf',height = 6,width = 9)


####
risktype.col=c('#CCCC00',"#CC6699")

risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.subtype.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=risk.tcga)
#######
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(tcga.risktype.cli$Riskscore),'High','Low')

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,3,5))
tcga.roc
tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                  data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = T,pval = T,risk.table = T, 
                      fun = "pct",size = 1,surv.median.line = 'hv',
                      title='TCGA-LUAD',legend.title='Risktype',
                      legend.labs = c('High','Low'),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,
                      ylab='Overall Survival(OS)',
                      legend=c(0.85,0.8),#
                      ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(3,1),align = 'v')

tcga.km.OS

tcga.risktype.cli$Status=ifelse(tcga.risktype.cli$OS==0,'Alive','Dead')
tcga.model.p=my_riskplot(cli_dat = tcga.risktype.cli,cols =risktype.col,xlab = 'sample',
            a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = median(tcga.risktype.cli$Riskscore),labs = '')


tcga.km.DSS=ggsurvplot(fit=survfit(Surv(DSS.time/365, DSS) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',
                       title='TCGA-LUAD',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,
                       ylab='Disease-Specific Survival(DSS)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.DSS=mg_merge_plot(tcga.km.DSS$plot,tcga.km.DSS$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.DSS

tcga.km.DFI=ggsurvplot(fit=survfit(Surv(DFI.time/365, DFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',
                       title='TCGA-LUAD',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,
                       ylab='Disease-Free interval(DFI)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.DFI=mg_merge_plot(tcga.km.DFI$plot,tcga.km.DFI$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.DFI

tcga.km.PFI=ggsurvplot(fit=survfit(Surv(PFI.time/365, PFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',title='TCGA-LUAD',
                       legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       palette = risktype.col,
                       ylab='Progression-Free Interval(PFI)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.PFI=mg_merge_plot(tcga.km.PFI$plot,tcga.km.PFI$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.PFI


########
model.gene.df=data.frame(tcga.cli[,c('OS','OS.time')],t(tcga.exp[names(lan),tcga.cli$Samples]))
head(model.gene.df)
module.gene.km=list()
for (i in 1:length(names(lan))) {
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km[[i]]=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                             data = model.gene.df1),
                                 data=model.gene.df1,
                                 conf.int = F,pval = T,risk.table = T, 
                                 fun = "pct",size = 1,surv.median.line = 'hv',
                                 title='TCGA',legend.title=names(lan)[i],
                                 # legend.labs = c('High','Low'),
                                 linetype = c("solid", "dashed","strata")[1],
                                 palette = risktype.col,
                                 ylab='Overall Survival(OS)',
                                 legend=c(0.8,0.8),#
                                 ggtheme = theme_bw(base_size = 12))
  module.gene.km[[i]]=module.gene.km[[i]]$plot
}
tcga.module.km=mg_merge_plot(module.gene.km,ncol=4,nrow=1)



########4.2 ##########
GSE31210_model_data <- data.frame(GSE31210.cli[, c("OS.time", "OS")],
                             t(GSE31210.exp[intersect(names(lan),rownames(GSE31210.exp)), GSE31210.cli$Samples]))
colnames(GSE31210_model_data) <- gsub('-', '_', colnames(GSE31210_model_data))

risk.GSE31210=as.numeric(lan%*%as.matrix(t(GSE31210_model_data[GSE31210.cli$Samples,names(lan)])))

GSE31210.risktype.cli=data.frame(GSE31210.cli,Riskscore=risk.GSE31210)
GSE31210.risktype.cli$Risktype=ifelse(GSE31210.risktype.cli$Riskscore>median(risk.GSE31210),'High','Low')
GSE31210.roc=ggplotTimeROC(GSE31210.risktype.cli$OS.time,
                           GSE31210.risktype.cli$OS,
                           GSE31210.risktype.cli$Riskscore,mks = c(1,2,3))
GSE31210.roc
GSE31210.km=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                   data = GSE31210.risktype.cli),
                       data=GSE31210.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='GSE31210',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ylab='Overall Survival(OS)',
                       legend=c(0.85,0.25),#
                       ggtheme = theme_bw(base_size = 12))
GSE31210.km=mg_merge_plot(GSE31210.km$plot,GSE31210.km$table,nrow=2,heights = c(3,1),align = 'v')
GSE31210.km

gse31210.model.p=my_riskplot(cli_dat = GSE31210.risktype.cli,cols = risktype.col,xlab = 'sample',
            a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = median(GSE31210.risktype.cli$Riskscore),labs = '')



model.gene.df=data.frame(GSE31210.cli[,c('OS','OS.time')],t(GSE31210.exp[names(lan),GSE31210.cli$Samples]))
head(model.gene.df)
module.gene.km=list()
for (i in 1:length(names(lan))) {
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km[[i]]=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                             data = model.gene.df1),
                                 data=model.gene.df1,
                                 conf.int = F,pval = T,risk.table = T, 
                                 fun = "pct",size = 1,surv.median.line = 'hv',
                                 title='GSE31210',legend.title=names(lan)[i],
                                 # legend.labs = c('High','Low'),
                                 linetype = c("solid", "dashed","strata")[1],
                                 palette = risktype.col,
                                 ylab='Overall Survival(OS)',
                                 legend=c(0.7,0.35),#
                                 ggtheme = theme_bw(base_size = 12))
  module.gene.km[[i]]=mg_merge_plot(module.gene.km[[i]]$plot,module.gene.km[[i]]$table,nrow=2,heights = c(3,1),align = 'v')
}
GSE31210.module.km=mg_merge_plot(module.gene.km,ncol=4,nrow=1)

###########
fig4_1=mg_merge_plot(mg_merge_plot(tcga.lasso$plot,p1,labels = c('','C')),
                   mg_merge_plot(tcga.model.p,tcga.roc,tcga.km.OS,ncol=3,labels = LETTERS[4:6]),
                   mg_merge_plot(tcga.km.DFI,tcga.km.DSS,tcga.km.PFI,ncol=3,labels = LETTERS[7:9]),nrow=3)

fig4_2=mg_merge_plot(mg_merge_plot(gse31210.model.p,GSE31210.roc,GSE31210.km,ncol=3,labels = LETTERS[1:3]),
                     GSE31210.module.km,nrow=2,labels = c('','D'))
ggsave('04_model/Fig4-1.pdf',fig4_1,height = 18,width = 20)
ggsave('04_model/figs4-2.pdf',fig4_2,height = 12,width = 20)


####05.##########
dir.create('05_Risktype.immune')
###estimate 
tcga.estimate <- immu_estimate(tcga.exp)
fig5a=my_mutiboxplot(dat = tcga.estimate[tcga.risktype.cli$Samples,],
                     group = tcga.risktype.cli$Risktype,
                     group_cols = risktype.col,legend.position = 'none',
                     ylab = 'score')+ggtitle('ESTIMATE')
saveRDS(tcga.estimate,file ='05_Risktype.immune/tcga.estimate.RDS')
###timer
tcga.timer <- immu_timer(tcga.exp)
fig5b=my_mutiboxplot(dat = tcga.timer[tcga.risktype.cli$Samples,],
                     group = tcga.risktype.cli$Risktype,
                     group_cols = risktype.col,legend.position = 'none',
                     ylab = 'score')+ggtitle('TIMER')
saveRDS(tcga.timer,file ='05_Risktype.immune/tcga.timer.RDS')
###ssgsea
tcga.immu.ssgsea=immu_ssgsea(tcga.exp,isTCGA = T)
saveRDS(tcga.immu.ssgsea,file ='05_Risktype.immune/tcga.immu.ssgsea.RDS')
tcga.immu.ssgsea <- readRDS("05_Risktype.immune/tcga.immu.ssgsea.RDS")
fig5c=my_mutiboxplot(dat = tcga.immu.ssgsea[tcga.risktype.cli$Samples,],
                     group = tcga.risktype.cli$Risktype,
                     group_cols = risktype.col,
                     ylab = 'score')+ggtitle('ssGSEA')
###mcp
tcga.mcp <- immu_MCPcounter(tcga.exp)
fig5d=my_mutiboxplot(dat = tcga.mcp[tcga.risktype.cli$Samples,],
                     group = tcga.risktype.cli$Risktype,
                     group_cols = risktype.col,legend.position = 'none',
                     ylab = 'score')+ggtitle('MCPCounter')
saveRDS(tcga.mcp,file ='05_Risktype.immune/tcga.mcp.RDS')
##estimate
tcga.estimate <- immu_estimate(tcga.exp)
fig5ab=mg_merge_plot(fig5a,fig5b,labels = c('A','B'),ncol=2,widths = c(1,2))
fig5abcd=mg_merge_plot(fig5ab,
                      fig5c,fig5d,labels = c('','C',"D"),align = 'v', ncol =1,nrow = 3,heights = c(1,1.5,1.3)
                     )
savePDF('05_Risktype.immune/Fig5ABCD.pdf',fig5abcd,height = 15,width =15)



####06.############
dir.create('06_TIDE')

########TIDE
# tcga_tide_dat <- t(scale(t(tcga.exp),scale = F))
# dim(tcga_tide_dat)
# write.table(tcga_tide_dat,file = '06_TIDE/tcga_tide_dat.txt',quote = F, sep = '\t')

tcga_tide_res<-read.csv('06_TIDE/LUAD_TIDE_res.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)
tcga_tide_merge=cbind(tcga_tide_res[tcga.risktype.cli$Samples,],tcga.risktype.cli)
tide_sel=c('TIDE','IFNG','Exclusion','Dysfunction','MDSC')

tide.p1=my_violin(dat = tcga_tide_merge$TIDE,group = tcga_tide_merge$Risktype,
          ylab = 'TIDE',fill = 'Risktype',group_cols = risktype.col)
tide.p1

tide.p2=plotMutiBar(table(tcga_tide_merge$Responder,tcga_tide_merge$Risktype),legTitle = 'Responder',color = mycolor[c(1,4)])
tide.p2

tide.p=mg_merge_plot(tide.p1,tide.p2,ncol =2,labels = c('A','B'))


#########IC50
load('06_TIDE/tcga_durg_ic50.RData')
tcga_durg_ic50_res[1:5,1:5]
colnames(tcga_durg_ic50_res)[1]='Cisplatin'
colnames(tcga_durg_ic50_res)=gsub('-','_',colnames(tcga_durg_ic50_res))
colnames(tcga_durg_ic50_res)=gsub(' ','_',colnames(tcga_durg_ic50_res))
colnames(tcga_durg_ic50_res)
dim(tcga_durg_ic50_res)

###
tcga.ic50.diff<-diff_pathway(dat=t(tcga_durg_ic50_res[tcga.risktype.cli$Samples,]),group=tcga.risktype.cli$Risktype)
table(tcga.ic50.diff$p.value<0.01)
# FALSE  TRUE 
# 26    50 
###
library(ggcorrplot)
library(psych)
IC50_RS_cor <- corr.test(x =tcga.risktype.cli$Riskscore,
                         y = tcga_durg_ic50_res[tcga.risktype.cli$Samples,],
                         method = "spearman",adjust = "BH",ci = F)


IC50_RS_cor_res=data.frame(drugs=colnames(tcga_durg_ic50_res))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
IC50_RS_cor_res <- data.frame(IC50_RS_cor_res,diff.pvalue=tcga.ic50.diff$p.value)
head(IC50_RS_cor_res)
table((IC50_RS_cor_res$p.adj<0.05)&(IC50_RS_cor_res$cor<(0.25)&(IC50_RS_cor_res$diff.pvalue<(0.05))))
# FALSE  TRUE 
#  30    46 
IC50_RS_cor_res <- IC50_RS_cor_res[(IC50_RS_cor_res$p.adj<0.05)&(IC50_RS_cor_res$cor<(0.25)),]
IC50_RS_cor_res <-IC50_RS_cor_res [order(IC50_RS_cor_res$cor,decreasing = F),]
dim(IC50_RS_cor_res)
###
colnames(tcga_durg_ic50_res)[c(1,2,7,15,21,22,23,24,25,28)]#
drug.df=data.frame(tcga_durg_ic50_res[tcga.risktype.cli$Samples,c(1,2,7,15,21,22,23,24,25,28)],Risktype=tcga.risktype.cli$Risktype)
drug.df=melt(drug.df)
head(drug.df)
pdf('06_TIDE/Fig6c.pdf',height = 8,width = 12)
dodge_width <-1.5

ic50=ggplot(drug.df, aes(x=variable, y=value, fill=Risktype)) +
  geom_violin(position = position_dodge(dodge_width), trim = FALSE,show.legend = T) +
  geom_boxplot(width = 0.3, position = position_dodge(dodge_width), outlier.shape = NA,show.legend = F) +
  scale_fill_manual(values = risktype.col) +
  facet_wrap(~variable, scales = 'free', nrow = 2, ncol = 5) +
  ggpubr::stat_compare_means(aes(group=Risktype), label = 'p.signif', method = 'wilcox.test') +
  ylab('IC50') + xlab('') +
  theme(legend.position = 'top')
dev.off()
FIG6=mg_merge_plot(tide.p,ic50,heights = c(1,1.6),labels = c('','C'),ncol=1,nrow=2)
ggsave('06_TIDE/Fig6.pdf',FIG6,height =15,width =13)



####07.############
dir.create('07_nomogram')
######7.1 ####
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)
table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'

table(tcga_cox_datas$N.stage)
tcga_cox_datas$N.stage[tcga_cox_datas$N.stage=='N1'|tcga_cox_datas$N.stage=='N2'|tcga_cox_datas$N.stage=='N3']<-'N1+N2+N3'

table(tcga_cox_datas$M.stage)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'

mg_compare_uni_muti_cox_use=function(dat,event,os){
  # dat=crbind2DataFrame(dat)
  sig.clini=colnames(dat)
  dat$time=os
  dat$status=event
  dat=dat[which(!is.na(os)&!is.na(event)),]
  all.cox=rbind()
  rnames=c()
  for(s in sig.clini){
    fmla <- as.formula(paste0("Surv(time, status) ~",s))
    cox <- coxph(fmla, data = dat)
    #summary(cox)[[7]]
    #print(summary(cox))
    re=cbind(summary(cox)[[7]][,5],summary(cox)[[7]][,2],summary(cox)[[8]][,3],summary(cox)[[8]][,4])
    if(nrow(re)==1){
      rnames=c(rnames,s)
    }else{
      rnames=c(rnames,row.names(summary(cox)[[7]]))
    }
    all.cox=rbind(all.cox,re)
  }
  row.names(all.cox)=rnames
  colnames(all.cox)=c('p.value','HR','Low 95%CI','High 95%CI')
  
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(sig.clini,collapse = '+')))
  cox <- coxph(fmla, data = dat)
  muti.re=cbind(summary(cox)[[7]][,5],summary(cox)[[7]][,2],summary(cox)[[8]][,3],summary(cox)[[8]][,4])
  row.names(muti.re)=row.names(summary(cox)[[7]])
  colnames(muti.re)=c('p.value','HR','Low 95%CI','High 95%CI')
  return(list(muti=crbind2DataFrame(muti.re),uni=crbind2DataFrame(all.cox)))
}
getForestplotData=function(res_sig){
  res_sig<-signif(res_sig,digits=2)
  res_sig$CI_for_HR=paste0(" (",res_sig$`Low 95%CI`, "-", res_sig$`High 95%CI`, ")")
  colnames(res_sig)=c("p.value","HR","CI_lower","CI_upper", "(95%_CI_for_HR)")
  
  forest_table<-data.frame(Features=rownames(res_sig),HR=res_sig$HR,`(95%CI)`=res_sig$`(95%_CI_for_HR)`, pvalue=res_sig$p.value,check.names = F,stringsAsFactors = F)
  forest_table$sig<-mg_format_p_values(forest_table$pvalue)
  
  forest_table2<-data.frame(Features="Features",HR="HR",`(95%CI)`="(95%CI)", pvalue="p-value",sig='Significant',check.names = F,stringsAsFactors=F)
  tabletext<-rbind(forest_table2,forest_table)
  
  forest_stastic<-data.frame(mean=as.numeric(as.character(res_sig$HR)),lower=as.numeric(as.character(res_sig$CI_lower)),upper=as.numeric(as.character(res_sig$CI_upper)))
  forest_stastic1<-data.frame(mean=NA,lower=NA,upper=NA)
  cochrane_from_rNAM<-rbind(forest_stastic1,forest_stastic)
  return(list(tabletext,cochrane_from_rNAM))
}
colnames(tcga_cox_datas)
cm.cox=mg_compare_uni_muti_cox_use(tcga_cox_datas[,c(18,4,6:9,19)],
                                   os=tcga_cox_datas$OS.time,
                                   event = tcga_cox_datas$OS)

cm.cox.uni=signif(cm.cox$uni,digits=3)
cm.cox.uni$`Hazard Ratio(95%CI)`=paste0(cm.cox.uni$HR,"(",cm.cox.uni$`Low 95%CI`, "-", cm.cox.uni$`High 95%CI`, ")")
cm.cox.uni
table(cm.cox.uni$p.value<0.05)
cm.cox.uni[which(cm.cox.uni$p.value<0.05),]

# ####### 
cm.cox2=mg_compare_uni_muti_cox_use(tcga_cox_datas[,c('T.stage','N.stage','M.stage','Stage','Riskscore')],
                                    os=tcga_cox_datas$OS.time,
                                    event = tcga_cox_datas$OS)
cm.cox.muti=signif(cm.cox2$muti,digits=3)
cm.cox.muti
cm.cox.muti$`Hazard Ratio(95%CI)`=paste0(cm.cox.muti$HR,"(",cm.cox.muti$`Low 95%CI`, "-", cm.cox.muti$`High 95%CI`, ")")

writeMatrix(cm.cox.uni,outpath = '07_nomogram/tcga.cli.single.cox.txt')
writeMatrix(cm.cox.muti,outpath = '07_nomogram/tcga.cli.multi.cox.txt')
library(forestplot)
cm.cox.uni=getForestplotData(cm.cox$uni)
cm.cox.muti=getForestplotData(cm.cox2$muti)

pdf('07_nomogram/unicox.pdf',height = 5,width = 8,onefile = F)
tabletext=cm.cox.uni[[1]]
cochrane_from_rNAM=cm.cox.uni[[2]]
forestplot(tabletext, 
           mean=cochrane_from_rNAM$mean,
           lower=cochrane_from_rNAM$lower,
           upper=cochrane_from_rNAM$upper,
           graph.pos=2,#is.summary = c(T,F,F,F,F,F,F),
           hrzl_lines=list('2'=gpar(lty=1,col="black"),
                           '9'=gpar(lty=1,col="black")),
           zero = 1,xlog=F,
           fn.ci_norm = fpDrawDiamondCI,
           col=fpColors(line = "#669933", box="#669933",zero = '#669933',summary = 'black'),
           lty.ci = 3,lwd.ci =1,
           ci.vertices.height = 0.05,
           # lineheight = "auto", 
           xlab="Hazard ratio" )

dev.off()

pdf('07_nomogram/muticox.pdf',height = 5,width = 8,onefile = F)
tabletext=cm.cox.muti[[1]]
cochrane_from_rNAM=cm.cox.muti[[2]]
forestplot(tabletext, 
           mean=cochrane_from_rNAM$mean,
           lower=cochrane_from_rNAM$lower,
           upper=cochrane_from_rNAM$upper,
           graph.pos=2,is.summary = c(T,F,F,F,F,F),
           hrzl_lines=list('2'=gpar(lty=1,col="black"),
                           '7'=gpar(lty=1,col="black")),
           zero = 1,xlog=F,
           fn.ci_norm = fpDrawDiamondCI,
           col=fpColors(line = "#669933", box="#669933",zero = '#669933',summary = 'black'),
           lty.ci = 3,lwd.ci =1,
           ci.vertices.height = 0.05,
           # lineheight = "auto", 
           xlab="Hazard ratio" )
dev.off()

#######7.2 #######
mg_nomogram=function(clinical_riskscore,os,status,title='Nomogram',
                     quick=T,mks = c(1,3,5)){
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
  }
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
    }
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],
         lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}
pdf('07_nomogram/Nomogram.pdf', width = 12, height = 5)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                T.stage=tcga_cox_datas$T.stage,
                                N.stage=tcga_cox_datas$N.stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5))
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))

######C-INDEX#######
sum1=summary(coxph(formula=Surv(OS.time, OS)~T.stage+N.stage+Riskscore, data=tcga_cox_datas))
sum1$concordance[1]
c_index=as.numeric(sum1$concordance[1])

clinical.feature=c('Age','Gender','T.stage','N.stage','M.stage','Stage','Riskscore')
for (i in 1:length(clinical.feature)) {
  # i=2
  as.formula(paste0("Surv(OS.time, OS) ~",clinical.feature[i]))
  sum_res=summary(coxph(as.formula(paste0("Surv(OS.time, OS) ~",clinical.feature[i])), data=tcga_cox_datas))
  index=as.numeric(sum_res$concordance[1])
  c_index=append(c_index,index)
}
names(c_index)=c('Nomogram',clinical.feature)
c_index[order(c_index)]
pdf('07_nomogram/c_index_barplot.pdf',height = 5,width = 10,onefile = F)
barplot(c_index[order(c_index)],xlab="Clinical Features", ylab="C-index",ylim = c(0,0.9),main = 'Prediction effect',col = '#669933')
dev.off()

write.csv(data.frame(feature=names(c_index[order(c_index)]),c_index=as.numeric(c_index[order(c_index)])),
          'results/08.nomogram/c_index.csv',row.names = F)

#######7.3 ###########
####ROC
tcga_model_data_all=data.frame(type=ifelse(tcga_type$type=='Normal',0,1),t(tcga.exp.all[names(lan),tcga_type$Samples]))
tcga_model_data_all$type=factor(tcga_model_data_all$type)
head(tcga_model_data_all)

LR_model <- glm(type ~ ., data = tcga_model_data_all,  family=binomial(link = "logit"))
summary(LR_model)
pre_rate<-predict(LR_model )
modelroc <- roc(tcga_model_data_all$type,pre_rate)

pdf('07_nomogram/diag.roc.pdf',height = 4,width = 4)
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), grid.col=c("green", "red"), 
     max.auc.polygon=T, auc.polygon.col="#669933", print.thres=TRUE)
dev.off()

save.image('LUAD_NAM.RData')
