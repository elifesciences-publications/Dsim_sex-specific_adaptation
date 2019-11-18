#full analytical pipeline florida fullset data final version (publication)
#SKHsu_20191117

rm(list=ls())
gc()
library(limma)
library(edgeR)
library(pheatmap)
library(biomaRt)
library(topGO)
library(RcisTarget)
library(VennDiagram)
setwd("/Volumes/Temp1/shengkai/remap_run/florida_result_fullset_v7/publication/")
####import basic function####
cont_table=function(query,background,classifyer){
  p1=length(Reduce(intersect,list(query,background,classifyer)))
  q1=length(intersect(query,background))-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}
ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}
pseudo_chr=function(x){
  for (i in 1:length(unique(x$CHR))){
    x$pseudo_CHR[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=i
  }
  return(x)
}
pseudo_pos=function(x){
  for (i in 1:length(unique(x$CHR))){
    if (i==1) x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=x$BP[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]
    else x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]=x$BP[x$CHR%in%c("X","2L","2R","3L","3R","4")[i]]+max(x$pseudo_POS[x$CHR%in%c("X","2L","2R","3L","3R","4")[i-1]])
  }
  return(x)
}
p_resolution=function(x){
  x$P=10^(-as.numeric(strsplit2(x$"-logp",split = "=")[,2]))
  return(x)
}
SNPID_gen=function(x){
  x$SNPID=paste(x$CHR,x$BP,sep = "_")
  return(x)
}
snp_identifier=function(query,snpset,threshold1,threshold2){
  subset1=snpset[snpset$CHR%in%query$Chr&snpset$CHR=="X"&snpset$P<threshold1,]
  subset2=snpset[snpset$CHR%in%query$Chr&snpset$CHR!="X"&snpset$P<threshold2,]
  subset=rbind(subset1,subset2)
  count=sum(query$start<subset$BP&query$end>subset$BP,na.rm = T)
  SNPID=as.character(na.omit(subset$SNPID[query$start<subset$BP&query$end>subset$BP]))
  
  query$SNP_counts=count
  if(length(SNPID)>0) query$SNP_ID=paste(SNPID,collapse = ",")
  else query$SNP_ID=NA
  return(query)
}

#################################
####Part1: Data import and QC####
#################################
####Data input and flitering####
#RNASeq_CGE1
count_dat1=read.csv("./RNASeq_CGE1_count_table.csv",stringsAsFactors = F,header = T,row.names = 1)

#RNASeq_CGE2
count_dat2=read.csv("./RNASeq_CGE2_count_table.csv",stringsAsFactors = F,row.names = 1)
raw_counts_carc=count_dat2[,grep(".....c",colnames(count_dat2))]
raw_counts_gonad=count_dat2[,grep(".....g",colnames(count_dat2))]
raw_counts_body=count_dat2[,grep(".....b",colnames(count_dat2))]

#cpm filtration
raw_counts_dat1=count_dat1[apply(cpm(count_dat1),1,function(x) !sum(x<0.1)>=1),]
raw_counts_carc_filtered=raw_counts_carc[rownames(raw_counts_dat1),]
raw_counts_gonad_filtered=raw_counts_gonad[rownames(raw_counts_dat1),]
raw_counts_body_filtered=raw_counts_body[rownames(raw_counts_dat1),]

#DGEset and grouping
evo=substr(colnames(count_dat1),1,1)
sex=substr(colnames(count_dat1),11,11)
group=paste0(evo,sex)
y=DGEList(counts=raw_counts_dat1,group = group)

group_carc=paste0(substr(colnames(raw_counts_carc_filtered),1,1),substr(colnames(raw_counts_carc_filtered),5,5))
group_gonad=paste0(substr(colnames(raw_counts_gonad_filtered),1,1),substr(colnames(raw_counts_gonad_filtered),5,5))
group_body=paste0(substr(colnames(raw_counts_body_filtered),1,1),substr(colnames(raw_counts_body_filtered),5,5))

y_carc=DGEList(counts=raw_counts_carc_filtered,group = group_carc)
y_gonad=DGEList(counts=raw_counts_gonad_filtered,group = group_gonad)
y_body=DGEList(counts=raw_counts_body_filtered,group = group_body)

y_carc_f=y_carc[,substr(colnames(raw_counts_carc_filtered),5,5)%in%"f"]
y_carc_m=y_carc[,substr(colnames(raw_counts_carc_filtered),5,5)%in%"m"]
y_gonad_f=y_gonad[,substr(colnames(raw_counts_gonad_filtered),5,5)%in%"f"]
y_gonad_m=y_gonad[,substr(colnames(raw_counts_gonad_filtered),5,5)%in%"m"]
y_body_f=y_body[,substr(colnames(raw_counts_body_filtered),5,5)%in%"f"]
y_body_m=y_body[,substr(colnames(raw_counts_body_filtered),5,5)%in%"m"]


################################################
####Part2: Correction for allometric factors####
################################################
####estimation of the confounding effect from allometric evolution####
mean_body_f=t(apply(cpm(y_body_f),1,function(x) tapply(x,substr(names(x),1,1),mean)))
mean_carc_f=t(apply(cpm(y_carc_f),1,function(x) tapply(x,substr(names(x),1,1),mean)))
mean_gonad_f=t(apply(cpm(y_gonad_f),1,function(x) tapply(x,substr(names(x),1,1),mean)))
mean_body_m=t(apply(cpm(y_body_m),1,function(x) tapply(x,substr(names(x),1,1),mean)))
mean_carc_m=t(apply(cpm(y_carc_m),1,function(x) tapply(x,substr(names(x),1,1),mean)))
mean_gonad_m=t(apply(cpm(y_gonad_m),1,function(x) tapply(x,substr(names(x),1,1),mean)))

a.hat.anc=c()
for(i in 1:dim(y)[1]){
  a.hat.anc=c(a.hat.anc,(mean_body_f[i,1]-mean_carc_f[i,1])/(mean_gonad_f[i,1]-mean_carc_f[i,1]))
}
a.hat.anc.filter=a.hat.anc[!(a.hat.anc>1|a.hat.anc<0)]

a.hat=c()
for(i in 1:dim(y)[1]){
  a.hat=c(a.hat,(mean_body_f[i,2]-mean_carc_f[i,2])/(mean_gonad_f[i,2]-mean_carc_f[i,2]))
}
a.hat.filter=a.hat[!(a.hat>1|a.hat<0)]

b.hat.anc=c()
for(i in 1:dim(y)[1]){
  b.hat.anc=c(b.hat.anc,(mean_body_m[i,1]-mean_carc_m[i,1])/(mean_gonad_m[i,1]-mean_carc_m[i,1]))
}
b.hat.anc.filter=b.hat.anc[!(b.hat.anc>1|b.hat.anc<0)]

b.hat=c()
for(i in 1:dim(y)[1]){
  b.hat=c(b.hat,(mean_body_m[i,2]-mean_carc_m[i,2])/(mean_gonad_m[i,2]-mean_carc_m[i,2]))
}
b.hat.filter=b.hat[!(b.hat>1|b.hat<0)]


ks.test(a.hat.filter,a.hat.anc.filter)
ks.test(b.hat.filter,b.hat.anc.filter)

allo.est=sapply(list(a.hat.anc.filter,a.hat.filter,b.hat.anc.filter,b.hat.filter),median)

png("./FigureS6.png",width = 12,height = 12,units = "cm",res = 600, pointsize = 10)
plot_dat=data.frame(a=c(sort(a.hat.anc.filter),sort(a.hat.filter),sort(b.hat.anc.filter),sort(b.hat.filter)),
                    pop=as.factor(c(rep("female anc.",length(a.hat.anc.filter)),rep("female evo.",length(a.hat.filter)),
                                    rep("male anc.",length(b.hat.anc.filter)),rep("male evo.",length(b.hat.filter)))))
qplot(pop,a,data = plot_dat,geom = "violin",fill=pop,alpha=I(0.5),xlab = "",ylab="allometric estimates")+
  theme_classic()+
  geom_boxplot(width=0.1,fill="white")+
  scale_fill_manual(values=c("orange","salmon","cyan","royalblue"))
dev.off()


####recontruct pseudo WB samples####
index_f_allo=!(a.hat.anc>1|a.hat.anc<0)&!(a.hat>1|a.hat<0)
index_m_allo=!(b.hat.anc>1|b.hat.anc<0)&!(b.hat>1|b.hat<0)
y_f_evo_hat=data.frame(matrix(rep(mean_body_f[,1],6),9457,6))
y_f_evo=data.frame(matrix(rep(mean_body_f[,2],6),9457,6))
y_m_evo_hat=data.frame(matrix(rep(mean_body_m[,1],6),9457,6))
y_m_evo=data.frame(matrix(rep(mean_body_m[,2],6),9457,6))
y_f_evo_hat[index_f_allo,]=(1-a.hat[index_f_allo])*cpm(y_carc[index_f_allo,y_carc$samples$group%in%"Bf"])+
  a.hat[index_f_allo]*cpm(y_gonad[index_f_allo,y_gonad$samples$group%in%"Bf"])
y_f_evo[index_f_allo,]=(1-a.hat[index_f_allo])*cpm(y_carc[index_f_allo,y_carc$samples$group%in%"Hf"])+
  a.hat[index_f_allo]*cpm(y_gonad[index_f_allo,y_gonad$samples$group%in%"Hf"])

y_m_evo_hat[index_m_allo,]=(1-b.hat[index_m_allo])*cpm(y_carc[index_m_allo,y_carc$samples$group%in%"Bm"])+
  b.hat[index_m_allo]*cpm(y_gonad[index_m_allo,y_gonad$samples$group%in%"Bm"])
y_m_evo[index_m_allo,]=(1-b.hat[index_m_allo])*cpm(y_carc[index_m_allo,y_carc$samples$group%in%"Hm"])+
  b.hat[index_m_allo]*cpm(y_gonad[index_m_allo,y_gonad$samples$group%in%"Hm"])

####calculating the correcting factor and apply it to monster CGE####
corr_factor_f=apply(y_f_evo_hat,1,mean)/mean_body_f[,1]
corr_factor_m=apply(y_m_evo_hat,1,mean)/mean_body_m[,1]
corrected_anc_f=y$counts[,y$samples$group%in%"Bf"]*corr_factor_f
corrected_anc_m=y$counts[,y$samples$group%in%"Bm"]*corr_factor_m
##############################################
####Part3: Linear modeling and DE analysis####
##############################################
####specify experimental design and known factors####
group_corrected = c(rep("Bf",5),rep("Bm",5),group[-1:-10])
sex_corrected=substr(group_corrected,2,2)
evo_corrected=substr(group_corrected,1,1)
count_corrected=cbind(corrected_anc_f,corrected_anc_m,y$counts[,-1:-10])
y_corrected=DGEList(counts=count_corrected,
                    group = group_corrected)
y_corrected=DGEList(counts=cbind(round(cbind(corrected_anc_f,corrected_anc_m),0),y$counts[,-1:-10]),
                    group = group_corrected)
y_corrected=calcNormFactors(y_corrected)
colnames(y_corrected)[1:10]=c("B_pesudo1_f","B_pesudo2_f","B_pesudo3_f","B_pesudo4_f","B_pesudo5_f",
                              "B_pesudo1_m","B_pesudo2_m","B_pesudo3_m","B_pesudo4_m","B_pesudo5_m")
colnames(y_corrected)=gsub("2818","",colnames(y_corrected))
#write.csv(y_corrected$counts,"./corrected_readcount.csv",quote = F)
#count_corrected=read.csv("./corrected_readcount.csv",header = T,row.names = 1)
####estimate the variance explained by known variables (PCA)####
pca=prcomp(t(log(cpm(y_corrected))))
plot(pca$sdev/sum(pca$sdev))
plot(pca$x[,c(1,2)],col=ifelse(substr(group_corrected,1,1)=="H","salmon","seagreen"),
     pch=ifelse(substr(group_corrected,2,2)=="m",21,22))
plot(pca$x[,c(1,3)],col=ifelse(substr(group_corrected,1,1)=="H","salmon","seagreen"),
     pch=ifelse(substr(group_corrected,2,2)=="m",21,22))
pca_m=prcomp(t(log(cpm(y_corrected[,sex_corrected%in%"m"]))))
pca_f=prcomp(t(log(cpm(y_corrected[,sex_corrected%in%"f"]))))

layout(matrix(c(1:3,3),2,2,byrow = T),width=c(4,4),height=c(4,1),respect=T)
layout.show(3)
plot(pca_m$x,col=ifelse(evo_corrected[sex_corrected%in%"m"]=="H","salmon","seagreen"),pch=19,asp=1,
     main="male",xlab="PC1 (16.35%)",ylab="PC2 (6.99%)")
plot(pca_f$x,col=ifelse(evo_corrected[sex_corrected%in%"f"]=="H","salmon","seagreen"),pch=19,asp=1,
     main="female",xlab="PC1 (23.28%)",ylab="PC2 (10.22%)")
par(mar=c(0,4.1,0,2.1))
plot(0,0,type = 'n',axes = F,xlab = "",ylab="")
legend("top",c("anc. corrected","evo. ori."),col = c("seagreen","salmon"),pch=19,cex=1,horiz = T)


####modeling and DE analysis####
ModelDesign=model.matrix(~0+group_corrected)
DGE=estimateDisp(y_corrected,design = ModelDesign,robust = T)
GLM=glmFit(DGE,design = ModelDesign)
my.contrasts=makeContrasts(sex_b=group_correctedBf-group_correctedBm,
                           sex_h=group_correctedHf-group_correctedHm,
                           sex=(group_correctedBf+group_correctedHf)/2-(group_correctedBm+group_correctedHm)/2,
                           hotevolF=group_correctedHf-group_correctedBf,
                           hotevolM=group_correctedHm-group_correctedBm,
                           interH=(group_correctedHf-group_correctedBf)-(group_correctedHm-group_correctedBm),
                           levels=ModelDesign)
LRT_hotevoF=glmLRT(GLM,contrast = my.contrasts[,"hotevolF"])
LRT_hotevoM=glmLRT(GLM,contrast = my.contrasts[,"hotevolM"])
LRT_sexbiasedH=glmLRT(GLM,contrast = my.contrasts[,"sex_h"])
LRT_sexbiasedB=glmLRT(GLM,contrast = my.contrasts[,"sex_b"])
LRT_sexbiased=glmLRT(GLM,contrast = my.contrasts[,c("sex")])
LRT_interH=glmLRT(GLM,contrast = my.contrasts[,"interH"])

res_hotevoF=LRT_hotevoF$table
res_hotevoF$padj=p.adjust(res_hotevoF$PValue,method = "BH")
res_hotevoM=LRT_hotevoM$table
res_hotevoM$padj=p.adjust(res_hotevoM$PValue,method = "BH")
res_sexbiasedH=LRT_sexbiasedH$table
res_sexbiasedH$padj=p.adjust(res_sexbiasedH$PValue,method = "BH")
res_sexbiasedB=LRT_sexbiasedB$table
res_sexbiasedB$padj=p.adjust(res_sexbiasedB$PValue,method = "BH")
res_sexbiased=LRT_sexbiased$table
res_sexbiased$padj=p.adjust(res_sexbiased$PValue,method = "BH")

res_interH=LRT_interH$table
res_interH$padj=p.adjust(res_interH$PValue,method = "BH")
corrected_res=list(res_hotevoF,res_hotevoM,res_interH,res_sexbiasedB,res_sexbiasedH,res_sexbiased)

write.csv(res_hotevoF,"./LRT_res_hotevoF.csv",quote = F)
write.csv(res_hotevoM,"./LRT_res_hotevoM.csv",quote = F)
write.csv(res_sexbiasedB,"./LRT_res_sexbiasedB.csv",quote = F)
write.csv(res_sexbiasedH,"./LRT_res_sexbiasedH.csv",quote = F)
write.csv(res_sexbiased,"./LRT_res_sexbiased.csv",quote = F)
write.csv(res_interH,"./LRT_res_interH.csv",quote = F)

####output####
write.table(rownames(corrected_res[[1]])[corrected_res[[2]]$padj<0.05&corrected_res[[2]]$logFC>0],
            "./goi/m_up.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[2]]$padj<0.05&corrected_res[[2]]$logFC<0],
            "./goi/m_dn.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[1]]$padj<0.05&corrected_res[[1]]$logFC>0],
            "./goi/f_up.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[1]]$padj<0.05&corrected_res[[1]]$logFC<0],
            "./goi/f_dn.txt",
            quote = F,col.names = F,row.names = F)

write.table(rownames(corrected_res[[1]])[corrected_res[[2]]$padj<0.05&corrected_res[[2]]$logFC>0&!corrected_res[[1]]$padj<0.05],
            "./goi/m_s_up.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[2]]$padj<0.05&corrected_res[[2]]$logFC<0&!corrected_res[[1]]$padj<0.05],
            "./goi/m_s_dn.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[1]]$padj<0.05&corrected_res[[1]]$logFC>0&!corrected_res[[2]]$padj<0.05],
            "./goi/f_s_up.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[1]]$padj<0.05&corrected_res[[1]]$logFC<0&!corrected_res[[2]]$padj<0.05],
            "./goi/f_s_dn.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[1]]$padj<0.05&corrected_res[[2]]$padj<0.05&corrected_res[[2]]$logFC>0&corrected_res[[1]]$logFC>0],
            "./goi/conc_up.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[1]]$padj<0.05&corrected_res[[2]]$padj<0.05&corrected_res[[2]]$logFC<0&corrected_res[[1]]$logFC<0],
            "./goi/conc_dn.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[1]]$padj<0.05&corrected_res[[2]]$padj<0.05&corrected_res[[2]]$logFC>0&corrected_res[[1]]$logFC<0],
            "./goi/anta_mf.txt",
            quote = F,col.names = F,row.names = F)
write.table(rownames(corrected_res[[1]])[corrected_res[[1]]$padj<0.05&corrected_res[[2]]$padj<0.05&corrected_res[[2]]$logFC<0&corrected_res[[1]]$logFC>0],
            "./goi/anta_fm.txt",
            quote = F,col.names = F,row.names = F)
background=rownames(y_corrected)
write.table(background,"./goi/background.txt",quote = F,row.names = F,col.names = F)


####scatter plots####
mycol=rep(NA,dim(corrected_res[[1]])[1])
mycol[corrected_res[[2]]$padj<0.05&corrected_res[[1]]$padj>0.05]="royalblue"
mycol[corrected_res[[1]]$padj<0.05&corrected_res[[2]]$padj>0.05]="salmon"
mycol[corrected_res[[1]]$padj<0.05&corrected_res[[2]]$padj<0.05]="seagreen"
mycol[corrected_res[[2]]$logFC*corrected_res[[1]]$logFC<0&corrected_res[[1]]$padj<0.05&corrected_res[[2]]$padj<0.05]="gold"

png("./Figure1a.png",width = 12,height = 12,units = "cm",res = 600,pointsize = 10)
par(mar=c(5,5,2,2))
plot(corrected_res[[1]]$logFC,corrected_res[[2]]$logFC,xlab="Evolutionary response in female",ylab="Evolutionary response in male",
     xlim = c(-4,4),ylim = c(-4,4),pch=19,col="grey",asp=1,cex=0.75)
abline(h=0,v=0,lty=3,col="grey70")
abline(a=0,b=1,lty=3,col="grey70")
abline(a=0,b=-1,lty=3,col="grey70")
points(corrected_res[[1]]$logFC,corrected_res[[2]]$logFC,pch=19,col=mycol,asp=1,cex=0.75)
legend("topleft",legend = c("male-specific (1295)","female-specific (3080)","concordant (760)","antagonistic (311)","n.s. (4011)"),
       col =c("royalblue","salmon","seagreen","gold","grey"),pch = 19 ,bty = "n",cex = 1)
dev.off()

####sex-biased genes and the evolution of sex-bias####
f_bias_base=background[corrected_res[[4]]$logFC>0&corrected_res[[4]]$padj<0.05]
m_bias_base=background[corrected_res[[4]]$logFC< 0&corrected_res[[4]]$padj<0.05]
f_bias_hot=background[corrected_res[[5]]$logFC>0&corrected_res[[5]]$padj<0.05]
m_bias_hot=background[corrected_res[[5]]$logFC< 0&corrected_res[[5]]$padj<0.05]
f_bias=background[corrected_res[[6]]$logFC>0&corrected_res[[6]]$padj<0.05]
m_bias=background[corrected_res[[6]]$logFC< 0&corrected_res[[6]]$padj<0.05]
sex_bias_ID=list(f_bias_base,m_bias_base,f_bias_hot,m_bias_hot,f_bias,m_bias)
names(sex_bias_ID)=c("f_biased_base","m_biased_base","f_biased_evo","m_biased_evo","f_biased","m_biased")

png("./FigureS5.png",width = 12,height=8,units = "cm",res = 600,pointsize = 8)
par(mfrow=c(1,2),mar=c(5,5,4,2))
boxplot(abs(res_sexbiasedB[setdiff(Reduce(union,sex_bias_ID[3:4]),Reduce(union,sex_bias_ID[1:2])),]$logFC),
        abs(res_sexbiasedH[setdiff(Reduce(union,sex_bias_ID[3:4]),Reduce(union,sex_bias_ID[1:2])),]$logFC),
        ylab=expression(paste("Magnitude of dimorphism ",abs(log[2](FC)))),xlab="Population",names=c("Anc.","Evo."),main="Gain")
boxplot(abs(res_sexbiasedB[setdiff(Reduce(union,sex_bias_ID[1:2]),Reduce(union,sex_bias_ID[3:4])),]$logFC),
        abs(res_sexbiasedH[setdiff(Reduce(union,sex_bias_ID[1:2]),Reduce(union,sex_bias_ID[3:4])),]$logFC),
        ylab=expression(paste("Magnitude of dimorphism ",abs(log[2](FC)))),xlab="Population",names=c("Anc.","Evo."),main="Loss")
dev.off()
####consistency between replicates####
expr_sig_m=log2(cpm(y_corrected)[corrected_res[[2]]$padj<0.05,sex_corrected%in%"m"])
expr_sig_f=log2(cpm(y_corrected)[corrected_res[[1]]$padj<0.05,sex_corrected%in%"f"])
colnames(expr_sig_m)=c("B01","B02","B03","B04","B05",paste0("H",rep(sprintf("%02d",1:10),each=3)))
colnames(expr_sig_f)=c("B01","B02","B03","B04","B05",paste0("H",rep(sprintf("%02d",1:10),each=2)))

pheatmap(expr_sig_m,scale = "row",show_rownames = F,filename = "./FigureS1a.png",res=600)
pheatmap(expr_sig_f,scale = "row",show_rownames = F,filename = "./FigureS1b.png",res=600)
###############################################################################
####Part4: Gene set enrichment analysis and related phenotyping experiments####
###############################################################################
####data input####
flyatlas=read.table("./fly_atlas_enrichment.table",header = T,stringsAsFactors = F)
query_path="./goi/"
query_files=list.files(path=query_path,pattern=".txt")[-3]
query_ID=lapply(query_files,function(x) read.table(paste0(query_path,x),header=F,stringsAsFactors = F)[,1])
names(query_ID)=strsplit2(query_files,".txt")[,1]
sapply(query_ID,length)

venn.diagram(query_ID[c(5,8,9,12)],
             "./Figure1b.png",
             imagetype = "png",height = 12,width = 12,units = "cm",resolution = 600,
             fill=c("violet","salmon","skyblue","Royalblue"),main = "DE genes",cex=2,cat.cex=1.5,main.cex = 2,
             category.names = paste0(rep(c("F.","M."),each=2),c("down","up")))

gtf_transform=read.table("./gtf_resolved_table.txt",stringsAsFactors = F,header = F)
gene_chr=tapply(gtf_transform[,1],INDEX = gtf_transform$V5,unique)
gene_start=tapply(gtf_transform[,2],INDEX = gtf_transform$V5,min)
gene_end=tapply(gtf_transform[,3],INDEX = gtf_transform$V5,max)
gene_strand=tapply(gtf_transform[,4],INDEX = gtf_transform$V5,unique)
gene_summary=data.frame(gene_chr,gene_start,gene_end,gene_strand,sort(unique(gtf_transform$V5)))
colnames(gene_summary)=c("Chr","start","end","strand","FB_ID")
gene_summary_filtered=gene_summary[background,]

####candidate gene ID conversion####
ensembl=useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
background_conv=background
latest_FB_ID=read.delim("./fbgn_annotation_ID_fb_2018_01_R.txt",header = T,stringsAsFactors = F,sep="\t")
for(i in 1:dim(latest_FB_ID)[1]){
  background_conv[which(background%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
}
query_ID_new=lapply(query_ID,function(x) {
  for(i in 1:dim(latest_FB_ID)[1]){
    x[which(x%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
  }
  return(x)
})
sex_bias_ID_new=lapply(sex_bias_ID,function(x) {
  for(i in 1:dim(latest_FB_ID)[1]){
    x[which(x%in%strsplit2(latest_FB_ID[i,4],","))]=latest_FB_ID[i,3]
  }
  return(x)
})
conv_query=lapply(query_ID_new,function(x) ID_converter(ID = x,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1])
conv_background=ID_converter(ID=background_conv,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1]
sex_bias_ID_conv=lapply(sex_bias_ID_new,function(x) ID_converter(ID = x,db = ensembl,attributes = "external_gene_name",filters = "flybase_gene_id")[,1])

gene_summary_filtered$FB_ID_new=c(background_conv)
gene_summary_filtered=gene_summary_filtered[gene_summary_filtered[,1]%in%c("2L","2R","3L","3R","4","X"),]
rownames(gene_summary_filtered)=gene_summary_filtered$FB_ID

####GO enrichment analysis####
GO_res_table=list()
for (i in names(query_ID_new)){
  tmp=factor(as.integer(background_conv%in%query_ID_new[[i]]))
  names(tmp)=background_conv
  tgd=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  GO_res_table[[i]]=tmp_res
}
GO_res_table=lapply(GO_res_table,function(x) {
  x$Fisher.weight01[x$Fisher.weight01=="< 1e-30"]=1e-30
  return(x)})

sig_GO=sapply(GO_res_table,function(x) x$GO.ID[x$Fisher.weight01<0.05])
venn.diagram(sig_GO[c(5,8,9,12)],
             "./Figure1c.png",
             imagetype = "png",height = 12,width = 12,units = "cm",resolution = 600,
             fill=c("violet","salmon","skyblue","Royalblue"),main = "GO: BP",cex=2,cat.cex=1.5,main.cex = 2,
             category.names = paste0(rep(c("F.","M."),each=2),c("down","up")))

####tissue_enrichment####
p.val=c()
odds=c()
exp.num=c()
tissue_enriched_ID=list()
for(i in seq(4,36,2)){
  tissue_specific_ID=flyatlas$FB[flyatlas[,i]>2]
  p=c()
  o=c()
  e=c()
  g=list()
  for (j in names(query_ID_new)){
    p=c(p,fisher.test(cont_table(query_ID_new[[j]],background_conv,tissue_specific_ID),alternative = "greater")$p.value)
    o=c(o,fisher.test(cont_table(query_ID_new[[j]],background_conv,tissue_specific_ID),alternative = "greater")$estimate)
    g[[j]]=intersect(query_ID_new[[j]],tissue_specific_ID)
    e=c(e,length(query_ID_new[[j]])*sum(tissue_specific_ID%in%background_conv)/length(background_conv))
  }  
  p.val=rbind(p.val,p)
  odds=rbind(odds,o)
  tissue_enriched_ID[[(i-2)/2]]=g
  exp.num=rbind(exp.num,e)
}
odds=odds[-c(7:9,13:14),]
exp.num=exp.num[-c(7:9,13:14),]
exp.num=round(exp.num,0)
obs.num= t(sapply(tissue_enriched_ID,function(x) sapply(x,length)))[-c(7:9,13:14),]
padj=matrix(p.adjust(p.val[-c(7:9,13:14),],method = "BH"),12,12)
names(tissue_enriched_ID)=strsplit2(colnames(flyatlas)[seq(4,36,2)],split = "_")[,1]
rownames(padj)=c("Br","Hd","Cr","Mg","Hg","Tb","Tg","Cs","Sg","Fb","Ey","Hr")
colnames(padj)=paste(names(query_ID),sapply(query_ID,length))
rownames(odds)=c("Br","Hd","Cr","Mg","Hg","Tb","Tg","Cs","Sg","Fb","Ey","Hr")
colnames(odds)=paste(names(query_ID),sapply(query_ID,length))
rownames(exp.num)=c("Br","Hd","Cr","Mg","Hg","Tb","Tg","Cs","Sg","Fb","Ey","Hr")
rownames(obs.num)=c("Br","Hd","Cr","Mg","Hg","Tb","Tg","Cs","Sg","Fb","Ey","Hr")

breaks=-log10(c(1,0.05,0.01,0.001,1e-1000))
color=c("lightyellow","gold","orange","firebrick")
pheatmap(-log10(padj),cluster_cols = F,cluster_rows = F,breaks = breaks,color = color,legend = F)
pheatmap(t(odds[,c(8,5,12,9)]),cluster_cols = F,cluster_rows = F,color = brewer.pal(5,"OrRd"),cex=1.5,display_numbers = T,
        filename = "./Figure1d.png",height = 8/2.54,width = 24/2.54,legend = F,
        labels_row = c("F.up","F.down","M.up","M.down"))

####expression of GOI####
#fatty acid metabolism
sapply(corrected_res,function(x) x["FBgn0012034",])#AcCoAS
sapply(corrected_res,function(x) x["FBgn0033246",])#ACC
sapply(corrected_res,function(x) x["FBgn0027571",])#FASN1
sapply(corrected_res,function(x) x["FBgn0042627",])#FASN2
sapply(corrected_res,function(x) x["FBgn0040001",])#FASN3
sapply(corrected_res,function(x) x["FBgn0040064",])#YIP2(ACAA2)
sapply(corrected_res,function(x) x["FBgn0025352",])#Thiolase
sapply(corrected_res,function(x) x["FBgn0028479",])#Mtpalpha (HADHA)
sapply(corrected_res,function(x) x["FBgn0033879",])#ECHS1
sapply(corrected_res,function(x) x["FBgn0033883",])#MECR
sapply(corrected_res,function(x) x["FBgn0263120",])#Acsl
sapply(corrected_res,function(x) x["FBgn0261862",])#CPT1
sapply(corrected_res,function(x) x["FBgn0035383",])#CPT2

FA_GOI=rev(c("FBgn0012034","FBgn0033246","FBgn0027571","FBgn0042627","FBgn0040001",
             "FBgn0040064","FBgn0025352","FBgn0028479","FBgn0033879","FBgn0033883","FBgn0263120",
             "FBgn0261862","FBgn0035383"))
FA_GOI_s=rev(c("AcCoAS","ACC","FASN1","FASN2","FASN3","yip2","Thiolase","Mtpalpha","Echs1","CG16935","Acsl","CPT1","CPT2"))
FC_FA_GOI=sapply(corrected_res[1:2],function(x) x[FA_GOI,]$logFC)
p_FA_GOI=sapply(corrected_res[1:2],function(x) x[FA_GOI,]$padj)
rownames(FC_FA_GOI)=FA_GOI_s
colnames(FC_FA_GOI)=c("females","males")
png("./Figure2a.png",width = 12,height = 12,units = "cm",res=600,pointsize = 10)
par(mar=c(5.1,5.1,2.1,2.1))
barplot(t(FC_FA_GOI),beside = T,horiz = T,las=1,font.axis=3,col=c("salmon","royalblue"),
        xlim=c(-1,1)*1.1,axes = F,names.arg = NA)
rect(rep(-1.1,2),c(0,24.5),c(1.1,1),c(6.5,40),col = "grey90",border = NA)
barplot(t(FC_FA_GOI),beside = T,horiz = T,las=1,font.axis=3,col=c("salmon","royalblue"),
        xlim=c(-1,1)*1.1,xlab=expression(log[2](FC)),add=T)
points(ifelse(as.vector(FC_FA_GOI)>0,as.vector(FC_FA_GOI)+0.025,as.vector(FC_FA_GOI)-0.025),
       c(seq(from=2,by = 3,length.out = length(FA_GOI))-0.5,seq(from=2,by = 3,length.out = length(FA_GOI))+0.5),
       pch=ifelse(as.vector(p_FA_GOI)<0.05,"*",""),cex=1.5)
text(c(0.75,0.75,0.75),c(4,15.5,32.25),labels = rev(c("Synthesis","Elongation","Degradation")),font = 4)
dev.off()

#dopaminergic neuronal activity
sapply(corrected_res,function(x) x["FBgn0014859",])#Hr38
sapply(corrected_res,function(x) x["FBgn0005626",])#ple
sapply(corrected_res,function(x) x["FBgn0000422",])#Ddc
sapply(corrected_res,function(x) x["FBgn0019643",])#Dat
sapply(corrected_res,function(x) x["FBgn0260964",])#Vmat
sapply(corrected_res,function(x) x["FBgn0034136",])#DAT
sapply(corrected_res,function(x) x["FBgn0011582",])#Dop1R1 filtered out (lowly expressed)
sapply(corrected_res,function(x) x["FBgn0015129",])#Dop1R2
sapply(corrected_res,function(x) x["FBgn0053517",])#Dop2R
sapply(corrected_res,function(x) x["FBgn0035538",])#DopEcR

dop_GOI=rev(c("FBgn0014859","FBgn0005626","FBgn0000422","FBgn0019643","FBgn0260964","FBgn0034136",
              "FBgn0015129","FBgn0053517","FBgn0035538"))
dop_GOI_s=rev(c("Hr38","Ple","Ddc","Dat","Vmat","DAT","Dop1R2","Dop2R","DopEcR"))
FC_dop_GOI=sapply(corrected_res[1:2],function(x) x[dop_GOI,]$logFC)
p_dop_GOI=sapply(corrected_res[1:2],function(x) x[dop_GOI,]$padj)
rownames(FC_dop_GOI)=dop_GOI_s
colnames(FC_dop_GOI)=c("females","males")

#octopaminergic
sapply(corrected_res,function(x) x["FBgn0259977",])#Tdc1
sapply(corrected_res,function(x) x["FBgn0050446",])#Tdc2
sapply(corrected_res,function(x) x["FBgn0010329",])#Tbh not annotated...
sapply(corrected_res,function(x) x["FBgn0260964",])#Vmat
sapply(corrected_res,function(x) x["FBgn0024944",])#Oamb
sapply(corrected_res,function(x) x["FBgn0038653",])#Octa2R
sapply(corrected_res,function(x) x["FBgn0038980",])#Octb1R
sapply(corrected_res,function(x) x["FBgn0038063",])#Octb2R
sapply(corrected_res,function(x) x["FBgn0250910",])#OCtb3R
sapply(corrected_res,function(x) x["FBgn0004514",])#Oct-TyrR

oct_GOI=rev(c("FBgn0259977","FBgn0050446","FBgn0260964","FBgn0024944","FBgn0038653","FBgn0038980",
              "FBgn0038063","FBgn0250910","FBgn0004514"))
oct_GOI_s=rev(c("Tdc1","Tdc2","Vmat","Oamb","Octa2R","Octb1R","Octb2R","Octb3R","Oct-TyrR"))
FC_oct_GOI=sapply(corrected_res[1:2],function(x) x[oct_GOI,]$logFC)
p_oct_GOI=sapply(corrected_res[1:2],function(x) x[oct_GOI,]$padj)
rownames(FC_oct_GOI)=oct_GOI_s
colnames(FC_oct_GOI)=c("females","males")

FC_dopoct_GOI=rbind(FC_dop_GOI,FC_oct_GOI)
p_dopoct_GOI=rbind(p_dop_GOI,p_oct_GOI)
png("./Figure2c.png",width = 12,height = 12,units = "cm",res=600,pointsize = 10)
par(mar=c(5.1,5.1,2.1,2.1))
barplot(t(FC_dopoct_GOI),beside = T,horiz = T,las=1,font.axis=3,col=c("salmon","royalblue"),
        xlim=c(-1,1),xlab=expression(log[2](FC)),main="",names.arg = NA)
rect(-1.1,0,1.1,27.5,col = "grey90",border = NA)
barplot(t(FC_dopoct_GOI),beside = T,horiz = T,las=1,font.axis=3,col=c("salmon","royalblue"),
        xlim=c(-1,1),xlab=expression(log[2](FC)),main="",add=T)
points(ifelse(as.vector(FC_dopoct_GOI)>0,as.vector(FC_dopoct_GOI)+0.025,as.vector(FC_dopoct_GOI)-0.025),
       c(seq(from=2,by = 3,length.out = dim(FC_dopoct_GOI)[1])-0.5,seq(from=2,by = 3,length.out = dim(FC_dopoct_GOI)[1])+0.5),
       pch=ifelse(as.vector(p_dopoct_GOI)<0.05,"*",""),cex=1.5)
text(c(0.75,0.75)-0.1,c(13,42.5),labels = c("Dopaminergic","Octopaminergic"),font = 4)
dev.off()



####phenotyping experiments####
#TAG content barghi et al., 2019
fat=read.delim("./bodyFatContent.txt",header = T,stringsAsFactors = F)
anova(lm(fat$LperFly~fat$Sex+fat$Pop))
hsd_res=HSD.test(fat$LperFly,paste(fat$Sex,fat$Pop),DFerror = 100,MSerror = 57)


png("./Figure2b.png",width = 12,height = 12,units = "cm",pointsize = 10,res = 600)
par(mar=c(4.1,4.1,2,2),las=1)
boxplot(fat$LperFly~fat$Pop+fat$Sex,col=c("pink","salmon","skyblue","royalblue"),names=c("Anc. F.","Evo. F.","Anc. M.","Evo. M."),
        ylab=expression(paste("TAG equivalents (",mu,"g) / fly")),ylim=c(20,95))
text(1:4,95,labels = c("a","b","c","c"))
dev.off()

#male reproductive activity
rel_chasing_duration=read.table("./rel_chasing_activity.txt",header = T)
rel_attcop_duration=read.table("./rel_attcop_duration.txt",header = T)

wilcox.test(rel_chasing_duration$rel_chasing_duration[rel_chasing_duration$evo%in%"B"],
            rel_chasing_duration$rel_chasing_duration[rel_chasing_duration$evo%in%"H"])
wilcox.test(rel_attcop_duration$rel_attcop_duration[rel_attcop_duration$evo%in%"B"],
            rel_attcop_duration$rel_attcop_duration[rel_attcop_duration$evo%in%"H"])

png("./Figure2e.png",width = 6,height = 12,units = "cm",res = 600,pointsize = 10)
par(mar=c(4.1,4.1,2,1),las=1)
boxplot((rel_chasing_duration$rel_chasing_duration*900)~as.character(rel_chasing_duration$evo),
        ylab="Chasing activity (sec. spent on chasing females)",names=c("Anc.M.","Evo.M."),
        col=c("skyblue","royalblue"))
dev.off()
png("./FigureS2.png",width = 6,height = 12,units = "cm",res = 600,pointsize = 10)
par(mar=c(4.1,4.1,2,1),las=1)
boxplot((rel_attcop_duration$rel_attcop_duration)~as.character(rel_attcop_duration$evo),
        ylab="Activity on attempted copulation",names=c("anc.\nmales","hot evo.\nmales"),
        col=alpha(c("forestgreen","maroon"),0.75),ylim=c(0,0.31))
dev.off()

#female reproductive dormancy
dorm10=read.delim("./Dormancy_Pop_10.txt",header = T,stringsAsFactors = F)
dorm12=read.delim("./Dormancy_Pop_12.txt",header = T,stringsAsFactors = F)
wilcox.test(dorm10$Freq[1:10],dorm10$Freq[11:13])
wilcox.test(dorm12$Freq[1:10],dorm12$Freq[11:13])

png("./Figure2d.png",width = 6,height = 12,units = "cm",pointsize = 10,res = 600)
par(mar=c(4.1,4.1,2,1),las=1)
boxplot(dorm10$Freq~c(rep("H",10),rep("B",3)),col=c("pink","salmon"),names=c("Anc. F.","Evo. F."),
        ylab="ovarian dormancy level at 10°C")
dev.off()
png("./FigureS3.png",width = 8,height = 12,units = "cm",pointsize = 10,res = 600)
par(mar=c(4.1,4.1,2,1),las=1)
boxplot(dorm12$Freq~c(rep("H",10),rep("B",3)),col=c("pink","salmon"),names=c("Anc. F.","Evo. F."),
        ylab="ovarian dormancy level at 12°C")
dev.off()

####in silico candidate TFs search using RcisTarget####
#download.file("https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc8nr/gene_based/dm6-5kb-upstream-full-tx-11species.mc8nr.feather",
#              destfile = "./dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
motifranking=importRankings("./dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
all_TFs=intersect(unique(motifAnnotations_dmel_v8$TF),conv_background)
all_TFs_tmp=ID_converter(all_TFs,db = ensembl,attributes = c("flybase_gene_id","external_gene_name","description"),
                         filters = "external_gene_name")
all_TFs_tmp=all_TFs_tmp[all_TFs_tmp$external_gene_name%in%all_TFs,]
gene_summary_filtered_tmp=gene_summary_filtered[!duplicated(gene_summary_filtered$FB_ID_new),]
rownames(gene_summary_filtered_tmp)=gene_summary_filtered_tmp$FB_ID_new
all_TFs_out=cbind(all_TFs_tmp[,1:2],gene_summary_filtered_tmp[all_TFs_tmp$flybase_gene_id,1:4],
                  sex_bias=apply(sapply(sex_bias_ID_new[c(1,2)],function(x) ifelse(all_TFs_tmp$flybase_gene_id%in%x,1,0)),1,function(x) ifelse(sum(x)==0,NA,names(which(x==1)))),
                  evolution=apply(sapply(query_ID_new[-c(5,8,9,12)],function(x) ifelse(all_TFs_tmp$flybase_gene_id%in%x,1,0)),1,function(x) ifelse(sum(x)==0,NA,names(which(x==1)))),
                  description=all_TFs_tmp[,3])

#criterion 1: TFs regulating sex-biased genes
motifEnrichmentTable_sex_biased=lapply(sex_bias_ID_conv[5:6],function(x) {
  cisTarget(x,motifRankings = motifranking,motifAnnot = motifAnnotations_dmel_v8,nesThreshold = 3,nCores = 8,
            geneErnMethod = "iCisTarget",aucMaxRank = 0.01*ncol(motifranking))
})
TFs_sex_biased=lapply(motifEnrichmentTable_sex_biased,function(x) {
  tmp=c(x$TF_highConf,x$TF_lowConf)
  genes <- gsub(" \\(.*\\). ", "; ", tmp, fixed=FALSE)
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  genesSplit <- unique(strsplit2(genesSplit, " ")[,1])
  return(genesSplit)
})

#criterion 2: TFs regulating significantly evolving genes
motifEnrichmentTable=lapply(conv_query,function(x) {
  cisTarget(x,motifRankings = motifranking,motifAnnot = motifAnnotations_dmel_v8,nesThreshold = 3,nCores = 8,
            geneErnMethod = "iCisTarget",aucMaxRank = 0.01*ncol(motifranking))
})
TFs=lapply(motifEnrichmentTable,function(x) {
  tmp=c(x$TF_highConf,x$TF_lowConf)
  genes <- gsub(" \\(.*\\). ", "; ", tmp, fixed=FALSE)
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  genesSplit <- unique(strsplit2(genesSplit, " ")[,1])
  return(genesSplit)
})

#criterion 3: TFs evolving in a direction compatible with the changes of their target genes 
index_TFs=Reduce(union, list(Reduce(intersect, list(unlist(TFs_sex_biased),unlist(conv_query[c(6,7)]),unlist(TFs[c(1,2,6,7)]))),
                             Reduce(intersect, list(unlist(TFs_sex_biased),unlist(conv_query[c(10,11)]),unlist(TFs[c(1,2,10,11)]))),
                             Reduce(intersect, list(unlist(TFs_sex_biased),unlist(conv_query[c(1,2)]),unlist(TFs[c(1,2,6,7,10,11)])))))

candidate_TF=all_TFs_out[all_TFs_out$external_gene_name%in%index_TFs,]

######################################
####Part5: simulated data analysis####
######################################
####data input####
dat_files=list.files(path="./simulation1/",pattern=".gpf")[1:700]
dat=lapply(dat_files,function(x) read.table(paste0("./simulation1/",x),header=F,stringsAsFactors = F))
dat_files2=list.files(path="./simulation2/",pattern=".gpf")
dat2=lapply(dat_files2,function(x) read.table(paste0("./simulation2/",x),header=F,stringsAsFactors = F))
names(dat)=strsplit2(dat_files,split = ".gpf")[,1]
names(dat2)=strsplit2(dat_files2,split = ".gpf")[,1]

avg=sapply(dat,function(x) tapply(x[,5],paste0(x[,3],x[,1],x[,2]),mean))
avg_2=sapply(dat2,function(x) tapply(x[,5],paste0(x[,3],x[,1],x[,2]),mean))
dev=sapply(dat,function(x) tapply(x[,5],paste0(x[,3],x[,1],x[,2]),sd))
dev_2=sapply(dat2,function(x) tapply(x[,5],paste0(x[,3],x[,1],x[,2]),sd))
scen=paste0(unique(apply(strsplit2(colnames(avg),"_")[,1:2],1,function(x) paste(x,collapse = "_"))),"_")
scen2=paste0(unique(apply(strsplit2(colnames(avg_2),"_")[,1:2],1,function(x) paste(x,collapse = "_"))),"_")

out_tab=data.frame(Male=(avg["M1100",]-avg["M10",])/dev["M10",]^2,
                   Female=(avg["F1100",]-avg["F10",])/dev["F10",]^2)
out_tab2=data.frame(Male=(avg_2["M1100",]-avg_2["M10",])/dev_2["M10",]^2,
                    Female=(avg_2["F1100",]-avg_2["F10",])/dev_2["F10",]^2)

####statistical testing####
prop_exp1=c()
p.sim=c()
for (j in scen[1:6]){
  p.sim=c(p.sim,t.test(out_tab[grep(j,rownames(out_tab)),1],out_tab[grep(j,rownames(out_tab)),2])$p.value)
  prop_exp1=rbind(prop_exp1,c(sum(out_tab[grep(j,rownames(out_tab)),1]>0&out_tab[grep(j,rownames(out_tab)),2]<0),
                              sum(out_tab[grep(j,rownames(out_tab)),1]>0&out_tab[grep(j,rownames(out_tab)),2]>0),
                              sum(out_tab[grep(j,rownames(out_tab)),1]<0&out_tab[grep(j,rownames(out_tab)),2]>0),
                              sum(out_tab[grep(j,rownames(out_tab)),1]<0&out_tab[grep(j,rownames(out_tab)),2]<0)))
}
p.adjust(p.sim,"bonferroni")
p.adjust(apply(prop_exp1,1,function(x) prop.test(x[1],100,0.34)$p.value),"bonferroni")

prop_exp2=c()
p.sim2=c()
for (j in scen2[1:5]){
  p.sim2=c(p.sim2,t.test(out_tab2[grep(j,rownames(out_tab2)),1],out_tab2[grep(j,rownames(out_tab2)),2])$p.value)
  prop_exp2=rbind(prop_exp2,c(sum(out_tab2[grep(j,rownames(out_tab2)),1]>0&out_tab2[grep(j,rownames(out_tab2)),2]<0),
                              sum(out_tab2[grep(j,rownames(out_tab2)),1]>0&out_tab2[grep(j,rownames(out_tab2)),2]>0),
                              sum(out_tab2[grep(j,rownames(out_tab2)),1]<0&out_tab2[grep(j,rownames(out_tab2)),2]>0),
                              sum(out_tab2[grep(j,rownames(out_tab2)),1]<0&out_tab2[grep(j,rownames(out_tab2)),2]<0)))
}
p.adjust(p.sim2,"bonferroni")
p.adjust(apply(prop_exp2,1,function(x) prop.test(x[1],100,0.34)$p.value),"bonferroni")

####visualization####
png("./Figure4.png",width = 16,height = 8,units = "cm",res=600,pointsize = 8)
par(mfrow=c(1,2),mar=c(5,5,2,2))
plot(NA,xlim=c(0.5,6.5),ylim=c(-2,2),type="n",xaxt="n",xlab="Numbers of sex-specific loci for each sex (# out of 50)",ylab="Normalized phenotypic responses",las=2)
for (j in 1:6){
  boxplot(out_tab[grep(scen[j],rownames(out_tab)),],axes=F,boxwex = 0.4,
          at=c(j-0.25,j+0.25),col=c("royalblue","salmon"),add=T)
}
#text(1,2,"n.s.")
#text(2:6,2,paste0("p = ",scientific(p.adjust(p.sim,"bonferroni")[2:6],2)))
axis(1,at = 1:6,labels = c(0,1,2,5,10,20))
legend("bottomleft",fill = c("salmon","royalblue"),legend = c("Female","Male"),bty="n")
bp=barplot(prop_exp1[,1],names.arg = c(0,1,2,5,10,20),col="darkblue",ylim=c(0,105),
           xlab="Numbers of sex-specific loci for each sex (# out of 50)",ylab = "Fraction of expected evolution")
text(bp[2:6],prop_exp1[2:6,1]+1.5,c("","***","***","***","***"))
dev.off()

png("./FigureS4.png",width = 16,height = 8,units = "cm",res=600,pointsize = 8)
par(mfrow=c(1,2),mar=c(5,5,2,2))
plot(NA,xlim=c(0.5,6.5),ylim=c(-1,1),type="n",xaxt="n",xlab="Numbers of sex-biased loci for each sex (# out of 50)",ylab="Normalized phenotypic responses",las=2)
boxplot(out_tab[grep(scen[1],rownames(out_tab)),],axes=F,boxwex = 0.4,
        at=c(1-0.25,1+0.25),col=c("royalblue","salmon"),add=T)
#text(1,1,"n.s.")
for (j in 1:5){
  boxplot(out_tab2[grep(scen2[j],rownames(out_tab2)),],axes=F,boxwex = 0.4,
          at=c(j+1-0.25,j+1+0.25),col=c("royalblue","salmon"),add=T)
  #  text(j+1,1,paste0("p = ",scientific(p.adjust(p.sim2,"bonferroni")[j],2)))
}
axis(1,at = 1:6,labels = c(0,1,2,5,10,20))
legend("bottomleft",fill = c("salmon","royalblue"),legend = c("Female","Male"),bty="n")
bp2=barplot(c(34,prop_exp2[,1]),names.arg = c(0,1,2,5,10,20),col="darkblue",ylim=c(0,105),
            xlab="Numbers of sex-biased loci for each sex (# out of 50)",ylab = "Fraction of expected evolution")
text(bp2[-1],prop_exp2[,1]+1.5,c("","","***","***","***"))
dev.off()
