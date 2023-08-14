###APRGs R code
##Author：Wang et al. 
##Date:2023-5-10
#KEGG_ARGININE_AND_PROLINE_METABOLISM
#https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/KEGG_ARGININE_AND_PROLINE_METABOLISM
#The Univariate cox file includes information on the 54 APRGs and corresponding samples downloaded from the MsigDB and TCGA database
#Fig1A-Univariate cox regression
library("survival")
library("survminer")
td=read.table("Univariate cox.txt",header = T,row.names = 1,sep = "\t")
pFilter=0.05
outResult=data.frame()
sigGenes=c("status","time") 
for(i in colnames(td[,3:ncol(td)])){ 
  tdcox <- coxph(Surv(time, status) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox) 
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] 
  if(pvalue<pFilter){ 
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outResult,file="UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)
UniCoxSurSigGeneExp=td[,sigGenes] 
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)
write.table(UniCoxSurSigGeneExp,file="UniCoxSurSigGeneExp.txt",sep="\t",row.names=F,quote=F)
tducs <- read.table("UniCoxSurvival.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(tducs)
hr <- sprintf("%.3f",tducs$"HR")
hrLow  <- sprintf("%.3f",tducs$"L95CI")
hrHigh <- sprintf("%.3f",tducs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))

pdf(file="UniCoxSurForestPlot.pdf", width = 6,height = 10)
n <- nrow(tducs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))

xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex);text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)

dev.off()

#pam clustering-FigS1A-B
library(ConsensusClusterPlus)
title <-"C:/Users/WZQ/Desktop/APRGs"
rt=read.table("Univariate cox results.txt", header=T, sep="\t", check.names=F, row.names=1)
dataset=as.matrix(rt)
mads <- apply(dataset, 1, mad)
dataset <- dataset[rev(order(mads))[1:12],]
dim(dataset)
results <- ConsensusClusterPlus(dataset, maxK = 6,
                                reps = 1000, pItem = 0.8,
                                pFeature = 1,  
                                clusterAlg = "pam", 
                                distance = "pearson",
                                title = title,
                                plot = "png")
sample_cluster <- results[[2]]$consensusClass
sample_cluster_df <- data.frame(sample = names(sample_cluster),
                                cluster = sample_cluster)
head(sample_cluster_df)
write.table(sample_cluster_df, file="cluster2.txt", sep="\t", quote=F, col.names=T)

#Fig1B-The survival curve
library(dplyr)
library(survival)
library(survminer)
score_2 <- read.csv("Clinical.csv")  
Sscore <- score_2
sfit <- surv_fit(Surv(time,status==1)~score,data = Sscore)#0 is alive, 1 is dead
sfit
sfit2 <- survfit(Surv(time,status==1)~score,data=Sscore)
sfit2
summary(sfit2,times = 365.25)
survdiff(Surv(time,status==1)~score ,data = Sscore)#log rank test
pdf("LUAD-APRGs.pdf", height = 5.5,width = 6,onefile = FALSE)
ggsurvplot(sfit2,conf.int = F,pval = T,
           risk.table = T,surv.median.line = "hv",
           legend = c(0.8, 0.83),
           legend.labs=c("cluster1","cluster2"),break.time.by=1000,
           legend.title="Group",palette = c("#2E9FDF","#EFC000"),
           title="TCGA-LUAD",ggtheme=theme_classic()+
             theme(plot.title = element_text(size=12,hjust=0.5)))+
  xlab("Time(days)")
dev.off()

#Fig1C
library(umap)
library(ggplot2)
df <- read.csv(file = "PCA.csv") 
df_umap <- df[,colnames(df)!='label']
df_umap <- data.frame(t(apply(df_umap,1,function(v){(v-mean(v,na.rm=T))/sd(v,na.rm=T)})),stringsAsFactors=F)
df_umap[is.na(df_umap)] <- 0
umap <- umap(df_umap,method='naive',n_neighbors = 10)
head(umap$layout)
df1 <- data.frame(umap$layout)
df1$label <- df$label 
colnames(df1) <- c('X','Y','label') 
p <- ggplot(df1, aes(x=X, y=Y, colour=label)) + geom_point(size=1.2)+
  xlab("umap_1")+ 
  ylab("umap_2")+
  theme(  panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title=element_blank(), 
          panel.border = element_blank(),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          panel.background = element_blank())+
  scale_colour_manual(values=c( "#2E9FDF","#EFC000"))
p
ggsave(p,filename = "umap.pdf",width = 5,height = 4)

#Fig1D,Fig2A,Fig2C,Fig3B-pheatmap
library(pheatmap)          
inputFile="exp.txt"      
groupFile="group.txt" 
outFile="pheatmap1.pdf"      
rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)    
ann=read.table(groupFile,header=T,sep="\t",row.names=1,check.names=F)   
rt <- as.matrix(rt)
dat <- t(scale(t(rt)))
dat[dat>4]=4
dat[dat<(-4)]=-4
Clustercolor<-c("#2E9FDF","#EFC000")#"#4682B4","#B22222"
names(Clustercolor)<-c("Cluster1","Cluster2")
ann_colors <- list(Group=Clustercolor) 
pdf(file=outFile,width=6.5,height=7)
pheatmap(dat,
         annotation_col = ann,
         cluster_cols = F,
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         show_colnames = F,
         fontsize = 8,
         fontsize_row=7.5,
         fontsize_col=6,
         cluster_rows = T,
         show_cluster=F,
         annotation_colors = ann_colors,
         treeheight_row = 4.5,
         )
dev.off()

#Fig2B-somatic cell mutation
library(maftools)
options(stringsAsFactors = F)
laml = read.maf(maf = 'combined_maf_value_final_1.txt')#Matrix of merged somatic mutation data
pdf(file="maf_summary1.pdf",width =12,height=7)
plotmafSummary(maf = laml,addStat = 'median')
dev.off()
col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')#Modify the colour scheme
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit')
genes <- c("TP53","KRAS","KEAP1","STK11","EGFR")
pdf(file="TMB_cluster1.pdf",width =8,height=6)
oncoplot(maf = laml, genes = genes,colors = col)
dev.off()
#Calculating the TMB score
tmb_table=tmb(maf = laml)
tmb_table=tmb(maf = laml,logScale = F)
write.csv(tmb_table,"tmb_score_results.csv")

#limma-Variance analysis
library(limma)
library(dplyr)
df <- read.table("TPM.normalize.txt", header = T, sep = "\t", row.names = 1, check.names = F)
list <- c(rep("cluster1",400),rep("cluster2",100)) %>% factor(., levels = c("cluster1", "cluster2"), ordered = F)
list <- model.matrix(~factor(list)+0)
colnames(list) <- c("cluster1", "cluster2")
df.fit <- lmFit(df, list) 
df.matrix <- makeContrasts(cluster1 - cluster2 , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
nrDEG = na.omit(tempOutput)
diffsig <- nrDEG  
write.csv(diffsig, "all.limmaOut.csv")
foldChange = 1
padj = 0.05
All_diffSig <-diffsig[(diffsig$adj.P.Val < padj & (diffsig$logFC>foldChange | diffsig$logFC < (-foldChange))),]
write.csv(All_diffSig, "all.diffsig.csv")
dim(All_diffSig)

#Fig2D-Go enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
hsa <- org.Hs.eg.db
keytypes(hsa)
head(keys(hsa,keytype = "SYMBOL"))
gene <- read.csv(file ="all.diffsig.csv")
go.all <- enrichGO(gene=gene$symbol,OrgDb = org.Hs.eg.db,ont = 'ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,keyType ='SYMBOL' )
go.all
dim(go.all[go.all$ONTOLOGY=='BP',]);dim(go.all[go.all$ONTOLOGY=='CC',]);dim(go.all[go.all$ONTOLOGY=='MF',])
write.csv(go.all@result,'gene_go.all.result.csv',row.names = F)
pdf("GO-cluster.pdf",width = 10,height =8 )
dotplot(go.all, showCategory = 10,size=NULL,font.size = 8,title = "GO enrichment",split="ONTOLOGY",label_format =100)+facet_grid(ONTOLOGY~.,scales = "free")
dev.off()

#Fig3A-Estiamte 
library(estimate)
in.file <- 'TPM.normalize.txt' 
outfile2E <- 'ESTIMATE_input.gct' 
outputGCT(in.file, outfile2E)
filterCommonGenes(input.f= in.file, output.f= outfile2E, id="GeneSymbol")
estimateScore("ESTIMATE_input.gct", "ESTIMATE_score.gct")
plotPurity(scores="ESTIMATE_score.gct", samples="s516")
ESTIMATE_score <- read.table("ESTIMATE_score.gct", skip = 2,
                             header = TRUE,row.names = 1) 
ESTIMATE_score <- ESTIMATE_score[,2:ncol(ESTIMATE_score)] 
ESTIMATE_score 
write.table(ESTIMATE_score,file = "ESTIMATE_score.txt",quote = F,sep = "\t")
library(ggpubr)
library(cowplot)
df <- read.table("estimate results.txt",header = T,sep = "\t",row.names = 1,check.names = F)
my_comparisons <- list( c("cluster1", "cluster2") )
p1 <- ggviolin(df,x = "group", y = "StromalScore",fill= "group", palette = c("#2E9FDF","#EFC000"),
               add = "boxplot", add.params = list(fill = "white"), legend="none",
               order = c("cluster1", "cluster2"),
               error.plot = "errorbar") + xlab("")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
p2 <- ggviolin(df,x = "group", y = "ImmuneScore",fill= "group", palette = c("#2E9FDF","#EFC000"),
               add = "boxplot", add.params = list(fill = "white"),  legend="none",
               order = c("cluster1", "cluster2"),
               error.plot = "errorbar") + xlab("")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
p3 <- ggviolin(df,x = "group", y = "ESTIMATEScore",fill= "group", palette = c("#2E9FDF","#EFC000"),
               add = "boxplot", add.params = list(fill = "white"),  legend="none",
               order = c("cluster1", "cluster2"),
               error.plot = "errorbar") + xlab("")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
p4 <- ggviolin(df,x = "group", y = "TumorPurity",fill= "group", palette = c("#2E9FDF","#EFC000"),
               add = "boxplot", add.params = list(fill = "white"),  legend="none",
               order = c("cluster1", "cluster2"),
               error.plot = "errorbar") + xlab("")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")
estimate <- plot_grid(p1,p2,p3,p4,ncol = 2)
estimate
ggsave("ESTIMATE.pdf",plot=estimate,height = 6.5,width = 8.5)

#Fig3C
#Xcell
library(xCell)
exprMatrix = read.table("TPM_tumor_matrix.txt",header=TRUE,sep = "\t",row.names = 1,check.names = F)
exp <- as.matrix(xCellAnalysis(exprMatrix))
write.table(file="Xcell results.txt",exp,sep = "\t",row.names = F)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggsci)
df <- read_tsv("Xcell DC.txt") %>% pivot_longer(- group)
a=ggboxplot(
  df, x = "name", y = "value",
  color = "group", palette = c("aaas"),#"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
  add = "jitter" )+ylab("TIDE score")+xlab("")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
  label ="p.format")
ggsave("XCELL-DC.pdf",plot = a,width = 5.5,height = 5)

#Fig3D-Cibersort
Source(Cibersort.R)
LM22.file <- "LM22.txt"
TCGA_exp.file<- "cibersort,exp.txt"
TCGA_TME.results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 1000, QN = F)  
write.csv(TCGA_TME.results, "TCGA_CIBERSORT_Results.csv")
df <- read_tsv("Cibersort.txt") %>% pivot_longer(- group)
location <- df %>% group_by(name) %>% slice_max(value)
location$x <- seq(1,22,by=1)
head(location,3)
ggplot(df,aes(x = name, y = value,fill=group))+
  geom_violin(scale = "width",alpha=0.8,width=0.5,size=0.5)+ 
  scale_fill_manual(values = c("#2E9FDF","#EFC000"))+        
  stat_compare_means(aes(group=group),                       
                     method = "wilcox.test",
                     paired = F,                             
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     label.y = location$value+0.02,      
                     size=4.5)+ 
  geom_segment(data=location,                                 
               aes(x=x,y=value,xend=x+0.2,yend=value),size=1)+
  xlab("")+                                                   
  ylab("Fraction")+                                           
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(size=1.2),                  
        axis.text.x = element_text(angle=60,size=10,vjust = 1,hjust =1,color = "black",face = "bold"),                                                   
        axis.text.y = element_text(size =10,face = "bold"),
        legend.position = c(0.9,0.85) )+
  geom_line(position=position_dodge(0.5))+
  stat_summary(fun = "median", geom = "point", position=position_dodge(0.5),colour="white",size=3)

#fig3E-ssGSEA
library(GSVA)
gs = read.csv("cellmarker.csv", stringsAsFactors = FALSE, check.names = FALSE)
a = read.table("TPM.normalize.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1, sep = "\t")
a = as.matrix(a)
gs = as.list(gs)
gs = lapply(gs, function(x) x[!is.na(x)])
ssgsea = gsva(a, gs, method = "ssgsea", kcdf='Gaussian',abs.ranking=TRUE) 
ssgsea=t(scale(t(ssgsea)))
write.csv(ssgsea, "ssgsea.csv")
library(ggplot2)
library(tidyverse)
df <- read_tsv("ssGSEA.txt") %>% pivot_longer(- group)
ggplot(df, aes(x = name, y = value))+ 
  labs(y="score",x= "",title = "group")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c( "#51C4D3","#EC7696"))+
  theme_classic() + 
  stat_compare_means(aes(group =  group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = F)+
  theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
        axis.title = element_text(size = 12,color ="black"), 
        axis.text = element_text(size= 12,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12))

#Fig4A-B,Immunological escape characteristics
df <- read_tsv("ssGSEA.txt") %>% pivot_longer(- group)
ggboxplot(df, x="name", y="value", color = "group", 
          ylab="Gene expression",
          xlab="",
          legend.title="Antigen presentation",
          palette = c("#3170A7","#F3BF4C"),
          width=0.6, add = "none")+
  theme_bw()+
  rotate_x_text(60)+
  stat_compare_means(aes(group=group),
                     method="wilcox.test",
                     label = "p.signif")+
  theme(legend.position = "top")

#Fig4C,FigS4A-D
d<-read.csv("all.limmaOut.csv")  
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
gene<-str_trim(d$symbol,"both") 
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=d$logFC, 
                      SYMBOL = d$symbol) 
gene_df <- merge(gene_df,gene,by="SYMBOL")
geneList<-gene_df $logFC 
names(geneList)=gene_df $ENTREZID 
geneList=sort(geneList,decreasing = T) 
geneList
kegmt<-read.gmt("c2.cp.kegg.v2022.1.Hs.entrez.gmt") 
KEGG<-GSEA(geneList,TERM2GENE = kegmt)
library(GseaVis)
dotplotGsea(data = KEGG,topn = 10,
            order.by = 'NES',
            add.seg = T)
a=gseaNb(object = KEGG,
         geneSetID = 'KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY',addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0)
ggsave("GSEA-Tcell receptor.pdf",plot = a,height = 6,width = 8.4)

#Fig4D
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
A  <- read.table("Immune escape.txt",header = T,sep = "\t",row.names = 1)
dat <- t(scale(t(A)))
dat[dat>2]=2
dat[dat<(-2)]=-2
dat <- as.matrix(dat)
samples <- rep(c('Cluster1', 'Cluster2'), c(400, 100)) #定义样本分组信息  
Group = factor(rep(c("Cluster1","Cluster2"),times = c(400, 100)))#分组信息，用于热图分割
Group = factor(Group,levels = c("Cluster1","Cluster2"))
top_annotation = HeatmapAnnotation(group=Group,#列注释
                                   border = T,
                                   show_annotation_name = T,
                                   col = list(group=c('Cluster1'="#7097A8",
                                                      'Cluster2'="#AE6378")))
df = data.frame(type = c(rep("Immune ihibtor", 16), rep("Immune stimulator", 22),rep("Chemokines and receptors", 51),rep("Interferons", 8),
                         rep("Interluekins and receptors", 46),rep("Others", 22))
)
split = rep(1:6, c(16,22,51,8,46,22))
col_fun <- colorRamp2(c(-2,0,2), c("navy", "#FFFFFF", "firebrick3"))  
pdf("finalheatmap1.pdf",height = 20,width = 6)
Heatmap(dat,cluster_rows = F,
        col = colorRampPalette(c("navy", "#FFFFFF", "firebrick3"))(100),      
        show_row_names = T,
        top_annotation = top_annotation,
        column_split = Group,
        column_title = NULL,
        show_column_names = F,
        border = T,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 8),
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill =2:2),
                                                         labels = c("Immune ihibtor","Immune stimulator",
                                                                    "Chemokines and receptors","Interferons","Interluekins and receptors","Others"), 
                                                         labels_gp = gpar(col = "white", fontsize = 10))),
        
        row_split = split,
        show_heatmap_legend = FALSE
)
lgd <- Legend(col_fun = col_fun, 
              at = c(-2,-1, 0,1,2),
              # labels = c("low", "zero", "high"),
              title = "Expression level",
              legend_height = unit(2, "cm"),
              title_position = "topcenter",
              title_gp = gpar(fontsize = 8),
              labels_gp = gpar(fontsize = 8),
              direction = "horizontal",
              grid_height = unit(4, "mm")
)
draw(lgd, x = unit(0.9, "npc"), y = unit(0.54, "npc"))
dev.off()

#Fig4E
library(ggplot2)
library(ggrepel)
df <- read.csv("all.limmaOut.csv")
p <- ggplot(data = df,aes(x = logFC, y = FDR)) + 
  geom_point(aes(colour = cd,size = abs(logFC)),alpha=0.5) +
  scale_color_manual(values=c("down" = "#80B1D3","else" = "grey","up" = "#FB8072"))+
  scale_x_continuous("logFC",limits = c(-0.5,0.5),breaks = seq(-0.5,0.5,0.2),labels = seq(-0.5,0.5,0.2)) +
  scale_y_continuous("-log10 (FDR)")+
  geom_vline(xintercept=c(-0.1,0.1),lty=2,col="black",lwd=1) +
  theme_bw()+
  theme(
    legend.background=element_blank(), legend.key=element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )
p
df1=df%>%filter(abs(logFC) > 0.1)
p+geom_text(data=df1,mapping = aes(label=X))
p+geom_label(aes(label=X),df1,alpha=0,nudge_y = 3)
p+ggrepel::geom_text_repel(
  aes(label=X),df1
)
p+ggrepel::geom_text_repel(
  aes(label=X,color=cd),df1,
  size = 4, 
  box.padding = 0.5, 
  point.padding = 0.8, 
  min.segment.length = 0.5, 
  segment.color = "black",
  show.legend = F)
p+ggrepel::geom_label_repel(
  aes(label=X),df1
)+theme(plot.title = element_text(size = 12,color="black",hjust = 0.5))+
  labs(title = "Cluster 1 vs. Cluster 2")

#Fig5A
library(ggpubr)
df <- read.table("GEP_score.txt",header = T,sep = "\t")
my_comparisons <- list( c("cluster1", "cluster2") )
ggviolin(df,x = "group", y = "GEP_score",fill= "group", palette = c("#2E9FDF","#EFC000"),
         add = "boxplot", add.params = list(fill = "white"), 
         order = c("cluster1", "cluster2"),
         error.plot = "errorbar") + xlab("")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")

#Fig5B
ggboxplot(
  df, x = "group", y = "PD-L1",
  ###颜色分类以subgroup列为准，后面是颜色参数，拿不准可以用帮助文件查看 
  color = "group", palette = c("aaas"),#"npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
  ###在箱式图上加点
  add = "jitter" )+ylab("PD-L1 Expression")+xlab("")+
  ####后面是增加比较数据并给标尺
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test",
                     label ="p.format")

#Fig5D-E
library(ggplot2)
library(reshape2)
library(ggfittext)
library(tidyverse)
df <- read.table("TIDE.txt",header = T,sep = "\t",row.names = 1)
df = melt(df) 
colnames(df) <- c("X","Tumor_Type", "value", "label")
df$Responder <- gsub("Responder","R",df$Responder)
df$Responder <- gsub("No.R","NR",df$Responder)
ggplot( df, aes( x = X, weight = value, fill = Responder))+
  geom_bar( position = "stack")+xlab("X-squared = 71.147, df = 4, p-value = 1.3e-14")+
  ylab("Percentage")+
  annotate("text",x="cluster1",y=19.625,label="39.25%")+
  annotate("text",x="cluster1",y=69.625,label="60.75%")+
  annotate("text",x="cluster2",y=15,label="30%")+
  annotate("text",x="cluster2",y=65,label="70%")+
  theme(panel.grid = element_blank())+theme_bw()+
  scale_fill_manual( values = c("#51C4D3","#EC7696"))
ggplot( df, aes( x = X, weight = value, fill = Tumor_Type))+
  geom_bar( position = "stack")+xlab("X-squared = 136.12, df = 6, p-value < 2.2e-16")+
  ylab("Percentage")+
  annotate("text",x="cluster1",y=4.375,label="8.75%")+
  annotate("text",x="cluster1",y=18.375,label="19.25%")+
  annotate("text",x="cluster1",y=64,label="72%")+
  annotate("text",x="cluster2",y=9.5,label="19%")+
  annotate("text",x="cluster2",y=31.5,label="25%")+
  annotate("text",x="cluster2",y=72,label="56%")+
  theme(panel.grid = element_blank())+theme_bw()+
  scale_fill_manual( values = c("#51C4D3","#EC7696","#F9D770"))

#Fig6A-B
library(survival)
library(rms)
LXT <- read.table("line graph.txt",header = T,sep = "\t",row.names = 1)
LXT$gender <- factor(LXT$gender,
                     levels = c(0,1),
                     labels = c("male", "female"))

dd <- datadist(LXT )
options(datadist = "dd")
coxph <- coxph(Surv(time, status) ~ age + gender+stage+riskscore, data = LXT )
coxph
survival <- Survival(cph)
survival1 <- function(x) survival(365, x)
survival2 <- function(x) survival(730, x)
survival3 <- function(x) survival(1825, x)
library(regplot)
regplot(coxph,
        observation=LXT[6,], 
        points=TRUE,
        plots=c("density","no plot"),
        failtime = c(365,1095,1825), 
        odds=F,
        droplines=F,
        leftlabel=T,
        prfail = TRUE, 
        showP = T,
        rank="range", 
        interval="confidence",
        title="Cox regression"
)
rcorrcens(Surv(time,status) ~ predict(coxph), data =  LXT)
psm <- psm(Surv(time, status) ~ age + gender+stage+riskscore, data = LXT , x = TRUE,
           y = TRUE, dist='lognormal')
psm
cal1 <- calibrate(psm, 
                  cmethod='KM', 
                  method="boot", 
                  u=365, 
                  m=164, 
                  B=1000) 
plot(cal1,lwd=2,lty=1,
     conf.int=T,
     errbar.col="blue",
     col="red", 
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 1-Year OS",
     ylab="Actual 1-Year OS (proportion)",
     subtitles = F)
cox1 <- cph(Surv(OS.time,OS) ~ age + gender +  tumor_stage + riskScore2,
            surv=T,x=T, y=T,time.inc = 1*365,data=riskScore_cli2)
cal <- calibrate(psm, cmethod="KM", method="boot", u=1*365,
                 m= 100, B=1000)
#3-year
cox3 <- cph(Surv(OS.time,OS) ~ age + gender +  tumor_stage + riskScore2,
            surv=T,x=T, y=T,time.inc = 1*365*3,data=riskScore_cli2)
ca3 <- calibrate(psm, cmethod="KM", method="boot", u=1*365*3,
                 m= 100, B=1000)
#5-year
cox5 <- cph(Surv(OS.time,OS) ~ age + gender +  tumor_stage + riskScore2,
            surv=T,x=T, y=T,time.inc = 1*365*5,data=riskScore_cli2)
ca5 <- calibrate(psm, cmethod="KM", method="boot", u=1*365*5,
                 m= 100, B=1000)
pdf("Correction curve.pdf",width = 6.5,height =5.5 )
plot(cal,lwd=2,lty=1,errbar.col="black",
     xlim = c(0.1,1) ,
     ylim = c(0.1,1), 
     xlab ="Nomogram-Predicted Probability of 1,3,5Year Survival",
     ylab="Actual 1,3,5Year Survival",col="#229453",sub=F)
plot(ca3,
     add = T ,
     lwd=2,lty=1,errbar.col="black",col='#FF1493',sub=F)
plot(ca5,
     add = T ,
     lwd=2,lty=1,errbar.col="black",col="#2775B6",sub=F)
legend("bottomright",
       c(paste0(" 1 year "), 
         paste0(" 3 year "), 
         paste0(" 5 year ")),
       col=c( "#229453",'#FF1493', "#2775B6"),
       lty=1, lwd=2,bty = "n")
dev.off()

#Fig6F
library(rms)
library(ggrisk)
library(ggplot2)
YHMX <- read.table("Clinical.txt",header = T,sep = "\t",row.names = 1)
fit <- cph(Surv(time,status)~CPS1+SMS, YHMX )
ggrisk(fit,new.data = YHMX ,
       cutoff.value = 'cut',
       cutoff.x = 150,
       cutoff.y = -1,
       code.0 = 'Alive',
       code.1 = 'Dead',
       code.highrisk = 'High Risk',
       code.lowrisk = 'Low Risk',
       title.A.ylab =  'Risk Score',
       title.B.ylab = 'Survival Time(year)',
       title.A.legend='Risk Group',
       title.B.legend='Status',
       title.C.legend='Expression')

#Fig6G
rm(list = ls())
library(timeROC)
library(survival)
df <- read.table("ROC.txt",header = T,sep = "\t",row.names = 1)
ROC <- timeROC(T=df$time,   
               delta=df$status,   
               marker=df$riskscore,   
               cause=1,                
               weighting="marginal",   
               times=c(1, 2, 3),       
               iid=TRUE)
pdf("ROC-risk.pdf",width = 6,height = 6)
plot(ROC, 
     time=1, col="#229453", lwd=4, title = "")   
plot(ROC,
     time=2, col='#FF1493', add=TRUE, lwd=4)    
plot(ROC,
     time=3, col="#2775B6", add=TRUE, lwd=4)
legend("bottomright",
       c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],4)), 
         paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],4)), 
         paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],4))),
       col=c("#229453",'#FF1493', "#2775B6"),
       lty=1, lwd=4,bty = "n")
dev.off()

#Multivariate cox regression
meta <- read.table("Multivariate cox.txt",header=T,sep="\t",row.names=1)
res.cox <- coxph(Surv(time, status) ~ CPS1+SMS , data =  meta)
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
multi_res
write.table(file="multivariate_cox_result.txt",multi_res,quote=F,sep="\t")
outResult=data.frame() 
tdcoxSummary = summary(res.cox)
outResult=rbind(outResult,
                cbind(
                  HR=tdcoxSummary$conf.int[,"exp(coef)"],
                  L95CI=tdcoxSummary$conf.int[,"lower .95"],
                  H95CI=tdcoxSummary$conf.int[,"upper .95"],
                  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
)
outFile="cluster-Multivariate.pdf" 
rt=outResult
gene=rownames(rt)
hr=sprintf("%.3f",rt$"HR")
hrLow=sprintf("%.3f",rt$"L95CI")
hrHigh=sprintf("%.3f",rt$"H95CI")
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file=outFile, width = 6, height =4.5)
n=nrow(rt)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

xlim = c(0,3)
par(mar=c(4,2,1.5,1.5))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#Fig7A
library(gghalves)
library(gglayer)
library(introdataviz)
df <- read.table("Gene Expression.txt")
color=c("#F18D91", "#639CCF")
my_comparisons <- list( c("Normal", "Tumor") )
ggplot(df, aes(x=group, y=ACLY)) +
  geom_flat_violin(aes(fill=group),position=position_nudge(x=.25),color="black") +
  scale_color_manual(values=rev(color))+
  scale_fill_manual(values=rev(color))+
  geom_jitter(aes(color=group), width=0.1,size=0.5) +
  geom_boxplot(width=.1,position=position_nudge(x=0.25),fill="white",size=0.5) +
  theme_bw()+
  xlab("")+
  ylab("The expression of SMS")+
  stat_compare_means(comparisons = my_comparisons,aes(group=group),
                     method="wilcox.test",bracket.size = 0.6,size=5,
                     label = "p.signif")

#FigS6-drug sensitivity
library(oncoPredict)
library(data.table)
library(ggplot2)
library(ggpubr)
load("GDSC2_Expr.Rdata")
load("GDSC2_Res.Rdata")
group = c(rep('cluster1',400),rep('cluster2',100))
testExpr <- read.table("TPM.normalize.500.txt",header = T,sep = "\t",row.names = 1,check.names = F)
testExpr <- as.matrix(testExpr)
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
res <- fread("calcPhenotype_Output\\DrugPredictions.csv")
res <- as.data.frame(res)
rownames(res) <- res$V1
res <- res[,-1]
for (i in 1:198) {
  Cam <- as.data.frame(res[,i])
  colnames(Cam) <- "senstivity"
  Cam$Risk <- group
  boxplot=ggboxplot(Cam, x="Risk", y="senstivity", fill="Risk",
                    xlab="",
                    ylab=paste0(colnames(res)[i], " senstivity (IC50)"),
                    legend.title="",
                    palette=c("#4995c6","#EE4431")
  )+
    stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("cluster1","cluster2")),method="wilcox.test")+#???Ӽ???
    pdf(file=paste0(colnames(res)[i], ".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}

#FigS7B
library(ggbeeswarm)
library(scales)
library(ggplot2)
df <- read.table("Gene Expression2.txt")
comparisons <- list(c("Normal","Tumor"))
ggplot(df,aes(group,ACLY,color=group))+
  stat_boxplot(geom = "errorbar",width=.2,color="black")+
  geom_boxplot(color="black",outlier.alpha = 0)+
  geom_beeswarm(priority='density',size=1.0,cex = 1.5)+
  scale_color_manual(values=c(  "blue","red","#999999" ))+
  theme_bw()+
  theme(legend.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid = element_blank())+
  xlab("")+
  ylab("SMS expression")+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")
















