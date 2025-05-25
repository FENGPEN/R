####差异和WCGNA####

#引用包
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)
library(randomcoloR)
if (!require('R.utils')) install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )

#参数设置
GSE="normalize"    #表达矩阵文件名称，不用后缀
C="Low"              #正常控制组名称
P="High"              #疾病实验组的名称
Ccol="blue2"        #热图注释条正常组颜色
Pcol="red2"         #热图注释条疾病组颜色
lowcol="blue"     #热图格子低表达量颜色
midcol="white"     #热图格子中表达量颜色
highcol="red"      #热图格子高表达量颜色
num=40             #热图分别展示的上下调基因数目

rt=read.table(paste0(GSE,".txt"),sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

######################################################################################################################1.数据准备
#分组
sample=read.table("分组.txt",sep="\t",header=F,check.names=F,row.names = 1)
rt=rt[,rownames(sample)]
afcon=sum(sample[,1]==C)
#判断原始数据是否去了log
max(rt)
if(max(rt)>50) rt=log2(rt+1)     #rt最大值大于50则取log

#使用normalizeBetweenArrays进行矫正，矫正后赋值为rt1
rt1=normalizeBetweenArrays(as.matrix(rt))

#未标准化,mar调整画布范围下左上右，自行调整哈
cols=rainbow(ncol(rt)) ###针对24个样本，设置颜色，整体呈现彩虹色
pdf(file = "1.raw.pdf",width=5,height = 4)
par(cex = 0.7,mar=c(8,8,8,8))
if(ncol(rt)>40) par(cex = 0.5,mar=c(8,8,8,8))   ###设置字体大小
boxplot(rt,las=2,col =cols ) ###绘图
dev.off()

#标准化
cols=rainbow(ncol(rt1)) ###针对24个样本，设置颜色，整体呈现彩虹色
pdf(file = "1.nor.pdf",width=5,height = 4.5)
par(cex = 0.5,mar=c(8,8,8,8))
if(ncol(rt1)>40) par(cex = 0.5,mar=c(8,8,8,8))   ###设置字体大小
boxplot(rt1,las=2,col =cols ) ###绘图
dev.off()

#保存标准化后结果
rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file=paste0("1.","norexp_",GSE,".txt"),sep="\t",quote=F,col.names = F)

#保留原始结果
rt3=rbind(ID=colnames(rt),rt)
write.table(rt3,file=paste0("1.","rawexp_",GSE,".txt"),sep="\t",quote=F,col.names = F)

#######################################################################################################################2.差异分析
data=rt
#data=rt1

conData=data[,as.vector(colnames(data)[1:afcon])]
aftreat=afcon+1
treatData=data[,as.vector(colnames(data)[aftreat:ncol(data)])]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#limma差异标准流程
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Diff=topTable(fit2,adjust='fdr',number=length(rownames(data)))
#保存所有基因的差异结果
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file=paste0("2.","DIFF_all.xls"),sep="\t",quote=F,col.names=F)
diffSig=Diff[with(Diff, (abs(logFC)>1 & adj.P.Val < 0.05 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file=paste0("2.","DIFF_less.xls"),sep="\t",quote=F,col.names=F)

#展示差异最大的前num个基因
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(2*num)){
  afGene=diffGene[c(1:num,(diffLength-num+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=rt[afGene,]
#分组标签
Type=c(rep(C,conNum),rep(P,treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
#分组标签的注释颜色
anncolor=list(Type=c(C=Ccol,P=Pcol))
names(anncolor[[1]])=c(C,P)

pdf(file=paste0("3.", "DIFF_heatmap.pdf"),height=7,width=6)
pheatmap(afExp,                                                                      #热图数据
         annotation=Type,                                                            #分组
         color = colorRampPalette(c(lowcol,midcol,highcol))(50),     #热图颜色
         cluster_cols =F,                                                           #不添加列聚类树
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 4,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=anncolor
)
dev.off()

#火山图差异标记颜色标准设置
adjP=0.05
aflogFC=0.5
Significant=ifelse((Diff$P.Value<adjP & abs(Diff$logFC)>aflogFC), ifelse(Diff$logFC>aflogFC,"Up","Down"), "Not")
#开始绘制
p = ggplot(Diff, aes(logFC, -log10(P.Value)))+
  geom_point(aes(col=Significant),size=3)+
  scale_color_manual(values=c(pal_npg()(2)[2], "#838B8B", pal_npg()(1)))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  geom_hline(aes(yintercept=-log10(adjP)), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=aflogFC), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=-aflogFC), colour="gray", linetype="twodash",size=1)
#查看，不添加标记可以直接保存
p
#添加基因点名称标记，按照,可自行根据差异分析的结果进行标记
point.Pvalue=0.01
point.logFc=2
#继续绘制
Diff$symbol=rownames(Diff)
pdf(paste0("3.", "DIFF_vol.pdf"),width=6.5,height=6)
p=p+theme_bw()
for_label <- Diff %>% 
  filter(abs(logFC) >point.logFc & P.Value< point.Pvalue )
p+geom_point(size = 1.5, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black",
    label.size =0.1
  )
dev.off()

###3.GSEA分析
deg=Diff
logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

#开始GSEA富集分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "fdr" )
#保存GSEA结果
GSEAOUT=as.data.frame(kk2@result)
write.table(GSEAOUT,file="4.GSEAOUT.xls",sep="\t",quote=F,col.names=T)
#看一看有没有富集到结果
head(GSEAOUT)

#排序后分别取GSEA结果的前5个和后5个
col=distinctColorPalette(100)        #随机颜色，不满意可以多跑几遍
num=5
pdf(paste0("4.","down_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)],color = col[1:(num)])
dev.off()
pdf(paste0("4.","up_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)],color = col[1:(num)])
dev.off()
#排序后取前5个和后5个一起展示
num=5
pdf(paste0("4.","all_GSEA.pdf"),width = 10,height = 10)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))],color = col[1:(num*2)])
dev.off()

#生成目录存放批量GSEA单图，取前5个和后5个分别绘制
num=5
afname=rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))]
afdDir <- paste0(getwd(),"/4.GSEA")           #路径必须中文
dir.create(afdDir)
for (j in afname) {
  #使用GseaVi包进行单个gsea的绘制，https://github.com/junjunlab/GseaVis
  if (!require('GseaVis')) devtools::install_github("junjunlab/GseaVis")
  library(GseaVis)
  dd=gseaNb(object = kk2,
            geneSetID = j,
            addPval = T,
            pvalSize = 6,            #p与nes字体大小
            pvalX = 0.8,pvalY = 0.7,      #p与nes坐标
            pCol = 'black',)
  pdf(paste0(afdDir,"/",j,".pdf"),width = 10,height = 10)
  print(dd)
  dev.off()
}

#山脊图，展示全部
library(stringr)
pdf(paste0("4.","all_ridgeplot_GSEA.pdf"),width = 6,height = 20)
ridgeplot(kk2)
dev.off()

#############################################################################################################5.WGCNA

###准备
afdir <- paste0(getwd(),"/5.WGCNA")           #路径必须中文
dir.create(afdir)

traitData=sample
traitData[,2]=traitData[,1]
traitData[,1]=ifelse(traitData[,1]==C,1,0)
traitData[,2]=ifelse(traitData[,2]==P,1,0)
#修改性状名称
colnames(traitData)=c(C,P)


###############正式分析
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set
fpkm = read.table(paste0("1.rawexp_",GSE,".txt"),header=T,comment.char = "",check.names=F)#########file name can be changed#####数据文件名，根据实际修改，如果工作路径不是实际数据路径，需要添加正确的数据路径
rownames(fpkm)=fpkm[,1]
dim(fpkm)
names(fpkm)
datExpr0 = as.data.frame(t(fpkm[,-1]))
names(datExpr0) = fpkm[,1];
rownames(datExpr0) = names(fpkm[,-1])

datExpr0

##check missing value
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

##filter
meanFPKM=0.5  ####the threshold can be changed---过滤标准
n=nrow(datExpr0)
datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)
datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM]  # for meanFpkm in row n+1 and it must be above what you set--select meanFpkm>opt$meanFpkm(by rp)


filtered_fpkm=t(datExpr0)
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]="sample"
head(filtered_fpkm)
write.table(filtered_fpkm, file=paste0(afdir,"/FPKM_filter.xls"),row.names=F, col.names=T,quote=FALSE,sep="\t")

sampleTree = hclust(dist(datExpr0), method = "average")


pdf(file =paste0(afdir,"/1_sampleClustering.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

### Plot a line to show the cut
##abline(h = 15, col = "red")##是否选择剪切
dev.off()


### Determine cluster under the line
##clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
##table(clust)


### clust 1 contains the samples we want to keep.
##keepSamples = (clust==1)
##datExpr0 = datExpr0[keepSamples, ]


####载入性状数据
#Loading clinical trait data
for (df in colnames(traitData)) {
  traitData[,df]=traitData[,df]/max(traitData[,df])
  print(sd(traitData[,df]))
}
max(traitData)
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits) 
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.

#sizeGrWindow(12,12)
pdf(file=paste0(afdir,"/2_Sample dendrogram and trait heatmap.pdf"),width=12,height=11)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()







#############################network constr########################################
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
allowWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(1:30)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
#sizeGrWindow(9, 5)
pdf(file=paste0(afdir,"/3_Scale independence.pdf"),width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######chose the softPower
softPower =sft$powerEstimate
#显示软阈值
print(softPower)

adjacency = adjacency(datExpr0, power = softPower)

##### Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file=paste0(afdir,"/4_Gene clustering on TOM-based dissimilarity.pdf"),width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file=paste0(afdir,"/5_Dynamic Tree Cut.pdf"),width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file=paste0(afdir,"/6_Clustering of module eigengenes.pdf"),width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25######模块剪切高度
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#sizeGrWindow(12, 9)
pdf(file=paste0(afdir,"/7_merged dynamic.pdf"), width = 9, height = 6.5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Save module colors and labels for use in subsequent parts
#save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = "networkConstruction-stepByStep.RData")




##############################relate modules to external clinical triats######################################
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file=paste0(afdir,"/8_Module-trait relationships.pdf"),width=7,height=7.5)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


######## Define variable weight containing all column of datTraits

###MM and GS


# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#names of those trait
traitNames=names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


####plot MM vs GS for each trait vs each module


##########example:royalblue and CK
#module="royalblue"
#column = match(module, modNames)
#moduleGenes = moduleColors==module

#trait="CK"
#traitColumn=match(trait,traitNames)

#sizeGrWindow(7, 7)

######

for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      
      #sizeGrWindow(7, 7)
      pdf(file=paste(afdir,"/9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}

#####
names(datExpr0)
probes = names(datExpr0)


#################export GS and MM############### 

geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = paste0(afdir,"/10_GS_and_MM.xls"),sep="\t",row.names=F)



####################################################Visualizing the gene network#######################################################


nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# ls() ###查看工作空间
# rm(list=ls()) ###清空工作空间

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA



#随机选择基因数展示TOM热图
nSelect = 1000
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

# Open a graphical window
#sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA

pdf(file=paste0(afdir,"/13_Network heatmap plot_selected genes.pdf"),width=9, height=9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
dev.off()



####################################################Visualizing the gene network of eigengenes####################################################


#sizeGrWindow(5,7.5)
pdf(file=paste0(afdir,"/14_Eigengene dendrogram and Eigengene adjacency heatmap.pdf"), width=5, height=7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()

#####101种机器学习组合算法####
# 设置工作路径
getwd()
work.path <- "/Users/feifei/Desktop/文章投稿/跑程序性死亡/13.101种预后模型/"
setwd(work.path) 

# 设置其他路径
code.path <- file.path(work.path, "Codes") # 存放脚本
data.path <- file.path(work.path, "InputData") # 存在输入数据（需用户修改）
res.path <- file.path(work.path, "Results") # 存放输出结果
fig.path <- file.path(work.path, "Figures") # 存放输出图片

# 如不存在这些路径则创建路径
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# BiocManager::install("mixOmics")
# BiocManager::install("survcomp")
# devtools::install_github("binderh/CoxBoost")
# install.packages("randomForestSRC")
# install.packages("snowfall")

# 加载需要使用的R包
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

# 加载模型训练以及模型评估的脚本
source(file.path(code.path, "ML.R"))

# 选择最后生成的模型类型：panML代表生成由不同算法构建的模型； multiCox表示抽取其他模型所用到的变量并建立多变量cox模型
FinalModel <- c("panML", "multiCox")[2]

## Training Cohort ---------------------------------------------------------
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，表达谱需有一定变异性，以免建模过程报错）
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 训练集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型）
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 测试集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有StepCox和RunEnet支持输入算法参数
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------
min.selected.var <- 5 # 筛选变量数目的最小阈值
timeVar = "OS.time"; statusVar = "OS" # 定义需要考虑的结局事件，必须出现在Train_surv以及Test_surv中
model <- list() # 初始化模型结果列表
set.seed(seed = 123) # 设置建模种子，使得结果可重复
for (method in methods){ # 循环每一种方法组合
  # method <- "CoxBoost+plsRcox" # [举例]若遇到报错，请勿直接重头运行，可给method赋值为当前报错的算法来debug
  cat(match(method, methods), ":", method, "\n") # 输出当前方法
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
  
  Variable = colnames(Train_expr) # 最后用于构建模型的变量
  for (i in 1:length(method)){
    if (i < length(method)){
      selected.var <- RunML(method = method[i], # 机器学习方法
                            Train_expr = Train_expr, # 训练集有潜在预测价值的变量
                            Train_surv = Train_surv, # 训练集生存数据
                            mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                            timeVar = timeVar, statusVar = statusVar) # 用于训练的生存变量，必须出现在Train_surv中
      if (length(selected.var) > min.selected.var) Variable <- intersect(Variable, selected.var)
    } else {
      model[[method_name]] <- RunML(method = method[i],
                                    Train_expr = Train_expr[, Variable],
                                    Train_surv = Train_surv,
                                    mode = "Model",
                                    timeVar = timeVar, statusVar = statusVar)
    }
  }
  # 如果筛选出的变量小于阈值，则该算法组合无意义，置空（尤其针对以RSF筛选变量的情况，需在ML脚本中尝试调参）
  if (length(selected.var) <= min.selected.var) model[[method_name]] <- NULL
}
saveRDS(model, file.path(res.path, "model.rds")) # 保存所有模型输出

# 当要求最终模型为多变量cox时，对模型进行更新
if (FinalModel == "multiCox"){
  coxmodel <- lapply(model, function(fit){ # 根据各算法最终获得的变量，构建多变量cox模型，从而以cox回归系数和特征表达计算单样本风险得分
    tmp <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                 data = as.data.frame(Train_expr[, ExtractVar(fit)]))
    tmp$subFeature <- fit$subFeature
    return(tmp)
  })
}
saveRDS(coxmodel, file.path(res.path, "coxmodel.rds")) # 保存最终以多变量cox拟合所筛选变量的模型

## Evaluate the model -----------------------------------------------------

## Evaluate the model -----------------------------------------------------
gc()
# 读取已保存的模型列表（请根据需要调整）
model <- readRDS(file.path(res.path, "model.rds")) # 若希望使用各自模型的线性组合函数计算得分，请运行此行
# model <- readRDS(file.path(res.path, "coxmodel.rds")) # 若希望使用多变量cox模型计算得分，请运行此行

methodsValid <- names(model) # 取出有效的模型（变量数目小于阈值的模型视为无效）

# 根据给定表达量计算样本风险评分
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = Test_expr,
                                    type = "lp") # 同原文，使用linear Predictor计算得分
  
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件

# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- model[[method]]$subFeature
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){ data.frame(fit$subFeature) })
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"  # 数据框有两列，包含算法以及算法所筛选出的变量
write.table(fea_df, file.path(res.path, "fea_df.txt"),sep = "\t", row.names = F, col.names = T, quote = F)

# 对各模型计算C-index
Cindexlist <- list()
for (method in methodsValid){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                  Test_expr = Test_expr, # 测试集预后变量，应当包含训练集中所有的变量，否则会报错
                                  Test_surv = Test_surv, # 训练集生存数据，应当包含训练集中所有的变量，否则会报错
                                  Train_expr = Train_expr, # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                  Train_surv = Train_surv, # 若需要同时评估训练集，则给出训练集生存数据，否则置NULL
                                  Train_name = "TCGA", # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                  # Train_expr = NULL,
                                  # Train_surv = NULL, 
                                  cohortVar = "batch", # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                  # timeVar = OS.time, # 用于评估的生存时间，必须出现在Test_surv中；这里是OS.time
                                  # statusVar = statusVar
  ) # 用于评估的生存状态，必须出现在Test_surv中；这里是OS
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

Cindex_mat <- read.table(file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] # 最优模型（即测试集[或者训练集+测试集]C指数均值最大）所筛选的特征

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置绘图时的队列颜色
names(CohortCol) <- colnames(Cindex_mat)

# 调用简易绘图函数
# 极浅蓝
# 浅蓝

hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "#81C784", # 右侧柱状图颜色
                    col = c("#E1F5FE",  "#B3E5FC" ,"#4FC3F7"), # 热图颜色
                    cellwidth = 1, cellheight = 0.5, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(fig.path, "heatmap of cindex1.pdf"), width = 2 * ncol(Cindex_mat) + 3, height = 0.5 * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())


####保存其他文件，选择性跑####
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = Train_expr,
                                    type = "lp") # 同原文，使用linear Predictor计算得分
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, "Train_exprmat.txt",sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件


RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = Test_set,
                                    type = "lp") # 同原文，使用linear Predictor计算得分
  
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, "RS_mat.txt",sep = "\t", row.names = T, col.names = NA, quote = F)

# 确保加载了必要的包
library(ComplexHeatmap)
# circlize包也需要被加载，因为它包含colorRamp2函数
library(circlize)

# 创建一个颜色函数
col_func <- colorRamp2(c(min(Cindex_mat), 0, max(Cindex_mat)), c("#1CB8B2", "#FFFFFF", "#EEB849"))

# 调用SimpleHeatmap函数
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat,
                    avg_Cindex = avg_Cindex,
                    CohortCol = CohortCol,
                    barCol = "steelblue",
                    col = col_func, # 使用颜色函数而不是颜色向量
                    cellwidth = 1, cellheight = 0.5,
                    cluster_columns = FALSE, cluster_rows = FALSE)

####肿瘤的分型我们参考了以下官方教程####


####加载包####
library("MOVICS")

#####加载数据###
# 加载乳腺癌示例数据
load(system.file("extdata", "brca.tcga.RData", package = "MOVICS", mustWork = TRUE))
load(system.file("extdata", "brca.yau.RData", package = "MOVICS", mustWork = TRUE))
# 示例数据的打印名称
names(brca.tcga)
names(brca.yau)

mo.data <- brca.tcga[1:4]
count <- brca.tcga$count
fpkm <- brca.tcga$fpkm
maf <- brca.tcga$maf
segment <- brca.tcga$segment
surv.info <- brca.tcga$clin.info
#### scenario 1:####
# 考虑到我们正在处理具有 2 行且具有 NA 值的表达式数据
tmp <- brca.tcga$mRNA.expr # 获取表达式数据
dim(tmp) # 检查数据维度
tmp[1,1] <- tmp[2,2] <- NA # 使用 NA 值设置 2 行
tmp[1:3,1:3] # 检查数据
elite.tmp <- getElites(dat = tmp,
                       method = "mad",
                       na.action = "rm", # NA 值将被删除
                       elite.pct = 1) # 将选择所有 （100%） 要素
#> --2 具有 NA 值的要素将被移除。
#> 缺少 elite.num 然后使用 elite.pct
dim(elite.tmp$elite.dat) # 我们删除了 2 行包含 NA 数据
elite.tmp <- getElites(dat = tmp,
                       method = "mad",
                       na.action = "impute", # 将估算NA值
                       elite.pct = 1)

#> 缺少 elite.num 然后使用 elite.pct
dim(elite.tmp$elite.dat) # 保留所有数据
elite.tmp$elite.dat[1:3,1:3] #NA 值已估算

# scenario 2:
# 考虑到我们正在处理连续数据并使用MAD或SD来选择精英
tmp <- brca.tcga$mRNA.expr # 获取包含 500 个特征的表达式数据
elite.tmp <- getElites(dat = tmp,
                       method = "mad",
                       elite.pct = 0.1) # 保留前 10% 的功能
dim(elite.tmp$elite.dat) # 剩余 50 名精英

elite.tmp <- getElites(dat = tmp,
                       method = "sd",
                       elite.num = 100, # 保留前 100 个功能
                       elite.pct = 0.1) # 由于指示的精英而禁用。

dim(elite.tmp$elite.dat) # 剩下 100 名精英

# scenario 3:
# 考虑到我们正在处理连续数据并使用PCA来选择精英
tmp <- brca.tcga$mRNA.expr # 获取包含 500 个特征的表达式数据
elite.tmp <- getElites(dat = tmp,
                       method = "pca",
                       pca.ratio = 0.95) # 选择pcs比例
dim(elite.tmp$elite.dat) # 剩余 204 名精英 （PC）

# scenario 4:
# 考虑到我们正在处理数据并使用 Cox 选择精英
tmp <- brca.tcga$mRNA.expr # 获取表达式数据
elite.tmp <- getElites(dat = tmp,
                       method = "cox",
                       surv.info = surv.info, # 获取表达式数据
                       p.cutoff = 0.05,
                       elite.num = 100) # 禁用，因为 Cox 引用 P.cutoff
dim(elite.tmp$elite.dat) # 获得 125 名精英

table(elite.tmp$unicox$pvalue < 0.05) # 125个基因在Cox中具有pvalue<0.05

tmp <- brca.tcga$mut.status # 获取突变数据
elite.tmp <- getElites(dat = tmp,
                       method = "cox",
                       surv.info = surv.info,
                       p.cutoff = 0.05,
                       elite.num = 100)

dim(elite.tmp$elite.dat) # get 3 elites
table(elite.tmp$unicox$pvalue < 0.05) # 3 mutations have pvalue<0.05

tmp <- brca.tcga$mut.status # get mutation data
rowSums(tmp) # 检查突变频率
elite.tmp <- getElites(dat = tmp,
                       method = "freq", # 必须设置为“频率”
                       elite.num = 80, # 现在参考突变频率
                       elite.pct = 0.1) # 由于指示的精英而禁用。
rowSums(elite.tmp$elite.dat) # 保留在 >80 个样本中突变的基因
elite.tmp <- getElites(dat = tmp,
                       method = "freq",
                       elite.pct = 0.2) # 现在参考突变/样本量
rowSums(elite.tmp$elite.dat) # 保留>0.2*643=128.6样本突变的基因

####GET Module####

# 确定最佳聚类数（可能需要一段时间）
optk.brca <- getClustNum(data = mo.data,
                         is.binary = c(F,F,F,T), # 第 4 个数据是二进制矩阵
                         try.N.clust = 2:8, # 尝试从 2 到 8 分组
                         fig.name = "CLUSTER NUMBER OF TCGA-BRCA")

# 执行 iClusterBayes（可能需要一段时间）
iClusterBayes.res <-
  getiClusterBayes(data = mo.data,
                   N.clust = 5,
                   type = c("gaussian","gaussian","gaussian","binomial"),
                   n.burnin = 1800,
                   n.draw = 1200,
                   prior.gamma = c(0.5, 0.5, 0.5, 0.5),
                   sdev = 0.05,
                   thin = 3)

iClusterBayes.res <-
  getMOIC(data = mo.data,
          N.clust = 5,
          methodslist = "iClusterBayes", # 在此处仅指定一种算法
          type = c("gaussian","gaussian","gaussian","binomial"), # 数据类型
          n.burnin = 1800,
          n.draw = 1200,
          prior.gamma = c(0.5, 0.5, 0.5, 0.5),
          sdev = 0.05,
          thin = 3)

# 使用其余 9 种算法执行多组学综合聚类
moic.res.list <-
  getMOIC(data = mo.data,
          methodslist = list("SNF", "PINSPlus", "NEMO",
                             "COCA", "LRAcluster", "ConsensusClustering",
                             "IntNMF", "CIMLR", "MoCluster"),
          N.clust = 5,
          type = c("gaussian", "gaussian", "gaussian", "binomial"))
moic.res.list <- append(moic.res.list,
                        list("iClusterBayes" = iClusterBayes.res))

# 将 moic.res.list 保存到本地路径
save(moic.res.list, file = "moic.res.list.rda")
load("moic.res.list.rda")
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name = "CONSENSUS HEATMAP",
                               distance = "euclidean",
                               linkage = "average")
getSilhouette(sil = cmoic.brca$sil, # 一个由 getConsensusMOIC（） 返回的 sil 对象
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height = 5.5,
              width = 5)

# 将 beta 值转换为 M 值以获得更强的信号
indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# 热图的数据规范化
plotdata <- getStdiz(data = indata,
                     halfwidth = c(2,2,2,NA), # 无突变截断
                     centerFlag = c(T,T,T,F), # 无突变中心
                     scaleFlag = c(T,T,T,F)) # 无突变量表
feat <- iClusterBayes.res$feat.res
feat1 <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"]
feat2 <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3 <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat4 <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)

#为每个组学数据设置颜色
mRNA.col <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white" , "#FF3C38")
meth.col <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col <- c("grey90" , "black")
col.list <- list(mRNA.col, lncRNA.col, meth.col, mut.col)
# 综合热图（可能需要一段时间）
getMoHeatmap(data = plotdata,
             row.title = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary = c(F,F,F,T), # 第 4 个数据是二进制矩阵
             legend.name = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res = iClusterBayes.res$clust.res, # 群集结果
             clust.dend = NULL, # 无树状图
             show.rownames = c(F,F,F,F), # 为每个组学数据指定
             show.colnames = FALSE, # 不显示示例名称
             annRow = annRow, # 标记所选要素
             color = col.list,
             annCol = NULL, # 没有示例注释
             annColors = NULL, # 无注释颜色
             width = 10, # 每个子热图的宽度
             height = 5, # 每个子热图的高度
             fig.name = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")
getMoHeatmap(data = plotdata,
             row.title = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary = c(F,F,F,T),
             legend.name = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res = moic.res.list$COCA$clust.res, # 群集结果
             clust.dend = moic.res.list$COCA$clust.dend, # 显示样本树状图
             color = col.list,
             width = 10,
             height = 5,
             fig.name = "COMPREHENSIVE HEATMAP OF COCA")
# 提取PAM50，病理分期和年龄进行样本注释
annCol <- surv.info[,c("PAM50", "pstage", "age"), drop = FALSE]
# 为示例注释生成相应的颜色
annColors <- list(age =
                    circlize::colorRamp2(breaks
                                         = c(min(annCol$age),
                                             median(annCol$age),
                                             max(annCol$age)),
                                         colors
                                         = c("#0000AA", "#555555", "#AAAA00")),
                  PAM50 = c("Basal" = "blue",
                            "Her2" = "red",
                            "LumA" = "yellow",
                            "LumB" = "green",
                            "Normal" = "black"),
                  pstage = c("T1" = "green",
                             "T2" = "blue",
                             "T3" = "red",
                             "T4" = "yellow",
                             "TX" = "black"))
# 综合热图
getMoHeatmap(data = plotdata,
             row.title = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary = c(F,F,F,T),
             legend.name = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res = cmoic.brca$clust.res, # 共识MOIC结果
             clust.dend = NULL,
             show.rownames = c(F,F,F,F),
             show.colnames = FALSE,
             show.row.dend = c(F,F,F,F), # 不显示要素的树状图
             annRow = NULL, # 未选择要素
             color = col.list,
             annCol = annCol, # 示例注释
             annColors = annColors, #批注颜色
             width = 10,
             height = 5,
             fig.name = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")
####COMP Module ####
# 生存比较
surv.brca <-
  compSurv(moic.res = cmoic.brca,
           surv.info = surv.info,
           convt.time = "m", # 将日单位转换为月
           surv.median.line = "h", # 在中位生存期绘制水平线
           xyrs.est = c(5,10), # 估计 5 年和 10 年生存率
           fig.name = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")

print(surv.brca)


clin.brca <-
  compClinvar(moic.res = cmoic.brca,
              var2comp = surv.info, # 数据帧需要总结
              strata = NULL, # 分层变量（NULL 表示“亚型”）
              factorVars = c("PAM50","pstage","fustat"), # 分类变量
              nonnormalVars = "futime", # 使用非参数检验的特征
              exactVars = "pstage", # 使用精确测试的功能
              doWord = TRUE, # 在本地路径中生成.docx文件
              tab.name = "SUMMARIZATION OF CLINICAL FEATURES")
print(clin.brca$compTab)


# 突变频率比较
mut.brca <-
  compMut(moic.res = cmoic.brca,
          mut.matrix = brca.tcga$mut.status, # 二元体细胞突变基质
          doWord = TRUE, # 生成.docx格式的表格
          doPlot = TRUE, # 绘制OncoPrint
          freq.cutoff = 0.05, # 保留至少5%样本中突变的基因
          p.adj.cutoff = 0.05, # 在OncoPrint中保留具有调整p值<0.05的基因
          innerclust = TRUE, # 在每个子类型内执行聚类
          annCol = annCol, # 热图的相同注释
          annColors = annColors, # 热图的注释颜色相同
          width = 6,
          height = 2,
          fig.name = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
          tab.name = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
print(mut.brca) 



head(maf)
# 比较 TMB
tmb.brca <- compTMB(moic.res = cmoic.brca,
                    maf = maf,
                    rmDup = TRUE, # 删除每个样本的重复变体
                    rmFLAGS = FALSE, #保持标志突变
                    exome.size = 38, # 估计外显子组大小
                    test.method = "nonparametric", # 统计测试方法
                    fig.name = "DISTRIBUTION OF TMB AND TITV")
head(tmb.brca$TMB.dat)

# 更改区段数据的列名称
colnames(segment) <- c("sample","chrom","start","end","value")
head(segment)

# 比较 FGA、FGG 和 FGL(拷贝数变异)
fga.brca <-
  compFGA(moic.res = cmoic.brca,
          segment = segment,
          iscopynumber = FALSE, # 这是一个分段拷贝编号文件
          cnathreshold = 0.2, # 确定 CNA 收益或损失的阈值
          test.method = "nonparametric", # 统计测试方法
          fig.name = "BARPLOT OF FGA")
head(fga.brca$summary)


# 药物敏感性比较
drug.brca <-compDrugsen(moic.res = cmoic.brca,
                        norm.expr = fpkm[,cmoic.brca$clust.res$samID],
                        drugs = c("Cisplatin", "Paclitaxel"), # 药物名称的载体
                        tissueType = "breast", # 选择特定的组织类型
                        test.method = "nonparametric", # 统计测试方法
                        prefix = "BOXVIOLIN OF ESTIMATED IC50")
head(drug.brca$Cisplatin)


# 自定义 PSTAGE 的因子水平
surv.info$pstage <- factor(surv.info$pstage, levels = c("TX","T1","T2","T3","T4"))

# 协议比较（最多支持 6 种分类，包括当前子类型）
agree.brca <-
  compAgree(moic.res = cmoic.brca,
            subt2comp = surv.info[,c("PAM50","pstage")],
            doPlot = TRUE,
            box.width = 0.2,
            fig.name = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")
print(agree.brca)
# 使用 edgeR 运行 DEA
runDEA(dea.method = "edger",
       expr = count, # 原始计数数据
       moic.res = cmoic.brca,
       prefix = "TCGA-BRCA") #图形名称前缀

# 使用 DESeq2 运行 DEA
runDEA(dea.method = "deseq2",
       expr = count,
       moic.res = cmoic.brca,
       prefix = "TCGA-BRCA")
# 与limma一起运行DEA
runDEA(dea.method = "limma",
       expr = fpkm, # 标准化表达数据
       moic.res = cmoic.brca,
       prefix = "TCGA-BRCA")

# 选择 edgeR 结果以鉴定亚型特异性上调生物标志物
marker.up <-
  runMarker(moic.res = cmoic.brca,
            dea.method = "edger", # DEA 方法的名称
            prefix = "TCGA-BRCA", # 必须与 runDEA（） 中的参数相同
            dat.path = getwd(), # DEA 文件的路径
            res.path = getwd(), #保存标记文件的路径
            p.cutoff = 0.05, # p 截止值以识别显著 DEG
            p.adj.cutoff = 0.05, # 用于识别显著 DEG 的 PAJ 截止值
            dirct = "up", # 表达失调的方向
            n.marker = 100, # 每种亚型的生物标志物数量
            doplot = TRUE, # 生成对角线热图
            norm.expr = fpkm, # 使用规范化表达式作为热图输入
            annCol = annCol, # 热图中的示例注释
            annColors = annColors, # 示例注释的颜色
            show_rownames = FALSE, # 不显示行名（生物标志物名称）
            fig.name = "UPREGULATED BIOMARKER HEATMAP")

# 检查上调的生物标志物
head(marker.up$templates)
# 选择 LIMMA 结果以鉴定亚型特异性下调生物标志物
marker.dn <-
  runMarker(moic.res = cmoic.brca,
            dea.method = "limma",
            prefix = "TCGA-BRCA",
            dirct = "down",
            n.marker = 50, # switch to 50
            doplot = TRUE,
            norm.expr = fpkm,
            annCol = annCol,
            annColors = annColors,
            fig.name = "DOWNREGULATED BIOMARKER HEATMAP")

# 必须找到 msigdb 文件的绝对路径
MSIGDB.FILE <-
  system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

# 运行 GSEA 以使用 edgeR 的结果识别上调的 GO 通路
gsea.up <-
  runGSEA(moic.res = cmoic.brca,
          dea.method = "edger", # DEA 方法的名称
          prefix = "TCGA-BRCA", # 必须与 runDEA（） 中的参数相同
          dat.path = getwd(), # DEA 文件的路径
          res.path = getwd(), # 保存 GSEA 文件的路径
          msigdb.path = MSIGDB.FILE, # 必须是 msigdb 文件的绝对路径
          norm.expr = fpkm, # 使用规范化表达式计算扩充分数
          dirct = "up", # 通路失调的方向
          p.cutoff = 0.05, # p 截止以确定重要途径
          p.adj.cutoff = 0.25, # PADJ截止以确定重要途径
          gsva.method = "gsva", # 单样品富集分数的计算方法
          norm.method = "mean", # 计算亚型特异性富集分数的方法
          fig.name = "UPREGULATED PATHWAY HEATMAP")
print(gsea.up$gsea.list$CS1[1:6,3:6])
head(round(gsea.up$grouped.es,3))

# 运行 GSEA 以使用 DESeq2 的结果识别下调的 GO 途径
gsea.dn <- runGSEA(moic.res = cmoic.brca,
                   dea.method = "deseq2",
                   prefix = "TCGA-BRCA",
                   msigdb.path = MSIGDB.FILE,
                   norm.expr = fpkm,
                   dirct = "down",
                   p.cutoff = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method = "ssgsea", # 切换到SSGSEA
                   norm.method = "median", # 切换到中位数
                   fig.name = "DOWNREGULATED PATHWAY HEATMAP")

# 必须找到基因集文件的绝对路径
GSET.FILE <-
  system.file("extdata", "gene sets of interest.gmt", package = "MOVICS", mustWork = TRUE)

# 运行 GSVA 以根据给定的目标基因集估计单样本富集分数
gsva.res <-
  runGSVA(moic.res = cmoic.brca,
          norm.expr = fpkm,
          gset.gmt.path = GSET.FILE, # 基因集文件的绝对路径
          gsva.method = "gsva", # 单样品富集分数的计算方法
          annCol = annCol,
          annColors = annColors,
          fig.path = getwd(),
          fig.name = "GENE SETS OF INTEREST HEATMAP",
          height = 5,
          width = 8)

# 检查原始扩充分数
print(gsva.res$raw.es[1:3,1:3])

#检查 z 评分和截断扩充分数
print(gsva.res$scaled.es[1:3,1:3])

#使用上调的生物标志物在丘队列中运行NTP
yau.ntp.pred <-
  runNTP(expr = brca.yau$mRNA.expr,
         templates = marker.up$templates, # 模板已在 runMarker（） 中准备
         scaleFlag = TRUE, # 缩放输入数据
         centerFlag = TRUE, # 居中输入数据
         doPlot = TRUE, # 生成热图
         fig.name = "NTP HEATMAP FOR YAU")
head(yau.ntp.pred$ntp.res)

# 比较丘成桐队列的生存结局
surv.yau <-
  compSurv(moic.res = yau.ntp.pred,
           surv.info = brca.yau$clin.info,
           convt.time = "m", # 切换到月份
           surv.median.line = "hv", # 切换到两者
           fig.name = "KAPLAN-MEIER CURVE OF NTP FOR YAU")
print(surv.yau)

# 比较丘成桐队列中的一致性
agree.yau <-
  compAgree(moic.res = yau.ntp.pred,
            subt2comp = brca.yau$clin.info[, "PAM50", drop = FALSE],
            doPlot = TRUE,
            fig.name = "YAU PREDICTEDMOIC WITH PAM50")
yau.pam.pred <- runPAM(train.expr = fpkm,
                       moic.res = cmoic.brca,
                       test.expr = brca.yau$mRNA.expr)
print(yau.pam.pred$IGP)

# 使用 NTP 预测发现队列中的亚型
tcga.ntp.pred <- runNTP(expr = fpkm,
                        templates = marker.up$templates,
                        doPlot = FALSE)

#使用 PAM 预测发现队列中的子类型
tcga.pam.pred <- runPAM(train.expr = fpkm,
                        moic.res = cmoic.brca,
                        test.expr = fpkm)

# 检查发现TCGA-BRCA中当前和NTP预测亚型之间的一致性
runKappa(subt1 = cmoic.brca$clust.res$clust,
         subt2 = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")
# 检查发现TCGA-BRCA中当前和PAM预测的亚型之间的一致性
runKappa(subt1 = cmoic.brca$clust.res$clust,
         subt2 = tcga.pam.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM",
         fig.name = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and PAM")
# 在验证Yau-BRCA中检查NTP和PAM预测子类型之间的一致性
runKappa(subt1 = yau.ntp.pred$clust.res$clust,
         subt2 = yau.pam.pred$clust.res$clust,
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.name = "CONSISTENCY HEATMAP FOR YAU")
# 原始临床信息列表为“Clust.res”和“Mo.method”
pseudo.moic.res <- list("clust.res" = surv.info,
                        "mo.method" = "PAM50")
#制作伪 samID
pseudo.moic.res$clust.res$samID <- rownames(pseudo.moic.res$clust.res)

#使用映射关系创建伪克吕斯特（make pseudo clust using a mapping relationship）
pseudo.moic.res$clust.res$clust <- sapply(pseudo.moic.res$clust.res$PAM50,
                                          switch,
                                          "Basal" = 1, # relabel Basal as 1
                                          "Her2" = 2, # relabel Her2 as 2
                                          "LumA" = 3, # relabel LumA as 3
                                          "LumB" = 4, # relabel LumnB as 4
                                          "Normal" = 5) # relabel Normal as 5
head(pseudo.moic.res$clust.res)

# 生存比较
pam50.brca <-
  compSurv(moic.res = pseudo.moic.res,
           surv.info = surv.info,
           convt.time = "y", # 将日单位转换为年
           surv.median.line = "h", # 在中位生存期绘制水平线
           fig.name = "KAPLAN-MEIER CURVE OF PAM50 BY PSEUDO")