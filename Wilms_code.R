
setwd("X:\\Wilms")
load("Wilms.RData")
#rm(list = ls())

tumor.path <- "X:\\Wilms";setwd(tumor.path) #create dir
data.path   <- file.path(tumor.path, "InputData")
fig1.path    <- file.path(tumor.path, "Figure1")
fig2.path    <- file.path(tumor.path, "Figure2")
fig3.path    <- file.path(tumor.path, "Figure3")
fig4.path    <- file.path(tumor.path, "Figure4")
fig5.path    <- file.path(tumor.path, "Figure5")
fig6.path    <- file.path(tumor.path, "Figure6")
fig7.path    <- file.path(tumor.path, "Figure7")
Scripts    <- file.path(tumor.path, "Scripts")

library("tidyverse")
if (!file.exists(tumor.path)) { dir.create(tumor.path) }
if (!file.exists(data.path)) { dir.create(data.path) }
if (!file.exists(fig1.path)) { dir.create(fig1.path) }
if (!file.exists(fig2.path)) { dir.create(fig2.path) }
if (!file.exists(fig3.path)) { dir.create(fig3.path) }
if (!file.exists(fig4.path)) { dir.create(fig4.path) }
if (!file.exists(fig5.path)) { dir.create(fig5.path) }
if (!file.exists(fig6.path)) { dir.create(fig6.path) }
if (!file.exists(fig7.path)) { dir.create(fig7.path) }

# load R package
library(sva)
library(ConsensusClusterPlus)
library(pheatmap)
library(corrplot)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Boruta)
library(org.Hs.eg.db)
library(enrichplot)
library(ggalluvial)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(forcats)
library(maftools)
library(GSVA)
#library(ComplexHeatmap)
library(gplots)
library(clusterProfiler)
library(tidyr)
library(ggplot2)
library(estimate)
library(ggpubr)
library('progress')
source(file.path(Scripts,"twoclasslimma.R"))
# load R package
library(NMF)
library(survival)
library(survminer)
library(sva)
library(Rtsne)
library(ComplexHeatmap)
library(gplots)
source(file.path(Scripts,"batchPCA.R"))
# set color


# 主调色板（适合 5 组线条区分）
jco <- c("#E76F51", "#2A9D8F", "#E9C46A", "#264653", "#8E24AA")  
# 橙红-蓝绿-沙黄-深蓝-紫罗兰

# 三色调板1（对比鲜明）
jco2 <- c("#E63946", "#3399CC", "#EF8A09")  
# 鲜红-亮黄-蓝灰

# 三色调板2（相对柔和但依然区分明显）
jco3 <- c("#F4A261", "#2A9D8F", "#9B5DE5","blue")  
# 橙-蓝绿-薰衣草紫

# 三色调板3（深色背景上适用）
jco4 <- c("#1D3557", "#E63946", "#06D6A0")  
# 深蓝-玫红-薄荷绿

# 三色调板4（灰调中性色）
jco5 <- c("#3A86FF", "#FFBE0B", "#8338EC")  
# 蓝-亮黄-紫

# yjp（温和渐变型，适合浅背景线条）
yjp <- c("#F4A6B9", "#FFB347", "#FF6F61", "#2A9D8F", "#118AB2")  
# 粉-橘-红-蓝绿-蔚蓝
yjp2 <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")


tt<-jco4
pie(rep(1, length(tt)), col=tt)

#load R package

########## Figure 1 combat

#####GEO去批次

####combat of GEO datasets######
TARGET.expr<-read.table(file.path(data.path,"TARGET.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
TARGET.expr<-round(log2(TARGET.expr+1),2)
TARGET.expr<-TARGET.expr[!apply(TARGET.expr,1,function(x) length(x[x<1]))>(0.9*ncol(TARGET.expr)),]
TARGET.clin<-read.table(file.path(data.path,"TARGET.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
merge<-intersect(colnames(TARGET.expr),rownames(TARGET.clin))
TARGET.expr<-TARGET.expr[,merge]
TARGET.clin<-TARGET.clin[merge,]
TARGET.clin$OS.time<-TARGET.clin$OS.time/30.5
TARGET.clin$Age<-TARGET.clin$Age/30.5
range(TARGET.expr)

GSE60850.expr<-read.csv(file.path(data.path,"GSE60850.expr.csv"),header = T,row.names = 1,check.names = F,stringsAsFactors = F)
GSE60850.expr<-GSE60850.expr %>% 
  drop_na()
GSE60850.clin<-read.table(file = file.path(data.path,"GSE60850.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
merge<-intersect(colnames(GSE60850.expr),rownames(GSE60850.clin))
GSE60850.expr<-GSE60850.expr[,merge]
GSE60850.clin<-GSE60850.clin[merge,]


range(GSE60850.expr)

exp <- GSE60850.expr
exp=exp-apply(exp,1,min)
exp=normalizeBetweenArrays(exp)#校正
range(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)#画图


GSE60850.expr<-exp


library(sva)
library(cluster)
library(oompaBase)

# 
comgene <- intersect(rownames(TARGET.expr),  rownames(GSE60850.expr))
# 
combined.expr <- cbind.data.frame(TARGET.expr[comgene,],
                                  GSE60850.expr[comgene,])
# 
batchPCA(indata = t(scale(t(combined.expr))),
         batch = rep(c("TARGET","GSE60850"), times = c(ncol(TARGET.expr),ncol(GSE60850.expr))),
         fig.dir = fig1.path,
         PCA.fig.title = "Raw PCA for combined expression profile",
         cols = jco4[1:2],
         showID = F,
         cex = 0.7,
         showLegend = T) # 

range(combined.expr)

# 
batch <- data.frame(batch = rep(c("TARGET","GSE60850"), c(ncol(TARGET.expr),ncol(GSE60850.expr))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))


#
batchPCA(indata = t(scale(t(combined.expr.combat))),
         batch = rep(c("TARGET","GSE60850"), times = c(ncol(TARGET.expr),ncol(GSE60850.expr))),
         fig.dir = fig1.path,
         PCA.fig.title = "Combat PCA for combined expression profile",
         cols = jco4[1:2],
         showID = F,
         cex = 0.7,
         showLegend = T) #


TARGET.clin$cohort<-"TARGET"
GSE60850.clin$cohort<-"GSE60850"

GSE60850.expr<-combined.expr.combat[,rownames(GSE60850.clin)]
TARGET.expr<-combined.expr.combat[,rownames(TARGET.clin)]

#########find the key pathways

library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#START GSVA
load(file.path(data.path, "hallmark.gs.RData"))


gsym.expr <- TARGET.expr
head(gsym.expr)

# 这一句就完成了GSVA分析
gsva_es <- gsva(gsvaParam(as.matrix(gsym.expr), gs))

head(gsva_es)

write.csv(gsva_es, file = file.path(data.path,"gsva_output.csv"), quote = F)

group<-TARGET.clin
group<- group[order(group$OS, decreasing = F), ]
table(group$OS)

group_list <- data.frame(sample = rownames(group), group = c(rep("a", 74), rep("b", 50)))
head(group_list)

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design

# 构建差异比较矩阵
contrast.matrix <- makeContrasts(b-a, levels = design)

# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

write.csv(x, file = file.path(data.path,"gsva_limma.csv"), quote = F)

#输出t值，用做FigureYa39bar的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
write.csv(df, file = file.path(data.path,"easy_input2_for39bar.csv"), quote = F, row.names = F)

#开始画图
df <- read.csv(file = file.path(data.path,"easy_input2_for39bar.csv"))
head(df)

cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

# 添加一个新列 label_group，控制文字颜色
sortdf$label_group <- ifelse(sortdf$group == 2, "ns", "sig")

ggplot(sortdf, aes(x = ID, y = score)) +
  #coord_flip() +  # 横图
  
  # 棒棒糖“棒”部分（线），按 group 映射颜色
  geom_segment(aes(x = ID, xend = ID, y = 0, yend = score, color = group), size = 2) +
  
  # 棒棒糖“糖果”部分（点），按 group 映射颜色
  geom_point(aes(color = group), shape = 19, size = 5) +
  
  # 两条显著性阈值虚线
  geom_hline(yintercept = c(-cutoff, cutoff),
             color = "grey60", linetype = 2, size = 0.3) +
  
  # 显著负向标签
  geom_text(data = subset(sortdf, score < 0),
            aes(x=ID, y= 0.1, label= paste0(" ", ID), color = label_group),
            size = 3, hjust = "inward") +
  
  # 显著正向标签
  geom_text(data = subset(sortdf, score > 0),
            aes(x=ID, y= -0.1, label= paste0(" ", ID), color = label_group),
            size = 3, hjust = "outward") +
  
  
  # 设置线条/点的颜色
  scale_color_manual(
    name = "2",
    values = c(
      "1"="#00B072",  # green
      "2"="#BBBBBB",  # gray
      "3"="#FDB530",  # orange
      "sig" = "black",
      "ns" = "snow3"
    ),
    guide = FALSE
  ) +
  
  xlab("") +
  ylab("t value of GSVA score, Dead vs. Alive") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 0.6),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

ggsave(file.path(fig1.path,"gsva.pdf"), width = 6, height = 8)

##############横着

ggplot(sortdf, aes(x = ID, y = score)) +
  
  # 棒棒糖“棒”部分（线）
  geom_segment(aes(x = ID, xend = ID, y = 0, yend = score, color = group), size = 2) +
  
  # 棒棒糖“糖果”部分（点）
  geom_point(aes(color = group), shape = 19, size = 5) +
  
  # 显著性阈值线
  geom_hline(yintercept = c(-cutoff, cutoff), color = "grey60", linetype = 2, size = 0.3) +
  
  # 添加标签：统一放在点上方或下方并旋转 45 度
  
  # 显著负向标签
  geom_text(data = subset(sortdf, score < 0),
            aes(x=ID, y= 0.1, label= paste0(" ", ID), color = label_group),
            size = 2.2,angle = 60, hjust = "outward") +
  
  # 显著正向标签
  geom_text(data = subset(sortdf, score > 0),
            aes(x=ID, y= -0.1, label= paste0(" ", ID), color = label_group),
            size = 2.2, angle = 60,hjust = "inward") +
  
  # 手动设置颜色
  scale_color_manual(
    values = c(
      "1" = "#00B072",   # green
      "2" = "#BBBBBB",   # gray
      "3" = "#FDB530",   # orange
      "sig" = "black",
      "ns" = "snow3"
    ),
    guide = FALSE
  ) +
  
  ylab("") +
  xlab("t value of GSVA score, Dead vs. Alive") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.6),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),axis.text.x = element_blank()
  )
ggsave(file.path(fig1.path,"gsva2.pdf"), width = 8, height = 3.8)

##############GSE60850横着
gsym.expr <- GSE60850.expr
head(gsym.expr)

# 这一句就完成了GSVA分析
gsva_es <- gsva(gsvaParam(as.matrix(gsym.expr), gs))

head(gsva_es)

write.csv(gsva_es, file = file.path(data.path,"gsva_output.csv"), quote = F)

group<-GSE60850.clin
group<- group[order(group$OS, decreasing = F), ]
table(group$OS)

group_list <- data.frame(sample = rownames(group), group = c(rep("a", 23), rep("b", 12)))
head(group_list)

design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design

# 构建差异比较矩阵
contrast.matrix <- makeContrasts(b-a, levels = design)

# 差异分析，b vs. a
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

write.csv(x, file = file.path(data.path,"gsva_limma.csv"), quote = F)

#输出t值，用做FigureYa39bar的输入数据
pathway <- str_replace(row.names(x), "HALLMARK_", "")
df <- data.frame(ID = pathway, score = x$t)
write.csv(df, file = file.path(data.path,"easy_input2_for39bar.csv"), quote = F, row.names = F)

#开始画图
df <- read.csv(file = file.path(data.path,"easy_input2_for39bar.csv"))
head(df)

cutoff <- 1
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)
sortdf$label_group <- ifelse(sortdf$group == 2, "ns", "sig")

ggplot(sortdf, aes(x = ID, y = score)) +
  
  # 棒棒糖“棒”部分（线）
  geom_segment(aes(x = ID, xend = ID, y = 0, yend = score, color = group), size = 2) +
  
  # 棒棒糖“糖果”部分（点）
  geom_point(aes(color = group), shape = 19, size = 5) +
  
  # 显著性阈值线
  geom_hline(yintercept = c(-cutoff, cutoff), color = "grey60", linetype = 2, size = 0.3) +
  
  # 添加标签：统一放在点上方或下方并旋转 45 度
  
  # 显著负向标签
  geom_text(data = subset(sortdf, score < 0),
            aes(x=ID, y= 0.1, label= paste0(" ", ID), color = label_group),
            size = 2.2,angle = 60, hjust = "outward") +
  
  # 显著正向标签
  geom_text(data = subset(sortdf, score > 0),
            aes(x=ID, y= -0.1, label= paste0(" ", ID), color = label_group),
            size = 2.2, angle = 60,hjust = "inward") +
  
  # 手动设置颜色
  scale_color_manual(
    values = c(
      "1" = "#57A7B9",   # green
      "2" = "#BBBBBB",   # gray
      "3" = "#CB0066",   # orange
      "sig" = "black",
      "ns" = "snow3"
    ),
    guide = FALSE
  ) +
  
  ylab("") +
  xlab("t value of GSVA score, Dead vs. Alive") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(size = 0.6),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),axis.text.x = element_blank()
  )
ggsave(file.path(fig1.path,"gsva GSE60850.pdf"), width = 8, height = 3.8)


###########select cellcycle gene
cellcycle<-read.table(file = file.path(data.path,"pathways.gmt"),header=T,sep="\t",row.names=1,check.names=F)
cellcycle<-cellcycle[,-c(1:1)]

test<-t(cellcycle)

tt<-matrix(test,dimnames=list(t(outer(colnames(test),rownames(test),FUN=paste)),NULL))

tt[tt==""]<-NA
tt<-as.data.frame(tt)
tt<-tt %>% 
  drop_na()
select<-tt$V1[!duplicated(tt$V1)]


#################UpsetVenn plot 筛选可用的Histon基因#################
library("ggvenn")
# Default plot

TARGET.expr2<-read.table(file.path(data.path,"TARGET.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
TARGET.expr2<-TARGET.expr2[!apply(TARGET.expr2,1,function(x) length(x[x<1]))>(0.9*ncol(TARGET.expr2)),]
TARGET.expr2<-log2(TARGET.expr2+1)
range(TARGET.expr2)

x <- list(
  TCGA.TARGET=rownames(TARGET.expr2),  
  GSE60850=rownames(read.table(file = file.path(data.path,"GSE60850.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)),
  cellcycle_gene=select
)

ggvenn(x)

library(UpSetR)
upset(fromList(x),order.by = "freq",nsets = 4,
      queries = list(list(query = intersects, params = list("TCGA.TARGET","GSE60850","cellcycle_gene"), active = T)))

dd2<-fromList(x)

#devtools::install_github("PhDMeiwp/basicPackages@master", force = TRUE) 安装下面的包时使用
library(basicPackages)
#basicPackages::install.yyplot()安装下面的包时使用

library(yyplot)
p2 <- yyplot::ggvenn(dd2,fill_color = c("#4dbbd599", "#e8924299", "#e762d799"))+
  theme_void()+
  theme(legend.position = "right")
p2

library(ggplotify) #把别的图转为ggplot2
library(ggimage) # 组合图片
p1<-upset(fromList(x),order.by = "freq",nsets = 4,
          queries = list(list(query = intersects, params = list("TCGA.TARGET","GSE60850","cellcycle_gene"), active = T)))
g1<-as.ggplot(p1) # 转换为ggplot2图片

library(yyplot)

g5<-g1 + geom_subview(subview = p2 + theme_void(), x=.7, y=.7, w=.5, h=.5)
g5

ggsave(file.path(fig1.path,"upsetR Venn.pdf"),width = 10,height = 5)

gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are TARGETepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}

# load meta signature
cellcycle.sig <- gmt2list(file.path(data.path,"pathways.gmt"))
cellcycle.sig <- sapply(cellcycle.sig, function(x) setdiff(x,""))
cellcycle.class <- NULL
for (i in names(cellcycle.sig)) {
  tmp <- cellcycle.sig[i]
  for (j in tmp) {
    cellcycle.class <- rbind.data.frame(cellcycle.class,
                                      data.frame(gene = j,
                                                 path = i,
                                                 stringsAsFactors = F),
                                      stringsAsFactors = F)
  }
}

# calculate GSVA enrichment score计算通路分数
cellcycle.score <- gsva(gsvaParam(as.matrix(TARGET.expr2[comgene,]), cellcycle.sig))
  
###保存计算出来的分数，备用，投稿的时候可能会要原始数据
write.table(cellcycle.score,file.path(fig1.path,"gsva enrichment score of cellcycle signature in tcga cohort.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
range(cellcycle.score)

library(ClassDiscovery) #用于外部聚类

annCol <-TARGET.clin[,c("OS","Age","Gender","Stage","Histologic")]
#annCol$PSA<-ifelse(annCol$PSA=="unknow","unknow",
 #                  ifelse(annCol$PSA<10,"<=10",">10"))

head(annCol)
#annCol[is.na(annCol) | annCol == ""] <- "N/A"
##设定每个因子的颜色
annCol$OS<-ifelse(annCol$OS=="0","Alive","Dead")
table(annCol$Race)
annColors <- list(OS  = c("Alive" = "grey",
                              "Dead"   = "black"),
                  Age    =  c("white", "#990000"),
                  Gender  = c("Male" = "#009999",
                           "Female"   = "#CC0066"),
                  Race  = c("unknow" = "grey",
                           "Asian"   = "#34A12E",
                           "White"   = "#FF7E02",
                           "BoAA"   = "#53207C"),
                  laterality = c(
                    "Left"    = "#F6BF02",
                    "Right"    = "#00A16D",
                    "Bilateral"    = "#E24A19"),
                  Stage = c(
                    "Stage I"    = yjp[1],
                    "Stage II"    = yjp[2],
                    "Stage III"    = yjp[3],
                    "Stage IV"    = yjp[4],
                    "[Discrepancy]"    = "grey"))

annColors <- list(
  OS = c(
    "Alive" = "#B0B0B0",   # 柔和灰
    "Dead"  = "#2C2C2C"    # 深灰近黑
  ),
  Age = c("#E8F5E9", "#B71C1C"),     # 深红
  Gender = c(
    "Male"   = "#1E88E5",  # 蓝
    "Female" = "#D81B60"   # 洋红
  ),
  Histologic = c(
    "FHWT"      = "#29B6F6",  # 淡蓝
    "DAWT" = "#EF5350"   # 红
  ),
  Stage = c(
    "I"         = "#A5D6A7",   # 绿
    "II"        = "#FFF176",   # 黄
    "III"       = "#FFB74D",   # 橙
    "IV"        = "#E57373"   # 红
  )
)


annColors
#画图看一下基本情况
#pdf(file.path(fig1.path, "heatmap in tcga.pdf"), width = 10,height = 12)
pheatmap(cellcycle.score,
         scale = "row",
         color = c(
           "#2421BA","#63BFA5","#F8FCB4","#EC6146","#D41714"),
         annotation_col = annCol,
         #cutree_cols = 2,
         #kmeans_k  = 2,
         annotation_colors = annColors,
         show_rownames = T, show_colnames = F,
         filename = "raw_heatmap.pdf")
#dev.off()


#############批量尝试不同聚类分组方法组合。

library(pheatmap)
library(survival)
library(survminer)
library(dendextend)
cellcycle.score <- gsva(gsvaParam(as.matrix(TARGET.expr2[comgene,]), cellcycle.sig))
# 距离与连接方法组合
distance_methods <- c("euclidean", "maximum", "manhattan", "canberra",  "minkowski")
linkage_methods <- c("ward.D", "ward.D2",  "complete", "average")

mat<-cellcycle.score
# 设置调色板
jco3 <- c("#FFCC33", "#009999", "#CC0066", "#146FB5")

# 创建输出文件夹
dir.create("survival_plots", showWarnings = FALSE)

# 循环每一种组合
for (dist in distance_methods) {
  for (link in linkage_methods) {
    message(sprintf("Processing: %s + %s", dist, link))
    
    #dist<-"euclidean";link<-"ward.D"
    # 聚类
    hcs <- hclust(distanceMatrix(as.matrix(mat), metric = dist), method = link)
    group <- cutree(hcs, k = 4)
    
    # 添加分组信息
    TARGET.clin$Clust <- paste0("C", group[rownames(TARGET.clin)])
    
    # 生存分析
    fit <- survfit(Surv(OS.time, OS) ~ Clust, data = TARGET.clin)
    
    # 保存生存曲线
    plot_file <- file.path("survival_plots", paste0("survplot_", dist, "_", link, ".pdf"))
    pdf(plot_file, width = 6, height = 5)
    print(ggsurvplot(fit, pval = TRUE, risk.table = FALSE,
                     ggtheme = theme_minimal(),
                     title = paste(dist, "+", link)))
    dev.off()
  }
}






##聚类分组，调整，选取差别最明显的分组。可以尝试不同的k值
mat<-cellcycle.score

hcs <- hclust(distanceMatrix(as.matrix(mat), "manhattan"), "ward.D") # 请阅读distanceMatrix()以及hclust()，了解更多distance测度和linkage方法
#the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
# one of "ward.D", "ward.D2", "single", "complete", "average" 
#hcg <- hclust(distanceMatrix(t(as.matrix(mat)), "maximum"), "ward.D") # 注意距离函数是针对列的，所以对行聚类要转置
group <- cutree(hcs,k=4)

# 增加一行annotation，及其配色
annCol$Clust2 <- paste0("C",group[rownames(annCol)])
annColors[["Clust2"]] <- c("C1"=jco3[1],"C2"=jco3[2],"C3"=jco3[3],"C4"=jco3[4])
#head(annCol)
pdf(file.path(fig1.path, "heatmap in tcga.pdf"), width = 10,height = 6)
pheatmap(mat,
         scale = "row",
         color = c("#4BFDF8","black", "#F40000"),
         cutree_cols = 4,
         #cluster_rows = hcs,
         cluster_cols = hcs,
         annotation_col = annCol,
         annotation_colors = annColors,
         show_rownames = T,show_colnames = F)
dev.off()

####给分组加一个前缀，123变成C1C2C3
TARGET.clin$Clust<-paste0("C",group[rownames(annCol)])
table(TARGET.clin$Clust,TARGET.clin$OS)

library(survival)
library("survminer")

#TARGET.clin$OS.time<-TARGET.clin$OS.time/30.5
outTab=data.frame()
fit <- survfit(Surv(OS.time, OS) ~ Clust, data = TARGET.clin)
ggsurvplot(fit,pval = T)
pval_info <- surv_pvalue(fit, data = TARGET.clin)
pval1 <- pval_info$pval
pval1

library(survival)

# 保留两组数据比较（例如 C1 vs C2）
subset_data <- TARGET.clin[TARGET.clin$Clust %in% c("C1", "C4"), ]
#subset_data$Clust <- droplevels(subset_data$Clust)
# 将 Clust 转为因子，并设置 C2 为 reference（基线）
subset_data$Clust <- factor(subset_data$Clust, levels = c("C4", "C1"))

cox_model <- coxph(Surv(OS.time, OS) ~ Clust, data = subset_data)

coxSummary = summary(cox_model)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
outTab<-data.frame()
outTab=rbind(outTab,
             cbind(
               HR=coxSummary$conf.int[,"exp(coef)"],
               HR.95L=coxSummary$conf.int[,"lower .95"],
               HR.95H=coxSummary$conf.int[,"upper .95"],
               pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
)
outTab

HR <- paste("HR = ", round(outTab$HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(outTab$HR.95L,3), round(outTab$HR.95H,3), sep = " - "), sep = "")


names(fit$strata) <- gsub("Clust=", "", names(fit$strata))

pdf(file=file.path(fig1.path,"TARGET cellcycle subgroup survival.pdf"), width=5,height=7,onefile = FALSE)
ggsurvplot(fit, data = TARGET.clin,
           risk.table = T, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           #surv.median.line = "hv", # Specify median survival
           ggtheme = theme_classic2(), 
           conf.int = F, 
           conf.int.style = "step",
           censor = T, 
           palette = jco3, #
           ylim = c(0,1),
           xlab = 'Time in months',
           #risk.table.y.text.col = T, 
           #risk.table.y.text = T, 
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"),
           #fun = "event",
           legend.title = NULL,
           pval = paste(pval1 = ifelse(pval1 < 0.001, "p < 0.001", 
                                     paste("Overall P = ",round(pval1,4), sep = "")),
                        pval2 = ifelse(outTab$pvalue < 0.001, "p < 0.001", 
                                      paste("C1 vs. C4: P = ",round(outTab$pvalue,4), sep = "")),
                        HR, CI, sep = "\n"),
           pval.coord = c(0, 0.15)
)

dev.off()


########用log-rank 两两比较（可选）
fit <- survfit(Surv(OS.time, OS) ~ Clust, data = TARGET.clin)
ggsurvplot(fit,pval = T)

surv.obj <- survfit(Surv(OS.time, OS) ~ Clust, data = TARGET.clin)
surv.obj
#分组比较生存曲线
p1 <- ggsurvplot(surv.obj,pval = F,conf.int = F,
                 risk.table = T,
                 surv.median.line = "hv",
                 risk.table.col = "strata",
                 ggtheme = theme_bw(),
                 palette = jco3)
p1
restest <- pairwise_survdiff(Surv(OS.time, OS) ~ Clust, data = TARGET.clin)
restest

restest[["p.value"]]
pvalues <- as.data.frame(restest[["p.value"]])
pvalues <- round(pvalues,3)
pvalues[pvalues<0.001] <- '<0.001'
pvalues[is.na(pvalues)] <- '-'
pvalues <- rownames_to_column(pvalues,var = '  ')
pvalues


############twoclass limma 寻找两组间差异基因，即组蛋白修饰最强和最弱组间差异基因
source(file.path(Scripts,"twoclasslimma.R"))
condition<-as.matrix(TARGET.clin$Clust)
colnames(condition)<-c("condition")

subt <- data.frame(condition = condition,
                   row.names = rownames(TARGET.clin))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = TARGET.expr2[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "C1", # name of treatment group
              ctrlVar  = "C4", # name of control group
              prefix   = "TARGET", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path = fig1.path) # path for result

# extract group specific pathways
tmp1 <- read.table(file.path(fig1.path,"TARGET_limma_test_result.C1_vs_C4.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

###########筛选差异基因，筛选与之可调整

deg1up <- rownames(tmp1[which(tmp1$log2fc > 0.4 & tmp1$padj < 0.05),])
deg1dw <- rownames(tmp1[which(tmp1$log2fc < -0.4& tmp1$padj < 0.05),])

genes<-c(deg1up,deg1dw)
genes<-intersect(genes,comgene)

##保存
write.table(genes,file=file.path(fig1.path,"intersectgene.txt"),sep="\t",row.names=T,quote=F)

##########画火山图###########

library(ggplot2)
library(ggrepel) 

# 分组
tmp1$change <- ifelse(tmp1$pvalue < 0.01 & abs(tmp1$log2fc) >= 0.5, 
                      ifelse(tmp1$log2fc > 0.5, 'Up', 'Down'),
                      'Stable')

# 加 symbol 列
tmp1$symbol <- rownames(tmp1)

# 提取上下调 top10（按最小 p 值）
top10_up <- tmp1[tmp1$change == "Up", ]
top10_up <- top10_up[order(top10_up$pvalue), ][1:10, ]

top10_down <- tmp1[tmp1$change == "Down", ]
top10_down <- top10_down[order(top10_down$pvalue), ][1:10, ]

# 合并为一个数据框
top_genes <- rbind(top10_up, top10_down)


# volcano plot
pDEG <- ggplot(data = tmp1, 
            aes(x = log2fc, 
                y = -log10(pvalue), 
                colour = change,
                shape = change,
                label = symbol)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = c("Down" = "#3399CC", 
                                "Stable" = "lightgrey", 
                                "Up" = "#EF8A09")) +
  scale_shape_manual(values = c("Stable" = 16,  # 圆
                                "Down" = 18, # 菱形
                                "Up" = 17)) +  # 三角形
  geom_vline(xintercept = c(-0.5, 0.5), lty = 4, col = 'black', lwd = 0.8) +
  geom_hline(yintercept = -log10(0.01), lty = 4, col = 'black', lwd = 0.8) +
  geom_text_repel(data = top_genes,
                  aes(label = symbol),
                  max.overlaps = 1000,
                  size = 3.5,
                  box.padding = 0.4,
                  segment.color = "grey50") +
  labs(x = 'log2(fold change)',
       y = '-log10 (p-value)',
       title = 'cellcycle modification associated genes') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.1, 0.8), 
        legend.title = element_blank())

# 显示图
print(pDEG)

# 保存图
ggsave(file.path(fig1.path, "cellcycle volcano.pdf"), 
       pDEG, width = 5, height = 5)

#--------新型火山图---------
diff_express2<-tmp1
diff_express2$name<-rownames(diff_express2)
diff_express2$log2FoldChange<-diff_express2$log2fc
# Default plot
ggmaplot(diff_express2, main = expression("C1" %->% "C3"),
         fdr = 0.05, fc = 1.5157, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(diff_express2$name),
         legend = "top", top = 20,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())
ggsave(file.path(fig1.path,"cellcycle volcano new.pdf"),width = 6,height = 5)


#####差异基因通路富集
library("clusterProfiler")
options(connectionObserver = NULL)
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

#
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
gene<-entrezIDs <- as.character(entrezIDs)

genelist_input<-tmp1
genelist_input$enzid<-mget(genelist_input$symbol, org.Hs.egSYMBOL2EG, ifnotfound=NA)

genelist_input<-subset(genelist_input,genelist_input$enzid!="NA")

geneList = genelist_input[,2]
names(geneList) = as.character(genelist_input[,7])
geneList = sort(geneList, decreasing = TRUE)
edo2 <- gseDO(geneList)
ridgeplot(edo2)

#---Go通路富集
Go<- enrichGO(gene = gene,
              keyType = "ENTREZID",
              OrgDb = org.Hs.eg.db, 
              pvalueCutoff =0.05, 
              qvalueCutoff = 0.5,
              ont="all",
              readable =T)
write.table(Go,file.path(fig1.path,file="GO-up.txt"),sep="\t",quote=F,row.names = F)

upsetplot(Go)
library(DOSE)

###保存点图
pdf(file=file.path(fig1.path,"Go-up.pdf"),width = 10,height = 6)
dotplot(Go,showCategory = 10,label_format = 70,split="ONTOLOGY") + ##label_format 改前面名称显示的长度
  facet_grid(ONTOLOGY~., scale='free') + #是否根据BP，CC，MF分类
  scale_color_continuous(low='#009966', high='#FF0034')+ #设置颜色
  aes(shape=GeneRatio > 0.04)#设置三角形圆形
dev.off()

###保存调控网络图
pdf(file=file.path(fig2.path,"Go-tree.pdf"),width = 10,height = 6)
Go2<- pairwise_termsim(Go)
#treeplot(Go2)
p <- treeplot(Go2)
p + scale_color_viridis()+ scale_fill_manual(values = c("#FFA07A", "#87CEEB", "#98FB98",'#FF0034','#009966'))
dev.off()

library(viridis)

###保存通路联系、关键基因图

pdf(file=file.path(fig2.path,"Go-corelation.pdf"),width = 6,height = 6)
#cnetplot(Go2, showCategory = 5,categorySize="count", colorEdge = TRUE,foldChange=d)

cnetplot(Go, categorySize="geneNum",
         showCategory = 10, 
         foldChange=geneList)+scale_color_continuous(low='#009999', high='#53207C')

dev.off()


#---KEGG分析
Kegg <- enrichKEGG(gene = gene, organism = "hsa",
                   pvalueCutoff =0.5, qvalueCutoff =0.5)
write.table(Kegg,file.path(fig2.path,"KEGG-up.txt"),sep="\t",quote=F,row.names = F)


#设定你只想展示的通路
intest<-c("Nucleocytoplasmic transport","Cell cycle")
#          "Signaling pathways regulating pluripotency of stem cells",
#          "Axon guidance","Cellular senescence",
#          "mTOR signaling pathway","Proteoglycans in cancer",
#          "Ubiquitin mediated proteolysis","Phagosome")

pdf(file=file.path(fig2.path,"KEGG-up.pdf"),width = 8,height = 4)
dotplot(Kegg, showCategory = 10)+ scale_fill_gradientn(colors = yjp2)
dev.off()

Kegg2<- pairwise_termsim(Kegg)
Kegg2<-setReadable(Kegg2, 'org.Hs.eg.db', 'ENTREZID')

pdf(file=file.path(fig2.path,"KEGG-net.pdf"),width = 6,height = 6)
#test<-c("ECM-receptor interaction","Focal adhesion","PI3K-Akt signaling pathway")
#cnetplot(Kegg2, showCategory = 5,categorySize="count",colorEdge = TRUE)
cnetplot(Kegg2, categorySize="geneNum",
         showCategory = 10, 
         foldChange=geneList)+scale_color_continuous(low='#009999', high='#53207C')
dev.off()


pdf(file="Kegg-net.pdf",width = 10,height = 8)
emapplot(Kegg2, showCategory = 20,color = "p.adjust")+ scale_color_continuous(low='yellow', high='red')
dev.off()

####保存富集通路树图
pdf(file.path(fig2.path,"Kegg-tree.pdf"),width = 10,height = 5)
treeplot(Kegg2,showCategory = 30)


p <- treeplot(Kegg2,showCategory = 20)
p + scale_color_viridis()+ scale_fill_manual(values = c("#FFA07A", "#87CEEB", "#98FB98",'#FF0034','#009966'))
dev.off()

upsetplot(Kegg2)

heatplot(Kegg2, showCategory=5)

#---HALLMARK分析
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
gene<-entrezIDs <- as.character(entrezIDs)
library(msigdbr)
msigdbr_species() #支持的物种
#Hu_msigdbr <- msigdbr(species="Homo sapiens")
#head(Hu_msigdbr, 2) %>% as.data.frame

HALLMARK<- msigdbr(species="Homo sapiens",category="H") %>% 
  dplyr::select(gs_name, entrez_gene, gene_symbol)
head(HALLMARK)

HALLMARKresult<- enricher(gene,TERM2GENE=HALLMARK[,c(1,2)])

write.table(HALLMARKresult,file=file.path(fig1.path,"HALLMARKresult.txt"),sep="\t",quote=F,row.names = F)

pdf(file=file.path(fig1.path,"HALLMARKresult.pdf"),width = 6,height = 4)
dotplot(HALLMARKresult, showCategory = 10)+ scale_color_continuous(low='#009966', high='#FF0034')+ aes(shape=GeneRatio > 0.1)
dev.off()

##################################
#################
###计算每个队列中，2620个组蛋白差异基因的预后预测能力
library(survival)
rt=cbind(TARGET.clin[,c("OS.time","OS")],t(TARGET.expr2))
#rt=cbind(TARGET.clin[,c("OS.time","OS")],t(TARGET.expr))
TARGEToutTab=data.frame()
pb <- progress_bar$new(total=ncol(rt)-2)
for(i in colnames(rt[,3:ncol(rt)])){
  rt[,i]<- factor(ifelse(rt[,i]>median(rt[,i]),"high"," low"))
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  TARGEToutTab=rbind(TARGEToutTab,
                     cbind(id=i,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  pb$tick()
  Sys.sleep(0.05)
}
rownames(TARGEToutTab)<-TARGEToutTab$id
write.table(TARGEToutTab,file.path(fig2.path,"TARGETuniCox.txt"),sep="\t",row.names=F,quote=F)


rt=cbind(GSE60850.clin[,c("OS.time","OS")],t(GSE60850.expr[genes,]))
#rt=cbind(GSE60850.clin[,c("OS.time","OS")],t(GSE60850.expr))
GSE60850outTab=data.frame()
pb <- progress_bar$new(total=ncol(rt)-2)
for(i in colnames(rt[,3:ncol(rt)])){
  rt[,i]<- factor(ifelse(rt[,i]>median(rt[,i]),"high"," low"))
  cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  GSE60850outTab=rbind(GSE60850outTab,
                       cbind(id=i,
                             HR=coxSummary$conf.int[,"exp(coef)"],
                             HR.95L=coxSummary$conf.int[,"lower .95"],
                             HR.95H=coxSummary$conf.int[,"upper .95"],
                             pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  pb$tick()
  Sys.sleep(0.05)
}

rownames(GSE60850outTab)<-GSE60850outTab$id
write.table(GSE60850outTab,file.path(fig2.path,"GSE60850uniCox.txt"),sep="\t",row.names=F,quote=F)




Risk1<-rownames(TARGEToutTab[which(TARGEToutTab$pvalue<0.05&TARGEToutTab$HR>1),])
Prot1<-rownames(TARGEToutTab[which(TARGEToutTab$pvalue<0.05&TARGEToutTab$HR<1),])


Risk4<-rownames(GSE60850outTab[which(GSE60850outTab$pvalue<0.05&GSE60850outTab$HR>1),])
Prot4<-rownames(GSE60850outTab[which(GSE60850outTab$pvalue<0.05&GSE60850outTab$HR<1),])

survgene<-c(intersect(Risk4,Risk1),intersect(Prot4,Prot1))


deg1up <- rownames(tmp1[which(tmp1$log2fc > 0.5 & tmp1$padj < 0.05),])
deg1dw <- rownames(tmp1[which(tmp1$log2fc < -0.5& tmp1$padj < 0.05),])

genes<-c(deg1up,deg1dw)

outgene1<-intersect(Risk1,deg1up)
outgene2<-intersect(Prot1,deg1dw)

survgene<-c(outgene1,outgene2)

########################################################Venn图
library("ggVennDiagram")
# Default plot
venn1<-ggVennDiagram(list(Risk1,deg1up),
                     label_alpha = 0,
                     category.names = c("TARGET cohort","upregulate group"))+
  ggplot2::scale_fill_gradient(low="#2A9D8F",high = "#FFBE0B")
venn1
ggsave(file.path(fig2.path,"Risky gene Venn.pdf"),width = 5,height = 3)

venn2<-ggVennDiagram(list(Prot1,deg1dw),
                     label_alpha = 0,
                     category.names = c("TCGA-TARGET cohort","GSE60850 cohort"))+
  ggplot2::scale_fill_gradient(low="#3A86FF",high = "#8E24AA")
venn2
ggsave(file.path(fig2.path,"Protective gene Venn.pdf"),width = 5,height = 3)

"#E76F51" "#2A9D8F" "#E9C46A" "#264653" "#8E24AA""#3A86FF" "#FFBE0B" "#8338EC"
##############################
library(yyplot)
library(ggplot2)

TARGEToutTab$symbol<-rownames(TARGEToutTab)
TARGEToutTab$pvalue<-as.numeric(TARGEToutTab$pvalue)
TARGEToutTab$HR<-as.numeric(TARGEToutTab$HR)
TARGEToutTab$change = ifelse(TARGEToutTab$pvalue < 0.05, 
                           ifelse(TARGEToutTab$HR> 1 ,'Riskey','Protective'),
                           'Stable')
TARGETvolcano <- ggplot(data = TARGEToutTab, 
                      aes(x = TARGEToutTab$HR, 
                          y = -log10(TARGEToutTab$pvalue), 
                          colour=change,
                          label = TARGEToutTab$symbol)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c('#00B072', '#FDB530','#07346C'))+
  xlim(c(0, 3)) +
  geom_vline(xintercept=c(1),lty=4,col='black',lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col='black',lwd=0.8) +
  labs(x='Hazard Ratio',
       y='-log10 (p-value)',
       title='TCGA-TARGET cohort')  +
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5), 
    legend.position=c(0.8,0.3), 
    legend.title = element_blank())
TARGETvolcano

GSE60850outTab$symbol<-rownames(GSE60850outTab)
GSE60850outTab$pvalue<-as.numeric(GSE60850outTab$pvalue)
GSE60850outTab$HR<-as.numeric(GSE60850outTab$HR)
GSE60850outTab$change = ifelse(GSE60850outTab$pvalue < 0.05, 
                               ifelse(GSE60850outTab$HR> 1 ,'Riskey','Protective'),
                               'Stable')
GSE60850volcano <- ggplot(data = GSE60850outTab, 
                          aes(x = GSE60850outTab$HR, 
                              y = -log10(GSE60850outTab$pvalue), 
                              colour=change,
                              label = GSE60850outTab$symbol)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c('#E8863D','#732A7C', '#299D92'))+
  xlim(c(0, 4)) +
  geom_vline(xintercept=c(1),lty=4,col='black',lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col='black',lwd=0.8) +
  labs(x='Hazard Ratio',
       y='-log10 (p-value)',
       title='GSE60850 cohort')  +
  theme_classic()+
  theme(#plot.title = element_text(hjust = 0.5), 
    legend.position=c(0.8,0.3), 
    legend.title = element_blank())
GSE60850volcano


####拼图
ggarrange(TARGETvolcano, pDEG,venn1,venn2,
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

ggsave(file.path(fig2.path,"cohort HR volcano.pdf"),width = 8,height = 8)

###########
rt=cbind(annCol,t(TARGET.expr2[survgene,rownames(annCol)]))
#----------热图---------
library(dplyr)
library(pheatmap)                   #引用包
head(rt)
rt<-rt[order(rt$Clust2,decreasing = F),]

rt2 <- rt %>% 
  dplyr::select(-c(1:6)) 
rt2<-log2(rt2+1)

rt2 <- scale(t(rt2)) #scale标化

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(rt2), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

rt3 <- standarize.fun(rt2,halfwidth =0.5)

head(rt3[1:3,1:3])


cluster<- rt %>% 
  dplyr::select(c(1:6)) 


#绘制热图
pdf(file.path(fig2.path,"TARGET_heatmap_79genes.pdf"),height=5,width=10)

pheatmap(rt3,
         scale = "row",
         color = c("#23B1EB","black", "#FBFE00"),
         cluster_rows = T,
         cluster_cols = F,
         annotation_col = cluster,
         annotation_colors = annColors,
         show_rownames = F,show_colnames = F)
dev.off()


###########GSE60850
rt=cbind(GSE60850.clin,t(GSE60850.expr[survgene,rownames(GSE60850.clin)]))
#----------热图---------
library(dplyr)
library(pheatmap)                   #引用包
head(rt)
rt<-rt[order(rt$OS,decreasing = F),]

rt2 <- rt %>% 
  dplyr::select(-c(1:5)) 

rt2 <- scale(t(rt2)) #scale标化

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(rt2), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

rt3 <- standarize.fun(rt2,halfwidth =1)

head(rt3[1:3,1:3])


cluster<- rt %>% 
  dplyr::select(c(2:4)) 
table(cluster$Type)
annColors2 <- list(OS  = c("0" = "grey",
                          "1"   = "black"),
                  Gender  = c("male" = "#009999",
                              "female"   = "#CC0066"),
                  Type = c(
                    "metastasis"    = "#F6BF02",
                    "primary"    = "#00A16D"))
#绘制热图
pdf(file.path(fig2.path,"GSE60850_heatmap_192genes.pdf"),height=8,width=5)

pheatmap(rt3,
         scale = "row",
         color = c("#32AEEC","black", "#FDE805"),
         cluster_rows = T,
         cluster_cols = T,
         annotation_col = cluster,
         annotation_colors = annColors2,
         show_rownames = F,show_colnames = F)
dev.off()





########################
###在TCGA中训练预后预测模型
TARGET.expr3<-log2(TARGET.expr2+1)
LASSO.input<-cbind(TARGET.clin[,1:2],t(TARGET.expr3[survgene,]))

#-------------LASSO----------

df<-LASSO.input
dim(df)
head(df)
mydata<-df[,3:ncol(df)]
mydata<-as.matrix(mydata)

#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")
#####set.seed是为了保证每次重复结果都一样。如果一次结果不理想，可以调整seed

df<-LASSO.input
dim(df)
head(df)
mydata<-df[,3:ncol(df)]
mydata<-as.matrix(mydata)

#install.packages("glmnet")
#install.packages("survival")
library("glmnet")
library("survival")
#####set.seed是为了保证每次重复结果都一样。如果一次结果不理想，可以调整seed

#n<-round(runif(1, min = 0, max = 1000000),0)

set.seed(510181)
#做10倍交叉验证，算出lambda值
cvfit = cv.glmnet(mydata, Surv(df$OS.time,df$OS), 
                  family = "cox",
                  nfold=20) #10倍交叉验证
cvfit$lambda.min
plot(cvfit)


coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
geneids <- colnames(mydata)[active.min]
geneids
index.min = coef.min[active.min]
index.min

fit <- glmnet(mydata, Surv(df$OS.time,df$OS), 
              family = "cox",
              nfold=100) #10倍交叉验证


####fit plot 美化图片

library(survival)
library(glmnet)
library(ggplot2)
library(ggsci)

x <- coef(fit)  
tmp <- as.data.frame(as.matrix(x))
tmp$coef <- row.names(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1] # extract the lambda values
tmp$norm <- apply(abs(x[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm  


ggplot(tmp,aes(log(lambda),value,color = coef)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_line(size=1) + 
  xlab("Lambda (log scale)") + 
  #xlab("L1 norm")+
  ylab('Coefficients')+
  theme_bw(base_rect_size = 2)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12),
                                pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12),
                                pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0.01,0.01))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+
  #annotate('text',x = -3.3,y=0.26,label='Optimal Lambda = 0.012',color='black')+
  guides(col=guide_legend(ncol = 2))

ggsave(filename = file.path(fig2.path,"LASSO.fit.pdf"), width = 8,height = 5)

## 准备数据LASSO.cvfit.
xx <- data.frame(lambda=cvfit[["lambda"]],cvm=cvfit[["cvm"]],cvsd=cvfit[["cvsd"]],
                 cvup=cvfit[["cvup"]],cvlo=cvfit[["cvlo"]],nozezo=cvfit[["nzero"]])
xx$ll <- log(xx$lambda)
xx$NZERO <- paste0(xx$nozezo,' vars')

ggplot(xx,aes(ll,cvm,color=NZERO))+
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),width=0.05,size=0.5)+
  geom_vline(xintercept = xx$ll[which.min(xx$cvm)],size=0.8,color='grey60',alpha=0.8,linetype=2)+
  geom_point(size=2)+
  xlab("Log Lambda")+ylab('Partial Likelihood Deviance')+
  theme_bw(base_rect_size = 1.5)+ 
  scale_color_manual(values = c(pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12),
                                pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12),
                                pal_npg()(10),pal_d3()(10),pal_jco()(7),
                                pal_jama()(7),pal_uchicago()(9),
                                pal_gsea()(12),pal_ucscgb()(12)))+
  scale_x_continuous(expand = c(0.02,0.02))+
  scale_y_continuous(expand = c(0.02,0.02))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15,color='black'),
        axis.text = element_text(size=12,color='black'),
        legend.title = element_blank(),
        legend.text = element_text(size=12,color='black'),
        legend.position = 'right')+ #记得修改下方的最佳lambda值
  annotate('text',x = -6.8,y=13.2,label=paste0('Optimal Lambda = ',round(cvfit$lambda.min,3)),color='black')+
  guides(col=guide_legend(ncol = 2))
ggsave(filename = file.path(fig2.path,"LASSO.cvfit.pdf"), width = 7.5,height = 5)

cvfit$lambda.min #最佳lambda值

cvfit$lambda.1se #一倍SE内的更简洁的模型

# 输出基因顺序
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
geneids <- colnames(mydata)[active.min]
geneids
index.min = coef.min[active.min]

combine<-cbind(geneids, index.min)
write.csv(combine,file = file.path(fig3.path,"gene_index.csv"))

# 将纳入signature的变量拟合成一个变量，作为nomogram的输入

geneids
index.min

signature <- as.matrix(df[, geneids]) %*% as.matrix(index.min) 
write.csv(signature,file = file.path(fig3.path,"TARGET-signature.csv"))

total_signature<-signature1<-signature

TARGET.os<-LASSO.input[,1:2]

TARGET_KM_input<-cbind(TARGET.os,signature)
TARGET_KM_input1<-cbind(TARGET.os,signature1)

#------------K-M

library(survival)
library("survminer")

TARGET_KM_input$signature<- factor(ifelse(TARGET_KM_input$signature>median(TARGET_KM_input$signature),"high"," low"))
#TARGET_KM_input$OS.time<-TARGET_KM_input$OS.time/30.5
outTab=data.frame()
fit <- survfit(Surv(OS.time, OS) ~ signature, data = TARGET_KM_input)

cox <- coxph(Surv(OS.time, OS) ~ signature, data = TARGET_KM_input)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
outTab=rbind(outTab,
             cbind(
               HR=coxSummary$conf.int[,"exp(coef)"],
               HR.95L=coxSummary$conf.int[,"lower .95"],
               HR.95H=coxSummary$conf.int[,"upper .95"],
               pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
)
outTab

HR <- paste("Hazard Ratio = ", round(outTab$HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(outTab$HR.95L,3), round(outTab$HR.95H,3), sep = " - "), sep = "")


pdf(file=file.path(fig3.path,"TARGET survival.pdf"), width=4,height=4.2,onefile = FALSE)
ggsurvplot(fit, data = TARGET_KM_input,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), 
           conf.int = T, 
           conf.int.style = "ribbon",
           censor = T, 
           palette = jco5, #
           ylim = c(0,1),
           
           xlab = 'Time in months',
           legend.title='riskscore', 
           legend.labs=c('Low','High'), 
           
           risk.table.y.text.col = T, 
           risk.table.y.text = T, 
           
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"),
           #fun = "event",
           
           pval = paste(pval = ifelse(outTab$pvalue < 0.001, "p < 0.001", 
                                      paste("P = ",round(outTab$pvalue,3), sep = "")),
                        HR, CI, sep = "\n"),
           pval.coord = c(0, 0.15)
)

dev.off()

outTab


#----------ROC————————————

library(timeROC)
library(survival)

pdf(file=file.path(fig3.path,"TCGA HS ROC 1 3 5.pdf"), width=4,height=4.5)

ROC.DSST<-timeROC(T=TARGET_KM_input1$OS.time,#结局时间
                  delta=TARGET_KM_input1$OS,#生存结局
                  marker=TARGET_KM_input1$signature1,#预测变量
                  cause=1,#阳性结局赋值，比如死亡，复发的赋值
                  weighting="marginal",# 权重计算方法，marginal是默认值，采用km计算删失分布
                  times=c(12,36,30),# 时间点，选取10年和20年生存率
                  ROC = TRUE,
                  iid = TRUE
)
plot(ROC.DSST,time=12,col=jco2[1],title=FALSE,lwd=2)
plot(ROC.DSST,time=36,col=jco2[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=30,col=jco2[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC.DSST$AUC[1],3)),
         paste0('AUC at 3 years: ',round(ROC.DSST$AUC[2],3)),
         paste0('AUC at 5 years: ',round(ROC.DSST$AUC[3],3))),
       col=jco2,lwd=2,bty = 'n')
dev.off()

#####################Risk table
library(reshape2)
library(ggplot2)
library(scales)
library(cowplot)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

data <- cbind(LASSO.input,signature)
data[1:2, 1:4]

bestvars <-geneids
bestvars

# risk score，用于画顶部散点图
rs <- data$signature
names(rs) <- rownames(data)
rs_data <- data.frame(x=1:length(rs),rs=as.numeric(sort(rs)))
# 用中值分组
rs_data$Risk <- ifelse(rs_data$rs>=median(rs_data$rs), "High-risk", "Low-risk")
head(rs_data)


# follow-up，用于画中间B图
surv_data <- data.frame(x=1:length(rs),
                        t=data[names(sort(rs)),'OS.time']/12*12,
                        s=data[names(sort(rs)),'OS']) 
surv_data$Status <- as.factor(ifelse(surv_data$s==0,'Alive','Death'))
head(surv_data)

# 提取signature对应的data，并按risk score排序，用于画底部热图
exp_data <- data[names(sort(rs)),which(colnames(data) %in% bestvars)]
exp_data[1:2,1:4]

plot.A <- ggplot(rs_data, aes(x=x,y=rs))+
  geom_point(aes(col=Risk),size=1.5)+
  scale_color_manual(labels=c("High-risk","Low-risk"), 
                     #guide_legend(guide = NULL), #如果不想画图例就删掉#
                     name="Risk score", values =c("#67008F","#63C16A")) + 
  
  geom_segment(aes(x = sum(rs_data$Risk=="Low-risk"),
                   y = min(rs_data$rs), 
                   xend = sum(rs_data$Risk=="Low-risk"), 
                   yend = max(rs_data$rs)), linetype="dashed", size = 0.6)+
  # 画横线
  geom_segment(aes(x=0,y=median(rs_data$rs),
                   xend=nrow(rs_data),
                   yend=median(rs_data$rs)),linetype="dashed", size = 0.3)+
  
  # 写文字Cutoff:
  #geom_text(aes(x=sum(rs_data$Risk=="Low-risk")/2,
  #              y=median(rs_data$rs)+8,
  #              label=paste0("Cutoff: ",round(median(rs_data$rs),3))),
  #          col ="black",size = 4,alpha=0.8)+
  
  theme(axis.title.x=element_blank()) +
  scale_x_continuous(limits = c(0,NA),expand = c(0,0)) +
  labs(y="Risk score",x="",fill="Risk") +
  #scale_colour_discrete(name="Risk scores") +
  theme_classic() +
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样画坐标轴，就删掉这行
        axis.text.x=element_blank())

plot.A


plot.B <- ggplot(surv_data,aes(x=x,y=t))+
  geom_point(aes(col=Status),size=1.5)+
  geom_vline(aes(xintercept=sum(rs_data$Risk=="Low-risk")),size=1,linetype="dashed")+
  scale_x_continuous(limits = c(0,NA),expand = c(0,0))+
  scale_color_manual(labels=c("Alive","Dead"),
                     values =c("#63C16A","#67008F"))+
  labs(y="RFS(months)",x="")+
  theme_classic()+
  theme(axis.ticks.x=element_blank(),
        axis.line = element_blank(), #如果想像example2那样不画坐标轴，就删掉前面的#
        axis.text.x=element_blank())

plot.B

tmp <- t(scale(exp_data))
tmp[tmp > 1] = 1
tmp[tmp < -1] = -1
reorder_cormat <- function(cormat){
  dd <- dist(cormat)
  hc <- hclust(dd,method = "average")
  cormat <-cormat[hc$order,]
}
tmp1 <- reorder_cormat(tmp)
tmp1 <- rbind(tmp1,ifelse(rs_data$Risk=="Low-risk",-1.5,1.5))
tmp.m <- reshape2::melt(tmp1)

p2 <-ggplot(tmp.m, aes(Var2, Var1),size=1) + 
  geom_tile(aes(fill = value)) 

plot.C <- p2 + scale_fill_gradient2(name="Genes\nexpression", low="#63C16A", 
                                    high="#67008F", mid="black") +
  labs(x = "", y = "")+
  theme_classic()+
  theme(legend.title = element_text(size = 12), legend.position = "right",
        axis.line = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank())

plot.C

TARGETriskmap<-plot_grid(plot.A, plot.B, plot.C,
          labels = c("", "",""), # 或者按顺序标注ABC
          rel_heights = c(1,1,2), # 3个图的比例
          #label_x=0,
          #label_y=1,
          align = 'v',ncol = 1, axis="lr", scale = c(1,1,1), greedy = F)
TARGETriskmap

ggsave(file.path(fig3.path,paste0("TARGET riskmap figures.pdf")), width = 6, height = 10)



#-----------TARGET-临床信息森林图-----------
  library("tidyverse")
#TARGET.sig<-signature
signature<-TARGET_KM_input1[rownames(TARGET.clin),"signature1"]
TARGETmulti<-cbind(TARGET.clin,signature)

TARGETmulti$Risk<- factor(ifelse(TARGETmulti$signature>median(TARGETmulti$signature),"high-risk"," low-risk"))
TARGETmulti<-TARGETmulti[,c("OS.time","OS","Age","Gender","Stage","Histologic","Risk")]
#TARGETmulti$Grade<-ifelse(TARGETmulti$Grade=="G1","G1+G2",
 #                       ifelse(TARGETmulti$Grade=="G2","G1+G2",TARGETmulti$Grade))

library(survival)
library(survminer)
TARGETmulti[TARGETmulti=="unknow"]<-NA
TARGETmulti[TARGETmulti==""]<-NA
TARGETmulti[TARGETmulti=="[Discrepancy]"]<-NA

TARGETmulti<-TARGETmulti %>% 
  drop_na()
head(TARGETmulti)



#构建模型
model <- coxph( Surv(OS.time, OS) ~., data =TARGETmulti,na.action=na.exclude )
coxSummary = summary(model)
coxSummary


outTab=data.frame()
outTab=cbind(
  HR=coxSummary$conf.int[,"exp(coef)"],
  HR.95L=coxSummary$conf.int[,"lower .95"],
  HR.95H=coxSummary$conf.int[,"upper .95"],
  pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

write.table(outTab,file=file.path(fig3.path,"TCGAmultiCox2.xls"),sep="\t",row.names=F,quote=F)

pdf(file=file.path(fig3.path,"TARGET TOTAL clinical forest.pdf"),onefile = FALSE,
    width = 6,             #图片的宽度
    height = 6,            #图片的高度
)

ggforest(model,
         data=TARGETmulti,
         main = "TARGET Cohort",
         cpositions = c(0.01,0.14,0.36), 
         fontsize = 1, 
         refLabel = "reference", 
         noDigits = 3)
dev.off()

model


####
library(ggpubr)
library(rstatix)
library(survminer)

boxtt<-TARGETmulti[,c("signature","Stage")]
boxtt$Stage<-ifelse(boxtt$Stage=="I","I+II",
                    ifelse(boxtt$Stage=="II","I+II","III+IV"))
boxtt[boxtt=="unknow"]<-NA
boxtt[boxtt==""]<-NA
boxtt[boxtt=="[Discrepancy]"]<-NA

boxtt<-boxtt %>% 
  drop_na()

boxtt<-boxtt[order(boxtt$Stage, decreasing = F), ]


compare_means(signature ~ Stage,  data = boxtt,
              ref.group = ".all.", method = "t.test")

ggboxplot(boxtt, x = "Stage", y = "signature", 
          fill = "Stage", 
          #color = "Stage",
          legend = "none") +
  geom_jitter(width=0.15, alpha=0.6)+
  ##rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(boxtt$signature[ which(!is.na(boxtt$Stage))]), linetype = 2)+
  stat_compare_means(method = "anova", label.y = 0.8)+
  scale_fill_brewer(palette="PuBuGn")+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")
# 画小提琴图 + 箱线图 + jitter
ggplot(boxtt, aes(x = Stage, y = signature, fill = Stage)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  geom_hline(yintercept = mean(boxtt$signature[!is.na(boxtt$Stage)]), 
             linetype = 2, color = "black") +
  stat_compare_means(method = "t.test", label.y = 0.8) +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  scale_fill_manual(values = c("#63C16A", "#67008F")) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggviolin(boxtt, x = "Stage", y = "signature", fill = "Stage",
         add = "boxplot", add.params = list(fill = "white")) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  geom_hline(yintercept = mean(boxtt$signature[!is.na(boxtt$Stage)]), 
             linetype = 2) +
  stat_compare_means(method = "anova", label.y = 0.8) +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  scale_fill_brewer(palette = "PuBuGn") +
  theme(legend.position = "none")

ggsave(file.path(fig4.path,paste0("Stage signature subgroup.pdf")), width = 4, height = 4)

#
library(ggpubr)
library(rstatix)
library(survminer)

boxtt<-TARGETmulti[,c("signature","Clust")]
boxtt[boxtt=="GX"]<-NA
boxtt[boxtt==""]<-NA
boxtt[boxtt=="[Discrepancy]"]<-NA

boxtt<-boxtt %>% 
  drop_na()

boxtt<-boxtt[order(boxtt$Clust, decreasing = F), ]


compare_means(signature ~ Clust,  data = boxtt,
              ref.group = ".all.", method = "t.test")

ggboxplot(boxtt, x = "Clust", y = "signature", 
          fill = "Clust", 
          #color = "Grade",
          legend = "none") +
  geom_jitter(width=0.15, alpha=0.6)+
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(boxtt$signature[ which(!is.na(boxtt$Clust))]), linetype = 2)+
  stat_compare_means(method = "anova", label.y = -1)+
  scale_fill_brewer(palette="YlGnBu")+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")

ggsave(file.path(fig4.path,paste0("Clust signature subgroup.pdf")), width = 6, height = 4)

boxtt<-TARGETmulti[,c("signature","Histologic")]

boxtt<-boxtt %>% 
  drop_na()

boxtt<-boxtt[order(boxtt$Histologic, decreasing = F), ]


compare_means(signature ~ Histologic,  data = boxtt,
              ref.group = ".all.", method = "t.test")

# 画小提琴图 + 箱线图 + jitter
ggplot(boxtt, aes(x = Histologic, y = signature, fill = Histologic)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  geom_hline(yintercept = mean(boxtt$signature[!is.na(boxtt$Histologic)]), 
             linetype = 2, color = "black") +
  stat_compare_means(method = "t.test", label.y = 0.8) +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.") +
  scale_fill_manual(values = c( "#FDE820","#0B4FA2")) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig4.path,paste0("Histologic signature subgroup.pdf")), width = 4, height = 4)


#------森力图后multiROC
library("survival")
library("survminer")
library("riskRegression")
library(dplyr)

write.table(TARGETmulti,file=file.path(fig6.path,"TARGETmulti.txt"),sep="\t",row.names=T,quote=F)

TARGETmulti<-cbind(TARGET.clin,signature)
all_expdata<-cbind(TARGET.clin,signature)%>% 
  select(c(signature))

all_expdata<-t(all_expdata)

mark_expdata<-all_expdata

clindata<-dplyr::select(TARGETmulti,c("Stage","Gender","Age","Histologic"))
traitData<-dplyr::select(TARGETmulti,c("OS.time","OS"))

expdata = all_expdata[rownames(mark_expdata),]
rownames(clindata) = gsub("-", ".", rownames(clindata))
rownames(traitData) = gsub("-", ".", rownames(traitData))
colnames(traitData) = c("time", "status")
clindata2 = clindata
traitData2 = traitData
exp_sig = expdata

#Nomogram
predict_mat_all = na.omit(cbind(traitData2, clindata2, exp_sig))
cox_form = as.formula(paste("Surv(time, status)", paste("`", paste(colnames(predict_mat_all)[3:ncol(predict_mat_all)], collapse = "` + `"), "`", sep = ""), sep = " ~ "))
res.cox.all = coxph(cox_form, data = predict_mat_all, x = TRUE)
res.cox.lst = list(res.cox.all)

#Classifier
predict_mat_classifier = cbind(traitData2, exp_sig)[rownames(predict_mat_all),]
cox_form = as.formula(paste("Surv(time, status)", paste("`", paste(colnames(predict_mat_classifier)[3:ncol(predict_mat_classifier)], collapse = "` + `"), "`", sep = ""), sep = " ~ "))
res.cox.classifier = coxph(cox_form, data = predict_mat_classifier, x = TRUE)
res.cox.lst = c(res.cox.lst, list(res.cox.classifier))

# Single variate Cox regression ROC comparation
predict_mat_tmp = cbind(traitData2, clindata2)[rownames(predict_mat_all),]
predict_mat_tmp$Nomogram = predict(res.cox.all,predict_mat_all)
predict_mat_tmp$Classifier = predict(res.cox.classifier,predict_mat_classifier)
predict_mat_for_single = predict_mat_tmp[,c(1:2,ncol(predict_mat_tmp)-1,ncol(predict_mat_tmp),3:(ncol(predict_mat_tmp)-2))]
res.cox.lst = NULL
for (i in 3:(ncol(predict_mat_for_single))){
  test_mat_single = predict_mat_for_single[,c(1,2,i)]
  cox_form = as.formula(paste("Surv(time, status)", paste("`", colnames(predict_mat_for_single)[i], "`", sep = ""), sep = " ~ "))
  cox.res.test = coxph(cox_form, data = predict_mat_for_single, x = TRUE)
  if(is.null(res.cox.lst)){
    res.cox.lst = list(cox.res.test)
  }else{
    res.cox.lst = c(res.cox.lst, list(cox.res.test))
  }
}
names(res.cox.lst) = c(colnames(predict_mat_for_single)[3:ncol(predict_mat_for_single)])
xs <- Score(res.cox.lst, Hist(time,status)~1,data=predict_mat_for_single, plots="roc",metrics="auc")
pdf(file=file.path(fig4.path,"/total TARGET ROC combined.pdf", sep = ""), height = 6.5, width = 6)
plotROC(xs, xlab = "False negative rate", ylab = "Ture negative rate", legend = TRUE, auc.in.legend = TRUE,col=brewer.pal(6,"Set1"))
dev.off()

library(RColorBrewer)
n <- length(xs$models)  # 或 length(res.cox.lst)，总共几条ROC曲线
colors <- brewer.pal(n, "Dark2")



##################################
############Figure 5 Nomogram#####
##################################

#########nomogram

library(survival)

pbc<-TARGETmulti[,c("OS.time","OS","Age","Stage","Gender","signature")]
pbc[pbc=="[Discrepancy]"]<-NA
pbc[pbc==""]<-NA
pbc<-pbc %>% 
  drop_na()

pbc$Stage<-ifelse(pbc$Stage=="I","I+II",
                  ifelse(pbc$Stage=="II","I+II","III+IV"))
pbc$Age<-pbc$Age/12

pbccox <- coxph(formula = Surv(OS.time,OS) ~ ., data = pbc)
pbccox

library(regplot)

regplot(pbccox,
        #对观测2的六个指标在列线图上进行计分展示
        observation=pbc[c("TCGA-VS-A9UJ-01"),], #也可以不展示
        points = T, #If FALSE the regression scores of each βx contribution are shown. Otherwise contributions are represented by a 0-100 "points" scale.
        plots = c("bean", #可选"no plot" "density" "boxes" "ecdf" "bars" "boxplot" "violin" "bean" "spikes"
                  "bars"), #可选"no plot" "boxes" "bars" "spikes"
        subticks = TRUE, dencol="#BA6EAE",boxcol="#00B5B7",obscol="red",spkcol="brown",
        #预测3年和5年的死亡风险，此处单位是day
        center = T,
        failtime = c(36,12),
        prfail = TRUE, #cox回归中需要TRUE
        rank="sd", #rank="range" is by the range of the βx's, and rank="sd" is by the standard deviation of the βx's. 
        showP = T, #是否展示统计学差异
        droplines = T,#观测2示例计分是否画线
        #colors = mycol, #用前面自己定义的颜色
        #rank=NULL, #根据统计学差异的显著性进行变量的排序
        #interval="confidence"
) #展示观测的可信区间

dev.copy2pdf(file=file.path(fig5.path, "nomogram_new.pdf"), width = 9,height = 6)

#------------K-M

library(survival)
library("survminer")
points<-pbccox$linear.predictors
KM_input<-cbind(pbc,points)

KM_input$points<- factor(ifelse(KM_input$points>median(KM_input$points),"high"," low"))
outTab=data.frame()
fit <- survfit(Surv(OS.time, OS) ~ points, data = KM_input)

cox <- coxph(Surv(OS.time, OS) ~ points, data = KM_input)
coxSummary = summary(cox)
coxP=coxSummary$coefficients[,"Pr(>|z|)"]
outTab=rbind(outTab,
             cbind(
               HR=coxSummary$conf.int[,"exp(coef)"],
               HR.95L=coxSummary$conf.int[,"lower .95"],
               HR.95H=coxSummary$conf.int[,"upper .95"],
               pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
)
outTab

HR <- paste("Hazard Ratio = ", round(outTab$HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(outTab$HR.95L,3), round(outTab$HR.95H,3), sep = " - "), sep = "")


pdf(file=file.path(fig5.path,"nomogram survival.pdf"), width=4,height=4.2,onefile = FALSE)
ggsurvplot(fit, data = KM_input,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), 
           #fun = "event",
           conf.int = T, 
           #conf.int.style = "",
           censor = T, 
           palette = c( "#57A7B9","#F05C52"), #
           ylim = c(0,1),
           
           xlab = 'Time in months',
           legend.title='riskscore', 
           legend.labs=c('Low points','High points'), 
           
           risk.table.y.text.col = T, 
           risk.table.y.text = T, 
           
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"),
           
           
           
           
           pval = paste(pval = ifelse(outTab$pvalue < 0.001, "p < 0.001", 
                                      paste("P = ",round(outTab$pvalue,3), sep = "")),
                        HR, CI, sep = "\n"),
           pval.coord = c(20, 0.2)
)

dev.off()

outTab

#----------ROC————————————

library(timeROC)
library(survival)
KM_input1<-cbind(pbc,points)
head(KM_input1)

pdf(file.path(fig5.path,"nomograme ROC2.pdf"), 5, 5.5)
ROC.DSST<-timeROC(T=KM_input1$OS.time,#结局时间
                  delta=KM_input1$OS,#生存结局
                  marker=KM_input1$points,#预测变量
                  cause=1,#阳性结局赋值，比如死亡，复发的赋值
                  weighting="marginal",# 权重计算方法，marginal是默认值，采用km计算删失分布
                  times=c(12,36,30),# 时间点，选取10年和20年生存率
                  ROC = TRUE,
                  iid = TRUE
)
plot(ROC.DSST,time=12,col=jco3[1],title=FALSE,lwd=2)
plot(ROC.DSST,time=36,col=jco3[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC.DSST,time=30,col=jco3[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC.DSST$AUC[1],3)),
         paste0('AUC at 3 years: ',round(ROC.DSST$AUC[2],3)),
         paste0('AUC at 5 years: ',round(ROC.DSST$AUC[3],3))),
       col=jco3,lwd=2,bty = 'n')
dev.off()

#------------------Calibration---------------
#------------------Calibration---------------

library(ResourceSelection)
pbccox <- coxph(formula = Surv(OS.time,OS) ~Age+Stage+Gender+signature, data = pbc)
library(survival)
library(regplot)
library(rms)

f1<-cph(formula = Surv(OS.time,OS) ~Age+Stage+Gender+signature, data = pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 12) 
#参数m=50表示每组50个样本进行重复计算
cal1<-calibrate(f1, cmethod="KM", method="boot",u=12,m=30,B=1000) 

pbc$OS1<-ifelse(pbc$OS.time>12,"0",pbc$OS)
pbc$OS1<-as.numeric(pbc$OS1)

dat1 <- as.data.frame(cbind(pbc$OS1,f1$linear.predictors))
colnames(dat1)<-c("OS","model")

fullmodel_glm1 <- glm(OS ~ model, 
                      data = dat1, 
                      family = "binomial", 
                      control = list(maxit = 50))

p.hoslem1 <- hoslem.test(fullmodel_glm1$y, fitted(fullmodel_glm1), g=10)$p.value


f3<-cph(formula = Surv(OS.time,OS) ~Age+Stage+Gender+signature, data = pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 36) 
#参数m=50表示每组50个样本进行重复计算
cal3<-calibrate(f3, cmethod="KM", method="boot",u=36,m=30,B=1000) 

pbc$OS3<-ifelse(pbc$OS.time>36,"0",pbc$OS)
pbc$OS3<-as.numeric(pbc$OS3)

dat3 <- as.data.frame(cbind(pbc$OS3,f3$linear.predictors))
colnames(dat3)<-c("OS","model")

fullmodel_glm3 <- glm(OS ~ model, 
                      data = dat3, 
                      family = "binomial", 
                      control = list(maxit = 50))

p.hoslem3 <- hoslem.test(fullmodel_glm3$y, fitted(fullmodel_glm3), g=10)$p.value

f5<-cph(formula = Surv(OS.time,OS) ~Age+Stage+Gender+signature, data = pbc,x=T,y=T,surv = T,na.action=na.delete,time.inc = 60) 
#参数m=50表示每组50个样本进行重复计算
cal5<-calibrate(f5, cmethod="KM", method="boot",u=60,m=30,B=1000) 

pbc$OS5<-ifelse(pbc$OS.time>60,"0",pbc$OS)
pbc$OS5<-as.numeric(pbc$OS5)

dat5 <- as.data.frame(cbind(pbc$OS5,f5$linear.predictors))
colnames(dat5)<-c("OS","model")

fullmodel_glm5 <- glm(OS ~ model, 
                      data = dat5, 
                      family = "binomial", 
                      control = list(maxit = 50))

p.hoslem5 <- hoslem.test(fullmodel_glm5$y, fitted(fullmodel_glm5), g=10)$p.value


pdf(file.path(fig5.path,"calibration_compare.pdf"),width = 6,height = 6)
plot(cal1,lwd = 2,lty = 0,errbar.col = jco2[1],
     bty = "o", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = jco2[1],
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = jco2[1], pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 0,errbar.col = jco2[2],
     xlim = c(0,1),ylim= c(0,1),col = jco2[2],add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = jco2[2], pch = 16)

plot(cal5,lwd = 2,lty = 0,errbar.col = jco2[3],
     xlim = c(0,1),ylim= c(0,1),col = jco2[3],add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = jco2[3], pch = 16)

abline(0,1, lwd = 1, lty = 3, col = c("lightgrey"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =jco2[1:3], #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框

text(0.25,0.205,bquote("Hosmer-Lemeshow 1-year"~italic(P)~" = "~.(round(p.hoslem1,3))),adj = 0)
text(0.25,0.135,bquote("Hosmer-Lemeshow 3-year"~italic(P)~" = "~.(round(p.hoslem3,3))),adj = 0)
text(0.25,0.065,bquote("Hosmer-Lemeshow 5-year"~italic(P)~" = "~.(round(p.hoslem5,3))),adj = 0)

dev.off()


pdf(file=file.path(fig5.path,"calibration_compare 1-year.pdf"),width = 5,height = 5)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#14B96A"),
     bty = "o", #只画左边和下边框
     xlim = c(0.6,1),ylim= c(0.6,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#14B96A"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#14B96A"), pch = 16)
mtext("")


abline(0,1, lwd = 1, lty = 3, col = c("lightgrey"))

text(0.60,0.63,bquote("Hosmer-Lemeshow 1-year"~italic(P)~" = "~.(round(p.hoslem1,3))),adj = 0)

dev.off()


pdf(file=file.path(fig5.path,"calibration_compare 3-year.pdf"),width = 5,height = 5)
plot(cal3,lwd = 2,lty = 0,errbar.col = c("#FFA900"),
     bty = "o", #只画左边和下边框
     xlim = c(0.35,1),ylim= c(0.35,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#FFA900"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#FFA900"), pch = 16)
mtext("")


abline(0,1, lwd = 1, lty = 3, col = c("lightgrey"))

text(0.45,0.40,bquote("Hosmer-Lemeshow 3-year"~italic(P)~" = "~.(round(p.hoslem3,3))),adj = 0)

dev.off()

pdf(file=file.path(fig5.path,"calibration_compare 5-year.pdf"),width = 5,height = 5)
plot(cal5,lwd = 2,lty = 0,errbar.col = c("#5C0A98"),
     bty = "o", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#5C0A98"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#5C0A98"), pch = 16)
mtext("")


abline(0,1, lwd = 1, lty = 3, col = c("lightgrey"))

text(0.15,0.12,bquote("Hosmer-Lemeshow 5-year"~italic(P)~" = "~.(round(p.hoslem5,3))),adj = 0)

dev.off()

##-------------------

for (i in 3:5) {
  fit <- coxph(Surv(lenfol, fstat)~rcspline.eval(whas500$bmi,nk=i,inclx = T)+gender,data=whas500, x=TRUE)
  tmp <- extractAIC(fit)
  if(i == 3) {AIC = tmp[2]; nk = 3}
  if(tmp[2] < AIC) {AIC = tmp[2]; nk = i} #nk保存最优的knots数目，且样条数目为nk-2。具体见参考文献“理论1.pdf"第2页2.3节公式1
}


#-------------DCA-------------------


library(rms)
library(rmda)

modul<- decision_curve(OS ~ Age+Stage+Gender+signature,
                       data= pbc,
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals = 0.95)
modul1<- decision_curve(OS ~ Age,
                        data= pbc,
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals = 0.95)
modul2<-decision_curve(OS ~ Stage,
                       data= pbc,
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals = 0.95)
modul3<- decision_curve(OS ~ Gender,
                        data= pbc,
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals = 0.95)

pbccox <- glm(OS ~ ., data = pbc)

list<-list(modul,modul1,modul2,modul3)
jcot <- c("#411050","#97B959","#C72E24","#0088A2","#F16A2F","#80A7DE")

pdf(file=file.path(fig5.path,"DCA.pdf"),width = 4,height = 4.5)
plot_decision_curve(list,
                    curve.names=c("Nomogram","Age","Stage","Tumor_status"),
                    xlab="Threshold probability",
                    cost.benefit.axis =FALSE,col= jcot,
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()

modul


#-------------------clinical impact curve nomogram--------------

pdf(file=file.path(fig5.path,"Clinical_impact_curve.pdf"),width = 4,height = 4.5)
plot_clinical_impact(modul,xlim = c(0, 1),legend.position = "topright",
                     col = c("#411050", "#97B959"))
dev.off()
#--------------------points correlation

KM_input1$points2<-(54.58+KM_input1$signature*41.25)+
  ifelse(KM_input1$Gender=="Female",50,69)+
  ifelse(KM_input1$Stage=="I+II",50,78)+
  (-2.29*KM_input1$Age+60)

library(ggstatsplot)
ggscatterstats(
  data = KM_input1,
  x = points2,
  y = OS.time,
  bf.message = FALSE
)
ggsave(file=file.path(fig5.path,"Correlation points OS time.pdf"),width = 5,height = 5)

KM_input1$OS2<-as.factor(KM_input1$OS)
ggscatter(KM_input1, x = "points2", y = "OS.time",
          add = "reg.line",               # Add regression line
          conf.int = TRUE,                # Add confidence interval
          color = "OS2", palette = "Set1", # Color by groups "cyl"
          #shape = "OS"                   # Change point shape by groups "cyl"
)+
  stat_cor(aes(color = OS2), label.x = 180)       # Add correlation coefficient
ggsave(file=file.path(fig5.path,"Correlation points OS time subgroup.pdf"),width = 5,height = 5)

KM_input1$OS2<-ifelse(KM_input1$OS2==1,"Death","Alive")

p <- ggboxplot(KM_input1, x = "OS2", y = "points2",notch = T,
                fill = "OS2",alpha = 0.6,
               add = "jitter")+
  scale_fill_brewer(palette = "Set1", direction = -1)

p + stat_compare_means(method = "t.test")+
  labs(x = '',
       y = 'Nomogram Points')

ggsave(file=file.path(fig5.path,"points subgroup.pdf"),width = 5,height = 5)

##计算预测C-index 并绘图

library(dynpred)
cindex.sig<-cindex(Surv(OS.time,OS) ~Age+Stage+Gender+signature, data = pbc)
cindex.age<-cindex(Surv(OS.time,OS) ~Age, data = pbc)
cindex.stage<-cindex(Surv(OS.time,OS) ~Stage, data = pbc)
cindex.gender<-cindex(Surv(OS.time,OS) ~Gender, data = pbc)

df <- data.frame(dose=c(" Nomogram", "Age", "Stage","Gender"),
                 len=c(round(cindex.sig$cindex,3), round(cindex.age$cindex,3),
                       round(cindex.stage$cindex,3),round(cindex.gender$cindex,3)))

# Change barplot fill colors by groups
p<-ggplot(data=df, aes(x=dose, y=len,fill=dose)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=len), vjust=-0.3, size=4,color = "black")+
  theme_bw()+
  ylim(0,0.8)+
  labs(x = '',
       y = 'C-index')+
  scale_fill_brewer(palette="Set1")+ 
  theme(
    legend.position = "none",
    
    # 设置主标题字体（font.main = c(14, "bold", "darkblue")）
    plot.title = element_text(size = 12, color = "darkblue", hjust = 0.5),
    
    # 设置x轴标题和刻度字体（font.x = c(14, "bold", "black")）
    axis.title.x = element_text(size = 12, color = "black"),
    axis.text.x  = element_text(size = 12, color = "black"),
    
    # 设置y轴标题和刻度字体（font.y = c(14, "bold", "black")）
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.y  = element_text(size = 12, color = "black")
  )
p2<-p+scale_fill_brewer(palette="Dark2")+ theme(legend.position="none")

p2
ggsave(file=file.path(fig5.path,"C-index nomogram others.pdf"),width = 4,height = 4)


###############singlecell part##################

# Seurat v4.4
library(Seurat)
library(Matrix)
library(data.table)

library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

expr_dt <- fread("GSM5344367_WTTOT.txt", sep = "\t", header = TRUE, check.names = FALSE)
cell_ids <- expr_dt[[1]]
expr_dt <-expr_dt[,-c(1:11)]
expr <- as.matrix(expr_dt)
rownames(expr) <- cell_ids   





mat <- t(expr)
mat <- Matrix(mat, sparse = TRUE)
mat[!is.finite(mat)] <- 0                  
mat <- mat[Matrix::rowSums(mat > 0) > 0, ] 
colnames(mat) <- make.names(colnames(mat)) 
rownames(mat) <- make.names(rownames(mat)) 
obj <- CreateSeuratObject(counts = mat, assay = "RNA", project = "wlims")
saveRDS(obj,'obj.rds')


# ==== 人类 QC 指标 ====
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")   # 可选
obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]") # 可选

# ==== 自适应阈值（MAD），同时设经验上限 ====
set.seed(10086)
mad_cut <- function(x, k=3) { m <- median(x); c(m - k*mad(x), m + k*mad(x)) }
nf_lim <- mad_cut(obj$nFeature_RNA); nf_lim[1] <- max(nf_lim[1], 200); nf_lim[2] <- min(nf_lim[2], 6000)
mt_max <- min(mad_cut(obj$percent.mt)[2], 20)

obj <- subset(obj, subset = nFeature_RNA >= nf_lim[1] &
                nFeature_RNA <= nf_lim[2] &
                percent.mt   <= mt_max)


# ==== 聚类配色 ====
cluster_cols <- c("#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72",
                  "#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
                  "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D",
                  "#E6AB02","#7570B3","#BEAED4","#666666","#999999",
                  "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000",
                  "#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00")
# ==== 标准流程（非 SCT） ====
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)
set.seed(123)
obj <- RunTSNE(
  obj,
  dims = 1:30,          # 视你的数据/PC拐点而定
  reduction = "pca"     # 用 PCA 作为输入
)

# ==== 绘图（聚类配色） ====
n.clust <- length(levels(Idents(obj)))
stopifnot(n.clust <= length(cluster_cols))
p_umap <- DimPlot(obj, reduction = "tsne", label = TRUE,
                  cols = cluster_cols[seq_len(n.clust)])
FeaturePlot(obj,features = 'CD3D')
# 可选 QC 图
p_qc <- VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                pt.size = 0, ncol = 3)

# ==== 保存 ====
ggsave(paste0("UMAP_clusters.pdf"), p_qc, width = 12, height = 4)
saveRDS(obj, paste0(SAMPLE_NAME, "_seurat_v44_human_noSCT.rds"))


#计算细胞周期，seurat包自带了
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, 
                        s.features = s.genes, 
                        g2m.features = g2m.genes,
                        set.ident = TRUE)
VlnPlot(obj,features = c("S.Score","G2M.Score"))



#如果你还想要其他的打分
#aucell
#addumodulescore
library(Seurat)
library(AUCell)

stem_genes <- c(
  "SOX2","POU5F1","NANOG","KLF4","MYC","LIN28A","LIN28B","ESRRB","PRDM14",
  "DPPA4","ZFP42","UTF1","TERT","PROM1","ALDH1A1"
)
expr <- GetAssayData(obj, assay = "RNA", slot = "data")  # gene x cell

# AUCell 需要 gene×cell 的矩阵
cells_rankings <- AUCell_buildRankings(expr, nCores = 1, plotStats = FALSE, verbose = FALSE)

geneSets <- list(Stemness = stem_genes)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# === 3) 写回 Seurat meta.data ===
obj$Stemness_AUC <- as.numeric(getAUC(cells_AUC)["Stemness", colnames(obj)])
head(obj)
# 可视化
VlnPlot(obj, features = "Stemness_AUC", pt.size = 0) +
  ggtitle("Stemness (AUCell)")
FeaturePlot(obj, features = "Stemness_AUC",order = T)



#########计算不同区分度clustree
obj <- FindClusters(obj, 
                    resolution = seq(from = 0.1, 
                                     to = 1.0, 
                                     by = 0.2))

library(clustree)

clustree(obj)
clustree(obj@meta.data, prefix = "RNA_snn_res.")
ggsave("clustree.pdf", width = 8, height = 6)

########根据PCA计算 umap、tsne
obj <- RunUMAP(obj, reduction = "pca", dims = 1:10)
obj <- RunTSNE(obj, reduction = "pca", dims = 1:10)

p_umap1 <- DimPlot(obj, reduction = "umap", label = TRUE
                   # cols = cluster_cols[seq_len(n.clust)]
)
p_umap2 <- DimPlot(obj, reduction = "tsne", label = TRUE
                   # cols = cluster_cols[seq_len(n.clust)]
)

p_umap1+p_umap2

# 想画的分辨率
res_vec   <- c("0.1","0.3","0.5","0.7","0.9")
meta_cols <- paste0("RNA_snn_res.", res_vec)

# 保险起见：检查这些列都在 meta.data 中
stopifnot(all(meta_cols %in% colnames(obj@meta.data)))

# 如果没有 tSNE，则计算一个（基于 pca 或当前默认降维）
if (!"tsne" %in% Reductions(obj)) {
  base_red <- if ("pca" %in% Reductions(obj)) "pca" else DefaultReduction(obj)
  set.seed(123)
  obj <- RunTSNE(obj, dims = 1:30, reduction = base_red, check_duplicates = FALSE)
}

# 颜色：给足够多的离散色
my_cols<-cols_all <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
  "#8C564B", "#17BECF", "#BCBD22", "#9467BD", "#FF9896",
  "#C5B0D5", "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5",
  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
  "#D6616B", "#E7969C", "#DE9ED6", "#9C9EDE", "#6BAED6",
  "#9ECAE1", "#C6DBEF", "#FDD0A2", "#FDAE6B", "#F16913"
)

# 一个小工具函数：在指定降维上画某个分辨率
plot_one <- function(seu, red, meta_col) {
  DimPlot(
    seu, reduction = red, group.by = meta_col,
    label = TRUE, repel = TRUE, cols = cols_all
  ) + ggtitle(paste0(toupper(red), " | res = ", sub("RNA_snn_res\\.", "", meta_col))) +
    theme(plot.title = element_text(size = 12, face = "bold"))
}

# 逐分辨率作图：TSNE 与 UMAP
plt_tsne <- lapply(meta_cols, function(m) plot_one(obj, "tsne", m))
plt_umap <- lapply(meta_cols, function(m) plot_one(obj, "umap", m))

# 拼图（每行2列，可按需调整）
p_tsne_grid <- wrap_plots(plt_tsne, ncol = 3)
p_umap_grid <- wrap_plots(plt_umap, ncol = 3)

# 展示
p_tsne_grid
p_umap_grid

# 如需保存到文件（可选）
ggsave("tsne_all_res.pdf", p_tsne_grid, width = 10, height = 6)
# ggsave("umap_all_res.png", p_umap_grid, width = 10, height = 12, dpi = 300)


# 根据你的数据类型调整：这里给一套免疫常用面板
DotPlot_grouped <- function(obj, panel, group.by = "seurat_clusters",
                            alpha_bg = 0.25,
                            palette = "Set3") {
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(ggh4x)
  require(RColorBrewer)
  
  # 1) 去重，保证一个基因只出现一次
  panel <- lapply(panel, unique)
  .seen <- character(0)
  panel <- lapply(panel, function(g) { keep <- setdiff(g, .seen); .seen <<- c(.seen, keep); keep })
  
  # 2) 展开 features 与分组名
  feat      <- unlist(panel, use.names = FALSE)
  grp_names <- rep(names(panel), lengths(panel))
  
  # 3) 基础 DotPlot
  p <- DotPlot(obj, features = feat, group.by = group.by) + RotatedAxis()
  
  # 4) 加入 gene_groups 并固定顺序
  map_df <- data.frame(features.plot = feat, gene_groups = grp_names)
  p$data <- left_join(p$data, map_df, by = "features.plot")
  facet_levels <- names(panel)
  p$data$gene_groups <- factor(p$data$gene_groups, levels = facet_levels)
  
  # 5) 统一配色（上下同色）
  base_cols <- setNames(
    colorRampPalette(brewer.pal(12, palette))(length(facet_levels)),
    facet_levels
  )
  bg_cols <- sapply(facet_levels,
                    function(g) grDevices::adjustcolor(base_cols[g], alpha.f = alpha_bg))
  names(bg_cols) <- facet_levels
  
  # 6) strip 主题
  strip_order <- levels(p$data$gene_groups)
  strip_theme <- ggh4x::strip_themed(
    background_x = lapply(strip_order,
                          function(g) element_rect(fill = bg_cols[g], color = NA)),
    text_x       = lapply(strip_order,
                          function(g) element_text(face = "bold", size = 10))
  )
  
  # 7) 每个分面背景矩形
  bg_df <- data.frame(
    gene_groups = factor(strip_order, levels = strip_order),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
  
  # 8) 拼接图层
  p2 <- p +
    geom_rect(
      data = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = gene_groups),
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = bg_cols, guide = "none") +
    ggh4x::facet_nested(~ gene_groups,
                        scales = "free_x", space = "free_x", strip = strip_theme) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.spacing.x  = unit(0.6, "lines"),
      strip.background = element_blank(),
      strip.placement  = "outside",
      plot.margin      = margin(10, 30, 10, 10)
    ) +
    coord_cartesian(clip = "off")
  
  return(p2)
}

panel <- list(
  Epithelial  = c("EPCAM", "KRT17", "KRT8", "KRT18"),
  Fibro       = c("DCN", "COL1A2", "COL1A1",'ACTA2'),
  Endo        =c("PECAM1", "VWF", "RAMP2"),
  Plasma      = c("MZB1", "JCHAIN", "IGHG1"),
  B           = c("CD79A", "MS4A1"),
  CD4CD8      = c("CD3D", "CD3E", "CD4", "CD8A"),
  Mono        = c("CD14", "S100A8", "S100A9"),
  Macrophage  = c("C1QA", "C1QB", "CD68", "CD86"),
  Myeloid     = c("CD74", "CST3"),
  Dendritic   = c("ITGAX", "BST2", "CLEC9A", "SIGLEC1")  ,
  Prof        = c("MKI67", "STMN1", "PCNA"),
  NK          = c("NKG7", "KLRD1",'GNLY','KLRF1'),
  Naive       = c("TCF7", "SELL", "LEF1"),
  Neutrophils = c("CSF3R", "CXCL3"),
  Pericyte   =c("RGS5", "CSPG4", "ABCC9", "KCNJ8"),
  Mast  = c("TPSAB1", "TPSB2", "MS4A2", "CPA3")
)

# 直接调用函数
DotPlot_grouped(obj, panel, group.by = "seurat_clusters")

# 如需保存到文件（可选）
ggsave("group identify dotplot.pdf",  width = 16, height = 8)


#############画图展示整理后的分簇

# 建立映射关系：0,1,2 -> 0; 3 -> 1; 4 -> 2; 5 -> 3
mapping <- c("0" = "0", "1" = "0", "2" = "0","3" = "0","4" = "0","5" = "0",
             "6" = "1", "7" = "2", "8" = "3")

# 应用映射到 RNA_snn_res.0.5
new_clusters <- mapping[ as.character(obj@meta.data[["RNA_snn_res.0.9"]]) ]

# 写入 seurat_clusters
obj@meta.data[["seurat_clusters"]] <- factor(new_clusters,
                                             levels = c("0","1","2","3"))

# 检查结果
table(obj@meta.data[["seurat_clusters"]])

Idents(obj) <- "seurat_clusters"

new.cluster.ids <- c("0"="Epi", 
                     "1"="T_NK", 
                     "2"="Myeloid", 
                     "3"="Prof")
obj <- RenameIdents(obj, new.cluster.ids)                        
obj$celltype <- obj@active.ident
DimPlot(obj, label = T,reduction = 'tsne',pt.size = 1,group.by = "celltype")

library("scplotter")

p1<-CellDimPlot(obj, group_by = "RNA_snn_res.0.9", reduction = "tsne", theme= "theme_blank")+
  scale_color_manual(values=my_cols)

p2<-CellDimPlot(obj, group_by = "celltype", reduction = "tsne", theme= "theme_blank")+
  scale_color_manual(values=my_cols)

# 不用 collect（容易触发 bug），先用 keep
p1 + p2 

ggsave("CellDimPlot for groups.pdf",  width = 10, height = 5)

#--------------------------------------#
#########识别每组差异基因，画图#########
#--------------------------------------#

all_genes <- rownames(obj)
genes_to_keep <- all_genes[!grepl("^MT-|^RP", all_genes)]
uterus <- subset(obj, features = genes_to_keep)
all.markers  <- FindAllMarkers(uterus, 
                               only.pos = FALSE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.75)
significant.markers  <- all.markers [all.markers $p_val_adj < 0.2, ]
write.csv(significant.markers, file = "significant.markers.csv")

myColor_sorted <- c(
  "#D87F32", "#48A75A", "#E3BC06", "#E41B1B", "#FA9B93", "#E9358B",
  "#A0094E", "#BD5E95", "#D690C6", "#B17A7D", "#87638F",
  "#737690", "#999999", "#204B75", "#4376AC", "#4285BF",
  "#6FCDDC", "#B6DB7B", "#588257", "#847A74"
)

library(scRNAtoolVis)
markerVolcano(markers = all.markers,
              topn = 5,
              labelCol = my_cols)

ggsave("DEGs for groups.pdf",  width = 10, height = 6)


DEGs       = c("TCIM","KLRC3","C1QC","TNNT3")
FeaturePlot(obj,features = DEGs,order = T,reduction = "tsne",cols = c("lightgrey", "#E41B1B"))

ggsave("DEGs ditribution expression.pdf",  width = 7, height = 6)

#--------------------------------------#
#---画出细胞周期分数在不同细胞的表达---#
#--------------------------------------#

# 画出细胞周期分布
# Seurat 内置的 cell cycle 基因集
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

# 打分 & 写入 obj@meta.data$cellcycle
obj <- CellCycleScoring(obj,
                        s.features = s.genes,
                        g2m.features = g2m.genes,
                        set.ident = FALSE)

# 查看新加的列
head(obj@meta.data[, c("S.Score", "G2M.Score", "Phase")])

p <- DimPlot(obj, reduction = "tsne", group.by = "Phase")

# ggplot 对象的数据在 p$data
p$data <- p$data[order(p$data$Phase, decreasing = TRUE), ]  # 顺序可调

p + scale_color_manual(values = c("G1"  = "#1f77b4",
                                  "G2M" = "#d62728",
                                  "S"   = "#2ca02c"))

ggsave("cellcycle phase ditribution expression.pdf",  width = 3.5, height = 3) 


#---画出细胞周期评分小提琴图---#
#--------------------------------------#

# 1) 指定你要看的四种细胞类型（名字必须与 obj$celltype 的水平一致）
#cells_of_interest <- c("Epithelial", "Fibro", "Immune", "Endothelial")  # 示例

# 2) 取数：分簇、细胞类型、S/G2M 打分
df <- FetchData(obj, vars = c("seurat_clusters", "celltype", "S.Score", "G2M.Score")) %>%
  mutate(
    seurat_clusters = factor(seurat_clusters),       # 确保按因子处理
    celltype = factor(celltype)
  ) %>%
  droplevels() %>%
  pivot_longer(cols = c("S.Score", "G2M.Score"),
               names_to = "Score", values_to = "value")

## 4) 计算 Epithelial 组在每个 Score 下的均值（用于画横线）
baseline_group <- "Epi"   # ← 这里改成你想作为基准的组名

# 计算每个面板中 Epi 的均值
epi_means <- df %>%
  dplyr::filter(celltype == baseline_group) %>%
  dplyr::group_by(Score) %>%
  dplyr::summarise(yint = mean(value, na.rm = TRUE), .groups = "drop")

## 5) 画图：并排小提琴 + 中位数点 + k检验 + Epi均值虚线
p <- ggplot(df, aes(x = celltype, y = value, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width",
              width = 0.9, color = NA, alpha = 0.85) +
  # 中位数点
  stat_summary(fun = mean, geom = "point",
               position = position_dodge(width = 0.9),
               size = 0.8, color = "black") +
  # k检验（整体比较）
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.y = 0.8) +
  facet_wrap(~ Score, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = my_cols, guide = "none") +
  labs(x = "细胞类型", y = "Score") +
  theme_classic(base_size = 12) +
  theme(
    strip.text  = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 画图（p 是你的小提琴图主体）
p +
  geom_hline(
    data = epi_means,
    aes(yintercept = yint),
    linetype = 2, linewidth = 0.5,
    inherit.aes = FALSE
  )

ggsave("cellcycle score in different cells.pdf",  width = 6, height = 7) 


#--------------------------------------#
#---画出一组基因在tsne图的表达---#
#--------------------------------------#

featureplot_to_pdf <- function(obj, genes, file = "cellcycle_featureplot.pdf",
                               nrow = 3, ncol = 4,
                               reduction = "tsne",
                               cols = c("lightgrey", "#E41B1B")) {
  stopifnot(requireNamespace("patchwork", quietly = TRUE))
  library(ggplot2)
  
  # 只保留对象中存在的基因，并提示缺失
  present <- intersect(genes, rownames(obj))
  missing <- setdiff(genes, rownames(obj))
  if (length(missing) > 0) {
    message("These genes are not in the object and will be skipped: ",
            paste(missing, collapse = ", "))
  }
  if (length(present) == 0) stop("No valid genes to plot.")
  
  # 打开 PDF 设备；确保出错也能关闭
  grDevices::pdf(file, width = 10, height = 7)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  # 分页，每页 nrow*ncol 个图
  page_size <- nrow * ncol
  idx <- seq(1, length(present), by = page_size)
  for (i in idx) {
    subset_genes <- present[i:min(i + page_size - 1, length(present))]
    
    plots <- lapply(subset_genes, function(g) {
      p <- FeaturePlot(obj, features = g, order = TRUE,
                       reduction = reduction, cols = cols) +
        ggtitle(g) +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              legend.position = "none")
      # 避免某些版本返回列表
      if (inherits(p, "list")) p <- p[[1]]
      p
    })
    
    # 如果最后一页不足，填充空白占位，避免 patchwork 报错
    if (length(plots) < page_size) {
      n_blank <- page_size - length(plots)
      blanks <- replicate(n_blank, ggplot() + theme_void(), simplify = FALSE)
      plots <- c(plots, blanks)
    }
    
    pg <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow)
    print(pg)  # 关键：把拼好的对象“打印”到当前 pdf 设备
  }
  
  invisible(NULL)
}


cellcyclegene <- c(
  "AURKA","AURKB",
  "CDC6","CDC45","CDC25A","CDC25B","CDC25C","FZR1",
  "CDK1","CDK2","CDK3","CDK4","CDK5","CDK6","CDK7","CDK9",
  "CDT1","CHEK1","CHEK2","MYC",
  "CCNA2","CCNB1","CCNB2","CCND1","CCND2","CCND3","CCNE1","CCNE2","CCNF","CCNH",
  "E2F1","FOXM1","GMNN","MKI67",
  "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7",
  "PKMYT1",
  "ORC1","ORC2","ORC3","ORC4","ORC5","ORC6",
  "CDKN2B","CDKN2A","CDKN2C","CDKN2D","CDKN1A","CDKN1B","CDKN1C",
  "TP53","PCNA","RB1","RBL1","RBL2","SKP2","WEE1","PLK1"
)


featureplot_to_pdf(
  obj,
  genes = cellcyclegene,          # 你的基因向量
  file  = "cellcycle_featureplot.pdf",
  nrow  = 3, ncol = 4,
  reduction = "tsne",
  cols = c("lightgrey", "#E41B1B")
)


geneids <- c(
  "VRK1", "MYCN", "SPIN4", "PAG1", "OLFML3", "FMOD", "BICC1",
  "MATN2", "SCD5", "DDIT4L", "SPOCK1", "FZD1", "METTL7A", "COL4A5"
)

featureplot_to_pdf(
  obj,
  genes = geneids,          # 你的基因向量
  file  = "geneids_featureplot.pdf",
  nrow  = 3, ncol = 4,
  reduction = "tsne",
  cols = c("lightgrey", "#1f77b4")
)


obj <- AddModuleScore(obj, features = list(geneids), name = "GeneSet")
# 结果存储在 obj@meta.data$GeneSet1

FeaturePlot(
  obj,
  features  = "GeneSet1",
  reduction = "tsne",
  cols      = c("lightgrey", "#2ca02c")
)



# 以簇为身份（通常默认就是 seurat_clusters）
Idents(obj) <- "seurat_clusters"

# 找正向marker（在本簇表达更高），常用稳妥参数
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,          # 只要上调基因
  test.use = "wilcox",      # 也可用 "MAST"/"DESeq2"（慢）
  min.pct = 0.25,           # 至少25%细胞表达
  logfc.threshold = 0.25    # 至少0.25对数倍数变化
)

# 多重检验后过滤显著基因
markers_sig <- subset(markers, p_val_adj < 0.05)

# 每簇取前10个（按 avg_log2FC 排序）
library(dplyr)
top10 <- markers_sig %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

# 可视化：热图 & 点图
DoHeatmap(obj, features = unique(top10$gene), group.by = "seurat_clusters")
DotPlot(obj, features = unique(top10$gene), group.by = "seurat_clusters") + RotatedAxis()

# 保存
write.csv(markers_sig, "markers_per_cluster.csv", row.names = FALSE)




# 3. 绘制 DotPlot 图
p_all_markers <- DotPlot(
  obj,
  features = feat,
  scale = TRUE,
  assay = 'RNA',
  group.by = "RNA_snn_res.0.5"  # 根据聚类进行分组
) +
  theme_bw() + 
  # scale_color_gradientn(colors = brewer.pal(9, "BuPu")) +  # 使用 BuPu 调色板
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 0, unit = 'cm'),
    axis.text.x = element_text(color = "black", size = 12, angle = 45, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 12),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black")
  )
# 显示图形
print(p_all_markers)


#############画图展示整理后的分簇

# 建立映射关系：0,1,2 -> 0; 3 -> 1; 4 -> 2; 5 -> 3
mapping <- c("0" = "0", "1" = "0", "2" = "0","3" = "0","4" = "0","5" = "0",
             "6" = "1", "7" = "2", "8" = "3")

# 应用映射到 RNA_snn_res.0.5
new_clusters <- mapping[ as.character(obj@meta.data[["RNA_snn_res.0.9"]]) ]

# 写入 seurat_clusters
obj@meta.data[["seurat_clusters"]] <- factor(new_clusters,
                                             levels = c("0","1","2","3"))

# 检查结果
table(obj@meta.data[["seurat_clusters"]])

Idents(obj) <- "seurat_clusters"

new.cluster.ids <- c("0"="Epi", 
                     "1"="T_NK", 
                     "2"="Myeloid", 
                     "3"="Prof")
obj <- RenameIdents(obj, new.cluster.ids)                        
obj$celltype <- obj@active.ident
DimPlot(obj, label = T,reduction = 'tsne',pt.size = 1,group.by = "celltype")

library("scplotter")

p1<-CellDimPlot(obj, group_by = "RNA_snn_res.0.9", reduction = "tsne", theme= "theme_blank")+
  scale_color_manual(values=my_cols)

p2<-CellDimPlot(obj, group_by = "celltype", reduction = "tsne", theme= "theme_blank")+
  scale_color_manual(values=my_cols)

# 不用 collect（容易触发 bug），先用 keep
p1 + p2 






cellcyclegene <- c(
  "AURKA","AURKB",
  "CDC6","CDC45","CDC25A","CDC25B","CDC25C","FZR1",
  "CDK1","CDK2","CDK3","CDK4","CDK5","CDK6","CDK7","CDK9",
  "CDT1","CHEK1","CHEK2","MYC",
  "CCNA2","CCNB1","CCNB2","CCND1","CCND2","CCND3","CCNE1","CCNE2","CCNF","CCNH",
  "E2F1","FOXM1","GMNN","MKI67",
  "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7",
  "PKMYT1",
  "ORC1","ORC2","ORC3","ORC4","ORC5","ORC6",
  "CDKN2B","CDKN2A","CDKN2C","CDKN2D","CDKN1A","CDKN1B","CDKN1C",
  "TP53","PCNA","RB1","RBL1","RBL2","SKP2","WEE1","PLK1"
)

FeaturePlot(obj,features = cellcyclegene,order = T,reduction = "tsne",cols = c("lightgrey", "#E41B1B"))
head(obj)
s.genes


# 定义基因集
genes <- cellcyclegene

# PDF 输出
pdf("cellcycle_featureplot.pdf", width = 14, height = 10)  # A4横版, 3x4 合适

# 分页循环，每页 12 个基因
for (i in seq(1, length(genes), by = 12)) {
  # 当前页的基因
  #i=1
  subset_genes <- genes[i:min(i+11, length(genes))]
  
  # 画图并组合成 3x4
  plots <- lapply(subset_genes, function(g) {
    FeaturePlot(obj, features = g, order = TRUE, reduction = "tsne",
                cols = c("lightgrey", "#E41B1B")) + 
      ggtitle(g) + 
      theme(plot.title = element_text(hjust = 0.5, size = 12))
  })
  
  # patchwork 拼图
  wrap_plots(plots, ncol = 4, nrow = 3)
}

dev.off()


# Seurat 内置的 cell cycle 基因集
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

# 打分 & 写入 obj@meta.data$cellcycle
obj <- CellCycleScoring(obj,
                        s.features = s.genes,
                        g2m.features = g2m.genes,
                        set.ident = FALSE)

# 查看新加的列
head(obj@meta.data[, c("S.Score", "G2M.Score", "Phase")])

FeaturePlot(obj, features = "S.Score", reduction = "tsne",
            cols = c("lightgrey", "#1f77b4"), order = TRUE)

FeaturePlot(obj, features = "G2M.Score", reduction = "tsne",
            cols = c("lightgrey", "#d62728"), order = TRUE)



FeaturePlot(obj,features = S.Score,order = T,reduction = "tsne",cols = c("lightgrey", "#E41B1B"))




#--------------------------------------#
#---画出一组基因在tsne图的表达---#
#--------------------------------------#

featureplot_to_pdf <- function(obj, genes, file = "cellcycle_featureplot.pdf",
                               nrow = 3, ncol = 4,
                               reduction = "tsne",
                               cols = c("lightgrey", "#E41B1B")) {
  stopifnot(requireNamespace("patchwork", quietly = TRUE))
  library(ggplot2)
  
  # 只保留对象中存在的基因，并提示缺失
  present <- intersect(genes, rownames(obj))
  missing <- setdiff(genes, rownames(obj))
  if (length(missing) > 0) {
    message("These genes are not in the object and will be skipped: ",
            paste(missing, collapse = ", "))
  }
  if (length(present) == 0) stop("No valid genes to plot.")
  
  # 打开 PDF 设备；确保出错也能关闭
  grDevices::pdf(file, width = 14, height = 10)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  # 分页，每页 nrow*ncol 个图
  page_size <- nrow * ncol
  idx <- seq(1, length(present), by = page_size)
  for (i in idx) {
    subset_genes <- present[i:min(i + page_size - 1, length(present))]
    
    plots <- lapply(subset_genes, function(g) {
      p <- FeaturePlot(obj, features = g, order = TRUE,
                       reduction = reduction, cols = cols) +
        ggtitle(g) +
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              legend.position = "none")
      # 避免某些版本返回列表
      if (inherits(p, "list")) p <- p[[1]]
      p
    })
    
    # 如果最后一页不足，填充空白占位，避免 patchwork 报错
    if (length(plots) < page_size) {
      n_blank <- page_size - length(plots)
      blanks <- replicate(n_blank, ggplot() + theme_void(), simplify = FALSE)
      plots <- c(plots, blanks)
    }
    
    pg <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow)
    print(pg)  # 关键：把拼好的对象“打印”到当前 pdf 设备
  }
  
  invisible(NULL)
}


cellcyclegene <- c(
  "AURKA","AURKB",
  "CDC6","CDC45","CDC25A","CDC25B","CDC25C","FZR1",
  "CDK1","CDK2","CDK3","CDK4","CDK5","CDK6","CDK7","CDK9",
  "CDT1","CHEK1","CHEK2","MYC",
  "CCNA2","CCNB1","CCNB2","CCND1","CCND2","CCND3","CCNE1","CCNE2","CCNF","CCNH",
  "E2F1","FOXM1","GMNN","MKI67",
  "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7",
  "PKMYT1",
  "ORC1","ORC2","ORC3","ORC4","ORC5","ORC6",
  "CDKN2B","CDKN2A","CDKN2C","CDKN2D","CDKN1A","CDKN1B","CDKN1C",
  "TP53","PCNA","RB1","RBL1","RBL2","SKP2","WEE1","PLK1"
)


featureplot_to_pdf(
  obj,
  genes = cellcyclegene,          # 你的基因向量
  file  = "cellcycle_featureplot.pdf",
  nrow  = 3, ncol = 4,
  reduction = "tsne",
  cols = c("lightgrey", "#E41B1B")
)



# 提取表达矩阵（通常是 RNA assay 的 count 矩阵）
expr <- GetAssayData(obj, assay = "RNA", slot = "counts")

# 提取元数据
meta <- obj@meta.data
matr.filter <- function(mat, min.cells = 10, min.genes = 10) {
  # 过滤掉表达少于 min.cells 个细胞的基因
  gene_filter <- Matrix::rowSums(mat > 0) >= min.cells
  # 过滤掉检测到的基因数少于 min.genes 的细胞
  cell_filter <- Matrix::colSums(mat > 0) >= min.genes
  mat[gene_filter, cell_filter]
}
library(ROGUE)
# 应用过滤函数
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)
#计算表达熵模型，这是后续的基础：
ent.res <- SE_fun(expr)#表达异质性评估函数。
SEplot(ent.res)

# 计算 ROGUE 值
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.res.sample <- rogue(
  expr = expr,
  labels = as.character(meta$RNA_snn_res.0.3),
  samples = as.character(meta$group),
  platform = "UMI",
  span = 0.6
)


write.csv(rogue.res.sample, file = "rogue.res.sample.csv")
rogue_sample <- read.csv("rogue.res.sample.csv", header = TRUE, row.names = 1)
head(rogue_sample)
myColor <- c(
  "#E41B1B", "#4376AC", "#48A75A", "#87638F", "#D87F32",
  "#737690", "#D690C6", "#B17A7D", "#847A74", "#4285BF",
  "#204B75", "#588257", "#B6DB7B", "#E3BC06", "#FA9B93",
  "#E9358B", "#A0094E", "#999999", "#6FCDDC", "#BD5E95"
)

# 整理数据用于绘图
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrastr)  # 如果你使用了 geom_jitter_rast

plot_rogue_sample <- rogue_sample %>%
  tidyr::gather(key = clusters, value = ROGUE) %>%
  dplyr::filter(!is.na(ROGUE))

# 绘制 boxplot + jitter 图
ggplot(data = plot_rogue_sample, aes(clusters, ROGUE, color = clusters)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter_rast(shape = 16, position = position_jitter(0.2)) +
  scale_color_manual(values = myColor_sorted) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, colour = "black",angle = 45),
    axis.title = element_text(size = 13, colour = "black")
  ) +
  labs(x = "", y = "ROGUE index") +
  ylim(0, 1)





