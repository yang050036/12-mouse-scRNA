library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)
library(clusterProfiler)
library(CellChat)
library(harmony)

###### standard analysis #####
out = "result"
if (!dir.exists(out)) dir.create(out)
mtx_path = "matrix/"
files = list.files(mtx_path)
data.list = lapply(files, function(x){
  obj <- CreateSeuratObject(Read10X(file.path(mtx_path, x)),
                            project = x, 
                            assay = "RNA", 
                            min.cells = 3, 
                            min.feaatures = 200)
  obj[["percent.mt"]] = PercentageFeatureSet(object = obj, pattern = "^mt")
  return(obj)
})
data = merge(data.list[[1]], data.list[-1])
qcout = file.path(out, "QC")
if (!dir.exists(qcout)) dir.create(qcout)
p = VlnPlot(data, features = c("nFeature_RNA"))
ggsave(file.path(qcout, "vlnplot.nGene.png"), p)
ggsave(file.path(qcout, "vlnplot.nGene.pdf"), p)
p = VlnPlot(data, features = c("nCount_RNA"))
ggsave(file.path(qcout, "vlnplot.nUMI.png"), p)
ggsave(file.path(qcout, "vlnplot.nUMI.pdf"), p)
p = VlnPlot(data, features = c("percent.mt"))
ggsave(file.path(qcout, "vlnplot.mt.png"), p)
ggsave(file.path(qcout, "vlnplot.mt.pdf"), p)
data = subset(data, nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
#去除双细胞
library(DoubletFinder)
data.list = SplitObject(data, split.by = "orig.ident")
doublet_num = 0
data.list = lapply(data.list, function(obj){
  obj = NormalizeData(obj)
  obj = FindVariableFeatures(obj, selection.method  = "vst", nfeatures = 2000)
  obj = ScaleData(obj)
  obj = RunPCA(obj)
  obj = RunUMAP(obj, dims = 1:20)
  # obj = FindNeighbors(obj, dims = 1:20)
  # obj = FindClusters(obj, resolution = 0.8)
  if (ncol(obj) > 7000 & ncol(obj) < 8700){
    rate = 0.031
  }else if(ncol(obj) > 8700 & ncol(obj) < 10500){
    rate = 0.046
  }else if(ncol(obj) > 10500 & ncol(obj) < 12200){
    rate = 0.054
  }else if(ncol(obj) > 12200 & ncol(obj) < 14000){
    rate = 0.061
  }else if(ncol(obj) > 14000 & ncol(obj) < 15700){
    rate = 0.069
  }else if(ncol(obj) > 15700 & ncol(obj) < 17400){
    rate = 0.076
  }
  nExp_poi <- round(rate*nrow(obj@meta.data)) 
  obj = doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  DF.name = colnames(obj@meta.data)[grepl("DF.classification", colnames(obj@meta.data))]
  obj = obj[, obj@meta.data[, DF.name] == "Singlet"]
  return(obj)
})

## 剔除BL03重新分析
data = subset(data, orig.ident != "BL03")
DefaultAssay(data) = "RNA"

##改用harmony去除批次效应
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
s.genes[which(s.genes == "Mlf1ip")] = "Cenpu"
g2m.genes[which(g2m.genes == "Fam64a")] = "Pimreg"
g2m.genes[which(g2m.genes == "Hn1")] = "Jpt1"
data = CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
data = FindVariableFeatures(data, nfeatures = 3000)
options(future.globals.maxSize = 3000 * 1024^2)
data = ScaleData(data, vars.to.regress = c("S.Score", "G2M.Score"))
data = RunPCA(data, npcs = 50)
png("tmp/convergence.png", width = 1500, height = 1500, res = 300)
data = RunHarmony(data, group.by.vars = "orig.ident", plot_convergence = TRUE, reduction = "pca")
dev.off()
#qc
Idents(data) = "orig.ident"
for (sample in names(data.list)){
  sub = data.list[[sample]]
  p = FeatureScatter(sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend() +
    labs("Correlative of nUMI and nGene", x = "nUMI perl cell", y = "nGene perl cell")
  ggsave(file.path(qcout, paste("correlative", sample, "png", sep = ".")), p, bg = "white")
  ggsave(file.path(qcout, paste("correlative", sample, "pdf", sep = ".")), p, bg = "white")
}
Idents(data) = "orig.ident"
p = VlnPlot(data, features = c("nFeature_RNA"), pt.size = 0) + labs("nGene")
ggsave(file.path(qcout, "vlnplot.nGene.png"), p, bg = "white")
ggsave(file.path(qcout, "vlnplot.nGene.pdf"), p, bg = "white")
p = VlnPlot(data, features = c("nCount_RNA"), pt.size = 0) + labs("nUMI")
ggsave(file.path(qcout, "vlnplot.nUMI.png"), p, bg = "white")
ggsave(file.path(qcout, "vlnplot.nUMI.pdf"), p, bg = "white")
p = VlnPlot(data, features = c("percent.mt"), pt.size = 0) + labs("percent.mito")
ggsave(file.path(qcout, "vlnplot.mt.png"), p, bg = "white")
ggsave(file.path(qcout, "vlnplot.mt.pdf"), p, bg = "white")

p = ElbowPlot(data, ndims = 20, reduction = "pca")
ggsave("tmp/elbowplot.png", p)
data = RunUMAP(data, reduction = "harmony", dims = 1:8, n.neighbors = 50L)
data = RunTSNE(data, reduction = "harmony", dims = 1:8)
data = FindNeighbors(data, reduction = "harmony", dims = 1:8, k.param = 110)
data = FindClusters(data, resolution = 1:15/10)
p = DimPlot(data, group.by = "orig.ident", reduction = "umap", raster=FALSE)
ggsave("tmp/tsne.sample.png", p)
p = clustree(data, prefix = "RNA_snn_res.")
ggsave("tmp/clustree.png", p, height = 10, width = 18)
Idents(data) = "RNA_snn_res.0.2"
data$cluster = Idents(data)
p = DimPlot(data, group.by = "cluster", reduction = "umap", label = TRUE)
ggsave("tmp/tsne.cluster.png", p)

###### group infomation #####
group = read.table("tmp/group.txt", sep = "\t", header = TRUE)
data$group = group$Group[match(data$orig.ident, group$Sample)]
###### DE of cluster #####
DefaultAssay(data) = "RNA"
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 1500 * 1024^2)
Idents(data) = "cluster"
deg = FindAllMarkers(data)
deout = file.path(out, "DEG")
if (!dir.exists(deout)) dir.create(deout)
write.csv(deg, file.path(deout, "DEG.cluster.csv"), row.names = F, quote = F)
#top 20 vlnplot
vout = file.path(deout, "DEG_cluster_top20")
if (!dir.exists(vout)) dir.create(vout)
subdeg = deg %>% group_by(cluster) %>% top_n(20, wt=avg_log2FC)
for (name in unique(subdeg$cluster)){
  subout = file.path(vout, name)
  if (!dir.exists(subout)) dir.create(subout)
  genes = subdeg %>% filter(cluster == name) %>% .$gene
  for (gene in genes){
    p = VlnPlot(data, features = gene, group.by = "cluster", pt.size = 0) + NoLegend()
    ggsave(file.path(subout, paste("vln", gene, "png", sep = ".")), p, height = 4, bg = "white")
    ggsave(file.path(subout, paste("vln", gene, "pdf", sep = ".")), p, height = 4, bg = "white")
  }
}

###### define cluster #####
library(SingleR)
library(BiocParallel)
ref = readRDS("/mnt/beegfs/Research/Database/Seurat/Mouse_RNA-seq.rds")
DefaultAssay(data) = "RNA"
sce = as.SingleCellExperiment(data)
# pred <- SingleR(sce, ref=ref, assay.type.test=1, labels=ref$label.fine,
#                   BPPARAM=MulticoreParam(8))
pred <- SingleR(sce, ref=ref, assay.type.test=1, labels=ref$label.fine, method = "cluster",
                BPPARAM=MulticoreParam(8))
data$celltype_singler=pred$labels
p = DimPlot(data, group.by = "celltype_singler", reduction = "umap", label = TRUE)
ggsave("tmp/tsne.celltype.png", p)

Idents(data) = data$cluster
data = RenameIdents(data,
                    "0" = "Stage II neutrophil",
                    "1" = "Stage I neutrophil",
                    "2" = "T Cell",
                    "3" = "Monocyte",
                    "4" = "B Cell",
                    "5" = "NKT",
                    "6" = "NK",
                    "7" = "Monocyte",
                    "8" = "Stage I neutrophil",
                    "9" = "B Cell")
data$celltype = Idents(data)
p = DimPlot(data, group.by = "celltype", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap") + NoLegend()
ggsave("tmp/tsne.celltype.png", p)

p = DimPlot(data, group.by = "group", reduction = "umap")
ggsave("tmp/tsne.group.png", p)

###### boxplot #####
library(cowplot)
df <- data.frame(data$orig.ident,data$cluster,data$celltype,data$nCount_RNA)
colnames(df) <- c('times','cluster','celltype','nUMI')
p = ggplot(df,aes(x=celltype,y=nUMI)) + 
  geom_boxplot() + coord_flip() +
  xlab("Celltypes") + ylab("nUMI") + scale_fill_discrete(guide=FALSE) +
  #  scale_y_continuous(expand=c(0,0),limits=c(0,17000),breaks=c(0,5000,10000,15000),labels=c('0','5000','10000','>=15000')) + 
  theme_cowplot()
bout = file.path(out, "Stat_of_per_cluster")
if (!dir.exists(bout)) dir.create(bout)
ggsave(file.path(bout, "boxplot_of_per_celltype.png"),p,bg="white")
ggsave(file.path(bout, "boxplot_of_per_celltype.pdf"),p,bg="white")
p = ggplot(df,aes(x=cluster,y=nUMI)) + 
  geom_boxplot() + coord_flip() +
  xlab("Clusters") + ylab("nUMI") + scale_fill_discrete(guide=FALSE) +
  #  scale_y_continuous(expand=c(0,0),limits=c(0,17000),breaks=c(0,5000,10000,15000),labels=c('0','5000','10000','>=15000')) + 
  theme_cowplot()
ggsave(file.path(bout, "boxplot_of_per_cluster.png"),p,bg="white")
ggsave(file.path(bout, "boxplot_of_per_cluster.pdf"),p,bg="white")

###### umap #####
uout = file.path(out, "UMPA")
if (!dir.exists(uout)) dir.create(uout)
p = DimPlot(data, group.by = "celltype", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap", raster = FALSE)
ggsave(file.path(uout, "umap.celltype.png"), p, bg = "white", width = 9)
ggsave(file.path(uout, "umap.celltype.pdf"), p, bg = "white", width = 9)
p = DimPlot(data, group.by = "orig.ident", reduction = "umap", raster = FALSE)
ggsave(file.path(uout, "umap.sample.png"), p, bg = "white")
ggsave(file.path(uout, "umap.sample.pdf"), p, bg = "white")
p = DimPlot(data, group.by = "cluster", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap", raster = FALSE)
ggsave(file.path(uout, "umap.cluster.png"), p, bg = "white")
ggsave(file.path(uout, "umap.cluster.pdf"), p, bg = "white")

mdata = subset(data, celltype == "Monocyte")
p = DimPlot(mdata, group.by = "celltype", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap")
ggsave(file.path(uout, "Monocyte_umap.celltype.png"), p, bg = "white")
ggsave(file.path(uout, "Monocyte_umap.celltype.pdf"), p, bg = "white")
p = DimPlot(mdata, group.by = "orig.ident", reduction = "umap")
ggsave(file.path(uout, "Monocyte_umap.sample.png"), p, bg = "white")
ggsave(file.path(uout, "Monocyte_umap.sample.pdf"), p, bg = "white")
p = DimPlot(mdata, group.by = "cluster", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap")
ggsave(file.path(uout, "Monocyte_umap.cluster.png"), p, bg = "white")
ggsave(file.path(uout, "Monocyte_umap.cluster.pdf"), p, bg = "white")


###### DE of celltype #####
Idents(data) = "celltype"
deg = FindAllMarkers(data)
if (!dir.exists(deout)) dir.create(deout)
write.csv(deg, file.path(deout, "DEG.celltype.csv"), row.names = F, quote = F)
#top 10 vlnplot
vout = file.path(deout, "DEG_celltype_top10")
if (!dir.exists(vout)) dir.create(vout)
subdeg = deg %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC)
for (name in unique(subdeg$cluster)){
  subout = file.path(vout, name)
  if (!dir.exists(subout)) dir.create(subout)
  genes = subdeg %>% filter(cluster == name) %>% .$gene
  for (gene in genes){
    p = VlnPlot(data, features = gene, group.by = "celltype", pt.size = 0) + NoLegend()
    ggsave(file.path(subout, paste("vln", gene, "png", sep = ".")), p, height = 4, bg = "white")
    ggsave(file.path(subout, paste("vln", gene, "pdf", sep = ".")), p, height = 4, bg = "white")
  }
}

###### Monocyte stat #####
msout = file.path(out, "Monocyte_stat")
if (!dir.exists(msout)) dir.create(msout)
tb = table(mdata$celltype, mdata$orig.ident)
tb = cbind(celltype = rownames(tb), tb)
tb = as.data.frame(tb)
write.csv(tb, file.path(msout, "Monocyte_stat.csv"), row.names = F, quote = F)

exp = GetAssayData(mdata, assay = "RNA", slot = "data")
write.csv(exp, file.path(msout, "Monocyte_Expression.csv"), row.names = T, quote = F)

###### celltype stat #####
csout = file.path(out, "Celltype_stat")
if (!dir.exists(csout)) dir.create(csout)
tb = table(data$celltype)
tb = cbind(celltype = rownames(tb), cell_num = tb)
tb = as.data.frame(tb)
write.csv(tb, file.path(csout, "celltype_stat.csv"), row.names = F, quote = F)

###### Monocyte repair gene #####
library(ggpubr)
library(stringr)
mrout = file.path(out, "Monocyte_repair_gene")
if (!dir.exists(mrout)) dir.create(mrout)
##Tgfb改成Tgfb1
genes = "Arg1/Vegfa/Il10/Tgfb1/Cd68/Cd163/Azin1/ORM2/ACTG1/NAGLU/Serpinc1/Suco/F5/Olfml2b/Dusp10/Lamb3/Creb3l1/Tgm2/Ppp1r16b/Mmp9/Sytl4/Il12a/Mef2d/S100a1/Adamtsl4/Zc3h12a/Ctnnbip1/Mmp23/Emilin1/Arhgap24/Hspb1/Flt1/Prss2/Gp9/P3h3/Fgfr1op2/Gp6/Pepd/Hps5/Smpd1/Rnh1/Cd151/Vsir/Adamts14/Egr2/Tbxa2r/Elk3/Lemd3/Mmp19/Col4a2/F10/Vegfc/Rab3a/Colgalt1/Hmox1/Exoc8/Mmp27/Angptl6/Ets1/Mcam/Mpzl3/Jaml/Il10ra/Tgfbr2/Anxa6/Clec10a/Myo1c/Plxdc1/Timp2/Cbx8/E2f3/Nqo2/F13a1/Rreb1/Tgfbi/Rhoj/Ahnak2/Ptp4a3/Mkl1/Irak4/Ano6/Gp1bb/Ccdc80/Il10rb/Mmp25/Pdpk1/Jmjd8/Adamts10/Ager/Vegfa/Treml1/Lrg1/Egr1/Sema6a/Pdgfrb/Vegfb/Ric1/Myof/Angptl2/Col27a1/Mmrn1/Gfod2/Efemp2/Timp1/Dcbld2/Ptpn14/Procr/F9/Adamtsl1/Plpp3/Epha2/Mmp17/Tfpi2/Anpep/Adamtsl3/Serpinh1/Mmp11/Adamtsl5/Gas6/Cysltr2/Col5a3/Col23a1/Sparc/Kazald1/Tns2/Col17a1/F8/Creb1/Dstyk/Qsox1/Itga8/Atf2/Sdc4/Atp6ap2/Sdcbp/Prkd2/Furin/Lrrc32/Wnt11/Tgfb1i1/Shcbp1/Ier2/Carmil2/Nrp1/Cep57/Rnf111/Kif9/Tvp23b/F12/Sdc1/Nrros/Noxo1/Ddr1/Ltbp1/Tcf4/Fgfbp3/Adamts3/Apold1/Col14a1/Osm/Anxa1"
genes = unlist(strsplit(genes, "/"))
DefaultAssay(data) = "RNA"
# mdata = subset(data, celltype == "Monocyte")
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
  }
  sub = subset(mdata, orig.ident %in% samples)
  compare_list = combn(as.character(unique(sub$group)), 2, simplify = FALSE)
  subout = file.path(mrout, source)
  if (!dir.exists(subout)) dir.create(subout)
  lapply(genes, function(gene){
    gene = str_to_title(gene)
    pt.size = AutoPointSize(sub)
    vln = VlnPlot(sub, features = gene, group.by = "group", slot = "count")
    df = vln$data
    vln = ggplot(df, aes_string("ident", gene, fill = "ident")) + geom_violin(trim = TRUE, scale = "width")
    vln = vln + geom_jitter(height = 0, size = pt.size, show.legend = FALSE)
    vln = vln + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
    vln = vln + xlab("Identity") + ylab("Expression Level") + labs(title = gene) + NoLegend()
    vln = vln + stat_compare_means(
      comparisons = compare_list, method = "wilcox.test", 
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                         symbols = c("****", "***", "**", "*", "ns")))
    ggsave(file.path(subout, paste0("vln.", gene, ".png")), vln, bg = "white")
    ggsave(file.path(subout, paste0("vln.", gene, ".pdf")), vln, bg = "white")
  })
}

###### group deg GO and GSVA #####
library(GSVA)
library(limma)
library(msigdbr)
library(clusterProfiler)
library(KEGGREST)
library(org.Mm.eg.db)
library(stringr)
library(GSEABase)
library(ggthemes)
library(ggprism)
## function ##
gsva_de = function(gsva_result, design, compare, source, gsvaout){
  fit = lmFit(gsva_result, design)
  fit2 = contrasts.fit(fit, compare)
  fit3 = eBayes(fit2)
  Diff = topTable(fit3, coef = 2, number = Inf, adjust.method='BH')
  ##barplot
  Diff$threshold = ifelse(Diff$t > -2, ifelse(Diff$t >=2, "Up", "NoSig"), "Down")
  Diff = Diff %>% arrange(t)
  Diff$id = rownames(Diff)
  Diff <- t(apply(Diff,1,function(x){
    if(nchar(x[ncol(Diff)]) > 45){
      x[ncol(Diff)] <- substr(x[ncol(Diff)],1,45)
      x[ncol(Diff)] <- paste0(x[ncol(Diff)],'...')
    }
    return(x)
  }))
  Diff = as.data.frame(Diff)
  Diff$t = as.numeric(Diff$t)
  Diff$id <- factor(Diff$id,levels = Diff$id)
  lim = ceiling(max(abs(Diff$t)))
  p <- ggplot(data = Diff,aes(x = id,y = t,fill = threshold)) +
    geom_col()+
    coord_flip() +
    ylim(-lim, lim) +
    scale_fill_manual(values = c('Up'= '#36638a','NoSig'='#cccccc','Down'='#7bcd7b')) +
    geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
    xlab('') + 
    ylab('t value of GSVA score') + #注意坐标轴旋转了
    guides(fill=F)+ # 不显示图例
    theme_prism(border = T) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  # 小于-2的数量
  low1 <- Diff %>% filter(t < -2) %>% nrow()
  # 小于0总数量
  low0 <- Diff %>% filter( t < 0) %>% nrow()
  # 小于2总数量
  high0 <- Diff %>% filter(t < 2) %>% nrow()
  # 总的柱子数量
  high1 <- nrow(Diff)
  p <- p + geom_text(data = Diff[1:low1,],aes(x = id,y = 0.1,label = id),
                     hjust = 0,color = 'black') + # 小于-1的为黑色标签
    geom_text(data = Diff[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
              hjust = 0,color = 'grey') + # 灰色标签
    geom_text(data = Diff[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
              hjust = 1,color = 'grey') + # 灰色标签
    geom_text(data = Diff[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
              hjust = 1,color = 'black') # 大于1的为黑色标签
  ggsave(file.path(gsvaout, paste("gsva", source, "png", sep = ".")), p, bg = "white", width = 8)
  ggsave(file.path(gsvaout, paste("gsva", source, "pdf", sep = ".")), p, bg = "white", width = 8)
  write.csv(Diff, file.path(gsvaout, paste("gsva", source, "csv", sep = ".")), row.names = T, quote = F)
}
## main ##
eout = file.path(out, "Enrichment")
if (!dir.exists(eout)) dir.create(eout)
Idents(mdata) = "group"
geneSet = getGmt("/mnt/beegfs/Research/Database/MSigDB/c5.go.v7.5.1.symbols.gmt")
get_GO_data = getFromNamespace("get_GO_data", "clusterProfiler")
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    ident.1 = "A-2"
    ident.2 = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    ident.1 = "B-2"
    ident.2 = "B-1"
  }
  subout = file.path(eout, source)
  goout = file.path(subout, "GO")
  gsvaout = file.path(subout, "GSAV")
  if (!dir.exists(subout)) dir.create(subout)
  if (!dir.exists(goout)) dir.create(goout)
  if (!dir.exists(gsvaout)) dir.create(gsvaout)
  sub = subset(mdata, orig.ident %in% samples)
  compare_list = combn(as.character(unique(sub$group)), 2, simplify = FALSE)
  deg = FindMarkers(sub, ident.1 = ident.1, ident.2 = ident.2, test.use = "DESeq2")
  #deg = deg %>% filter(avg_log2FC > 0.5, p_val_adj < 0.01)
  write.csv(deg, file.path(subout, paste("DEG", source, "csv", sep = ".")), row.names = T, quote = F)
  genes = rownames(deg)
  eg <- bitr(genes, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb='org.Mm.eg.db')
  BP_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="BP",pvalueCutoff = 0.05,readable = TRUE)
  MF_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="MF",pvalueCutoff = 0.05,readable = TRUE)
  CC_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="CC",pvalueCutoff = 0.05,readable = TRUE)
  bp <- as.data.frame(BP_go)
  cc <- as.data.frame(CC_go)
  mf <- as.data.frame(MF_go)
  write.csv(bp, file.path(goout, paste("GO", "BP", "csv", sep = ".")), row.names = F, quote = F)
  write.csv(cc, file.path(goout, paste("GO", "CC", "csv", sep = ".")), row.names = F, quote = F)
  write.csv(mf, file.path(goout, paste("GO", "MF", "csv", sep = ".")), row.names = F, quote = F)
  
  ## go_list
  bp = subset(bp, Count > 2)
  cc = subset(cc, Count > 2)
  mf = subset(mf, Count > 2)
  bpterms <- bp$ID
  ccterms <- cc$ID
  mfterms <- mf$ID
  go_list = c()
  if (length(bpterms) != 0){
    BP_GO <- get_GO_data("org.Mm.eg.db", "BP", "SYMBOL")
    bp_list <- BP_GO$PATHID2EXTID[bpterms]
    names(bp_list) <- paste(names(BP_GO$PATHID2NAME[bpterms]), "BP",BP_GO$PATHID2NAME[bpterms],sep=':')
    go_list = c(go_list, bp_list)
  }
  if (length(ccterms) != 0){
    CC_GO <- get_GO_data("org.Mm.eg.db", "CC", "SYMBOL") 
    cc_list <- CC_GO$PATHID2EXTID[ccterms]
    names(cc_list) <-  paste(names(CC_GO$PATHID2NAME[ccterms]), "CC",CC_GO$PATHID2NAME[ccterms],sep=":")
    go_list = c(go_list, bp_list)
  }
  if (length(mfterms) != 0){
    MF_GO <- get_GO_data("org.Mm.eg.db", "MF", "SYMBOL")
    mf_list <- MF_GO$PATHID2EXTID[mfterms]
    names(mf_list) <-  paste(names(MF_GO$PATHID2NAME[mfterms]), "MF",MF_GO$PATHID2NAME[mfterms],sep=":")
    go_list = c(go_list, bp_list)
  }
  
  ## GSVA
  exp = GetAssayData(sub, assay = "RNA", slot = "data")
  gsva_result = gsva(as.matrix(exp), method = "gsva", kcdf = "Gaussian", mx.diff=T, gset.idx.list=go_list, parallel.sz=10)
  if (source == "BL"){
    group <- factor(sub$group,levels = c("A-2", "A-1"))
    design <- model.matrix(~0+group)
    colnames(design) = c("A2", "A1")
    compare = makeContrasts("A2", "A1", levels = design)
    colnames(compare) = c("A-2", "A-1")
    rownames(compare) = c("A-2", "A-1")
    colnames(design) = c("A-2", "A-1")
  }else{
    group <- factor(sub$group,levels = c("B-2", "B-1"))
    design <- model.matrix(~0+group)
    colnames(design) = c("B2", "B1")
    compare = makeContrasts("B2", "B1", levels = design)
    colnames(compare) = c("B-2", "B-1")
    rownames(compare) = c("B-2", "B-1")
    colnames(design) = c("B-2", "B-1")
  }
  gsva_de(gsva_result, design, compare, source, gsvaout)
}

###### Monocyte Ly6C2#####
lout = file.path(out,"Ly6c2_anlysis")
if (!dir.exists(lout)) dir.create(lout)
DefaultAssay(mdata) = "RNA"
exp = FetchData(mdata,vars = "Ly6c2", slot = "count")
mdata$Ly6c2 = ifelse(exp$Ly6c2 > 0, "Ly6c2+", "Ly6c2-")
mdata$Ly6c2 = factor(mdata$Ly6c2, levels = c("Ly6c2+", "Ly6c2-"))
p = DimPlot(mdata, group.by = "Ly6c2", split.by = "Ly6c2", label = TRUE, label.box = TRUE)
ggsave(file.path(lout, "umap.Ly6c2.png"), p, width = 15)
ggsave(file.path(lout, "umap.Ly6c2.pdf"), p, width = 15)
p = DimPlot(mdata, group.by = "orig.ident", split.by = "Ly6c2")
ggsave(file.path(lout, "umap.sample.png"), p, width = 15)
ggsave(file.path(lout, "umap.sample.pdf"), p, width = 15)
p = DimPlot(mdata, group.by = "cluster", split.by = "Ly6c2", label = TRUE, label.box = TRUE)
ggsave(file.path(lout, "umap.cluster.png"), p, width = 15)
ggsave(file.path(lout, "umap.cluster.pdf"), p, width = 15)
p = DimPlot(mdata, group.by = "celltype", split.by = "Ly6c2", label = TRUE, label.box = TRUE)
ggsave(file.path(lout, "umap.celltype.png"), p, width = 15)
ggsave(file.path(lout, "umap.celltype.pdf"), p, width = 15)

## Ly6c2在单核细胞中的细胞数量
tb = table(mdata$group, mdata$Ly6c2)
tb = cbind(group = rownames(tb), tb)
tb = as.data.frame(tb)
write.csv(tb, file.path(lout, "stat.Ly6c2.csv"), row.names = F, quote = F)

## Ly6c2差异分析
deout = file.path(lout, "DEG")
if (!dir.exists(deout)) dir.create(deout)
for (i in c("Ly6c2+", "Ly6c2-")){
  subdata = subset(mdata, Ly6c2 == i)
  for(source in c("BL", "BM")){
    if (source == "BL"){
      samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
      ident.1 = "A-2"
      ident.2 = "A-1"
    }else{
      samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
      ident.1 = "B-2"
      ident.2 = "B-1"
    }
    subout = file.path(deout, paste0(i, "_", source))
    goout = file.path(subout, "GO")
    keggout = file.path(subout, "KEGG")
    if (!dir.exists(subout)) dir.create(subout)
    if (!dir.exists(goout)) dir.create(goout)
    if (!dir.exists(keggout)) dir.create(keggout)
    sub = subset(subdata, orig.ident %in% samples)
    Idents(sub) = "group"
    deg = FindMarkers(sub, ident.1 = ident.1, ident.2 = ident.2)
    # deg = deg %>% filter(abs(avg_log2FC) > 0.5, p_val_adj < 0.01)
    write.csv(deg, file.path(subout, paste("DEG", i, source, "csv", sep = ".")), row.names = T, quote = F)
    genes = rownames(deg)
    eg <- bitr(genes, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb='org.Mm.eg.db')
    BP_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="BP",pvalueCutoff = 0.05,readable = TRUE)
    MF_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="MF",pvalueCutoff = 0.05,readable = TRUE)
    CC_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="CC",pvalueCutoff = 0.05,readable = TRUE)
    KEGG <- enrichKEGG(gene = unique(eg$ENTREZID),keyType = "ncbi-geneid",organism = 'mmu',use_internal_data = TRUE,pvalueCutoff = 0.05)
    kegg <- as.data.frame(KEGG)
    kegg$genename <- unlist(lapply(kegg$geneID,  FUN =function(y){
      paste( unlist(lapply(unlist(strsplit(y,"/")),FUN=function(y){eg$SYMBOL[which(eg$ENTREZID ==y)]})) , collapse = "/")} ))
    bp <- as.data.frame(BP_go)
    cc <- as.data.frame(CC_go)
    mf <- as.data.frame(MF_go)
    write.csv(bp, file.path(goout, paste("GO", "BP", "csv", sep = ".")), row.names = F, quote = F)
    write.csv(cc, file.path(goout, paste("GO", "CC", "csv", sep = ".")), row.names = F, quote = F)
    write.csv(mf, file.path(goout, paste("GO", "MF", "csv", sep = ".")), row.names = F, quote = F)
    write.csv(kegg, file.path(keggout, paste("KEGG", "csv", sep = ".")), row.names = F, quote = F)
    
    ## deg expression
    exp = FetchData(sub, vars = genes, slot = "data")
    write.csv(exp, file.path(subout, paste("DEG_expression", i, source, "csv", sep = ".")), row.names = F, quote = F)
  }
}

## gene of repair
vout = file.path(lout,"Repair_gene_vln")
if (!dir.exists(vout)) dir.create(vout)
genes = "Arg1/Vegfa/Il10/Tgfb1/Cd68/Cd163/Azin1/ORM2/ACTG1/NAGLU/Serpinc1/Suco/F5/Olfml2b/Dusp10/Lamb3/Creb3l1/Tgm2/Ppp1r16b/Mmp9/Sytl4/Il12a/Mef2d/S100a1/Adamtsl4/Zc3h12a/Ctnnbip1/Mmp23/Emilin1/Arhgap24/Hspb1/Flt1/Prss2/Gp9/P3h3/Fgfr1op2/Gp6/Pepd/Hps5/Smpd1/Rnh1/Cd151/Vsir/Adamts14/Egr2/Tbxa2r/Elk3/Lemd3/Mmp19/Col4a2/F10/Vegfc/Rab3a/Colgalt1/Hmox1/Exoc8/Mmp27/Angptl6/Ets1/Mcam/Mpzl3/Jaml/Il10ra/Tgfbr2/Anxa6/Clec10a/Myo1c/Plxdc1/Timp2/Cbx8/E2f3/Nqo2/F13a1/Rreb1/Tgfbi/Rhoj/Ahnak2/Ptp4a3/Mkl1/Irak4/Ano6/Gp1bb/Ccdc80/Il10rb/Mmp25/Pdpk1/Jmjd8/Adamts10/Ager/Vegfa/Treml1/Lrg1/Egr1/Sema6a/Pdgfrb/Vegfb/Ric1/Myof/Angptl2/Col27a1/Mmrn1/Gfod2/Efemp2/Timp1/Dcbld2/Ptpn14/Procr/F9/Adamtsl1/Plpp3/Epha2/Mmp17/Tfpi2/Anpep/Adamtsl3/Serpinh1/Mmp11/Adamtsl5/Gas6/Cysltr2/Col5a3/Col23a1/Sparc/Kazald1/Tns2/Col17a1/F8/Creb1/Dstyk/Qsox1/Itga8/Atf2/Sdc4/Atp6ap2/Sdcbp/Prkd2/Furin/Lrrc32/Wnt11/Tgfb1i1/Shcbp1/Ier2/Carmil2/Nrp1/Cep57/Rnf111/Kif9/Tvp23b/F12/Sdc1/Nrros/Noxo1/Ddr1/Ltbp1/Tcf4/Fgfbp3/Adamts3/Apold1/Col14a1/Osm/Anxa1"
genes = unlist(strsplit(genes, "/"))
DefaultAssay(data) = "RNA"
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
  }
  sub = subset(mdata, orig.ident %in% samples)
  compare_list = combn(as.character(unique(sub$group)), 2, simplify = FALSE)
  subout = file.path(vout, source)
  if (!dir.exists(subout)) dir.create(subout)
  lapply(genes, function(gene){
    gene = str_to_title(gene)
    pt.size = AutoPointSize(sub)
    vln = VlnPlot(sub, features = gene, group.by = "group", split.by = "Ly6c2", slot = "count")
    df = vln$data
    vln = ggplot(df, aes_string("ident", gene, fill = "ident")) + geom_violin(trim = TRUE, scale = "width")
    vln = vln + geom_jitter(height = 0, size = pt.size, show.legend = FALSE)
    vln = vln + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~split)
    vln = vln + xlab("Identity") + ylab("Expression Level") + labs(title = gene) + NoLegend()
    vln = vln + stat_compare_means(
      comparisons = compare_list, method = "wilcox.test", 
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                         symbols = c("****", "***", "**", "*", "ns")))
    ggsave(file.path(subout, paste0("vln.", gene, ".png")), vln, bg = "white")
    ggsave(file.path(subout, paste0("vln.", gene, ".pdf")), vln, bg = "white")
  })
}


vout = file.path(out, "Gene_vln")
if (!dir.exists(vout)) dir.create(vout)
egenes = c()
for (gene in genes){
  if (gene %in% rownames(mdata)){
    p = VlnPlot(mdata, features = gene, pt.size = 0, group.by = "orig.ident")
    ggsave(file.path(vout, paste("vln", gene, "png", sep = ".")), p)
    ggsave(file.path(vout, paste("vln", gene, "pdf", sep = ".")), p)
  }else{
    egenes = c(egenes, gene)
  }
}

##### 其他 ####
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    case = "A-2"
    control = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    case = "B-2"
    control = "B-1"
  }
  sub = subset(mdata, orig.ident %in% samples)
  Idents(sub) = "orig.ident"
  p = StackedVlnPlot(sub, features = gene)
  ggsave(file.path("result-2022-07-12", paste("stackedvlnplot", source, "png", sep = ".")), p, height = 10, width = 5)
  ggsave(file.path("result-2022-07-12", paste("stackedvlnplot", source, "pdf", sep = ".")), p, height = 7, width = 2)
  
  subout = file.path("result-2022-07-12",  "Repair_gene_vln")
  if (!dir.exists(subout)) dir.create(subout)
  lapply(gene, function(g){
    g = str_to_title(g)
    pt.size = AutoPointSize(sub)
    vln = VlnPlot(sub, features = g, group.by = "group", slot = "data")
    df = vln$data
    df$ident = ifelse(df$ident == case, "MI", "Sham")
    vln = ggplot(df, aes_string("ident", g, fill = "ident")) + geom_violin(trim = TRUE, scale = "width")
    vln = vln + geom_jitter(height = 0, size = pt.size, show.legend = FALSE)
    vln = vln + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
    vln = vln + xlab("Group") + ylab("Expression Level") + labs(title = g) + NoLegend()
    vln = vln + stat_compare_means(
      comparisons = list(c("MI", "Sham")), method = "wilcox.test", 
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                         symbols = c("****", "***", "**", "*", "ns")))
    ggsave(file.path(subout, paste0(source, g, ".png")), vln, bg = "white", width = 2, height = 4)
    ggsave(file.path(subout, paste0(source, g, ".pdf")), vln, bg = "white", width = 2, height = 4)
  })
}

hout= file.path(out, "heatmap")
if (!dir.exists(hout)) dir.create(hout)
library(pheatmap)
for (source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL04","BL05","BL06")
    case = "A-2"
    control = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    case = "B-2"
    control = "B-1"
  }
  sub = subset(mdata, orig.ident %in% samples)
  DefaultAssay(sub) = "RNA"
  # df = AverageExpression(sub, vars = gene, slot = "data", group.by = "orig.ident")
  Idents(sub) = "group"
  df = lapply(samples, function(x){
    sample_data = subset(sub, orig.ident == x)
    Idents(sample_data) = "group"
    subdeg = FindMarkers(sample_data, ident.1 = case, ident.2 = control, logfc.threshold = 0)
    subdeg = subdeg[rownames(subdeg) %in% gene,]
    subdeg = subdeg %>% select(avg_log2FC)
    colnames(subdeg) = x
  })
  df = do.call(cbind, df)
  p = pheatmap(df, cellwidth = 20, cellheight = 20, cluster_row = FALSE, 
               display_numbers = matrix(ifelse(df > 0.1, ifelse(df > 0.25, "**", "*"), ""), nrow(df)))
  ggsave(file.path(hout, paste("heatmap", source, "png", sep = ".")), p)
}


##### 2022-07-14 ######
out = "/mnt/beegfs/Research/Users/yangyupeng/scRNA/KY202205250064/2022-07-04/03.seurat/KY202205250064/result-2022-07-14"
if (!dir.exists(out)) dir.create(out)
uout = file.path(out, "UMAP")
if (!dir.exists(uout)) dir.create(uout)
p = DimPlot(data, group.by = "celltype", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap", raster=FALSE)
ggsave(file.path(uout, "umap.celltype.png"), p, width = 9)
ggsave(file.path(uout, "umap.celltype.pdf"), p, width = 9)
p = DimPlot(data, group.by = "cluster", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap", raster=FALSE)
ggsave(file.path(uout, "umap.cluster.png"), p, width = 9)
ggsave(file.path(uout, "umap.cluster.pdf"), p, width = 9)
p = DimPlot(data, group.by = "orig.ident", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap", raster=FALSE)
ggsave(file.path(uout, "umap.sample.png"), p, width = 9)
ggsave(file.path(uout, "umap.sample.pdf"), p, width = 9)

p = DimPlot(mdata, group.by = "cluster", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap", raster=FALSE)
ggsave(file.path(uout, "umap.Monocyte_cluster.png"), p, width = 9)
ggsave(file.path(uout, "umap.Monocyte_cluster.pdf"), p, width = 9)
p = DimPlot(mdata, group.by = "celltype", label = TRUE, repel = TRUE,label.box = TRUE, reduction = "umap", raster=FALSE)
ggsave(file.path(uout, "umap.Monocyte_celltype.png"), p, width = 9)
ggsave(file.path(uout, "umap.Monocyte_celltype.pdf"), p, width = 9)
p = DimPlot(mdata, group.by = "orig.ident", reduction = "umap", raster=FALSE)
ggsave(file.path(uout, "umap.Monocyte_sample.png"), p, width = 9)
ggsave(file.path(uout, "umap.Monocyte_sample.pdf"), p, width = 9)

mout = file.path(out, "Marker_Plot")
if (!dir.exists(mout)) dir.create(mout)
levels = "4, 9, 3, 7, 6, 5, 2, 1, 8, 0"
levels = unlist(strsplit(levels, ", "))
data$cluster = factor(data$cluster, levels = levels)
markers = "Cd74, Cd79a, Cd79b, Igkc, Csf1r, Ctsc, Ctss, S100a4, Ms4a6c, Gzma, Klra8, Ncr1, Klrd1, Ccl5, Nkg7, Cd3d, Cd3g, Cd3e, Trbc2, Ltf, Ngp, Camp, Lcn2, Ifitm6, Cd177, Anxa1, Adpgk, Ly6g, Syne1, Dstn, Stfa2l1, Il1b, Ifitm1"
markers = unlist(strsplit(markers, ", "))
p = DotPlot(data, assay = "RNA", group.by = "cluster", features = rev(markers)) + coord_flip()
ggsave(file.path(mout, "dotplot.png"), p, bg = "white", height = 10, width = 10)
ggsave(file.path(mout, "dotplot.pdf"), p, bg = "white", height = 10, width = 10)
p = StackedVlnPlot(data, assay = "RNA", group.by = "cluster", features = markers, 
                   color.use = c(ggsci::pal_simpsons()(16), ggsci::pal_igv()(50)))
ggsave(file.path(mout, "stackedvlnplot.png"), p, bg = "white", height = 10)
ggsave(file.path(mout, "stackedvlnplot.pdf"), p, bg = "white", height = 10)

DefaultAssay(data) = "RNA"
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 1500 * 1024^2)
Idents(data) = "cluster"
deg = FindAllMarkers(data, )
deout = file.path(out, "DEG")
if (!dir.exists(deout)) dir.create(deout)
write.csv(deg, file.path(deout, "DEG.cluster.csv"), row.names = F, quote = F)
#top 20 vlnplot
vout = file.path(deout, "DEG_cluster_top20")
if (!dir.exists(vout)) dir.create(vout)
subdeg = deg %>% group_by(cluster) %>% top_n(20, wt=avg_log2FC)
for (name in unique(subdeg$cluster)){
  subout = file.path(vout, name)
  if (!dir.exists(subout)) dir.create(subout)
  genes = subdeg %>% filter(cluster == name) %>% .$gene
  for (gene in genes){
    p = VlnPlot(data, features = gene, group.by = "cluster", pt.size = 0) + NoLegend()
    ggsave(file.path(subout, paste("vln", gene, "png", sep = ".")), p, height = 4, bg = "white")
    ggsave(file.path(subout, paste("vln", gene, "pdf", sep = ".")), p, height = 4, bg = "white")
  }
}

Idents(data) = "celltype"
deg = FindAllMarkers(data)
if (!dir.exists(deout)) dir.create(deout)
write.csv(deg, file.path(deout, "DEG.celltype.csv"), row.names = F, quote = F)
#top 10 vlnplot
vout = file.path(deout, "DEG_celltype_top10")
if (!dir.exists(vout)) dir.create(vout)
subdeg = deg %>% group_by(cluster) %>% top_n(10, wt=avg_log2FC)
for (name in unique(subdeg$cluster)){
  subout = file.path(vout, name)
  if (!dir.exists(subout)) dir.create(subout)
  genes = subdeg %>% filter(cluster == name) %>% .$gene
  for (gene in genes){
    p = VlnPlot(data, features = gene, group.by = "celltype", pt.size = 0) + NoLegend()
    ggsave(file.path(subout, paste("vln", gene, "png", sep = ".")), p, height = 4, bg = "white")
    ggsave(file.path(subout, paste("vln", gene, "pdf", sep = ".")), p, height = 4, bg = "white")
  }
}

library(ggpubr)
library(stringr)
mrout = file.path(out, "Monocyte_repair_gene")
if (!dir.exists(mrout)) dir.create(mrout)
##Tgfb改成Tgfb1
genes = "Cdkn1a, Fosl2, Isg15, Kmt2c, Krit1, Myh9, Ncl, Nr4a1, Phip, Pim1, Prkx, Rnf213, Stat3, Syk, Tnfaip3, Wasf2, SOD1, YWHAZ"
genes = unlist(strsplit(genes, ", "))
DefaultAssay(data) = "RNA"
# mdata = subset(data, celltype == "Monocyte")
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    case = "A-2"
    control = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    case = "B-2"
    control = "B-1"
  }
  sub = subset(mdata, orig.ident %in% samples)
  subout = file.path(mrout, source)
  if (!dir.exists(subout)) dir.create(subout)
  lapply(genes, function(gene){
    gene = str_to_title(gene)
    pt.size = AutoPointSize(sub)
    vln = VlnPlot(sub, features = gene, group.by = "group", slot = "data")
    df = vln$data
    df$ident = ifelse(df$ident == case, "MI", "Sham")
    df$ident = factor(df$ident, levels = c("Sham", "MI"))
    compare_list = combn(as.character(unique(df$ident)), 2, simplify = FALSE)
    vln = ggplot(df, aes_string("ident", gene, fill = "ident")) + geom_violin(trim = TRUE, scale = "width")
    vln = vln + geom_jitter(height = 0, size = pt.size, show.legend = FALSE)
    vln = vln + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
    vln = vln + xlab("Identity") + ylab("Expression Level") + labs(title = gene) + NoLegend()
    vln = vln + stat_compare_means(
      comparisons = compare_list, method = "wilcox.test", 
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                         symbols = c("****", "***", "**", "*", "ns")))
    ggsave(file.path(subout, paste0("vln.", gene, ".png")), vln, bg = "white", width = 3)
    ggsave(file.path(subout, paste0("vln.", gene, ".pdf")), vln, bg = "white", width = 3)
  })
}

lout = file.path(out,"Ly6c2_anlysis")
if (!dir.exists(lout)) dir.create(lout)
DefaultAssay(mdata) = "RNA"
exp = mdata@assays$RNA@counts["Ly6c2",]
mdata$Ly6c2 = ifelse(exp > 0, "Ly6c2+", "Ly6c2-")
mdata$Ly6c2 = factor(mdata$Ly6c2, levels = c("Ly6c2+", "Ly6c2-"))
p = DimPlot(mdata, group.by = "Ly6c2", split.by = "Ly6c2", label = TRUE, label.box = TRUE)
ggsave(file.path(lout, "umap.Ly6c2.png"), p, width = 14)
ggsave(file.path(lout, "umap.Ly6c2.pdf"), p, width = 14)
p = DimPlot(mdata, group.by = "orig.ident", split.by = "Ly6c2")
ggsave(file.path(lout, "umap.sample.png"), p, width = 14)
ggsave(file.path(lout, "umap.sample.pdf"), p, width = 14)
p = DimPlot(mdata, group.by = "cluster", split.by = "Ly6c2", label = TRUE, label.box = TRUE)
ggsave(file.path(lout, "umap.cluster.png"), p, width = 14)
ggsave(file.path(lout, "umap.cluster.pdf"), p, width = 14)
p = DimPlot(mdata, group.by = "celltype", split.by = "Ly6c2", label = TRUE, label.box = TRUE)
ggsave(file.path(lout, "umap.celltype.png"), p, width = 14)
ggsave(file.path(lout, "umap.celltype.pdf"), p, width = 14)
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    case = "A-2"
    control = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    case = "B-2"
    control = "B-1"
  }
  sub = subset(mdata, orig.ident %in% samples)
  subout = file.path(lout, source)
  if (!dir.exists(subout)) dir.create(subout)
  lapply(genes, function(gene){
    gene = str_to_title(gene)
    pt.size = AutoPointSize(sub)
    vln = VlnPlot(sub, features = gene, group.by = "group", slot = "data", split.by = "Ly6c2")
    df = vln$data
    df$ident = ifelse(df$ident == case, "MI", "Sham")
    df$ident = factor(df$ident, levels = c("Sham", "MI"))
    compare_list = combn(as.character(unique(df$ident)), 2, simplify = FALSE)
    vln = ggplot(df, aes_string("ident", gene, fill = "ident")) + geom_violin(trim = TRUE, scale = "width")
    vln = vln + geom_jitter(height = 0, size = pt.size, show.legend = FALSE)
    vln = vln + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~split)
    vln = vln + xlab("Identity") + ylab("Expression Level") + labs(title = gene) + NoLegend()
    vln = vln + stat_compare_means(
      comparisons = compare_list, method = "wilcox.test", 
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                         symbols = c("****", "***", "**", "*", "ns")))
    ggsave(file.path(subout, paste0("vln.", gene, ".png")), vln, bg = "white", width = 7)
    ggsave(file.path(subout, paste0("vln.", gene, ".pdf")), vln, bg = "white", width = 7)
  })
}

eout = file.path(out, "Group_de")
if (!dir.exists(eout)) dir.create(eout)
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    ident.1 = "A-2"
    ident.2 = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    ident.1 = "B-2"
    ident.2 = "B-1"
  }
  subout = file.path(eout, source)
  if (!dir.exists(subout)) dir.create(subout)
  sub = subset(mdata, orig.ident %in% samples)
  compare_list = combn(as.character(unique(sub$group)), 2, simplify = FALSE)
  Idents(sub) = "group"
  deg = FindMarkers(sub, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0)
  deg = deg %>% filter(p_val_adj < 0.05)
  write.csv(deg, file.path(subout, paste("DEG", source, "csv", sep = ".")), row.names = T, quote = F)
}

eout = file.path(out, "Ly6c2_group_de")
if (!dir.exists(eout)) dir.create(eout)
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    ident.1 = "A-2"
    ident.2 = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    ident.1 = "B-2"
    ident.2 = "B-1"
  }
  subout = file.path(eout, source)
  if (!dir.exists(subout)) dir.create(subout)
  sub = subset(mdata, orig.ident %in% samples)
  compare_list = combn(as.character(unique(sub$group)), 2, simplify = FALSE)
  sub1 = subset(mdata, Ly6c2 == "Ly6c2+")
  Idents(sub1) = "group"
  deg = FindMarkers(sub1, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0)
  deg = deg %>% filter(p_val_adj < 0.05)
  write.csv(deg, file.path(subout, paste("DEG", source, "Ly6c2+", "csv", sep = ".")), row.names = T, quote = F)
  sub2 = subset(mdata, Ly6c2 == "Ly6c2-")
  Idents(sub2) = "group"
  deg = FindMarkers(sub2, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0)
  deg = deg %>% filter(p_val_adj < 0.05)
  write.csv(deg, file.path(subout, paste("DEG", source, "Ly6c2-", "csv", sep = ".")), row.names = T, quote = F)
}

###### 新的修复基因 #####
out = "/mnt/beegfs/Research/Users/yangyupeng/scRNA/KY202205250064/2022-07-04/03.seurat/KY202205250064/result-2022-07-16"
mrout = file.path(out, "Monocyte_repair_gene")
if (!dir.exists(out)) dir.create(out)
if (!dir.exists(mrout)) dir.create(mrout)
genes = "AGO2, BCL2L11, BTG1, CCL5, CD274, CD74, CDKN1A, CFLAR, CXCL10, CYBB, DHX9, FLNA, FOSL2, FOXP1, GRK2, HGF, IL1B, ITGB1, ITGB2, ITPKB, KMT2C, KRIT1, LILRA5, MAPK14, MYH9, NFE2L2, PDE3B, PHIP, PML, PRKX, PTEN, PTGS2, RB1, RUNX1, SOD1, STAT3, SYK, TNFRSF1A, VEZF1, WAPL, WASF2, YWHAZ, ZMIZ1"
genes = unlist(strsplit(genes, ", "))
# mdata = subset(data, celltype == "Monocyte")
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    case = "A-2"
    control = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    case = "B-2"
    control = "B-1"
  }
  sub = subset(mdata, orig.ident %in% samples)
  subout = file.path(mrout, source)
  if (!dir.exists(subout)) dir.create(subout)
  lapply(genes, function(gene){
    gene = str_to_title(gene)
    pt.size = AutoPointSize(sub)
    vln = VlnPlot(sub, features = gene, group.by = "group", slot = "data")
    df = vln$data
    df$ident = ifelse(df$ident == case, "MI", "Sham")
    df$ident = factor(df$ident, levels = c("Sham", "MI"))
    compare_list = combn(as.character(unique(df$ident)), 2, simplify = FALSE)
    vln = ggplot(df, aes_string("ident", gene, fill = "ident")) + geom_violin(trim = TRUE, scale = "width")
    vln = vln + geom_jitter(height = 0, size = pt.size, show.legend = FALSE)
    vln = vln + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
    vln = vln + xlab("Identity") + ylab("Expression Level") + labs(title = gene) + NoLegend()
    vln = vln + stat_compare_means(
      comparisons = compare_list, method = "wilcox.test", 
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                         symbols = c("****", "***", "**", "*", "ns")))
    ggsave(file.path(subout, paste0("vln.", gene, ".png")), vln, bg = "white", width = 3)
    ggsave(file.path(subout, paste0("vln.", gene, ".pdf")), vln, bg = "white", width = 3)
  })
}

lout = file.path(out, "Ly6C2_repair_gene")
if (!dir.exists(lout)) dir.create(lout)
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    case = "A-2"
    control = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    case = "B-2"
    control = "B-1"
  }
  sub = subset(mdata, orig.ident %in% samples)
  subout = file.path(lout, source)
  if (!dir.exists(subout)) dir.create(subout)
  lapply(genes, function(gene){
    gene = str_to_title(gene)
    pt.size = AutoPointSize(sub)
    vln = VlnPlot(sub, features = gene, group.by = "group", slot = "data", split.by = "Ly6c2")
    df = vln$data
    df$ident = ifelse(df$ident == case, "MI", "Sham")
    df$ident = factor(df$ident, levels = c("Sham", "MI"))
    compare_list = combn(as.character(unique(df$ident)), 2, simplify = FALSE)
    vln = ggplot(df, aes_string("ident", gene, fill = "ident")) + geom_violin(trim = TRUE, scale = "width")
    vln = vln + geom_jitter(height = 0, size = pt.size, show.legend = FALSE)
    vln = vln + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~split)
    vln = vln + xlab("Identity") + ylab("Expression Level") + labs(title = gene) + NoLegend()
    vln = vln + stat_compare_means(
      comparisons = compare_list, method = "wilcox.test", 
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                         symbols = c("****", "***", "**", "*", "ns")))
    ggsave(file.path(subout, paste0("vln.", gene, ".png")), vln, bg = "white", width = 7)
    ggsave(file.path(subout, paste0("vln.", gene, ".pdf")), vln, bg = "white", width = 7)
  })
}

##qc图##
qcout = file.path(out, "QC")
if (!dir.exists(qcout)) dir.create(qcout)
p = VlnPlot(data, features = c("nFeature_RNA")) + labs("nGene") + NoLegend()
ggsave(file.path(qcout, "vlnplot.nGene.png"), p)
ggsave(file.path(qcout, "vlnplot.nGene.pdf"), p)
p = VlnPlot(data, features = c("nCount_RNA")) + labs("nUMI") + NoLegend()
ggsave(file.path(qcout, "vlnplot.nUMI.png"), p)
ggsave(file.path(qcout, "vlnplot.nUMI.pdf"), p)
p = VlnPlot(data, features = c("percent.mt")) + NoLegend()
ggsave(file.path(qcout, "vlnplot.mt.png"), p)
ggsave(file.path(qcout, "vlnplot.mt.pdf"), p)
for(source in c("BL", "BM")){
  if (source == "BL"){
    samples = c("BL01","BL02","BL03","BL04","BL05","BL06")
    case = "A-2"
    control = "A-1"
  }else{
    samples = c("BM01","BM02","BM03","BM04","BM05","BM06")
    case = "B-2"
    control = "B-1"
  }
  sub = subset(data, orig.ident %in% samples)
  sub$group = ifelse(sub$group == case, "MI", "Sham")
  sub$group = factor(sub$group, levels = c("Sham", "MI"))
  p = VlnPlot(sub, features = c("nFeature_RNA"), group.by = "group") + labs("nGene") + NoLegend()
  ggsave(file.path(qcout, "vlnplot.nGene_by_group.png"), p, width = 3)
  ggsave(file.path(qcout, "vlnplot.nGene_by_group.pdf"), p, width = 3)
  p = VlnPlot(sub, features = c("nCount_RNA"), group.by = "group") + labs("nUMI") + NoLegend()
  ggsave(file.path(qcout, "vlnplot.nUMI_by_group.png"), p, width = 3)
  ggsave(file.path(qcout, "vlnplot.nUMI_by_group.pdf"), p, width = 3)
  p = VlnPlot(sub, features = c("percent.mt"), group.by = "group") + NoLegend()
  ggsave(file.path(qcout, "vlnplot.mt_by_group.png"), p, width = 3)
  ggsave(file.path(qcout, "vlnplot.mt_by_group.pdf"), p, width = 3)
}

## boxplot ##
library(cowplot)
df <- data.frame(data$orig.ident,data$cluster,data$celltype,data$nCount_RNA)
colnames(df) <- c('times','cluster','celltype','nUMI')
p = ggplot(df,aes(x=celltype,y=nUMI)) + 
  geom_boxplot() + coord_flip() +
  xlab("Celltypes") + ylab("nUMI") + scale_fill_discrete(guide=FALSE) +
  #  scale_y_continuous(expand=c(0,0),limits=c(0,17000),breaks=c(0,5000,10000,15000),labels=c('0','5000','10000','>=15000')) + 
  theme_cowplot()
bout = file.path(out, "Stat_of_per_cluster")
if (!dir.exists(bout)) dir.create(bout)
ggsave(file.path(bout, "boxplot_of_per_celltype.png"),p,bg="white")
ggsave(file.path(bout, "boxplot_of_per_celltype.pdf"),p,bg="white")
p = ggplot(df,aes(x=cluster,y=nUMI)) + 
  geom_boxplot() + coord_flip() +
  xlab("Clusters") + ylab("nUMI") + scale_fill_discrete(guide=FALSE) +
  #  scale_y_continuous(expand=c(0,0),limits=c(0,17000),breaks=c(0,5000,10000,15000),labels=c('0','5000','10000','>=15000')) + 
  theme_cowplot()
ggsave(file.path(bout, "boxplot_of_per_cluster.png"),p,bg="white")
ggsave(file.path(bout, "boxplot_of_per_cluster.pdf"),p,bg="white")

## 单核细胞的颜色和总的uamp相同 ##
library(scales)
uout = file.path(out, "UMAP")
if (!dir.exists(uout)) dir.create(uout)
p = DimPlot(mdata, group.by = "cluster", reduction = "umap", raster=FALSE) +
  scale_color_manual(values = hue_pal()(9)[c(4, 8)])
ggsave(file.path(uout, "umap.Monocyte_cluster.png"), p, width = 9)
ggsave(file.path(uout, "umap.Monocyte_cluster.pdf"), p, width = 9)
p = DimPlot(mdata, group.by = "celltype", reduction = "umap", raster=FALSE) +
  scale_color_manual(values = hue_pal()(7)[4])
ggsave(file.path(uout, "umap.Monocyte_celltype.png"), p, width = 9)
ggsave(file.path(uout, "umap.Monocyte_celltype.pdf"), p, width = 9)
p = DimPlot(mdata, group.by = "orig.ident", reduction = "umap", raster=FALSE)
ggsave(file.path(uout, "umap.Monocyte_sample.png"), p, width = 9)
ggsave(file.path(uout, "umap.Monocyte_sample.pdf"), p, width = 9)

## celltype stat ##
sout = file.path(out, "Celltype_stat")
if (!dir.exists(sout)) dir.create(sout)
tb = table(data$celltype, data$orig.ident)
tatol = rowSums(tb)
tb = cbind(Celltype = rownames(tb), tb, Total = tatol)
tb = as.data.frame(tb)
write.csv(tb, file.path(sout, "celltype_stat.csv"), row.names = F, quote = F)
