rm(list=ls());gc()
library(Seurat)
library(tidyverse)
library(export)
library(parallel)
library(parallelMap)
parallelStartSocket(cpus = detectCores())
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*12)
memory.limit(size = size_in_megabytes * 1024^12)

sce_dong <- readRDS("Dong/1.dong_seurat/data_dong_harmony.rds")


#############################################################
sce_dong_tumor <- subset(sce_dong,celltype=="tumor")
rm(sce_dong);gc()


sce_dong_tumor <- NormalizeData(sce_dong_tumor) %>%
  FindVariableFeatures() %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

ElbowPlot(sce_dong_tumor,ndims = 50)

sce_dong_tumor <- RunUMAP(sce_dong_tumor, reduction = "harmony", dims = 1:30)
sce_dong_tumor <- RunTSNE(sce_dong_tumor, reduction = "harmony", dims = 1:30)
sce_dong_tumor <- FindNeighbors(sce_dong_tumor, reduction = "harmony", dims = 1:30)
sce_dong_tumor <- FindClusters(sce_dong_tumor, resolution = .1)

DimPlot(sce_dong_tumor,reduction = "umap",group.by = "sample")

library(Scissor)
load("bulk data/cluster_clin.Rda")
load("bulk data//GSE49710_norm.Rda")
clin_tumor <- clin[,"Cluster",drop=F] %>%
  mutate(Cluster=ifelse(.$Cluster=="C3","1","0"))

expr_tumor <- as.data.frame(GSE49710_norm) %>% .[,rownames(clin_tumor)] %>%as.matrix()
identical(rownames(clin_tumor),colnames(expr_tumor))


phenotype = as.vector(clin_tumor$Cluster)


info_tumor <- Scissor(bulk_dataset=expr_tumor,sc_dataset =sce_dong_tumor,
                      phenotype=phenotype,
                      alpha =,
                      # cutoff=0.9,
                      tag=c("nonC3","C3"),
                      family = "binomial",
                      Load_file = "Dong/2.tumor cluster/2.1.scissor/Scissor_inputs.Rdata")

names(info_tumor)
length(info_tumor$Scissor_pos);info_tumor$Scissor_pos[1:4]
length(info_tumor$Scissor_neg);info_tumor$Scissor_neg[1:4]

Scissor_select <- rep("Scissor", ncol(sce_dong_tumor))
names(Scissor_select) <- colnames(sce_dong_tumor)
Scissor_select[info_tumor$Scissor_pos] <- "Scissor+"
Scissor_select[info_tumor$Scissor_neg] <- "Scissor-"
sce_dong_tumor <- AddMetaData(sce_dong_tumor, metadata = Scissor_select, col.name = "scissor")
DimPlot(sce_dong_tumor, reduction = 'umap',
        group.by = 'scissor',raster=FALSE,
        cols = c('grey','royalblue','indianred1'),
        pt.size = .5#, order = c("Scissor+","Scissor-")
)

table(sce_dong_tumor$scissor)


x <- sce_dong_tumor@active.ident %>%as.data.frame()
names(x) <- 'Type'
x_data <- sce_dong_tumor@reductions[["umap"]]@cell.embeddings
y <- as.data.frame(x_data)
mydata <- merge(x,y,by = 0,all.x = T)

celltype_position <- mydata %>%
  group_by(Type) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2))

ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,col= Type))+      
  geom_point(size=1)+
  labs(x = 'UMAP1',y = 'UMAP2')+
  theme_bw(base_rect_size = 1)+
  guides(colour=guide_legend(override.aes = list(size=2,alpha=1)))+
  theme(axis.text = element_text(size = 12,colour = 'black',angle = 0),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        # axis.title.x = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.text = element_text(size = 10),
        legend.position = "bottom",legend.title = element_blank())

graph2pdf(file="Dong/2.tumor cluster/2.1.scissor/Scissor UMAP.pdf",width = 5,height = 5.3)






library(clusterProfiler)
library(ReactomePA)
Idents(sce_dong_tumor) <- sce_dong_tumor$scissor
sce_dong_tumor.markers <- FindAllMarkers(object = sce_dong_tumor,
                                   only.pos = T,
                                   min.pct = 0.25,
                                   logfc.threshold = 0.25
)

marker.sig <- sce_dong_tumor.markers %>% filter(p_val_adj < 0.05, avg_log2FC >.5)

p_kegg <- data.frame()
for (i in  c("Scissor+")) {#,"Scissor-","Scissor"
  C <- subset(marker.sig ,cluster==i)
  de <- bitr(C$gene,"SYMBOL","ENTREZID",'org.Hs.eg.db') %>% merge(.,C,by.x="SYMBOL",by.y="gene")
  enrich_kegg <- enrichKEGG(de$ENTREZID,pvalueCutoff = 0.05 )
  p <- data.frame(enrich_kegg@result)  %>% data.frame(.,cluster=rep(i,nrow(.))) 
  p_kegg <- rbind(p_kegg,p)
  
}

dotplot(enrich_kegg,showCategory = 20)+
  scale_color_gradient(low ="#C91212" ,high = "#1F78B4")+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'darkred',face='bold'),
        # axis.title.x = element_blank(),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'right')
graph2pdf(file="Dong/2.tumor cluster/2.1.scissor/Scissor+ kegg.pdf",width = 6,height = 8.5)

sce_dong_tumor_scissor <- sce_dong_tumor
save(sce_dong_tumor_scissor,file = "Dong/2.tumor cluster/2.1.scissor/sce_dong_tumor_scissor.Rda")



#Scissor CytoTRACE----
rm(list=ls());gc()
library(CytoTRACE)
library(ggpubr)
library(ggsci)
load("Dong/2.tumor cluster/2.1.scissor/sce_dong_tumor_scissor.Rda")

mat_3k <- as_matrix(sce_dong_tumor_scissor@assays$RNA@counts)
mat_3k[1:4,1:4]
results <- CytoTRACE(mat = mat_3k,ncores = 1)
sce_dong_tumor_scissor_data <- AddMetaData(sce_dong_tumor_scissor,
                                     results$CytoTRACE,
                                     col.name = "CytoTRACE")
metadata <- sce_dong_tumor_scissor_data@meta.data
metadata$scissor <- factor(metadata$scissor)
colnames(metadata)
library(PupillometryR) 
library(ggpubr)
library(gghalves)
ggplot(metadata, aes(x = scissor, y = CytoTRACE, fill = scissor)) +
  geom_half_violin(aes(fill = scissor),position = position_nudge(x = -0.1, y = 0),  trim = TRUE, alpha = .5, colour = NA)+
  # geom_point(aes(x = .55, y = CytoTRACE, colour = scissor),position = position_jitter(width = .05), size = .5, shape = 20)+
  geom_boxplot(aes(x = scissor, y = CytoTRACE, fill = scissor),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_manual(values =  c("#868686FF","#0072B599","#BC3C2999" ))+
  scale_fill_manual(values = c("#868686FF","#0072B599","#BC3C2999" ))+
  stat_compare_means(aes(label=..p.signif..),#label.x.npc = .5,
                     tip.length = 0,size=5,vjust=.4,
                     comparisons = list(c("Scissor+","Scissor-")))+
  # geom_signif(map_signif_level = TRUE,test = "wilcox.test")+
  theme(axis.title = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12))+
  # ylim(0,1)+
  ylab('CytoTRACE score')+
  theme_bw(base_rect_size = 1)+
  theme(axis.text.x = element_text(size = 12,face = "bold",colour = 'black'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'none',
        axis.title.x = element_blank()
  )

graph2pdf(file="Dong/2.tumor cluster/2.2.CytoTRACE/Scissor UMAP CytoTRACE.pdf",width = 4.3,height = 4)


#######
rm(list=ls());gc()
load("Dong/2.tumor cluster/2.1.scissor/sce_dong_tumor_scissor.Rda")

library(beyondcell)
library(Seurat)
library(ggsci)
library(RColorBrewer)
library(scales)
library(tidyverse)

set.seed(1)
sc <- sce_dong_tumor_scissor
rm(sce_dong_tumor_scissor);gc()


gs <- GetCollection(SSc)
nopath <- GetCollection(PSc, include.pathways = FALSE)

bc <- bcScore(sc, gs, expr.thres = 0) 
bc <- bcUMAP(bc,pc=10, k.neighbors = 20, res = 0.1)
bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.1", pt.size = 1.5)

bc <- bcRanks(bc, idents = "bc_clusters_res.0.1", extended = FALSE)
bc <- bcRanks(bc, idents = "scissor", extended = FALSE)

head(bc@ranks$scissor)

bc <- readRDS("Dong/2.tumor cluster/2.4.beyondcell/bc.RDS")
bcClusters(bc, UMAP = "beyondcell", idents = "bc_clusters_res.0.1", pt.size = 1)+
  theme_bw(base_rect_size = 1)+
  xlab("UMAP1")+ylab("UMAP2")+
  scale_color_manual(values = pal_nejm("default",0.8)(8)[3:8])+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text(size = 10)
        # axis.title.x = element_blank()
  )
graph2pdf(file="Dong/2.tumor cluster/2.4.beyondcell/Scissor UMAP bc_clusters_res.0.1.pdf",width = 4.7,height = 4)



FindDrugs(bc, "CRIZOTINIB")
Genes <-bcSignatures(bc, UMAP = "beyondcell", genes = list(values = c("ALK","BIRC5")), pt.size = 1)
Genes[[1]]  +theme_bw(base_rect_size = 1)+
  scale_color_gradient(low = "gray90",high = "#d53f4d")+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'right',
        # legend.title = element_blank(),
        legend.text = element_text(size = 10)
        # axis.title.x = element_blank()
  )

graph2pdf(file="Dong/2.tumor cluster/2.4.beyondcell/ALK UMAP Beyondcell.pdf",width = 4.7,height = 4.5)




bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = c("sig-20955")), 
             pt.size = .5)+theme_bw()+theme(title = element_text(size=10))+
  theme_bw(base_rect_size = 1)+
  labs(x="UMAP1",y="UMAP2")+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=13,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'right',
        # legend.title = element_blank(),
        legend.text = element_text(size = 10)
        # axis.title.x = element_blank()
  )

graph2pdf(file="Dong/2.tumor cluster/2.4.beyondcell/sig-20955 UMAP Beyondcell.pdf",width = 4.7,height = 4.7)
bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = c("sig-21335")), 
             pt.size = .5)+theme_bw()+theme(title = element_text(size=10))+
  theme_bw(base_rect_size = 1)+
  labs(x="UMAP1",y="UMAP2")+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=13,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'right',
        # legend.title = element_blank(),
        legend.text = element_text(size = 10)
        # axis.title.x = element_blank()
  )

graph2pdf(file="Dong/2.tumor cluster/2.4.beyondcell/sig-21335 UMAP Beyondcell.pdf",width = 4.7,height = 4.5)


bcSignatures(bc, UMAP = "beyondcell", signatures = list(values = c("sig-20955","sig-21335")), 
             merged = "direct",pt.size = .5)+
  theme_bw(base_rect_size = 1)+
  labs(x="UMAP1",y="UMAP2")+
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=13,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = 'right',
        # legend.title = element_blank(),
        legend.text = element_text(size = 10)
        # axis.title.x = element_blank()
  )

graph2pdf(file="Dong/2.tumor cluster/2.4.beyondcell/sig-20955 sig-21335 UMAP Beyondcell.pdf",width = 4.7,height = 4.5)


