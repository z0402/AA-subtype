rm(list=ls());gc()
library(Seurat)
library(dorothea)
library(tidyverse)
library(viper)
library(parallel)
library(parallelMap)
parallelStartSocket(cpus = detectCores())
load("2.tumor cluster/2.1.scissor/sce_dong_tumor_scissor.Rda")

seu <- subset(sce_dong_tumor_scissor,scissor %in% c("Scissor-","Scissor+"))
rm(sce_dong_tumor_scissor);gc()


## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

regulon <- dorothea_regulon_human %>% 
  dplyr::filter(confidence %in% c("A","B"))

seu <- run_viper(seu, regulon,options = list(method = "scale", minsize = 4, 
                                               eset.filter = FALSE, cores = 1, verbose = FALSE))
Assays(seu)
saveRDS(seu,file="Dong/2.tumor cluster/2.5.dorothea/seu_dorothea.RDS")

seu <- readRDS("Dong/2.tumor cluster/2.5.dorothea/seu_dorothea.RDS")
DefaultAssay(seu) <- "dorothea"
Idents(seu) <- seu$scissor
table(Idents(seu))


seu.markers <- FindAllMarkers(object = seu, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               thresh.use = 0.25)                                                          
DT::datatable(seu.markers)
write.csv(seu.markers,file='Dong/2.tumor cluster/2.5.dorothea/seu.markers.csv')


library(dplyr) 
seu.markers$fc = seu.markers$pct.1 - seu.markers$pct.2
top10 <- seu.markers %>% group_by(cluster) %>% top_n(10, fc)
seu@assays$dorothea@data[1:4,1:4]
top10 =top10[top10$fc > 0.3,] 
seu <- ScaleData(seu)
DoHeatmap(seu,top10$gene,size=3,slot = 'scale.data')


viper_scores_df <- GetAssayData(seu, slot = "scale.data", assay = "dorothea") %>%
  data.frame(check.names = F) %>% t()
viper_scores_df[1:4,1:4]


## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(seu)), 
                            cell_type = as.character(Idents(seu)),
                            check.names = F)
head(CellsClusters)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)
head(viper_scores_clusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

head(summarized_viper_scores)


## We select the 20 most variable TFs.
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(20, var) %>%
  distinct(tf)
highly_variable_tfs

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) %>% t() %>%
  as.data.frame()


ann_colors = list(
  CellType = c(`Scissor-` = "#0072B599",
               `Scissor+` =  "#BC3C2999"
  )
)
celltype <- data.frame(CellType=c("Scissor-","Scissor+"))
rownames(celltype) <- c("Scissor-","Scissor+")

library(pheatmap)
pdf("Dong/2.tumor cluster/2.5.dorothea/Scissor+ TFs.pdf",width = 5,height =4)

pheatmap(summarized_viper_scores_df,fontsize=12, 
        fontsize_row = 8, show_colnames = F,
        annotation_col=celltype,
        annotation_colors = ann_colors ,
        color= colorRampPalette( c('gray90','darkred') )(10),# breaks = my_breaks, 
        # main = "DoRothEA (ABC)", 
        # angle_col = 45,
        treeheight_col = 0, treeheight_row = 0, border_color = "white") 
dev.off()







