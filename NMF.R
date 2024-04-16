#NMF共聚类----
library(tidyverse)
library(parallel)
library(parallelMap)
parallelStartSocket(cpus = 16)

load("./GSE49710/0.data/GSE49710_norm.Rda")
load("./GSE49710/0.data/clinical.Rda")
load("GSE49710/1.NMF聚类/aa_genes.rda")
nmf.input <- GSE49710_norm[aa_genes,]
ranks <- 2:9
library(NMF)
estim <- lapply(ranks, function(r){
  fit <- nmf(nmf.input, r, nrun = 10, seed = 4) #nrun设置为5以免运行时间过长
  list(fit = fit, consensus = consensus(fit), .opt = "vp",coph = cophcor(fit))
})
names(estim) <- paste('rank', ranks)

estim.r <- nmf(as.matrix(nmf.input), 2:9, nrun=10,seed=4)

pdf("GSE49710/1.NMF聚类/cluster.pdf",width = 8,height =6.5,onefile=F)
plot(estim.r)+theme_bw(base_rect_size = 1.5)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = "right",legend.title = element_blank())
dev.off()

save(estim.r,file = "GSE49710/1.NMF聚类/NMF estim.r.Rda")
load("GSE49710/1.NMF聚类/NMF estim.r.Rda")

rank <- 3
seed <- 2
mut.nmf <- nmf(nmf.input,
               rank = rank,
               seed = seed)
group <- predict(mut.nmf) # 提出亚型

# consensus heatmap
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', 
               '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', 
               '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')
jco <- c("#0072B599","#E1872799","#BC3C2999",'#E95C59','#58A4C3',"#E31A1C","#1F78B4",'#B0997F','#93a8ac','#ffc857','#61a5c2','#119da4','#FF6666')#设置颜色
Cluster <- data.frame("cluster"=group[colnames(nmf.input)])
colnames(Cluster) <- 'Cluster'
Cluster$Cluster <- paste0('C',Cluster$Cluster)##前面加上C
table(Cluster$Cluster)
Cluster <- Cluster[,1,drop=F]

pdf("GSE49710/1.NMF聚类/cluster heatmap.pdf",width = 6,height =5,onefile=F)
consensusmap(mut.nmf,color = colorRampPalette(c('#CCE0F5', "white", "#B24745B2"))(50),
             labRow = NA,
             labCol = NA,
             annCol = Cluster,
             annColors = list(Cluster=c("C1"=jco[1],"C2"=jco[2],"C3"=jco[3]))
)
dev.off()


library(tinyarray)
dp = nmf.input[,order(group)]
draw_heatmap(dp,sort(group),
             color_an = c("#0072B599","#E1872799","#BC3C2999","#868686"),
             annotation_legend = T,
             cluster_cols = F,
             show_rownames = T)

library(survival)
library(survminer)
clin <- merge(clinical,Cluster,by=0)
colnames(clin)[4:5] <- c("OS_Time","OS_Event")
colnames(clin)[12] <- "Cluster"
colnames(clin)[2:3] <- c("EFS_Time","EFS_Event")
fit <- survfit(Surv(EFS_Time,EFS_Event) ~ Cluster,data = clin)
surv_pvalue(fit)$pval
summary(coxph(Surv(EFS_Time,EFS_Event)~Cluster,data=clin))

pdf("GSE49710/1.NMF聚类/cluster KM EFS.pdf",width = 4,height =4.3,onefile=F)
ggsurvplot(fit,data = clin,
           title="GSE49710",
           conf.int = T,
           pval = paste0('Log-rank\np = ',signif(surv_pvalue(fit)$pval,3)),
           pval.coord = c(1200,0.1), pval.size = 4,
           risk.table = F, risk.table.title = element_blank(),
           risk.table.fontsize =3,tables.height = 0.2,# risk.table.posS = "in",
           font.tickslab = 8,
           censor.shape = 124, censor.size = 2,
           legend.labs = c("C1","C2","C3"),
           palette = c("#0072B599","#E1872799","#BC3C2999"),
           # surv.median.line = "hv",
           ggtheme = theme_survminer(base_size = 10),
           tables.theme = theme_light(),
           # risk.table.y.col = T,
           legend.title=element_blank(),
           legend = c(0.1,0.15),
           font.main = c(14, "bold.italic"),
           font.x = c(12, "bold.italic", "darkred"),
           font.y = c(12, "bold.italic", "darkred")
)+xlab('Time (Days)')+ylab("Event Free Survival")#Event Free
dev.off()
# ggsave(filename ="GSE49710/1.NMF聚类/cluster KM OS.pdf")


save(clin,file="./GSE49710/0.氨基酸代谢/cluster_clin.Rda")




#临床特征----

library(ggsci)
load("./GSE49710/1.NMF聚类//cluster_clin.Rda")
clin$MYCN_status <- ifelse(clin$MYCN_status==1,"Amp","Non-Amp")
clin_p <- table(clin$MYCN_status,clin$Cluster)
chisq.test(clin_p)$p.value

clin_p <- table(clin$MYCN_status,clin$Cluster) %>% as.data.frame() %>% 
  set_names(.,c("MYCN_status","Cluster","Freq"))%>%
  mutate(labels=paste0("N = ",.$Freq))

ggplot(clin_p,aes(Cluster,Freq,fill=MYCN_status)) +
  geom_bar(stat = "identity",position = "fill",width = 0.4,alpha=1) +
  theme_bw()+ 
  scale_fill_manual(values = c("#BC3C2999","#62b3e9"))+
  # geom_text(aes(x=Cluster,label = labels),
  #           # position =position_stack(vjust = 0.5),size=2.5,
  #           show.legend = F,angle=90,color="black")+
  geom_text(aes(x=2,y=1.05,label="****",size=4),show.legend = F) +
  xlab("")+ylab("Percentage (%)")+
  ggtitle("GSE49710")+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = "right",legend.title = element_blank())

ggsave("GSE49710/1.NMF聚类/MYCN.pdf",width = 5.2,height = 4.3)




clin_p <- table(clin$INSS,clin$Cluster)
chisq.test(clin_p)$p.value

clin_p <- table(clin$INSS,clin$Cluster) %>% as.data.frame() %>% 
  set_names(.,c("INSS","Cluster","Freq"))%>%
  mutate(labels=paste0("N = ",.$Freq))

ggplot(clin_p,aes(Cluster,Freq,fill=INSS)) +
  geom_bar(stat = "identity",position = "fill",width = 0.4,alpha=1) +
  theme_bw()+ 
  scale_fill_manual(values = c("#62b3e9", "#829fcc", "#eed4b5","#BC3C2999",  "#727fae"))+
  # geom_text(aes(x=Cluster,label = labels),
  #           # position =position_stack(vjust = 0.5),size=2.5,
  #           show.legend = F,angle=90,color="black")+
  geom_text(aes(x=2,y=1.05,label="****",size=4),show.legend = F) +
  xlab("")+ylab("Percentage (%)")+
  ggtitle("GSE49710")+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = "right",legend.title = element_blank())

ggsave("GSE49710/1.NMF聚类/INSS.pdf",width = 5,height = 4.3)




