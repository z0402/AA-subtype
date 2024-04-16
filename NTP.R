
rm(list=ls())
load("E-MTAB-8248/data/E_MTAB_8248_expr_cl.Rda")
load("/DEG.Rda")

gene_C1 <- rownames(subset(deg1,logFC> 1 & adj.P.Val < 0.05))  
gene_C2 <- rownames(subset(deg2,logFC> 1 & adj.P.Val < 0.05)) 
gene_C3 <- rownames(subset(deg3,logFC> 1 & adj.P.Val < 0.05)) 

Cluster <- rep(c('C1','C2','C3'),c(376,138,217))
mydata <- data.frame(Cluster,c(gene_C1,gene_C2,gene_C3))
colnames(mydata) <- c('Cluster','Gene')  

col <- c('#B0997F','#93a8ac','#ffc857','#61a5c2','#119da4','#FF6666')


w <- intersect(rownames(E_MTAB_8248_exp),mydata$Gene)
ee <- E_MTAB_8248_exp[w,]
emat <- t(scale(t(ee)))
templates <-mydata[mydata$Gene%in%w,]
colnames(templates) <- c('class','probe')

library(CMScaller)
res <- ntp(emat      = emat,
           templates = templates,
           doPlot    = F,
           nPerm     = 1000,
           distance  = "cosine",
           nCores    = 1,
           seed      = 2020104,
           verbose   = T)
subHeatmap(emat = emat, res = res, templates = templates,
           classCol = c("#E31A1C","#1F78B4",'#119da4'),heatCol = c('seagreen','orange'))

pdf("E-MTAB-8248/cluster heatmap.pdf",width = 5,height =5.5,onefile=F)
subHeatmap(emat = emat, res = res, templates = templates,
           classCol =c("#0072B599","#E1872799","#BC3C2999","#868686"),
           heatCol = c('seagreen','orange'))+title("E-MTAB-8248")
dev.off()


table(res$FDR<0.05)        
res2 <- res[res$FDR<0.05,]   
table(res2$prediction) 

cl0 <- merge(E_MTAB_8248_cl,res2,by=0)
colnames(cl0)[14] <- 'Cluster'


library(survival)
library(survminer)
fit <- survfit(Surv(OS_Time,OS_Event) ~ Cluster,data = cl0)
surv_pvalue(fit)$pval
summary(coxph(Surv(OS_Time,OS_Event)~Cluster,data=cl0))

mytheme <- theme_survminer(font.legend = c(14,"plain", "black"),
                           font.x = c(14,"plain", "black"),
                           font.y = c(14,"plain", "black"))
pdf("E-MTAB-8248/cluster KM OS.pdf",width = 4,height =4.3,onefile=F)

ggsurvplot(fit,data = cl0,
           title="E-MTAB-8248",
           conf.int = T,
           pval = paste0('Log-rank\np = ',signif(surv_pvalue(fit)$pval,3)),
           pval.coord = c(1100,0.1), pval.size = 4,
           risk.table = F, risk.table.title = element_blank(),
           risk.table.fontsize =3,tables.height = 0.2,# risk.table.pEFS = "in",
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
)+xlab('Time (Days)')+ylab("Overall Survival")#Event Free  Overall
dev.off()


#
cl0$MYCN <- ifelse(cl0$MYCN =="amplified","Amp",
                   ifelse(cl0$MYCN =="non amplified","Non-Amp",""))
cl0$`1p` <- ifelse(cl0$`1p` =="del/im","1p deletion",
                   ifelse(cl0$`1p` =="normal","1p non-deletion",""))
cl0$ALT <- ifelse(cl0$ALT =="negative","Negative","Positive")


library(ggsci)
library(tidyverse)
cl0_p <- table(cl0$MYCN,cl0$Cluster)
fisher.test(cl0_p)$p.value
cl0_1 <- cl0[cl0$MYCN !="",]
cl0_p <- table(cl0_1$MYCN,cl0_1$Cluster) %>% as.data.frame() %>% 
  set_names(.,c("MYCN","Cluster","Freq"))%>%
  mutate(labels=paste0("N = ",.$Freq))

ggplot(cl0_p,aes(Cluster,Freq,fill=MYCN)) +
  geom_bar(stat = "identity",position = "fill",width = 0.4,alpha=1) +
  theme_bw()+ 
  scale_fill_manual(values =c("#BC3C2999","#62b3e9"))+

  geom_text(aes(x=2,y=1.05,label="****",size=5),show.legend = F) +
  xlab("")+ylab("Percentage (%)")+ggtitle("E-MTAB-8248")+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = "right",legend.title = element_blank())

ggsave("E-MTAB-8248//MYCN.pdf",width = 5.2,height = 4.3)




cl0_p <- table(cl0$INSS,cl0$Cluster)
fisher.test(cl0_p,simulate.p.value = T)

cl0_p <- table(cl0$INSS,cl0$Cluster) %>% as.data.frame() %>% 
  set_names(.,c("INSS","Cluster","Freq"))%>%
  mutate(labels=paste0("N = ",.$Freq))

ggplot(cl0_p,aes(Cluster,Freq,fill=INSS)) +
  geom_bar(stat = "identity",position = "fill",width = 0.4,alpha=1) +
  theme_bw()+ 
  scale_fill_manual(values =c("#62b3e9", "#829fcc", "#eed4b5","#BC3C2999",  "#727fae"))+
  # geom_text(aes(x=Cluster,label = labels),
  #           # position =position_stack(vjust = 0.5),size=2.5,
  #           show.legend = F,angle=90,color="black")+
  geom_text(aes(x=2,y=1.05,label="***",size=5),show.legend = F) +
  xlab("")+ylab("Percentage (%)")+ggtitle("E-MTAB-8248")+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = "right",legend.title = element_blank())

ggsave("E-MTAB-8248/INSS.pdf",width = 5,height = 4.3)



cl0_2 <- cl0[cl0$`1p` !="",]
cl0_p <- table(cl0_2$`1p`,cl0_2$Cluster)
chisq.test(cl0_p)$p.value

cl0_p <- table(cl0_2$`1p`,cl0_2$Cluster) %>% as.data.frame() %>% 
  set_names(.,c("1p","Cluster","Freq"))%>%
  mutate(labels=paste0("N = ",.$Freq))

ggplot(cl0_p,aes(Cluster,Freq,fill=`1p`)) +
  geom_bar(stat = "identity",position = "fill",width = 0.4,alpha=1) +
  theme_bw()+ 
  scale_fill_manual(values =c("#BC3C2999","#62b3e9"))+
  # geom_text(aes(x=Cluster,label = labels),
  #           # position =position_stack(vjust = 0.5),size=2.5,
  #           show.legend = F,angle=90,color="black")+
  geom_text(aes(x=2,y=1.05,label="****",size=5),show.legend = F) +
  xlab("")+ylab("Percentage (%)")+
  ggtitle("E-MTAB-8248")+
  theme_bw(base_rect_size = 1.5)+
  theme(axis.text.y = element_text(size = 12,colour = 'black'),
        axis.text.x = element_text(size = 12,colour = 'darkred',face='bold'),
        axis.title = element_text(size = 14,colour = 'darkred',face='bold'),
        plot.title = element_text(size=15,hjust=0.5,face='bold'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"),
        panel.background = element_rect(fill='#f3f6f6'),
        legend.position = "right",legend.title = element_blank())

ggsave("E-MTAB-8248//1p.pdf",width = 5.5,height = 4.3)







