#
setwd("F:/Dropbox/Scripts/16S_Copy_Number/Version3/")
source("fitPE_functions.R")
setwd("F:/Dropbox/Scripts/16S_Copy_Number/Version4/")
source("RasperGade_reconstruction.R")
source("RasperGade_diagnosis.R")
setwd("F:/Dropbox/Scripts/16S_Copy_Number/20210302/")
library(ggplot2)
library(ggpubr)
#
panelA.phy.df = data.frame(x=c(-1,-0.5,0,0.5,1),y=c(-1,-0.5,0,-0.5,-1),
                           label=c("x[1]*\"  \"*v[1]","l[1]","x[0]","l[2]","x[2]*\"  \"*v[2]"),
                           color=c("b","b","r","b","b"),
                           x.adjust=c(0,-0.12,0,0.12,0),y.adjust=c(-0.12,0.12,0.12,0.12,-0.12))
panelA.plot = ggplot()+
  geom_line(mapping = aes(x=x,y=y),data = panelA.phy.df,size=2)+
  geom_text(mapping = aes(x=x+x.adjust,y=y+y.adjust,label=label,color=color),
            data = panelA.phy.df,size=4,parse = TRUE)+
  coord_cartesian(xlim = c(-2,2),ylim = c(0.5,-1.5))+
  scale_color_manual(values = c(b="black",r="red"))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#
panelB.phy.df = data.frame(x=c(-1.5,-0.5,0.5,1.5,-1,1),
                           xend=c(-1,-1,1,1,0,0),
                           y=c(-2,-2,-2,-2,-1,-1),
                           yend=c(-1,-1,-1,-1,0,0),
                           x.adjust = c(0,0,0,0,-0.12,0.12),
                           y.adjust = c(-0.06,-0.06,-0.06,-0.06,0,0),
                           label=c("x[3]","x[4]","x[5]","x[6]","x[1]","x[2]"))
panelB.plot = ggplot()+
  geom_segment(mapping = aes(x=x,xend=xend,y=y,yend=yend),data = panelB.phy.df,size=2)+
  geom_point(mapping = aes(x=xend,y=yend),data=panelB.phy.df,size=1.25)+
  geom_text(mapping = aes(x=x+3*x.adjust,y=y+3*y.adjust,label=label),data = panelB.phy.df,size=4,parse = TRUE)+
  annotate(geom = "text",x = 0,y = 0.1,label="x[0]",size=4,parse=TRUE,color="red")+
  coord_cartesian(xlim = c(-2.5,2.5),ylim = c(0.25,-2.25))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#
panelC.phy.df = data.frame(x=c(-1.5,-0.5,0.5,-1,0),
                           xend=c(-1,-1,0,0,0),
                           y=c(-2,-2,0,-1,1),
                           yend=c(-1,-1,1,1,2),
                           x.adjust=c(0,0,0,0.12,-0.12),
                           y.adjust=c(-0.06,-0.06,-0.06,0,0),
                           label=c("x[3]","x[4]","x[5]","x[1]","x[2]"))
panelC.plot = ggplot()+
  geom_segment(mapping = aes(x=x,xend=xend,y=y,yend=yend),data = panelC.phy.df,size=2)+
  geom_point(mapping = aes(x=xend,y=yend),data=panelC.phy.df,size=1.25)+
  geom_text(mapping = aes(x=x+3*x.adjust,y=y+3*y.adjust,label=label),data = panelC.phy.df,size=4,parse = TRUE)+
  annotate(geom = "text",x = 0,y = 2.25,label="x[6]",size=4,parse=TRUE,color="red")+
  coord_cartesian(xlim = c(-2.5,2.5),ylim = c(2.25,-2.25))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#
combined.plot = ggarrange(plotlist = list(panelA.plot,panelB.plot,panelC.plot),nrow = 1,ncol = 3,
                          labels = c("A","B","C"))
#
ggsave(filename = "FigS4.pdf",device = "pdf",plot = combined.plot,path = "./",
       scale = 2.5,width = 3.25,height = 1.25,units = "in")
ggsave(filename = "FigS4.png",device = "png",plot = combined.plot,path = "./",
       scale = 2.5,width = 3.25,height = 1.25,units = "in",dpi = "print",type="cairo")
#