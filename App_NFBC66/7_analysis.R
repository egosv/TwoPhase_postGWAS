### Note that to generate these tables/plots all previous analyses must have been performed across all genes/phase 2 sample sizes
setwd("/out/path/processed_data/")

library(data.table)
library(ggplot2)
library(scales) 
library(RColorBrewer)
# install.packages("https://cran.r-project.org/src/contrib/UpSetR_1.4.0.tar.gz", type="source")
library(UpSetR)
library(reshape)
library(vcd)
library(gridExtra)

### Load files ----
resOpt_tot <- data.frame()
for (file in grep(paste0("NFBC66_TG_twophase_(LPL|GCKR|APOA5)_.+[.]RData"),dir(),value=T) ) { 
  load(file=paste0(file))
  GENE <- strsplit(file,split = "_|=")[[1]][4]
  sf <- strsplit(file,split = "_|=|[.]")[[1]][6]
  ss = c(449,1123,2246)[sf]
  resOpt_tot <- rbind(resOpt_tot,data.frame(GENE=GENE,ss=ss,scanresults))
}
rm(scanresults)
summary(resOpt_tot)
resOpt_tot$p.value <- pchisq(resOpt_tot$S,1,lower.tail = FALSE)

table(resOpt_tot$Alloc,resOpt_tot$GENE)

head(resOpt_tot)

###* RV ----
resRV_tot <- data.frame()
for (file in grep(paste0("NFBC66_TG_twophase_RVanalysis_.+[.]RData"),dir(),value=T) ) {
  load(file=paste0("./",file))
  str <- strsplit(file,split = "_|=")[[1]]
  ss <- str[6]
  resRV_tot <- rbind(resRV_tot, data.frame(ss=ss,resRV))
}
rm(resRV)
resRV_tot$p.value <- pchisq(resRV_tot$LR,1,lower.tail = FALSE)
head(resRV_tot)
summary(resRV_tot)

## Table with main results ----
cols1 <- c("GENE","chromosome","position","ss","snp.name","Alloc","beta1","p.value")
tabres <- subset(resOpt_tot,snp.name%in%c("rs268","rs2266788","rs3135506"),select = cols1)
formt <- function(x,ndec){format(round(x, ndec), nsmall = ndec)}
tabres$res <- paste0(formt(tabres$beta1,3)," (",format(tabres$p.value,digits = 3,scientific = T),")") 
tabres$GENE <- factor(tabres$GENE,c("LPL","APOA5"))
tabres$ss <- as.numeric(as.character(tabres$ss))
tabres$Alloc <- factor(tabres$Alloc,c("Complete","LM","GA","RDS","TZL"))

tabres1 <- cast(data = tabres,formula = ss+GENE+chromosome+position+snp.name ~ Alloc,value = "res")
uniq_val <- function(val){
  tmp <- val
  for(i in 2:length(tmp)){
    if(val[i]==val[i-1]) tmp[i] <- NA
  }
  return( tmp)
}
tabres1$ss <- uniq_val(tabres1$ss)
tabres1$GENE <- uniq_val(tabres1$GENE)
tabres1$chromosome <- uniq_val(tabres1$chromosome)

tabres1


## RVs table ----
cols2 <- c("GENE","ss","Weight","RVthres","corSNP1vsRV","corSNP2vsRV","nRV","Alloc","beta1","var_beta1","p.value")
tabRVres <- resRV_tot[, cols2]
tabRVres$res <- paste0(formt(tabRVres$beta1,4)," (",formt(sqrt(tabRVres$var_beta1),3),") p=",format(tabRVres$p.value,digits = 2,scientific = F)) 
tabRVres$ss <- as.numeric(as.character(tabRVres$ss))
tabRVres$GENE <- factor(tabRVres$GENE,c("GCKR","LPL","APOA5"))
tabRVres$Alloc <- factor(tabRVres$Alloc,c("Complete","LM","GA","RDS","TZL"))

tabRVres1 <- cast(data = subset(tabRVres, Alloc!="LM"), formula = Weight+ss+RVthres+nRV+GENE ~ Alloc, value = "res")
tabRVres2 <- subset(tabRVres, Alloc!="LM" & !is.na(corSNP1vsRV))[,c("Weight","ss","RVthres","nRV","GENE","corSNP1vsRV","corSNP2vsRV")]
tabRVres1 <- merge(tabRVres1, tabRVres2, sort=F)
tabRVres1$ss <- uniq_val(tabRVres1$ss)
tabRVres1$GENE <- uniq_val(tabRVres1$GENE)
tabRVres1$Weight <- uniq_val(tabRVres1$Weight)
tabRVres1$RVthres <- uniq_val(tabRVres1$RVthres)

tabRVres1


###* Beta-Beta plots ----
{
  loccols <- c(brewer.pal(n = 4, name = 'Set1')[c(2:4)])
  
  merged_data =  merge( subset(resOpt_tot, !Alloc %in% c("Complete")), subset(resOpt_tot, Alloc=="Complete"), by=c("GENE","ss","p","snp.name","position","allele.1", "allele.2") )
  merged_data$ss <- as.numeric(as.character(merged_data$ss))
  merged_data$Alloc.x <- factor(merged_data$Alloc.x, levels=c("LM","GA","RDS","TZL"))
  merged_data$GENE <- factor(merged_data$GENE, levels=c("GCKR","LPL","APOA5"))
  
  bb_plot <- ggplot(merged_data, aes(x = beta1.x, y=beta1.y)) + 
    geom_abline(intercept=0, slope=1, colour="red", size=0.9) +
    geom_point(aes(col=factor(ss)), size=2) + 
    facet_grid(GENE ~ Alloc.x) + 
    scale_colour_manual(values=loccols) + 
    guides(color=guide_legend(title="n")) + 
    xlab(expression(SPML)) + ylab(expression(Complete)) + 
    theme_bw(base_size = 25) + 
    theme(axis.text.x = element_text(vjust=0.5,angle=0,size=14), 
          axis.text.y  = element_text(vjust=0.5, size=14), 
          axis.title.x = element_text(face="bold",size=17,vjust=-0.5), 
          axis.title.y = element_text(face="bold",size=17), 
          strip.text.x = element_text(size=17, face="bold",color="#FFFFFF"), 
          strip.text.y = element_text(size=17, face="bold",color="#FFFFFF",angle=270), 
          strip.background = element_rect(fill="black"), 
          legend.position="bottom", 
          legend.title=element_text(size=17), legend.text=element_text(size=20),
          legend.key = element_rect(color="black",size=0.5)) 
  
  ggsave(bb_plot, file=paste0("App_beta-betaplot.png"), width=11, height=8.5, units="in", bg="transparent", dpi=500)

}

###* Region plots ----
{
  ### select data
  regplot_data =  subset(resOpt_tot)
  regplot_data$Alloc = factor(regplot_data$Alloc, levels=c("Complete","LM","GA","RDS","TZL"))
  regplot_data$ss <- as.numeric(as.character(regplot_data$ss))
  
  cbPalette <- c("black",brewer.pal(n = 4, name = 'RdYlGn'))
  
  regplot_data$mlog10p <- -log10(regplot_data$p.value)
  regplot_data$Mbp <- regplot_data$position/1000000
  regplot_data$GENE1 <- factor(regplot_data$GENE,levels = c("GCKR","LPL","APOA5"))
  
  regplot <- ggplot(regplot_data, aes(y=mlog10p,x=Mbp)) + 
    geom_point(aes(color=Alloc,shape=Alloc),size=2.5,stroke=1.2) + 
    facet_grid(ss ~ GENE1, scales = "free_x") + 
    scale_color_manual(values=cbPalette,labels=expression("Complete","LM","GA","RDS","TZL")) + 
    scale_shape_manual(values=c(1,0,5,2,6),labels=expression("Complete","LM","GA","RDS","TZL")) + 
    guides(color=guide_legend(title=""),shape=guide_legend(title=""))  + 
    theme_bw(base_size = 11) + ylab(expression(-log["10"]*(p))) + xlab("Genomic position in Mbps (hg19)") + theme(axis.text.x = element_text(vjust=0.5,angle = 0, size=8), axis.text.y  = element_text(vjust=0.5, size=10), axis.title.x = element_text(face="bold",size=15,vjust=-0.5), axis.title.y = element_text(face="bold",size=15), strip.text.x = element_text(size=12, face="bold",color="#FFFFFF"), strip.text.y = element_text(size=14, face="bold",color="#FFFFFF"), strip.background = element_rect(fill="black"), legend.position="bottom", legend.title=element_text(size=17), legend.text=element_text(size=13),legend.key = element_rect(color="black",fill="transparent"),panel.background=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.background=element_blank(),legend.background=element_blank(),plot.title = element_text(lineheight=.8, face="bold",size=17))
  
  ggsave(regplot, file=paste0("RegionPlot_App.png"), width=8, height=4, units="in", bg="transparent", dpi=500)
  
} 

###* Mosaic plots ----

### load data ---
ph2ind <- data.frame()
for (file in grep(paste0("Phase2_indicators_sf=.+[.]RData"),dir(),value=T) ) { 
  load(file=file)
  data_phase2 <- cbind(Phase1_data_nonmiss_wCTS,R1=R1,R2=R2,R3=R3,R4=R4)
  str <- strsplit(file,split = "_|=")[[1]]
  sf <- str[4]
  ss = c(449,1123,2246)[as.numeric(sf)]
  ph2ind <- rbind(ph2ind,data.frame(ss=ss, data_phase2))
}

df.mosaic = reshape(ph2ind, varying=list(paste0("R",1:4)), times=c("RDS","TZL","LM","GA"), timevar = "design", direction="long")
names(df.mosaic)[names(df.mosaic)=="R1"] <- "R"
df.mosaic$design <- factor(df.mosaic$design,levels=c("LM","GA","RDS","TZL"))
head(df.mosaic[df.mosaic$R==1,])
summary(df.mosaic)


myColors = grey.colors(9)

vnames <- list(set_varnames = c(design="",S="Yst"))
lnames <- list(design = c("LM","GA","RDS","TZL"), S=c(expression(T["1"]),expression(T["2"]),expression(T["3"])))
lnames0 <- list(design = rep("",4), StrataZ=rep("",9), S=c(expression(T["1"]),expression(T["2"]),expression(T["3"])))
vnames0 <- list(set_varnames = c(design="",S="Yst",StrataZ=""))

### Mosaic all
p0 = grid.grabExpr(mosaic(StrataZ~S, data=subset(df.mosaic,ss==2246), labeling_args=vnames, set_labels=lnames, zero_size = 0, main = "", margins=unit(c(2.5,3,0.5,3), "lines"), highlighting_fill = myColors,split_vertical=TRUE, main_gp = gpar(fontsize = 15)))

png(paste0("Mosaics_App_ph1.png"),width = 7, height = 4, units = "in", bg="transparent", res=300)
grid.draw(p0)
dev.off()

### Plot comparing designs across ph2 sample sizes
p2 = grid.grabExpr(mosaic(StrataZ~S|design, data=xtabs(R ~ StrataZ + S + design, data=droplevels(subset(df.mosaic,ss==1123 & design!="LM"))), labeling_args=vnames0, set_labels=lnames, zero_size = 0, main = "1123", margins=unit(c(2.5,3,0.5,3), "lines"), highlighting_fill = myColors, main_gp = gpar(fontsize = 15), rot_labels=c(0,0,90,0)))

p3 = grid.grabExpr(mosaic(StrataZ~S|design, data=xtabs(R ~ StrataZ + S + design, data=droplevels(subset(df.mosaic,ss==2246 & design!="LM"))), labeling_args=vnames, set_labels=lnames, zero_size = 0, main = "2246", margins=unit(c(2.5,3,0.5,3), "lines"), highlighting_fill = myColors, main_gp = gpar(fontsize = 15), rot_labels=c(0,0,90,0)))


obj <- grid.arrange(p2, p3, ncol=2)

png(paste0("Mosaics_App.png"),width = 7, height = 4, units = "in", bg="transparent", res=300)
grid.draw(obj)
dev.off()
  
###* UpSet plot ----
{
  # subset(inter_tot_sim2,it==1 & ID==10)
  desgs <- c("RDS","TZL","GA")
  
  dt1 <- subset(ph2ind, ss==1123)[,c(2,13,14,16)]
  names(dt1)[2:4] <- paste(desgs,"1123", sep="\t")
  
  dt2 <- subset(ph2ind, ss==2246)[,c(2,13,14,16)]
  names(dt2)[2:4] <- paste(desgs,"2246", sep="\t")
  
  # check if rows match
  # all.equal(dt1$FID,dt2$FID)
  
  inter_app <- cbind(dt1,dt2[,-1])
  
  png(paste0("Upset_App.png"),width = 8, height = 5, units = "in", bg="transparent", res=300)
  print(upset(inter_app, nintersects=NA, sets = c(paste(rep(desgs,2), rep(c("1123","2246"),each=3),  sep="\t")), order.by = "freq", group.by = "degree", set_size.show=FALSE, show.numbers=TRUE))
  dev.off()
  
  postscript(paste0("Upset_App.eps"), width = 8, height = 5, bg="transparent", horizontal=FALSE, onefile = FALSE, paper="special")
  print(upset(inter_app, nintersects=NA, sets = c(paste(rep(desgs,2), rep(c("1123","2246"),each=3),  sep="\t")), order.by = "freq", group.by = "degree", set_size.show=FALSE, show.numbers=TRUE))
  dev.off() 
  
}

