setwd("~/Downloads/")

#import data from mothur
otus<-as.data.frame(import_mothur(mothur_shared_file ="biofilm.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"))
taxa <- as.data.frame(import_mothur(mothur_constaxonomy_file =("biofilm.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy")))
meta <-read.csv("biometa.csv", header=TRUE, row.names=1, check.names = FALSE)


#Format
MET=sample_data(meta)
OTU=as.matrix(otus)
colnames(taxa) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX=tax_table(as.matrix(taxa))

otus <-read.csv("otus.csv", header=TRUE, row.names=1, check.names = FALSE)
meta_table <-read.csv("biometa.csv", header=TRUE, row.names=1, check.names = FALSE)

#Adapted script from https://link.springer.com/chapter/10.1007/978-981-13-1534-3_10

rdata<-otus
#transpose
rdata_t<-t(rdata)
#replace 0s with count multiplicative method w/ zCompositions package
rdata_r<-t(cmultRepl((rdata_t), method="CZM", output= "p-counts"))
#convert to proportions
rdata_p<-apply(rdata_r, 2, function(x){x/sum(x)})
#clr transform
names_add<- rownames(rdata_pf)[ + order(apply(rdata_pf, 1, sum), decreasing=T)]
rdata_pr<-rdata_pf[names_add,]
rdata_clr<-t(apply(rdata_pr, 2, function(x){log(x)-mean(log(x))}))
#Singular Value Decomposition
rdata_PCX<- prcomp(rdata_clr)
sum(rdata_PCX$sdev[1:2]^2)/mvar(rdata_clr)
#create distance matrix and cluster data
dist<-dist(rdata_clr, method="euclidian")
hc<-hclust(dist, method="ward.D2")
re_order<-rdata_pr[,hc$order]
abund_PCX<- prcomp(rdata_clr)

#Permanova looking at differences in the microbiota composition of samples from 3 different facilities
#experimental endpoint, and by treatment
rdata_dist<-dist(rdata_clr, method='euclidean')

permanovaoverall<-pairwise.adonis(rdata_dist, factors=meta_table$Facility, perm = 999, p.adjust.m = 'bonferroni')

permanovaodays<-pairwise.adonis(rdata_dist, factors=meta_table$Days, perm = 999, p.adjust.m = 'bonferroni')

permanovaotreat<-pairwise.adonis(rdata_dist, factors=meta_table$Treatment, perm = 999, p.adjust.m = 'bonferroni')


#differentially abundant taxa when comparing 3D to 5D biofilms
# meta table and reformatting
D<-read.csv("3Daldex.csv", row.names=1, check.names = FALSE)
D.row<-row.names(D)
alde35D<-as.data.frame(subset((t(rdata)), rownames(t(rdata)) %in% D.row))
alde35D_t<-t(alde35D)
#if 3d then 3d, if not then label 5d
groups <- with(D,ifelse(as.factor(Days)%in% c("3D"), c("3D"), c("5D")))

#clr on transposed data
vdr <- aldex.clr(alde35D_t, groups, mc.samples=128, verbose=FALSE)
#welchs t and wilcox rank sums test
vdr_t<-aldex.ttest(vdr, groups, paired.test = FALSE)
#calculate effect size
vdr_effect <- aldex.effect(vdr, groups, verbose=FALSE)
#merge the data into one data frame
vdr_all <- data.frame(vdr_t, vdr_effect)
# Identify significant taxa
sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
#reach significance when the p-values are adjusted for multiple testing corrections using the Benjamini-Hochberg’s method.
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)

#export data
D3v5D <-xtable(vdr_all[sig_by_both_fdr, c(1:11)], caption="Table of significant taxa", digits=3,label="sig.table", align=c("l",rep("r",11) ))
print.xtable(D3v5D, type="html", file="D3v5D.html")

#differentially abundant taxa when comparing 3D to 15D biofilms

D15<-read.csv("D315D.csv", row.names=1, check.names = FALSE)
D15.row<-row.names(D15)
alde35D<-as.data.frame(subset((t(rdata)), rownames(t(rdata)) %in% D15.row))
alde35D_t<-t(alde35D)

groups <- with(D15,ifelse(as.factor(Days)%in% c("3D"), c("3D"), c("15D")))

vdr <- aldex.clr(alde35D_t, groups, mc.samples=128, verbose=FALSE)
vdr_t<-aldex.ttest(vdr, groups, paired.test = FALSE)
vdr_effect <- aldex.effect(vdr, groups, verbose=FALSE)
vdr_all <- data.frame(vdr_t, vdr_effect)

sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)

D3v15D <-xtable(vdr_all[sig_by_both_fdr, c(1:11)], caption="Table of significant taxa", digits=3,label="sig.table", align=c("l",rep("r",11) ))
print.xtable(D3v15D, type="html", file="D3v15D.html")

#differentially abundant taxa when comparing 5D to 15D biofilms

D515<-read.csv("D515D.csv", row.names=1, check.names = FALSE)
D515.row<-row.names(D515)
alde55D<-as.data.frame(subset((t(rdata)), rownames(t(rdata)) %in% D515.row))
alde55D_t<-t(alde55D)

groups <- with(D515,ifelse(as.factor(Days)%in% c("5D"), c("5D"), c("15D")))

vdr <- aldex.clr(alde55D_t, groups, mc.samples=128, verbose=FALSE)
vdr_t<-aldex.ttest(vdr, groups, paired.test = FALSE)
vdr_effect <- aldex.effect(vdr, groups, verbose=FALSE)
vdr_all <- data.frame(vdr_t, vdr_effect)

sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)

D5v15D <-xtable(vdr_all[sig_by_both_fdr, c(1:11)], caption="Table of significant taxa", digits=3,label="sig.table", align=c("l",rep("r",11) ))
print.xtable(D5v15D, type="html", file="D5v15D.html")

# aldex for F1 v F2 3D
stomach <- subset_samples(carbom, Days=="3D")
D3F<-otu_table(stomach)
D3F1<-read.csv("D3F1.csv", row.names=1, check.names = FALSE)
D3F1.row<-row.names(D3F1)
alde3DF1<-as.data.frame(subset((t(D3F)), rownames(t(D3F)) %in% D3F1.row))
alde3DF1_t<-t(alde3DF1)

groups <- with(D3F1,ifelse(as.factor(Facility)%in% c("F1"), c("F1"), c("F2")))

vdr <- aldex.clr(alde3DF1_t, groups, mc.samples=128, verbose=FALSE)
vdr_t<-aldex.ttest(vdr, groups, paired.test = FALSE)
vdr_effect <- aldex.effect(vdr, groups, verbose=FALSE)
vdr_all <- data.frame(vdr_t, vdr_effect)

sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)

D3v15D <-xtable(vdr_all[sig_by_both_fdr, c(1:11)], caption="Table of significant taxa", digits=3,label="sig.table", align=c("l",rep("r",11) ))
print.xtable(D3v15D, type="html", file="D3v15D.html")


library(compositions)

#adapte from Laura's Script for overall pca

mvar_clr<-mvar(rdata_clr)

row_16sY1<-rownames(rdata_clr) #Make vector with sample names
pc_out_16sY1<-as.data.frame(rdata_PCX$x[,1:2]) #Get PC1 and PC2
pc_out_meta_16sY1<-as.data.frame(bind_cols(pc_out_16sY1,meta$Treatment,meta$Days)) #Add metadata information
row.names(pc_out_meta_16sY1)<-row_16sY1



# test out all categories but will use treatment as shape and days as color for most distinct clusters

PCA <- ggplot(pc_out_16sY1, aes(x=PC1,y=PC2, color=meta$Days, shape=meta$Treatment))+
  geom_point(size=3)+
  geom_text(aes(label=row.names(pc_out_meta_16sY1)),hjust=0, vjust=0) +
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(rdata_PCX$sdev[1]^2/mvar_clr*100, digits=2), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(rdata_PCX$sdev[2]^2/mvar_clr*100, digits=2), "%", sep="")) +
  ggtitle("Microbiota of 3 Apple Packing House Facility Attached Biomass")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='viridis')

#subset otu table in excel to seperate facilities

otus <-read.csv("F1otus.csv", header=TRUE, row.names=1, check.names = FALSE)
meta <-read.csv("biometaF1.csv", header=TRUE, row.names=1, check.names = FALSE)

rdata<-otus
rdata_t<-t(rdata)
rdata_r<-t(cmultRepl((rdata_t), method="CZM", output= "p-counts"))
rdata_p<-apply(rdata_r, 2, function(x){x/sum(x)})
rdata_pf<- rdata_r[apply(rdata_p, 1, min) > 0.00001,]
names_add<- rownames(rdata_pf)[ + order(apply(rdata_pf, 1, sum), decreasing=T)]
rdata_pr<-rdata_pf[names_add,]
rdata_clr<-t(apply(rdata_pr, 2, function(x){log(x)-mean(log(x))}))
rdata_PCX<- prcomp(rdata_clr)

row_16sY1<-rownames(rdata_clr) #Make vector with sample names
pc_out_16sY1<-as.data.frame(rdata_PCX$x[,1:2]) #Get PC1 and PC2
pc_out_meta_16sY1<-as.data.frame(bind_cols(pc_out_16sY1,meta$Treatment,meta$Days)) #Add metadata information
row.names(pc_out_meta_16sY1)<-row_16sY1

#transform into multivariate

mvar_clr<-mvar(rdata_clr)

PCA <- ggplot(pc_out_16sY1, aes(x=PC1,y=PC2, color=meta$Days, shape=meta$Treatment))+
  geom_point(size=3)+
  geom_text(aes(label=row.names(pc_out_meta_16sY1)),hjust=0, vjust=0) +
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(rdata_PCX$sdev[1]^2/mvar_clr*100, digits=2), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(rdata_PCX$sdev[2]^2/mvar_clr*100, digits=2), "%", sep="")) +
  ggtitle("Food Post Harvest Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='viridis')

#aldex2 for each treatment in facility and day


#F1 by day
groups <- with(meta,ifelse(as.factor(Days)%in% c("3D"), c("3D"), c("OtherGrowth")))
vdr <- aldex.clr(rdata, groups, mc.samples=128, verbose=FALSE)
vdr_t<-aldex.ttest(vdr, groups, paired.test = FALSE)
vdr_effect <- aldex.effect(vdr, groups, verbose=FALSE)
vdr_all <- data.frame(vdr_t, vdr_effect)

sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)
F13Dtable <-xtable(vdr_all[sig_by_both, c(1:11)], caption="Table of significant taxa", digits=3,label="sig.table", align=c("l",rep("r",11) ))
print.xtable(F13Dtable, type="html", file="F13D_Table.html")

groups <- with(meta,ifelse(as.factor(Days)%in% c("5D"), c("5D"), c("OtherGrowth")))
vdr <- aldex.clr(rdata, groups, mc.samples=128, verbose=FALSE)
vdr_t<-aldex.ttest(vdr, groups, paired.test = FALSE)
vdr_effect <- aldex.effect(vdr, groups, verbose=FALSE)
vdr_all <- data.frame(vdr_t, vdr_effect)

sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)
F15Dtable <-xtable(vdr_all[sig_by_both, c(1:11)], caption="Table of significant taxa", digits=3,label="sig.table", align=c("l",rep("r",11) ))
print.xtable(F15Dtable, type="html", file="F15D_Table.html")

groups <- with(meta,ifelse(as.factor(Days)%in% c("15D"), c("15D"), c("OtherGrowth")))
vdr <- aldex.clr(rdata, groups, mc.samples=128, verbose=FALSE)
vdr_t<-aldex.ttest(vdr, groups, paired.test = FALSE)
vdr_effect <- aldex.effect(vdr, groups, verbose=FALSE)
vdr_all <- data.frame(vdr_t, vdr_effect)

sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)
F15Dtable <-xtable(vdr_all[sig_by_both, c(1:11)], caption="Table of significant taxa", digits=3,label="sig.table", align=c("l",rep("r",11) ))
print.xtable(F15Dtable, type="html", file="F15D_Table.html")

m<-read.csv("D3v5Dlda.csv", check.names = FALSE)
ggplot(m, aes(x = m$Family, y = m$effect)) +
       geom_bar(stat = "identity") +  # color by class
       coord_flip() +
        expand_limits(y=c(-3,3)) +
  ggtitle("3 day vs 5 day ALDEx2 Results")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

n<-read.csv("D5v15Dsu.csv", check.names = FALSE)

ggplot(n, aes(x = n$OTU, y = n$Effect)) +
       geom_bar(stat = "identity") +  # color by class
     coord_flip() +
      expand_limits(y=c(-3,3)) +
      ggtitle("5 day vs 15 day ALDEx2 Results")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

y<-read.csv("D3v15Dsum.csv", check.names = FALSE)

ggplot(y, aes(x = y$OTU, y = y$Effect)) +
  geom_bar(stat = "identity") +  # color by class
  coord_flip() +
  expand_limits(y=c(-3,3)) +
  ggtitle("5 day vs 15 day ALDEx2 Results")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))


par(mfrow=c(1,3)) +
  x +
  y +
  z


grid.arrange(
  tableGrob(otu),
  x,
  y,
  ncol = 2,
  widths = c(1.5, 1),
  clip = FALSE
)

#identify outliers and remove in excel

otus <-read.csv("F1otusrefine.csv", header=TRUE, row.names=1, check.names = FALSE)
meta <-read.csv("biometaF1refine.csv", header=TRUE, row.names=1, check.names = FALSE)

rdata<-otus
rdata_t<-t(rdata)
rdata_r<-t(cmultRepl((rdata_t), method="CZM", output= "p-counts"))
rdata_p<-apply(rdata_r, 2, function(x){x/sum(x)})
rdata_pf<- rdata_r[apply(rdata_p, 1, min) > 0.000001,]
names_add<- rownames(rdata_pf)[ + order(apply(rdata_pf, 1, sum), decreasing=T)]
rdata_pr<-rdata_pf[names_add,]
rdata_clr<-t(apply(rdata_pr, 2, function(x){log(x)-mean(log(x))}))
rdata_PCX<- prcomp(rdata_clr)

rdata_dist<-dist(rdata_clr, method='euclidean')

permanovaF1days<-pairwise.adonis(rdata_dist, factors=meta$Days, perm = 999, p.adjust.m = 'bonferroni')

permanovaF1treatment<-pairwise.adonis(rdata_dist, factors=meta$Treatment, perm = 999, p.adjust.m = 'bonferroni')

row_16sY1<-rownames(rdata_clr) #Make vector with sample names
pc_out_16sY1<-as.data.frame(rdata_PCX$x[,1:2]) #Get PC1 and PC2
pc_out_meta_16sY1<-as.data.frame(bind_cols(pc_out_16sY1,meta$Treatment,meta$Days)) #Add metadata information
row.names(pc_out_meta_16sY1)<-row_16sY1

#transform into multivariate

mvar_clr<-mvar(rdata_clr)
Days<-(meta$Days)
Treatment<-(meta$Treatment)
dotfill= c("blue","red","green")

F1PCA <- ggplot(pc_out_16sY1, aes(x=PC1,y=PC2, color=Days, shape=Treatment))+
  geom_point(size=3)+
  geom_text(aes(label=row.names(pc_out_meta_16sY1)),hjust=0, vjust=0) +
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(rdata_PCX$sdev[1]^2/mvar_clr*100, digits=2), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(rdata_PCX$sdev[2]^2/mvar_clr*100, digits=2), "%", sep="")) +
  ggtitle("F1 Attached Biomass Microbiota")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_fill_manual(values=dotfill,
                    guide="none")
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='viridis')

#F2 PCA

otus <-read.csv("F2otus.csv", header=TRUE, row.names=1, check.names = FALSE)
meta <-read.csv("F2biometa.csv", header=TRUE, row.names=1, check.names = FALSE)

rdata<-otus
rdata_t<-t(rdata)

# Bayesian multiplicative replacement of count zeros, using the count zero multiplicative method and output counts
rdata_r<-t(cmultRepl((rdata_t), method="CZM", output= "p-counts"))
#convert to proportion
rdata_p<-apply(rdata_r, 2, function(x){x/sum(x)})
#filter
rdata_pf<- rdata_r[apply(rdata_p, 1, min) > 0.000001,]
names_add<- rownames(rdata_pf)[ + order(apply(rdata_pf, 1, sum), decreasing=T)]
rdata_pr<-rdata_pf[names_add,]
rdata_clr<-t(apply(rdata_pr, 2, function(x){log(x)-mean(log(x))}))
rdata_PCX<- prcomp(rdata_clr)

rdata_dist<-dist(rdata_clr, method='euclidean')

permanovaF2days<-pairwise.adonis(rdata_dist, factors=meta$Days, perm = 999, p.adjust.m = 'bonferroni')

permanovaF2treatment<-pairwise.adonis(rdata_dist, factors=meta$Treatment, perm = 999, p.adjust.m = 'bonferroni')


row_16sY1<-rownames(rdata_clr) #Make vector with sample names
pc_out_16sY1<-as.data.frame(rdata_PCX$x[,1:2]) #Get PC1 and PC2
pc_out_meta_16sY1<-as.data.frame(bind_cols(pc_out_16sY1,meta$Treatment,meta$Days)) #Add metadata information
row.names(pc_out_meta_16sY1)<-row_16sY1

#transform into multivariate

mvar_clr<-mvar(rdata_clr)
Days<-(meta$Days)
Treatment<-(meta$Treatment)
dotfill= c("blue","red","green")

F2PCA <- ggplot(pc_out_16sY1, aes(x=PC1,y=PC2, color=Days, shape=Treatment))+
  geom_point(size=3)+
  geom_text(aes(label=row.names(pc_out_meta_16sY1)),hjust=0, vjust=0) +
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(rdata_PCX$sdev[1]^2/mvar_clr*100, digits=2), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(rdata_PCX$sdev[2]^2/mvar_clr*100, digits=2), "%", sep="")) +
  ggtitle("F2 Attached Biomass Microbiota")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_fill_manual(values=dotfill,
                    guide="none")
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='viridis')

#F3 PCA

otus <-read.csv("F3otus.csv", header=TRUE, row.names=1, check.names = FALSE)
meta <-read.csv("F3biometa.csv", header=TRUE, row.names=1, check.names = FALSE)

rdata<-otus
rdata_t<-t(rdata)
rdata_r<-t(cmultRepl((rdata_t), method="CZM", output= "p-counts"))
rdata_p<-apply(rdata_r, 2, function(x){x/sum(x)})
rdata_pf<- rdata_r[apply(rdata_p, 1, min) > 0.00001,]
names_add<- rownames(rdata_pf)[ + order(apply(rdata_pf, 1, sum), decreasing=T)]
rdata_pr<-rdata_pf[names_add,]
rdata_clr<-t(apply(rdata_pr, 2, function(x){log(x)-mean(log(x))}))
rdata_PCX<- prcomp(rdata_clr)

rdata_dist<-dist(rdata_clr, method='euclidean')

permanovaF3days<-pairwise.adonis(rdata_dist, factors=meta$Days, perm = 999, p.adjust.m = 'bonferroni')

permanovaF3treatment<-pairwise.adonis(rdata_dist, factors=meta$Treatment, perm = 999, p.adjust.m = 'bonferroni')


row_16sY1<-rownames(rdata_clr) #Make vector with sample names
pc_out_16sY1<-as.data.frame(rdata_PCX$x[,1:2]) #Get PC1 and PC2
pc_out_meta_16sY1<-as.data.frame(bind_cols(pc_out_16sY1,meta$Treatment,meta$Days)) #Add metadata information
row.names(pc_out_meta_16sY1)<-row_16sY1

#transform into multivariate

mvar_clr<-mvar(rdata_clr)

Days<-(meta$Days)
Treatment<-(meta$Treatment)
dotfill= c("blue","red","green")

F3PCA <- ggplot(pc_out_16sY1, aes(x=PC1,y=PC2, color=Days, shape=Treatment))+
  geom_point(size=3)+
  geom_text(aes(label=row.names(pc_out_meta_16sY1)),hjust=0, vjust=0) +
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(rdata_PCX$sdev[1]^2/mvar_clr*100, digits=2), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(rdata_PCX$sdev[2]^2/mvar_clr*100, digits=2), "%", sep="")) +
  ggtitle("F3 Attached Biomass Microbiota")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_fill_manual(values=dotfill,
                    guide="none")
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='viridis')

#Combine PCAs- cant ge to work- combine n adobe illustrator

AttachedBiomass= plot_grid(F1PCA+ theme(legend.position = "none"),
                        F2PCA+ theme(legend.position = "none"),
                        F3PCA+ theme(legend.position = "none"),
                        ncol=1, nrow=3, labels = c("A","B","C"), label_size = 20)



#make counts compositional

compbio<-as.data.frame(acomp((rdata_r)), total=1)

otu=otu_table(compbio, taxa_are_rows=TRUE)
tax=tax_table(TAX)
samples=sample_data(MET)

carbom <- phyloseq(otu, tax, samples)

#permanova for 3D

otus <-read.csv("otus.csv", header=TRUE, row.names=1, check.names = FALSE)
rdata<-otus
otu=otu_table(rdata, taxa_are_rows=TRUE)
carbom <- phyloseq(otu, tax, samples)
stomach <- subset_samples(carbom, Days=="3D")
t<-as.data.frame(otu_table(stomach))
rdata_t<-t(t)
rdata_r<-t(cmultRepl((rdata_t), method="CZM", output= "p-counts"))
rdata_p<-apply(rdata_r, 2, function(x){x/sum(x)})
rdata_pf<- rdata_r[apply(rdata_p, 1, min) > 0.000001,]
names_add<- rownames(rdata_pf)[ + order(apply(rdata_pf, 1, sum), decreasing=T)]
rdata_pr<-rdata_pf[names_add,]
rdata_clr<-t(apply(rdata_pr, 2, function(x){log(x)-mean(log(x))}))
rdata_PCX<- prcomp(rdata_clr)

metttt<-as.data.frame(sample_data(stomach))

rdata_dist<-dist(rdata_clr, method='euclidean')

permanova3D<-pairwise.adonis(rdata_dist, factors=metttt$Facility, perm = 999, p.adjust.m = 'bonferroni')



f <- carbom %>%
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

#https://github.com/joey711/phyloseq/issues/901 abundance graphs

compbio<-as.data.frame(acomp((rdata_r)), total=1)

otu=otu_table(compbio, taxa_are_rows=TRUE)
tax=tax_table(TAX)
samples=sample_data(MET)

carbom <- phyloseq(otu, tax, samples)

#transform into relative abundance and filter to a mean threshold

physeq2 = filter_taxa(carbom, function(x) mean(x) > 0.1, TRUE)
physeq3 = transform_sample_counts(carbom, function(x) x / sum(x) )

#subset for incubation period and Facility
stomach <- subset_samples(physeq3, Facility=="F3")

#create dataframe w/o condensed phyla

stom<-psmelt(stomach)

#Turn all OTU into Family counts

glom <- tax_glom(stomach, taxrank = 'Family')
glom # should list # taxa as # phyla
data <- psmelt(glom) # create dataframe from phyloseq object
data$Family <- as.character(data$Family) #convert to character

#simple way to rename phyla with < 1% abundance
data$Family[data$Abundance < 0.01] <- "< 1% abund."

#
mach <- subset_samples(physeq3, Facility=="F1")
#create dataframe w/o condensed phyla

mache<-psmelt(mach)

#Turn all OTU into Family counts

glomm <- tax_glom(mach, taxrank = 'Family')
glomm # should list # taxa as # phyla
dataa <- psmelt(glomm) # create dataframe from phyloseq object
dataa$Family <- as.character(dataa$Family) #convert to character

#simple way to rename phyla with < 1% abundance
dataa$Family[dataa$Abundance < 0.01] <- "< 1% abund."

#
ach <- subset_samples(physeq3, Facility=="F2")
#create dataframe w/o condensed phyla

ach<-psmelt(ach)

#Turn all OTU into Family counts

glommm <- tax_glom(ach, taxrank = 'Family')
glommm # should list # taxa as # phyla
dataaa <- psmelt(glommm) # create dataframe from phyloseq object
dataaa$Family <- as.character(dataaa$Family) #convert to character

#simple way to rename phyla with < 1% abundance
dataaa$Family[dataa$Abundance < 0.01] <- "< 1% abund."


#bubble plot for days overall

bubbleplot_biofilm<-ggplot(data, aes(x=Category, y=Family, size=Abundance, color=Category))+
  geom_point() +
  scale_size(range=c(.1,6), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'),
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='viridis')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_biofilm

#F1 bubble
physeq3 = transform_sample_counts(carbom, function(x) x / sum(x) )

#subset for incubation period and Facility
stomach <- subset_samples(physeq3, Facility=="F1")
#create dataframe w/o condensed phyla

stom<-psmelt(stomach)

#Turn all OTU into Family counts

glom <- tax_glom(stomach, taxrank = 'Family')
glom # should list # taxa as # phyla
data <- psmelt(glom) # create dataframe from phyloseq object
data$Family <- as.character(data$Family)

data$Family[data$Abundance < 0.01] <- "< 1% abund."

data <-read.csv("F3data.csv", header=TRUE, check.names = FALSE)

#bubble plot for days overall

bubbleplot_16sY1<-ggplot(data, aes(x=Days, y=Family, size=Abundance, color=Days))+
  geom_point() +
  scale_size(range=c(.1,6), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("F3 Attached Biomass Relative Abundance at Different Growth Periods")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='viridis')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_16sY1

#plot with condensed phyla into "< 1% abund" category
p <- ggplot(data=data, aes(x=Days, y=Abundance, fill=Family))
p + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                               "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "pink", "red")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

#Permanova from Laura Rolon
#Directly from Phyloseq
F1_16sY1_count <- subset_samples(phyloseq16sY1_family, Facility == "F1")

#From OTU table, example based on Aldex output
#Subset table with relative abundance per sample, only for significant families
Aldex_family_16sY1_F1F3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_family_16sY1_F1F3.sig.row))

Here is the code for two-way permanova:
#Permanova by Lm, Facility and interaction. Use Aitchinson distance matrix

#Permanova
rdata_dist<-dist(rdata_clr, method='euclidean')

permanovaoverall<-pairwise.adonis2(rdata_dist, factors=meta_table$Facility, perm = 999, p.adjust.m = 'bonferroni')

#insignificant differences

#Correct p-values for multiple comparison using bonferroni's method
p.adjust(permanova_16sY1_Lm$`+_vs_-`$`Pr(>F)`, method = 'bonferroni')




