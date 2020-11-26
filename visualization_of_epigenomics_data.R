library(data.table)
library(ggplot2)
library(reshape2)
library(ggdendro)
library(FactoMineR)
library(edgeR)
library(dynamicTreeCut)
library(pheatmap)

# 2.3.2 load annotation
annot.features<-fread("data/cro_annotation.gff3")
setkey(annot.features, gene_id, start, end)
comparisons<-scan("data/comparisons", what="character")
comparisons2<-gsub("Lea", "leaf", gsub("ARoo", "Aroot", gsub("LBas", "Basal", gsub("LDis", "Distal", comparisons))))
names(comparisons2)<-comparisons
samples<-c("Root", "Flow", "MLea", "MLeacat", "YLeacat", "LBas", "LDis", "ARoo", "YLea")
samples2<-c("Root", "Flow", "Mleaf", "Mleafcat", "Yleafcat", "Basal", "Distal", "Aroot", "Yleaf")

annot<-fread("data/cro_ap_asm_v2.annotation.xls")
annot.uniprot<-fread("data/cro_uniprot.tab")
blastx.cro<-fread("data/cro_ap_asm_v2.transcripts.fasta.blastx.out")
blastx.cro<-blastx.cro[!duplicated(V1)]
colnames(annot.uniprot)[3]<-"V2"
annot.uniprot<-merge(blastx.cro, annot.uniprot)
annot.uniprot$V1<-gsub("CRO_T", "CRO_", annot.uniprot$V1)
go.bp.universe<-lapply(annot.uniprot$V1, function(x)unlist(strsplit(annot.uniprot[V1==x, `Gene ontology (biological process)`], "; ")))
names(go.bp.universe)<-annot.uniprot$V1
go.bp.universe.table<-table(unlist(go.bp.universe))

# 3.5.1 combine comparisons into a single table.
list.dss.dml.C<-readRDS("data/dml_C.RDS")
list.dss.dml.CHG<-readRDS("data/dml_CHG.RDS")
list.dss.dml.CHH<-readRDS("data/dml_CHH.RDS")
list.dss.dmr.C<-readRDS("data/dmr_C.RDS")
list.dss.dmr.CHG<-readRDS("data/dmr_CHG.RDS")
list.dss.dmr.CHH<-readRDS("data/dmr_CHH.RDS")

list.dss.dml<-rbindlist(list("CG"=rbindlist(list.dss.dml.C[comparisons], idcol = "Comp"),
                             "CHG"=rbindlist(list.dss.dml.CHG[comparisons], idcol = "Comp"),
                             "CHH"=rbindlist(list.dss.dml.CHH[comparisons], idcol = "Comp")), idcol="Context")
colnames(list.dss.dml)[4]<-"start"

# 3.5.2 annotate C sites
list.dss.dml[,end:=start]
list.dss.dml<-foverlaps(list.dss.dml, annot.features[,.(start, end, feature, gene_id)], by.x=c("chr", "start", "end"))
list.dss.dml[is.na(start), feature:="intron"]
list.dss.dml$feature<-factor(list.dss.dml$feature, levels=c("upstream", "five_prime_UTR", "exon", "intron", "three_prime_UTR", "downstream"),
                             labels=c("Up", "5'UTR", "Exon", "Intron", "3'UTR", "Down"))


# 3.5.3 inspect methylation profiles per gene features and context
list.dss.dml[,RelPos:=findInterval((i.start-start)/(end-start), seq(0,1, by=0.1)), by=1:nrow(list.dss.dml)]
draw<-list.dss.dml[, .N, by=.(feature, Context, RelPos)]
ggplot(draw[!is.na(RelPos), ], aes(x=RelPos, y=N, colour=Context)) + geom_line() + facet_grid(~feature) + theme_light() + theme(panel.spacing = unit(0, "lines"), axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position="top") + ylab("Number of DMC")

# 3.5.4 clusters samples according to profiles
#combine numbers of sites for each context and per comparison
list.dss.dml.counts<-rbindlist(list("CG"=data.table(melt(lapply(list.dss.dml.C[comparisons], nrow))),
                                    "CHG"=data.table(melt(lapply(list.dss.dml.CHG[comparisons], nrow))),
                                    "CHH"=data.table(melt(lapply(list.dss.dml.CHH[comparisons], nrow)))), idcol="Context")
list.dss.dml.counts[,L2:=comparisons2[L1]]
list.dss.dml.counts[,c("From", "To"):=as.list(strsplit(L2, "_")[[1]]), by=1:nrow(list.dss.dml.counts)]
list.dss.dml.counts$From<-factor(list.dss.dml.counts$From, levels=rev(samples2))
list.dss.dml.counts$To<-factor(list.dss.dml.counts$To, levels=rev(samples2))

#use ggdendro to visualize a dendrogram
#the distance is based on pairwise differences in DML numbers
list.labels<-list.coord<-vector(mode="list", length=3)
names(list.labels)<-names(list.coord)<-c("CG", "CHG", "CHH")
for (x in c("CG", "CHG", "CHH")){
  dist.tmp<-data.frame(dcast(list.dss.dml.counts[Context==x], From~To, value.var = "value"), row.names=1)
  dist.tmp2<-as.dist(cbind(rbind("YLeaf"=rep(NA, 8), dist.tmp), "Root"=rep(NA, 9)))
  h.tmp<-hclust(dist.tmp2, method="ward.D2")
  h.tmp$labels<-gsub("C_count_", "", h.tmp$labels)
  ddata <- dendro_data(h.tmp, type = "rectangle")
  list.coord[[x]]<-ddata$segments
  list.labels[[x]]<-ddata$labels
}
table.coord<-rbindlist(list.coord, idcol=T)
table.label<-rbindlist(list.labels, idcol=T)
ggplot(table.coord) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand = c(1, 0)) + theme_dendro() + geom_text(data = table.label, aes(x = x, y = y, label = label), size = 2.5, hjust=0) + facet_wrap(".id", scales="free_x") + theme(panel.spacing = unit(0,"line"))

# 3.5.5 visualize the number of DMC and DMR
list.dss.dmr<-rbindlist(list("CG"=rbindlist(list.dss.dmr.C[comparisons], idcol = "Comp"),
                             "CHG"=rbindlist(list.dss.dmr.CHG[comparisons], idcol = "Comp"),
                             "CHH"=rbindlist(list.dss.dmr.CHH[comparisons], idcol = "Comp")), idcol="Context")
list.dss.dmr<-foverlaps(list.dss.dmr, annot.features[,.(start, end, feature, gene_id)], by.x=c("chr", "start", "end"))
list.dss.dmr[is.na(start), feature:="intron"]
list.dss.dmr$feature<-factor(list.dss.dmr$feature, levels=c("upstream", "five_prime_UTR", "exon", "intron", "three_prime_UTR", "downstream"),
                             labels=c("Up", "5'UTR", "Exon", "Intron", "3'UTR", "Down"))

dm.count<-rbindlist(list("DMC"=list.dss.dml[, .N, by=.(Comp, Context, feature)],
                         "DMR"=list.dss.dmr[, .N, by=.(Comp, Context, feature)]), idcol="Level")
dm.count[,Comp2:=comparisons2[Comp]]
dm.count$Comp2<-factor(dm.count$Comp2, levels=dm.count[, sum(N), by=Comp2][order(V1)][,Comp2])
dm.count.dcast<-dcast.data.table(dm.count, Comp2~feature+Context+Level, value.var = "N")
dm.count.dcast[is.na(dm.count.dcast)]=0

#cluster comparisons
tr<-hclust(dist(data.frame(dm.count.dcast, row.names=1)), method="ward.D2")
ddata <- dendro_data(tr, type = "rectangle")
ggplot(ddata$segments) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand = c(1, 0)) + theme_dendro() + geom_text(data = ddata$labels, aes(x = x, y = y, label = label), size = 3, hjust=0)

#plot the number of DMC/DMR per context, feature and comparison
dm.count$Comp2<-factor(dm.count$Comp2, levels=tr$labels[tr$order])
ggplot(dm.count, aes(x=feature, y=Comp2, size=N, colour=Context, shape=Level)) + geom_point(alpha=0.7) + facet_grid(~Context, scales="free") + theme_light() + theme(panel.spacing = unit(0, "lines"), axis.text.x = element_text(angle=270, hjust=0, vjust=0.5), axis.title = element_blank(), axis.text.y=element_text(hjust=0)) + scale_shape_manual(values=c(1,4))

#3.5.6 Principal component analysis
res.pca<-PCA(data.frame(dm.count.dcast, row.names=1), graph=F)
par(mar=rep(0,4))

#plot comparisons (individuals)
plot(res.pca, choix="ind", autoLab="auto", graph.type="ggplot", cex=0.75, select="contrib 15") + theme(plot.title = element_blank(), plot.margin = margin(0, 0, 0, 0, "cm"))
#plot context/features (variables)
plot(res.pca, choix="var", autoLab="auto")

#decreasing contribution of context/features
dim.dt<-rbindlist(list("Dim1"=data.table("Feature"=names(sort(res.pca$var$contrib[,1], decreasing=T)[1:10]), "Value"=sort(res.pca$var$contrib[,1], decreasing=T)[1:10]),
                       "Dim2"=data.table("Feature"=names(sort(res.pca$var$contrib[,2], decreasing=T)[1:10]), "Value"=sort(res.pca$var$contrib[,2], decreasing=T)[1:10])), idcol="Dim")
dim.dt$Feature<-factor(dim.dt$Feature, levels=dim.dt[,max(Value), by=Feature][order(V1),Feature])
ggplot(dim.dt, aes(x=Feature, y=Value)) + geom_bar(stat="identity") + facet_wrap(~Dim, scales="free",strip.position = "right") + coord_flip() + theme(axis.title.y = element_blank()) +ylab("Contributions") + theme_light()

#decreasing contributions of comparisons
sort(res.pca$ind$contrib[,1], decreasing=T)[1:10]
sort(res.pca$ind$contrib[,2], decreasing=T)[1:10]

# 3.5.7 correlation between logFC and methylation level
#load expression matrix
exp<-fread("data/200220_cro.numreads")
exp$Name<-gsub("CRO_T", "CRO_", exp$Name)

#detect DEG in comparisons of interest
list.comp<-t(combn(samples, 2))
list.comp<-sapply(1:nrow(list.comp), function(x)list(list.comp[x,]))
samples.exp<-c("SRR1271858", "TDDB_CAF2_ATCACG_L008", "C4_TGACCA_L003", "ERR1512373", "ERR1512376", "TDDB_CAB19_TAGCTT_L008", "TDDB_CAA15_ACTTGA_L008", "ERR2112587", "ERR1512377")

exp.full<-lapply(1:length(list.comp), function(x){
  print(x)
  sample1<-list.comp[[x]][1]
  sample2<-list.comp[[x]][2]
  sample.exp.A=samples.exp[which(samples==sample1)]
  sample.exp.B=samples.exp[which(samples==sample2)]
  y.exp<-DGEList(counts=exp[,.SD, .SDcols=c(sample.exp.A,sample.exp.B)], group=c(2,1))
  rownames(y.exp$counts)<-exp[,Name]
  y.exp$common.dispersion=0.1
  et.exp<-exactTest(y.exp)
  et.exp.dt<-data.table(et.exp$table, keep.rownames = T)
  et.exp.dt[PValue<0.05,]
})
names(exp.full)<-unlist(lapply(list.comp, function(x)paste(x, collapse="_")))
exp.table<-rbindlist(exp.full, idcol="Comp")
exp.table<-exp.table[,Comp2:=comparisons2[Comp]]

#cluster comparisons according to DEG number SupplFig2
list.exp<-exp.table[, .N, by=Comp]
list.exp[,Comp2:=comparisons2[Comp]]
list.exp[,c("From", "To"):=as.list(strsplit(Comp2, "_")[[1]]), by=1:nrow(list.exp)]
list.exp$From<-factor(list.exp$From, levels=rev(samples2))
list.exp$To<-factor(list.exp$To, levels=rev(samples2))
dist.tmp<-data.frame(dcast(list.exp, From~To, value.var = "N"), row.names=1)
dist.tmp2<-as.dist(cbind(rbind("YLeaf"=rep(NA, 8), dist.tmp), "Root"=rep(NA, 9)))

#cluster samples according to the number of DEGs
h.tmp<-hclust(dist.tmp2, method="ward.D2")
ddata <- dendro_data(h.tmp, type = "rectangle")
ggplot(ddata$segments) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand = c(1, 0)) + theme_dendro() + geom_text(data = ddata$labels, aes(x = x, y = y, label = label), size = 2.5, hjust=0)

#merge DEG and DMR
dmr.exp<-merge(list.dss.dmr, exp.table, by.x=c("chr", "Comp"), by.y=c("rn", "Comp"))
dmr.exp[, length(unique(chr)), by=.(Context, Comp)][, sum(V1), by=Context]$V1/list.dss.dmr[, length(unique(chr)), by=.(Context, Comp)][, sum(V1), by=Context]$V1*100

#calculate correlations between logFC and methylation difference
dmr.exp[, c("pcc", "pval"):=list(cor.test(logFC, diff.Methy)$estimate, cor.test(logFC, diff.Methy)$p.value), by=.(feature, Context)]
dmr.exp[,"sym":=symnum(pval, corr = FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1),symbols = c("***","**","*","."," "))]
dmr.exp$feature<-factor(dmr.exp$feature, levels=rev(levels(dmr.exp$feature)))
dmr.exp[,Comp2:=comparisons2[Comp]]
ggplot(dmr.exp[, .SD[1],by=.(feature, Context)], aes(x=Context, y=feature, fill=pcc)) + geom_tile() + scale_fill_distiller(palette="PuOr", breaks = c(-1, 0, 1), limits=c(-1, 1)) + theme_light() + 
  theme(panel.grid = element_blank(), axis.title = element_blank(), legend.position="top", axis.text.x = element_text(angle=270, hjust=0, vjust=0.5)) + geom_text(aes(label=sym), size=2, nudge_y = 0.1, fontface="bold") + geom_text(data=dmr.exp[,.N, by=.(feature, Context)], aes(label=N, x=Context, y=feature), nudge_y=-0.1, inherit.aes = F, size=2, fontface="bold")

#plot relationships logFC=f(Methylation difference)
ggplot(dmr.exp, aes(x=diff.Methy, y=logFC)) + geom_point(aes(colour=Context), size=1, alpha=0.5) + facet_wrap(~feature, strip.position = "right") + theme_light() + theme(panel.grid = element_blank(), axis.text.x=element_text(angle=270, hjust=0, vjust=0.5)) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +xlab("Methylation difference") + ylab("Expression logFC") + guides(colour = guide_legend(override.aes = list(size=4)))

# 3.5.8 Methylation state of DEG
#table with methylation difference
dmr.exp.table<-dcast.data.table(dmr.exp, chr~Comp2, value.var="diff.Methy", fun.aggregate = mean)
dmr.exp.table[is.na(dmr.exp.table)]<-0

#cluster genes into groups of similar methylation profile
tr<-hclust(as.dist(1-cor(t(data.frame(dmr.exp.table, row.names=1)))), method="ward.D2")
#use dynamic tree cut to improve clustering
cuttr.all<-cutreeDynamic(tr, method="tree", minClusterSize =30)
cuttr.all<-cuttr.all+1
# cuttr.all<-cutreeHybrid(tr, distM=as.matrix(dist(data.frame(dmr.exp.table, row.names=1))), minClusterSize =50)
names(cuttr.all)<-tr$labels
#generate annotation table for pheatmap
cuttr.df.all.diffexp<-data.frame("Group"=cuttr.all)
cuttr.df.all.diffexp$Group<-as.character(cuttr.df.all.diffexp$Group)
vec.col<-list("Group"=c("#a6cee3",
                        "#1f78b4",
                        "#b2df8a",
                        "#33a02c",
                        "#fb9a99",
                        "#e31a1c",
                        "#fdbf6f",
                        "#ff7f00",
                        "#cab2d6",
                        "#6a3d9a",
                        "#ffff99",
                        "#b15928"))
names(vec.col[[1]])<-as.character(1:12)
pheatmap(data.frame(dmr.exp.table, row.names=1), cluster_rows = tr, show_rownames = F, clustering_method = "ward.D2", annotation_row = cuttr.df.all.diffexp, annotation_colors = vec.col, fontsize=8, treeheight_row=10, treeheight_col = 10)

#the following code is used to sum up
#methylation and expression trends
rn<-rownames(data.frame(dmr.exp.table, row.names=1)) 
cn<-colnames(data.frame(dmr.exp.table, row.names=1))
#table with logFC
exp.dmr.table<-dcast.data.table(dmr.exp, chr~Comp2, value.var = "logFC", fun.aggregate = function(x)x[1])
exp.dmr.table<-data.frame(exp.dmr.table, row.names = 1)[rn,cn] 
exp.dmr.table<-apply(exp.dmr.table, 2, function(x)ifelse(x>0, "+", "-"))
exp.dmr.table[is.na(exp.dmr.table)]<-""
exp.dmr.table[data.frame(dmr.exp.table, row.names=1)==0]<-""

m.status<-dmr.exp.table[, lapply(.SD, function(x)ifelse(x==0, "", ifelse(x>0, "M+", "M-"))), .SDcols=2:ncol(dmr.exp.table)]
exp.dmr.table<-data.table(exp.dmr.table)
e.status<-exp.dmr.table[, lapply(.SD, function(x)ifelse(x=="", "", ifelse(x=="+", "E+", "E-"))), .SDcols=1:ncol(exp.dmr.table)]
status.em<-do.call("cbind", lapply(1:ncol(m.status), function(x)paste(unlist(m.status[,..x]), unlist(e.status[,..x]), sep="")))
colnames(status.em)<-colnames(e.status)
status.em<-cbind("Group"=paste("K", cuttr.all, sep=""), "Genes"=dmr.exp.table$chr, status.em)
status.em.m<-melt.data.table(data.table(status.em), id.vars=c("Group", "Genes"))
status.em.m<-status.em.m[, .N, by=.(value, Group)]
status.em.m$Group<-factor(status.em.m$Group, levels=rev(paste("K", unique(cuttr.all[tr$order]), sep="")))
status.em.m[, c("Methyl", "logFC"):=as.list(strsplit(as.vector(value), "E")[[1]]), by=1:nrow(status.em.m)]
status.em.m$Methyl[is.na(status.em.m$Methyl)]<-"nM"
ggplot(status.em.m[Methyl!="nM"], aes(x=N, y=Group, shape=logFC, fill=Methyl, colour=Methyl)) + geom_point(size=2) + theme_light() +theme(legend.position = "top", axis.title.y = element_blank(), axis.text.x = element_text(angle=270, hjust=0, vjust=0), legend.box="vertical") + xlab("events per group") + scale_shape_manual(values=c(25, 24)) + scale_color_manual(values=c("#4575b4", "#d73027")) + scale_fill_manual(values=c("#74add1", "#f46d43"))

#repartition of DMR & DMC and differentially expressed genes
counts.demg<-rbindlist(list("DMC"=merge(merge(exp.table[,.N, by=Comp], dml.exp[,length(unique(chr)), by=Comp], by="Comp"), list.dss.dml[, length(unique(chr)), by=Comp], by="Comp"),
                            "DMR"=merge(merge(exp.table[,.N, by=Comp], dmr.exp[,length(unique(chr)), by=Comp], by="Comp"), list.dss.dmr[, length(unique(chr)), by=Comp], by="Comp")), idcol="Level")
counts.demg.deg<-counts.demg[Level=="DMC"]
counts.demg.deg[,Level:="DEG"]
counts.demg.deg[,c("V1.x", "V1.y"):=0]
counts.demg<-rbind(counts.demg, counts.demg.deg)
counts.demg[Level!="DEG", N:=0]
counts.demg.m<-melt(counts.demg)
counts.demg.m[,Comp2:=comparisons2[Comp]]
counts.demg.m$Comp2<-factor(counts.demg.m$Comp2, levels=counts.demg.m[Level=="DEG" & variable=="N"][order(value), Comp2])
counts.demg.m$variable<-factor(counts.demg.m$variable, labels=c("DEG", "DEM locus", "DM locus"))
counts.demg.text<-counts.demg[Level!="DEG", list("ratio"=V1.x/V1.y*100, "value"=(V1.x+V1.y)), by=.(Comp, Level)]
counts.demg.text[,Comp2:=comparisons2[Comp]]
ggplot(counts.demg.m, aes(x=Comp2, y=value/1000, fill=variable)) + geom_bar(stat="identity") + facet_grid(~Level, scales="free_x") + coord_flip() + theme_light() + scale_fill_brewer(palette="Set2") + theme(axis.title.y = element_blank(), legend.position="top", legend.title = element_blank(), axis.text.x = element_text(angle=315, hjust=0)) +ylab("Number of gene models (x1000)") + geom_text(data=counts.demg.text, aes(label=round(ratio, 0),x=Comp2, y=value/1000), inherit.aes = F, nudge_y=0.1, size=2) 

###1 comp vs all other: summarize biological functions with uniprot keywords
list.all.comp<-lapply(samples, function(x)paste(x, samples, sep="_"))
list.all.comps<-lapply(list.all.comp, function(x){
  list.tmp<-lapply(x, function(y){
    if ( y %in% dmr.exp$Comp) {
      dmr.exp[Comp==y] } else {
        tmp<-strsplit(y, "_")[[1]]
        y.new<-paste(tmp[2], tmp[1], sep="_")
        dmr.tmp<-dmr.exp[Comp==y.new]
        dmr.tmp[, Comp:=y]
        dmr.tmp[,logFC:=-logFC]
        dmr.tmp[,diff.Methy:=-diff.Methy]
      }
  })
  names(list.tmp)<-x
  rbindlist(list.tmp, idcol="Comp3")
})
names(list.all.comps)<-samples

very.large.kw<-table(unlist(strsplit(as.vector(annot.uniprot[,Keywords]), ";")))
keyw.all.prot<-lapply(list.all.comps, function(x){
  set.Mmore.Emore<-x[diff.Methy>0 & logFC >0, chr]
  set.Mmore.Eless<-x[diff.Methy>0 & logFC<0, chr]
  set.Mless.EMore<-x[diff.Methy<0 & logFC>0, chr]
  set.Mless.Eless<-x[diff.Methy<0 & logFC<0, chr]
  list.tmp<-lapply(list(set.Mmore.Emore, set.Mmore.Eless, set.Mless.EMore, set.Mless.Eless), function(y){
  table.tmp<-table(unlist(strsplit(as.vector(annot.uniprot[V1 %in% y, Keywords]), ";")))
  sort(table.tmp[!names(table.tmp) %in% c("3D-structure", "Reference proteome", "Repeat", "Direct protein sequencing")], decreasing=T)[1:20]
  })
  names(list.tmp)<-c("M+E+", "M+E-", "M-E+", "M-E-")
  list.tmp
})
keyw.all.prot.m<-data.table(melt(keyw.all.prot))
to.rm<-as.vector(keyw.all.prot.m[, .N, by=.(Var1, L1)][, .N, by=Var1][N==9, Var1])[-18]
keyw.all.prot.m$L1<-factor(keyw.all.prot.m$L1, levels=c("ARoo", "Root", "Flow", "MLea", "YLea", "MLeacat", "YLeacat", "LBas", "LDis"), labels=c("Aroot", "Root", "Flow", "Mleaf", "Yleaf", "Mleafcat", "Yleafcat", "Basal", "Distal"))
class.sample<-c(rep("Devel", 5), rep(c("Environ", "Onto"), each=2))
names(class.sample)<-levels(keyw.all.prot.m$L1)
keyw.all.prot.m[, "Class":=class.sample[L1]]
keyw.all.prot.m$Var1<-factor(keyw.all.prot.m$Var1, levels=keyw.all.prot.m[, sum(value), by=Var1][order(V1),Var1])
keyw.all.prot.m$L2<-factor(keyw.all.prot.m$L2, levels=c("M-E-", "M-E+", "M+E-", "M+E+"))
ggplot(keyw.all.prot.m[!(Var1 %in% to.rm) & value>2], aes(x=L1, y=Var1, size=value, colour=Class)) + geom_point() + facet_wrap(~L2, ncol=4) + theme_light(base_size = 10) + theme(axis.text.x=element_text(angle=315, hjust=0), axis.title = element_blank(), legend.position = "right") + scale_size_continuous(name="Number of entries") + scale_colour_brewer(palette="Dark2")

ggplot(keyw.all.prot.m[!(Var1 %in% to.rm) & value>2], aes(x=L1, y=Var1, size=value, colour=Class)) + 
  geom_point() + facet_grid(~L2, scales="free") + theme_light(base_size = 10) + 
  theme(axis.text.x=element_text(angle=315, hjust=0), axis.title = element_blank(), axis.text.y = element_text(angle=45, vjust=0), legend.position = "top") + scale_size_continuous(name="Number of entries") + scale_colour_brewer(palette="Dark2") + coord_flip()


