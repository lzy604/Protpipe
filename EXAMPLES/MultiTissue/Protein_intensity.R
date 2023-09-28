#!/usr/bin/env Rscript
# R/4
#proteomics analysis for DIA-NN and Spectronaut quantity estimates

#### ARG PARSING ###################################################################################
package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 'magick', 'ggdendro', 'ecodist',
                 'ggbeeswarm', 'ggrepel', 'ggthemes', 'foreach','reshape2','org.Hs.eg.db',
                 'clusterProfiler','STRINGdb','eulerr','pheatmap')
all((lapply(package_list, require, character.only=TRUE)))

#Working 
setwd('/Users/liz36/Documents/pro_pip/cw_protpip/ProtPipe_de/MultiTissue/')
source('/Users/liz36/Documents/GitHub/DIA_MS_pro_pip/src/functions.R')
outdir='AD_Fasta_new/'
##data and design matrix
dat <- fread('Report_AD_add_muscle.csv')
#dat <- standardize_format(dat)

setnames(dat, trim_colnames(dat))

dat.long <- melt_intensity_table(dat)
dat.long <- dat.long[! is.na(Intensity)][Intensity != 0]



#### MAKE DIRS #####################################################################################

QC_dir <- paste0(outdir, '/QC/')
if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
}

cluster_dir <- paste0(outdir, '/clustering/')
if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
}

DI_dir <- paste0(outdir, '/differential_intensity/')
if(! dir.exists(DI_dir)){
    dir.create(DI_dir, recursive = T)
}

EA_dir <- paste0(outdir, '/Enrichiment_analysis/')
if(! dir.exists(EA_dir)){
  dir.create(EA_dir, recursive = T)
}


#### QC ############################################################################################

#intensity distribution
plot_pg_intensities(dat.long, QC_dir, 'intensities.pdf', plot_title='Un-normalized intensities')

# pgcounts represents the distribution of Protein Groups with Intensity > 0
pgcounts <- dat.long[, .N, by=Sample]
ezwrite(pgcounts, QC_dir, 'all_protein_group_nonzero_counts.tsv')
plot_pg_counts(pgcounts, QC_dir, 'all_protein_group_nonzero_counts.pdf')

#correlations
dat.correlations <- get_spearman(dat)
ezwrite(dat.correlations, QC_dir, 'sample_correlation.tsv')
plot_correlation_heatmap(dat.correlations, QC_dir, 'sample_correlation.png')

##normalization
dat_nor=dat
dat_nor[,3:ncol(dat_nor)]=as.data.frame(apply(dat_nor[,3:ncol(dat_nor)],2,function(x){x*10e4/median(x,na.rm=TRUE)}))
dat.long_nor=melt_intensity_table(dat_nor)
dat.long_nor <- dat.long_nor[! is.na(Intensity)][Intensity != 0]


plot_pg_intensities(dat.long_nor, QC_dir, 'intensities_normalized.pdf', plot_title='normalized intensities')


#### CLUSTERING ####################################################################################

# PCA
pca <- get_PCs(dat_nor)
ezwrite(pca$components, cluster_dir, 'PCA.tsv')
ezwrite(pca$summary, cluster_dir, 'PCA_summary.tsv')
pca$components$Condition=gsub("[0-9]+P$","Plasma",pca$components$Condition)
pca$components$Condition=gsub("[0-9]+C$","CSF",pca$components$Condition)

ezwrite(pca$components, cluster_dir, 'PCA.tsv')

plot_PCs(pca, cluster_dir, 'all_gene_PCA.pdf')

pca$components$Condition=gsub("Brain.*","Brain",pca$components$Condition)
pca$components$Condition=gsub("Muscle.*","Muscle",pca$components$Condition)
plot_PCs(pca, cluster_dir, 'all_gene_PCA_1.pdf')

pca <- get_PCs(na.omit(dat_nor))
pca$components$Condition=gsub("[0-9]+P$","Plasma",pca$components$Condition)
pca$components$Condition=gsub("[0-9]+C$","CSF",pca$components$Condition)
ezwrite(pca$components, cluster_dir, 'commen_gene_PCA.tsv')
ezwrite(pca$summary, cluster_dir, 'commen_gene_PCA_summary.tsv')
plot_PCs(pca, cluster_dir, 'commen_gene_PCA.pdf')

############HC#######################
cluster_data <- DT[,-c(1:2)]
cluster_data[is.na(cluster_data)]=0
cluster_data=cluster_data[which(rowSums(cluster_data)>0),]
log2_cluster_data=log2(cluster_data+1)

dist_mat <- dist(t(log2_cluster_data)) #
hc_cluster <- hclust(dist_mat,method = "complete")
ggdendrogram(hc_cluster,rotate=TRUE) + labs(title='Hierarchical clustering')
cat(paste0('   -> ', output_dir, '\n'))
ggsave(g, filename=paste0(output_dir, 'hc_cluster_log2.pdf'))




####umap
umap=get_umap(dat_nor,15)
umap$condition=gsub("[0-9]+P$","Plasma",umap$condition)
umap$condition=gsub("[0-9]+C$","CSF",umap$condition)
umap$condition=gsub("Brain.*","Brain",umap$condition)
umap$condition=gsub("Muscle.*","Muscle",umap$condition)

plot_umap(umap, cluster_dir, 'all_gene_umap.pdf')


####pheatmap########################################################################

sample_anno <- data.frame(condition = colnames(dat_nor)[3:ncol(dat_nor)])
row.names(sample_anno) <- sample_anno$condition
sample_anno$condition=gsub("_[0-9]*$","",sample_anno$condition)
sample_anno$condition=gsub("[0-9][0-9]*P$","Plasma",sample_anno$condition)
sample_anno$condition=gsub("[0-9]*C$","CSF",sample_anno$condition)

sample_anno$condition=gsub("iNeuron_D28.*","iNeuron_D28",sample_anno$condition)
sample_anno$condition=factor(sample_anno$condition, levels = c("Brain_C9FTLD_TDP","Brain_Control","iMicroglia_D45","Astrocyte","iNeuron_D28","Muscle_IBM","Muscle_IBM_CTRL","CSF" ,"Plasma"    ))

annoCol<-list(condition=c(Brain_C9FTLD_TDP ="#8dd3c7", Brain_Control="#ffffb3",
                          iMicroglia_D45="#bebada", Astrocyte="#fb8072",
                          iNeuron_D28='#80b1d3',
                          Muscle_IBM='#b3de69',Muscle_IBM_CTRL='#fccde5',
                          CSF='#d9d9d9',Plasma='#bc80bd'))

write.csv(sample_anno,"sample_anno.csv")
p=pheatmap(log2(na.omit(dat[,3:ncol(dat_nor)])+1),
         scale='none',
         show_rownames = F,
         annotation_col = sample_anno,
         annotation_colors = annoCol,
         clustering_distance_cols = 'correlation')

p=pheatmap(log2(na.omit(dat[,3:ncol(dat_nor)])+1),
           scale='none',
           show_rownames = F,
           cluster_cols = F,
           annotation_col = sample_anno,
           annotation_colors = annoCol,
           clustering_distance_cols = 'correlation')

save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p,paste0(outdir,'commen_gene_pheatmap.pdf'),12,10)

dat_0=dat_nor
dat_0[is.na(dat_0)]=0

p=pheatmap(log2(dat_0[which(rowSums(dat_0[,3:ncol(dat_0)])>0),3:ncol(dat_0)]+1),
         scale='none',
         clustering_distance_cols = "correlation",
         show_rownames = F,
         annotation_col = sample_anno,
         annotation_colors = annoCol,)
save_pheatmap_pdf(p,paste0(outdir,'all_gene_pheatmap.pdf'),12,10)

#correlation
Brain_dat=fread('Brain_CSF_Plasma.csv')
setnames(Brain_dat, trim_colnames(Brain_dat))
dat.correlations <- get_spearman(Brain_dat)
plot_correlation_heatmap(dat.correlations, QC_dir, 'Brain_sample_correlation.pdf')


fread('Report_AD_add_muscle.csv')
#dat <- standardize_format(dat)
gsub("[0-9][0-9]*P$","Plasma",sample_anno$condition)





p=pheatmap(log2(na.omit(dat[,3:ncol(dat_nor)])+1),
           scale='none',
           show_rownames = F,
           annotation_col = sample_anno,
           annotation_colors = annoCol,
           clustering_distance_cols = 'correlation',
           cluster_cols = F)







save_pheatmap_pdf(p,"pheatmap_commen_gene_nor.pdf")


#upset
upset_Data=dat[,3:ncol(dat)]
upset_Data[is.na(upset_Data)]=0
upset_Data[upset_Data!=0]=1
library(UpSetR)
pdf(paste0(outdir,"upset.pdf"),width = 20,height = 15)
upset(upset_Data,nsets = 74,order.by = "freq")
dev.off()


#90% gene

condition=unique(sample_anno$condition)
for (i in condition){
  samplename=rownames(sample_anno)[grep(i,sample_anno$condition)]
  tmp=dat[,grep(i,sample_anno$condition)+2]
  rownames(tmp)=dat$PG.ProteinGroups
  tmp[is.na(tmp)] <- 0
  tmp$missing= apply(tmp, 1, function(x) sum(x>0))
  tmp_gene=rownames(tmp)[tmp$missing>length(samplename)*0.9]
  gene_name=union(gene_name,tmp_gene)
}



pheatmap(log2(na.omit(dat[,3:ncol(dat_nor)])+1),
         scale='none',
         clustering_distance_cols = "correlation",
         show_rownames = F,
         annotation_col = sample_anno,
         annotation_colors = ann_colors)
sample_anno$condition=c(rep('A',40),rep('B',33))
ann_colors = list(condition = c(A="black",B="orange"))
















ggplot(melt_intensity_table(dat), aes(x=log2(Intensity),color=Sample)) + 
  geom_density()+
  theme_classic()+
  theme(axis.text.x = element_text( angle=90))
ggsave(filename = "intenstiy_unnor.pdf",width = 12,height = 6)

ggplot(melt_intensity_table(dat_nor), aes(x=log10(Intensity),color=Sample)) + 
  geom_density()+
  theme_classic()+
  theme(axis.text.x = element_text( angle=90))
ggsave(filename = "intenstiy_nor.pdf",width = 12,height = 6)

ggplot(melt_intensity_table(na.omit(dat_nor)), aes(x=log10(Intensity),color=Sample)) + 
  geom_density()+
  theme_classic()+
  labs(fill = "",x="",y='Log10 Protein Intensity')+
  theme(axis.text.x = element_text( angle=90))
ggsave(filename = "intenstiy_nor_commen.pdf",width = 12,height = 6)



###common gene normalization
data_common_gene=na.omit(dat)

dat_comm_nor=data_common_gene
dat_comm_nor[,3:ncol(dat_comm_nor)]=as.data.frame(apply(dat_comm_nor[,3:ncol(dat_comm_nor)],2,function(x){x*10e4/mean(x)}))

ggplot(melt_intensity_table(dat_comm_nor), aes(x=log2(Intensity),color=Sample)) + 
  geom_density()+
  theme_classic()+
  theme(axis.text.x = element_text( angle=90))
ggsave("AD_Fasta/QC/intenstiy_commen_nor.pdf",width = 12,height = 6)

pca <- get_PCs(na.omit(dat_comm_nor))
pca$components$Condition=gsub("LP.*","Brain",pca$components$Condition)
pca$components$Condition=gsub("[0-9]*P$","Plasma",pca$components$Condition)
pca$components$Condition=gsub("[0-9]*C$","CSF",pca$components$Condition)
plot_PCs(pca, cluster_dir, 'commen_gene_nor_PCA.pdf')


