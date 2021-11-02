library(DESeq2)
library(tidyverse)
library(cowplot)
library(fgsea)
library(readxl)
library("RColorBrewer")
library(ComplexHeatmap)
library(biomaRt)
library(dendsort)
#library(dendextend)

set.seed(2021)

#read in
gem <- readRDS('data/gem.rds')
gencode <- read.csv('data/gencode.vM23.ids_names_types.csv')
md <- readxl::read_excel('data/metadata.xlsx')



#### plot the lib size of all samples, including failed #####
gencode <- gencode[gencode$gene_type == 'protein_coding',]
gem <- gem[match(gencode$gene_id, rownames(gem)),]

pdf <- data.frame(samp = colnames(gem), 
                  numreadsaligned = colSums(gem),
                  condition = md$Condition,
                  batch = md$Batch,
                  color = md$Color,
                  stringsAsFactors = F
)


#order them from hi to low
pdf$samp <- factor(pdf$samp, levels = pdf[order(pdf$numreadsaligned, decreasing = T),"samp"])
pdf$condition <- factor(pdf$condition, levels = unique(pdf[order(pdf$numreadsaligned, decreasing = T),"condition"]))
cols <- unique(pdf[order(pdf$numreadsaligned, decreasing = T),"color"])

data.frame(cols = cols, condition = levels(pdf$condition))

libsize_rawall <- ggplot(pdf, aes(x = samp, y = numreadsaligned, fill = condition))+
  geom_bar(stat = 'identity')+
  scale_fill_brewer(palette = 'Set1',direction = -1)+
  scale_y_continuous(labels = scales::comma)+
  theme_light()+
  #scale_fill_brewer(palette = 'Set2')+
  scale_fill_manual(values = cols)+
  scale_y_continuous(limits = c(0,25000000), labels = scales::comma)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = 'Number of reads aligned to protein-coding genes', 
       subtitle = 'Non-normalized',
       y = 'Number of reads aligned', x = 'Sample')

libsize_rawall

ggsave(libsize_rawall, filename = 'results/allsamples/libsize-proteincoding.jpg', height = 5, width = 5, dpi = 300)





### filtering ###

#remove empty genes
gem <- gem[rowSums(gem) > 0,]

#filter lowly exp (under 10 counts) genes
dim(gem[rowSums(gem) >= 10,])
gem <- gem[rowSums(gem) >= 10,]
dim(gem)




#### plot the lib size #####



#for plotting only, see how size factor norm affect lib size

sizefactors <- DESeq2::estimateSizeFactorsForMatrix(counts = gem)


gemnorm <- as.data.frame(t( t(gem) / sizefactors ))



#### plot the lib size #####
pdf <- data.frame(samp = colnames(gemnorm), 
                  numreadsaligned = colSums(gemnorm),
                  condition = md$Condition,
                  batch = md$Batch,
                  color = md$Color,
                  stringsAsFactors = F
)


#order them from hi to low
pdf$samp <- factor(pdf$samp, levels = pdf[order(pdf$numreadsaligned, decreasing = T),"samp"])
pdf$condition <- factor(pdf$condition, levels = unique(pdf[order(pdf$numreadsaligned, decreasing = T),"condition"]))
cols <- unique(pdf[order(pdf$numreadsaligned, decreasing = T),"color"])

data.frame(cols = cols, condition = levels(pdf$condition))

libsize_norm <- ggplot(pdf, aes(x = samp, y = numreadsaligned, fill = condition))+
  geom_bar(stat = 'identity')+
  scale_fill_brewer(palette = 'Set1',direction = -1)+
  scale_y_continuous(labels = scales::comma)+
  theme_light()+
  #scale_fill_brewer(palette = 'Set2')+
  scale_fill_manual(values = cols)+
  scale_y_continuous(limits = c(0,25000000), labels = scales::comma)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = 'Number of reads aligned to protein-coding genes', 
       subtitle = 'With DEseq2 size factor normalization',
       y = 'Number of reads aligned', x = 'Sample')

libsize_norm

ggsave(libsize_norm, filename = 'results/allsamples/libsize-proteincoding-normalized.jpg', height = 5, width = 5, dpi = 300)




##### PCA #####
#Novogene PCA showed one sample (DJ208F0) as a major outlier...

##use DESEQ2: size factor norm, and rlog ##

#create dds obj
dds <- DESeqDataSetFromMatrix(gem, md,
                              design = ~ Condition)


#run DESeq2 
dds <- DESeq(dds)

#rlog
#setting blind to T --> makes an extremely tiny difference (tested at ntop = 500, 2000)
rlog <- rlog(dds, blind = F)
#rlog <- vst(dds, blind = F)

#run PCA


#default deseq pca
defaultDEseqpca = F

if(defaultDEseqpca == T){
  
  ntop = 2000
  
  pcadata <- plotPCA(rlog, returnData=T, ntop=ntop, intgroup = c('Condition'))
  
  
  pca_defaultdeseq2 <- ggplot(pcadata, aes(PC1, PC2, col = Condition))+
    geom_point(size = 3)+
    ggrepel::geom_text_repel(aes(label = name), box.padding = 0.4)+
    scale_color_brewer(palette = 'Set1', direction = -1)+
    theme_light()
}


## recode pca myself... ##

#get the regularized log transformed data
rlmat <- assay(rlog)

#beep, since the previos will take a while
beepr::beep()


#get the top row variance genes, as performed by DESeq2 plotPCA function
rvs <- rowVars(rlmat)
rlmat <- rlmat[order(rvs, decreasing = T),]

#set the number of variable genes
ntop = 2000

#see how ntop selects genes
plot(sort(rvs, decreasing = T)) ; abline(v = ntop)

#see how ntop selects genes, with log y
plot(sort(log(rvs), decreasing = T)) ; abline(v = ntop)


#subset rlmat by ntop
rlmat <- rlmat[1:ntop,]

#instead of scaling later, subtract by mean
# https://www.biostars.org/p/387863/
rlmat <- rlmat - rowMeans(rlmat)

#transpose the mat
rlmat <- t(rlmat)

#scale the matrix; z-transformation. scales each column.
#in this case, columns are genes; thus, this is genewise scaling across samples
#rlmat <- scale(rlmat)


# perform a PCA on the data in assay(x) for the selected genes
pcaobj <- prcomp( rlmat )

pdf <- as.data.frame(pcaobj$x[,1:2])
pdf$name <- dds$Sample
pdf$batch <- dds$Batch
pdf$age <- dds$`Age_at_sack (week`
pdf$tumor_location <- dds$Tumor_location
pdf$condition <- factor(dds$Condition)
pdf$color <- factor(md$Color, levels = unique(md$Color)[order(unique(pdf$condition))] )



data.frame(levels(pdf$condition), levels(pdf$color))


pca_recoded <- ggplot(pdf, aes(PC1, PC2, col = condition, shape = batch))+
  geom_point(size = 3)+
  ggrepel::geom_text_repel(aes(label = name), box.padding = 0.4)+
  #scale_color_brewer(palette = 'Set2', direction = -1)+
  scale_color_manual(values = levels(pdf$color))+
  theme_light()

pca_recoded
ggsave('results/allsamples/pca.jpg', pca_recoded, height = 5, width = 6, dpi = 300)


#try to cluster with the HVGs
#HC on correlation matrix instead of distance

# Pairwise correlation between samples (columns)
cols.cor <- cor(t(rlmat), use = "pairwise.complete.obs", method = "pearson")
cols.cor <- as.dist(1-cols.cor)

col_dend <- dendsort(as.dendrogram(hclust(cols.cor)))

colors <- rev(colorRampPalette( brewer.pal(9, "RdYlBu"))(255))

colors <- circlize::colorRamp2(seq(min(-1), max(1), length = length(colors)), colors)

cormat <- cor(t(rlmat))



annotdf <- data.frame(sample = colnames(cormat),
                      sampcheck = md[match(colnames(cormat), md$Sample), "Sample"],
                      cond = md[match(colnames(cormat), md$Sample), "Condition"],
                      color = md[match(colnames(cormat), md$Sample), "Color"])

hacol <- list(Condition = annotdf$Color) ; names(hacol[[1]]) <- annotdf$Condition
ha <- HeatmapAnnotation(Condition = annotdf[,3], col = hacol, name = 'Condition')
hal <- rowAnnotation(Condition = annotdf[,3], col = hacol, name = 'Condition')

hm <- Heatmap(cormat, name = "Pearson r", 
              col=colors,
              cluster_rows = col_dend, cluster_columns = col_dend, 
              row_names_side = 'left',
              column_names_side = 'top',
              column_dend_height = unit(5, "cm"),
              row_dend_width =  unit(5, "cm"),
              #width = unit(12, "cm"), height = unit(12, "cm"),
              top_annotation = ha, left_annotation = hal
)


hm


jpeg('results/allsamples/hierarchicalclustering.jpg', height = 8, width = 10, units = 'in', res = 300)

print(hm)

dev.off()



