library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(FerrenaBulkRNAseq)

library(mMCPcounter)

set.seed(2021)


#read in
gem <- readRDS('data/gem.rds')
gencode <- read.csv('data/gencode.vM23.ids_names_types.csv')
md <- readxl::read_excel('data/metadata.xlsx')



#keep only protein coding
gencode <- gencode[gencode$gene_type == 'protein_coding',]
gem <- gem[match(gencode$gene_id, rownames(gem)),]


# for easier time, remove id version info
gencode$id_with_version <- gencode$gene_id
split <- str_split_fixed(string = gencode$gene_id, '\\.', 2)
gencode$gene_id <- split[,1] ; rm(split)

rownames(gem) <- gencode$gene_id





### prep for mMCP-Counter ###

# input is normalized, log-transformed counts


#normalize
sizefactors <- DESeq2::estimateSizeFactorsForMatrix(gem)
gemnorm <- as.data.frame(t( t(gem) / sizefactors ))

#log transform
gemlog <- log2(gemnorm + 1)
rm(gemnorm)


#estimate deconvolved results...
dc <- mMCPcounter.estimate(gemlog, features = "ENSEMBL.ID")





#try to plot...
#transpose, samples = rows
dct <- as.data.frame(t(dc))

#merge with metadata
mdx <- cbind(md, dct)

#plot
celltypes <- colnames(dct)
ctlist <- list()
for(ct in celltypes){
   
  pdf <- data.frame(Sample = mdx$Sample,
                       Condition = mdx$Condition,
                    estimate = mdx[,ct],
                       celltype = ct
  )
  
  pdf[pdf$mMCP_estimate==-Inf,"estimate"] <- 0
  
  p <- t.test(pdf[pdf$Condition=='TKO',"estimate"],
              pdf[pdf$Condition=='DKO',"estimate"])$p.value
  
  p <- signif(p, digits = 2)
  
  pdf$celltype_pvalue <- paste0(ct, '\nP = ', p)
  
  ctlist[[ct]] <- pdf
  

}
ctres <- dplyr::bind_rows(ctlist)

#plot them

#sort celltypes by global mean?
ctag <- aggregate(estimate ~ celltype, ctres, mean)
ctres$celltype <- factor(ctres$celltype, levels =   ctag[order(ctag$estimate, decreasing = T),"celltype"])

ctag <- aggregate(estimate ~ celltype_pvalue, ctres, mean)
ctres$celltype_pvalue <- factor(ctres$celltype_pvalue, levels =   ctag[order(ctag$estimate, decreasing = T),"celltype_pvalue"])



#rename estimate column based on tool
colnames(ctres)[3] <- 'mMCP_estimate'


#boxplots comparing conditions

ggplot(ctres, aes(Condition, mMCP_estimate, fill = Condition))+
  geom_point()+
  geom_boxplot()+
  facet_wrap(~ celltype_pvalue)+
  scale_fill_manual(values = rev(unique(md$Color)))+
  labs(caption = 'Statistics via T test')



#sample barplot

#better color scale...
pal <- c( RColorBrewer::brewer.pal(Inf, "Set1"),
          RColorBrewer::brewer.pal(Inf, "Set2")
)

pal <- sample(pal, size = length(levels(ctres$celltype)), replace = F)


ggplot(ctres, aes(Sample, mMCP_estimate, fill = celltype))+
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))+
  scale_fill_manual(values = pal)










