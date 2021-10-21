library(tidyverse)


md <- readxl::read_excel('data/metadata.xlsx')
files <- list.files('data/counts/', recursive = T, full.names = T)

samplist <- list()
for(file in files){
  
  sampname <- str_split_fixed(file, '/', Inf)[,4]
  basename <- str_sub(sampname, 10)
  
  message(sampname)
  
  
  samp <- read.table(file, sep = '\t', skip = 4)
  samp <- data.frame(row.names = samp[,1],
                     counts = samp[,3])
  
  colnames(samp)[1] <- basename
  
  
  
  samplist[[sampname]] <- samp
  
  rm(samp, sampname, basename)

  
}

gem <- dplyr::bind_cols(samplist)
rm(samplist)


#match with metadata
gem <- gem[,match(md$Sample, colnames(gem))]


geneIDs_noversion <- gsub("\\..*","", rownames(gem))
rownames(gem) <- geneIDs_noversion


### check sex ###

#normalize
sizefactors <- DESeq2::estimateSizeFactorsForMatrix(counts = gem)
gemnorm <- as.data.frame(t( t(gem) / sizefactors ))




#read in md, with sex labels
xist <- t(gemnorm['ENSMUSG00000086503',])

df <- data.frame(xist = xist,
                 sample = rownames(xist),
                 sex = md$Sex)

xistplot <- ggplot(df, aes(sample, xist, color = sex))+
  geom_point()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab('Normalized Xist expression')

xistplot

ggsave(xistplot, filename = 'results/allsamples/xist-expression.jpg', height = 5, width = 5, dpi = 300)


rm(gemnorm, df, xist, xistplot)






#### plot the lib size of all samples, including failed #####
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
  scale_y_continuous(limits = c(0,15000000), labels = scales::comma)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(title = 'Number of aligned reads', 
       subtitle = 'Non-normalized',
       y = 'Number of reads aligned', x = 'Sample')

libsize_rawall

ggsave(libsize_rawall, filename = 'results/allsamples/libsize-whollelib.jpg', height = 5, width = 5, dpi = 300)




#select the protein-coding genes; need to use ensembl
library(biomaRt)
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")




bmlist <- NULL

while(is.null(bmlist) ){
  try(
    bmlist <- getBM(attributes=c("ensembl_gene_id","mgi_symbol","transcript_biotype"),
                    filters = "ensembl_gene_id", values=geneIDs_noversion, mart=mouse)
    
  )
  
}

rm(geneIDs_noversion, mouse)

saveRDS(file = 'data/biomart.rds', bmlist)




#select non-duplicate coding genes...
coding <- bmlist[bmlist$transcript_biotype == 'protein_coding',]
rm(bmlist)

gem <- gem[match(coding$ensembl_gene_id, rownames(gem)),]

saveRDS(file = 'data/coding.rds', coding)

saveRDS(file = 'data/gem.rds', gem)

rm(coding, gem)




