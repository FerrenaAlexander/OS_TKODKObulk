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

ggplot(df, aes(sample, xist, color = sex))+
  geom_point()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab('Normalized Xist expression')


rm(gemnorm, df, xist)




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

saveRDS(file = 'data/coding,rds', coding)

saveRDS(file = 'data/gem.rds', gem)

rm(coding, gem)
