## Download required files and copy into working directory

## Link to R environment with GTEx data
## https://drive.google.com/file/d/16Wh-6nagVIBXwrFINCoFKBB5zuUqGFeL/view?usp=drive_link

### Link to metadata including sex datatable
## https://drive.google.com/file/d/13Yg1a0uXagMffX0fTimLUhGhSgw-04I1/view?usp=drive_link

## Link to human/mouse orthologs table
## https://drive.google.com/file/d/1-6R13c546yWgHMlbunW5snDlkflZGicp/view?usp=drive_link

# set working directory
setwd('')

library(ggVennDiagram)
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(WGCNA)
library(limma)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)

## load GTEx environment
load('GTEx NA included env.RData')

###############################################################
# DEG results, from "Differential expression and integration" section
# copy to working directory
de_results = read.csv('results from limma on KO over WT.csv')
res1 = de_results
res1$tissue =gsub(".*_", "", res1$ID) 
res1$gene_symbol = gsub("\\_.*","",res1$ID)


##########################
#GTEx sex specific integration
working_dataset=GTEx_subfiltered
GTEx_subfiltered = NULL
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))
working_dataset[1:5,1:5]

sex_table = read.delim('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
table(sex_table$sexMF)
new_trts = sex_table[sex_table$GTEx_ID %in% row.names(working_dataset),]
table(new_trts$sexMF)
#  F   M 
#327 653

males = new_trts[new_trts$sexMF=='M',]
females = new_trts[!new_trts$sexMF=='M',]

#Subset two datasets for expression based on sex
working_datasetM = working_dataset[row.names(working_dataset) %in% males$GTEx_ID,]
working_datasetF = working_dataset[row.names(working_dataset) %in% females$GTEx_ID,]
working_dataset = NULL
full_sig_degs = res1[res1$P.Value<0.001,]

## Integration with human orthologs
orths = read.delim('Mouse Gene info with Human Orthologues.txt')
full_sig_degs$human_orth = orths$human_orth[match(full_sig_degs$gene_symbol, orths$Symbol)]
table(full_sig_degs$tissue)
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='Liver', 'Liver', '')
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='iWAT', 'Adipose - Subcutaneous', paste0(full_sig_degs$human_tissue))
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='gWAT', 'Adipose - Visceral (Omentum)', paste0(full_sig_degs$human_tissue))
full_sig_degs$human_tissue = ifelse(full_sig_degs$tissue=='Muscle', 'Muscle - Skeletal', paste0(full_sig_degs$human_tissue))
full_sig_degs = na.omit(full_sig_degs)
full_sig_degs$gene_tissueH = paste0(full_sig_degs$human_orth, '_', full_sig_degs$human_tissue)
full_sig_degs = full_sig_degs[!full_sig_degs$gene_tissueH=='ZBTB9_Adipose ??? Subcutaneous',]
gene_set = c(full_sig_degs$gene_tissueH, 'C1QL2_Adipose - Visceral (Omentum)', 'C1QL2_Adipose - Subcutaneous')

isogenes = working_datasetM[,colnames(working_datasetM) %in% gene_set]
ii = na.omit(isogenes)

cc1 = bicorAndPvalue(ii, ii, use = 'p')
cc3 = cc1$bicor
cc3[is.na(cc3)] = 0
cc4 = ifelse(cc1$p<0.01, '*', '')

colnames(ii)[1:10]

anno = data.frame(row.names(cc3), Group=gsub(".*_","", row.names(cc3)))
row.names(anno) = row.names(cc3)

anno$row.names.cc3.=NULL
pdf(file = 'global cor structure of DEGS - Males.pdf')
breaksList = seq(-1, 1, by = .1)
pheatmap::pheatmap(cc3, annotation_row = anno, display_numbers = cc4,fontsize_row = 1, fontsize_col = 1,  fontsize_number = 5, labels_col = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(length(breaksList)) )
dev.off()

full_cors_table = reshape2::melt(cc1$bicor)
head(full_cors_table)
colnames(full_cors_table) = c('gene1', 'gene2', 'bicor') 
full_cors_table$pvalue = reshape2::melt(cc1$p)$value
full_cors_table$sex = paste0('Males')
ff1 = full_cors_table
#####Now in Females
isogenes = working_datasetF[,colnames(working_datasetF) %in% gene_set]
ii = na.omit(isogenes)

cc1 = bicorAndPvalue(ii, ii, use = 'p')
cc3 = cc1$bicor
cc3[is.na(cc3)] = 0
cc3 = na.omit(cc3)
cc4 = ifelse(cc1$p<0.01, '*', '')

colnames(ii)[1:10]

anno = data.frame(row.names(cc3), Group=gsub(".*_","", row.names(cc3)))
row.names(anno) = row.names(cc3)

anno$row.names.cc3.=NULL
pdf(file = 'global cor structure of DEGS - Females.pdf')
breaksList = seq(-1, 1, by = .1)
pheatmap::pheatmap(cc3, annotation_row = anno, display_numbers = cc4,fontsize_row = 1, fontsize_col = 1,  fontsize_number = 5, labels_col = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(length(breaksList)) )
dev.off()

full_cors_table = reshape2::melt(cc1$bicor)
head(full_cors_table)
colnames(full_cors_table) = c('gene1', 'gene2', 'bicor') 
full_cors_table$pvalue = reshape2::melt(cc1$p)$value
full_cors_table$sex = paste0('Females')

full_cors_table = as.data.frame(rbind(full_cors_table, ff1))
table(full_cors_table$sex)

ff2 = na.omit(full_cors_table)
ff2 = ff2[!ff2$gene1==ff2$gene2,]
ff2$logp = -log(ff2$pvalue)
ff2$absbic = abs(ff2$bicor)
ff2$tissue1 =gsub(".*_", "", ff2$gene1) 
ff2$tissue2 =gsub(".*_", "", ff2$gene2) 
ff2$gene_symbol1 = gsub("\\_.*","",ff2$gene1)
ff2$gene_symbol2 = gsub("\\_.*","",ff2$gene2)

head(ff2)
ff3 = ff2[ff2$tissue1==ff2$tissue2,]
pdf(file = 'male vs female degs bicor comparison within tissue.pdf')
ggboxplot(ff3, x = "sex", y = "absbic",
          color = "sex", palette = "jco",
          add = "jitter") + stat_compare_means() + facet_wrap(~tissue1)
dev.off()

ff3 = ff2[!ff2$tissue1==ff2$tissue2,]
ff3$t1t2 = paste0(ff3$tissue1, ' - ', ff3$tissue2)
pdf(file = 'male vs female degs bicor comparison between tissue.pdf')
ggboxplot(ff3, x = "sex", y = "absbic",
          color = "sex", palette = c('darkslategray4', 'chartreuse3'),
          add = "jitter") + stat_compare_means() + facet_wrap(~t1t2, ncol = 3) + theme_classic() + 
  theme( axis.text = element_text( size = 8 ),
         axis.text.x = element_text( size = 5 ),
         axis.title = element_text( size = 5, face = "bold" ),
         legend.position="none",
         # The new stuff
         strip.text = element_text(size = 4))
dev.off()




ff2$gene1_gene2 = paste0(ff2$gene1, '_', ff2$gene2)
table(ff2$sex)
new_delta = ff2[ff2$sex=='Males',]
new_delta1 = ff2[!ff2$sex=='Males',]
new_delta$female_bicor = new_delta1$bicor[match(new_delta$gene1_gene2, new_delta1$gene1_gene2)]
new_delta$female_pvale = new_delta1$pvalue[match(new_delta$gene1_gene2, new_delta1$gene1_gene2)]
new_delta$delta_bicor = new_delta$bicor - new_delta$female_bicor
new_delta = new_delta[order(abs(new_delta$delta_bicor), decreasing = T),]
write.csv(new_delta, file = 'GTEx DEGS  correlation shift by sex.csv', row.names = F)

nn3 = new_delta[grepl('AGT_Adipose', new_delta$gene1),]
#C1QL2_Adipose - Visceral (Omentum)
#GPAT3_Liver
#SLC41A3_Liver
#ISCA1_Liver

## function to plot gene 1_tissue1 x gene 2_tissue2
gene1 = 'GPAT3_Liver'
gene2 = 'C1QL2_Adipose - Visceral (Omentum)'

plot_gene_cor_bysex = function(gene1, gene2){
  gg1 = as.data.frame(working_datasetM[,colnames(working_datasetM)==gene1])
  row.names(gg1) = row.names(working_datasetM)
  colnames(gg1) = 'gene1'
  gg2 = as.data.frame(working_datasetM[,colnames(working_datasetM)==gene2])
  row.names(gg2) = row.names(working_datasetM)
  gg1$gene2 = gg2$`working_datasetM[, colnames(working_datasetM) == gene2]`[match(row.names(gg1), row.names(gg2))]
  gg1$sex = paste0('Males')
  
  cors_plot = gg1
  
  gg1 = as.data.frame(working_datasetF[,colnames(working_datasetF)==gene1])
  row.names(gg1) = row.names(working_datasetF)
  colnames(gg1) = 'gene1'
  gg2 = as.data.frame(working_datasetF[,colnames(working_datasetF)==gene2])
  row.names(gg2) = row.names(working_datasetF)
  gg1$gene2 = gg2$`working_datasetF[, colnames(working_datasetF) == gene2]`[match(row.names(gg1), row.names(gg2))]
  gg1$sex = paste0('Females')
  
  cors_plot = as.data.frame(rbind(gg1, cors_plot))
  
  pdf(file = paste0(gene1, ' x ', gene2, ' sex-specific correlations.pdf'))
  g1 =ggplot(cors_plot, aes(x=gene1, y=gene2, color = sex, fill = sex)) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_cor(method = "pearson") +
    scale_color_manual(values = c(Males = "red",  Females = 'dodgerblue1'), aesthetics = c("color", "fill")) +
    theme_pubr() +
    labs(x = paste0('Normalized ',gene1,   ' expression'), y = paste0('Normalized ',gene2,   ' expression')) 
  print(g1)
  dev.off()
}

#C1QL2_Adipose - Visceral (Omentum)
#GPAT3_Liver
#SLC41A3_Liver
#ISCA1_Liver

#example
plot_gene_cor_bysex('GPAT3_Liver', 'C1QL2_Adipose - Subcutaneous')
