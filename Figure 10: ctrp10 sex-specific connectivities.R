#########################################
#                                       #
#             For Figure 10             #
#    Sex specific connectivities        #
#                                       #
#                                       #
#########################################

## Set working directory (folder where the file in the provided link is downloaded)
# i.e.: setwd('C:/My Drive/Data/')

# All the plots generated in these scripts will be saved in that directory

setwd('')


## Download required files and copy into working directory

## Link to R environment with GTEx data
## https://drive.google.com/file/d/104d0lBOFHlt_iwUjzeOT4POIZRV7nY_l/view?usp=drive_link

### Link to metadata including sex datatable
## https://drive.google.com/file/d/13Yg1a0uXagMffX0fTimLUhGhSgw-04I1/view?usp=drive_link

## Link to human/mouse orthologs table
## https://drive.google.com/file/d/1-6R13c546yWgHMlbunW5snDlkflZGicp/view?usp=drive_link




## library packages required

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
## Also here: https://drive.google.com/file/d/1RkEXp5QF1R30YiFea1vOxj3yAGmcqtL7/view?usp=drive_link
# copy to working directory

de_results = read.csv('results from limma on KO over WT.csv')
res1 = de_results
res1$tissue =gsub(".*_", "", res1$ID) 
res1$gene_symbol = gsub("\\_.*","",res1$ID)


##########################
#GTEx sex specific integration

working_dataset=GTEx_subfiltered
GTEx_subfiltered = NULL
GTEx_full = NULL
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



##################################################
#                                                #
#     Human orthologs DEG correlation - Males    #
#                                                #
##################################################

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

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_1', 'gene_2', 'bicor')
set2 = melt(cc1$p)
set1$pvalue = set2$value
set1 = set1[order(set1$pvalue),]
set1 = set1[set1$gene_1 != set1$gene_2,]
set2=NULL

### Human orthologs DEG correlation - Males table, with bicorrelation and p values
## Write table in working directory
#write.csv(set1, 'correlated human DEG orthologs Males.csv')

#############################################################################################################
#                                                                                                           #
# Or download here: https://drive.google.com/file/d/16Drfq8_C2WjYILYE0e4BmAuqcBg2rmTX/view?usp=drive_link   #
#                                                                                                           #
#############################################################################################################


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



####################################################
#                                                  #
#     Human orthologs DEG correlation - Females    #
#                                                  #
####################################################

isogenes = working_datasetF[,colnames(working_datasetF) %in% gene_set]
ii = na.omit(isogenes)

cc1 = bicorAndPvalue(ii, ii, use = 'p')
cc3 = cc1$bicor
cc3[is.na(cc3)] = 0
cc3 = na.omit(cc3)
cc4 = ifelse(cc1$p<0.01, '*', '')

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_1', 'gene_2', 'bicor')
set2 = melt(cc1$p)
set1$pvalue = set2$value
set1 = set1[order(set1$pvalue),]
set1 = set1[set1$gene_1 != set1$gene_2,]
set2=NULL

### Human orthologs DEG correlation - Females table, with bicorrelation and p values

## Write table in working directory
#write.csv(set1, 'correlated human DEG orthologs Females.csv')

##########################################################################################################
#                                                                                                        #
# Or download here: https://drive.google.com/file/d/16IYNeCuW69YoX1_OJObu2ub5ShUz1u3Z/view?usp=sharing   #
#                                                                                                        #
##########################################################################################################


colnames(ii)[1:10]

anno = data.frame(row.names(cc3), Group=gsub(".*_","", row.names(cc3)))
row.names(anno) = row.names(cc3)

anno$row.names.cc3.=NULL
pdf(file = 'global cor structure of DEGS - Females.pdf')
breaksList = seq(-1, 1, by = .1)
pheatmap::pheatmap(cc3, annotation_row = anno, display_numbers = cc4,fontsize_row = 1, fontsize_col = 1,  fontsize_number = 5, labels_col = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(length(breaksList)) )
dev.off()



####################################################
#                                                  #
#     Within tissue gene-gene connectivity         #
#                                                  #
####################################################


full_cors_table = reshape2::melt(cc1$bicor)
head(full_cors_table)
colnames(full_cors_table) = c('gene1', 'gene2', 'bicor') 
full_cors_table$pvalue = reshape2::melt(cc1$p)$value
full_cors_table$sex = paste0('Females')

full_cors_table = as.data.frame(rbind(full_cors_table, ff1))
table(full_cors_table$sex)

ff2 = na.omit(full_cors_table)
ff2 = ff2[!ff2$gene1==ff2$gene2,]

### Human orthologs DEG correlation - Full table (males and females), with bicorrelation and p values

## Write table in working directory
#write.csv(ff2, 'correlated human DEG orthologs Full Table.csv')

##########################################################################################################
#                                                                                                        #
# Or download here: https://drive.google.com/file/d/16PIwlsQw33m7AEPNHo9-FtmIws96SX1C/view?usp=sharing   #
#                                                                                                        #
##########################################################################################################


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


#################################################
#                                               #
#     Across tissues gene-gene connectivity     #
#                                               #
#################################################


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
