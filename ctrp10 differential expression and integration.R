setwd('G:/My Drive/lab files/Will Wong/CTRP10/')
library(ggVennDiagram)
library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(WGCNA)
library(limma)
library(ggrepel)

full_melted_cnts = read.csv('full melted kallisto tpm CTRP10.csv')
head(full_melted_cnts)

full_melted_cnts$gen_tissue = paste0(full_melted_cnts$gene_symbol, '_', full_melted_cnts$tissue)
full_melted_cnts$genotype = ifelse(grepl('WT', full_melted_cnts$sample), 'WT', 'KO')

full_melted_cnts$mouse_ID = str_sub(full_melted_cnts$sample,-1,-1)
full_melted_cnts$ms_geno = paste0(full_melted_cnts$mouse_ID, '_', full_melted_cnts$genotype)
full_melted_cnts$logcnts = log2(full_melted_cnts$est_counts + 1)

##########################################################

new_cnts = dcast(full_melted_cnts, ms_geno ~ gen_tissue, value.var = 'logcnts', fun.aggregate = mean, na.rm=T)
row.names(new_cnts) = new_cnts$ms_geno
new_cnts$ms_geno=NULL
new_cnts[1:10,1:10]
rowSums(new_cnts[1,])
cnts_mat = new_cnts[,colSums(new_cnts)>5]
cnts_mat$dm = ifelse(grepl('WT', row.names(cnts_mat)), 'WT', 'KO')
cnts_mat$dm = factor(cnts_mat$dm, levels=c('WT', 'KO'))
table(cnts_mat$dm)
design = model.matrix(~dm, data=cnts_mat)
head(design)
table(cnts_mat$dm)
dim(design)
new_cnts1 = as.data.frame(t(cnts_mat[, !colnames(cnts_mat)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)
row.names(fit)[1:10]

res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                     sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
head(res_table)
test1 = as.data.frame(cnts_mat[,colnames(cnts_mat) == 'Zbtb9_iWAT'])
row.names(test1) = row.names(cnts_mat)
test1$cond = ifelse(grepl('HF', row.names(test1)),  'HF', 'Chow') 

write.csv(res_table, file = 'results from limma on KO over WT.csv', row.names = F)

res1 = res_table
head(res1)
#need to play around to assessing proper thresholds.  
test_genes = as.vector(res1$ID[res1$P.Value<0.0001])
label_key = res1$ID[res1$P.Value<0.0001]
res1$label2 = ifelse(res1$ID %in% label_key, paste0(res1$ID), '')
table(res1$label2)
color_key_table = as.data.frame(res1$ID)
colnames(color_key_table) = 'ID'
color_key_table$tissue =gsub(".*_", "", color_key_table$ID)
color_sets = as.data.frame(table(color_key_table$tissue))
colnames(color_sets) = c('tissue_ID', 'Freq')
color_sets$colors = c( 'deepskyblue4',  'darkorange',  'deeppink4', 'forestgreen')

color_key_table$color = color_sets$colors[match(color_key_table$tissue, color_sets$tissue_ID)]

res1$label_col1 = color_key_table$color[match(res1$ID, color_key_table$ID)]

res1$label_col2 = ifelse(res1$P.Value<0.01, paste0(res1$label_col1), 'gray74')
#Number of genes which will be labelled
#Volcano plot
pdf(file = 'Volcano Plot of KO over WT.pdf')
ggplot(res1, aes(x=logFC, y=-log10(P.Value))) + theme_classic() +
  geom_point(aes(x=logFC, y=-log10(P.Value)), color=res1$label_col2)  + geom_label_repel(aes(x=logFC, y=-log10(P.Value), label = res1$label2), color = res1$label_col2, size = 2, label.size=NA, box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf, segment.color = 'grey50')  +   ggtitle('Volcano plot over KO over WT')
dev.off()

new_pal = color_sets$colors
names(new_pal) = color_sets$tissue_ID
pdf(file = paste0('TISSUE LEGEND for Volcano Plot of KO over WT.pdf'))
plot.new()
legend("center",
       legend = names(new_pal),
       fill = new_pal,       # Color of the squares
       border = "black")
dev.off()


proportions_plot = function(pval_threshold){
  cof1_table = res1[res1$P.Value<pval_threshold,]
  cof1_table$FC_cat = ifelse(cof1_table$logFC>1.2, 'up-regulated', 'down-regulated')
  cof1_table$tissue =gsub(".*_", "", cof1_table$ID)
  cof1_table$class_FC = paste0(cof1_table$tissue, '_', cof1_table$FC_cat)
  
  pp1 = cof1_table %>%
    group_by(class_FC) %>%
    summarise(cnt = n()) %>%
    mutate(freq = round(cnt / sum(cnt), 3)) %>% 
    arrange(desc(freq))
  pp1$direction = gsub(".*_", "", pp1$class_FC)
  pp1$tissue = gsub("\\_.*","",pp1$class_FC)
  pdf(file = paste0('Proportions of genotype changes P_less ', pval_threshold, '.pdf'))
  g1 = ggplot(pp1, aes(x=tissue, y=freq, fill=direction)) + geom_bar(stat="identity", width=.5, position = "dodge") + theme_classic() + scale_fill_hue(l=20, c=60) + ylab('Frequency by datatype') + xlab('')+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ggtitle(paste0('Genotype changes P < ', pval_threshold))
  print(g1)
  dev.off()
  }

proportions_plot(0.001)
proportions_plot(0.01)
proportions_plot(0.05)
plot_intersects = function(pval_threshold){
  cof1_table = res1[res1$P.Value<pval_threshold,]
  cof1_table$tissue =gsub(".*_", "", cof1_table$ID)
  table(cof1_table$tissue)
  gene_list <- list(Liver= cof1_table$ID[cof1_table$tissue=='Liver'], gWAT= cof1_table$ID[cof1_table$tissue=='gWAT'], iWAT= cof1_table$ID[cof1_table$tissue=='iWAT'], Muscle= cof1_table$ID[cof1_table$tissue=='Muscle'])
  pdf(file = paste0('Intersection of Genotype changes for genes P_less ', pval_threshold, '.pdf'))
  g2=ggVennDiagram(gene_list, label_alpha = 0, label = "count") + scale_fill_distiller(palette = "RdYlBu") + ggtitle('Intersection of genotype changes for genes P_less ', pval_threshold)
  print(g2)
  dev.off()
} 
plot_intersects(0.05)
plot_intersects(0.01)
plot_intersects(0.001)
plot_intersects(0.0001)


#pathways
library(enrichR)


setEnrichrSite("Enrichr")
dbs1 <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")
dbs = listEnrichrDbs()
cof1_table = res1
cof1_table$FC_cat = ifelse(cof1_table$logFC>1.2, 'up-regulated', 'down-regulated')
cof1_table$tissue =gsub(".*_", "", cof1_table$ID)
table(cof1_table$tissue)

tissue1 = 'gWAT'
get_paths = function(tissue1){
pp2 = cof1_table[cof1_table$tissue %in% tissue1,]
pp1 = pp2[pp2$FC_cat== 'down-regulated'& pp2$P.Value<0.01,]
#first for positive
head(pp1)
name_D = paste0('_', tissue1)
pp1$gene_symbol = gsub(name_D, '', pp1$ID, fixed = T) 
gg1 = pp1$gene_symbol
enriched <- enrichr(gg1, dbs1)
g1 = plotEnrich(enriched[[1]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0('downregulated DEGs ', names(enriched[1]))) +theme(axis.text=element_text(size=5))
                                                                                                                          

g2 = plotEnrich(enriched[[2]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0('downregulated DEGs ', names(enriched[2])))+theme(axis.text=element_text(size=5))

g3 = plotEnrich(enriched[[3]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0('downregulated DEGs ', names(enriched[3])))+theme(axis.text=element_text(size=5))


pp1 = pp2[pp2$FC_cat== 'up-regulated'& pp2$P.Value<0.01,]
#first for positive
head(pp1)
name_D = paste0('_', tissue1)
pp1$gene_symbol = gsub(name_D, '', pp1$ID, fixed = T) 
gg1 = pp1$gene_symbol
enriched <- enrichr(gg1, dbs1)
g4 = plotEnrich(enriched[[1]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0('upregulated DEGs ', names(enriched[1])))+theme(axis.text=element_text(size=5))

g5 = plotEnrich(enriched[[2]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0('upregulated DEGs ', names(enriched[2])))+theme(axis.text=element_text(size=5))

g6 = plotEnrich(enriched[[3]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0('upregulated DEGs ', names(enriched[3])))+theme(axis.text=element_text(size=5))


pdf(file = paste0(tissue1, ' full_paths.pdf'))
gg3 = gridExtra::grid.arrange(g1, g2, g3, g4, g5, g6)
print(gg3)
dev.off()
}
table(cof1_table$tissue)
get_paths('gWAT')
get_paths('Liver')

get_paths('iWAT')
get_paths('Muscle')
