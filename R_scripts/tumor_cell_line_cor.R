#compare expression profiles of TCGA primary tumors and CCLE cell lines; adapted from https://github.com/Bin-Chen-Lab/HCC_NEN/blob/master/tumor_cell_line/compute_tumor_cell_line_cor_update.R

library(data.table)
library(pheatmap)
library(gplots)
library(ggplot2)
library(pheatmap)
library(viridis)
library(edgeR)


cancer = "PAAD" #select tumor type

comparison_gene_set = "TCGA_5k" # select number of most variable genes (other options: "TCGA_10k", "all_genes")

############
# Functions
############

#compare correlations between tumors and dz related cell lines to non dz related cell lines

is_outlier <- function(cancer, cell_line_one_tumor_anno){
  ccle_tcga_target_celllines = CCLE.meta$CCLE.name[(CCLE.meta$disease == cancer)]
  cor_same = cell_line_one_tumor_anno$cor[cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines]
  cor_others = cell_line_one_tumor_anno$cor[!(cell_line_one_tumor_anno$CCLE.name  %in% ccle_tcga_target_celllines)]
  p = wilcox.test(cor_same, cor_others, alternative = 'greater')  
  return(p$p.value)
}

#function to check if the cell line most correlated with the TCGA sample is from the same tumor type

is_outlier3 <- function(cancer, cell_line_one_tumor_anno){
  cell_line_one_tumor_anno.ordered = cell_line_one_tumor_anno[order(-cell_line_one_tumor_anno$cor),]
  if(as.character(cell_line_one_tumor_anno.ordered$disease[1]) == cancer){
    return("YES")
  } else {
    return("NO")
  }
}

###############
#    MAIN 
###############

#read in UQ normalized + log transformed CCLE and TCGA expression data

CCLE= fread("data/CCLE_normalized_expression.txt", data.table = F)
rownames(CCLE) = CCLE$gene_ids
CCLE = CCLE[,-1]

TCGA = fread(paste("data/", cancer, "_normalized_expression.txt", sep = ""), data.table = F)
rownames(TCGA) = TCGA$gene_ids
TCGA = TCGA[,-1]

#read in CCLE mapping info and TCGA purity info

CCLE.meta = read.table("data/CCLE_meta.txt", sep = "\t", header = T, stringsAsFactors = F)
purity_info = read.table("data/TCGA_Purity_Estimates.txt", header = T)

##adjust for tumor purity if purity information available (solid tumors only)
if (cancer == "LAML" | cancer == "DLBC"){
  #purity estimates not available for blood cancers, skip purity adjustment
  TCGA_adjusted = TCGA
}else{
  #remove genes that are significantly correlated with tumor purity
  TCGA = TCGA[, intersect(purity_info$sample, colnames(TCGA))]
  tumors_purity = purity_info[purity_info$sample %in% colnames(TCGA), ]
  tmp=lapply(rownames(TCGA), function(x)cor.test(as.numeric(TCGA[x,]),tumors_purity$purity, method = "s"))
  tmp.df = data.frame(sapply(tmp,function(x)x$p.value), sapply(tmp,function(x)x$estimate))
  rownames(tmp.df) = rownames(TCGA)
  colnames(tmp.df) = c("pvalue", "rho")
  tmp.df$padj = p.adjust(tmp.df$pvalue)
  tmp.df.sig = tmp.df[tmp.df$padj < 0.01 & tmp.df$rho < -0.4 & complete.cases(tmp.df),]
  purity_genes = rownames(tmp.df.sig)
  TCGA = TCGA[!(rownames(TCGA) %in% purity_genes),]
  
  #adjust for tumor purity with linear regression
  tumors_infiltrate = 1 - tumors_purity$purity
  design = model.matrix(~ tumors_infiltrate)
  fit <- lmFit(TCGA,design)
  beta <- fit$coefficients[,2,drop=FALSE]
  TCGA_adjusted = as.matrix(TCGA) - beta %*% t(tumors_infiltrate)
}

## select most variable genes
if (comparison_gene_set == "TCGA_5k"){
  iqr_gene = apply(TCGA_adjusted, 1, IQR)
  varying_genes = iqr_gene[order(iqr_gene, decreasing = T)]
  varying_genes = names(varying_genes[1:5000])
}else if (comparison_gene_set == "TCGA_10k"){
  iqr_gene = apply(TCGA2, 1, IQR)
  varying_genes = iqr_gene[order(iqr_gene, decreasing = T)]
  varying_genes = names(varying_genes[1:10000])
}else if (comparison_gene_set == "all_genes"){
  varying_genes = rownames(TCGA2)
}else{
  stop("error: need comparison gene set")
}

TCGA_CCLE = cbind(TCGA_adjusted[varying_genes,],CCLE[varying_genes,])

##make cancer specific directory if does not exist already
ifelse(!dir.exists(file.path(paste("results/", cancer, sep = ""))), dir.create(file.path(paste("results/", cancer, sep = ""))), FALSE)

TCGA_CCLE_cor = cor(TCGA_CCLE, method = "s")

cell_lines <- colnames(CCLE)
tumor_cell_all = data.frame()
for (tumor in colnames(TCGA)){
  cell_line_one_tumor_cor = TCGA_CCLE_cor[tumor,cell_lines]
  cell_line_one_tumor_cor = data.frame(sample = names(cell_line_one_tumor_cor), cor=cell_line_one_tumor_cor)
  cell_line_one_tumor_anno = merge(CCLE.meta, cell_line_one_tumor_cor, by.y="sample", by.x="barcode")
  cell_line_one_tumor_anno$barcode = tumor
  cell_line_one_tumor_anno$patient_id = substr(tumor,1,12)
  cell_line_one_tumor_anno$sample_id = substr(tumor,1,16)
  cell_line_one_tumor_anno$outlier = is_outlier(cancer, cell_line_one_tumor_anno)
  cell_line_one_tumor_anno$top_cancer_line_matching = is_outlier3(cancer, cell_line_one_tumor_anno) 
  tumor_cell_all = rbind(tumor_cell_all, cell_line_one_tumor_anno)
}

saveRDS(tumor_cell_all, paste("results/", cancer, "/tumor_cell_all_", cancer, ".rds", sep=""))

########################
# visualize results
########################  

readRDS(tumor_cell_all, paste("results/", cancer, "/tumor_cell_all_", cancer, ".rds", sep=""))

tumor_cell_all$Cell.Line = gsub("_.*", "", tumor_cell_all$CCLE.name)
cancer_specific_cell_lines = unique(as.character(tumor_cell_all$Cell.Line[tumor_cell_all$disease == cancer]))

cell_lines = unique(tumor_cell_all$Cell.Line)
samples = unique(tumor_cell_all$barcode)

#violin plot of tumor vs cell line correlations separated by cell line primary tissue sites
pdf(paste("results/", cancer, "/", cancer, "_by_primary_sites.pdf", sep=""))
cell_line_type_order <- aggregate(cor ~ Site.Primary, tumor_cell_all, median)
cell_line_type_ordered <- cell_line_type_order$Site.Primary[order(cell_line_type_order$cor, decreasing=T)]
tumor_cell_all$Site.Primary = factor(tumor_cell_all$Site.Primary, levels = cell_line_type_ordered)
p <- ggplot(tumor_cell_all, aes(Site.Primary, cor))
print(p + geom_violin()  +  ylab("correlation") + geom_boxplot(width=0.2) +
        xlab("") +
        stat_summary(geom = "crossbar", width=0.2, fatten=2, color="red", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +  theme(panel.background = element_rect(color = 'white'), axis.text.x = element_text(angle = 70, hjust = 1, size=8)) 
)
dev.off()

#violin plot of cell line correlations with primary tumors
tumor_cell_all_subset = tumor_cell_all[tumor_cell_all$disease == cancer,]
pdf(paste("results/", cancer, "/", cancer, "_specific_cell_lines.pdf", sep=""))
cell_line_order <- aggregate(cor ~ Cell.Line, tumor_cell_all_subset, median)
cell_line_ordered <- cell_line_order$Cell.Line[order(cell_line_order$cor, decreasing=T)]
tumor_cell_all_subset$Cell.Line = factor(tumor_cell_all_subset$Cell.Line, levels = cell_line_ordered)
p <- ggplot(tumor_cell_all_subset, aes(Cell.Line, cor))
print(p + geom_violin()  +  ylab("correlation") + geom_boxplot(width=0.2) +
        xlab("") +
        stat_summary(geom = "crossbar", width=0.2, fatten=2, color="red", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +  theme(panel.background = element_rect(color = 'white'), axis.text.x = element_text(angle = 70, hjust = 1, size=8))
)
dev.off()


#heatmap of tumor vs cell line correlations

cell_lines = unique(tumor_cell_all_subset$Cell.Line)
samples = unique(tumor_cell_all_subset$barcode)

cell_tumor_matrix = matrix(NA, nrow=length(cell_lines), ncol=length(samples), dimnames=list(cell_lines, samples))
for(i in 1:nrow(tumor_cell_all_subset)){
  cell_tumor_matrix[tumor_cell_all_subset$Cell.Line[i],tumor_cell_all_subset$barcode[i] ] = tumor_cell_all_subset$cor[i]
}

pdf(paste("results/", cancer, "/", cancer, "_heatmap_TCGA_CCLE.pdf", sep=""))
print(
  
  pheatmap(
    mat               = t(cell_tumor_matrix),
    color             = inferno(50),
    border_color      = NA,
    show_colnames     = T,
    show_rownames     = F,
    fontsize          = 8,
    drop_levels = F,
    main = paste(cancer, " correlation matrix", sep = "")
  )
)

dev.off()






