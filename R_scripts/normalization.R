library(data.table)
library(edgeR)
library(sva)

#perform upper-quartile normalization and log transformation on RNA-seq counts data

TCGA_CCLE = fread("TCGA_CCLE_counts.txt")
colnames(TCGA_CCLE) = gsub("\\.", "-", colnames(TCGA_CCLE))
rownames(TCGA_CCLE) = TCGA_CCLE$gene_id
TCGA_CCLE$gene_id = NULL
my_DGEList = DGEList(counts=TCGA_CCLE)
my_DGEList_norm = calcNormFactors(my_DGEList, method="upperquartile")
TCGA_CCLE = cpm(my_DGEList_norm, log=TRUE, prior.count=1)
saveRDS(TCGA_CCLE, "TCGA_CCLE_normalized_expression.rds")

#correct for sequencing platform (GAII vs Hiseq) differences using ComBat
cancers = c("COAD", "READ", "UCEC", "STAD", "LAML")

for (cancer in cancers){
  batches = readRDS(paste("data/batch_info_", cancer, ".rds", sep = ""))
  samples = intersect(batches$Extract.Name, colnames(TCGA_CCLE))
  TCGA = TCGA_CCLE[, samples]
  batches = batches[match(colnames(TCGA), batches$Extract.Name),]
  if (length(unique(batches$Type)) > 1){
    modcombat = model.matrix(~Type, data=batches)
  }else {
    modcombat = model.matrix(~1, data=batches)
  }
  combat <- ComBat(dat=TCGA, batch=batches$Platform, mod=modcombat, par.prior = T)
  saveRDS(combat, paste("data/", cancer, "normalized_expression.rds", sep = ""))
}
