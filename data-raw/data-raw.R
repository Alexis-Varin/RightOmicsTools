# Raw data (all non-coding RNA symbols in Homo sapiens) were downloaded from https://www.genenames.org/data/genegroup/#!/group/475 as a txt file
ncRNA_human = data.frame(data.table::fread("genegroup_475.txt"))
ncRNA_human = c(ncRNA_human$Approved.symbol, ncRNA_human$Alias.symbols)
ncRNA_human = unique(ncRNA_human)
ncRNA_human = strsplit(ncRNA_human, ",")
ncRNA_human = unlist(ncRNA_human)
usethis::use_data(ncRNA_human)

# Raw data (all pseudogenes from all tissues in Mus musculus) were downloaded from https://rna.sysu.edu.cn/dreamBase2/scrna.php?SClade=mammal&SOrganism=mm10&SDataId=0&SProteinID=0 as a txt file
pseudogenes_mouse = data.frame(data.table::fread("'dreamBase-Expressions_of_psedogenes_2024_05_10_am_28_51.txt"), row.names = 1)
pseudogenes_mouse = rownames(pseudogenes_mouse)
usethis::use_data(pseudogenes_mouse)

# 62 cells in Seurat's pbmc3k do not have annotations (because these cells were removed in Seurat tutorial due to high mito percent), randomly assigning cell types to these cells so that various functions used as examples do not throw an error with NA
library(Seurat)
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")
pbmc3k = UpdateSeuratObject(pbmc3k)
pbmc3k.anno = pbmc3k$seurat_annotations
pbmc3k.anno[is.na(pbmc3k.anno)] = sample(levels(pbmc3k.anno), 1)
usethis::use_data(pbmc3k.anno, internal = TRUE)
