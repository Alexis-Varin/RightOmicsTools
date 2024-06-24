#data (all non-coding RNA symbols in Homo sapiens) were downloaded from https://www.genenames.org/data/genegroup/#!/group/475 as a txt file
ncRNA_human = data.frame(data.table::fread("genegroup_475.txt"))
ncRNA_human = c(ncRNA_human$Approved.symbol, ncRNA_human$Alias.symbols)
ncRNA_human = unique(ncRNA_human)
ncRNA_human = strsplit(ncRNA_human, ",")
ncRNA_human = unlist(ncRNA_human)
usethis::use_data(ncRNA_human)

#data (all pseudogenes from all tissues in Mus musculus) were downloaded from https://rna.sysu.edu.cn/dreamBase2/scrna.php?SClade=mammal&SOrganism=mm10&SDataId=0&SProteinID=0 as a txt file
pseudogenes_mouse = data.frame(data.table::fread("'dreamBase-Expressions_of_psedogenes_2024_05_10_am_28_51.txt"), row.names = 1)
pseudogenes_mouse = rownames(pseudogenes_mouse)
usethis::use_data(pseudogenes_mouse)
