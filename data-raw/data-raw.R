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

# Data preparation for pbmc3k, to speed up data loading with Right_Data
utils::download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                     destfile = "pbmc3k_filtered_gene_bc_matrices.tar.gz", quiet = TRUE)
utils::untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
pbmc3k.mat = Seurat::Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
colnames(pbmc3k.mat) = gsub("-1$","",colnames(pbmc3k.mat))
pbmc3k = Seurat::CreateSeuratObject(counts = pbmc3k.mat, project = "pbmc3k", min.cells = 3, min.features = 200)
unlink(c("pbmc3k_filtered_gene_bc_matrices.tar.gz","filtered_gene_bc_matrices"), recursive = TRUE)
pbmc3k = Seurat::NormalizeData(pbmc3k)
pbmc3k = Seurat::FindVariableFeatures(pbmc3k)
pbmc3k = Seurat::ScaleData(pbmc3k, features = rownames(pbmc3k))
pbmc3k = Seurat::RunPCA(pbmc3k)
pbmc3k = Seurat::FindNeighbors(pbmc3k, dims = 1:10)
pbmc3k = Seurat::FindClusters(pbmc3k, resolution = 0.5)
pbmc3k = Seurat::RunUMAP(pbmc3k, dims = 1:10)
pbmc3k.markers = RightOmicsTools::Find_Annotation_Markers(pbmc3k)
RightOmicsTools::DotPlot_Heatmap(pbmc3k, features = pbmc3k.markers, cluster.idents = FALSE, cluster.features = FALSE)
new.cluster.ids = c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T",
                     "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelets")
names(new.cluster.ids) = levels(pbmc3k)
pbmc3k = Seurat::RenameIdents(pbmc3k, new.cluster.ids)
pbmc3k.data = list(anno = pbmc3k@active.ident,
                   hvg = Seurat::VariableFeatures(pbmc3k),
                   umap = pbmc3k[["umap"]])
usethis::use_data(pbmc3k.data, internal = TRUE)
