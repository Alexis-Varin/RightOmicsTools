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
pbmc3k$orig.ident = factor(sample(c("Donor_1","Donor_2","Donor_3"), ncol(pbmc3k), replace = TRUE))
pbmc3k@meta.data$treatment = factor(sample(c("Control","Stimulated"), ncol(pbmc3k), replace = TRUE))
pbmc3k = Seurat::NormalizeData(pbmc3k)
pbmc3k = Seurat::FindVariableFeatures(pbmc3k)
pbmc3k = Seurat::ScaleData(pbmc3k, features = rownames(pbmc3k))
pbmc3k = Seurat::RunPCA(pbmc3k)
pbmc3k = Seurat::FindNeighbors(pbmc3k, dims = 1:10)
pbmc3k = Seurat::FindClusters(pbmc3k, resolution = 0.5)
pbmc3k = Seurat::RunUMAP(pbmc3k, dims = 1:10)
pbmc3k.markers = RightOmicsTools::Find_Annotation_Markers(pbmc3k)
Seurat::DimPlot(pbmc3k, label = T, pt.size = 1) + Seurat::NoLegend()
Seurat::DimPlot(pbmc3k, label = T, pt.size = 1, group.by = "orig.ident") + Seurat::NoLegend()
RightOmicsTools::DotPlot_Heatmap(pbmc3k, features = pbmc3k.markers, cluster.idents = FALSE, cluster.features = FALSE)
new.cluster.ids = c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T",
                     "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelets")
names(new.cluster.ids) = levels(pbmc3k)
pbmc3k = Seurat::RenameIdents(pbmc3k, new.cluster.ids)
Seurat::DimPlot(pbmc3k, label = T, pt.size = 1) + Seurat::NoLegend()
pbmc3k.data = list(orig.ident = pbmc3k$orig.ident,
                   seurat_annotations = pbmc3k@active.ident,
                   treatment = pbmc3k$treatment,
                   hvg = Seurat::VariableFeatures(pbmc3k),
                   umap = pbmc3k[["umap"]])

# Data preparation for monolps, to speed up data loading with Right_Data
utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7077nnn/GSM7077865/suppl/GSM7077865%5FD1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.tar.gz",
                     destfile = "pbmc_resting.tar.gz", quiet = TRUE)
utils::untar("pbmc_resting.tar.gz", exdir = "pbmc_resting")
pbmc.resting.mat = Seurat::Read10X(data.dir = "pbmc_resting")
pbmc.resting = Seurat::CreateSeuratObject(counts = pbmc.resting.mat, project = "resting", min.cells = 3, min.features = 200)
unlink(c("pbmc_resting.tar.gz","pbmc_resting"), recursive = TRUE)
utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7077nnn/GSM7077866/suppl/GSM7077866%5FG1%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.tar.gz",
                     destfile = "pbmc_stimulated.tar.gz", quiet = TRUE)
utils::untar("pbmc_stimulated.tar.gz", exdir = "pbmc_stimulated")
pbmc.stimulated.mat = Seurat::Read10X(data.dir = "pbmc_stimulated")
pbmc.stimulated = Seurat::CreateSeuratObject(counts = pbmc.stimulated.mat, project = "stimulated", min.cells = 3, min.features = 200)
unlink(c("pbmc_stimulated.tar.gz","pbmc_stimulated"), recursive = TRUE)
pbmc.all = merge(pbmc.resting, pbmc.stimulated)
pbmc.all@project.name = "Monocytes resting or LPS stimulated"
pbmc.all = Seurat::NormalizeData(pbmc.all)
pbmc.all = Seurat::FindVariableFeatures(pbmc.all)
pbmc.all = Seurat::ScaleData(pbmc.all, features = rownames(pbmc.all))
pbmc.all = Seurat::RunPCA(pbmc.all)
pbmc.all = Seurat::IntegrateLayers(pbmc.all, method = Seurat::CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
pbmc.all[["RNA"]] = SeuratObject::JoinLayers(pbmc.all[["RNA"]])
pbmc.all = Seurat::FindNeighbors(pbmc.all, dims = 1:20, reduction = "integrated.cca")
pbmc.all = Seurat::FindClusters(pbmc.all, resolution = 0.5)
pbmc.all = Seurat::RunUMAP(pbmc.all, dims = 1:30, reduction = "integrated.cca")
pbmc.all.markers = RightOmicsTools::Find_Annotation_Markers(pbmc.all)
Seurat::DimPlot(pbmc.all, label = T, pt.size = 1) + Seurat::NoLegend()
Seurat::DimPlot(pbmc.all, label = T, pt.size = 1, group.by = "orig.ident") + Seurat::NoLegend()
RightOmicsTools::DotPlot_Heatmap(pbmc.all, features = pbmc.all.markers, cluster.idents = FALSE, cluster.features = FALSE)
mono = subset(pbmc.all, idents = c(1, 8))
mono = Seurat::FindVariableFeatures(mono)
mono = Seurat::FindNeighbors(mono, dims = 1:20, reduction = "integrated.cca")
mono = Seurat::FindClusters(mono, resolution = 0.2)
mono = Seurat::RunUMAP(mono, dims = 1:30, reduction = "integrated.cca")
mono.markers = RightOmicsTools::Find_Annotation_Markers(mono)
Seurat::DimPlot(mono, label = T, pt.size = 1) + Seurat::NoLegend()
Seurat::DimPlot(mono, label = T, pt.size = 1, group.by = "orig.ident") + Seurat::NoLegend()
RightOmicsTools::DotPlot_Heatmap(mono, features = mono.markers, cluster.idents = FALSE, cluster.features = FALSE)
new.cluster.ids = c("Intermediate Mono", "Classical Mono", "Activated Mono", "Non-classical Mono")
names(new.cluster.ids) = levels(mono)
mono = Seurat::RenameIdents(mono, new.cluster.ids)
Seurat::DimPlot(mono, label = T, pt.size = 1) + Seurat::NoLegend()
mono.data = list(cells = colnames(mono),
                 seurat_annotations = mono@active.ident,
                 hvg = Seurat::VariableFeatures(mono),
                 umap = mono[["umap"]])

# Saving all objects
all.data = list(pbmc3k = pbmc3k.data, mono = mono.data)
usethis::use_data(all.data, internal = TRUE)
