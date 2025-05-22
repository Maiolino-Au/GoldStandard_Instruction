# this script

## Load needed libraries

text
```R
# LOAD LIBRARIES
library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)
library(presto)
library(cowplot)
```

## Set up

text
```R
# SET UP NAMES
path_to_data <- "/sharedFolder/Data/"
name_of_the_data <- ""

# Make folders for partials and results
name_new_dir_results <- paste(getwd(), "/Results", sep = "")
if (!dir.exists(name_new_dir_results)) {dir.create(name_new_dir_results)}
name_new_dir_partial <- paste(getwd(), "/Partial", sep = "")
if (!dir.exists(name_new_dir_partial)) {dir.create(name_new_dir_partial)}
```
# Load, normalize, find variable features and scale the Data
## Load the data

paste() is the command to unite two strings; sep = "" specifies which caracter to put between the strings, with "" it puts nothing, by default it ads a space. paste0() is a paste() but with sep = "" as default. \
\
Some dataset have double names for the genes, both the Gene cards symbol (ex.  SRCIN1) and the Ensembl code (ex. ENSG00000277363). With gene.column you select wich one to use, in this case i used gene card. 
```R
# Load the data
sc_data <- Read10X(data.dir = paste(path_to_data, name_of_the_data, sep = ""), 
                   gene.column = 1) 
```

## Create Seurat object

text
```R
# Create Seurat object
sc_data <- CreateSeuratObject(counts = sc_data, 
                              min.cells = 3, 
                              min.features = 500, 
                              project = name_of_the_data, 
                              names.delim = "-", 
                              names.field = 2)
```

## Normalize and Find variable features

text
```R
# Normalize the data
sc_data <- NormalizeData(sc_data, 
                         normalization.method = "LogNormalize", 
                         scale.factor = 1e6)

# Find variable features
sc_data <- FindVariableFeatures(sc_data, 
                                selection.method = "mvp", 
                                nfeatures = 2000)
```

## Scale the data

text
```R
# Scale the data
sc_data <- ScaleData(sc_data)
```


## Save
```R
# Save the Scaled data
print(paste("Saving PCA for time point", name_of_the_data, "in", name_new_dir))
save(sc_data, file = paste(name_new_dir_partial, "/Scaled_", name_of_the_data, ".Robj", sep = ""))
```

# PCA and Clusterization
## Principal Component Analisys

text
```R
# PCA
sc_data <- RunPCA(sc_data, 
                  npcs = 50, # Number of "principal components" (dimensions)
                  verbose = FALSE) # Removes output 
```

This plot shows how much of the diversity is covered by a determinate number of dimensions
```R
print(ElbowPlot(object = sc_data, ndims = 50)) 
```

## Cluster the cells

More dimension mean more precision, but it will aslo take longer. It is unnecessary to go above 50, many use 40, some 30. If you are giust doing excercise 10 is enough
```R
# Cluster the cells
sc_data <- FindNeighbors(sc_data, dims = 1:40)
```

The higher the resolution the more cluster will be found
```R
sc_data <- FindClusters(sc_data, resolution = 1)
```

view stuff
```R
print(table(Idents(sc_data)))
```

## Save
```R
# Save the PCA plot
print(paste("Saving PCA for time point", name_of_the_data, "in", name_new_dir))
save(sc_data, file = paste(name_new_dir_partial, "/PCA_res_", name_of_the_data, ".Robj", sep = ""))
```

# UMAP 

text
```R
# Make and plot UMAP
sc_UMAP <- RunUMAP(sc_data, dims = 1:10)
DimPlot(sc_UMAP, reduction = "umap", label = TRUE, pt.size = 1)
```

# Find all Markers
Analyze the clusters to find which genes are differentially expressed 

text
```R
# Find all markers for every cluster compared to all remaining cells
cluster_markers <- FindAllMarkers(sc_data, 
                                  only.pos = TRUE, # Considera solo i marker espressi positivamente
                                  min.pct = 0.25, # Percentuale minima di espressione nelle cellule del cluster
                                  logfc.threshold = 0.25) # Soglia minima di LogFC
```


Save
```R
# Save the markers
print(paste("Saving cluster markers for time point", name_of_the_data, "in", name_new_dir))
save(cluster_markers, file = paste(name_new_dir_partial, "/cluster_markers_", name_of_the_data, ".Robj", sep = ""))
```

## Search for a particular gene

text
```R
de.genes <- function(genes_oi) {
    print(paste("Finding differentially expressed genes for time point:", name_of_the_data))

    # Find differentially expressed genes
    de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
    print(de_genes)

    # Save the DE genes
    print(paste("Saving differentially expressed genes for time point", name_of_the_data, "in", name_new_dir))
    write.csv(de_genes, file = paste(name_new_dir_results, "/de_genes_", name_of_the_data, ".csv", sep = ""))

    return(de_genes)
}
```

text
```R
de.genes("name of the gene")
```

# Reload the data

For convenience I prepared functions to load the saved files
```R
load.sc_data <- function(name_of_the_data) {
    load(paste(name_new_dir_partial, "/Scaled_", name_of_the_data, ".Robj", sep = ""))
    return(sc_data)
}

load.clusters <- function(name_of_the_data) {
    load(paste(name_new_dir_partial, "/PCA_", name_of_the_data, ".Robj", sep = ""))
    return(sc_data)
}

load.markers <- function(name_of_the_data) {
    load(paste(name_new_dir_partial, "/cluster_markers_", name_of_the_data, ".Robj", sep = ""))
    return(cluster_markers)
}
```
For the scaled data
```R
sc_data <- load.sc_data("")
```
For the clusterized data
```R
sc_data <- load.clusters("")
```
For the markers
```R
cluster_markers <- load.markers("")
```
