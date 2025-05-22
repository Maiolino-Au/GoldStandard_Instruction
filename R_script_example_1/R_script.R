# LOAD LIBRARIES
library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)
library(dplyr)
library(presto)
library(cowplot)


# SET UP NAMES
housekeeping_genes <- c("")
genes_of_interest <- c("")
path_to_data <- "/sharedFolder/Data/"

name_new_dir_results <- paste(getwd(), "/Results", sep = "")
if (!dir.exists(name_new_dir_results)) {
    dir.create(name_new_dir_results)
}

name_new_dir_partial <- paste(getwd(), "/Partial", sep = "")
if (!dir.exists(name_new_dir_partial)) {
    dir.create(name_new_dir_partial)
}


# LOAD DATA, NORMALIZE, FIND VARIABLE FEATURES, SCALE DATA
load.data <- function(name_of_the_data) {
    # Load the data
    sc_data <- Read10X(data.dir = paste(path_to_data, #    paste() is the command to unite two strings
                                        name_of_the_data, 
                                        sep = ""), #       This specifies which caracter to put between the strings, with "" it puts nothing, by default it ads a space. paste0() is a paste() but with sep = "" as default
                       gene.column = 1) #    Some dataset have double names for the genes, both the Gene cards symbol (ex.  SRCIN1) and the Ensembl code (ex. ENSG00000277363). With gene.column you select wich one to use, in this case i used gene card 

    # Create Seurat object
    sc_data <- CreateSeuratObject(counts = sc_data, 
                                  min.cells = 3, 
                                  min.features = 500, 
                                  project = timepoints[time_point], 
                                  names.delim = "-", 
                                  names.field = 2)

    # Normalize the data
    sc_data <- NormalizeData(sc_data, 
                             normalization.method = "LogNormalize", 
                             scale.factor = 1e6)

    # Find variable features
    sc_data <- FindVariableFeatures(sc_data, 
                                    selection.method = "mvp", 
                                    nfeatures = 2000)

    # Scale the data
    sc_data <- ScaleData(sc_data)

    # Save the Scaled data
    name_new_dir <- paste(name_new_dir_partial, "/", name_of_the_data, sep = "")
    if (!dir.exists(name_new_dir)) {
        dir.create(name_new_dir)
    }
    print(paste("Saving PCA for time point", name_of_the_data, "in", name_new_dir))
    save(sc_data, 
         file = paste(name_new_dir, 
                      "/Scaled_", 
                      name_of_the_data, 
                      ".Robj", 
                      sep = ""))

    return(sc_data)
}


# PCA&Clusterization
PCA.cluster <- function(res) {
    print(paste("Running PCA and clustering for time point:", timepoints[time_point]))

    # PCA
    sc_data <- RunPCA(sc_data, npcs = 50, verbose = FALSE)
    # print(ElbowPlot(object = sc_data, ndims = 50))

    # Cluster the cells
    sc_data <- FindNeighbors(sc_data, dims = 1:40)
    sc_data <- FindClusters(sc_data, resolution = res)

    print(table(Idents(sc_data)))

    # Save the PCA plot
    name_new_dir <- paste(name_new_dir_partial, "/", timepoints[time_point], sep = "")
    if (!dir.exists(name_new_dir)) {
        dir.create(name_new_dir)
    }
    print(paste("Saving PCA for time point", timepoints[time_point], "in", name_new_dir))
    save(sc_data, file = paste(name_new_dir, "/PCA_res_", res, "_", timepoints[time_point], ".Robj", sep = ""))

    return(sc_data)
}


# FIND ALL MARKERS
cluster.markers <- function() {
    print(paste("Finding all markers for time point:", timepoints[time_point]))

    # Find all markers for every cluster compared to all remaining cells
    cluster_markers <- FindAllMarkers(sc_data,
        only.pos = TRUE, # Considera solo i marker espressi positivamente
        min.pct = 0.25, # Percentuale minima di espressione nelle cellule del cluster
        logfc.threshold = 0.25
    ) # Soglia minima di LogFC

    # Save the markers
    name_new_dir <- paste(name_new_dir_partial, "/", timepoints[time_point], sep = "")
    if (!dir.exists(name_new_dir)) {
        dir.create(name_new_dir)
    }
    print(paste("Saving cluster markers for time point", timepoints[time_point], "in", name_new_dir))
    save(cluster_markers, file = paste(name_new_dir, "/cluster_markers_", timepoints[time_point], ".Robj", sep = ""))

    return(cluster_markers)
}

# FIND CONSERVED MARKERS____________________________________________________________________________________________________________

# RELOAD DATA
load.sc_data <- function(time_point) {
    name_new_dir <- paste(name_new_dir_partial, "/", timepoints[time_point], sep = "")
    load(paste(name_new_dir, "/cluster_markers_", timepoints[time_point], ".Robj", sep = ""))
    return(sc_data)
}

load.clusters <- function(time_point, res) {
    name_new_dir <- paste(name_new_dir_partial, "/", timepoints[time_point], sep = "")
    load(paste(name_new_dir, "/PCA_res_", res, "_", timepoints[time_point], ".Robj", sep = ""))
    return(sc_data)
}

load.markers <- function(time_point) {
    name_new_dir <- paste(name_new_dir_partial, "/", timepoints[time_point], sep = "")
    load(paste(name_new_dir, "/cluster_markers_", timepoints[time_point], ".Robj", sep = ""))
    return(cluster_markers)
}


# FIND DIFFERENTIALLY EXPRESSED GENES
de.genes <- function(genes_oi) {
    print(paste("Finding differentially expressed genes for time point:", timepoints[time_point]))

    # Find differentially expressed genes
    de_genes <- cluster_markers %>% filter(gene %in% genes_of_interest)
    print(de_genes)

    # Save the DE genes
    name_new_dir <- paste(name_new_dir_results, "/", timepoints[time_point], sep = "")
    if (!dir.exists(name_new_dir)) {
        dir.create(name_new_dir)
    }
    print(paste("Saving differentially expressed genes for time point", timepoints[time_point], "in", name_new_dir))
    write.csv(de_genes, file = paste(name_new_dir, "/de_genes_", timepoints[time_point], ".csv", sep = ""))

    return(de_genes)
}


# ANNOTATION OF THE CLUSTERS
# EnrichR


# GPTCelltype
# Paper:            https://doi.org/10.1038/s41592-024-02235-4
# GitHub:           https://github.com/Winnie09/GPTCelltype
# To install:       install.packages("openai")
#                   remotes::install_github("Winnie09/GPTCelltype")


# Monocle


# PLOTS!!!!!!!!!!!!!!!!!!!!!!!!!!!

#   distribution of the genes of interest in the cells of the various clusters
#   compare with the distribution of the housekeeping genes

























# MAIN FUNCTION

# Load the data
sc_data <- load.data(1)
time_point <- sc_data$time_point
sc_data <- sc_data$sc_data

# Run PCA and clustering
sc_data <- PCA.cluster(sc_data, res = 1)
sc_data <- PCA.cluster(res = 1)

# Find all markers
cluster_markers <- cluster.markers(sc_data)
sluster_markers <- cluster.markers()

# Find differentially expressed genes
de_genes <- de.genes(genes_of_interest)































# Plotting
plotting <- function(sc_data, genes_oi) {
    # Plot the expression of the genes of interest
    for (gene in genes_oi) {
        p1 <- FeaturePlot(sc_data, features = gene, cols = c("lightgrey", "blue"), pt.size = 0.5) +
            ggtitle(paste("Expression of", gene)) +
            theme_minimal()

        # Save the plot
        ggsave(filename = paste(name_new_dir_partial, "/FeaturePlot_", gene, "_", timepoints[time_point], ".png", sep = ""), plot = p1)
    }
}


do_iteration(1)
# Plot the expression of the housekeeping genes
plotting(sc_data, housekeeping_genes)
# Plot the expression of the genes of interest
plotting(sc_data, genes_of_interest)
# Plot the expression of the housekeeping genes and genes of interest together
plotting(sc_data, c(housekeeping_genes, genes_of_interest))
