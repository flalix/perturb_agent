
#!/usr/bin/env Rscript
#
# read_scAtlas2025.R
#
# @Author:
#    Flavio Lichtenstein
#    Claude: Opus 4.8
# @Date: 2026/09/07

# conda deactivate

# libpng pulls in the headers the png R package wants; zlib provides libz for the -lz link
# conda install -n renv -c conda-forge zlib libpng

# conda activate renv
# export LDFLAGS="-L${CONDA_PREFIX}/lib"
# export CPPFLAGS="-I${CONDA_PREFIX}/include"
# export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

# A faster and usually less painful alternative — install the heavy binary deps from conda-forge instead of compiling them, then let remotes build only Seurat:

# conda install -n renv -c conda-forge \
#   r-reticulate r-leiden r-png r-rcppeigen r-rcppprogress \
#   r-matrix r-seuratobject

# conda install -n renv -c conda-forge r-seuratobject

# Option 1 — enlarge swap, load once, extract, never load again (fastest to try)

# sudo fallocate -l 100G /swapfile2
# sudo chmod 600 /swapfile2
# sudo mkswap /swapfile2
# sudo swapon /swapfile2
# free -h

# Setting up swapspace version 1, size = 100 GiB (107374178304 bytes)
# no label, UUID=c02cab8e-a335-4a7a-a456-93f80456030a
#                total        used        free      shared  buff/cache   available
# Mem:            62Gi        18Gi        43Gi       905Mi       2,1Gi        43Gi
# Swap:          101Gi       1,9Gi       100Gi


# -------------------
# R version 4.5.3

library(SeuratObject)
library(Matrix)
## installed / loaded versions
packageVersion("SeuratObject")  # ‘5.4.0’
packageVersion("Matrix")        # ‘1.7.5’

setwd("/home/flavio/uv/perturb_agent/data/multi_progs/PAAD/scAtlas2025/prism")

obj <- readRDS("scAtlas.rds")            # after gunzip

## the version the object was created under (the one that matters)
obj@version                       # Seurat/SeuratObject version at save time

sessionInfo()   

print("\n-------------------\n")

print(obj)

# An object of class Seurat 
# 36601 features across 726107 samples within 1 assay 
# Active assay: RNA (36601 features, 5000 variable features)
#  3 layers present: counts, data, scale.data
#  3 dimensional reductions calculated: pca, harmony, umap

print("\n-------------------\n")

## 1. find the columns before assuming names
md <- obj@meta.data
print(str(md))   # dump all columns + example values
print("\n-------------------\n")

cnt <- GetAssayData(obj, layer = "counts")
print(str(cnt))   # counts
print("\n-------------------\n")

rm(obj); gc()                        # release the giant object immediately

writeMM(cnt, "atlas_counts.mtx")
writeLines(rownames(cnt), "genes.txt")
writeLines(colnames(cnt), "cells.txt")
write.csv(md, "atlas_meta.csv")


sapply(md, function(x) if (is.character(x) || is.factor(x)) length(unique(x)) else NA)   # which cols are categorical

if (3==5) {

    STUDY <- "orig.study"      # <- replace with the real column name
    STATE <- "disease_state"   # <- replace with the real column name

    table(md[[STUDY]], md[[STATE]])           # the counts you asked for
    sort(unique(md[[STUDY]]))                 # confirm how Peng is spelled
    sort(unique(md[[STATE]]))                 # confirm how Met is spelled


    keep <- !(md[[STUDY]] %in% peng_labels) & !(md[[STATE]] %in% met_labels)
    cat(sum(keep), "cells kept of", nrow(md), "\n")
    print(table(md[[STUDY]][keep], md[[STATE]][keep]))

    sub <- subset(obj, cells = rownames(md)[keep])

    ## Seurat v5 gotcha: counts may be split into per-sample layers.
    ## Join them before extracting, or GetAssayData returns empty/partial.
    if (inherits(sub[["RNA"]], "Assay5")) sub <- JoinLayers(sub)

    cnt <- GetAssayData(sub, assay = "RNA", layer = "counts")   # genes x cells, raw UMIs
    stopifnot(all(cnt@x == floor(cnt@x)))                       # verify integer counts

    writeMM(cnt, "atlas_counts_noPeng_noMet.mtx")
    writeLines(rownames(cnt), "genes.txt")
    writeLines(colnames(cnt), "cells.txt")
    write.csv(sub@meta.data, "atlas_meta_noPeng_noMet.csv")
    saveRDS(sub, "scAtlas_noPeng_noMet.rds")   # optional, much smaller
}


print("\n------------------ end -----------------")
