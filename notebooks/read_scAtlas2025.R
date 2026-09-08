
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

writeMM(cnt, "atlas_counts.mtx")
writeLines(rownames(cnt), "genes.txt")
writeLines(colnames(cnt), "cells.txt")
write.csv(md, "atlas_meta.csv")

rm(obj); gc()                        # release the giant object immediately

stopifnot(all(cnt@x == floor(cnt@x)))    # raw UMIs, not normalized

sapply(md, function(x) if (is.character(x) || is.factor(x)) length(unique(x)) else NA)   # which cols are categorical

#             nCount_RNA            nFeature_RNA              percent.mt 
#                     NA                      NA                      NA 
#        seurat_clusters                   Count  Study..Citation..PMID. 
#                     43                      NA                      12 
#        GSE.SRA..Study.                    Name If.metastatic..location 
#                     12                     234                       7 
#               Clusters               Treatment            DiseaseState 
#                     14                       4                       4 
#          TreatmentType 
#                     14 


# There are your columns. The study lives in Study..Citation..PMID. or GSE.SRA..Study. (both 12 levels — one is citation/PMID, 
# the other the accession),
# and state is DiseaseState (4 levels: expect PT/Met/AdjN/Donor). 
# Name (234) is the per-sample ID — that's your donor count column. If.metastatic..location (7) confirms the met sites.

# Print the actual level strings so you use the exact spelling:

sort(unique(md$GSE.SRA..Study.)) 

'''
 [1] "EGAS00001002543" "GSE154778"       "GSE155698"       "GSE156405"      
 [5] "GSE158356"       "GSE194247"       "GSE202051"       "GSE205013"      
 [9] "GSE211644"       "GSE229413"       "phs001840.v1.p1" "PRJCA001063" 
'''

sort(unique(md$Study..Citation..PMID.))
'''
 [1] "Carpenter et al Cancer Discovery"                                                                                                                                                                                                                                                                                                                                                                                                                                                         
 [2] "Chan-Seng-Yue, M., Kim, J.C., Wilson, G.W. et al. Transcription phenotypes of pancreatic cancer are driven by genomic events during tumor evolution. Nat Genet 52, 231–240 (2020). https://doi.org/10.1038/s41588-019-0566-9 PMID: 31932696"                                                                                                                                                                                                                                              
 [3] "Elyada E, Bolisetty M, Laise P, Flynn WF, Courtois ET, Burkhart RA, Teinor JA, Belleau P, Biffi G, Lucito MS, Sivajothi S, Armstrong TD, Engle DD, Yu KH, Hao Y, Wolfgang CL, Park Y, Preall J, Jaffee EM, Califano A, Robson P, Tuveson DA. Cross-Species Single-Cell Analysis of Pancreatic Ductal Adenocarcinoma Reveals Antigen-Presenting Cancer-Associated Fibroblasts. Cancer Discov. 2019 Aug;9(8):1102-1123. doi: 10.1158/2159-8290.CD-19-0094. Epub 2019 Jun 13. PMID: 31197017"
 [4] "Kemp et al., Life Sci Alliance, 2020"                                                                                                                                                                                                                                                                                                                                                                                                                                                     
 [5] "Lee et al Clinical Cancer Research"                                                                                                                                                                                                                                                                                                                                                                                                                                                       
 [6] "Lin W, Noel P, Borazanci EH, Lee J et al. Single-cell transcriptome analysis of tumor and stromal compartments of pancreatic ductal adenocarcinoma primary tumors and metastatic lesions. Genome Med 2020 Sep 29;12(1):80."                                                                                                                                                                                                                                                               
 [7] "Not published yet"                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
 [8] "Peng J, Sun BF, Chen CY, Zhou JY, Chen YS, Chen H, Liu L, Huang D, Jiang J, Cui GS, Yang Y, Wang W, Guo D, Dai M, Guo J, Zhang T, Liao Q, Liu Y, Zhao YL, Han DL, Zhao Y, Yang YG, Wu W. Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma. Cell Res. 2019 Sep;29(9):725-738. doi: 10.1038/s41422-019-0195-y. Epub 2019 Jul 4. Erratum in: Cell Res. 2019 Aug 13;: PMID: 31273297"                                 
 [9] "Schalck A, Sakellariou-Thompson D, Forget MA, Sei E et al. Single-Cell Sequencing Reveals Trajectory of Tumor-Infiltrating Lymphocyte States in Pancreatic Cancer. Cancer Discov 2022 Oct 5;12(10):2330-2349. PMID: 35849783"                                                                                                                                                                                                                                                             
[10] "Shiau C, Su J, Guo JA, Hong TS et al. Treatment-associated remodeling of the pancreatic cancer endothelium at single-cell resolution. Front Oncol 2022;12:929950. PMID: 36185212"                                                                                                                                                                                                                                                                                                         
[11] "Steele NG, Carpenter ES, Kemp SB, Sirihorachai VR et al. Multimodal Mapping of the Tumor and Peripheral Blood Immune Landscape in Human Pancreatic Cancer. Nat Cancer 2020 Nov;1(11):1097-1112. PMID: 34296197"                                                                                                                                                                                                                                                                           
[12] "Werba G, Weissinger D, Kawaler EA, Zhao E et al. Single-cell RNA sequencing reveals the effects of chemotherapy on human pancreatic adenocarcinoma and its tumor microenvironment. Nat Commun 2023 Feb 13;14(1):797. PMID: 36781852" 
'''

sort(unique(md$DiseaseState))
'''
[1] Donor             Adjacent normal   Primary tumor     Metastatic lesion
Levels: Donor Adjacent normal Primary tumor Metastatic lesion
'''

table(md$GSE.SRA..Study., md$DiseaseState)

'''
                
                   Donor Adjacent normal Primary tumor Metastatic lesion
  EGAS00001002543      0               0         76094                 0
  GSE154778            0               0         11806              8558
  GSE155698            0            4502         19489                 0
  GSE156405            0               0         10011              7798
  GSE158356            0               0             0              2746
  GSE194247            0               0         29622                 0
  GSE202051            0           10660        136983                 0
  GSE205013            0               0        130026             37340
  GSE211644            0               0         40971                 0
  GSE229413        33309               0             0                 0
  phs001840.v1.p1      0            4181         14104                 0
  PRJCA001063          0           54692         93215                 0
'''

'''
Numbers check out — the table sums to exactly 726,107, so nothings miscounted. Dropping Peng (PRJCA001063) and all Metastatic lesion cells:

Peng removes 147,907 (54,692 AdjN + 93,215 PT)
Mets remove 56,442 (matches the papers Met total)
Remaining: 521,758 cells across 10 studies

Two studies vanish entirely: PRJCA001063 (Peng) and GSE158356, which is 100% metastatic (Schalck TIL study, 2,746 cells).

Accession	Study	Kept cells	States kept
EGAS00001002543	Chan-Seng-Yue 2020	76,094	PT
GSE154778	Lin 2020	11,806	PT (8,558 met dropped)
GSE155698	Steele 2020	23,991	AdjN + PT
GSE156405	Lee 2021	10,011	PT (7,798 met dropped)
GSE194247	(verify)	29,622	PT
GSE202051	(verify — see below)	147,643	AdjN + PT
GSE205013	Werba 2023	130,026	PT (37,340 met dropped)
GSE211644	(verify)	40,971	PT
GSE229413	not published	33,309	Donor
phs001840.v1.p1	Elyada 2019	18,285	AdjN + PT

Correction to something I told you earlier. 
Id said GSE202051 was Hwang snRNA and to drop it. 
The citation list here has no Hwang entry, and Lovelesss methods state every included study is 10x scRNA with raw fastqs re-aligned uniformly. 
So my snRNA flag was probably wrong — and its your single largest remaining contributor (147,643 cells, 28% of the reference), 
so dont drop it on my earlier say-so. Resolve its identity before deciding:

'''

# map accession -> citation -> confirm platform
table(md$GSE.SRA..Study., md$Study..Citation..PMID.)

STUDY <- "GSE.SRA..Study."
STATE <- "DiseaseState"

drop_study <- c("PRJCA001063")           # Peng
drop_state <- c("Metastatic lesion")     # exact factor level

keep <- !(md[[STUDY]] %in% drop_study) & !(md[[STATE]] %in% drop_state)
cat(sum(keep), "cells kept of", nrow(md), "\n")     # expect 521758
print(table(droplevels(md[[STUDY]][keep]), droplevels(md[[STATE]][keep])))

# 521758 cells kept of 726107 

# donors remaining per study — the number that governs reference quality
tapply(md$Name[keep], md[[STUDY]][keep], function(x) length(unique(x)))
'''
EGAS00001002543       GSE154778       GSE155698       GSE156405       GSE194247 
             13               9              20               5               5 
      GSE202051       GSE205013       GSE211644       GSE229413 phs001840.v1.p1 
             71              17              14              11               9 
'''
print(table(md[[STUDY]][keep], droplevels(md[[STATE]][keep])))
'''        
                   Donor Adjacent normal Primary tumor
  EGAS00001002543      0               0         76094
  GSE154778            0               0         11806
  GSE155698            0            4502         19489
  GSE156405            0               0         10011
  GSE194247            0               0         29622
  GSE202051            0           10660        136983
  GSE205013            0               0        130026
  GSE211644            0               0         40971
  GSE229413        33309               0             0
  phs001840.v1.p1      0            4181         14104
'''

print("\n------------------ end -----------------")


sort(table(md$Clusters), decreasing = TRUE)
grep("type|ident|label|celltype|annotation|compartment",
     colnames(md), ignore.case = TRUE, value = TRUE)

'''
          DUCTAL              TNK          MYELOID      FIBROBLASTS 
          277301           102067            76253            66192 
     ENDOTHELIAL        PERICYTES           ACINAR   CYCLING DUCTAL 
           52388            35774            35109            23276 
       ENDOCRINE      CYCLING TNK          B CELLS           PLASMA 
           17833            14607            12880             6117 
            MAST CYCLING. MYELOID 
            3688             2622 
[1] "TreatmentType"

'''


#---------------- Final save ------------------------

# cnt = genes x cells, raw counts (extract if not already):
# cnt <- GetAssayData(obj, layer = "counts")   # or slot = "counts" on pre-v5

stopifnot(identical(colnames(cnt), rownames(md)))   # matrix/meta aligned before subsetting

keep <- !(md$GSE.SRA..Study. %in% "PRJCA001063") &
        !(md$DiseaseState %in% "Metastatic lesion")
cat(sum(keep), "cells kept\n")                       # expect 521758

# 521758 cells kept

cnt_sub <- cnt[, keep]                               # sparse column slice — cheap
md_sub  <- md[keep, , drop = FALSE]

stopifnot(identical(colnames(cnt_sub), rownames(md_sub)))          # alignment holds
stopifnot(all(cnt_sub@x == floor(cnt_sub@x)))                      # raw integer UMIs

writeMM(cnt_sub,               "atlas_counts_noPeng_noMet.mtx")
writeLines(rownames(cnt_sub),  "genes_noPeng_noMet.txt")
writeLines(colnames(cnt_sub),  "cells_noPeng_noMet.txt")
write.csv(md_sub,              "atlas_meta_noPeng_noMet.csv")
