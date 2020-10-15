packages.cran <- c("devtools","data.table","NanoStringNorm","nanostringr","ggrepel")
packages.bioC <- c("limma","GEOquery","DESeq2")

installed <- rownames(installed.packages())
packages.cran <- packages.cran[!packages.cran%in%installed]
packages.bioC <- packages.bioC[!packages.cran%in%installed]

if (length(packages.cran)>0){
  install.packages(packages.cran, dependencies = TRUE)
}

if (!requireNamespace("nanoR", quietly = TRUE))
  install.packages("nanoR/",repos=NULL,type="source")


if (length(packages.cran)>0){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(packages.bioC)
}

if (!requireNamespace("ggbiplot", quietly = TRUE)) {
  library(devtools)
  install_github("vqv/ggbiplot")
}

