# go to the directory you downloaded
setwd("the.path.to.dir")
# define path to usefull directories
data.dir <- "./data"
RCC.dir <- file.path(data.dir,"GSE146204_RAW")


# install packages
source(file.path("./scripts","setup.R"))

# source functions
source(file.path("./scripts","functions.R"))

# use GEOQuery to upload annotation data
GOIDs <- c("GSE146204")
GEOdata <- getGEO(GOIDs,GSEMatrix=T, AnnotGPL=F)
sample.annots <- get.annotations(GEOdata)
sample.annots <- sample.annots[,c("title","tp53 status")]


# get the names of the RCC files containing the raw data
rcc.files <- system(paste("ls",RCC.dir),intern=T)
# map file names with sample IDs
samp.GSM <- unlist(sapply(strsplit(rcc.files,"_"),"[[",1))
m <- match(samp.GSM,row.names(sample.annots))
sample.annots$fileName <- sub(".RCC","",rcc.files)
colnames(sample.annots) <- sub(" ",".",colnames(sample.annots))
sample.annots$tp53.status <- sub("Wild type/NA/ND","WildType",sample.annots$tp53.status)
# read RCC files from your computer
raw.data <- import.raw.rcc(RCC.dir,sample.annots)


# options for tools :
# nappa.NS, nappa.param1, nappa.param2,nanostringnorm,desq2,nanoR.top100,nanostringR,nanoR.total

# compute all comparisons: each tool will normalize the data and compute differentation expression
# The function returns the data.frame with 2 columns : the log Fold Change and the SYMBOL of the gene
tools <- c("nappa.NS","nappa.param1", "nappa.param2","nappa.param3","nanostringnorm.default","nanostringnorm.param1","nanostringnorm.param2","desq2","nanoR.top100","nanoR.total","nanostringR")
data.to.comp <- tools.inspect(raw.data,tool="nappa.default",nanoR=F)
for (tool in tools){
  print(tool)
  nanoR=F
  raw <- raw.data
  if (tool%in%c("nanoR.top100","nanoR.total")){
    nanoR <- T
    raw <- RCC.dir
  }
  tmp <- tools.inspect(raw,tool,nanoR)
  data.to.comp <- merge(data.to.comp,tmp,by="SYMBOL",all=T)
}

row.names(data.to.comp) <- data.to.comp$SYMBOL
data.to.comp <- data.to.comp[,-1]

# if you need to remove genes not present in all analyses: NA (missing) data
data.to.comp <- na.omit(data.to.comp)

# transpose the matrix to use the genes as descriptive values to compare the tools
data.to.comp <- as.data.frame(t(data.to.comp))

################################################################################
################################################################################
# compute PCA

################################################################################
################################################################################
# compute intersection with upsetR






