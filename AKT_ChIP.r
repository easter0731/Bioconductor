## ChIP seq analysis of AKT project using DiffBind Package

# Function to check bioconductor packages
installpkg <- function(pkg){
        if(!require(pkg,character.only=T)){
                source("http://bioconductor.org/biocLite.R")
                biocLite(pkg)
        }
        require(pkg,character.only=T)
}

# Setting R file folder as a working directory
setwd("F:/Chip - and RNA-seq/AKT Chip-Seq/#Bin Anaylsis")

# Backup information
tmp =  tempfile(as.character(Sys.getpid()))
pdf(tmp)
savewarn = options("warn")
options(warn=-1)

# Loading Package
installpkg("DiffBind")

# Get sample sheet
samples = read.csv("AKT_ChIP_samplesheet.csv")
samples
    


