# Function to check bioconductor packages
installpkg <- function(pkg){
        if(!require(pkg,character.only=T)){
                source("http://bioconductor.org/biocLite.R")
                biocLite(pkg)
        }
        require(pkg,character.only=T)
}

installpkg("DiffBind")

source("http://bioconductor.org/biocLite.R")
biocLite("GenomicAlignments")



library("GenomicAlignments")
choosebank("genbank")
