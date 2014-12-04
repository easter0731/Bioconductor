# Function to check bioconductor packages
installpkg <- function(pkg){
        if(!require(pkg,character.only=T)){
                source("http://bioconductor.org/biocLite.R")
                biocLite(pkg)
        }else{
                require(pkg,character.only=T)
        }
}

installpkg('seqinr')
?require