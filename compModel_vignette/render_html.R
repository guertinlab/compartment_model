
# UNCOMMENT if you don't have these packages installed
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager", repos="https://cloud.r-project.org")
#if(!require("rmarkdown", quietly = TRUE))
#    install.packages("rmarkdown", repos="https://cloud.r-project.org")

#BiocManager::install(ask = F)
#BiocManager::install("BiocStyle")
# pandoc needs to be installed. Follow this page: https://pandoc.org/installing.html

#install.packages("xfun")
library(BiocManager)
library(rmarkdown)
rmarkdown::render("compartmentModel_vignette.Rmd", BiocStyle::html_document())

