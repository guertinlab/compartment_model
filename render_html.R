if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
if(!require("rmarkdown", quietly = TRUE))
    install.packages("rmarkdown", repos="https://cloud.r-project.org")

BiocManager::install(version = "3.15")
BiocManager::install("BiocStyle")
# pandoc needs to be installed. Follow this page: https://pandoc.org/installing.html

rmarkdown::render("compartmentModel_efferocytosis.Rmd", BiocStyle::html_document())

