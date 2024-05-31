# Move to the vignettes directory
old_dir <- setwd("vignettes")

# Install the latest version of the package
devtools::install()

# Generate static Rmds
knitr::knit("source/_ehymet-01-introduction.Rmd", "ehymet-01-introduction.Rmd")
knitr::knit("source/_ehymet-02-clustering-epi-hypo.Rmd", "ehymet-02-clustering-epi-hypo.Rmd")

# Restore
setwd(old_dir)
