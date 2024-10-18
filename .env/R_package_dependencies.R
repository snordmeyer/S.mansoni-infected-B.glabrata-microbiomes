#System 
options("repos" = c(CRAN = "https://cloud.r-project.org"))

if (packageVersion("cli") < "3.1.0") { install.packages("cli") }

library("devtools")
library("rlang")

myinstall_version <- function(x, ...) { if(! is_installed(x)) install_version(x, ...) }
myinstall_bioc    <- function(x, ...) { if(! is_installed(basename(x))) install_bioc(x, ...) }
myinstall_github  <- function(x, ...) { if(! is_installed(basename(x))) install_github(x, ...) }

# Microbiome analysis
myinstall_version("AICcPermanova", version = "0.0.2", upgrade = "never")
myinstall_version("rjson", version = "0.2.21", upgrade = "never") # microViz dependency
myinstall_github("david-barnett/microViz", ref = "739f388", upgrade = "never")
myinstall_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", ref = "6e09713", upgrade = "never")
myinstall_github("mahendra-mariadassou/phyloseq-extended", ref = "093aac5", upgrade = "never")
myinstall_github("mikemc/speedyseq", ref = "0057652", upgrade = "never")
myinstall_github("gauravsk/ranacapa", ref = "58c0cab", upgrade = "never")
