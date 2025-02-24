install.packages("dplyr",dependencies=TRUE)
install.packages("rrBLUP")
install.packages("vcfR")

install.packages("shiny")
install.packages("NAM")
install.packages("bWGR")
install.packages("STPGA")
install.packages("BGLR")
install.packages("SOMMER")
install.packages("DT")

install.packages("BiocManager")
BiocManager::install(version = "3.16")
 
options(repos = BiocManager::repositories())
library(BiocManager)

install.packages("crayon")
install.packages("devtools")
library(devtools)

install.packages("rJava")

Sys.setenv(LD_LIBRARY_PATH='/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/java-17-openjdk-amd64/lib/server/')
Sys.getenv("LD_LIBRARY_PATH")
dyn.load("/usr/lib/jvm/java-17-openjdk-amd64/lib/server/libjvm.so")
library(rJava)


devtools::install_bitbucket(
		repo = "bucklerlab/rTASSEL",
		ref = "master",
		build_vignettes = FALSE
) 
library(rTASSEL)




