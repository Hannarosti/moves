pkglist = c("amt","bayesplot","bookdown","broom","CircStats","circular","coda",
            "crawl","dplyr","ggExtra","ggmap","ggplot2","gridExtra","here",
            "knitr","leaflet","lme4","lubridate","lwgeom","maptools","MASS",
            "Matrix","MHadaptive","momentuHMM","move","moveHMM","mvtnorm",
            "nlme","raster","Rcpp","readr","rstan","sp","spatstat","tibble",
            "tidyr","tidyverse")

inst.pkgs = rownames(installed.packages())
newpkgs <- pkglist[!pkglist %in% inst.pkgs]
if (length(newpkgs) > 0) {
  do.call("install.packages", list(pkglist))
}