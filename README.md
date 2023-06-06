# mND

mND uses **m**ulti-layer **N**etwork **D**iffusion to find gene networks that contain high scoring genes. It considers two or more "layers" of genome-wide scores and an interactome. 

<img src="vignettes/overview.jpeg">

image from: https://academic.oup.com/bioinformatics/article/36/3/865/5553095

## Citation

If you use this software, please cite:

Di Nanni N, Gnocchi M, Moscatelli M, Milanesi L and Mosca E, *Gene relevance based on multiple evidences in complex networks*, **Bioinformatics**, Volume 36, Issue 3, 1 February 2020, Pages 865–871, https://doi.org/10.1093/bioinformatics/btz652

## Installation
```{r, eval=FALSE}
install.packages(c("devtools", "igraph", "parallel"))
library(devtools)
install_github("emosca-cnr/mND", build_vignettes = TRUE)
```

## Documentation
Please look at the vignette included in the R package:
```{r, eval=FALSE}
library(mND)
vignette("mND")
```

