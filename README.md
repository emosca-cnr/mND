# mND

mND uses **m**ulti-layer **N**etwork **D**iffusion to find gene networks that contain high scoring genes. It considers two or more "layers" of genome-wide scores and an interactome. 

![Overview of mND](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/36/3/10.1093_bioinformatics_btz652/1/m_btz652f1.jpeg?Expires=1630659321&Signature=Bm7kxTRXzLzJgrT~lkODRftelkSmk9omb-kGJJPash5cL2NiA1TfVzfHnYYWT6Cq4PzBw7T7bvTSHdBWc5YnjEccVUnppU5CgdK6yFo1LQs8C8G1fkJ9HdOGFulAYRdIDG66Sqz0ngQRBVT5dLinh0TyEeO3HJ6gIwBQeWOopYx~bHYKIHKOXVmBICcJMSf9V5orHv3Q1aCMmJREWSPOi4zQWP5MOJb7z9Nu7X7vFaezxsb-DeawwVe52Zp8U9nftSr940QdNv5-ULDTjsNAosVXPlakXYmbmn4fayuCZpz6f8iv6PXTl7xJ9i~BgO4HyHQ5bMpzc2RTW80TWxEusA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

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

