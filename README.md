
# SCTCwhatateam: Predicting cell locations based on location-marker genes

### Supplementary folder: https://www.dropbox.com/sh/q7fq3y5yw1lmtrq/AAAnVDVAWzdBRS6FJor-Bk1Aa?dl=0

## 1. Using a package


The Python functions (linGen, linFwd, linRev) could be retrieved by running the linmethods.py provided in the Supplementary folder.


```{r}
library(reticulate)
# source_python("./linmethods.py")
```

The data example is provided by DREAM Single Cell Transcriptomics Challenge in the Supplementary folder.

## 2. Shiny app

Please see tutorial in the "Shiny app tutorial.docx" in the Supplementary folder.

To simplify we have developed an easy-to-use user interface at https://github.com/pvvhoang/SCTCwhatateam-ShinyApp

The input requirements are based on each method's preference. In general,

- Reference gene expression for locations (Locations x Genes): see bdtnp.csv

- Reference scRNA-Seq binarized file (Locations x Genes): see binarized_bdtnp.csv

- Raw scRNA-Seq expression (Cell x Genes): see dge_raw.csv

- scRNA-Seq normalized file (Cells x Genes): see dge_normalized.csv

- 3D coordinates for each locations (Locations x (x,y,z) coordinates): see geometry.csv

- Seed list (Locations x coordinates): seed_list.csv
