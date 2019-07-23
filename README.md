# SCTCwhatateam

## 1. Using package
Our tutorial is coming soon.

The Python functions (linGen, linFwd, linRev) could be retrieved by
```{r}
library(reticulate)
source_python("./linmethods.py")
```

The data example is provided by DREAM Single Cell Transcriptomics Challenge at https://www.dropbox.com/sh/q7fq3y5yw1lmtrq/AAAnVDVAWzdBRS6FJor-Bk1Aa?dl=0

## 2. Shiny app
To simplify we have developed an easy-to-use user interface at http://nugget.unisa.edu.au:3838/SCTCWhatATeam

What we need to input is dependent on which methods to use:

- Reference locations for gene markers (Locations x Genes): see bdtnp.csv

- Reference scRNA-Seq binarized file (Cells x Genes): see binarized_bdtnp.csv

- scRNA-Seq normalized file (Cells x Genes): see dge_normalized.csv

- 3D coordinates for each location: see geometry.csv

- Seed list: seed_list.csv
