
## ------------------- Normalizacion de los datos-------------------------------

## Cargar librerias
## Normalización de los datos
library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP095512, "counts"),
  genes = rowData(rse_gene_SRP095512)
)
dge <- calcNormFactors(dge)
