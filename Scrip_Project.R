## Cargar paquete Recount3 R
library("recount3")

## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()


## Identificador de mi proyecto de interes:SRP095512
##Obtencion del proyecto de interes
proj_info <- subset(
  human_projects,
  project == "SRP095512" & project_type == "data_sources"
)


## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP095512 <- create_rse(proj_info)

## Explora el objeto RSE
rse_gene_SRP095512

## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
assay(rse_gene_SRP095512, "counts") <- compute_read_counts(rse_gene_SRP095512)

## Visualizar la informacion de los atributos del objeto "rse_gene_SRP095512"

rse_gene_SRP095512 <- expand_sra_attributes(rse_gene_SRP095512)
colData(rse_gene_SRP095512)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP095512)))
]

