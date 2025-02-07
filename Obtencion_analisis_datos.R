## ------------------------Obtención y análisis de los datos -------------------
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

## Formato correcto en el que debemos tener la informacion
## Pasar de character a numeric o factor
rse_gene_SRP095512$sra_attribute.disease_state <- factor(tolower(rse_gene_SRP095512$sra_attribute.disease_state))

rse_gene_SRP095512$sra_attribute.gender <- factor(rse_gene_SRP095512$sra_attribute.gender)

## Resumen de las variables de interés
summary(as.data.frame(colData(rse_gene_SRP095512)[
  ,
  grepl("^sra_attribute.[gender|disease_state]", colnames(colData(rse_gene_SRP095512)))
]))


# Chechar la calidad de los datos
rse_gene_SRP095512$assigned_gene_prop <-
  rse_gene_SRP095512$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP095512$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP095512$assigned_gene_prop)

# Observar la diferencia entre muestras
with(colData(rse_gene_SRP095512),
     tapply(assigned_gene_prop,
            sra_attribute.disease_state,
            summary))

# Plot the distribution of the assigned gene proportion
# Grafica de la distribución de la proporción de genes asignados
hist(rse_gene_SRP095512$assigned_gene_prop, col = "orange")

# Guardad los datos sin filtrar, en caso de querer recuperarlos despúes
rse_gene_SRP095512_unfiltered <- rse_gene_SRP095512
