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


# Check quality of the data
rse_gene_SRP095512$assigned_gene_prop <-
  rse_gene_SRP095512$recount_qc.gene_fc_count_all.assigned /
  rse_gene_SRP095512$recount_qc.gene_fc_count_all.total
summary(rse_gene_SRP095512$assigned_gene_prop)

# Observe differences among samples injected with different substances
with(colData(rse_gene_SRP095512),
     tapply(assigned_gene_prop,
            sra_attribute.disease_state,
            summary))

# Plot the distribution of the assigned gene proportion
hist(rse_gene_SRP095512$assigned_gene_prop, col = "orange")

# Save the unfiltered data
rse_gene_SRP095512_unfiltered <- rse_gene_SRP095512

## ------------------- Normalizacion de los datos-------------------------------

library("edgeR") # BiocManager::install("edgeR", update = FALSE)
dge <- DGEList(
  counts = assay(rse_gene_SRP095512, "counts"),
  genes = rowData(rse_gene_SRP095512)
)
dge <- calcNormFactors(dge)


## ---------------- Expresion diferencial --------------------------------------

library("ggplot2")
ggplot(as.data.frame(colData(rse_gene_SRP095512)), aes(y = assigned_gene_prop, x = sra_attribute.disease_state)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Gender")

# Modelo estadístico

mod <- model.matrix(~ sra_attribute.gender +sra_attribute.disease_state + assigned_gene_prop,
                    data = colData(rse_gene_SRP095512)
)
colnames(mod)

## Usamos limma para realizar el análisis de expresión diferencial
library("limma")
vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP095512),
  sort.by = "none"
)
dim(de_results)

## Genes diferencialmente expresados entre pre y post natal con FDR < 5%
table(de_results$adj.P.Val < 0.05)

## Visualicemos los resultados estadísticos
plotMA(eb_results, coef = 3)

volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)
de_results[de_results$gene_name %in% c("TTY14", "RPS4Y1", "RP11-331"), ]


volcanoplot(eb_results, coef = 3, highlight = 3, names = de_results$gene_name)
de_results[de_results$gene_name %in% c("WNK2", "NUD111", "FAM83H"), ]


## Visualizando la expresion diferencial a traves de un Heatmap------

## Extraer valores de los genes de interés
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

## Creemos una tabla con información de las muestras
## y con nombres de columnas más amigables
df <- as.data.frame(colData(rse_gene_SRP095512)[, c("sra_attribute.gender", "sra_attribute.disease_state")])
colnames(df) <- c( "Gender", "disease_state")

## Hagamos un heatmap
library("pheatmap")
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df,
  scale = "row"
)

ComplexHeatmap::pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df,
  scale = "row"
)
## Para colores
library("RColorBrewer")

## Conviertiendo los valores de Gender a colores
col.gender <- df$sra_attribute.gender
levels(col.gender) <- brewer.pal(nlevels(col.gender), "Dark2")

col.gender <- as.character(col.gender)

## MDS por genero
plotMDS(vGene$E, labels = df$sra_attribute.gender, col = col.gender)

## Conviertiendo los valores de Disease a colores
col.disease <- df$sra_attribute.disease_state
levels(col.disease) <- brewer.pal(nlevels(col.disease), "Dark2")

col.disease <- as.character(col.disease)

## MDS por genero
plotMDS(vGene$E, labels = df$sra_attribute.disease_state, col = col.disease)



