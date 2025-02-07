## ---------------- Expresion diferencial --------------------------------------
## Cargar librerias
library("ggplot2")
library("limma")
library("pheatmap")
library("RColorBrewer")

## Creacion de un boxplot para analizar la distribución de los genes en muestras
## de pacientes enfermos y sanos

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

vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_SRP095512),
  sort.by = "none"
)
dim(de_results)

## Genes diferencialmente expresados con FDR < 5%
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

## MDS por estado de salud
plotMDS(vGene$E, labels = df$sra_attribute.disease_state, col = col.disease)
