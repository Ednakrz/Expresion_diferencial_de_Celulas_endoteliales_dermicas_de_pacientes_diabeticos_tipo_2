# Análisis de la expresión diferencial de células endoteliales dérmicas de pacientes diabéticos tipo 2

## Autor: 
### Edna Karen Rivera Zagal  
  ednakrz@lcg.unam.mx 

## Descripción del Proyecto

Este proyecto tiene como objetivo analizar la expresión diferencial de genes en células endoteliales de personas diabéticas y sanas, utilizando datos de RNA-seq obtenidos a través del paquete [Recount3](https://rna.recount.bio/) de [bioconductor](https://bioconductor.org/). Utilizamos herramientas bioinformáticas para identificar genes diferencialmente expresados y visualizamos los resultados con gráficos y heatmaps.

Para mayor información del proyecto, favor de leer el [reporte](https://github.com/Ednakrz/Expresion_diferencial_de_Celulas_endoteliales_dermicas_de_pacientes_diabeticos_tipo_2/blob/master/Reporte_Analisis_Diferencial.pdf) proporcionado en este repocitorio.  

## Datos

- Número de identificación de los datos en Recount3 : **SRP095512**
- [GEO accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92724) a los datos : **GSE92724**  

## Uso 
Para obtener los resultados del análisis por favor corra los scripts en el orden indicado: 

#### 1. Obtención de Datos: Ejecuta el script [obtencion_analisis_datos.R](https://github.com/Ednakrz/Expresion_diferencial_de_Celulas_endoteliales_dermicas_de_pacientes_diabeticos_tipo_2/blob/master/src/Obtencion_analisis_datos.R) para descargar los datos de RNA-seq, y depurarlos para su análisis. 

#### 2. Normalización de Datos: Ejecuta el script [Normalizacion_datos.R](https://github.com/Ednakrz/Expresion_diferencial_de_Celulas_endoteliales_dermicas_de_pacientes_diabeticos_tipo_2/blob/master/src/Normalizacion_datos.R) para normalizar los datos.

#### 3. Análisis de Expresión Diferencial: Ejecuta el script [Analisis_Expresion_diferencial.R](https://github.com/Ednakrz/Expresion_diferencial_de_Celulas_endoteliales_dermicas_de_pacientes_diabeticos_tipo_2/blob/master/src/Analisis_Expresion_diferencial.R) para realizar el análisis de expresión diferencial y obtener los gráficos resultantes.

## Licencia

Este proyecto está licenciado bajo la Licencia MIT. Para más detalles, consulta el archivo [LICENSE.txt](https://github.com/Ednakrz/Expresion_diferencial_de_Celulas_endoteliales_dermicas_de_pacientes_diabeticos_tipo_2/blob/master/LICENSE.txt).


## Citación 
 
Si tú usas este script en tu trabajo, por favor citalo como: Rivera-Zagal, E. Expresion_diferencial_de_Celulas_endoteliales_dermicas_de_pacientes_diabeticos_tipo_2 (v1.0.0) [GitHub Repository]. GitHub. https://github.com/Ednakrz/Expresion_diferencial_de_Celulas_endoteliales_dermicas_de_pacientes_diabeticos_tipo_2/tree/master

## Contactanos

Si tiene problemas o preguntas, favor de crear un issue en este repositorio o contactenos a través del correo: [ednakrz@lcg.unam.mx]





