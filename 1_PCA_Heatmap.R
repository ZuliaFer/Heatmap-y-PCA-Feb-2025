# Autor: M. C. Zulia Fernandina Nieves López 
# Visualización de datos de mRNA-seq mediante PCA (Análisis de componentes principales)

# Cargar librerías necesarias para nuestro PCA
library(DESeq2)
library(ggplot2)
library(magrittr)
library(ggplotify)
library(ggrepel)
library(pheatmap)
library(dplyr)

# Datos crudos: Raíz (dos librería control/tratamiento), Hojas (dos librería control/tratamiento)
head(datos)


# Creación del objeto datosm (metadatos) con las réplicas técnicas
datosm <- data.frame(sample.id = c("Control_1_Raiz", "Control_2_Raiz", 
                "Trat_1_Raiz", "Trat_2_Raiz", 
                "Control_1_Hojas", "Control_2_Hojas", 
                "Trat_1_Hojas", "Trat_2_Hojas"),
  # Las réplicas técnicas dentro de cada condición experimental
  num.tech.reps = rep(1:2, 4),  # 1 para la primera réplica, 2 para la segunda réplica
  # El factor de condición experimental (control vs. tratamiento)
  strain = c("Control", "Control", "Trat", "Trat", 
             "Control", "Control", "Trat", "Trat"),
  # Las condiciones experimentales (en este caso las raíces o las hojas)
  experiment.number = c(1, 2, 1, 2, 1, 2, 1, 2),  
  # El tipo de tratamiento (Mock o Tratamiento)
  Experiment = c("Mock", "Mock", "Tratamiento", "Tratamiento", 
                 "Mock", "Mock", "Tratamiento", "Tratamiento"), 
  stringsAsFactors = FALSE
)
# Agregar una nueva columna 'tissue' que distingue entre raíces y hojas
datosm$tissue <- ifelse(grepl("Raiz", datosm$sample.id), "Raiz", "Hoja")


# Verifica el objeto datosm


# Utilizando el objeto "metadata" asignamos como nombre de las filas la columna "sample.id"


# Confirmamos que tenemos todas las condiciones en nuestro objeto "metadata"

# Normalizar y transformar las cuentas de datos usando "default desing"
# Trabajaremos con una matriz

# Vamos a factorizar el contenido de la columna "strain" desde nuestro objeto "metadata"

############---Ahora vamos a crear un objeto DESeqDataSet que necesitaremos para crear nuestro PCA---#############
# Crear el objeto DESeqDataSet considerando las réplicas técnicas
dds <- DESeqDataSetFromMatrix(
  countData = round(counts.matrix),  # Tu matriz de conteos (filas=genes, columnas=muestras)
  colData = datosm,  # El data frame con los metadatos
  design = ~ strain + experiment.number  # Considera tanto la condición como las réplicas técnicas
)

# Crear el objeto dds para DESeq2, usando la columna 'tissue' como parte del diseño


# Transformar los datos de cuentas a una forma útil para su visualización
# VST para normalizar las cuentas. Es más rápido


# rlog para normalizar las cuentas


# Obtener la matriz de los valores transformados 


# PCA usando la funcion plotPCA de DESeq2 
# ntop=500 indica que solo se deben considerar los 500 genes con mayor varianza en los datos de expresión para realizar el PCA



##############################################################################
######-- Ahora vamos a construir un plot PCA mejor usando ggplot2 --##########
# Primero vamos a obtener los datos de PCA sin generar un gráfico

# Cambiar ntop a 1000


##### Calcular el porcentaje de varianza explicado por cada componente princial (PC)
percentVar<- round(100 * attr(pcaData, "percentVar"))

#### Crear un PCA básico con ggplot
g<- ggplot(pcaData,aes(PC1,PC2))
show(g)


# Cambiar el color de las muestras
library(RColorBrewer)
display.brewer.all()


g3 <- g + geom_point(size=4, aes(col=tissue,fill=sample.id), shape = 21, stroke = 1) +
  # Colores para "tissue"
  scale_color_manual(values = c("Raiz" = "#E31A1C", "Hoja" = "#FEB24C")) +  
  # Colores para cada muestra
  scale_fill_manual(values = c("Control_1_Raiz" = "#1D91C0", "Control_2_Raiz" ="#1D91C0", 
                               "Trat_1_Raiz" = "#081D58", "Trat_2_Raiz" = "#081D58",
                               "Control_1_Hojas" = "#B2DF8A", "Control_2_Hojas" = "#B2DF8A", 
                               "Trat_1_Hojas" = "#33A02C", "Trat_2_Hojas" = "#33A02C"))
g3 + theme_bw()

# Ahora le vamos a agregar el nombre a cada muestra dentro de nuestro PC2


# Añadir un titulo y cambiar las etiquetas del eje x y y
g5<- g4 + theme_bw() + ggtitle("PCA con todas las librerías (lecturas crudas)") +
  xlab(paste0("PC1:",percentVar[1], "% variance")) +
  ylab(paste0("PC2:",percentVar[2], "% variance"))
show(g5)

# Si queremos separar el gráfico por tipo de tejido



#######-- Construyendo un Heatmap con el paquete pheatmap --########
# Primero, vamos a obtener una matriz desde el objeto transformado DESeq2
myMatrix<- SummarizedExperiment::assay(vsd)

# Obtener las ubicaciones de los genes más variables (n=500) y hacer un subconjunto de la matriz
topVarGenes<- head(order(rowVars(myMatrix, useNames = TRUE), decreasing = TRUE),500)

# Ahora vamos a ordenar nuestras muestras 



# Ahora vamos a dibujar un heatmap utilizando los 500 genes más variables seleccionados en el paso anterior.
phm<- pheatmap(myMatrix_for_heatmap_order,
               #Escalar los datos por filas
               scale = "row",
               # Agrupar las filas y columnas 
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               # Mostrar/ocultar los genes o las muestras gene or sample names
               show_rownames = FALSE,
               show_colnames = TRUE)
show(phm)

# Podemos modifcar algunas caracteristicas del plot


# Como luciria un heatmap pero con 1000 genes?



###-- Añadir el nombre de las muestras a los puntos en el plot
# colours
datosm$Experiment<- c("Control_Raiz","Control_Raiz","Tratamiento_Raiz", "Tratamiento_Raiz", 
                      "Control_Hojas", "Control_Hojas", "Tratamiento_Hojas", "Tratamiento_Hojas")


ann_colors <- list(
  strain = c("Control" = "#807DBA",
                     "Trat" = "#DE77AE"), 
  Experiment = c("Control_Hojas" = "darkgreen",
            "Control_Raiz" = "#377EB8",
            "Tratamiento_Raiz"="#FF7F00",
            "Tratamiento_Hojas"="#E41A1C")
)


phm2<- pheatmap(myMatrix_for_heatmap_order,
               #scale the data by row
               scale = "row",
               main="Draw heatmap based on the most variable genes",
               legend_labels = c("Low", "Medium", "High"),
               #cluster rows and columns are use default order
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               
               # show/hide gene or sample names
               show_rownames = FALSE,
               show_colnames = TRUE,
               
               # column annotation
               annotation_col = datosm[,c("strain", "Experiment"), drop=FALSE],
               annotation_colors = ann_colors)
show(phm2)

###-- Extracción de la informacion del objeto pheatmap.
# Utilizaremos str para ver la estructura del objeto phm 
str(phm2, list.len = 2)
# El simbolo $ puede ser utilizado para extraer los elementos de la lista
phm2$tree_row$labels[phm2$tree_row$order] %>% head()
# También podemos obtener las columnas en el mismo orden en el plot 
phm2$tree_col$labels[phm2$tree_col$order] %>% head()







