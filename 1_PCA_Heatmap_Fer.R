# Autor: M. C. Zulia Fernandina Nieves López 
# Visualización de datos de mRNA-seq mediante PCA (Análisis de componentes principales)

# Cargar librerías
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(magrittr)


# Datos para trabajar: tabla con el número de cuentas por librería (sin ningún filtro)
# Librerías: Raíz (dos librería control/tratamiento), Hojas (dos librería control/tratamiento)
head(datos)
dim(datos)

# Creación del objeto datosm con las réplicas técnicas
datosm <- data.frame(
  sample.id = c("Control_1_Raiz", "Control_2_Raiz", 
                "Trat_1_Raiz", "Trat_2_Raiz", 
                "Control_1_Hojas", "Control_2_Hojas", 
                "Trat_1_Hojas", "Trat_2_Hojas"),
  # Las réplicas técnicas dentro de cada condición experimental
  num.tech.reps = rep(1:2, 4),  # 1 para la primera réplica, 2 para la segunda réplica
  # El factor de condición experimental (control vs. tratamiento)
  strain = c("Control", "Control", "Trat", "Trat", 
             "Control", "Control", "Trat", "Trat"),
  # Las condiciones experimentales (en este caso las raíces o las hojas)
  experiment.number = c(1, 2, 1, 2, 1, 2, 1, 2),  # Número de réplica dentro de cada condición
  # El tipo de tratamiento (Mock o Tratamiento)
  Experiment = c("Mock", "Mock", "Tratamiento", "Tratamiento", 
                 "Mock", "Mock", "Tratamiento", "Tratamiento"), 
  stringsAsFactors = FALSE
)
# Agregar una nueva columna 'tissue' que distingue entre raíces y hojas
datosm$tissue <- ifelse(grepl("Raiz", datosm$sample.id), "Raiz", "Hoja")


# Verifica el objeto datosm
head(datosm)

# Utilizando el objeto "metadata" asignamos como nombre de las filas la columna "sample.id"
rownames(datosm)<-datosm$sample.id
head(datosm)
# Confirmamos que tenemos todas las condiciones en nuestro objeto "metadata"
all(colnames(datos) == rownames(datosm))
# Normalizar y transformar las cuentas de datos usando "default desing"
# Trabajaremos con una matriz
counts.matrix<- as.matrix(datos)
# Vamos a factorizar el contenido de la columna "strain" desde nuestro objeto "metadata"
datosm$strain = factor(datosm$strain)
datosm$experiment.number <- factor(datosm$experiment.number)
datosm$tissue <- factor(datosm$tissue)
############---Ahora vamos a crear un objeto DESeqDataSet que necesitaremos para crear nuestro PCA---#############
# Para crear este objeto necesitaremos la matriz de conteos de genes head(counts.matrix) y el conjunto de 
# metadatos head(datosm)
# La funcion DESeqDataSet es una clase de objeto utilizada por DESeq2 para almacenar datos de expresion genica junto con 
# la informacion sobre los experimnetos y el diseño experimental. Por tanto, el objeto resultante del siguiente codigo
# estara listo para ser analizado con DESeq2.
# countsData=rount(counts.matrix) es la matriz donde las filas son los genes y las columnas las muestras, mientras que los valores son los conteos 
# de lecturas (numero de veces que una secuencia especifica ha sido leida durante el secuenciamiento de RNA-seq)
# rount(counts.matrix) redondea los valores de los conteos de lecturas. 
# colData= datosm es un data frame que contiene la informacion sobre las muestras, como las condiciones experiemntales o los factores biologicos
# (tratamiento, control, etc) este objeto debe tener las mismas filas que el numero de muestras en counts.matrix
# design = ~ strain + experiment.number en este parametro especificamos los factores que se utilizaran para realizar la comparacion entre los datos, entre muestras
# En este caso, strain podría referirse a diferentes cepas o tipos de organismos, y experiment.number podría ser un número que identifica diferentes lotes de experimentos o grupos experimentales.
# El signo ~ es utilizado para indicar una fórmula en R que describe el diseño del experimento. Esto significa que los efectos de las variables strain y experiment.number se considerarán en el modelo estadístico para evaluar la expresión génica diferencial. En este caso, el modelo evaluará cómo varían los genes en función de estas dos variables.
# Define el diseño experimental: El diseño experimental es crucial para especificar cómo se van a analizar los efectos de las condiciones experimentales sobre la expresión génica. Esto permite a DESeq2 realizar un análisis de expresión diferencial adecuado.
# Con experiment.number: Captura la variabilidad técnica dentro de cada condición experimental (útil si las réplicas técnicas podrían aportar variabilidad importante).
#  strain refleje la condición experimental (control o tratamiento) y que experiment.number indique las réplicas técnicas dentro de cada condición.

# Crear el objeto DESeqDataSet considerando las réplicas técnicas
dds <- DESeqDataSetFromMatrix(
  countData = round(counts.matrix),  # Tu matriz de conteos (filas=genes, columnas=muestras)
  colData = datosm,  # El data frame con los metadatos
  design = ~ strain + experiment.number  # Considera tanto la condición como las réplicas técnicas
)

# Crear el objeto dds para DESeq2, usando la columna 'tissue' como parte del diseño
dds <- DESeqDataSetFromMatrix(
  countData = round(counts.matrix),  # Tu matriz de conteos (de genes x muestras)
  colData = datosm,  # El data frame con los metadatos
  design = ~ tissue + strain + experiment.number  # Fórmula de diseño (tissue, strain y experiment.number)
)

# Transformar los datos de cuentas a una forma útil para su visualización
# Los datos de cuentas se utilizan directamente para el analisis de expresion diferencial, mientras que los datos
# transformados se utilizan para la visualización.
# Enseguida necesitamos utilizar un algoritmo que nos permita transformar los datos de expresion para hacer
# que la variabilidad de los datos sea mas uniforma y adecuada para las visualizaciones, como el PCA o los heatmaps.

# Vst (Variance Stabilizing transformation) estabiliza la varianza a lo largo de los niveles de expresion. El algoritmo realiza una 
# transformación para hacer que los datos sigan una distribución más homogénea, de modo que las variaciones en los datos se distribuyan de manera más uniforme, independientemente de la magnitud de la expresión.
# VST es ideal para visualizaciones como PCA o heatmaps, especialmente cuando tienes datos con una amplia gama de valores de expresión (es decir, genes que tienen tanto bajas como altas expresiones). 
# Esto es porque estabiliza la varianza, lo que permite una mejor representación visual.

# rlog también es un algoritmo de transformación, pero está diseñado para hacer regularización de los datos de conteo. 
# A diferencia de VST, que estabiliza la varianza, rlog tiene el propósito de reducir el efecto de la dispersión aleatoria 
# (ruido) en los datos. Esto lo logra aplicando una regularización para suavizar los valores extremos de alta expresión.
# rlog es ideal cuando se desea reducir el ruido en las muestras con baja expresión. Es útil cuando los datos contienen muchas 
# cuentas cercanas a cero y se requiere una mayor estabilidad en las relaciones entre las muestras.

# VST para normalizar las cuentas. Es más rápido
vsd<- vst(dds, blind = FALSE) # blind FALSE para considerar el diseño experimental

# rlog para normalizar las cuentas
rld<- rlog(dds, blind = FALSE)

# Obtener la matriz de los valores transformados 
mat<- assay(rld)
mat_v<- assay(vsd)

# PCA usando la funcion plotPCA de DESeq2 
# ntop=500 indica que solo se deben considerar los 500 genes con mayor varianza en los datos de expresión para realizar el PCA
p<- plotPCA(rld, intgroup=c("strain"))
show(p)
p<- plotPCA(rld, intgroup=c("strain"), ntop=1000)
show(p)

p_v<- plotPCA(vsd, intgroup=c("strain"))
show(p_v)
p_v<- plotPCA(vsd, intgroup=c("strain"), ntop=1000)
show(p_v)

##############################################################################
######-- Ahora vamos a construir un plot PCA mejor usando ggplot2 --##########
library(ggplot2)
library(ggplotify)
library(ggrepel)
# Primero vamos a obtener los datos de PCA sin generar un gráfico
pcaData<- plotPCA(vsd, intgroup=c("sample.id","tissue","experiment.number"), returnData=TRUE)
# Cambiar ntop a 1000
pcaData <- plotPCA(vsd, intgroup = c("sample.id", "tissue", "experiment.number"), returnData = TRUE, ntop = 1000)

##### Calcular el porcentaje de varianza explicado por cada componente princial (PC)
# Este paso calcula qué porcentaje de la variabilidad total de los datos es explicado por cada uno de los componentes principales (PC1, PC2, etc.).
# attr(pcaData, "percentVar"): Este código obtiene la varianza explicada por cada componente principal (PC) a partir de los resultados del PCA que has calculado previamente. Los resultados del PCA, como las componentes principales (PC1, PC2, etc.), incluyen una medida de cuánta variabilidad en los datos es explicada por cada componente.
# 100 * attr(pcaData, "percentVar"): Multiplica por 100 para convertir esa proporción en un porcentaje.
# round(...): Redondea los valores para hacerlos más fáciles de interpretar (por ejemplo, redondear 24.678 a 25).
percentVar<- round(100 * attr(pcaData, "percentVar"))

#### Crear un gráfico básico con ggplot
# ggplot(pcaData, aes(PC1, PC2)): Crea un gráfico con ggplot2, utilizando los datos de pcaData. En este gráfico:
# PC1: Será el eje X (primer componente principal).
# PC2: Será el eje Y (segundo componente principal).
# Los puntos en este gráfico representan las muestras o las observaciones en tu análisis PCA.
g<- ggplot(pcaData,aes(PC1,PC2))
show(g)
g2<-g + geom_point(size=4, aes(col=tissue,fill=sample.id), shape = 21, stroke = 1)
show(g2)

# Cambiar el color de las muestras
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9, "YlOrRd")
brewer.pal(n = 12, name = "Paired")

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
g4<- g3 + geom_text_repel(aes(label=sample.id), direction = "both", nudge_y = 1.1, point.padding = 0.6,
                          box.padding = 0.5, min.segment.length = unit(0.2,"lines"), size=2.5)
g4 + theme_bw()

# Añadir un titulo y cambiar las etiquetas del eje x y y
g5<- g4 + theme_bw() + ggtitle("PCA con todas las librerías (lecturas crudas)") +
  xlab(paste0("PC1:",percentVar[1], "% variance")) +
  ylab(paste0("PC2:",percentVar[2], "% variance"))
show(g5)
# Si queremos separar el gráfico por tipo de tejido
g6<- g5 + facet_wrap(~ tissue)
g6


#######-- Construyendo un Heatmap con el paquete pheatmap --########
# Los heatmaps son utiles para explorar patrones generales del cambio de expression entre muestras.
# En este heatmap, el color rojo indica un alto nivel de expression y un color azul un bajo nivel de expresion
# En las filas de nuestra matriz corresponderan a los genes
# En las columnas de nuestra matriz corresponderan a las muestras

# Cargar los paquetes que necesitaremos para hacer nuestro heatmap
library(pheatmap)
library(dplyr)

# Primero, vamos a obtener una matriz desde el objeto transformado DESeq2
myMatrix<- SummarizedExperiment::assay(vsd)

# Obtener las ubicaciones de los genes más variables (n=500) y hacer un subconjunto de la matriz
topVarGenes<- head(order(rowVars(myMatrix, useNames = TRUE), decreasing = TRUE),500)
myMatrix_for_heatmap<- myMatrix[topVarGenes,]
# orde the columns 
column_order <- c("Control_1_Raiz", "Control_2_Raiz","Trat_1_Raiz", "Trat_2_Raiz",
                  "Control_1_Hojas", "Control_2_Hojas", "Trat_1_Hojas", "Trat_2_Hojas")
myMatrix_for_heatmap_order <- myMatrix_for_heatmap[, column_order]

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
phm1<- pheatmap(myMatrix_for_heatmap_order,
               #Escalar los datos por filas
               scale = "row",
               # Agrupar las filas y columnas 
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               # Mostrar/ocultar los genes o las muestras gene or sample names
               show_rownames = FALSE,
               show_colnames = TRUE,
               main = "Heatmap con 500 mas variables",
               fontsize = 13, 
               fontsize_row = 13,
               fontsize_col = 10,
               angle_col = "45") #90 grados
show(phm1)

# Como luciria un heatmap pero con 1000 genes?


###-- Add the annotation
# colours
datosm$Experiment<- c("Control_Raiz","Control_Raiz","Tratamiento_Raiz", "Tratamiento_Raiz", 
                      "Control_Hojas", "Control_Hojas", "Tratamiento_Hojas", "Tratamiento_Hojas")

display.brewer.pal(9, "Set1")
brewer.pal(n = 9, name = "Set1")

ann_colors <- list(
  strain = c("Control" = "#807DBA",
                     "Trat" = "#DE77AE"), 
  Experiment = c("Control_Hojas" = "darkgreen",
            "Control_Raiz" = "#377EB8",
            "Tratamiento_Raiz"="#FF7F00",
            "Tratamiento_Hojas"="#E41A1C")
)


phm2<- pheatmap(myMatrix_for_heatmap_order,
               #sacle the data by row
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







