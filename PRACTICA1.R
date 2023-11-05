#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 1 OCTUBRE 23:59
## Se requiere la entrega con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #esto mide las dimensiones
head(data) #para ver los datos del encabezado
tail(data) #para ver los datos finales

# Hacemos un primer histograma para explorar los datos
hist(data)

# Transformamos los datos con un logaritmo 
log2(data)
data_log=log2(data)
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
hist(data_log)tiene forma de campana de Gauss, pasa de una distribución sesgada a una distribución normal, para visualizar mejor los datos

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
boxplot(data_log)
boxplot(data_log, col=c("blue","blue","blue","green","green","green"), main="GSE5583 - boxplots", las=2)
#main es para cambiar el título, col para el color y las2 para poner los nombres d elos ejes en vertical

# ¿Qué es un boxplot?
#Es un tipo de gráfico que muestra un resumen de una gran cantidad de datos en medidas descriptivas, además de intuir su morfología y simetría.
#Este tipo de gráficos nos permite identificar valores atípicos y comparar distribuciones.

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
hc=hclust(as.dist(1-cor(data_log)))
plot(hc, main="clustering")
#sí, es correcta la separación
#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
wt <-data[,1:3]
ko <-data[,4:6]
class(wt)
#hemos generado una matriz
#la coma es para escoger solo las columnas, los wt van de la columna 1 a la 3 y los ko de la 4 a la 6
head(wt) #para acceder mejor a los datos de arriba, visualizar las primeras líneas

# Calcula las medias de las muestras para cada condición. Usa apply
media_wt = apply(wt, 1, mean) #calcula la media de cada fila de la tabla (1 filas, 2 columnas)
media_ko = apply(ko, 1, mean)

# ¿Cuál es la media más alta?
max(media_wt)
max(media_ko)

# Ahora hacemos un scatter plot (gráfico de dispersión)
#es un gráfico de dispersión, compara las medias del eje x con las del eje y
plot(media_ko ~ media_wt, xlab = "WT", ylab = "KO", main = "dispersión")

# Añadir una línea diagonal con abline
abline(0, 1, col = "red")
abline(h=2,col="blue")
abline(v=5,col="green")
#abline hay que ejecutarlo con el plot abierto, h es horizontal y v vertical, hay que indicarlo porque si sólo pongo un número no sabe a qué me refiero

# ¿Eres capaz de añadirle un grid?
grid() #para añadir cuadrícula

# Calculamos la diferencia entre las medias de las condiciones
dif.media = media_wt - media_ko

# Hacemos un histograma de las diferencias de medias
hist(dif.media)

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values, si sale <0.05 es significativo
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué? 
	#Porque estaría manipulando los datos y no serían fiables.
# ¿Cuántas valores tiene cada condición? Tenemos 2 condiciones y 6 muestras, 3 muestras para cada condición (réplicas)
	#hacemos un t.test para cada fila, hay 12488 genes, tenemos un pvalue para cada gen
	#los primeros 6 genes head(pvalue) no tienen pvalues significativos

pvalue=NULL
tstat=NULL
#bucle
for(i in 1 : nrow(data)) { #para cada gen
	x = wt[i,] #gene wt numero i
	y = ko[i,] #gene ko numero i
	#hacemos en test
	t = t.test(x, y)
	#añadimos p-value a la lista
	pvalue[i] = t$p.value
	#añadimos las estadisticas a la lista
	tstat[i] = t$statistic
}

head(pvalue)
length(pvalue)
# Ahora comprobamos que hemos hecho TODOS los cálculos

#la transformación sólo se emplea para las gráficas

# Hacemos un histograma de los p-values.
hist(pvalue) #tenemos que transformar los datos porque este histograma es muy grande y no es representativo
# ¿Qué pasa si le ponemos con una transformación de -log10? Cambiamos la escala de proporción y la distribución para que sea más visual
hist(-log10(pvalue),col = "blue")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(dif.media, -log10(pvalue), main = "GSE5583-Volcano")
#mirar genes sobreexpresados (rojo) y deprimidos (azul)

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
dif.media_cutoff = 2
pvalue_cutoff = 0.01
abline(v = dif.media_cutoff, col = "blue", lwd = 3)
#abline(v = -dif.media_cutoff, col ="red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)



# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_dif.media = abs(dif.media) >= dif.media_cutoff
dim(data[filter_by_dif.media, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined = filter_by_dif.media & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(dif.media, -log10(pvalue), main= "volcano 2")
points(dif.media[filter_combined], -log10(pvalue[filter_combined]),col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés? Por como hemos calculado la diferencia de medias wt-ko, por lo que los ko dan valores mayores, los sobreexpresados están en negativo y por eso están en la gráfica volcano a la izquierda en rojo
plot(dif.media, -log10(pvalue), main = "volcano 3")
points(dif.media[filter_combined & dif.media < 0], -log10(pvalue[filter_combined & dif.media < 0]),col = "red")
points(dif.media[filter_combined & dif.media > 0], -log10(pvalue[filter_combined & dif.media > 0]),col = "blue")
#dif.media = wt.media - ko.media

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7, labRow=FALSE)

# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap? 
#hemos hecho agrupación por ko y wt y por expresión

# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap


# Lo guardamos en un archivo PDF


# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered,"GSE5583_DE.txt", sep= "\t", quote = FALSE)
