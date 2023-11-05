#Trabajo R por Ana Caro
#1. Carga los datos y examínalos
ruta_completa <- "C:/Users/Usuario/Desktop/Trabajo R/datos_trabajor.txt"
data <- read.table(ruta_completa, header = TRUE, sep = "\t")
head(data) #para ver los datos del encabezado
dim(data) #para medir las dimensiones
str(data) #para obtener información de la estructura: nombres de las columnas, tipos de datos de cada columna y algunas de las filas iniciales del data frame
summary(data) #para obtener un resumen estadístico
#¿Cuántas variables hay? ¿Cuántos tratamientos?
#Hay 2 variables y 5 tratamientos


#2. Haz un boxplot para nuestros datos. Uno para cada variable. Colorea Variable 1 y Variable 2 de forma diferente.
#primero hacemos una transformación logarítmica de los datos para visualizarlos mejor
log2(data)
data_log=log2(data)
boxplot(data_log) #boxplot total de los datos
#para crear dos boxplot de un data frame con dos variables: 
boxplot(Variable1 ~ Tratamiento, data = data, col='blue', main = 'Variable 1')
boxplot(Variable2 ~ Tratamiento, data = data, col='green', main = 'Variable 2')
#variable1 en azul y variable2 en verde, col para seleccionar color y main para el título


#3. Haz un gráfico de dispersión con las dos variables.
plot(data$Variable1, data$Variable2,
	col = data$Tratamiento,
	xlab = "Variable 1",
	ylab = "Variable 2",
	main = "Gráfico de Dispersión de Variable1 vs. Variable2")
	
#4. Ponle leyenda al gráfico del apartado anterior. En el margen inferior derecho.
legend('bottomright', 
	legend = unique(data$Tratamiento), 
	title = 'tratamientos', 
	fill = unique(data$Tratamiento)
)


#5. Haz un histograma para cada variable.
hist(datos$Variable1, main = "Histograma de Variable 1", xlab = "Valores de Variable 1", col = "blue")
hist(datos$Variable2, main = "Histograma de Variable 2", xlab = "Valores de Variable 2", col = "green")


#6. Haz un factor en la columna tratamiento y guárdalo en una variable.
data$Tratamiento_Factor <- factor(data$Tratamiento)
data$Tratamiento_Factor

#7. Calcula la media y la desviación estándar para cada tratamiento.
resultados_media <- aggregate(cbind(Variable1, Variable2) ~ Tratamiento, data = datos, FUN = mean)
resultados_desviacion <- aggregate(cbind(Variable1, Variable2) ~ Tratamiento, data = datos, FUN = sd)
resultados_media
resultados_desviacion
#FUN sirve para dar flexibilidad para especificar la operación que quieres aplicar a tus datos
#cbind para visualizar las columnas juntas
#mean para la media y sd para la desviación estándar

# Para renombrar las columnas:
colnames(resultados_media) <- c("Tratamiento", "Media_Variable1", "Media_Variable2")
colnames(resultados_desviacion) <- c("Tratamiento", "Desviacion_Variable1", "Desviacion_Variable2")

# Para combinar los resultados en un único data frame utilizamos merge:
resultados <- merge(resultados_media, resultados_desviacion, by = "Tratamiento")
resultados #también se pueden imprimir los resultados usando print()


#8. Averigua cuántos elementos tiene cada tratamiento.
#table nos sirve para el conteo de frecuencias
conteo_tratamientos <- table(data$Tratamiento)
conteo_tratamientos


#9. Extrae los datos para el tratamiento 1 y el tratamiento 4 y guárdalos cada uno en una variable diferente.
tratamiento_1 <- data[data$Tratamiento == 1, ]
tratamiento_1
tratamiento_4 <- data[data$Tratamiento == 4, ]
tratamiento_4


#10. Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la variable 1 son iguales. 
	#¿Puedes comprobarlo? Para ello, necesitarás comprobar
	#primero si los datos se distribuyen de forma normal. En función del resultado de la
	#prueba de normalidad, ¿qué test usarías? ** En general, asumimos que las muestras
	#son independientes, pero ¿son sus varianzas iguales? Actúa de acuerdo a tus resultados
#Shapiro.test se emplea para determinar si la distribución de una muestra es normal o no.
#La hipótesis nula para esta prueba es que los datos siguen una distribución normal.

# Prueba de normalidad para tratamiento 1
shapiro.test(tratamiento_1$Variable1)
# Prueba de normalidad para tratamiento 4
shapiro.test(tratamiento_4$Variable1)
#Tras realizar saphiro test no se puede descartar la hipótesis nula dado que los valores p en ambos casos no son significativamente bajos y están por encima del nivel de significancia común (alfa = 0.05).

#procedo a comprobar si las varianzas son iguales
# Prueba de igualdad de varianzas 
var.test(tratamiento_1$Variable1, tratamiento_4$Variable1)
#el valor p es muy bajo, lo que sugiere que hay evidencia significativa en contra de la hipótesis nula de igualdad de varianzas entre Tratamiento 1 y Tratamiento 4 para la Variable 1. 
#las varianzas de los dos grupos son estadísticamente diferentes

#en este caso usamos un test no paramétrico llamado Wilcoxon
wilcox.test(tratamiento_1$Variable1, tratamiento_4$Variable1)
#Hay evidencia estadística para afirmar que las medianas de Tratamiento 1 y Tratamiento 4 para la Variable 1 son diferentes.
#Esto respalda la hipótesis alternativa de que existe una diferencia significativa en la ubicación de las distribuciones entre los dos tratamientos.





