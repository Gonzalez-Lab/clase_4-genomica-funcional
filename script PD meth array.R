# BiocManager::install(c("minfi",
#                       "IlluminaHumanMethylation450kmanifest",
#                       "IlluminaHumanMethylation450kanno.ilmn12.hg19"))

library(minfi)
library(meffil)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

# -------------------------------------------------------
# 1. Lectura de muestras y archivos IDAT
# -------------------------------------------------------
path <- "data/IDATs"
list.files(path)

samples <- read.csv2(file.path(path, "sample_sheet.csv"))

#editamos la columna Basename para transformarla en el path a los archivos
samples$Basename <- file.path(path, samples$Basename)

#FUNCION DE minfi PARA LEER LOS DATOS DEL MICROARRAY read.metharray
#usamos la columna Basename con los path a los IDATs como argumento 
rgset <- read.metharray(samples$Basename, verbose = TRUE)

# Chequear que el set de muestras y las samples coincidan
ncol(rgset) == nrow(samples)

# Asignar identificadores GSM al set de muestras
colnames(rgset) <- samples$id

# Incorporar metadatos
rownames(samples) <- samples$id
pData(rgset) <- S4Vectors::DataFrame(samples)

# -------------------------------------------------------
# 2. Procesamiento de metilación
# -------------------------------------------------------

# FUNCION DE minfi PARA NORMALIZACION NOOB
mset <- preprocessNoob(rgset) #generamos un objeto GenomicMethylSet 

# Asociar coordenadas genómicas
mset <- mapToGenome(mset)

head(granges(mset))

#obtener matrices de valores beta y M

beta <- getBeta(mset)
M <- getM(mset)

dim(beta)
dim(M)

# Verificar el orden de las muestras
ncol(beta) == nrow(samples)
identical(colnames(beta), samples$id)

# -------------------------------------------------------
# 3. Control de calidad y sexo estimado
# -------------------------------------------------------
# estimacion de sexo por metilacion de sondas en chry
sex_pred <- getSex(mset)

colData(mset)$predictedSex <- sex_pred$predictedSex
colData(mset)$xMed <- sex_pred$xMed
colData(mset)$yMed <- sex_pred$yMed

plotSex(mset)

#para evaluar la senial met y unmet y detectar muetsras problematicas
plot(as.matrix(getQC(mset)))

# -------------------------------------------------------
# 4. Grupos clínicos
# -------------------------------------------------------

# creamos las etiquetas
samples$labels <- ifelse( samples$state == "PD-free control", "control","PD" )

# convertimos etiquetas a factor
samples$labels <- factor(samples$labels, levels = c("control", "PD"))

# vemos N de cada grupo
table(samples$labels)

#creamos indices para cada grupo
normalIndex <- which(samples$labels == "control")
PDIndex <- which(samples$labels == "PD")

# -------------------------------------------------------
# 5. Distribución de valores beta
# -------------------------------------------------------

# grafico de densidad de la muestra control 1
plot( density( beta[, normalIndex[1]] , na.rm = TRUE, from = 0, to = 1),
  main = "Distribución de metilación",
  xlab = "Valor beta",
  ylab = "Densidad",
  xlim = c(0, 1),
  col = "black",
  lwd = 2
)

#agregamos al grafico todas las muestras control
for (i in normalIndex) { 
    lines( density( beta[, i], na.rm = TRUE, from = 0, to = 1 ),
    col = adjustcolor("black", alpha.f = 0.3) )  }

#agregamos al grafico todas muestras PD
for (i in PDIndex) {
    lines( density( beta[, i], na.rm = TRUE, from = 0, to = 1),
    col = adjustcolor("red", alpha.f = 0.3)  ) }

legend( "top", legend = c("Control", "Parkinson"), col = c("black", "red"),
  lwd = 2, bty = "n" )

# -------------------------------------------------------
# 6. Estimación de composición celular con meffil
# -------------------------------------------------------

cell_counts_matrix <- meffil.estimate.cell.counts.from.betas(
                              beta, "blood gse35069 complete" )

cell_counts <- as.data.frame(cell_counts_matrix)

head(cell_counts)

identical(rownames(cell_counts), colnames(beta))

#Unimos la estimacion de poblaciones celulares a samples
samples <- cbind(samples, cell_counts)

# -------------------------------------------------------
# 7. Filtrado de sondas por detection p-value
# -------------------------------------------------------

# Calcular detection p-values a partir de los IDAT (TARDA)
detP <- detectionP(rgset)

# Reordenar las sondas de detP
detP <- detP[rownames(M), ]

# Verificar que ahora coincidan
identical(rownames(detP), rownames(M))
identical(rownames(M), rownames(beta))

# Calcular proporción de fallos de detección
failed_probes <- rowMeans( is.na(detP) | detP >= 0.01 )

# Conservar sondas detectadas en al menos el 95% de las muestras
keep_detection <- failed_probes <= 0.05

# Conservar sondas con valores M finitos
keep_finite <- rowSums(is.finite(M)) == ncol(M)

# Aplicar ambos filtros
keep <- keep_detection & keep_finite

# limpiamos las matrices beta y M
M_clean <- M[keep, , drop = FALSE]
beta_clean <- beta[keep, , drop = FALSE]

# Resumen
cat("Sondas originales:", nrow(M), "\n")
cat("Sondas después del filtrado:", nrow(M_clean), "\n")
cat("Sondas eliminadas:", sum(!keep), "\n")

dim(M_clean)
dim(beta_clean)

# -------------------------------------------------------
# 8. Diseño estadístico
#---------------------------------------------------------

#convertimos explicitamente las variables a tipo numerico
samples$age <- as.numeric(samples$age)

cell_vars <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Eos", "Neu")

samples[cell_vars] <- lapply(samples[cell_vars],as.numeric)

#definimos las variables del modelo
model_vars <- c("labels", "age", cell_vars)

# Modelo ajustado por edad y composición celular
#Se omite Neu para evitar redundancia entre las fracciones celulares

design <- model.matrix(  ~ labels + age + CD8T + CD4T + NK + Bcell + Mono + Eos,
                          data = samples )

colnames(design)

# -------------------------------------------------------
# 9. Análisis diferencial con limma
# -------------------------------------------------------

# Modelo ajustado por edad y composición celular
fit <- lmFit(M_clean, design)
fit <- eBayes(fit)

# Resultados Parkinson vs controles
results_PD <- topTable( fit, coef = "labelsPD", number = Inf, 
                                                adjust.method = "BH" )

head(results_PD)

# Número de CpG significativas
table(results_PD$adj.P.Val < 0.05)

# Seleccionar CpG significativas
sig_PD <- subset( results_PD, adj.P.Val < 0.05 )

sig_PD

#-------------------------------------------------------------------------
# 10. obtener anotacion de las cg diferenciales
#-------------------------------------------------------------------------

# Obtener anotación del array
anno <- as.data.frame(getAnnotation(mset))

# Ver qué columnas contiene
colnames(anno)

# Seleccionar las CpG significativas
cpg_sig <- rownames(sig_PD)

# Extraer anotación
anno_sig <- anno[cpg_sig, , drop = FALSE]

# Armar tabla con resultados estadísticos y anotación
tabla_CpG <- data.frame(CpG = cpg_sig,
                        sig_PD[cpg_sig, c("logFC", "P.Value", "adj.P.Val")],
                        chr = as.character(anno_sig$chr),
                        pos = anno_sig$pos,
                        Gene = anno_sig$UCSC_RefGene_Name,
                        Gene_Region = anno_sig$UCSC_RefGene_Group,
                        CpG_Island = anno_sig$Relation_to_Island,
                        row.names = NULL)

# Dirección de la diferencia en escala M
tabla_CpG$Direccion <- ifelse( tabla_CpG$logFC > 0, "hipermetilado en PD",
                                                    "hipometilado en PD"  )

print(tabla_CpG)

#guardar resultados finales
write.csv2( tabla_CpG,"results/CpG_significativas_anotadas.csv")

# grafico de cg diferenciales
ggplot(tabla_CpG, aes(x = reorder(CpG, logFC), y = logFC, fill = CpG_Island )) +
       geom_col() + coord_flip() +
  labs(title = "CpG asociadas con Parkinson",
       subtitle = "EWAS ajustado por edad y composición celular",
       x = "Sonda CpG",
       y = "Diferencia de metilación (escala M)",
       fill = "Contexto CpG" ) +
  theme_classic(base_size = 12)


# -------------------------------------------------------
# 11. PCA de metilación usando valores M normalizados
# -------------------------------------------------------

# Filtrar sondas con varianza cero
var_M <- apply(M_clean, 1, var)

M_pca <- M_clean[var_M > 0, ]

# Seleccionar las 500 sondas más variables para facilitar el cálculo del PCA
var_M <- var_M[var_M > 0]

top_var <- names(sort(var_M, decreasing = TRUE))[1:500]

M_pca <- M_pca[top_var, ]

# PCA: muestras en filas, sondas en columnas
pca <- prcomp(t(M_pca), center = TRUE, scale. = FALSE)

# Porcentaje de varianza explicada
var_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Gráfico
plot( pca$x[, 1], pca$x[, 2],
  col = as.numeric(samples$labels), pch = 19,
  xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
  ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
  main = "PCA - Metilación del ADN" )

legend( "topright", legend = levels(samples$labels),col = 1:length(levels(samples$labels)),
  pch = 19)

# -------------------------------------------------------
# 12. Volcano plot
# -------------------------------------------------------

EnhancedVolcano(results_PD,
                lab = rownames(results_PD),
                x = "logFC",
                y = "adj.P.Val",
                pCutoff = 0.05,
                FCcutoff = 0, # no queremos aplicar un umbral mínimo de tamaño de efecto
                title = "Parkinson vs Control - EWAS",
                subtitle = "Ajustado por edad y composición celular",
                xlab = "Diferencia de metilación (escala M)",
                ylab = "-log10(FDR)" )


# -------------------------------------------------------
# 13. Heatmap de CpG diferencialmente metiladas
# -------------------------------------------------------

# -------------------------------------------------------
# 13. Heatmap de CpG diferencialmente metiladas
# -------------------------------------------------------

mat.beta <- beta_clean[rownames(sig_PD),, drop = FALSE]

# Z-score por CpG
mat.beta <- t(scale(t(mat.beta)))

# Ordenar muestras por grupo
orden <- order(samples$labels)

mat.beta <- mat.beta[, orden, drop = FALSE]

# Anotación de grupos
annotation_col <- data.frame(
  Grupo = samples$labels[orden],
  row.names = colnames(mat.beta)
)

annotation_colors <- list(Grupo = c(control = "#77A6B6",
                                         PD = "#D98989"   ))

pheatmap(
  mat.beta,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_col = sum(samples$labels == "control"),
  fontsize_row = 10,
  main = "CpG diferencialmente metiladas - Parkinson"
)
