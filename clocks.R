# -------------------------------------------------------
# Estimación de edad epigenética con methylclock
# -------------------------------------------------------

# Instalar una sola vez, si es necesario:
# BiocManager::install("methylclock")

library(methylclock)
library(ggplot2)

# -------------------------------------------------------
# 1. Preparar matriz de metilación
# -------------------------------------------------------
# Usar beta original normalizada con NOOB
# NO usar beta_clean ni solo las 7 CpG significativas

# Comprobar cobertura de los relojes
clock_qc <- checkClocks(beta)

# -------------------------------------------------------
# 2. Calcular relojes epigenéticos
# -------------------------------------------------------

clock_results <- DNAmAge( beta,
                          clocks = c("Horvath", "Hannum", "Levine"),
                          age = samples$age,
                          cell.count = FALSE,
                          fastImp = TRUE,
                          normalize = FALSE )

clock_results <- as.data.frame(clock_results)

head(clock_results)
colnames(clock_results)

clock_results$labels <- samples$labels

# -------------------------------------------------------
# 6. Exportar resultados y graficar
# -------------------------------------------------------
write.csv2( clock_results, "results/relojes_epigeneticos_PD.csv")

#plots
plotDNAmAge(clock_results$Horvath,
            clock_results$age)

boxplot(clock_results$ageAcc2.Levine~clock_results$labels)

boxplot(clock_results$ageAcc.Horvath~clock_results$labels)

boxplot(clock_results$ageAcc2.Hannum~clock_results$labels)

plotDNAmAge(clock_results$Horvath, clock_results$age)

#analisis estadistico no parametrico
kruskal.test(clock_results$ageAcc2.Levine~clock_results$labels)
