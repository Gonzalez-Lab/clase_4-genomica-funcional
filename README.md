Análisis de metilación del ADN: EWAS y relojes epigenéticos

Curso de Epigenética y Genómica Funcional | ITBA

Este repositorio contiene el material práctico para el análisis bioinformático de datos de metilación del ADN obtenidos mediante microarrays Illumina HumanMethylation450 BeadChip (450K).

A partir de un conjunto de muestras de sangre periférica de pacientes con enfermedad de Parkinson y controles, se desarrollan dos análisis complementarios:

Estudio de asociación epigenómica (EWAS): identificación de sitios CpG diferencialmente metilados entre pacientes y controles.

Relojes epigenéticos: estimación de la edad epigenética y evaluación de su aceleración mediante los modelos de Horvath, Hannum y Levine PhenoAge.

Objetivos de aprendizaje

Comprender la estructura de los datos de metilación obtenidos mediante microarrays.

Importar y procesar archivos IDAT utilizando minfi.

Realizar controles de calidad y normalización NOOB.

Interpretar valores beta y valores M de metilación.

Estimar la composición celular de muestras de sangre periférica.

Identificar sitios CpG asociados con una condición clínica mediante modelos lineales.

Anotar las CpG significativas según su ubicación y contexto genómico.

Visualizar resultados mediante PCA, volcano plots y heatmaps.

Aplicar relojes epigenéticos y comparar las estimaciones de aceleración entre grupos.

Dataset

GSE111629 — Gene Expression Omnibus (GEO)

Se tomo un subconjunto de datos correspondientes a hombres caucasicos
Datos de metilación del ADN de muestras de sangre periférica de pacientes con enfermedad de Parkinson y controles.

Plataforma: Illumina HumanMethylation450 BeadChip.

Tipo de datos: intensidades de metilación en archivos IDAT.

Organismo: Homo sapiens.

Grupos de estudio: Parkinson (PD) y controles.

Los datos originales están disponibles en:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111629

Organización del análisis

Script 1. Análisis EWAS

Procesamiento de datos y análisis diferencial de metilación

El flujo de trabajo comprende:

Lectura de archivos IDAT y metadatos.

Normalización mediante NOOB.

Control de calidad y estimación del sexo mediante metilación.

Exploración de las distribuciones de valores beta.

Estimación de proporciones celulares con meffil.

Filtrado de sondas según sus valores p de detección.

Construcción de la matriz de diseño estadístico.

Análisis diferencial mediante limma, ajustando por edad y composición celular.

Anotación genómica de las CpG significativas.

Visualización mediante PCA, volcano plot y heatmap.

Se utiliza un umbral de significación de FDR < 0,05, aplicando la corrección de Benjamini-Hochberg.

El análisis diferencial se realiza sobre valores M, mientras que los valores beta se utilizan para facilitar la interpretación biológica de los niveles de metilación.

Script 2. Relojes epigenéticos

Estimación de edad epigenética y aceleración del envejecimiento

Se utiliza el paquete methylclock para calcular tres relojes:

Horvath: reloj multitejido basado en metilación del ADN.

Hannum: reloj desarrollado a partir de datos de sangre.

Levine PhenoAge: estimador epigenético relacionado con un fenotipo de envejecimiento.

El análisis incluye:

Evaluación de cobertura de las sondas requeridas por los relojes.

Estimación de edad epigenética a partir de valores beta normalizados.

Comparación entre edad cronológica y edad epigenética.

Visualización de las métricas de aceleración epigenética.

Comparación no paramétrica entre pacientes con Parkinson y controles.
