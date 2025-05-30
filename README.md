# Tesis de Licenciatura (Lucía Castelli, 2025)

Este repositorio contiene los archivos finales (`final_outputs`) y los scripts utilizados en el análisis de expresión genética diferencial y otros análisis transcriptómicos realizados para mi tesis de licenciatura en Biotecnología.

## 🧪 El proyecto

Identificación de firmas glicano-relacionadas en colitis ulcerosa y cáncer colorrectal asociado mediante análisis transcriptómicos y posterior validación in vivo.

## 📂 Estructura del Repositorio

├── final_outputs/ # Tablas finales del análisis (ej. genes diferencialmente expresados)<br>
├── scripts/ # Scripts en R utilizados en los análisis<br>
├── README.md # Este archivo



## 📊 Contenido

- final_outputs/: Contiene tablas procesadas listas para interpretación o presentación, como:
  - Listas de genes diferencialmente expresados (DEGs)
  - Resultados de análisis de enriquecimiento funcional (GO, KEGG)

- scripts/: Scripts escritos en **R** para reproducir los análisis, incluyendo:
  - Preprocesamiento y normalización de datos
  - Análisis de expresión diferencial con `limma`
  - Análisis de enriquecimiento funcional
  - Generación de gráficos (volcano plots, heatmaps, etc.)

## 🧬 Datos

Los análisis se realizaron utilizando datos públicos de expresión génica con los datasets de GEO explicitados en la tesis entregada y al principio de cada script.

## 💻 Requisitos

Los scripts fueron desarrollados en R (versión ≥ 4.2.0) y requieren al menos el uso de las siguientes librerías principales:

- `limma`
- `edgeR`
- `DESeq2`
- `biomaRt`
- `clusterProfiler`
- `ggplot2`
- `tidyverse`

🧾 Licencia

Este repositorio está bajo la Licencia MIT. Podés usar, modificar y compartir el contenido, citando adecuadamente.

👩‍🔬 Autoría

Lucía Castelli
Estudiante de Biotecnología
Tesis de Licenciatura – UADE
Año: 2025

📬 Contacto
Si tenés dudas o sugerencias, podés contactarme a través de GitHub o por email a [lucicastelli@uade.edu.ar].
