# MR.GARCH.VaR_GOBIXDR

Código fuente que replica la investigación: ["Volatilidad realizada y el riesgo de mercado: Aplicación del Value-at-Risk (VaR) en el índice GOBIXDR"](https://sb.gob.do/publicaciones/publicaciones-tecnicas/volatilidad-realizada-y-el-riesgo-de-mercado-aplicacion-del-value-at-risk-var-en-el-indice-gobixdr/), Superintendencia de Bancos de la República Dominicana (2025).

**Palabras clave:** portafolio de inversiones, riesgos financieros, VaR, GARCH.
**Clasificación JEL:** G00, G10, G17, G20, G21.

## Resumen

En este repositorio se aplica el Value-at-Risk (VaR) utilizando la metodología híbrida GARCH bajo el procedimiento de Simulación Histórica Filtrada (FHS) y el método de Monte-Carlo. En el proceso de búsqueda del modelo óptimo, realizamos 32,400 calibraciones de distintas familias GARCH variando los parámetros de manera iterativa dentro de un amplio universo de configuraciones. El grupo de modelos que minimiza el criterio de selección es sometido a las pruebas y, posteriormente, a la fase de validación. A partir de las estimaciones realizadas en el año 2022, aplicamos criterios para formular un escenario de Stress-VaR y estimar el potencial consumo de riesgo del portafolio benchmark.

## Contenido del repositorio

| Archivo | Descripción |
|---|---|
| `VR_y_RM_Aplicacion_VaR_en_GOBIXDR.R` | Código fuente completo: funciones, análisis y resultados de la investigación (punto de entrada principal). |
| `GARCH_msperlin_functions.R` | Funciones auxiliares para la búsqueda/calibración de modelos GARCH. Adaptado de [Perlin et al. (2021), "A GARCH Tutorial with R"](https://github.com/msperlin/GARCH-RAC). |
| `Best_GARCH_Gobix1.csv` | Resultado precalculado de la calibración: criterio de selección de cada uno de los 32,400 modelos probados. |
| `Top1PCT_model_validation_test.csv` | Resultado precalculado de la prueba/validación del top 1% de modelos/configuraciones. |
| `paper_VaR_GOBIXDR_final.pdf` | Versión final del documento publicado. |
| `refactor/` | Versión modular y limpia del script (`functions.R` + `main.R`), con la calibración costosa desactivada por defecto (`RUN_FULL_CALIBRATION <- FALSE`). Ver `refactor/README.md`. |

## Cómo ejecutar

1. **Requisitos:** R ≥ 4.0 y conexión a internet — el script descarga la serie del GOBIXDR desde la BVRD (`https://www.bvrd.com.do/indice/Data/GobixDataIRP.csv`) y los benchmarks internacionales (EMB, BND, IEI, TLT, HYG, AGG, IEF, entre otros) desde Yahoo Finance vía `quantmod`.
2. **Instalar los paquetes requeridos:**

```r
install.packages(c("xts", "zoo", "lubridate", "rugarch", "quantmod",
                   "PerformanceAnalytics", "roll", "pastecs", "stargazer",
                   "knitr", "fAssets", "Hmisc", "kableExtra", "stringr"))
```

3. **Fijar el directorio de trabajo en la raíz del repositorio** (el script hace `source("GARCH_msperlin_functions.R")` y lee los `.csv` por ruta relativa):

```r
setwd("ruta/al/repositorio/MR.GARCH.VaR_GOBIXDR")
source("VR_y_RM_Aplicacion_VaR_en_GOBIXDR.R")
```

También puede ejecutarse por secciones desde RStudio (el script está organizado con encabezados `# Sección ----`).

**Advertencia sobre tiempos de ejecución:** la calibración completa de los 32,400 modelos GARCH (sección *GARCH calibration* del script) toma varias horas. Por esa razón, sus resultados se distribuyen precalculados en `Best_GARCH_Gobix1.csv` y `Top1PCT_model_validation_test.csv`, que el script lee directamente — el análisis puede reproducirse sin re-calibrar.

## Principales funciones

* **`get.gobix()`**: realiza la descarga de la serie de tiempo del GOBIXDR.
* **`garch.ts.forecast()`**: realiza la estimación del cono de probabilidades en el nivel de precios aplicando la simulación GARCH diaria.
* **`garch.oos.sim()`**: realiza la prueba y validación del modelo (out-of-sample).

## Licencia

Este proyecto se distribuye bajo la licencia [GNU GPL v3](LICENSE).
