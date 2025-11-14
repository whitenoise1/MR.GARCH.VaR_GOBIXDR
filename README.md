# MR.GARCH.VaR_GOBIXDR
Código fuente que replica la investigacion: ["Volatilidad realizada y el riesgo de mercado: Aplicación del
Value-at-Risk (VaR) en el índice GOBIXDR"](https://sb.gob.do/publicaciones/publicaciones-tecnicas/volatilidad-realizada-y-el-riesgo-de-mercado-aplicacion-del-value-at-risk-var-en-el-indice-gobixdr/), (2025).

Palabras clave: portafolio de inversiones, riesgos financieros, VaR, GARCH.
Clasificación JEL: G00, G10, G17, G20, G21.

En este repositorio, se aplica el Value-at-Risk (VaR) utilizando la metodología híbrida GARCH bajo el procedimiento de Simulación Histórica Filtrada (FHS) y el método de Monte-Carlo. En el proceso de búsqueda del modelo óptimo, realizamos 32,400 calibraciones de distintas familias GARCH variando los parámetros de manera iterativa dentro de un amplio universo de configuraciones. El grupo de modelos que minimiza el criterio de selección es sometido a las pruebas y, posteriormente, a la fase de validación. A partir de las estimaciones realizadas en el año 2022, aplicamos criterios para formular un escenario de Stress-VaR y estimar el potencial consumo de riesgo del portafolio benchmark. 

* Best_GARCH_Gobix1.csv: Es resumen del resultado presenta el criterio de calibración de cada modelo probado.
* Top1PCT_model_validation_test.csv: Presenta el top 1% de modelos/configuraciones.
* GARCH_msperlin_functions.R: Función auxiliar que corre la busqueda de modelo.
* VR_y_RM_Aplicacion_VaR_en_GOBIXDR.R: El archivo contiene las funciones, el análisis y los resultados completos de la investigación.

# Principales funciones
* **get.gobix()**: Realiza la descarga de la serie de tiempo del GOBIXDR.
* **garch.ts.forecast()**: Realiza la estimación del cono de probabilidades de movimiento de precios.
* **garch.oos.sim()**: Realiza la prueba y validación del modelo.

# Librerias requeridas
```{r}
library(xts)
library(zoo)
library(lubridate)
library(rugarch)
library(quantmod)
library(PerformanceAnalytics)
library(roll)
library(pastecs)
library(stargazer)
library(knitr)
library(fAssets)
library(Hmisc)
library(kableExtra)
library(stringr)
```
