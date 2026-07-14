# Versión modular del análisis VaR GOBIXDR

Esta carpeta contiene una versión modular y depurada del script original
`VR_y_RM_Aplicacion_VaR_en_GOBIXDR.R`:

- **`functions.R`** — todas las funciones auxiliares, extraídas textualmente del original.
  Nota: `plot.quadrants` y `rolling.var.backtest` se redefinen intencionalmente a mitad
  de `main.R` (igual que en el original); `functions.R` contiene la primera versión de cada una.
- **`main.R`** — el flujo de análisis completo, en el mismo orden que el original,
  sin los gráficos no utilizados ni los bloques de código muerto comentados.

**El script original sigue siendo la referencia canónica para replicar el paper.**

## Ejecución

Ejecutar `main.R` con el directorio de trabajo en la **raíz del repositorio** (no en
`refactor/`), para que se encuentren `GARCH_msperlin_functions.R`, `Best_GARCH_Gobix1.csv`
y `Top1PCT_model_validation_test.csv`:

```r
# desde la raíz del repositorio
source("refactor/main.R")
```

## Bandera `RUN_FULL_CALIBRATION`

En el bloque `# Parámetros ----` de `main.R`, la bandera `RUN_FULL_CALIBRATION`
(por defecto `FALSE`) controla los dos bloques de cómputo pesado:

- la calibración completa de 32,400 modelos GARCH (`find_best_arch_model`), y
- el backtesting/validación rolling del Top 1% de modelos.

Con `FALSE` (recomendado) se usan los resultados precomputados de
`Best_GARCH_Gobix1.csv` y `Top1PCT_model_validation_test.csv`. Con `TRUE`
la calibración se re-ejecuta desde cero (tarda horas).

## Corrección de `source()`

El original hace `source("GARCH_msperlin_functions.r")` (extensión `.r` minúscula),
pero el archivo real es `.R`. En macOS/Windows funciona igual, pero en Linux el
sistema de archivos distingue mayúsculas; `main.R` usa la extensión correcta `.R`.
