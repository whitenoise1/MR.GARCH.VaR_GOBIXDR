# title: "RIESGO DE MERCADO Y VOLATILIDAD REALIZADA: Aplicacion del Value-at-Risk (VaR) en el indice GOBIXDR"
# author: Stefan Bolta, FRM.
# affiliation: Superintendencia de Bancos
# creado: 17/09/2023
# modificado: 18/07/2024
# subtitle: Superintendencia de Bancos, Departamento de Estudios Económicos
#
# Versión modular de VR_y_RM_Aplicacion_VaR_en_GOBIXDR.R — ejecutar con el
# directorio de trabajo en la raíz del repositorio.

# Carga las librerias utilizadas ----
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

# Parámetros ----
# TRUE re-ejecuta la calibración completa (32,400 fits + backtesting del Top 1%,
# tarda horas); FALSE usa los resultados precomputados en Best_GARCH_Gobix1.csv
# y Top1PCT_model_validation_test.csv.
RUN_FULL_CALIBRATION <- FALSE

fecha.inicio.datos <- '2014-01-01'   # inicio muestra (descargas Yahoo Finance)
fecha.fin.datos    <- "2022-12-31"   # fin muestra (descargas Yahoo Finance)
fecha.inicio.mc    <- "2019-01-01"   # ventana forecast Monte-Carlo (garch.ts.forecast)
fecha.fin.mc       <- "2024-12-31"

important.bond.tickers <- c('EMB','BND','IEI','TLT','HYG','AGG')
selected.bond.tickers <- c("EMB","HYG","AGG","BND",
                           "EMLC","EBND","LEMB","EMHY","ELD","VWOB","EMTL")

conf.niveles.rolling <- c(0.01,0.025,0.975,0.99)  # bandas VaR del rolling GARCH
seed.garch.roll <- 333   # semilla simulaciones rolling (garch.oos.sim)
seed.mc         <- 1     # semilla forecast Monte-Carlo (garch.ts.forecast)

# Helper functions ----
# GARCH model fitting functions
# Perlin, M. S., Mastella, M., Vancin, D. F., & Ramos, H. P. (2021). A GARCH tutorial with R. Revista de Administração Contemporânea, 25(1), e200088. https://doi.org/10.1590/1982-7849rac2021200088
# https://github.com/msperlin/GARCH-RAC
# NOTA (refactor): extensión corregida de .r a .R (el archivo real es .R;
# en Linux las rutas distinguen mayúsculas/minúsculas).
source("GARCH_msperlin_functions.R")

# Funciones auxiliares del análisis (extraídas del script original)
source("refactor/functions.R")

# Procesamiento de datos ----
# gobix daily
gobix <- read.csv(url("https://www.bvrd.com.do/indice/Data/GobixDataIRP.csv")) # descarga
gobix <- gobix[,1:2] # se queda con dos primeras columnas
colnames(gobix)[1] <- "fecha" # renombra la columna 1
gobix$fecha <- lubridate::mdy(gobix[,1]) # la convierte en formato fechas mes-dia-año
gobix <- xts::as.xts(gobix[,-1], order.by=gobix$fecha) # convierte el formato en XTS
colnames(gobix)[1] <- "Close.GOBIX" # renombra la serie (indice al cierre del dia)

gobix$var.diaria <- CalculateReturns(Cl(gobix), method = "discrete") # diferencia simple
gobix$daily.vol <- diff(log(Cl(gobix))) # diferencia logarítmica
gobix$rolling.sd <- roll_sd(gobix$daily.vol, width = 65)*sqrt(252) # stdev 3 meses anualizada
gobix$std.dev <- roll_sd(gobix$var.diaria, width = 90)
gobix <- na.omit(gobix)
gobix$drawdown <- Drawdowns(gobix[,'daily.vol']) # estimación caída histórica

gobix.completo <- gobix
gobix.completo <- as.data.frame(gobix.completo["2014::2022"])

# Período de Prueba de Modelo + Validación: Out of Sample.
oos <- gobix["2019::2022"] 

# Período de Entrenamiento: In the Sample.
gobix <- gobix["2014::2022"] 
gobix.train <- gobix["2014::2018"] 

# Análisis: Comparativo bonos internacionales ----
getSymbols(important.bond.tickers, src = 'yahoo', from = fecha.inicio.datos, auto.assign = TRUE)

# Bonos: Diferencias diarias RV
bond.vol <- merge(merge(merge(Cl(TLT),merge(Cl(IEI),Cl(AGG))),Cl(HYG)),Cl(EMB))
bond.vol <- apply(bond.vol, MARGIN = 2, FUN = function(x){ data.table::shift(x, n = 1, type = "lead") / x-1})
bond.vol <- as.xts(bond.vol, order.by = lubridate::ymd(rownames(bond.vol)))

bond.vol.df <- na.omit(merge(bond.vol, gobix$var.diaria))
colnames(bond.vol.df)[1:6] <- c("TLT: US Treasuty vencimiento 20A",
                                "IEI: US Treasury vencimiento 5-7A",
                                "AGG: US Corporativos Grado Inversion",
                                "HYG: US Especulativos ",
                                "EMB: JPM EMBI",
                                "GOBIX.diario") 

bond.vol.df <- as.data.frame(bond.vol.df)

# Universo ampliado de bonos

bond.data <- list()
bond.performance <- list()
for(i in 1:length(selected.bond.tickers)){
  bond.data[[i]] <- getSymbols(selected.bond.tickers[i], src = 'yahoo', from = fecha.inicio.datos, to = fecha.fin.datos, auto.assign = FALSE)
  bond.data[[i]]$daily.ret <- CalculateReturns(Cl(bond.data[[i]]), method = "discrete")
}

bond.performance <- do.call(cbind,lapply(bond.data, FUN = function(x){Cl(x)}))
bond.performance <- merge(gobix[,c('Close.GOBIX')], bond.performance)
bond.performance <- na.omit(bond.performance)
colnames(bond.performance) <- c('GOBIXDR',selected.bond.tickers)
index.day <- which(bond.performance$GOBIXDR == max(na.omit(bond.performance$GOBIXDR)))

df1 <- na.omit(bond.performance[index.day:nrow(bond.performance),])
df.indexed <- as.data.frame(apply(as.data.frame(df1), MARGIN = 2, FUN = function(x){x/head(x)}))
df.indexed$date <- lubridate::ymd(rownames(df.indexed))
names(df.indexed)[1] <- "GOBIXDR"

# Drawdown desde ATH
performance.2022 <- tail(df.indexed[,-ncol(df.indexed)],1) 
rownames(performance.2022) <- ""


bonds.df1 <- as.data.frame(apply(as.data.frame(bond.performance), MARGIN = 2, FUN = function(x){CalculateReturns(x, method="discrete")}))
bonds.df1$date <- lubridate::ymd(index(bond.performance))
bonds.df1 <- as.xts(bonds.df1[,-ncol(bonds.df1)], order.by = lubridate::ymd(index(bond.performance)))

# Tabla resumen
bond.vol.df <- bonds.df1

colnames(bond.vol.df)[2:ncol(bond.vol.df)] <- c("EMB: JPM EMBI",
                                                "HYG: US Especulativos ",
                                                "AGG: US Corporativos Grado Inversion",
                                                "BND: Bloomberg US Aggreate",
                                                "VanEck JPM EM Local Currency Bond",
                                                "SPDR Bloomberg EM Local Bond",
                                                "iShares JPM EM Local Currency",
                                                "iShares JPM EM High Yield Bond",
                                                "WisdomTree EM Local Debt Fund",
                                                "Vanguard EM Government Bond",
                                                "DoubleLine EM Fixed Income") 

bond.vol.df <- as.data.frame(bond.vol.df)


## Correlaciones ----
# computa correlaciones (con p-value) y las transforma en tabla
correlation.tb.daily <- rcorr(as.matrix(na.omit(bonds.df1)))
correlation.tb.daily <- prep.corr.matrix(correlation.tb.daily$r, correlation.tb.daily$P)
correlation.tb.daily <- correlation.tb.daily[correlation.tb.daily$row == 'GOBIXDR',]
correlation.tb.daily <- correlation.tb.daily[,-1]
correlation.tb.daily[,2:3] <- round(correlation.tb.daily[,2:3],4)
colnames(correlation.tb.daily) <- c('ETF','coeficiente','p.value')

## Características de riesgos ----
# Aquí se computan datos presentados en las págs. 10-13.
# Riesgos benchmarks internacionales----
bond.risk <- list()
for(i in 1:length(selected.bond.tickers)){
  bond.risk[[i]] <- risk.summary(ticker = selected.bond.tickers[i], benchmark_ticker = "IEF", start.date = fecha.inicio.datos, end.date = fecha.fin.datos)
}

risk.resume <- do.call(cbind, lapply(bond.risk, FUN = function(x){x[[3]]}))
relative_performance.resume <- do.call(cbind, lapply(bond.risk, FUN = function(x){x[[5]]}))

colnames(risk.resume) <- selected.bond.tickers
colnames(relative_performance.resume) <- selected.bond.tickers

### Riesgos GOBIX y 7-10Y UST benchmark ----
benchmark <- getSymbols("IEF", src = 'yahoo', from = fecha.inicio.datos, to = fecha.fin.datos, auto.assign = FALSE)
benchmark$benchmark.daily.rets <- diff(log(Cl(benchmark)))

gobix.risk <- merge(gobix,benchmark)
gobix.risk <- na.omit(gobix.risk)

return.distributions <- table.Distributions(gobix.risk[,"var.diaria", drop=FALSE])
table.DDs <- table.Drawdowns(gobix.risk[,"var.diaria", drop=FALSE], top = 25, digits = 4, geometric = FALSE)
downside.risk <- table.DownsideRisk(gobix.risk[,"var.diaria", drop=FALSE], Rf=.02/12, MAR =.05/12, p=.95)
calendar.rets <- calendar.ReturnTable(gobix.risk[,"var.diaria", drop=FALSE])
annualized.rets <- table.SFM(Ra = gobix.risk[,"var.diaria", drop=FALSE], Rb= gobix.risk[,"benchmark.daily.rets", drop=FALSE], Rf=0, digits=3)
names(annualized.rets)[1] <- 'var.diaria'

summary.gobix.risk <- list(return.distributions, table.DDs, downside.risk, calendar.rets, annualized.rets)
gobix.risk.resume <- rbind(downside.risk,annualized.rets)
names(gobix.risk.resume)[1] <- "GOBIX"

# consolidación tabla 6: PAG 13.
risk_performance <- cbind(t(risk.resume),t(relative_performance.resume))
risk_performance <- t(risk_performance)
risk_performance <- cbind(gobix.risk.resume,
                          risk_performance)


# par(mfrow=c(1,3), mar=c(4,4,4,4))
# plot.quadrants(t(risk_performance), y.var='Maximum Drawdown', x.var='Historical ES (95%)', description = "MaxDD vs c-VAR" )
# plot.quadrants(t(risk_performance), y.var='R-squared', x.var='Historical ES (95%)', description = "UST10Y R2 vs c-VaR" )
plot.quadrants(t(risk_performance), y.var='Beta', x.var='Maximum Drawdown', description = "UST10Y Beta vs MaxDD")
plot.quadrants(t(risk_performance), y.var='Active Premium', x.var='Correlation', description = "vs. UST10Y")

# par(mfrow=c(1,1), mar=c(3,3,3,3))


# Características ST GOBIXDR ----
# Aquí se computan datos presentados en las págs. 14-15.
## Tabla autocorrelación diaria ----
autocorrelation.full <- table.Autocorrelation(gobix.completo[,'daily.vol',drop=FALSE])
autocorrelation.training <- table.Autocorrelation(gobix["2014::2018",'daily.vol',drop=FALSE])
autocorrelation.test <- table.Autocorrelation(gobix["2019::2021",'daily.vol',drop=FALSE])
autocorrelation.validation <-  table.Autocorrelation(gobix["2022",'daily.vol',drop=FALSE])

table.autocorrelation <- cbind( autocorrelation.full,
                                autocorrelation.training,
                                autocorrelation.test,
                                autocorrelation.validation)

# consolidación tabla 7: PAG 14.
colnames(table.autocorrelation) <- c('periodo.completo',
                                     'entrenamiento',
                                     'prueba.de.modelo',
                                     'validacion')

## Tabla LB-test ----
LB.full <- round(Box.test(gobix[,'daily.vol',drop=FALSE], lag = 5, type = "Ljung-Box")$p.value,4)
LB.training <- round(Box.test(gobix["2014::2018",'daily.vol',drop=FALSE], lag = 5, type = "Ljung-Box")$p.value,4)
LB.test <- round(Box.test(gobix["2019::2021",'daily.vol',drop=FALSE], lag = 5, type = "Ljung-Box")$p.value,4)
LB.validation <- round(Box.test(gobix["2022",'daily.vol',drop=FALSE], lag = 5, type = "Ljung-Box")$p.value,4)

# consolidación tabla 8: PAG 15.
table.LB <- matrix(NA, nrow=1, ncol = 4, dimnames = list("p-value",
                                                         c('periodo.completo',
                                                           'entrenamiento',
                                                           'prueba.de.modelo',
                                                           'validacion')))
table.LB[1,1] <- LB.full
table.LB[1,2] <- LB.training
table.LB[1,3] <- LB.test
table.LB[1,4] <- LB.validation

## Tabla ARCH test ----
arch.test.full <- do_arch_test(x = gobix$daily.vol, max_lag = 5)
arch.test.training <- do_arch_test(x = gobix["2014::2018",'daily.vol',drop=FALSE], max_lag = 5)
arch.test.test <- do_arch_test(x = gobix["2019::2021",'daily.vol',drop=FALSE], max_lag = 5)
arch.test.validation <- do_arch_test(x = gobix["2022",'daily.vol',drop=FALSE], max_lag = 5)

arch.test.full <- as.data.frame(arch.test.full)
arch.test.training <- as.data.frame(arch.test.training)
arch.test.test <- as.data.frame(arch.test.test)
arch.test.validation <- as.data.frame(arch.test.validation)

# consolidación tabla 9: PAG 15.
table.arch <- cbind( arch.test.full[,c(1,2)],
                     arch.test.training[,2],
                     arch.test.test[,2],
                     arch.test.validation[,2])

colnames(table.arch)[2:ncol(table.arch)] <- c('periodo.completo',
                                              'entrenamiento',
                                              'prueba.de.modelo',
                                              'validacion')

# GARCH calibration ----
# Aquí se computan datos presentados en las págs. 16-17.
## Calibración GARCH ----

## 32,400 model fit. Se toma tiempo
# genera el output Best_GARCH_Gobix1.csv
# NOTA (refactor): bloque costoso (horas de cómputo). Por defecto NO se ejecuta;
# sus resultados precomputados se leen más abajo desde Best_GARCH_Gobix1.csv.
if (RUN_FULL_CALIBRATION) {
max_lag_AR <- 5 # m - parametro AR
max_lag_MA <- 5 # n - parametro MA
max_lag_ARCH <- 5 # p - parametro ARCH
max_lag_GARCH <- 5 # q - parametro GARCH
dist_to_use <- c('norm','snorm','ged','std','sstd','jsu') # ver rugarch::ugarchspecs
models_to_estimate <- c('sGARCH', 'eGARCH', 'iGARCH', 'gjrGARCH', 'apARCH', 'csGARCH') # see rugarch::rugarchspec help for more

out <- find_best_arch_model(x = gobix.train[,'daily.vol'], 
                            type_models = models_to_estimate,
                            dist_to_use = dist_to_use,
                            max_lag_AR = max_lag_AR,
                            max_lag_MA = max_lag_MA,
                            max_lag_ARCH = max_lag_ARCH,
                            max_lag_GARCH = max_lag_GARCH)

## tabla resumen con resultados
tab_out <- out$tab_out
#print(tab_out)

models_names <- unique(tab_out$model_name)
best_models <- c(tab_out$model_name[which.min(tab_out$AIC)],
                 tab_out$model_name[which.min(tab_out$BIC)])

# Define el modelo que minimiza el BIC
best_spec = ugarchspec(variance.model = list(model =  out$best_bic$type_model, 
                                             garchOrder = c(out$best_bic$lag_arch,
                                                            out$best_bic$lag_garch)),
                       mean.model = list(armaOrder = c(out$best_bic$lag_ar, 
                                                       out$best_bic$lag_ma)),
                       distribution = out$best_bic$type_dist)
}

# carga el resultado del proceso anterior
GOBIXDR.GARCH.best_model <- read.csv("Best_GARCH_Gobix1.csv")
GOBIXDR.GARCH.best_model <- na.omit(GOBIXDR.GARCH.best_model) # drop non-converged
GOBIXDR.GARCH.best_model$AIC_BIC_mean <- as.numeric(apply(GOBIXDR.GARCH.best_model[,c('AIC','BIC')], MARGIN = 1, FUN = mean))

GOBIXDR.GARCH.best_model <- GOBIXDR.GARCH.best_model %>% dplyr::arrange(AIC_BIC_mean)
GOBIXDR.GARCH.best_model <- head(GOBIXDR.GARCH.best_model, floor(nrow(GOBIXDR.GARCH.best_model)*0.01))

oos.df <- matrix(NA, nrow=nrow(GOBIXDR.GARCH.best_model), ncol = 25)
colnames(oos.df) <- c('model','training.AIC','training.BIC','uncond.mean','uncond.variance','persistence','half-life',
                      'excesos.esperados','excesos.registrados','H0', 'UC.LR.stat',
                      'UC.critico','UC.LRp','UC.decision','CC.H0','CC_LR.stat','CC.critico',
                      'CC.LRp','CC.decision','b','uLL','rLL','LRp','H0.CCduration','decision')
oos.df <- as.data.frame(oos.df)

## GARCH model test, validation ----
my_best_garch <- list()
best_spec <- list()
garch.roll <- list()
unconditioned.coverage.test <- list()
conditioned.coverage.test <- list()

## Corre el modelo, hace la prueba y validación
# NOTA (refactor): bloque costoso (horas de cómputo: rolling backtest del Top 1%
# de modelos). Por defecto NO se ejecuta; sus resultados precomputados se leen
# más abajo desde Top1PCT_model_validation_test.csv.
if (RUN_FULL_CALIBRATION) {
for(i in 1:nrow(GOBIXDR.GARCH.best_model)){
  
  print(paste("PROGRESS:", i / nrow(GOBIXDR.GARCH.best_model)))
  print(GOBIXDR.GARCH.best_model[i,])
  
  best_spec[[i]] = ugarchspec(variance.model = list(model =  GOBIXDR.GARCH.best_model$type_model[i], 
                                               garchOrder = c(GOBIXDR.GARCH.best_model$lag_arch[i],
                                                              GOBIXDR.GARCH.best_model$lag_garch[i])),
                              mean.model = list(armaOrder = c(GOBIXDR.GARCH.best_model$lag_ar[i],
                                                              GOBIXDR.GARCH.best_model$lag_ma[i])),
                              distribution = GOBIXDR.GARCH.best_model$type_dist[i])
  
  my_best_garch[[i]] <- ugarchfit(spec = best_spec[[i]], data = gobix["2014/2018",'daily.vol'])
  
  symbol = "GOBIXDR"
  n=262
  
  # CALIBRACION CORRIDA Y VALIDACION
  # Ejecuta la predicción de 1-día hacía adelante (rolling GARCH) con recalibración cada 30 días
  garch.roll[[i]] = ugarchroll(best_spec[[i]], gobix["2019/2022",'daily.vol'], n.ahead=1,
                               forecast.length = n, solver = "hybrid", refit.every=30, 
                               refit.window="moving", VaR.alpha=c(0.01,0.99))
  
  unconditioned.coverage.test[[i]] <- VaRTest(alpha=0.01, conf.level = 0.99, actual = garch.roll[[i]]@forecast$VaR$realized, VaR = garch.roll[[i]]@forecast$VaR$`alpha(1%)` )
  conditioned.coverage.test[[i]] <- VaRDurTest(alpha=0.01, conf.level = 0.99, actual = garch.roll[[i]]@forecast$VaR$realized, VaR = garch.roll[[i]]@forecast$VaR$`alpha(1%)` )
  
  oos.df$model[i] <- as.character(GOBIXDR.GARCH.best_model$model_name[i])
  oos.df$training.AIC[i] <- as.character(GOBIXDR.GARCH.best_model$AIC[i])
  oos.df$training.BIC[i] <- as.character(GOBIXDR.GARCH.best_model$BIC[i])
  oos.df$uncond.mean[i] <- uncmean(my_best_garch[[i]])
  oos.df$uncond.variance[i] <- uncvariance(my_best_garch[[i]])
  oos.df$persistence[i] <- persistence(my_best_garch[[i]])
  oos.df$`half-life`[i] <- halflife(my_best_garch[[i]])
  oos.df$excesos.esperados[i] <- as.numeric(unconditioned.coverage.test[[i]][1])
  oos.df$excesos.registrados[i] <- as.numeric(unconditioned.coverage.test[[i]][2])
  oos.df$H0[i] <- as.character(unconditioned.coverage.test[[i]][3])
  oos.df$UC.LR.stat[i] <- as.numeric(unconditioned.coverage.test[[i]][4])
  oos.df$UC.critico[i] <- as.numeric(unconditioned.coverage.test[[i]][5])
  oos.df$UC.LRp[i] <- as.numeric(unconditioned.coverage.test[[i]][6])
  oos.df$UC.decision[i] <- as.character(unconditioned.coverage.test[[i]][7])
  oos.df$CC.H0[i] <- as.character(unconditioned.coverage.test[[i]][8])
  oos.df$CC_LR.stat[i] <- as.numeric(unconditioned.coverage.test[[i]][9])
  oos.df$CC.critico[i] <- as.numeric(unconditioned.coverage.test[[i]][10])
  oos.df$CC.LRp[i] <- as.numeric(unconditioned.coverage.test[[i]][11])
  oos.df$CC.decision[i] <- as.character(unconditioned.coverage.test[[i]][12])
  oos.df$b[i] <- as.numeric(conditioned.coverage.test[[i]][1])
  oos.df$uLL[i] <- as.numeric(conditioned.coverage.test[[i]][2])
  oos.df$rLL[i] <- as.numeric(conditioned.coverage.test[[i]][3])
  oos.df$LRp[i] <- as.numeric(conditioned.coverage.test[[i]][4])
  oos.df$H0.CCduration[i] <- as.character(conditioned.coverage.test[[i]][5])
  oos.df$decision[i] <- as.character(conditioned.coverage.test[[i]][6])
}
}

best_model.resume <- GOBIXDR.GARCH.best_model

# RESUME CONFIGURACIONES (m,n)
best_model.resume$ar_ma <- paste(best_model.resume$lag_ar, best_model.resume$lag_ma)
  ar.resume <- table(best_model.resume$ar_ma)
  ar.resume_1 <- round(prop.table(ar.resume),3)
  ar.resume <- t(rbind(t(as.data.frame(ar.resume)),t(as.data.frame(ar.resume_1))))
  ar.resume <- as.data.frame(ar.resume[,-3])
  colnames(ar.resume) <- c("ar,ma","freq_ma","freq_ma%") 
  ar.resume[,'freq_ma'] <- as.numeric(ar.resume[,'freq_ma'])
  ar.resume[,'freq_ma%'] <- as.numeric(ar.resume[,'freq_ma%'])
  ar.resume <- ar.resume %>% arrange(freq_ma, decreasing=TRUE)
  ar.resume$ar.cumsum <- cumsum(ar.resume[,'freq_ma%'])
  
# RESUME CONFIGURACIONES (p,q)
best_model.resume$p_q <- paste(best_model.resume$lag_arch, best_model.resume$lag_garch)
  pq.resume <- table(best_model.resume$p_q)
  pq.resume_1 <- round(prop.table(pq.resume),3)
  pq.resume <- t(rbind(t(as.data.frame(pq.resume)),t(as.data.frame(pq.resume_1))))
  pq.resume <- as.data.frame(pq.resume[,-3])
  colnames(pq.resume) <- c("p,q","freq_pq","freq_pq%") 
  pq.resume[,'freq_pq'] <- as.numeric(pq.resume[,'freq_pq'])
  pq.resume[,'freq_pq%'] <- as.numeric(pq.resume[,'freq_pq%'])
  pq.resume <- pq.resume %>% arrange(freq_pq, decreasing=TRUE)
  pq.resume$pq_cumsum <- cumsum(pq.resume[,'freq_pq%'])

# RESUME TABLA CON TOP80% PARAMS
best_model_freq.resume <- cbind(ar.resume[1:20,],pq.resume[1:20,])

# RESUME LOS MEJORES GARCH
garch_type.resume <- table(best_model.resume$type_model)
  garch_type.resume_1 <- round(prop.table(garch_type.resume),3)
  garch_type.resume <- t(rbind(t(as.data.frame(garch_type.resume)),t(as.data.frame(garch_type.resume_1))))
  garch_type.resume <- as.data.frame(garch_type.resume[,-3])
  colnames(garch_type.resume) <- c("model_type","freq_model","freq_model%") 
  garch_type.resume[,'freq_model'] <- as.numeric(garch_type.resume[,'freq_model'])
  garch_type.resume[,'freq_model%'] <- as.numeric(garch_type.resume[,'freq_model%'])
  garch_type.resume <- garch_type.resume %>% arrange(freq_model, decreasing=TRUE)
  garch_type.resume$garch.cumsum <- cumsum(garch_type.resume[,'freq_model%'])
  
# RESUME LAS MEJORES DISTRIBUCINOES DEL TERMINO ERROR
dist_type.resume <- table(best_model.resume$type_dist)
  dist_type.resume_1 <- round(prop.table(dist_type.resume),3)
  dist_type.resume <- t(rbind(t(as.data.frame(dist_type.resume)),t(as.data.frame(dist_type.resume_1))))
  dist_type.resume <- as.data.frame(dist_type.resume[,-3])
  colnames(dist_type.resume) <- c("dist_type","freq_dist","freq_dist%") 
  dist_type.resume[,'freq_dist'] <- as.numeric(dist_type.resume[,'freq_dist'])
  dist_type.resume[,'freq_dist%'] <- as.numeric(dist_type.resume[,'freq_dist%'])
  dist_type.resume <- dist_type.resume %>% arrange(freq_dist, decreasing=TRUE)
  dist_type.resume$dist.cumsum <- cumsum(dist_type.resume[,'freq_dist%'])

## Resultado prueba de modelo ----  
# Toma el resultado del proceso anterior
oos.df <- read.csv("Top1PCT_model_validation_test.csv")
oos.df.clean <- na.omit(oos.df)
oos.df.passed <- oos.df.clean %>% dplyr::filter(UC.decision == 'Fail to Reject H0' &
                                                CC.decision == 'Fail to Reject H0' &
                                                decision == 'Fail to Reject H0')

oos.df.passed$training.AIC <- as.numeric(oos.df.passed$training.AIC)
oos.df.passed$training.BIC <- as.numeric(oos.df.passed$training.BIC)
oos.df.passed$avg_AIC.BIC <- apply(oos.df.passed[,c('training.AIC','training.BIC')], MARGIN = 1, FUN = mean)
oos.df.passed <- oos.df.passed %>% dplyr::select(model:training.BIC, avg_AIC.BIC, uncond.mean:decision)

# consolidación tabla 24: PAG 32 (anexos).
oos.df.resumed <- oos.df.passed %>% dplyr::select(model:training.BIC,
                                                 excesos.esperados,
                                                 excesos.registrados,
                                                 UC.LR.stat:UC.LRp,
                                                 CC_LR.stat:CC.LRp,
                                                 b:LRp)

## Best fit por familia GARCH ----
iGARCH_best <- oos.df.passed[str_detect(oos.df.passed$model, "iGARCH"), ] 
iGARCH_best <- iGARCH_best %>% arrange(avg_AIC.BIC)
iGARCH_best[1,] 

sGARCH_best <- oos.df.passed[str_detect(oos.df.passed$model, "sGARCH"), ] 
sGARCH_best <- sGARCH_best %>% arrange(avg_AIC.BIC)
sGARCH_best[1,]

csGARCH_best <- oos.df.passed[str_detect(oos.df.passed$model, "csGARCH"), ] 
csGARCH_best <- csGARCH_best %>% arrange(avg_AIC.BIC)
csGARCH_best[1,]

backtest_best.models <- rbind(rbind(iGARCH_best[1,],sGARCH_best[1,]),csGARCH_best[1,])

x <- intersect(backtest_best.models$model, GOBIXDR.GARCH.best_model$model_name)

final.models <- list()
final.models.row <- list()
for(i in 1:length(x)){
  final.models.row[[i]] <- which(GOBIXDR.GARCH.best_model[,'model_name'] == x[i])
  final.models[[i]] <- GOBIXDR.GARCH.best_model[which(GOBIXDR.GARCH.best_model[,'model_name'] == x[i]), ]
}

final.models.row <- do.call(rbind, final.models.row)
final.models <- as.data.frame(do.call(rbind, final.models))

## Entrenamiento Top 3 modelos ----
training.model <- list()
training.best_garch <- list()
for(i in 1:nrow(final.models)){
  training.model[[i]] = ugarchspec(variance.model = list(model =  final.models$type_model[i], 
                                                    garchOrder = c(final.models$lag_arch[i],
                                                                   final.models$lag_garch[i])),
                              mean.model = list(armaOrder = c(final.models$lag_ar[i],
                                                              final.models$lag_ma[i])),
                              distribution = final.models$type_dist[i])
  
  training.best_garch[[i]] <- ugarchfit(spec = training.model[[i]], data = gobix.train["2014/2018",'daily.vol'])
  
}

# GARCH diario prueba de modelo
# Contiene gráficos pág. 22
n = nrow(final.models)
GARCH.daily.marginal <- list()
for(i in 1:n){
  #print(paste0(i,":",final.models$model_name[i]))
  # ROLLING GARCH
  GARCH.daily.marginal[[i]] <- garch.oos.sim(model.df = final.models[i,],
                                             model.description = final.models$model_name[i],
                                             training.data = gobix["2014/2021",'daily.vol'],
                                             method = "rolling",
                                             confidence = conf.niveles.rolling,
                                             band.colors = c("red4","red","green","green4"),
                                             rolling.length = nrow(gobix["2019/2021",'daily.vol']),
                                             no.simulations = 10000,
                                             set.seed.n = seed.garch.roll,
                                             show.limit = 0.10,
                                             oos.realized = gobix["2019/2021",'daily.vol'])
  #print(GARCH.daily.marginal[[i]])
}

# GARCH diario validación
GARCH.daily.validation <- list()
for(i in 1:n){
  #print(paste0(i,":",final.models$model_name[i]))
  # ROLLING GARCH
  GARCH.daily.validation[[i]] <- garch.oos.sim(model.df = final.models[i,],
                                               model.description = final.models$model_name[i], 
                                               training.data = gobix["2014/2022",'daily.vol'],
                                               method = "rolling", 
                                               confidence=conf.niveles.rolling,
                                               band.colors = c("red4","red","green","green4"),
                                               rolling.length = nrow(gobix["2022",'daily.vol']), 
                                               no.simulations = 10000,
                                               set.seed.n = seed.garch.roll,
                                               show.limit = 0.10, 
                                               oos.realized = gobix["2022",'daily.vol']) 
  #print(GARCH.daily.validation[[i]])
}


# GARCH model stats ----

# Entrenamiento (ITS)
df.sigma.entrenamiento <- as.data.frame(cbind(sigma(training.best_garch[[1]]),
                                              sigma(training.best_garch[[2]]),
                                              sigma(training.best_garch[[3]])))
df.sigma.entrenamiento <- na.omit(df.sigma.entrenamiento)
colnames(df.sigma.entrenamiento) <- c("iGARCH","sGARCH","csGARCH")
df.sigma.entrenamiento$date <- lubridate::ymd(rownames(df.sigma.entrenamiento))

# Prueba de modelo (OOS)
df.sigma_test <- as.data.frame(cbind(GARCH.daily.marginal[[1]][[4]]$Sigma,
                                     GARCH.daily.marginal[[2]][[4]]$Sigma,
                                     GARCH.daily.marginal[[3]][[4]]$Sigma))

colnames(df.sigma_test) <- c("iGARCH","sGARCH","csGARCH")
date <- lubridate::ymd(rownames(GARCH.daily.marginal[[1]][[4]]))
df.sigma_test <- cbind(date,df.sigma_test)

# Validacion (OOS)
df.sigma_validacion<- as.data.frame(cbind(GARCH.daily.validation[[1]][[4]]$Sigma,
                                          GARCH.daily.validation[[2]][[4]]$Sigma,
                                          GARCH.daily.validation[[3]][[4]]$Sigma))

colnames(df.sigma_validacion) <- c("iGARCH","sGARCH","csGARCH")
date <- lubridate::ymd(rownames(GARCH.daily.validation[[1]][[4]]))
df.sigma_validacion <- cbind(date,df.sigma_validacion)

df.sigma <- rbind(df.sigma.entrenamiento,
                  df.sigma_test,
                  df.sigma_validacion)

# TABLA RESUMIDA
resume.sigma.entrenamiento <- melt(df.sigma.entrenamiento, id.vars = "date")
resume.sigma.entrenamiento$variable <- as.character(resume.sigma.entrenamiento$variable)
resume.sigma.entrenamiento <- resume.sigma.entrenamiento %>% 
                                      group_by(variable) %>% 
                                      summarise(entrenamiento_media = mean(value)*sqrt(250),
                                                entrenamiento_q.975 = quantile(value, 0.975)*sqrt(250),
                                                entrenamiento_q.99 = quantile(value, 0.99)*sqrt(250))

resume.sigma.entrenamiento[,2:ncol(resume.sigma.entrenamiento)] <- round(as.data.frame(resume.sigma.entrenamiento[,2:ncol(resume.sigma.entrenamiento)]),4)
resume.sigma.entrenamiento <- t(resume.sigma.entrenamiento)

resume.sigma.test <- melt(df.sigma_test, id.vars = "date")
resume.sigma.test$variable <- as.character(resume.sigma.test$variable)
resume.sigma.test <- resume.sigma.test %>% 
  group_by(variable) %>% 
  summarise(prueba.modelo_media = mean(value)*sqrt(250),
            prueba.modelo_q.975 = quantile(value, 0.975)*sqrt(250),
            prueba.modelo_q.99 = quantile(value, 0.99)*sqrt(250))

resume.sigma.test[,2:ncol(resume.sigma.test)] <- round(as.data.frame(resume.sigma.test[,2:ncol(resume.sigma.test)]),4)
resume.sigma.test <- t(resume.sigma.test)

resume.sigma.validacion<- melt(df.sigma_validacion, id.vars = "date")
resume.sigma.validacion$variable <- as.character(resume.sigma.validacion$variable)
resume.sigma.validacion <- resume.sigma.validacion %>% 
  group_by(variable) %>% 
  summarise(validacion_media = mean(value)*sqrt(250),
            validacion_q.975 = quantile(value, 0.975)*sqrt(250),
            validacion_q.99 = quantile(value, 0.99)*sqrt(250))

resume.sigma.validacion[,2:ncol(resume.sigma.validacion)] <- round(as.data.frame(resume.sigma.validacion[,2:ncol(resume.sigma.validacion)]),4)
resume.sigma.validacion <- t(resume.sigma.validacion)

resume.sigma <- rbind(as.data.frame(resume.sigma.entrenamiento),
                      as.data.frame(resume.sigma.test[2:nrow(resume.sigma.test),]),
                      as.data.frame(resume.sigma.validacion[2:nrow(resume.sigma.validacion),]))

colnames(resume.sigma) <- resume.sigma[1,]
resume.sigma <- resume.sigma[-1,]

# SIGMA
sigma.model.1 <- stat.desc(sigma(training.best_garch[[1]]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
sigma.model.2 <- stat.desc(sigma(training.best_garch[[2]]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
sigma.model.3 <- stat.desc(sigma(training.best_garch[[3]]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)

  sigma.model.1 <- round(sigma.model.1,5)
  sigma.model.2 <- round(sigma.model.2,5)
  sigma.model.3 <- round(sigma.model.3,5)
  
  sigma.df <- cbind(cbind(sigma.model.1,sigma.model.2), sigma.model.3)
  colnames(sigma.df) <- final.models$type_model
  #View(sigma.df)
  

# ERROR TERM
error.model.1 <- stat.desc(residuals(training.best_garch[[1]]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
error.model.2 <- stat.desc(residuals(training.best_garch[[2]]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
error.model.3 <- stat.desc(residuals(training.best_garch[[3]]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
  
error.model.1 <- round(error.model.1,5)
error.model.2 <- round(error.model.2,5)
error.model.3 <- round(error.model.3,5)
  
  residual.df <- cbind(cbind(error.model.1,error.model.2), error.model.3)
  colnames(residual.df) <- final.models$type_model
  #View(residual.df)

# INFOCRITERIA
info.model.1 <- infocriteria(training.best_garch[[1]])
info.model.2 <- infocriteria(training.best_garch[[2]])
info.model.3 <- infocriteria(training.best_garch[[3]])

  info.df <- cbind(cbind(info.model.1,info.model.2), info.model.3)
  colnames(info.df) <- final.models$type_model

  # INFOCRITERIA
  info.model.1 <- infocriteria(training.best_garch[[1]])
  info.model.2 <- infocriteria(training.best_garch[[2]])
  info.model.3 <- infocriteria(training.best_garch[[3]])
  
  info.df <- cbind(cbind(info.model.1,info.model.2), info.model.3)
  colnames(info.df) <- final.models$type_model
  # knitr::kable(info.df, caption="Criterios informativos")
  
  # LIKELiHOOD
  likelihood.df <- as.data.frame(cbind(likelihood(training.best_garch[[1]]),
                                       likelihood(training.best_garch[[2]]),
                                       likelihood(training.best_garch[[3]])))
  colnames(likelihood.df) <- final.models$type_model
  rownames(likelihood.df) <- "Likelihood"
  # knitr::kable(likelihood.df, caption="Likelihood")
  
  # JOINT STATS
  asymptotics.IC <- cbind(as.data.frame(nyblom(training.best_garch[[1]])$IndividualCritical),
                          as.data.frame(nyblom(training.best_garch[[2]])$IndividualCritical),
                          as.data.frame(nyblom(training.best_garch[[3]])$IndividualCritical))
  colnames(asymptotics.IC) <- colnames(likelihood.df)
  # kable(asymptotics.IC, caption="Hansen-Nyblom Stability Test: Individual Critical Values")
  
  asymptotics.JC <- cbind(as.data.frame(nyblom(training.best_garch[[1]])$JointCritical),
                          as.data.frame(nyblom(training.best_garch[[2]])$JointCritical),
                          as.data.frame(nyblom(training.best_garch[[3]])$JointCritical))
  colnames(asymptotics.JC) <- colnames(likelihood.df)
  # knitr::kable(asymptotics.JC, caption="Hansen-Nyblom Stability Test: Joint Critical Values")
  
  # SIGN-BIAS
  sign.bias.df <- cbind(as.data.frame(signbias(training.best_garch[[1]])[,-2]),
                        as.data.frame(signbias(training.best_garch[[2]])[,-2]),
                        as.data.frame(signbias(training.best_garch[[3]])[,-2]))
  sign.bias.df <- sign.bias.df[,-c(2,4,6)]
  colnames(sign.bias.df) <- colnames(likelihood.df)
  # knitr::kable(sign.bias.df, caption="Sesgo: t-value")
  
  # Goodness-of-Fit
  gof.test.df <- cbind(as.data.frame(gof(training.best_garch[[1]],c(10,20,30,40,50))[,c(1,3)]),
                       as.data.frame(gof(training.best_garch[[2]],c(10,20,30,40,50))[,c(3)]),
                       as.data.frame(gof(training.best_garch[[3]],c(10,20,30,40,50))[,c(3)]))
  colnames(gof.test.df)[2:ncol(gof.test.df)] <- colnames(likelihood.df)
  # knitr::kable(gof.test.df, caption="Pearson Goodness-of-Fit Ajustado: p-values")
  
  # merge together
  merged.df <- rbind(info.df,likelihood.df)
  merged.df <- rbind(merged.df,asymptotics.IC)
  rownames(merged.df)[6:8] <- c('Individual Critical Values.10%',
                                'Individual Critical Values.5%',
                                'Individual Critical Values.1%')
  merged.df <- rbind(merged.df,asymptotics.JC)
  rownames(merged.df)[9:11] <- c('Joint Critical Values.10%',
                                 'Joint Critical Values.5%',
                                 'Joint Critical Values.1%')
  merged.df <- rbind(merged.df,sign.bias.df)
  rownames(merged.df)[12:15] <- c('Sign Bias (t-value)',
                                  'Negative Sign Bias (t-value)',
                                  'Positive Sign Bias (t-value)',
                                  'Joint Effect (t-value)')
  merged.df <- rbind(merged.df,gof.test.df[,2:4])
  rownames(merged.df)[16:20] <- c('Goodness of Fit.10 (p-value)',
                                  'Goodness of Fit.20 (p-value)',
                                  'Goodness of Fit.30 (p-value)',
                                  'Goodness of Fit.40 (p-value)',
                                  'Goodness of Fit.50 (p-value)')

# Forecast + functions ----
# NOTA (refactor): se eliminaron aquí bloques comentados que duplicaban los bucles
# GARCH.daily.marginal / GARCH.daily.validation ya ejecutados arriba, y un bucle
# "PEND" abandonado (GARCH.MC.bootstrapped.sim).

# GARCH Monte-Carlo Bootstrapped OOS.
GARCH.MC.bootstrapped <- list()
for(i in 1:n){
  print(paste0(i,":",final.models$model_name[i]))
  GARCH.MC.bootstrapped[[i]] <- garch.ts.forecast(xts.object=gobix, label="GOBIXDR",
                                                  garch.family= final.models$type_model[i],
                                                  arma.params=c(final.models$lag_ar[i],
                                                                final.models$lag_ma[i]),
                                                  garch.params=c(final.models$lag_arch[i],
                                                                 final.models$lag_garch[i]),
                                                  distribution=final.models$type_dist[i],
                                                  set.seed.n = seed.mc,
                                                  start.date=fecha.inicio.mc, end.date=fecha.fin.mc,
                                                  days.forward=262, n.sim=10000, plot = TRUE)
  #print(GARCH.MC.bootstrapped[[i]])
}

length(GARCH.MC.bootstrapped[[3]])
GARCH.MC.bootstrapped[[2]][[3]]


cdf.1 <- na.omit(as.numeric(GARCH.MC.bootstrapped[[1]][[1]][nrow(GARCH.MC.bootstrapped[[1]][[1]]),]))
cdf.2 <- na.omit(as.numeric(GARCH.MC.bootstrapped[[2]][[1]][nrow(GARCH.MC.bootstrapped[[1]][[1]]),]))
cdf.3 <- na.omit(as.numeric(GARCH.MC.bootstrapped[[3]][[1]][nrow(GARCH.MC.bootstrapped[[1]][[1]]),]))

level.prob.1 <- quantInv(cdf.1,as.numeric(tail(Cl(oos),1)))
level.prob.2 <- quantInv(cdf.2,as.numeric(tail(Cl(oos),1)))
level.prob.3 <- quantInv(cdf.3,as.numeric(tail(Cl(oos),1)))

# Plot 3 CDFs - numeric
par(mfrow=c(1,1), mar=c(4,4,4,2))
plot(ecdf(cdf.1),
     xlab="Nivel Esperado",
     ylab="Probabilidad (Función Densidad Cumulativa)",
     xlim=c(80,110),
     ylim=c(0,0.05),
     main="Probabilidad teoríca de alcanzar el nivel de precio realizado",
     col="green", lwd=1.5)
lines(ecdf(cdf.2), col="red")
lines(ecdf(cdf.3), col="blue")
mtext("Validación (día T): Simulación vs. realización (out-of-sample)",  side=3)
abline(v=tail(Cl(oos),1), col="green4", lty=3)
legend("topleft", legend=c(paste0("iGARCH ",round(level.prob.1*100,2),"%"),
                           paste0("sGARCH ", round(level.prob.2*100,2),"%"),
                           paste0("csGARCH ",round(level.prob.3*100,2),"%")),
       col=c("green","red","blue"), lwd=1, lty=1, bty="n")
points(x = as.numeric(tail(Cl(oos),1)), y = level.prob.1, pch=21, cex=1.25, col="white", bg="green")
points(x = as.numeric(tail(Cl(oos),1)), y = level.prob.2, pch=21, cex=1.25, col="white", bg="red")
points(x = as.numeric(tail(Cl(oos),1)), y = level.prob.3, pch=21, cex=1.25, col="white", bg="blue")
box(col="grey")
par(mfrow=c(1,1), mar=c(4,4,4,4))


# Transform in % loss.
last.observation <- as.numeric(Cl(gobix["2021-12-31"]))

cdf.1_percent <- apply(as.data.frame(GARCH.MC.bootstrapped[[1]][[1]]), MARGIN = 2, FUN = function(x){x/last.observation-1})
cdf.1_percent <- na.omit(as.numeric(cdf.1_percent[nrow(cdf.1_percent),]))
cdf.2_percent <- apply(as.data.frame(GARCH.MC.bootstrapped[[2]][[1]]), MARGIN = 2, FUN = function(x){x/last.observation-1})
cdf.2_percent <- na.omit(as.numeric(cdf.2_percent[nrow(cdf.2_percent),]))
cdf.3_percent <- apply(as.data.frame(GARCH.MC.bootstrapped[[3]][[1]]), MARGIN = 2, FUN = function(x){x/last.observation-1})
cdf.3_percent <- na.omit(as.numeric(cdf.3_percent[nrow(cdf.3_percent),]))

realized.loss <- as.numeric(tail(Cl(oos),1)) / last.observation - 1
level.prob.1 <- quantInv(cdf.1_percent,realized.loss)
level.prob.2 <- quantInv(cdf.2_percent,realized.loss)
level.prob.3 <- quantInv(cdf.3_percent,realized.loss)

# Plot 3 CDFs - relative
par(mfrow=c(1,1), mar=c(4,4,4,2))
plot(ecdf(cdf.1_percent),
     xlab="Pérdida Esperada (Cumulativa)",
     ylab="Probabilidad (Función Densidad Cumulativa)",
     xlim=c(-0.3,-0.1),
     ylim=c(0,0.05),
     main="Probabilidad teoríca de alcanzar P(T)",
     col="white", lwd=1.5)
lines(ecdf(cdf.1_percent), col="green")
lines(ecdf(cdf.2_percent), col="red")
lines(ecdf(cdf.3_percent), col="blue")
mtext("Validación (T = 2022-12-31): Simulado vs. realizado",  side=3)
abline(v=realized.loss, col="green4", lty=3)
legend("topleft", legend=c(paste0("iGARCH ",round(level.prob.1*100,2),"%"),
                           paste0("sGARCH ", round(level.prob.2*100,2),"%"),
                           paste0("csGARCH ",round(level.prob.3*100,2),"%")),
       col=c("green","red","blue"), lwd=1, lty=1, pt.lwd = 0.5,  bty="n", xjust = 1)
points(x = realized.loss, y = level.prob.1, pch=21, cex=1.25, col="white", bg="green")
points(x = realized.loss, y = level.prob.2, pch=21, cex=1.25, col="white", bg="red")
points(x = realized.loss, y = level.prob.3, pch=21, cex=1.25, col="white", bg="blue")
box(col="grey")
par(mfrow=c(1,1), mar=c(4,4,4,4))


# stress from current levels
#VaR.bootstrap <- apply(GARCH.MC.bootstrapped[[1]], MARGIN = 1, FUN = function(x){quantile(na.omit(x),0.20)})
#as.numeric(tail(VaR.bootstrap,1)) / as.numeric(head(Cl(gobix["2022"]),1)) - 1

# Conditional
# This would be very interesting thing to do on many assets,
# and the perform a quantile regression.
condtional.excess <- GARCH.daily.marginal[[1]][[2]]$realized - GARCH.daily.marginal[[1]][[2]]$`alpha(1%)`
relative.rank <- GARCH.daily.marginal[[1]][[2]]$var01.rank
conditional.df <- as.data.frame(cbind(condtional.excess,relative.rank))

# NOTA (refactor): redefinición intencional de plot.quadrants (versión con
# label.position). Se mantiene aquí, en su posición original, para preservar
# el comportamiento del script (functions.R contiene la primera versión).
plot.quadrants <- function(data, x.var, y.var, label.position, description){
  
  require(scales)
  
  asset.data <- na.omit(data)
  
  a <- lm(coredata(asset.data[,y.var])~coredata(asset.data[,x.var]))
  spline <- smooth.spline(x=coredata(asset.data[,x.var]), y=coredata(asset.data[,y.var]), df=5)
  
  #par(mfrow=c(1,1), mar=c(4,4,4,4), xpd=FALSE) 
  par(xpd=FALSE) 
  plot(coredata(asset.data[,x.var]), coredata(asset.data[,y.var]), 
       main=description, ylab=y.var, xlab=x.var, pch=16, cex=0.85, col = scales::alpha("grey", 0.60))
  #mtext(description, side=3)
  points(as.numeric(tail(asset.data[,x.var],1)), as.numeric(tail(asset.data[,y.var],1)), col="red", pch=16, cex=1.25)
  abline(v=mean(coredata(asset.data[,x.var])), lty=2, lwd=0.75, col="grey")
  abline(h=mean(coredata(asset.data[,y.var])), lty=2, lwd=0.75, col="grey")
  text(x=min(asset.data[,x.var]), y=mean(asset.data[,y.var])-sd(asset.data[,y.var])/3, labels="Quadrant III", pos=4, cex=0.8, col="red")
  text(x=max(asset.data[,x.var]), y=mean(asset.data[,y.var])+sd(asset.data[,y.var])/3, labels="Quadrant I", pos=2, cex=0.8, col="green")
  text(x=max(asset.data[,x.var]), y=mean(asset.data[,y.var])-sd(asset.data[,y.var])/3, labels="Quadrant IV", pos=2, cex=0.8, col="grey40")
  text(x=min(asset.data[,x.var]), y=mean(asset.data[,y.var])+sd(asset.data[,y.var])/3, labels="Quadrant II", pos=4, cex=0.8, col="grey40")
  abline(a, col="red", lwd=2)
  #lines(spline, col="purple", lwd=2, lty=2)
  legend(label.position, legend=c(paste("y =", round(summary(a)$coefficients[1,1],3),"+",round(summary(a)$coefficients[2,1],3),"x"),
                             paste("y.t-value:", round(summary(a)$coefficients[1,3],3)),
                             paste("x.std error:", round(summary(a)$coefficients[2,2],3)),
                             paste("x.t-value:", round(summary(a)$coefficients[2,3],3)),
                             paste("r.squared:", round(summary(a)$r.squared,3))), bty="n")
  box(col = "grey")
}
plot.quadrants(data = conditional.df, 
               x.var = "relative.rank", 
               y.var = "condtional.excess", 
               label.position = "bottomright",
               description = "Conditional Excess on VaR Strength")

# VaR TABLES: Model TEST (2019-2021) ----

# oos_historical.VaR.01 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[1]]$daily.vol, 0.01)})
# oos_historical.cVaR.01 <- lapply(GARCH.daily.marginal, FUN = function(x){cVaR(x[[1]]$daily.vol, 0.01)})
# oos_historical.VaR.025 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[1]]$daily.vol, 0.025)})
# oos_historical.cVaR.025 <- lapply(GARCH.daily.marginal, FUN = function(x){cVaR(x[[1]]$daily.vol, 0.025)})

# REGULAR VaR
oos_historical.VaR.01 <- apply.fromstart(GARCH.daily.marginal[[1]][[1]][,'daily.vol'], FUN = function(x){quantile(x, 0.01)}, gap = 30)
oos_historical.VaR.025 <- apply.fromstart(GARCH.daily.marginal[[1]][[1]][,'daily.vol'], FUN = function(x){quantile(x, 0.025)}, gap = 30)

VaR.df_model.test <- matrix(nrow = 8, ncol = 6)
colnames(VaR.df_model.test) <- c('min','max','q.50','q.05','q.025','q.01')
rownames(VaR.df_model.test) <- c('historical.VaR.01',
                  'iGARCH.VaR.01',
                  'sGARCH.VaR.01',
                  'csGARCH.VaR.01',
                  'historical.VaR.025',
                  'iGARCH.VaR.025',
                  'sGARCH.VaR.025',
                  'csGARCH.VaR.025')

VaR.df_model.test[1,1] <- round(min(na.omit(oos_historical.VaR.01)),4)
VaR.df_model.test[1,2] <- round(max(na.omit(oos_historical.VaR.01)),4)
VaR.df_model.test[1,3] <- round(quantile(na.omit(oos_historical.VaR.01),0.5),4)
VaR.df_model.test[1,4] <- round(quantile(na.omit(oos_historical.VaR.01),0.05),4)
VaR.df_model.test[1,5] <- round(quantile(na.omit(oos_historical.VaR.01),0.025),4)
VaR.df_model.test[1,6] <- round(quantile(na.omit(oos_historical.VaR.01),0.01),4)

VaR.df_model.test[5,1] <- round(min(na.omit(oos_historical.VaR.025)),4)
VaR.df_model.test[5,2] <- round(max(na.omit(oos_historical.VaR.025)),4)
VaR.df_model.test[5,3] <- round(quantile(na.omit(oos_historical.VaR.025),0.5),4)
VaR.df_model.test[5,4] <- round(quantile(na.omit(oos_historical.VaR.025),0.05),4)
VaR.df_model.test[5,5] <- round(quantile(na.omit(oos_historical.VaR.025),0.025),4)
VaR.df_model.test[5,6] <- round(quantile(na.omit(oos_historical.VaR.025),0.01),4)

oos_GARCH.VaR.01_min <- lapply(GARCH.daily.marginal, FUN = function(x){min(x[[2]]$`alpha(1%)`)})
oos_GARCH.VaR.01_max <- lapply(GARCH.daily.marginal, FUN = function(x){max(x[[2]]$`alpha(1%)`)})
oos_GARCH.VaR.01_q50 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[2]]$`alpha(1%)`, 0.50)})
oos_GARCH.VaR.01_q05 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[2]]$`alpha(1%)`, 0.05)})
oos_GARCH.VaR.01_q025 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[2]]$`alpha(1%)`, 0.025)})
oos_GARCH.VaR.01_q01 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[2]]$`alpha(1%)`, 0.01)})

VaR.df_model.test[2:4,1] <- as.numeric(do.call(rbind,oos_GARCH.VaR.01_min))
VaR.df_model.test[2:4,2] <- as.numeric(do.call(rbind,oos_GARCH.VaR.01_max))
VaR.df_model.test[2:4,3] <- as.numeric(do.call(rbind,oos_GARCH.VaR.01_q50))
VaR.df_model.test[2:4,4] <- as.numeric(do.call(rbind,oos_GARCH.VaR.01_q05))
VaR.df_model.test[2:4,5] <- as.numeric(do.call(rbind,oos_GARCH.VaR.01_q025))
VaR.df_model.test[2:4,6] <- as.numeric(do.call(rbind,oos_GARCH.VaR.01_q01))

oos_GARCH.VaR.025_min <- lapply(GARCH.daily.marginal, FUN = function(x){min(x[[2]]$`alpha(3%)`)})
oos_GARCH.VaR.025_max <- lapply(GARCH.daily.marginal, FUN = function(x){max(x[[2]]$`alpha(3%)`)})
oos_GARCH.VaR.025_q50 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[2]]$`alpha(3%)`, 0.50)})
oos_GARCH.VaR.025_q05 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[2]]$`alpha(3%)`, 0.05)})
oos_GARCH.VaR.025_q025 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[2]]$`alpha(3%)`, 0.025)})
oos_GARCH.VaR.025_q01 <- lapply(GARCH.daily.marginal, FUN = function(x){quantile(x[[2]]$`alpha(3%)`, 0.01)})

VaR.df_model.test[6:8,1] <- as.numeric(do.call(rbind,oos_GARCH.VaR.025_min))
VaR.df_model.test[6:8,2] <- as.numeric(do.call(rbind,oos_GARCH.VaR.025_max))
VaR.df_model.test[6:8,3] <- as.numeric(do.call(rbind,oos_GARCH.VaR.025_q50))
VaR.df_model.test[6:8,4] <- as.numeric(do.call(rbind,oos_GARCH.VaR.025_q05))
VaR.df_model.test[6:8,5] <- as.numeric(do.call(rbind,oos_GARCH.VaR.025_q025))
VaR.df_model.test[6:8,6] <- as.numeric(do.call(rbind,oos_GARCH.VaR.025_q01))

VaR.df_model.test <- round(VaR.df_model.test,4)

# c-VaR
oos_historical.cVaR.01 <- apply.fromstart(GARCH.daily.marginal[[1]][[1]][,'daily.vol'], FUN = function(x){cVaR(x, 0.01)}, gap = 30)
oos_historical.cVaR.025 <- apply.fromstart(GARCH.daily.marginal[[1]][[1]][,'daily.vol'], FUN = function(x){cVaR(x, 0.025)}, gap = 30)
mean_historical.cVaR.01 <- quantile(na.omit(oos_historical.cVaR.01),0.5)
mean_historical.cVaR.025 <- quantile(na.omit(oos_historical.cVaR.025),0.5)

conditional.excess <- lapply(GARCH.daily.marginal, FUN = function(x){ lapply(x[[3]], FUN = function(x){ x[,'daily.vol'] - x[,2] })})
conditional.excess.VaR.01 <- lapply(conditional.excess, FUN = function(x){ x[[1]] })
conditional.excess.VaR.025 <- lapply(conditional.excess, FUN = function(x){ x[[2]] })
mean.excess.VaR.01 <- do.call(rbind,lapply(conditional.excess.VaR.01, FUN = function(x) {mean(x)}))
mean.excess.VaR.025 <- do.call(rbind,lapply(conditional.excess.VaR.025, FUN = function(x) {mean(x)}))

cVaR.01 <- apply(VaR.df_model.test[2:4,], MARGIN = 2, FUN = function(x){ x + mean.excess.VaR.01[1:3,1] })
rownames(cVaR.01) <- c('iGARCH.cVaR.01','sGARCH.cVaR.01','csGARCH.cVaR.01')

cVaR.025 <- apply(VaR.df_model.test[6:8,], MARGIN = 2, FUN = function(x){ x + mean.excess.VaR.025[1:3,1] })
rownames(cVaR.025) <- c('iGARCH.cVaR.025','sGARCH.cVaR.025','csGARCH.cVaR.025')

VaR.df_model.test <- rbind(VaR.df_model.test,rbind(cVaR.01,cVaR.025))
VaR.df_model.test <- round(VaR.df_model.test,4)
VaR.df_model.test <- VaR.df_model.test[c(-1,-5),]

VaR.df_model.test_trading <- round(apply(VaR.df_model.test[,4:6], MARGIN = 2, FUN = function(x){x*sqrt(10)}),4)
VaR.df_model.test_trading <- as.data.frame(VaR.df_model.test_trading)
colnames(VaR.df_model.test_trading) <- paste0("consumo.trading_",colnames(VaR.df_model.test_trading))
VaR.df_model.test_trading$consumo.sVaR <- round(as.numeric(VaR.df_model.test[,3]*sqrt(250)),4)

# para hacer la serie de tiempo del cVaR es necesario evitar el lookahead bias.
# hay que recalcular (walk forward) e imputar por fecha, con precision
#GARCH.daily.marginal[[1]][[2]]$`alpha(1%)` + mean.excess.VaR.01[[1]]
#lapply(GARCH.daily.marginal, FUN = function(x){ lapply(mean.excess.VaR.01, FUN = function(x){ x[,'daily.vol'] - x[,2] })})

#plot.pdf(GARCH.daily.marginal[[1]][[2]]$`alpha(1%)`)
#plot.pdf(GARCH.daily.marginal[[1]][[2]]$`alpha(3%)`)


# VaR TABLES: VALIDATION (2019-2022) ----

# REGULAR VaR
validation_historical.VaR.01 <- apply.fromstart(GARCH.daily.validation[[1]][[1]][,'daily.vol'], FUN = function(x){quantile(x, 0.01)}, gap = 30)
validation_historical.VaR.025 <- apply.fromstart(GARCH.daily.validation[[1]][[1]][,'daily.vol'], FUN = function(x){quantile(x, 0.025)}, gap = 30)

VaR.df_validation <- matrix(nrow = 8, ncol = 6)
colnames(VaR.df_validation) <- c('min','max','q.50','q.05','q.025','q.01')
rownames(VaR.df_validation) <- c('historical.VaR.01',
                                 'iGARCH.VaR.01',
                                 'sGARCH.VaR.01',
                                 'csGARCH.VaR.01',
                                 'historical.VaR.025',
                                 'iGARCH.VaR.025',
                                 'sGARCH.VaR.025',
                                 'csGARCH.VaR.025')

VaR.df_validation[1,1] <- round(min(na.omit(validation_historical.VaR.01)),4)
VaR.df_validation[1,2] <- round(max(na.omit(validation_historical.VaR.01)),4)
VaR.df_validation[1,3] <- round(quantile(na.omit(validation_historical.VaR.01),0.5),4)
VaR.df_validation[1,4] <- round(quantile(na.omit(validation_historical.VaR.01),0.05),4)
VaR.df_validation[1,5] <- round(quantile(na.omit(validation_historical.VaR.01),0.025),4)
VaR.df_validation[1,6] <- round(quantile(na.omit(validation_historical.VaR.01),0.01),4)

VaR.df_validation[5,1] <- round(min(na.omit(validation_historical.VaR.025)),4)
VaR.df_validation[5,2] <- round(max(na.omit(validation_historical.VaR.025)),4)
VaR.df_validation[5,3] <- round(quantile(na.omit(validation_historical.VaR.025),0.5),4)
VaR.df_validation[5,4] <- round(quantile(na.omit(validation_historical.VaR.025),0.05),4)
VaR.df_validation[5,5] <- round(quantile(na.omit(validation_historical.VaR.025),0.025),4)
VaR.df_validation[5,6] <- round(quantile(na.omit(validation_historical.VaR.025),0.01),4)

validation_GARCH.VaR.01_min <- lapply(GARCH.daily.validation, FUN = function(x){min(x[[2]]$`alpha(1%)`)})
validation_GARCH.VaR.01_max <- lapply(GARCH.daily.validation, FUN = function(x){max(x[[2]]$`alpha(1%)`)})
validation_GARCH.VaR.01_q50 <- lapply(GARCH.daily.validation, FUN = function(x){quantile(x[[2]]$`alpha(1%)`, 0.50)})
validation_GARCH.VaR.01_q05 <- lapply(GARCH.daily.validation, FUN = function(x){quantile(x[[2]]$`alpha(1%)`, 0.05)})
validation_GARCH.VaR.01_q025 <- lapply(GARCH.daily.validation, FUN = function(x){quantile(x[[2]]$`alpha(1%)`, 0.025)})
validation_GARCH.VaR.01_q01 <- lapply(GARCH.daily.validation, FUN = function(x){quantile(x[[2]]$`alpha(1%)`, 0.01)})

VaR.df_validation[2:4,1] <- as.numeric(do.call(rbind,validation_GARCH.VaR.01_min))
VaR.df_validation[2:4,2] <- as.numeric(do.call(rbind,validation_GARCH.VaR.01_max))
VaR.df_validation[2:4,3] <- as.numeric(do.call(rbind,validation_GARCH.VaR.01_q50))
VaR.df_validation[2:4,4] <- as.numeric(do.call(rbind,validation_GARCH.VaR.01_q05))
VaR.df_validation[2:4,5] <- as.numeric(do.call(rbind,validation_GARCH.VaR.01_q025))
VaR.df_validation[2:4,6] <- as.numeric(do.call(rbind,validation_GARCH.VaR.01_q01))

validation_GARCH.VaR.025_min <- lapply(GARCH.daily.validation, FUN = function(x){min(x[[2]]$`alpha(3%)`)})
validation_GARCH.VaR.025_max <- lapply(GARCH.daily.validation, FUN = function(x){max(x[[2]]$`alpha(3%)`)})
validation_GARCH.VaR.025_q50 <- lapply(GARCH.daily.validation, FUN = function(x){quantile(x[[2]]$`alpha(3%)`, 0.50)})
validation_GARCH.VaR.025_q05 <- lapply(GARCH.daily.validation, FUN = function(x){quantile(x[[2]]$`alpha(3%)`, 0.05)})
validation_GARCH.VaR.025_q025 <- lapply(GARCH.daily.validation, FUN = function(x){quantile(x[[2]]$`alpha(3%)`, 0.025)})
validation_GARCH.VaR.025_q01 <- lapply(GARCH.daily.validation, FUN = function(x){quantile(x[[2]]$`alpha(3%)`, 0.01)})

VaR.df_validation[6:8,1] <- as.numeric(do.call(rbind,validation_GARCH.VaR.025_min))
VaR.df_validation[6:8,2] <- as.numeric(do.call(rbind,validation_GARCH.VaR.025_max))
VaR.df_validation[6:8,3] <- as.numeric(do.call(rbind,validation_GARCH.VaR.025_q50))
VaR.df_validation[6:8,4] <- as.numeric(do.call(rbind,validation_GARCH.VaR.025_q05))
VaR.df_validation[6:8,5] <- as.numeric(do.call(rbind,validation_GARCH.VaR.025_q025))
VaR.df_validation[6:8,6] <- as.numeric(do.call(rbind,validation_GARCH.VaR.025_q01))

VaR.df_validation <- round(VaR.df_validation,4)

# c-VaR
validation_historical.cVaR.01 <- apply.fromstart(GARCH.daily.validation[[1]][[1]][,'daily.vol'], FUN = function(x){cVaR(x, 0.01)}, gap = 30)
validation_historical.cVaR.025 <- apply.fromstart(GARCH.daily.validation[[1]][[1]][,'daily.vol'], FUN = function(x){cVaR(x, 0.025)}, gap = 30)
mean_historical.cVaR.01 <- quantile(na.omit(validation_historical.cVaR.01),0.5)
mean_historical.cVaR.025 <- quantile(na.omit(validation_historical.cVaR.025),0.5)

conditional.excess <- lapply(GARCH.daily.validation, FUN = function(x){ lapply(x[[3]], FUN = function(x){ x[,'daily.vol'] - x[,2] })})
conditional.excess.VaR.01 <- lapply(conditional.excess, FUN = function(x){ x[[1]] })
conditional.excess.VaR.025 <- lapply(conditional.excess, FUN = function(x){ x[[2]] })
mean.excess.VaR.01 <- do.call(rbind,lapply(conditional.excess.VaR.01, FUN = function(x) {mean(x)}))
mean.excess.VaR.025 <- do.call(rbind,lapply(conditional.excess.VaR.025, FUN = function(x) {mean(x)}))

cVaR.01 <- apply(VaR.df_validation[2:4,], MARGIN = 2, FUN = function(x){ x + mean.excess.VaR.01[1:3,1] })
rownames(cVaR.01) <- c('iGARCH.cVaR.01','sGARCH.cVaR.01','csGARCH.cVaR.01')

cVaR.025 <- apply(VaR.df_validation[6:8,], MARGIN = 2, FUN = function(x){ x + mean.excess.VaR.025[1:3,1] })
rownames(cVaR.025) <- c('iGARCH.cVaR.025','sGARCH.cVaR.025','csGARCH.cVaR.025')

VaR.df_validation <- rbind(VaR.df_validation,rbind(cVaR.01,cVaR.025))
VaR.df_validation <- round(VaR.df_validation,4)
VaR.df_validation <- VaR.df_validation[c(-1,-5),]

VaR.df_validation_trading <- round(apply(VaR.df_validation[,4:6], MARGIN = 2, FUN = function(x){x*sqrt(10)}),4)
VaR.df_validation_trading <- as.data.frame(VaR.df_validation_trading)
colnames(VaR.df_validation_trading) <- paste0("consumo.trading_",colnames(VaR.df_validation_trading))
VaR.df_validation_trading$consumo.sVaR <- round(as.numeric(VaR.df_validation[,3]*sqrt(250)),4)

# para hacer la serie de tiempo del cVaR es necesario evitar el lookahead bias.
# hay que recalcular (walk forward) e imputar por fecha, con precision
#GARCH.daily.marginal[[1]][[2]]$`alpha(1%)` + mean.excess.VaR.01[[1]]
#lapply(GARCH.daily.marginal, FUN = function(x){ lapply(mean.excess.VaR.01, FUN = function(x){ x[,'daily.vol'] - x[,2] })})

#plot.pdf(GARCH.daily.marginal[[1]][[2]]$`alpha(1%)`)
#plot.pdf(GARCH.daily.marginal[[1]][[2]]$`alpha(3%)`)

# Out-of-Sample 2019-2021 ----
# Plot 1. VaR Only
# FUNCIONES

# sumariza el resultado del backtest. ADAPTADO solo para downside.
# NOTA (refactor): redefinición intencional de rolling.var.backtest (versión
# adaptada solo para downside, 3 bandas). Se mantiene en su posición original.
rolling.var.backtest <- function(rolling.garch.object, symbol, return.series, n){
  
  require(xts)
  require(dplyr)
  require(lubridate)
  
  garch.object.name <- convert.to.char(rolling.garch.object)
  VaR <- paste0(garch.object.name,"@forecast$VaR")
  VaR <- eval(parse(text = VaR))
  realized.var.name <- paste0(garch.object.name,"@forecast$VaR$realized")
  realized.var <- eval(parse(text = realized.var.name))
  
  var050.breach <- which(rolling.garch.object@forecast$VaR$realized < as.numeric(quantile(as.numeric(return.series), 0.050)))
  var025.breach <- which(rolling.garch.object@forecast$VaR$realized < as.numeric(quantile(as.numeric(return.series), 0.025)))
  var010.breach <- which(rolling.garch.object@forecast$VaR$realized < as.numeric(quantile(as.numeric(return.series), 0.010)))
  rows <- index(rolling.garch.object@forecast$VaR[var010.breach,'realized'])
  
  forecast.VaR <- as.data.frame(rolling.garch.object@forecast$VaR)
  forecast.VaR$var01.rank <- ntile(forecast.VaR$`alpha(1%)`,100)
  forecast.VaR$var025.rank <- ntile(forecast.VaR$`alpha(3%)`,100)
  forecast.VaR$var050.rank <- ntile(forecast.VaR$`alpha(5%)`,100)
  
  forecast.VaR$var050.breach <- ifelse(forecast.VaR$realized < rolling.garch.object@forecast$VaR$`alpha(5%)`,1,0)
  forecast.VaR$var025.breach <- ifelse(forecast.VaR$realized < rolling.garch.object@forecast$VaR$`alpha(3%)`,1,0)
  forecast.VaR$var01.breach <- ifelse(forecast.VaR$realized < rolling.garch.object@forecast$VaR$`alpha(1%)`,1,0)
  forecast.VaR$nextday.ret <- lead(forecast.VaR$realized)
  
  var.breaches <- list()
  # breaches are identified in cols: 8-11
  for(i in 1:3){
    breach.days <- forecast.VaR[forecast.VaR[,i+6] == 1,]
    breach.return <- forecast.VaR[forecast.VaR[,i+6] == 1,'nextday.ret']
    breach.days$return <- cumprod(1 + na.omit(breach.return))
    
    forecast.VaR <- as.xts(tail(forecast.VaR,n))
    forecast.density <- as.data.frame(rolling.garch.object@forecast$density)
    forecast.density <- as.xts(tail(forecast.density,n))
    forecast.density$price <- as.numeric(tail(return.series,n))
    
    # see consecutive breaches
    var.breaches[[i]] <- as.data.frame(breach.days)
    var.breaches[[i]]$date <- ymd(rownames(var.breaches[[i]]))
    #var.breaches[[i]]$days.between <- NA
    #for(j in 2:nrow(var.breaches[[i]])){
    #  var.breaches[[i]]$days.between[j] <- var.breaches[[i]]$date[j] - var.breaches[[i]]$date[j-1]
    #}
  }
  
  # upside breaches
  # upside.breaking.VaR <- forecast.VaR[forecast.VaR$realized > 0,]
  # upside.breaking.VaR$var99.breach.levels <- upside.breaking.VaR$realized - upside.breaking.VaR$`alpha(99%)`
  # mean.upside.breaking <- mean(upside.breaking.VaR[upside.breaking.VaR$var99.breach.levels > 0,'var99.breach.levels'])
  
  # downside breaches
  downside.breaking.VaR <- forecast.VaR[forecast.VaR$realized < 0,]
  downside.breaking.VaR$var01.breach.levels <- downside.breaking.VaR$realized - downside.breaking.VaR$`alpha(1%)`
  mean.downside.breaking <- mean(downside.breaking.VaR[downside.breaking.VaR$var01.breach.levels < 0,'var01.breach.levels'])
  
  VaR.backtest <- list(forecast.VaR,
                       var.breaches,
                       #upside.breaking.VaR,
                       downside.breaking.VaR)
  
  return(VaR.backtest)
  
}

i=1
model_spec = ugarchspec(variance.model = list(model = final.models[i,'type_model'],
                                              garchOrder = c(final.models[i,'lag_arch'],
                                                             final.models[i,'lag_garch'])),
                        mean.model = list(armaOrder = c(final.models[i,'lag_ar'],
                                                        final.models[i,'lag_ma'])),
                        distribution = final.models[i,'type_dist'])
symbol = "GOBIXDR"
n=1000
# CALIBRACION CORRIDA Y VALIDACION
# Idea: tal vez el modelo debería ser recalibrado de vez en cuando
# Ejecuta la predicción de 1-día hacía adelante (rolling GARCH) con recalibración cada 30 días
garch.roll = ugarchroll(model_spec, gobix["2014/2021-11-30",'daily.vol'], n.ahead=1,
                        forecast.length = n, solver = "hybrid",
                        refit.every=30, refit.window="moving", VaR.alpha=c(0.01, 0.025, 0.050))

# Prueba de modelo
report(garch.roll, type="VaR")
report(garch.roll, type="fpm")

VaRTest(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )
VaRloss(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )
VaRDurTest(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )

a <- rolling.var.backtest(rolling.garch.object=garch.roll, symbol=symbol, return.series=gobix["2014/2021",'daily.vol'], n=n)

forecast.VaR <- as.xts(a[[1]])
forecast.density <- as.data.frame(garch.roll@forecast$density)
forecast.density <- as.xts(forecast.density)
forecast.density$price <- as.numeric(tail(Cl(gobix["2014/2021-11-30"]),n))
forecast.density$Sigma.annualized <- SMA(forecast.density$Sigma,30) * sqrt(256)
forecast.density$Sigma.simple <- roll_sd(forecast.density$price, width = 60)*sqrt(256)
forecast.density <- na.omit(forecast.density)

#plot.quadrants(data=na.omit(forecast.density), x.var = 'Sigma.simple', y.var = 'gobix.return', description = "Retorno vs. RV")

# see consecutive breaches
forecast.VaR <- a[[1]]
forecast.VaR1 <- as.data.frame(a[[1]][a[[1]]$var025.breach==1,])
forecast.VaR1$date <- ymd(rownames(forecast.VaR1))
forecast.VaR1$days.between <- NA
for(i in 2:nrow(forecast.VaR1)){
  forecast.VaR1$days.between[i] <- forecast.VaR1$date[i] - forecast.VaR1$date[i-1]
}

forecast.VaR.model_test <- forecast.VaR["2019::2021-11-30"]
forecast.density <- forecast.density["2019::2021-11-30"]

# Display the VaR violations
VaR025.violations <- as.data.frame(forecast.VaR.model_test) %>% dplyr::filter(var025.breach == 1)
VaR025.violations$conditional.error <- VaR025.violations$realized - VaR025.violations$`alpha(3%)`
mean(VaR025.violations$conditional.error)
#VaR025.violations$cVaR <- VaR025.violations$`alpha(1%)` + VaR025.violations$conditional.error

VaR01.violations <- as.data.frame(forecast.VaR.model_test) %>% dplyr::filter(var01.breach == 1)
VaR01.violations$conditional.error <- VaR01.violations$realized - VaR01.violations$`alpha(1%)`
mean(VaR01.violations$conditional.error)
#VaR01.violations$cVaR <- VaR01.violations$`alpha(1%)` + VaR01.violations$conditional.error

# Basel 3
GARCH_VaR95.1d <- mean(forecast.VaR.model_test$`alpha(5%`)
GARCH_VaR975.1d <- mean(forecast.VaR.model_test$`alpha(3%`)
GARCH_VaR99.1d <- mean(forecast.VaR.model_test$`alpha(1%`)
GARCH_cVaR975.1d <- mean(forecast.VaR.model_test$`alpha(3%`) + mean(VaR025.violations$conditional.error)
GARCH_cVaR99.1d <- mean(forecast.VaR.model_test$`alpha(1%`) + mean(VaR01.violations$conditional.error)

GARCH_VaR95.10d <- mean(forecast.VaR.model_test$`alpha(5%`)*sqrt(10)
GARCH_VaR975.10d <- mean(forecast.VaR.model_test$`alpha(3%`)*sqrt(10)
GARCH_VaR99.10d <- mean(forecast.VaR.model_test$`alpha(1%`)*sqrt(10)
GARCH_cVaR975.10d <- (mean(forecast.VaR.model_test$`alpha(3%`) + mean(VaR025.violations$conditional.error))*sqrt(10)
GARCH_cVaR99.10d <- (mean(forecast.VaR.model_test$`alpha(1%`) + mean(VaR01.violations$conditional.error))*sqrt(10)

GARCH_VaR95.250d <- mean(forecast.VaR.model_test$`alpha(5%`)*sqrt(250)
GARCH_VaR975.250d <- mean(forecast.VaR.model_test$`alpha(3%`)*sqrt(250)
GARCH_VaR99.250d <- mean(forecast.VaR.model_test$`alpha(1%`)*sqrt(250)
GARCH_cVaR975.250d <- (mean(forecast.VaR.model_test$`alpha(3%`) + mean(VaR025.violations$conditional.error))*sqrt(250)
GARCH_cVaR99.250d <- (mean(forecast.VaR.model_test$`alpha(1%`) + mean(VaR01.violations$conditional.error))*sqrt(250)


df <- matrix(nrow = 5, ncol = 3)
rownames(df) <- c('VaR_95.0',
                  'VaR_97.5',
                  'VaR_99.0',
                  'cVaR_97.5',
                  'cVaR_99.0')

colnames(df) <- c('1d',
                  '10d',
                  '250d')

df[1,1] <- round(GARCH_VaR95.1d,4)
df[1,2] <- round(GARCH_VaR95.10d,4)
df[1,3] <- round(GARCH_VaR95.250d,4)

df[2,1] <- round(GARCH_VaR975.1d,4)
df[2,2] <- round(GARCH_VaR975.10d,4)
df[2,3] <- round(GARCH_VaR975.250d,4)

df[3,1] <- round(GARCH_VaR99.1d,4)
df[3,2] <- round(GARCH_VaR99.10d,4)
df[3,3] <- round(GARCH_VaR99.250d,4)

df[4,1] <- round(GARCH_cVaR975.1d,4)
df[4,2] <- round(GARCH_cVaR975.10d,4)
df[4,3] <- round(GARCH_cVaR975.250d,4)

df[5,1] <- round(GARCH_cVaR99.1d,4)
df[5,2] <- round(GARCH_cVaR99.10d,4)
df[5,3] <- round(GARCH_cVaR99.250d,4)

df

# only downside version
par(mfrow=c(1,1), mar=c(3,3,3,3))
rolling.plot <- plot(forecast.VaR.model_test$realized, ylim=c(min(forecast.VaR.model_test$`alpha(1%` -0.0025), 0), main="Prueba de Modelo (Mercado Normal)", type="h", col="black", lwd=0.50, grid.col = NA)
rolling.plot <- lines(forecast.VaR.model_test$`alpha(5%`, on=1, col="orange", lty=1, lwd=0.75)
rolling.plot <- lines(forecast.VaR.model_test$`alpha(3%`, on=1, col="red4", lty=1, lwd=1)
rolling.plot <- lines(forecast.VaR.model_test$`alpha(1%`, on=1, col="red", lty=1, lwd=1.25)
rolling.plot <- points(forecast.VaR.model_test[forecast.VaR.model_test$var01.breach==1,1], col="red", cex=1.25, pch=16)
rolling.plot <- addLegend("bottomleft", legend.names=c("Mean GARCH VaR:",
                                                       paste("VaR 95.0%:", round(GARCH_VaR95.1d*100,2),"%"),
                                                       paste("VaR 97.5%:", round(GARCH_VaR975.1d*100,2),"%"),
                                                       paste("VaR 99.0%:", round(GARCH_VaR99.1d*100,2),"%")), lty=1, lwd=c(0.75,1,1.25), col=c("white","orange","red4","red"), bty="n", y.intersp = 0.90)

plot(rolling.plot)

par(mfrow=c(1,2), mar=c(5,5,5,5), lend=2)
plot.pdf(data=as.numeric(forecast.VaR.model_test$`alpha(1%`), breaks = 32, freq=FALSE, title = "VaR.990", hist.col = "grey", line.col="red", cex.legend = 1)
plot.pdf(data=as.numeric(forecast.VaR.model_test$`alpha(3%`), breaks = 32, freq=FALSE, title = "VaR.975", hist.col = "grey", line.col="blue", cex.legend = 1)
par(mfrow=c(1,1), mar=c(5,5,5,5))

# Out-of-Sample 2022----
symbol = "GOBIXDR"
n=250
# CALIBRACION CORRIDA Y VALIDACION
# Idea: tal vez el modelo debería ser recalibrado de vez en cuando
# Ejecuta la predicción de 1-día hacía adelante (rolling GARCH) con recalibración cada 30 días

n = 250
# garch.roll = ugarchroll(model_spec, gobix["2014/2021-11-30",'daily.vol'], n.ahead=1,
#                         forecast.length = n, solver = "hybrid",
#                         refit.every=30, refit.window="moving", VaR.alpha=c(0.01, 0.025, 0.050))

garch.roll = ugarchroll(model_spec, oos[,'daily.vol'], n.ahead=1,
                        forecast.length = n, solver = "hybrid",
                        refit.every=30, refit.window="moving", VaR.alpha=c(0.01, 0.025, 0.05))

# Prueba de modelo
report(garch.roll, type="VaR")
report(garch.roll, type="fpm")

VaRTest(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )
VaRloss(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )
VaRDurTest(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )

symbol="GOBIXDR"
a <- rolling.var.backtest(rolling.garch.object=garch.roll, symbol=symbol, return.series=gobix[,'daily.vol'], n=n)

forecast.VaR <- as.xts(a[[1]])
forecast.density <- as.data.frame(garch.roll@forecast$density)
forecast.density <- as.xts(forecast.density)
forecast.density$price <- as.numeric(tail(Cl(oos),n))
forecast.density$Sigma.annualized <- SMA(forecast.density$Sigma,30) * sqrt(256)
forecast.density$Sigma.simple <- roll_sd(forecast.density$price, width = 60)*sqrt(256)
forecast.density <- na.omit(forecast.density)

#plot.quadrants(data=na.omit(forecast.density), x.var = 'Sigma.simple', y.var = 'gobix.return', description = "Retorno vs. RV")

# see consecutive breaches
forecast.VaR <- a[[1]]
forecast.VaR1 <- as.data.frame(a[[1]][a[[1]]$var025.breach==1,])
forecast.VaR1$date <- ymd(rownames(forecast.VaR1))
forecast.VaR1$days.between <- NA
for(i in 2:nrow(forecast.VaR1)){
  forecast.VaR1$days.between[i] <- forecast.VaR1$date[i] - forecast.VaR1$date[i-1]
}

forecast.VaR.validation <- forecast.VaR

# Display the VaR violations
VaR025.violations <- as.data.frame(forecast.VaR.validation) %>% dplyr::filter(var025.breach == 1)
VaR025.violations$conditional.error <- VaR025.violations$realized - VaR025.violations$`alpha(3%)`
mean(VaR025.violations$conditional.error)
#VaR025.violations$cVaR <- VaR025.violations$`alpha(1%)` + VaR025.violations$conditional.error

VaR01.violations <- as.data.frame(forecast.VaR.validation) %>% dplyr::filter(var01.breach == 1)
VaR01.violations$conditional.error <- VaR01.violations$realized - VaR01.violations$`alpha(1%)`
mean(VaR01.violations$conditional.error)
#VaR01.violations$cVaR <- VaR01.violations$`alpha(1%)` + VaR01.violations$conditional.error

# Basel 3
GARCH_VaR95.1d <- mean(forecast.VaR.validation$`alpha(5%`)
GARCH_VaR975.1d <- mean(forecast.VaR.validation$`alpha(3%`)
GARCH_VaR99.1d <- mean(forecast.VaR.validation$`alpha(1%`)
GARCH_cVaR975.1d <- mean(forecast.VaR.validation$`alpha(3%`) + mean(VaR025.violations$conditional.error)
GARCH_cVaR99.1d <- mean(forecast.VaR.validation$`alpha(1%`) + mean(VaR01.violations$conditional.error)

GARCH_VaR95.10d <- mean(forecast.VaR.validation$`alpha(5%`)*sqrt(10)
GARCH_VaR975.10d <- mean(forecast.VaR.validation$`alpha(3%`)*sqrt(10)
GARCH_VaR99.10d <- mean(forecast.VaR.validation$`alpha(1%`)*sqrt(10)
GARCH_cVaR975.10d <- (mean(forecast.VaR.validation$`alpha(3%`) + mean(VaR025.violations$conditional.error))*sqrt(10)
GARCH_cVaR99.10d <- (mean(forecast.VaR.validation$`alpha(1%`) + mean(VaR01.violations$conditional.error))*sqrt(10)

GARCH_VaR95.250d <- mean(forecast.VaR.validation$`alpha(5%`)*sqrt(250)
GARCH_VaR975.250d <- mean(forecast.VaR.validation$`alpha(3%`)*sqrt(250)
GARCH_VaR99.250d <- mean(forecast.VaR.validation$`alpha(1%`)*sqrt(250)
GARCH_cVaR975.250d <- (mean(forecast.VaR.validation$`alpha(3%`) + mean(VaR025.violations$conditional.error))*sqrt(250)
GARCH_cVaR99.250d <- (mean(forecast.VaR.validation$`alpha(1%`) + mean(VaR01.violations$conditional.error))*sqrt(250)

df <- matrix(nrow = 5, ncol = 3)
rownames(df) <- c('VaR_95.0',
                  'VaR_97.5',
                  'VaR_99.0',
                  'cVaR_97.5',
                  'cVaR_99.0')

colnames(df) <- c('1d',
                  '10d',
                  '250d')

df[1,1] <- round(GARCH_VaR95.1d,4)
df[1,2] <- round(GARCH_VaR95.10d,4)
df[1,3] <- round(GARCH_VaR95.250d,4)

df[2,1] <- round(GARCH_VaR975.1d,4)
df[2,2] <- round(GARCH_VaR975.10d,4)
df[2,3] <- round(GARCH_VaR975.250d,4)

df[3,1] <- round(GARCH_VaR99.1d,4)
df[3,2] <- round(GARCH_VaR99.10d,4)
df[3,3] <- round(GARCH_VaR99.250d,4)

df[4,1] <- round(GARCH_cVaR975.1d,4)
df[4,2] <- round(GARCH_cVaR975.10d,4)
df[4,3] <- round(GARCH_cVaR975.250d,4)

df[5,1] <- round(GARCH_cVaR99.1d,4)
df[5,2] <- round(GARCH_cVaR99.10d,4)
df[5,3] <- round(GARCH_cVaR99.250d,4)

df

# only downside version
par(mfrow=c(1,1), mar=c(3,3,3,3))
rolling.plot <- plot(forecast.VaR.validation$realized, ylim=c(min(forecast.VaR.validation$`alpha(1%` -0.0025), 0), main="Validacion de Modelo (Mercado en Estres)", type="h", col="black", lwd=0.50, grid.col = NA)
rolling.plot <- lines(forecast.VaR.validation$`alpha(5%`, on=1, col="orange", lty=1, lwd=0.75)
rolling.plot <- lines(forecast.VaR.validation$`alpha(3%`, on=1, col="red4", lty=1, lwd=1)
rolling.plot <- lines(forecast.VaR.validation$`alpha(1%`, on=1, col="red", lty=1, lwd=1.25)
rolling.plot <- points(forecast.VaR.validation[forecast.VaR.validation$var01.breach==1,1], col="red", cex=1.25, pch=16)
rolling.plot <- addLegend("bottomleft", legend.names=c("Mean GARCH VaR:",
                                                       paste("VaR 95.0%:", round(GARCH_VaR95.1d*100,2),"%"),
                                                       paste("VaR 97.5%:", round(GARCH_VaR975.1d*100,2),"%"),
                                                       paste("VaR 99.0%:", round(GARCH_VaR99.1d*100,2),"%")), lty=1, lwd=c(0.75,1,1.25), col=c("white","orange","red4","red"), bty="n", y.intersp = 0.90)

plot(rolling.plot)

par(mfrow=c(1,2), mar=c(5,5,5,5), lend=2)
plot.pdf(data=as.numeric(forecast.VaR.validation$`alpha(1%`), breaks = 32, freq=FALSE, title = "VaR.990", hist.col = "grey", line.col="red", cex.legend = 1)
plot.pdf(data=as.numeric(forecast.VaR.validation$`alpha(3%`), breaks = 32, freq=FALSE, title = "VaR.975", hist.col = "grey", line.col="blue", cex.legend = 1)
par(mfrow=c(1,1), mar=c(5,5,5,5))

# Out-of-Sample Simulacion 2022 ----
n.sim <- 1000
garch.sim <- matrix(nrow = nrow(oos), ncol=n.sim)
for(i in 1:n.sim){
  p.sim = ugarchsim(my_best_garch,n.sim=nrow(oos), startMethod="sample")
  garch.sim[,i] <- p.sim@simulation$seriesSim
}
garch.sim <- cbind(oos$daily.vol,garch.sim)

garch.sim$Q20 <- NA
garch.sim$Q05 <- NA
garch.sim$Q01 <- NA
garch.sim$Q80 <- NA
garch.sim$Q95 <- NA
garch.sim$Q99 <- NA

for(j in 2:nrow(garch.sim)){
  #upper band
  garch.sim$Q20 <- quantile(garch.sim[j,2:(ncol(garch.sim)-6)],0.80)
  garch.sim$Q05 <- quantile(garch.sim[j,2:(ncol(garch.sim)-6)],0.90)
  garch.sim$Q01 <- quantile(garch.sim[j,2:(ncol(garch.sim)-6)],0.99)
  #lower band
  garch.sim$Q80 <- quantile(garch.sim[j,2:(ncol(garch.sim)-6)],0.20)
  garch.sim$Q95 <- quantile(garch.sim[j,2:(ncol(garch.sim)-6)],0.05)
  garch.sim$Q99 <- quantile(garch.sim[j,2:(ncol(garch.sim)-6)],0.01)
}
garch.sim$mean.Q05 <- mean(na.omit(garch.sim$Q05))
garch.sim$mean.Q01 <- mean(na.omit(garch.sim$Q01))
garch.sim$mean.Q95 <- mean(na.omit(garch.sim$Q95))
garch.sim$mean.Q99 <- mean(na.omit(garch.sim$Q99))

sim.breaches.01 <- which(garch.sim[,1] > garch.sim$mean.Q01)
sim.breaches.05 <- which(garch.sim[,1] > garch.sim$mean.Q05)
sim.breaches.95 <- which(garch.sim[,1] < garch.sim$mean.Q95)
sim.breaches.99 <- which(garch.sim[,1] < garch.sim$mean.Q99)

exp.breaches.95 <- round(nrow(garch.sim)*0.05,0)
real.breaches.95 <- length(sim.breaches.95)
prob.breaches.95 <- real.breaches.95/nrow(garch.sim)

exp.breaches.99 <- round(nrow(garch.sim)*0.01,0)
real.breaches.99 <- length(sim.breaches.99)
prob.breaches.99 <- real.breaches.99/nrow(garch.sim)

sim.plot <- plot(garch.sim[,2:(ncol(garch.sim)-300)], ylim=c(-0.07,0.07),  type="p", cex=0.5, col=adjustcolor("grey", alpha=0.3),
                 main=paste0("Simulated Volatility: GARCH(1,1) + ARMA(0,0)"))
sim.plot <- lines(garch.sim[,1], on=1, col="black")
#sim.plot <- lines(garch.sim[,'Q20'], on=1, col="red", lty=3, lwd=0.5)
#sim.plot <- lines(garch.sim[,'Q05'], on=1, col="red4", lty=2, lwd=1)
sim.plot <- lines(garch.sim[,'Q01'], on=1, col="green4", lty=1, lwd=1)
#sim.plot <- lines(garch.sim[,'Q80'], on=1, col="red", lty=3, lwd=0.5)
#sim.plot <- lines(garch.sim[,'Q95'], on=1, col="red4", lty=2, lwd=0.5)
sim.plot <- lines(garch.sim[,'Q99'], on=1, col="red", lty=1, lwd=1)
sim.plot <- addLegend("topleft", legend.names=c("Realized",
                                                #paste0("20pct.Band.Simulated [",round(tail(garch.sim$Q20,1)*100,2),"%"," , ",round(tail(garch.sim$Q80,1)*100,2),"%","]"),
                                                paste0("5pct.Band.Simulated [",round(tail(garch.sim$mean.Q05,1)*100,2),"%"," , ",round(tail(garch.sim$mean.Q95,1)*100,2),"%","]"), 
                                                paste0("1pct.Band.Simulated [",round(tail(garch.sim$mean.Q01,1)*100,2),"%"," , ",round(tail(garch.sim$mean.Q99,1)*100,2),"%","]")), lty=c(1,2,1), col=c("black","red4","red"), bty="n", y.intersp = 0.75)
sim.plot <- addLegend("bottomleft", legend.names=c(paste0("VaR95.breaches [",round(prob.breaches.95*100,2),"%"," , ", exp.breaches.95,"/", real.breaches.95,"]"), 
                                                   paste0("VaR99.breaches [",round(prob.breaches.99*100,2),"%"," , ", exp.breaches.99,"/", real.breaches.99,"]")), pch=c(1,16), col=c("black","red"), bty="n", y.intersp = 0.75)
#sim.plot <- lines(garch.sim$mean.Q95, on=1, col="red4", lty=2)
#sim.plot <- lines(garch.sim$mean.Q99, on=1, col="red")
sim.plot <- points(garch.sim[sim.breaches.99,1], col="white", bg="red", cex=1, pch=21)
sim.plot <- points(garch.sim[sim.breaches.01,1], col="red", bg="white", cex=1, pch=21)

quantiles <- c(quantile(garch.sim,0.05),
               quantile(garch.sim,0.10),
               quantile(garch.sim,0.90),
               quantile(garch.sim,0.95))

median <- quantile(garch.sim,0.50)

par(mfrow=c(1,1), mar=c(3,3,3,3))
split.screen(rbind(c(0.55, 0.98, 0.08, 0.88), c(0.05, 0.52, 0.1, 0.95)))
screen(2)
plot(sim.plot)
screen(1)
hist(coredata(garch.sim[,300:1000]), breaks=128, xlim=c(-0.05, 0.05), ylim=c(0,180), col="grey", main="GARCH Simulated Vol", freq = FALSE, probability = TRUE, xlab="Daily Var (%)")
mtext("Fitted PDF vs. Realized", side=3)
lines(density(as.numeric(garch.sim[,990:1000]), bw=0.003), lwd=2)
lines(density(as.numeric(garch.sim[,1]), bw=0.003), lwd=2, col="red")
rug(garch.sim[,990:1000], ticksize = 0.03, lwd=1, col= "grey")
rug(as.numeric(quantiles), ticksize = 0.03, lwd=2, col= "red")
rug(as.numeric(median), ticksize = 0.03, lwd=2, col= "blue")
# lines(density(garch.sim, bw=sd(garch.sim)), lwd=2, col="blue")
# legend("topleft", legend=c(paste0("75 pctile: ", round(quantile(garch.sim,0.75),4)*100, "%"),
#                            paste0("90 pctile: ", round(quantile(garch.sim,0.90),4)*100, "%"),
#                            paste0("95 pctile: ", round(quantile(garch.sim,0.95),4)*100, "%"),
#                            paste0("Mean: ", round(mean(garch.sim),4)*100, "%"),
#                            paste0("25 pctile: ", round(quantile(garch.sim,0.25),4)*100, "%"),
#                            paste0("10 pctile: ", round(quantile(garch.sim,0.10),4)*100, "%"),
#                            paste0("5 pctile: ", round(quantile(garch.sim,0.05),4)*100, "%")), cex=1, y.intersp =0.75, bty="n")
legend("topleft", legend=c("GARCH Fit", "Empirical Density"), col=c("black","red"), lty=1, bty="n", xjust = 0, y.intersp = 0.9)
box(col="grey") 
close.screen(all.screens = TRUE)

# Out-of-Sample 2019-2022. FULL OOS ----
# Plot 1. VaR Only

# FUNCIONES
# funcion auxiliar para leer el nombre
# sumariza el resultado del backtest. UPSIDE + DOWNSIDE. REPLACED at the top...
# NOTA (refactor): redefinición intencional de rolling.var.backtest (versión
# upside + downside, 4 bandas). Se mantiene en su posición original.
rolling.var.backtest <- function(rolling.garch.object, symbol, return.series, n){
  
  require(xts)
  require(dplyr)
  require(lubridate)
  
  garch.object.name <- convert.to.char(rolling.garch.object)
  VaR <- paste0(garch.object.name,"@forecast$VaR")
  VaR <- eval(parse(text = VaR))
  #forecast.var <- paste0(garch.object.name,"@forecast$VaR")
  realized.var.name <- paste0(garch.object.name,"@forecast$VaR$realized")
  #assign("realized.var", realized.var.name)
  realized.var <- eval(parse(text = realized.var.name))
  
  var025.breach <- which(rolling.garch.object@forecast$VaR$realized < as.numeric(quantile(as.numeric(return.series), 0.025)))
  var010.breach <- which(rolling.garch.object@forecast$VaR$realized < as.numeric(quantile(as.numeric(return.series), 0.010)))
  var95.breach <- which(rolling.garch.object@forecast$VaR$realized > as.numeric(quantile(as.numeric(return.series), 0.95)))
  var99.breach <- which(rolling.garch.object@forecast$VaR$realized > as.numeric(quantile(as.numeric(return.series), 0.99)))
  rows <- index(rolling.garch.object@forecast$VaR[var010.breach,'realized'])
  
  forecast.VaR <- as.data.frame(rolling.garch.object@forecast$VaR)
  forecast.VaR$var01.rank <- ntile(forecast.VaR$`alpha(1%)`,100)
  forecast.VaR$var025.rank <- ntile(forecast.VaR$`alpha(3%)`,100)
  forecast.VaR$var95.rank <- ntile(forecast.VaR$`alpha(95%)`,100)
  forecast.VaR$var99.rank <- ntile(forecast.VaR$`alpha(99%)`,100)
  forecast.VaR$var025.breach <- ifelse(forecast.VaR$realized < rolling.garch.object@forecast$VaR$`alpha(3%)`,1,0)
  forecast.VaR$var01.breach <- ifelse(forecast.VaR$realized < rolling.garch.object@forecast$VaR$`alpha(1%)`,1,0)
  forecast.VaR$var95.breach <- ifelse(forecast.VaR$realized > rolling.garch.object@forecast$VaR$`alpha(95%)`,1,0)
  forecast.VaR$var99.breach <- ifelse(forecast.VaR$realized > rolling.garch.object@forecast$VaR$`alpha(99%)`,1,0)
  forecast.VaR$nextday.ret <- lead(forecast.VaR$realized)
  
  var.breaches <- list()
  # breaches are identified in cols: 8-11
  for(i in 1:4){
    breach.days <- forecast.VaR[forecast.VaR[,i+9] == 1,]
    breach.return <- forecast.VaR[forecast.VaR[,i+9] == 1,'nextday.ret']
    breach.days$return <- cumprod(1 + na.omit(breach.return))
    
    forecast.VaR <- as.xts(tail(forecast.VaR,n))
    forecast.density <- as.data.frame(rolling.garch.object@forecast$density)
    forecast.density <- as.xts(tail(forecast.density,n))
    forecast.density$price <- as.numeric(tail(return.series,n))
    
    # see consecutive breaches
    var.breaches[[i]] <- as.data.frame(breach.days)
    var.breaches[[i]]$date <- ymd(rownames(var.breaches[[i]]))
    var.breaches[[i]]$days.between <- NA
    for(j in 2:nrow(var.breaches[[i]])){
     var.breaches[[i]]$days.between[j] <- var.breaches[[i]]$date[j] - var.breaches[[i]]$date[j-1]
    }
  }
  
  # upside breaches
  upside.breaking.VaR <- forecast.VaR[forecast.VaR$realized > 0,]
  upside.breaking.VaR$var99.breach.levels <- upside.breaking.VaR$realized - upside.breaking.VaR$`alpha(99%)`
  mean.upside.breaking <- mean(upside.breaking.VaR[upside.breaking.VaR$var99.breach.levels > 0,'var99.breach.levels'])
  
  # downside breaches
  downside.breaking.VaR <- forecast.VaR[forecast.VaR$realized < 0,]
  downside.breaking.VaR$var01.breach.levels <- downside.breaking.VaR$realized - downside.breaking.VaR$`alpha(1%)`
  mean.downside.breaking <- mean(downside.breaking.VaR[downside.breaking.VaR$var01.breach.levels < 0,'var01.breach.levels'])
  
  VaR.backtest <- list(forecast.VaR,
                       var.breaches,
                       upside.breaking.VaR,
                       downside.breaking.VaR)
  
  return(VaR.backtest)
  
}


symbol = "GOBIXDR"
n=1000
# CALIBRACION CORRIDA Y VALIDACION
# Idea: tal vez el modelo debería ser recalibrado de vez en cuando
# Ejecuta la predicción de 1-día hacía adelante (rolling GARCH) con recalibración cada 30 días
garch.roll = ugarchroll(best_spec, gobix["2014/2022",'daily.vol'], n.ahead=1,
                        forecast.length = n, solver = "hybrid",
                        refit.every=30, refit.window="moving", VaR.alpha=c(0.01, 0.025, 0.95, 0.99))

# Prueba de modelo
report(garch.roll, type="VaR")
report(garch.roll, type="fpm")

VaRTest(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )
#VaRloss(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )
VaRDurTest(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )

a <- rolling.var.backtest(rolling.garch.object=garch.roll, symbol=symbol, return.series=gobix["2014/2022",'daily.vol'], n=n)

# static VaR. LT
asset1.ret <- gobix.train[,'daily.vol']
var025.breach <- which(garch.roll@forecast$VaR$realized < as.numeric(quantile(asset1.ret, 0.025)))
var010.breach <- which(garch.roll@forecast$VaR$realized < as.numeric(quantile(asset1.ret, 0.010)))
var95.breach <- which(garch.roll@forecast$VaR$realized > as.numeric(quantile(asset1.ret, 0.95)))
var99.breach <- which(garch.roll@forecast$VaR$realized > as.numeric(quantile(asset1.ret, 0.99)))
rows <- index(garch.roll@forecast$VaR[var025.breach,'realized'])

forecast.VaR <- as.data.frame(garch.roll@forecast$VaR)
forecast.VaR$var025.breach <- ifelse(forecast.VaR$realized < garch.roll@forecast$VaR$`alpha(3%)`,1,0)
forecast.VaR$var010.breach <- ifelse(forecast.VaR$realized < garch.roll@forecast$VaR$`alpha(1%)`,1,0)
# forecast.VaR$var95.breach <- ifelse(forecast.VaR$realized > garch.roll@forecast$VaR$`alpha(95%)`,1,0)
# forecast.VaR$var99.breach <- ifelse(forecast.VaR$realized > garch.roll@forecast$VaR$`alpha(99%)`,1,0)
forecast.VaR$nextday.ret <- lead(forecast.VaR$realized)

breach.days <- forecast.VaR[forecast.VaR$var025.breach == 1,]
breach.return <- forecast.VaR[forecast.VaR$var025.breach == 1,'nextday.ret']
breach.return <- cumprod(1 + na.omit(breach.return))

forecast.VaR <- as.xts(forecast.VaR)
forecast.density <- as.data.frame(garch.roll@forecast$density)
forecast.density <- as.xts(forecast.density)
forecast.density$return <- as.numeric(tail(asset1.ret,n))
#forecast.density$gobix.return <- as.numeric(tail(gobix.train$Close.GOBIX,n))-100
#forecast.density <- merge(forecast.density[,-ncol(forecast.density)], gobix.completo$Close.GOBIX / data.table::shift(gobix.completo$Close.GOBI, n = 30, type = "lag")-1)
forecast.density$Sigma.annualized <- SMA(forecast.density$Sigma,30) * sqrt(256)
forecast.density$Sigma.simple <- roll_sd(forecast.density$return, width = 60)*sqrt(256)
forecast.density <- na.omit(forecast.density)

plot.quadrants(data=na.omit(forecast.density), x.var = 'Sigma.simple', y.var = 'return', description = "Retorno vs. RV")

# see consecutive breaches
forecast.VaR <- a[[1]]
forecast.VaR1 <- as.data.frame(a[[1]][a[[1]]$var025.breach==1,])
forecast.VaR1$date <- ymd(rownames(forecast.VaR1))
forecast.VaR1$days.between <- NA
for(i in 2:nrow(forecast.VaR1)){
  forecast.VaR1$days.between[i] <- forecast.VaR1$date[i] - forecast.VaR1$date[i-1]
}

forecast.VaR.oos <- forecast.VaR["2019/2022"]

# VaR
oos_historical.VaR.025 <- quantile(forecast.VaR.oos$realized, 0.025)
oos_historical.cVaR.025 <- cVaR(forecast.VaR.oos$realized,0.025)
oos_GARCH.cVaR.025 <- cVaR(forecast.VaR.oos$`alpha(3%)`,0.025)

oos_historical.VaR.01 <- quantile(forecast.VaR.oos$realized, 0.010)
oos_historical.cVaR.01 <- cVaR(forecast.VaR.oos$realized,0.01)
oos_GARCH.cVaR.01 <- cVaR(forecast.VaR.oos$`alpha(1%)`,0.01)

forecast.VaR.oos$oos_GARCH.cVaR.025 <- oos_GARCH.cVaR.025
forecast.VaR.oos$oos_GARCH.cVaR.01 <- oos_GARCH.cVaR.01

oos.MR.consumption.historical <- oos_historical.cVaR.025*sqrt(10)
oos.MR.consumption.GARCH <- oos_GARCH.cVaR.025*sqrt(10)
#oos.MR.consumption.sVaR.025 <- oos_GARCH.sVaR.025*sqrt(10)

# Basel 3
#FY Stress
param.99.FY.Basel.sVaR <- -2.32*sd(forecast.VaR.oos$realized)*sqrt(256)
Historical_cVaR.975.Basel_sVaR <- oos_historical.cVaR.025*sqrt(256)
GARCH_cVaR.975.Basel_sVaR <- oos_GARCH.cVaR.025*sqrt(256)
#param.99.FY.sVaR <- -2.32*sd(regime.higher$var.diaria)*sqrt(256)
# regime.higher.sVaR.historical <- regime.higher_historical.cVaR.025*sqrt(256)
# regime.higher.sVaR.GARCH <- regime.higher_GARCH.cVaR.025*sqrt(256)
# 30 day stress
oos.sVaR_30d.historical <- oos_historical.cVaR.025*sqrt(30)
oos.sVaR_30d.GARCH <- oos_GARCH.cVaR.025*sqrt(30)
#baselIII_GARCH.sVaR.025 <- baselIII_GARCH.sVaR.025*sqrt(30)

df <- matrix(nrow = 1, ncol = 4)
colnames(df) <- c('Historical_cVaR.975_10days',
                  'GARCH_cVaR.975_10days',
                  #'GARCH_stress_cVaR.975',
                  #'HS_cVaR.975_30days',
                  #'GARCH_cVaR.975_30days',
                  'Historical_cVaR.975.Basel_sVaR',
                  'GARCH_cVaR.975.Basel_sVaR')
#'Parametric_S-VaR.99_250days'
#'Parametric_Stress_cVaR.975_FY',
#'GARCH_stress_cVaR.99_250days')

df[1,1] <- round(oos.MR.consumption.historical,4)
df[1,2] <- round(oos.MR.consumption.GARCH,4)
#df[1,3] <- round(oos.sVaR_30d.historical,4)
df[1,3] <- round(Historical_cVaR.975.Basel_sVaR,4)
df[1,4] <- round(GARCH_cVaR.975.Basel_sVaR,4)

df.training <- t(df)
df.training

# Display the VaR violations
as.data.frame(forecast.VaR.oos) %>% dplyr::filter(var01.breach == 1)

# Merge synthetic set with rolling forecast 
garch.sim.df <- as.data.frame(garch.sim)
garch.sim.df$date <- lubridate::ymd(rownames(garch.sim.df))
forecast.VaR.oos.df <- as.data.frame(forecast.VaR.oos)
forecast.VaR.oos.df$date <- lubridate::ymd(rownames(forecast.VaR.oos.df))

garch.complete <- merge(garch.sim.df,forecast.VaR.oos.df, by='date')
garch.complete <- as.xts(garch.complete[,-1], order.by = garch.complete$date)

sim.plot <- plot(garch.complete[,2:(ncol(garch.complete)-300)], ylim=c(-0.07,0.07),  type="p", cex=0.5, col=adjustcolor("grey", alpha=0.3),
                 main=paste0("Simulated Volatility: GARCH(1,1) + ARMA(0,0)"))
sim.plot <- lines(garch.complete[,1], on=1, col="black")
sim.plot <- lines(garch.complete$alpha.1.., on=1, col="red2", lty=1, lwd=1)
sim.plot <- lines(garch.complete$alpha.99.., on=1, col="green4", lty=1, lwd=1)
sim.plot <- addLegend("topleft", legend.names=c("Realized",
                                                #paste0("20pct.Band.Simulated [",round(tail(garch.complete$Q20,1)*100,2),"%"," , ",round(tail(garch.complete$Q80,1)*100,2),"%","]"),
                                                paste0("5pct.Band.Simulated [",round(tail(garch.complete$mean.Q05,1)*100,2),"%"," , ",round(tail(garch.complete$mean.Q95,1)*100,2),"%","]"), 
                                                paste0("1pct.Band.Simulated [",round(tail(garch.complete$mean.Q01,1)*100,2),"%"," , ",round(tail(garch.complete$mean.Q99,1)*100,2),"%","]")), lty=c(1,2,1), col=c("black","red4","red"), bty="n", y.intersp = 0.75)
sim.plot <- addLegend("bottomleft", legend.names=c(paste0("VaR95.breaches [",round(prob.breaches.95*100,2),"%"," , ", exp.breaches.95,"/", real.breaches.95,"]"), 
                                                   paste0("VaR99.breaches [",round(prob.breaches.99*100,2),"%"," , ", exp.breaches.99,"/", real.breaches.99,"]")), pch=c(1,16), col=c("black","red"), bty="n", y.intersp = 0.75)
#sim.plot <- lines(garch.complete$mean.Q95, on=1, col="red4", lty=2)
#sim.plot <- lines(garch.complete$mean.Q99, on=1, col="red")
sim.plot <- points(garch.complete[garch.complete$var01.breach == 1,1], col="white", bg="red", cex=1, pch=21)
sim.plot <- points(garch.complete[garch.complete$var99.breach == 1,1], col="red", bg="white", cex=1, pch=21)
sim.plot

# REGIME FILTERED VaR ----

symbol = "GOBIXDR"
n=250
# CALIBRACION CORRIDA Y VALIDACION
# Idea: tal vez el modelo debería ser recalibrado de vez en cuando
# Ejecuta la predicción de 1-día hacía adelante (rolling GARCH) con recalibración cada 30 días
garch.roll.filtered = ugarchroll(best_spec, gobix["2021-12-01/2022",'daily.vol'], n.ahead=1,
                        forecast.length = n, solver = "hybrid",
                        refit.every=30, refit.window="moving", VaR.alpha=c(0.01, 0.025))

# Prueba de modelo
report(garch.roll.filtered, type="VaR")
report(garch.roll.filtered, type="fpm")

#VaRTest(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )
#VaRloss(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )
#VaRDurTest(alpha=0.01, actual = garch.roll@forecast$VaR$realized, VaR = garch.roll@forecast$VaR$`alpha(1%)` )

a <- rolling.var.backtest(rolling.garch.object=garch.roll.filtered, symbol=symbol, return.series=gobix["2021-12-01/2022",'daily.vol'], n=n)

# static VaR. LT
asset1.ret <- tail(baselIII.12m.stress[,'daily.vol'],n)
var025.breach <- which(garch.roll.filtered@forecast$VaR$realized < as.numeric(quantile(asset1.ret, 0.025)))
var010.breach <- which(garch.roll.filtered@forecast$VaR$realized < as.numeric(quantile(asset1.ret, 0.010)))
# var95.breach <- which(garch.roll.filtered@forecast$VaR$realized > as.numeric(quantile(asset1.ret, 0.95)))
# var99.breach <- which(garch.roll.filtered@forecast$VaR$realized > as.numeric(quantile(asset1.ret, 0.99)))
rows <- index(garch.roll.filtered@forecast$VaR[var025.breach,'realized'])

forecast.VaR <- as.data.frame(garch.roll.filtered@forecast$VaR)
forecast.VaR$var025.breach <- ifelse(forecast.VaR$realized < garch.roll.filtered@forecast$VaR$`alpha(3%)`,1,0)
forecast.VaR$var010.breach <- ifelse(forecast.VaR$realized < garch.roll.filtered@forecast$VaR$`alpha(1%)`,1,0)
# forecast.VaR$var95.breach <- ifelse(forecast.VaR$realized > garch.roll.filtered@forecast$VaR$`alpha(95%)`,1,0)
# forecast.VaR$var99.breach <- ifelse(forecast.VaR$realized > garch.roll.filtered@forecast$VaR$`alpha(99%)`,1,0)
forecast.VaR$nextday.ret <- lead(forecast.VaR$realized)

breach.days <- forecast.VaR[forecast.VaR$var025.breach == 1,]
breach.return <- forecast.VaR[forecast.VaR$var025.breach == 1,'nextday.ret']
breach.return <- cumprod(1 + na.omit(breach.return))

forecast.VaR <- as.xts(forecast.VaR)
forecast.density <- as.data.frame(garch.roll.filtered@forecast$density)
forecast.density <- as.xts(forecast.density)
#forecast.density$price <- as.numeric(tail(asset1.ret,n))

# see consecutive breaches
forecast.VaR <- a[[1]]
forecast.VaR1 <- as.data.frame(a[[1]][a[[1]]$var025.breach==1,])
forecast.VaR1$date <- ymd(rownames(forecast.VaR1))
forecast.VaR1$days.between <- NA
for(i in 2:nrow(forecast.VaR1)){
  forecast.VaR1$days.between[i] <- forecast.VaR1$date[i] - forecast.VaR1$date[i-1]
}

# NOTA (refactor): se eliminó aquí un bloque comentado de alternativas
# abandonadas (construcción de VaR.higher vía merge de z/y).


VaR.higher <- tail(forecast.VaR,250)
names(VaR.higher)[1:2] <- c('alpha.1..','alpha.3..')
baselIII.filtered.stress.VaR <- VaR.higher["2021-11-30::2022-06-30"]
#names(baselIII.filtered.stress.VaR)[1:2] <- c('alpha.1..','alpha.3..')

# VaR
sample.1 = nrow(VaR.higher)/nrow(forecast.VaR)
regime.higher_historical.VaR.025 <- quantile(VaR.higher$realized, 0.025)
regime.higher_historical.cVaR.025 <- cVaR(VaR.higher$realized,0.025)
regime.higher_GARCH.cVaR.025 <- cVaR(VaR.higher$alpha.3..,0.025)
baselIII_GARCH.sVaR.025 <- cVaR(baselIII.filtered.stress.VaR$alpha.3..,0.025)

regime.higher_historical.VaR.01 <- quantile(VaR.higher$realized, 0.010)
regime.higher_historical.cVaR.01 <- cVaR(VaR.higher$realized,0.01)
regime.higher_GARCH.cVaR.01 <- cVaR(VaR.higher$alpha.1..,0.01)
baselIII_GARCH.sVaR.01 <- cVaR(baselIII.filtered.stress.VaR$alpha.1..,0.01)


# Capital charge
# Basel 1.5
regime.higher.MR.consumption.historical <- regime.higher_historical.cVaR.025*sqrt(10)
regime.higher.MR.consumption.GARCH <- regime.higher_GARCH.cVaR.025*sqrt(10)
baselIII.MR.consumption.sVaR.025 <- baselIII_GARCH.sVaR.025*sqrt(10)

# Basel 3
#FY Stress
param.99.FY.Basel.sVaR <- -2.32*sd(baselIII.12m.stress$var.diaria)*sqrt(256)
param.99.FY.sVaR <- -2.32*sd(regime.higher$var.diaria)*sqrt(256)
# regime.higher.sVaR.historical <- regime.higher_historical.cVaR.025*sqrt(256)
# regime.higher.sVaR.GARCH <- regime.higher_GARCH.cVaR.025*sqrt(256)
# 30 day stress
regime.higher.sVaR_30d.historical <- regime.higher_historical.cVaR.025*sqrt(30)
regime.higher.sVaR_30d.GARCH <- regime.higher_GARCH.cVaR.025*sqrt(30)
baselIII_GARCH.sVaR.025 <- baselIII_GARCH.sVaR.025*sqrt(30)

df <- matrix(nrow = 1, ncol = 5)
colnames(df) <- c('HS_cVaR.975',
                  'GARCH_cVaR.975',
                  'GARCH_stress_cVaR.975',
                  #'HS_cVaR.975_30d.sVaR',
                  #'GARCH_cVaR.975_30d.sVaR',
                  'Parametric_Stress_cVaR.975_FY',
                  'GARCH_stress_cVaR.99_FY')

df[1,1] <- round(regime.higher.MR.consumption.historical,4)
df[1,2] <- round(regime.higher.MR.consumption.GARCH,4)
df[1,3] <- round(baselIII.MR.consumption.sVaR.025,4)
df[1,4] <- round(param.99.FY.sVaR,4)
df[1,5] <- round(param.99.FY.Basel.sVaR,4)

df <- t(df)
df

# Plot modeled vol

baselIII.filtered.stress.VaR$HS_cVaR.975 <- regime.higher_historical.cVaR.025
baselIII.filtered.stress.VaR$GARCH_cVaR.975 <- regime.higher_GARCH.cVaR.025
baselIII.filtered.stress.VaR$GARCH_stress_cVaR.975 <- baselIII_GARCH.sVaR.025

par(mfrow=c(1,1), mar=c(3,3,3,3))
rolling.plot <- plot(baselIII.filtered.stress.VaR$realized, ylim=c(min(baselIII.filtered.stress.VaR$alpha.1.. -0.001), max(0)), main="Mercado Turbulento", type="h", col="black", lwd=0.50, grid.col = NA)
rolling.plot <- lines(baselIII.filtered.stress.VaR$alpha.1.., on=1, col="orange", lty=1, lwd=1)
rolling.plot <- points(baselIII.filtered.stress.VaR[baselIII.filtered.stress.VaR$var01.breach==1,1], col="red", cex=1.25, pch=16)

rolling.plot <- lines(baselIII.filtered.stress.VaR$HS_cVaR.975, on=1, col="red", lty=2, lwd=1.5)
rolling.plot <- lines(baselIII.filtered.stress.VaR$GARCH_cVaR.975, on=1, col="red4", lty=2, lwd=1.5)
#rolling.plot <- lines(baselIII.filtered.stress.VaR$GARCH_stress_cVaR.975, on=1, col="red4", lty=1, lwd=1.5)

rolling.plot <- addLegend("bottomleft", legend.names=c(paste("GARCH stress-cVaR 0.025:", round(last(baselIII.filtered.stress.VaR$alpha.3..)*100,2),"%"),
                                                       paste("GARCH stress-cVaR 0.010:", round(last(baselIII.filtered.stress.VaR$alpha.1..)*100,2),"%")), lty=c(1,1), col=c("orange","red"), bty="n", y.intersp = 0.75)
plot(rolling.plot)

# NOTA (refactor): se eliminó aquí un bloque comentado "Plot 2. VaR + Sigma + Mu.
# Pendiente ajustar" (gráfico abandonado).

forecast.density.plot <- plot(forecast.density[,2]*sqrt(256), main="Sigma", col="red", lwd=2, grid.col=NA)
forecast.density.plot <- lines(abs(forecast.density[,6]), on=NA, col="grey", lwd=1)
#forecast.density.plot <- lines(forecast.density[,1], on=NA, main="Mu", col="black", lwd=1)
plot(forecast.density.plot)

par(mfrow=c(1,2), mar=c(5,5,5,5))
plot.pdf(data=as.numeric(forecast.VaR$`alpha(1%)`), breaks = 32, freq=FALSE, title = "VaR.99", hist.col = "grey", line.col="red", cex.legend = 1)
plot.pdf(data=as.numeric(forecast.VaR$`alpha(3%)`), breaks = 32, freq=FALSE, title = "VaR.975", hist.col = "grey", line.col="blue", cex.legend = 1)
par(mfrow=c(1,1), mar=c(5,5,5,5))
#boxplot(forecast.VaR$`alpha(1%)`)

# TABLAS ----

# Tablas 4, 5, 6
kable(round(performance.2022-1,3), caption="Retorno cumulativo desde el máximo histórico del GOBIX (septiembre 9, 2021)")
kable(correlation.tb.daily, caption="Correlación Pearson (retornos diarios)") 
kable(risk_performance, caption="Performance y riesgo histórico - benchmark US Treasury 10A (IEF)")

# Tablas 7, 8, 9
knitr::kable(table.autocorrelation, caption="Tabla autocorrelación diaria")
knitr::kable(table.LB, caption="Estadística Ljung-Box")
knitr::kable(round(table.arch,4), caption="Tabla prueba ARCH (LM-stat)")

# Tabla 13
knitr::kable(round(merged.df,4), caption="Resumen pruebas relevantes")

# Resumen ejercicio fitting ARMA, ARCH
knitr::kable(best_model_freq.resume, caption="Distribución configuración Top 1% modelos")
knitr::kable(garch_type.resume, caption="Distribución modelo GARCH")
knitr::kable(dist_type.resume, caption="Distribución marginal término error")

# Tabla 23
knitr::kable(sigma.df, caption="Estimación ITS Sigma: Estadísticas Descriptivas")

# Residual
knitr::kable(residual.df, caption="Estimación ITS Residual: Estadísticas Descriptivas")


# Tabla 24: Prueba de modelo, validación
knitr::kable(oos.df.resumed, caption="Configuraciones GARCH")


# NOTA (refactor): las siguientes líneas usan `out` y `best_models`, objetos que
# solo existen si se ejecutó la calibración completa (RUN_FULL_CALIBRATION <- TRUE),
# igual que en el script original.
# Top 10 GARCH models
models.top20 <- na.omit(out$tab_out)
models.top20 <- models.top20 %>% dplyr::arrange(BIC)
models.top20 <- head(as.data.frame(models.top20[,1:8]),20)
knitr::kable(models.top20, caption="Mejor modelo")

# GARCH coeficientes y params
garch.params <- as.data.frame(t(round(coef(my_best_garch),3)))
rownames(garch.params) <- ''
knitr::kable(garch.params, caption=paste(best_models[1],"coeficientes"))

infocriteria(my_best_garch)
likelihood(my_best_garch)
nyblom(my_best_garch)
signbias(my_best_garch)
gof(my_best_garch,c(20,30,40,50))

ni <- newsimpact(z = NULL, my_best_garch)
par(mfrow=c(1,1), mar=c(5,5,5,5))
plot(ni$zx, ni$zy, ylab=ni$yexpr, xlab=ni$xexpr, type="l", main = "News Impact Curve")

#plot(sqrt(my_best_garch@fit$var)*sqrt(250), type="l")
plot(garch.roll@forecast$density$Sigma*sqrt(250), type="l")
lines(roll_mean(garch.roll@forecast$density$Sigma*sqrt(250), width = 60), col="red", lwd=2)



garch.criteria <- as.data.frame(round(cbind(uncmean(my_best_garch),
                                            uncvariance(my_best_garch),
                                            persistence(my_best_garch),
                                            halflife(my_best_garch)),3))

colnames(garch.criteria) <- c('media.incondicional',
                              'varianza.incondicional',
                              'persistencia',
                              'half-life')

garch.criteria <- t(garch.criteria)
colnames(garch.criteria) <- 'coeficiente'
knitr::kable(garch.criteria, caption="")



# Gráficos utilizados ----
## 001. GOBIXDR drawdown ----
par(mfrow=c(3,1), mar=c(0,4,2,4))
plot(gobix.completo$Close.GOBIX~gobix.completo$date, type="l", xaxt="n", ylab="", xlab="", main="GOBIXDR", cex.main=1.25, cex.sub=1, cex.lab=1, cex.axis=1)
legend("topleft", legend="Enero, 1 2014 = 100", cex=1, bty = "n")
box(col="grey")
#title(adj = 1, line = 2)
box(col="grey")

par(mar=c(0,4,0,4))
plot(gobix.completo$drawdown~gobix.completo$date, type="h", xaxt="n", col="grey", ylab="", xlab="", cex.main=1, cex.sub=1, cex.lab=1, cex.axis=1)
legend("bottomleft", legend="Pérdida desde el nivel máximo", cex=1, bty = "n")
box(col="grey")

par(mar=c(2,4,0,4))
plot(gobix.completo$rolling.sd~gobix.completo$date, type="h", ylab="", xlab="", cex.main=1, cex.sub=1, cex.lab=1, cex.axis=1)
legend("topleft", legend="Desviación Estándar (media móvil de 3 meses)", cex=1, bty = "n")
box(col="grey")
par(mfrow=c(1,1), mar=c(2,2,2,2))

## 002. RV: GOBIXDR vs. benchmarks comparables ----
par(mfrow=c(6,1), mar=c(0,4,4,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`EMB: JPM EMBI`, 
     type="l", main="Variación diaria alternativas de Renta Fija", xaxt="n", xlab="", ylab="")
legend("topleft", legend="EMB: JPM EMBI", cex=1, bty = "n")
mtext("enero 2014 - diciembre 2022", cex=0.90, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`HYG: US Especulativos `, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="HYG: US Especulativos", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`AGG: US Corporativos Grado Inversion`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="AGG: US Bonos Corporativos Grado Inversion", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`BND: Bloomberg US Aggreate`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="BND: Bloomberg US Aggreate", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`VanEck JPM EM Local Currency Bond`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="VanEck JPM EM Local Currency Bond", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`SPDR Bloomberg EM Local Bond`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="SPDR Bloomberg EM Local Bond", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`iShares JPM EM Local Currency`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="iShares JPM EM Local Currency", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`iShares JPM EM High Yield Bond`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="iShares JPM High Yield Bond", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`WisdomTree EM Local Debt Fund`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="WisdomTree EM Local Debt Fund", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`Vanguard EM Government Bond`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="Vanguard EM Government Bond", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`DoubleLine EM Fixed Income`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="DoubleLine EM Fixed Income", cex=1, bty = "n")
box(col="grey")

par(mar=c(4,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$GOBIXDR, 
     type="l", xlab="", ylab="")
legend("topleft", legend="GOBIXDR", cex=1, bty = "n")
box(col="grey")
par(mfrow=c(1,1), mar=c(4,4,4,4))
## 003. Distribución Probabilidad Marginal GOBIXDR ----
par(mfrow=c(1,1), mar=c(2,2,2,2))
hist(gobix[,'daily.vol',drop=FALSE], breaks=64, col="lightblue", main="Distribución Marginal", cex.main=0.75, freq=FALSE, xlab="", ylab="", axes = FALSE, probability = TRUE)
mtext("Retornos diarios", side=3, cex=0.65)
# abline(v=0, col="black", lwd=0.75)
axis(1,cex.axis=0.60)
axis(2,cex.axis=0.60)
rug(gobix[,'daily.vol',drop=FALSE], ticksize = 0.04, lwd=1, col= "grey")
# rug(as.numeric(quantiles), ticksize = 0.04, lwd=2, col= "red")
# rug(as.numeric(median), ticksize = 0.04, lwd=6, col= line.col)
# lines(density(gobix[,'daily.vol',drop=FALSE], bw=sd(gobix[,'daily.vol',drop=FALSE])), lwd=2, col="blue")
box(col="grey") 
## 004. QQplot ----
par(mfrow=c(1,3), mar=c(4,4,4,4))
chart.QQPlot(gobix["2014/2018",'daily.vol',drop=FALSE], 
             main = "Entrenamiento",
             xlab = "Cuartíles normales",
             ylab = "Cuartíles empíricos",
             line=c("robust"), 
             envelope=0.95)
mtext("(2014-2018)", cex = 0.7)

chart.QQPlot(gobix["2019/2021",'daily.vol',drop=FALSE], 
             main = "Prueba de modelo",
             xlab = "Cuartíles normales",
             ylab = "Cuartíles empíricos",
             line=c("robust"), 
             envelope=0.95)
mtext("(2019-2021)", cex = 0.7)

chart.QQPlot(gobix["2022",'daily.vol',drop=FALSE], 
             main = "Validación",
             xlab = "Cuartíles normales",
             ylab = "Cuartíles empíricos",
             line=c("robust"), 
             envelope=0.95)
mtext("(2022)", cex = 0.7)
par(mfrow=c(1,1))

## 005. Sigma todos los períodos ----
par(mfrow=c(1,1), mar=c(2,2,1,2))
plot(df.sigma$iGARCH~df.sigma$date, ylim=c(0,0.02), type="l", lwd=0.75, ylab="", main="Estimación sigma", col="white", cex.main=0.80, cex.sub=0.75, cex.lab=0.75, cex.axis=0.75)
rect(xleft=ymd("2019-01-03"), xright=ymd("2021-12-30"), ybottom=-0.0015, ytop=0.0215, col="white", border = "red", lty=2, lwd=1)
text(x = ymd("2020-06-30"), y = 0.0205, label = "Prueba de modelo", pos = 1, adj = 0, cex = 0.8)
rect(xleft=ymd("2022-01-03"), xright=ymd("2022-12-30"), ybottom=-0.0015, ytop=0.0215, col="white", border = "red", lty=2, lwd=1)
text(x = ymd("2022-06-30"), y = 0.0205, label = "Validación", pos = 1, adj = 0, cex = 0.8)
lines(df.sigma$iGARCH~df.sigma$date, col="green", lwd=0.5)
lines(df.sigma$sGARCH~df.sigma$date, col="red", lwd=0.5)
lines(df.sigma$csGARCH~df.sigma$date, col="blue", lwd=0.5)
legend("topleft", lty=1, lwd=1, legend=c("iGARCH","sGARCH","csGARCH"), horiz=TRUE, col=c("green","red","blue"), bty = "n", cex=0.80)
box(col="grey")
par(mfrow=c(1,1), mar=c(4,4,4,4))

## 009. ACF obs, ACF resids----
par(mfrow=c(2,2), mar=c(4,4,4,4))
plot(training.best_garch[[1]], which=4)
plot(training.best_garch[[1]], which=5)
plot(training.best_garch[[1]], which=6)
plot(training.best_garch[[1]], which=7)
par(mfrow=c(1,1))

par(mfrow=c(3,2), mar=c(4,4,4,4))
plot(training.best_garch[[1]], which=10)
plot(training.best_garch[[1]], which=11)
plot(training.best_garch[[2]], which=10)
plot(training.best_garch[[2]], which=11)
plot(training.best_garch[[3]], which=10)
plot(training.best_garch[[3]], which=11)
par(mfrow=c(1,1))
