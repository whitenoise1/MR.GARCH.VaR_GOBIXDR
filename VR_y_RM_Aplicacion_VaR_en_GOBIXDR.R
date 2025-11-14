# title: "RIESGO DE MERCADO Y VOLATILIDAD REALIZADA: Aplicacion del Value-at-Risk (VaR) en el indice GOBIXDR"
# author: Stefan Bolta, FRM.
# affiliation: Superintendencia de Bancos
# creado: 17/09/2023
# modificado: 18/07/2024
# subtitle: Superintendencia de Bancos, Departamento de Estudios Económicos

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

# Helper functions ----
# GARCH model fitting functions
# Perlin, M. S., Mastella, M., Vancin, D. F., & Ramos, H. P. (2021). A GARCH tutorial with R. Revista de Administração Contemporânea, 25(1), e200088. https://doi.org/10.1590/1982-7849rac2021200088
# https://github.com/msperlin/GARCH-RAC
source("GARCH_msperlin_functions.r")

cVaR <- function(x, conf){
  threshold <- quantile(na.omit(x),conf)
  excesos.VaR.historico <- which(x < threshold)
  cVaR.historico <- mean(x[excesos.VaR.historico])
  return(cVaR.historico)
}

plot.quadrants <- function(data, x.var, y.var, description){
  
  require(scales)
  
  asset.data <- na.omit(data)
  
  a <- lm(coredata(asset.data[,y.var])~coredata(asset.data[,x.var]))
  spline <- smooth.spline(x=coredata(asset.data[,x.var]), y=coredata(asset.data[,y.var]), df=4, lambda=0.03)
  
  #par(mfrow=c(1,1), mar=c(4,4,4,4), xpd=FALSE) 
  par(xpd=FALSE) 
  plot(coredata(asset.data[,x.var]), coredata(asset.data[,y.var]), 
       main=description, ylab=y.var, xlab=x.var, pch=21, cex=1, bg="grey80", col="blue", lwd=1.5)
  #mtext(description, side=3)
  #points(as.numeric(tail(asset.data[,x.var],1)), as.numeric(tail(asset.data[,y.var],1)), col="red", pch=16, cex=1.25)
  abline(v=mean(coredata(asset.data[,x.var])), lty=2, col="grey", lwd=0.75)
  abline(h=mean(coredata(asset.data[,y.var])), lty=2, col="grey", lwd=0.75)
  #text(x=min(asset.data[,x.var]), y=mean(asset.data[,y.var])-sd(asset.data[,y.var])/3, labels="Quadrant III", pos=4, cex=0.8, col="red")
  #text(x=max(asset.data[,x.var]), y=mean(asset.data[,y.var])+sd(asset.data[,y.var])/3, labels="Quadrant I", pos=2, cex=0.8, col="green")
  #text(x=max(asset.data[,x.var]), y=mean(asset.data[,y.var])-sd(asset.data[,y.var])/3, labels="Quadrant IV", pos=2, cex=0.8, col="grey40")
  #text(x=min(asset.data[,x.var]), y=mean(asset.data[,y.var])+sd(asset.data[,y.var])/3, labels="Quadrant II", pos=4, cex=0.8, col="grey40")
  abline(a, col="red", lwd=2)
  lines(spline, col="purple", lwd=1, lty=2)
  legend("topleft", legend=c(paste("y =", round(summary(a)$coefficients[1,1],3),"+",round(summary(a)$coefficients[2,1],3),"x"),
                             paste("y.t-value:", round(summary(a)$coefficients[1,3],3)),
                             paste("x.std error:", round(summary(a)$coefficients[2,2],3)),
                             paste("x.t-value:", round(summary(a)$coefficients[2,3],3)),
                             paste("r.squared:", round(summary(a)$r.squared,3))), bty="n")
  box(col = "grey")
}

prep.corr.matrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

return.q <- function(x){
  return(tail(PerformanceAnalytics::table.Autocorrelation(x),1))
}

convert.to.char <- function(v1) { deparse(substitute(v1)) }

do_arch_test <- function(x, max_lag = 5) {
  
  require(FinTS)
  require(tidyverse)
  
  do_single_arch <- function(x, used_lag)  {
    test_out <- FinTS::ArchTest(x, lags = used_lag)
    
    res_out <- tibble(Lag = used_lag,
                      `LMStatistic` = test_out$statistic, 
                      `pvalue` = test_out$p.value)
  }
  
  tab_out <- bind_rows(map(1:max_lag,.f = do_single_arch, x = x))
  
  return(tab_out)
}

plot.pdf <- function(data, breaks, freq, hist.col, line.col, title, cex.legend){
  
  data <- as.numeric(unlist(data))

  # PLOT3: PDF histogram
  quantiles <- c(quantile(data,0.05),
                 quantile(data,0.10),
                 quantile(data,0.90),
                 quantile(data,0.95))
  
  median <- quantile(data,0.50)
  
  hist(data, breaks=breaks, col=hist.col, main=title, freq=freq, xlab="", ylab="", axes = FALSE, border = "white")
  abline(v=0, col="black", lwd=0.75)
  axis(1,cex.axis=0.9)
  axis(2,cex.axis=0.9)
  rug(data, ticksize = 0.04, lwd=1, col= "grey")
  rug(as.numeric(quantiles), ticksize = 0.04, lwd=2, col= "red")
  rug(as.numeric(median), ticksize = 0.04, lwd=6, col= line.col)
  lines(density(data, bw=sd(data)), lwd=2, col=line.col)
  legend("topleft", legend=c(paste0("75 pctile: ", round(quantile(data,0.75),4)*100, "%"),
                             paste0("90 pctile: ", round(quantile(data,0.90),4)*100, "%"),
                             paste0("95 pctile: ", round(quantile(data,0.95),4)*100, "%"),
                             paste0("Mean: ", round(mean(data),4)*100, "%"),
                             paste0("25 pctile: ", round(quantile(data,0.25),4)*100, "%"),
                             paste0("10 pctile: ", round(quantile(data,0.10),4)*100, "%"),
                             paste0("5 pctile: ", round(quantile(data,0.05),4)*100, "%")), cex=cex.legend, y.intersp =0.75, x.intersp=0, bty="n", xjust=0, yjust=0)
  box(col="grey") 
}

rolling.var.backtest <- function(rolling.garch.object, symbol, return.series, n){
  
  require(xts)
  require(dplyr)
  require(lubridate)
  
  garch.object.name <- convert.to.char(rolling.garch.object)
  VaR <- paste0(garch.object.name,"@forecast$VaR")
  VaR <- eval(parse(text = VaR))
  realized.var.name <- paste0(garch.object.name,"@forecast$VaR$realized")
  realized.var <- eval(parse(text = realized.var.name))
  
  var010.breach <- which(rolling.garch.object@forecast$VaR$realized < as.numeric(quantile(as.numeric(return.series), 0.010)))
  var99.breach <- which(rolling.garch.object@forecast$VaR$realized > as.numeric(quantile(as.numeric(return.series), 0.99)))
  rows <- index(rolling.garch.object@forecast$VaR[var010.breach,'realized'])
  
  forecast.VaR <- as.data.frame(rolling.garch.object@forecast$VaR)
  forecast.VaR$var01.rank <- ntile(forecast.VaR$`alpha(1%)`,100)
  forecast.VaR$var99.rank <- ntile(forecast.VaR$`alpha(99%)`,100)
  forecast.VaR$var01.breach <- ifelse(forecast.VaR$realized < rolling.garch.object@forecast$VaR$`alpha(1%)`,1,0)
  forecast.VaR$var99.breach <- ifelse(forecast.VaR$realized > rolling.garch.object@forecast$VaR$`alpha(99%)`,1,0)
  forecast.VaR$nextday.ret <- lead(forecast.VaR$realized)
  
  var.breaches <- list()
  # i = number of bands
  # N is the number of rows to get to the breach variable. +5 if pct1/99.
  # breaches are identified in cols: 8-11.
  for(i in 1:2){
    breach.days <- forecast.VaR[forecast.VaR[,i+5] == 1,]
    breach.return <- forecast.VaR[forecast.VaR[,i+5] == 1,'nextday.ret']
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

get.gobix <- function(start.date, end.date){
  require(xts)
  require(PerformanceAnalytics)
  require(roll)
  
  gobix <- read.csv(url("https://www.bvrd.com.do/indice/Data/GobixDataIRP.csv")) # descarga
  gobix <- gobix[,1:2] # se queda con dos primeras columnas
  colnames(gobix)[1] <- "fecha" # renombra la columna 1
  gobix$fecha <- lubridate::mdy(gobix[,1]) # la convierte en formato fechas mes-dia-año
  gobix <- xts::as.xts(gobix[,-1], order.by=gobix$fecha) # convierte el formato en XTS
  colnames(gobix)[1] <- "Close.GOBIX" # renombra la serie (indice al cierre del dia)
  gobix <- na.omit(gobix)
  
  range <- paste0(start.date,"/",end.date)
  gobix.completo <- gobix[range]
  
  return(gobix.completo)
  
}

garch.oos.sim <- function(model.df, 
                          model.description, 
                          training.data, 
                          method, 
                          confidence,
                          band.colors,
                          rolling.length, 
                          no.simulations, 
                          set.seed.n,
                          show.limit, 
                          oos.realized){
  
  require(xts)
  require(dplyr)
  require(lubridate)
  
  #print("GARCH simulation")
  
  model_spec = ugarchspec(variance.model = list(model = model.df[,'type_model'],
                                                garchOrder = c(model.df[,'lag_arch'],
                                                               model.df[,'lag_garch'])),
                          mean.model = list(armaOrder = c(model.df[,'lag_ar'],
                                                          model.df[,'lag_ma'])),
                          distribution = model.df[,'type_dist'])
  
  fit = ugarchfit(data = training.data, spec = model_spec)
  
  n.sim <- no.simulations
  
  #print(paste0("simulating n=", n.sim, "... 1/3"))
  
  # code to fasten the function
  p.sim <- vector(mode='list', length=n.sim)
  
  garch.sim.extract <- function(fit, sim.length){
    p.sim <- ugarchsim(fit, n.sim=sim.length, startMethod="sample")
    p.sim <- p.sim@simulation$seriesSim
    return(p.sim)
  }
  
  # simulation
  garch.sim <- lapply(p.sim, FUN = function(x){garch.sim.extract(fit=fit, sim.length=nrow(oos.realized))})
  garch.sim <- as.data.frame(do.call(cbind,garch.sim))
  
  # determine quantiles to be estimated
  conf.levels <- confidence
  conf.names <- paste0("q_",conf.levels)
  conf.list <- list()
  
  if(method == "historical"){
    
    #print("historical simulation: identifying breaches... 3/3")
    
    # HISTORICAL
    mean.estimation <- list()
    for(i in 1:length(conf.levels)){
      conf.list[[i]] <- apply(as.data.frame(garch.sim), MARGIN = 1, FUN = function(x){quantile(x,conf.levels[i])})
      mean.estimation[[i]] <- mean(as.numeric(conf.list[[i]]))
    }
    
    conf.list <- as.data.frame(do.call(cbind, conf.list))
    names(conf.list) <- conf.names
    
    # estimate the mean
    mean.VaR <- apply(conf.list, MARGIN = 2, FUN = mean)
    mean.VaR <- as.data.frame(do.call("rbind", replicate(nrow(conf.list), mean.VaR, simplify = FALSE)))
    colnames(mean.VaR) <- paste0(colnames(mean.VaR),"_mean")
    
    # merge back
    conf.list <- cbind(conf.list,mean.VaR)
    garch.sim <- cbind(garch.sim,conf.list)
    garch.sim <- cbind(oos.realized, as.xts(garch.sim, order.by = ymd(index(oos.realized))))
    
    # historical sim
    sim.breaches <- list()
    expected.breaches <- list()
    real.breaches <- list()
    breach.prob <- list()
    for(i in 1:length(conf.levels)){
      #assign(paste0("VaR.breach_",conf.levels[[i]]), which(garch.sim[,1] < mean.estimation[[i]]))
      if(conf.levels[i] < 0.50){
        sim.breaches[[i]] <- which(garch.sim[,1] < mean.estimation[[i]])
      } else if(conf.levels[i] >= 0.50){
        sim.breaches[[i]] <- which(garch.sim[,1] > mean.estimation[[i]])
      }
      #assign(paste0("expected.breaches_",conf.levels[[i]]), round(conf.levels[[i]]*nrow(garch.sim),0))
      if(conf.levels[[i]] <= 0.50){
        expected.breaches[[i]] <- round(conf.levels[[i]]*nrow(garch.sim),0)
      } else if(conf.levels[i] > 0.50){
        expected.breaches[[i]] <- round((1-conf.levels[[i]])*nrow(garch.sim),0)
      }
      #expected.breaches[[i]] <- round(conf.levels[[i]]*nrow(garch.sim),0)
      real.breaches[[i]] <- length(sim.breaches[[i]])
      breach.prob[[i]] <- real.breaches[[i]] / nrow(garch.sim)
    }
    
    # plot
    cols.sampled <- sample(1:n.sim, show.limit*n.sim, replace = FALSE)
    
    colors <- c("red4","red","blue","blue4","orange","orange4","green","green4","purple","purple4")
    sim.plot <- plot(garch.sim[,cols.sampled], ylim=c(min(garch.sim[,1]) - 0.015,
                                                      max(garch.sim[,1]) + 0.015),  
                     type="p", cex=0.5, col=adjustcolor("grey", alpha=0.3),
                     grid.col=NA, main=model.description)
    sim.plot <- lines(garch.sim[,1], on=1, col="black", type="h")
    
    for(i in 1:length(conf.levels)){
      sim.plot <- lines(garch.sim[,colnames(mean.VaR)[i]], on=1, col=colors[i], lty=1, lwd=1)
    }
    
    # assign confidence labels
    legend.label <- vector(mode="character", length=length(conf.levels))
    legend.label[1] <- "GOBIXDR realizado"
    for(i in 2:(length(conf.levels)+1)){
      legend.label[i] <- paste0("VaR",conf.levels[i-1],"pct.banda [",round(tail(mean.VaR[1,i-1],1)*100,2),"%","]")
    }
    
    sim.plot <- addLegend("topleft", legend.names=legend.label, lty=rep(1,length(conf.levels)), lwd=rep(1,length(conf.levels)), col=c("black",colors), bty="n", y.intersp = 1.15)
    
    # assign breach label
    breaches.label <- vector(mode="character", length=length(conf.levels))
    breaches.label[1] <- "Retorno simulado"
    for(i in 2:(length(conf.levels)+1)){
      breaches.label[i] <- paste0("VaR",conf.levels[i-1],".brechas [",round(breach.prob[[i-1]][1]*100,2),"%"," , ", real.breaches[[i-1]][1],"/", expected.breaches[[i-1]][1],"]")
    }
    
    sim.plot <- addLegend("bottomleft", legend.names=breaches.label, pch=c(16,rep(16,length(conf.levels))), col=c("grey",colors), bty="n", y.intersp = 1.15)
    
    # shows only the first
    #sim.plot <- points(garch.sim[sim.breaches.99,1], col="white", bg="red", cex=1, pch=21)
    sim.plot <- points(garch.sim[sim.breaches[[1]],1], col="white", bg="red", cex=1, pch=21)
    sim.plot <- points(garch.sim[sim.breaches[[length(sim.breaches)]],1], col="red", bg="white", cex=1, pch=21)
    sim.plot
    
    results <- list(garch.sim, sim.plot)
    return(results)
    
  } else if(method == 'rolling'){
    
    convert.to.char <- function(v1) { deparse(substitute(v1)) }
    
    #print("rolling simulation: backtest... 3/3")
    
    # invert the order to match label colors
    conf.levels <- sort(conf.levels, decreasing = TRUE)
    
    garch.rolling = ugarchroll(model_spec, data = training.data, 
                               n.ahead=1, forecast.length = rolling.length, solver = "hybrid", 
                               refit.every=30, refit.window="moving", VaR.alpha=conf.levels)
    
    modeltest1 <- report(garch.rolling, type="VaR")
    modeltest2 <- report(garch.rolling, type="fpm")
    
    garch.object.name <- convert.to.char(garch.rolling)
    VaR <- paste0(garch.object.name,"@forecast$VaR")
    VaR <- eval(parse(text = VaR))
    realized.var.name <- paste0(garch.object.name,"@forecast$VaR$realized")
    realized.var <- eval(parse(text = realized.var.name))
    
    sigma <- garch.rolling@forecast$density
    
    # separate upper and lower
    split <-  as.numeric(gsub("\\D", "", colnames(VaR)))
    upper.cols <- which(split > 50)
    lower.cols <- which(split <= 50)
    
    upper.VaR <- VaR[,c(upper.cols,ncol(VaR))]
    lower.VaR <- VaR[,c(lower.cols,ncol(VaR))]
    
    upper.conf.levels <- conf.levels[conf.levels >= 0.5]
    lower.conf.levels <- conf.levels[conf.levels < 0.5]
    
    upper.breaches <- apply(as.data.frame(upper.VaR[,1:length(upper.conf.levels)]), MARGIN = 2, FUN = function(x){which(x < lower.VaR[,ncol(lower.VaR)])})
    upper.expected.breaches <- round((1-upper.conf.levels)*rolling.length,0)
    
    lower.breaches <- apply(as.data.frame(lower.VaR[,1:length(lower.conf.levels)]), MARGIN = 2, FUN = function(x){which(x > lower.VaR[,ncol(lower.VaR)])})
    lower.expected.breaches <- round(lower.conf.levels*rolling.length,0)
    
    sim.breaches <- c(upper.breaches,lower.breaches)
    expected.breaches <- c(upper.expected.breaches,lower.expected.breaches)
    
    real.breaches <- lapply(sim.breaches, FUN = length)
    breach.prob <- lapply(real.breaches, FUN=function(x){ x/rolling.length })
    
    cross <- na.omit(merge(training.data,oos.realized))
    
    garch.sim <- cbind(oos.realized, as.xts(garch.sim, order.by = ymd(index(cross))))
    VaR <- as.xts(VaR, order.by = ymd(rownames(VaR)))
    VaR <- VaR[,-ncol(VaR)]
    
    garch.sim <- na.omit(merge(garch.sim, VaR))
    
    # filter out just the realized vol and GARCH cutoff percentiles
    VaR.df <- garch.sim[,c(1,ncol(garch.sim):(ncol(garch.sim)-(length(conf.levels)-1)))]
    
    split.breaches <-  as.numeric(gsub("\\D", "", names(sim.breaches)))
    breaches.df <- data.frame(matrix(data = 0, ncol = length(split[1:length(conf.levels)]), nrow = nrow(VaR)))
    colnames(breaches.df) <- rep(paste0("VaR.",split.breaches,".breach"), 1)
    
    # impute the breaches
    for(i in 1:length(conf.levels)){
      breaches.df[sim.breaches[[i]],i] <- 1
    }
    # convert to XTS
    breaches.df <- as.xts(breaches.df, order.by = ymd(index(VaR)))
    # merge into VaR.df object
    VaR.df <- cbind(VaR.df, breaches.df)
    col.vars <- colnames(VaR.df)[ncol(VaR.df):(ncol(VaR.df)-length(conf.levels)+1)]
    
    #filter_all(as.data.frame(VaR.df), any_vars(.>0))
    a <- as.data.frame(VaR.df) %>% dplyr::select(col.vars)
    b <- filter_all(a, any_vars(.>0))
    
    breaches.list <- list()
    for(i in 1:ncol(b)){
      z <- b[b[,col.vars[i]] == 1,]
      #z %>% dplyr::select(col.vars[i])
      w <- z %>% dplyr::select(col.vars[i])
      w <- as.xts(w, order.by = ymd(rownames(w)))
      w <- na.omit(merge(VaR.df[,c(1,i+1)],w))
      w$days.between[2:nrow(w)] <- as.numeric(diff(ymd(index(w))))
      breaches.list[[i]] <- w
      #print(w)
    }
    
    #VaR$var01.rank <- ntile(VaR$`alpha(1%)`,100)
    # plot
    cols.sampled <- sample(1:n.sim, show.limit*n.sim, replace = FALSE)
    sim.breaches <- c(lower.breaches, upper.breaches)
    
    colors <- band.colors
    par(mfrow=c(1,1), mar=c(4,4,4,4))
    sim.plot <- plot(garch.sim[,cols.sampled], ylim=c(-0.06,0.06),  
                     type="p", cex=0.5, col=adjustcolor("grey", alpha=0.3),
                     grid.col=NA, main=model.description)
    sim.plot <- lines(garch.sim[,1], on=1, type="h", col="black")
    
    for(i in 1:length(conf.levels)){
      sim.plot <- lines(garch.sim[,ncol(garch.sim)-(i-1)], on=1, col=colors[i], lty=1, lwd=1)
    }
    
    # assign confidence labels
    legend.label <- vector(mode="character", length=length(conf.levels)+1)
    legend.label[1] <- "GOBIXDR realizado"
    for(i in 2:length(legend.label)){
      legend.label[i] <- paste0("VaR",conf.levels[i-1],"pct.banda")
    }
    legend.label <- legend.label[order(legend.label, decreasing = FALSE)]
    
    
    sim.plot <- addLegend("topleft", legend.names=legend.label, lty=rep(1,length(conf.levels)), lwd=rep(1,length(conf.levels)), col=c("black",colors), bty="n", y.intersp = 2)
    
    # assign breach label
    breaches.label <- vector(mode="character", length=length(conf.levels)+1)
    breaches.label[1] <- "Retorno simulado"
    for(i in 2:length(breaches.label)){
      breaches.label[i] <- paste0("VaR",conf.levels[i-1],".probabilidad_rebases [",round(breach.prob[[i-1]]*100,2),"%"," , ", real.breaches[[i-1]],"/", expected.breaches[[i-1]],"]")
    }
    breaches.label <- breaches.label[order(breaches.label, decreasing = FALSE)]
    
    
    sim.plot <- addLegend("bottomleft", legend.names=breaches.label, pch=c(16,rep(16,length(conf.levels))), col=c("grey",colors), bty="n", y.intersp = 2)
    sim.plot <- points(garch.sim[as.numeric(upper.breaches[[1]]),1], col="white", bg="red", cex=1.25, pch=21)
    sim.plot <- points(garch.sim[as.numeric(lower.breaches[[length(lower.breaches)]]),1], col="red", bg="white", cex=1.25, pch=21)
    sim.plot
    
    results <- list(garch.sim, VaR, breaches.list, sigma, modeltest1, modeltest2, sim.plot)
    return(results)
  }
}

garch.ts.forecast <- function(xts.object, 
                              label, 
                              garch.family, 
                              arma.params, 
                              garch.params, 
                              distribution,
                              start.date, 
                              end.date, 
                              days.forward, 
                              n.sim, 
                              set.seed.n,
                              show.limit,
                              plot){
  
  require(xts)
  require(dplyr)
  require(quantmod)
  require(rugarch)
  require(bizdays)
  
  symbol = label
  #sample = getSymbols(symbol, src = 'yahoo', from = start.date, to = end.date, auto.assign = FALSE)
  sample = xts.object
  sample$daily.vol <- diff(log(Cl(xts.object)))
  sample <- na.omit(sample)
  
  cut.date <- paste0(year(ymd(end.date))-1,"-12-31")
  
  ticker <- as.data.frame(sample) 
  ticker$date <- ymd(rownames(ticker))
  oos <- ticker %>% dplyr::filter(date > cut.date)
  ticker <- ticker %>% dplyr::filter(date < cut.date)
  
  oos <- as.xts(oos[,-ncol(oos)], order.by = oos[,'date'])
  ticker <- as.xts(ticker[,-ncol(ticker)], order.by = ticker[,'date'])
  
  # default GARCH(1,1)
  best_spec = ugarchspec(variance.model = list(model =  garch.family, 
                                               garchOrder = garch.params),
                         mean.model = list(armaOrder = arma.params),
                         distribution = distribution)
  
  my_best_garch <- ugarchfit(spec = best_spec, 
                             data = ticker$daily.vol)
  
  # basic GARCH criterias
  # infocriteria(my_best_garch)
  # coef(my_best_garch)
  # signbias(my_best_garch)
  # newsimpact(my_best_garch)
  # uncmean(my_best_garch)
  # uncvariance(my_best_garch)
  # persistence(my_best_garch)
  # halflife(my_best_garch)
  
  # SIMULACION
  asset1.ret <- ticker$daily.vol
  
  days.ahead = days.forward
  n.sim <- n.sim
  garch.sim <- matrix(nrow = days.ahead, ncol=n.sim)
  set.seed(set.seed.n)
  for(i in 1:n.sim){
    p.sim = ugarchsim(my_best_garch, n.sim=days.ahead, startMethod="sample")
    garch.sim[,i] <- p.sim@simulation$seriesSim
  }
  
  garch.sim <- as.data.frame(garch.sim)
  
  garch.sim$Q25 <- NA
  garch.sim$Q025 <- NA
  garch.sim$Q01 <- NA
  garch.sim$Q75 <- NA
  garch.sim$Q975 <- NA
  garch.sim$Q99 <- NA
  
  garch.sim$Q01 <- apply(garch.sim[,2:(ncol(garch.sim)-6)], FUN = function(x){quantile(na.omit(x),0.01)}, MARGIN = 1)
  garch.sim$Q025 <- apply(garch.sim[,2:(ncol(garch.sim)-6)], FUN = function(x){quantile(na.omit(x),0.025)}, MARGIN = 1)
  garch.sim$Q25 <- apply(garch.sim[,2:(ncol(garch.sim)-6)], FUN = function(x){quantile(na.omit(x),0.25)}, MARGIN = 1)
  garch.sim$Q75 <- apply(garch.sim[,2:(ncol(garch.sim)-6)], FUN = function(x){quantile(na.omit(x),0.75)}, MARGIN = 1)
  garch.sim$Q975 <- apply(garch.sim[,2:(ncol(garch.sim)-6)], FUN = function(x){quantile(na.omit(x),0.975)}, MARGIN = 1)
  garch.sim$Q99 <- apply(garch.sim[,2:(ncol(garch.sim)-6)], FUN = function(x){quantile(na.omit(x),0.99)}, MARGIN = 1)
  
  # adjust for business days
  business.calendar <- create.calendar('my_calendar', weekdays = c('saturday','sunday'))
  sim.dates <- bizdays::offset(ymd(last(index(ticker))+1), 1:days.ahead, cal = business.calendar)
  
  df <- cbind(as.data.frame(sim.dates),garch.sim)
  df <- as.xts(df[,-1], order.by=ymd(df$sim.dates))
  df <- as.numeric(tail(Cl(ticker),1)) * cumprod(1 + df)
  df$empty <- NA
  
  df.combined <- rbind(Cl(ticker),df$empty)
  df.combined <- merge(df.combined,df[,1:(ncol(df)-1)]) # do not merge the empty column
  
  last.day <- index(tail(Cl(ticker),1))
  last.day.level <- as.numeric(tail(Cl(ticker),1))
  last.price <- as.numeric(tail(Cl(oos),1))
  
  df.combined$Q01 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.01)}, MARGIN = 1)
  df.combined$Q025 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.025)}, MARGIN = 1)
  df.combined$Q25 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.25)}, MARGIN = 1)
  df.combined$Q75 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.75)}, MARGIN = 1)
  df.combined$Q975 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.975)}, MARGIN = 1)
  df.combined$Q99 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.99)}, MARGIN = 1)
  df.combined$mean <- as.numeric(apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){mean(na.omit(x))}, MARGIN = 1))
  
  quantInv <- function(distr, value) ecdf(distr)(value)
  level.prob <- quantInv(na.omit(as.numeric(tail(df.combined,1))),last.price) # has issues
  
  z <- as.data.frame(na.omit(merge(Cl(sample),df.combined[,2:ncol(df.combined)])))
  
  quantile.series <- apply(z[,2:ncol(df.combined)], MARGIN = 1, FUN = function(x) { quantInv(x,z[,1]) })
  quantile.series <- diag(quantile.series)
  
  
  if(plot == TRUE){
    
    cols_shown <- sample(2:(ncol(df.combined)-6), size=show.limit*(nrow(df.combined)), replace = FALSE)
    
    par(mfrow=c(1,1), mar=c(4,4,4,4))
    garch.sim.plot <- plot(tail(Cl(df.combined),days.forward*4), ylim=c(min(tail(Cl(oos),1), as.numeric(quantile(na.omit(df.combined[,2:ncol(df.combined)]),0.0025))), max(tail(Cl(oos),1), as.numeric(quantile(na.omit(df.combined[,2:ncol(df.combined)]),0.9975)))), 
                           main=paste0(symbol,": ARMA(", arma.params[1], ",",arma.params[2], ") + ", 
                                       garch.family, "(", garch.params[1], ",",garch.params[2],")"), grid.col=NA)
    garch.sim.plot <- lines(df.combined[,cols_shown], col=alpha("grey",0.5), on=1, lty=2, lwd=0.5)
    garch.sim.plot <- lines(df.combined[,'Q01'], col="red", on=1, lty=1, lwd=1.75)
    garch.sim.plot <- lines(df.combined[,'Q025'], col="red", on=1, lty=2, lwd=0.5)
    #garch.sim.plot <- lines(df.combined[,'Q25'], col="red", on=1, lty=3, lwd=1.5)
    garch.sim.plot <- lines(df.combined[,'mean'], col="grey40", on=1, lty=1, lwd=1.25)
    #garch.sim.plot <- lines(df.combined[,'Q75'], col="blue", on=1, lty=3, lwd=1.5)
    garch.sim.plot <- lines(df.combined[,'Q975'], col="blue", on=1, lty=2, lwd=0.5)
    garch.sim.plot <- lines(df.combined[,'Q99'], col="blue", on=1, lty=1, lwd=1.75)
    garch.sim.plot <- lines(Cl(oos), col="black", on=1, lty=1, lwd=1.5)
    garch.sim.plot <- points(tail(df.combined[,'mean'],1), col="red", pch=16)
    garch.sim.plot <- addLegend("bottomleft", c("Distribucion mostrada: min. 99.5%. Se muestran 1000 patrones simulados."), cex=0.8)
    garch.sim.plot <- addLegend("topleft", c(paste0("P01: ", round(tail(df.combined[,'Q01'],1),2), " (", round((tail(df.combined[,'Q01'],1)/last.day.level-1)*100,2), "%)"),
                                             paste0("P025: ", round(tail(df.combined[,'Q025'],1),2), " (", round((tail(df.combined[,'Q025'],1)/last.day.level-1)*100,2), "%)"),
                                             paste0("P975: ", round(tail(df.combined[,'Q975'],1),2), " (", round((tail(df.combined[,'Q975'],1)/last.day.level-1)*100,2), "%)"),
                                             paste0("P99: ", round(tail(df.combined[,'Q99'],1),2), " (", round((tail(df.combined[,'Q99'],1)/last.day.level-1)*100,2), "%)"),
                                             "realizado"), 
                                lty=c(1,2,2,1), lwd=c(1.75, 0.75, 0.75, 1.75), 
                                col=c("red","red","blue","blue","black"), cex=0.9, y.intersp = 1.75)
    garch.sim.plot
    
    # plot(ecdf(na.omit(as.numeric(df.combined[last.day]))),
    #      xlab="Nivel cumulativo esperado",
    #      ylab="Densidad de probabilidad cumulativa",
    #      main=paste(symbol),
    #      col="grey", lwd=2)
    # mtext("CDF of Simulated Price Level",  side=3)
    # points(x = last.price, y = level.prob, pch=16, col="red")
    # legend("topleft", legend=paste0("Al momento T: ", round(last.price,2) , " (" , round(level.prob*100,1),"%", ")"),
    #        col=c("red"), pch=c(16,16), bty="n")
    # box(col="grey")
    par(mfrow=c(1,1), mar=c(4,4,4,4))
    
    results <- list(df.combined, level.prob, quantile.series, garch.sim.plot)
    return(results)
    
  } else if(plot == FALSE){
    
    results <- list(df.combined, level.prob, quantile.series)
    return(results)
    
  }
}

risk.summary <- function(ticker, benchmark_ticker, start.date, end.date){
  data <- getSymbols(ticker, src = 'yahoo', from = start.date, to = end.date, auto.assign = FALSE)
  data$daily.rets <- diff(log(Cl(data)))
  
  benchmark <- getSymbols(benchmark_ticker, src = 'yahoo', from = start.date, to = end.date, auto.assign = FALSE)
  benchmark$benchmark.daily.rets <- diff(log(Cl(benchmark)))
  
  data <- merge(data,benchmark)
  data <- na.omit(data)
  
  return.distributions <- table.Distributions(data[,"daily.rets", drop=FALSE])
  table.DDs <- table.Drawdowns(data[,"daily.rets", drop=FALSE], top = 25, digits = 4, geometric = FALSE)
  downside.risk <- table.DownsideRisk(data[,"daily.rets", drop=FALSE], Rf=.02/12, MAR =.05/12, p=.95)
  calendar.rets <- calendar.ReturnTable(data[,"daily.rets", drop=FALSE])
  annualized.rets <- table.SFM(Ra = data[,"daily.rets", drop=FALSE], Rb= data[,"benchmark.daily.rets", drop=FALSE], Rf=0, digits=3)
  
  summary.data <- list(return.distributions, table.DDs, downside.risk, calendar.rets, annualized.rets)
  return(summary.data)
  
}

calendar.ReturnTable <- function(rets, digits = 3, percent = FALSE) {
  
  require(data.table)
  require(PerformanceAnalytics)
  
  pastePerc <- function(x) {return(paste0(x,"%"))}
  rowGsub <- function(x) {x <- gsub("NA%", "NA", x);x}
  
  dds <- apply.yearly(rets, maxDrawdown)
  rets <- apply.monthly(rets, Return.cumulative)

  dfRets <- cbind(year(index(rets)), month(index(rets)), coredata(rets))
  
  dfRets <- data.frame(dfRets)
  colnames(dfRets) <- c("Year", "Month", "Value")
  monthNames <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  for(i in 1:length(monthNames)) {
    dfRets$Month[dfRets$Month==i] <- monthNames[i]
  }
  dfRets <- data.table(dfRets)
  dfRets <- data.table::dcast(dfRets, Year~Month)
  
  dfRets <- data.frame(dfRets)
  yearNames <- dfRets$Year
  rownames(dfRets) <- yearNames; dfRets$Year <- NULL
  dfRets <- dfRets[,monthNames]
  
  yearlyRets <- apply.yearly(rets, Return.cumulative)
  dfRets$Annual <- yearlyRets
  dfRets$DD <- dds
  
  if(percent) {
    dfRets <- dfRets * 100
  }
  
  dfRets <- apply(dfRets, 2, round, digits)
  
  if(percent) {
    dfRets <- apply(dfRets, 2, pastePerc)
    dfRets <- apply(dfRets, 2, rowGsub)
    dfRets <- data.frame(dfRets)
    rownames(dfRets) <- yearNames
  }
  return(dfRets)
}


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
important.bond.tickers = c('EMB','BND','IEI','TLT','HYG','AGG')
getSymbols(important.bond.tickers, src = 'yahoo', from = '2014-01-01', auto.assign = TRUE)

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
selected.bond.tickers <- c("EMB","HYG","AGG","BND",
                           "EMLC","EBND","LEMB","EMHY","ELD","VWOB","EMTL")

bond.data <- list()
bond.performance <- list()
for(i in 1:length(selected.bond.tickers)){
  bond.data[[i]] <- getSymbols(selected.bond.tickers[i], src = 'yahoo', from = "2014-01-01", to = "2022-12-31", auto.assign = FALSE)
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
  bond.risk[[i]] <- risk.summary(ticker = selected.bond.tickers[i], benchmark_ticker = "IEF", start.date = "2014-01-01", end.date = "2022-12-31")
}

risk.resume <- do.call(cbind, lapply(bond.risk, FUN = function(x){x[[3]]}))
relative_performance.resume <- do.call(cbind, lapply(bond.risk, FUN = function(x){x[[5]]}))

colnames(risk.resume) <- selected.bond.tickers
colnames(relative_performance.resume) <- selected.bond.tickers

### Riesgos GOBIX y 7-10Y UST benchmark ----
benchmark <- getSymbols("IEF", src = 'yahoo', from = "2014-01-01", to = "2022-12-31", auto.assign = FALSE)
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


# Características ST GOBIXDR -----
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
# max_lag_AR <- 5 # m - parametro AR
# max_lag_MA <- 5 # n - parametro MA
# max_lag_ARCH <- 5 # p - parametro ARCH
# max_lag_GARCH <- 5 # q - parametro GARCH
# dist_to_use <- c('norm','snorm','ged','std','sstd','jsu') # ver rugarch::ugarchspecs
# models_to_estimate <- c('sGARCH', 'eGARCH', 'iGARCH', 'gjrGARCH', 'apARCH', 'csGARCH') # see rugarch::rugarchspec help for more
# 
# out <- find_best_arch_model(x = gobix.train[,'daily.vol'], 
#                             type_models = models_to_estimate,
#                             dist_to_use = dist_to_use,
#                             max_lag_AR = max_lag_AR,
#                             max_lag_MA = max_lag_MA,
#                             max_lag_ARCH = max_lag_ARCH,
#                             max_lag_GARCH = max_lag_GARCH)

## tabla resumen con resultados
# tab_out <- out$tab_out
#print(tab_out)

# models_names <- unique(tab_out$model_name)
# best_models <- c(tab_out$model_name[which.min(tab_out$AIC)],
#                  tab_out$model_name[which.min(tab_out$BIC)])

# Define el modelo que minimiza el BIC
# best_spec = ugarchspec(variance.model = list(model =  out$best_bic$type_model, 
#                                              garchOrder = c(out$best_bic$lag_arch,
#                                                             out$best_bic$lag_garch)),
#                        mean.model = list(armaOrder = c(out$best_bic$lag_ar, 
#                                                        out$best_bic$lag_ma)),
#                        distribution = out$best_bic$type_dist)

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
# for(i in 1:nrow(GOBIXDR.GARCH.best_model)){
#   
#   print(paste("PROGRESS:", i / nrow(GOBIXDR.GARCH.best_model)))
#   print(GOBIXDR.GARCH.best_model[i,])
#   
#   best_spec[[i]] = ugarchspec(variance.model = list(model =  GOBIXDR.GARCH.best_model$type_model[i], 
#                                                garchOrder = c(GOBIXDR.GARCH.best_model$lag_arch[i],
#                                                               GOBIXDR.GARCH.best_model$lag_garch[i])),
#                               mean.model = list(armaOrder = c(GOBIXDR.GARCH.best_model$lag_ar[i],
#                                                               GOBIXDR.GARCH.best_model$lag_ma[i])),
#                               distribution = GOBIXDR.GARCH.best_model$type_dist[i])
#   
#   my_best_garch[[i]] <- ugarchfit(spec = best_spec[[i]], data = gobix["2014/2018",'daily.vol'])
#   
#   symbol = "GOBIXDR"
#   n=262
#   
#   # CALIBRACION CORRIDA Y VALIDACION
#   # Ejecuta la predicción de 1-día hacía adelante (rolling GARCH) con recalibración cada 30 días
#   garch.roll[[i]] = ugarchroll(best_spec[[i]], gobix["2019/2022",'daily.vol'], n.ahead=1,
#                                forecast.length = n, solver = "hybrid", refit.every=30, 
#                                refit.window="moving", VaR.alpha=c(0.01,0.99))
#   
#   unconditioned.coverage.test[[i]] <- VaRTest(alpha=0.01, conf.level = 0.99, actual = garch.roll[[i]]@forecast$VaR$realized, VaR = garch.roll[[i]]@forecast$VaR$`alpha(1%)` )
#   conditioned.coverage.test[[i]] <- VaRDurTest(alpha=0.01, conf.level = 0.99, actual = garch.roll[[i]]@forecast$VaR$realized, VaR = garch.roll[[i]]@forecast$VaR$`alpha(1%)` )
#   
#   oos.df$model[i] <- as.character(GOBIXDR.GARCH.best_model$model_name[i])
#   oos.df$training.AIC[i] <- as.character(GOBIXDR.GARCH.best_model$AIC[i])
#   oos.df$training.BIC[i] <- as.character(GOBIXDR.GARCH.best_model$BIC[i])
#   oos.df$uncond.mean[i] <- uncmean(my_best_garch[[i]])
#   oos.df$uncond.variance[i] <- uncvariance(my_best_garch[[i]])
#   oos.df$persistence[i] <- persistence(my_best_garch[[i]])
#   oos.df$`half-life`[i] <- halflife(my_best_garch[[i]])
#   oos.df$excesos.esperados[i] <- as.numeric(unconditioned.coverage.test[[i]][1])
#   oos.df$excesos.registrados[i] <- as.numeric(unconditioned.coverage.test[[i]][2])
#   oos.df$H0[i] <- as.character(unconditioned.coverage.test[[i]][3])
#   oos.df$UC.LR.stat[i] <- as.numeric(unconditioned.coverage.test[[i]][4])
#   oos.df$UC.critico[i] <- as.numeric(unconditioned.coverage.test[[i]][5])
#   oos.df$UC.LRp[i] <- as.numeric(unconditioned.coverage.test[[i]][6])
#   oos.df$UC.decision[i] <- as.character(unconditioned.coverage.test[[i]][7])
#   oos.df$CC.H0[i] <- as.character(unconditioned.coverage.test[[i]][8])
#   oos.df$CC_LR.stat[i] <- as.numeric(unconditioned.coverage.test[[i]][9])
#   oos.df$CC.critico[i] <- as.numeric(unconditioned.coverage.test[[i]][10])
#   oos.df$CC.LRp[i] <- as.numeric(unconditioned.coverage.test[[i]][11])
#   oos.df$CC.decision[i] <- as.character(unconditioned.coverage.test[[i]][12])
#   oos.df$b[i] <- as.numeric(conditioned.coverage.test[[i]][1])
#   oos.df$uLL[i] <- as.numeric(conditioned.coverage.test[[i]][2])
#   oos.df$rLL[i] <- as.numeric(conditioned.coverage.test[[i]][3])
#   oos.df$LRp[i] <- as.numeric(conditioned.coverage.test[[i]][4])
#   oos.df$H0.CCduration[i] <- as.character(conditioned.coverage.test[[i]][5])
#   oos.df$decision[i] <- as.character(conditioned.coverage.test[[i]][6])
# }

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
                                             confidence = c(0.01,0.025,0.975,0.99),
                                             band.colors = c("red4","red","green","green4"),
                                             rolling.length = nrow(gobix["2019/2021",'daily.vol']),
                                             no.simulations = 10000,
                                             set.seed.n = 333,
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
                                               confidence=c(0.01,0.025,0.975,0.99),
                                               band.colors = c("red4","red","green","green4"),
                                               rolling.length = nrow(gobix["2022",'daily.vol']), 
                                               no.simulations = 10000,
                                               set.seed.n = 333,
                                               show.limit = 0.10, 
                                               oos.realized = gobix["2022",'daily.vol']) 
  #print(GARCH.daily.validation[[i]])
}


# GARCH model stats -----

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

# Forecast + functions -----

# n = nrow(final.models)

#  GARCH Daily Marginal OOS TEST
# GARCH.daily.marginal <- list()
# for(i in 1:n){
#   print(paste0(i,":",final.models$model_name[i]))
#   # ROLLING GARCH
#   GARCH.daily.marginal[[i]] <- garch.oos.sim(model.df = final.models[i,],
#                                              model.description = final.models$model_name[i],
#                                              training.data = gobix["2014/2021",'daily.vol'],
#                                              method = "rolling",
#                                              confidence = c(0.01,0.025,0.975,0.99),
#                                              band.colors = c("red4","red","green","green4"),
#                                              rolling.length = nrow(gobix["2019/2021",'daily.vol']),
#                                              no.simulations = 10000,
#                                              set.seed.n = 333,
#                                              show.limit = 0.10,
#                                              oos.realized = gobix["2019/2021",'daily.vol'])
#   print(GARCH.daily.marginal[[i]])
# }

#GARCH.daily.marginal[[1]][[3]]
#GARCH.daily.marginal[[1]][[3]][[1]]
#GARCH.daily.marginal[[1]][[3]][[4]]

#  GARCH Daily Marginal VALIDATION
# GARCH.daily.validation <- list()
# for(i in 1:n){
#   print(paste0(i,":",final.models$model_name[i]))
#   # ROLLING GARCH
#   GARCH.daily.validation[[i]] <- garch.oos.sim(model.df = final.models[i,],
#                                                model.description = final.models$model_name[i], 
#                                                training.data = gobix["2014/2022",'daily.vol'],
#                                                method = "rolling", 
#                                                confidence=c(0.01,0.025,0.975,0.99),
#                                                band.colors = c("red4","red","green","green4"),
#                                                rolling.length = nrow(gobix["2022",'daily.vol']), 
#                                                no.simulations = 10000,
#                                                set.seed.n = 333,
#                                                show.limit = 0.10, 
#                                                oos.realized = gobix["2022",'daily.vol']) 
#   #print(GARCH.daily.validation[[i]])
# }

# exesos
#GARCH.daily.validation[[1]][[3]]
#GARCH.daily.validation[[1]][[3]][[1]]
#GARCH.daily.validation[[1]][[3]][[4]]

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
                                                  set.seed.n = 1,
                                                  start.date="2019-01-01", end.date="2024-12-31",
                                                  days.forward=262, n.sim=10000, plot = TRUE)
  #print(GARCH.MC.bootstrapped[[i]])
}

# PEND
# GARCH.MC.bootstrapped.sim <- list()
# for(i in 1:n){
#   
#   name <- final.models$model_name[i]
#   modelo <- final.models$type_model[i]
#   ar <- final.models$lag_ar[i]
#   ma <- final.models$lag_ma[i]
#   arch <- final.models$lag_arch[i]
#   garch <- final.models$lag_garch[i]
#   distribution <- final.models$type_dist[i]
#   
#   GARCH.MC.bootstrapped.sim[[i]] <- list()
#   
#   for(j in 1:1000){
#     print(paste0("model:",j," sim:",i,":",name))
#     GARCH.MC.bootstrapped.sim[[i]][[j]] <- garch.ts.forecast(xts.object=gobix, label="GOBIXDR", 
#                                                     garch.family=modelo,
#                                                     arma.params=c(ar,ma),
#                                                     garch.params=c(arch,garch),
#                                                     distribution=distribution,
#                                                     set.seed.n = j,
#                                                     start.date="2019-01-01", end.date="2022-12-31",
#                                                     days.forward=262, n.sim=10000, plot = FALSE)
#     #print(GARCH.MC.bootstrapped[[i]])
#   }
# }

length(GARCH.MC.bootstrapped[[3]])
GARCH.MC.bootstrapped[[2]][[3]]

# aux function
quantInv <- function(distr, value) ecdf(distr)(value)

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

# VaR TABLES: Model TEST (2019-2021) -----

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


# VaR TABLES: VALIDATION (2019-2022) -----

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


# 01: Raising interest environment
#z <- as.data.frame(gobix[,1:2])
#z <- as.data.frame(regime.higher[,1:2])
#z <- as.data.frame(regime.flat[,1:2])
#z <- as.data.frame(forecast.VaR["2021-11-30/2022-06-30"])
#z <- as.data.frame(forecast.VaR["2021-11-30/2022"])
#z <- as.xts(z, order.by = ymd(rownames(z)))

#y <- as.data.frame(forecast.VaR)
#y <- as.xts(y, order.by = ymd(rownames(y)))

#VaR.higher <- na.omit(merge(z, y))
#names(VaR.higher)[3:4] <- c('alpha.1..','alpha.3..')
#baselIII.filtered.stress.VaR <- forecast.VaR["2022"]

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

# Plot 2. VaR + Sigma + Mu. Pendiente ajustar
# rolling.plot <- plot(forecast.VaR$realized, ylim=c(min(forecast.VaR$`alpha(1%)` -0.01), max(forecast.VaR$`alpha(99%)` +0.01)), main="", type="h", col="black", lwd=0.50, grid.col=NA)
# forecast.density.plot <- lines(forecast.VaR$`alpha(5%)`, on=1, col="orange", lty=1, lwd=0.75)
# forecast.density.plot <- lines(forecast.VaR$`alpha(1%)`, on=1, col="red", lty=1, lwd=1)
# forecast.density.plot <- lines(forecast.VaR$`alpha(95%)`, on=1, col="orange", lty=1, lwd=0.75)
# forecast.density.plot <- lines(forecast.VaR$`alpha(99%)`, on=1, col="red", lty=1, lwd=1)
# forecast.density.plot <- points(forecast.VaR[forecast.VaR$var99.breach==1,5], col="red", cex=1.25, pch=16)
# forecast.density.plot <- points(forecast.VaR[forecast.VaR$var01.breach==1,1], col="red", cex=1.25, pch=16)
# forecast.density.plot <- addLegend("topleft", legend.names=c("Banda percentil 99",
#                                                              "Banda percentil 95"), lty=c(1,1), col=c("red","orange"), bty="n", y.intersp = 0.75)
# 
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
# Gráficos no utilizados ----

## 001. RV: GOBIXDR vs. benchmarks importantes ----
par(mfrow=c(6,1), mar=c(0,4,4,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`TLT: US Treasuty vencimiento 20A`, 
     type="l", main="Variación diaria alternativas de Renta Fija", xaxt="n", xlab="", ylab="")
legend("topleft", legend="TLT: US Treasury vencimiento 20A", cex=1, bty = "n")
mtext("enero 2014 - diciembre 2022", cex=0.90, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`IEI: US Treasury vencimiento 5-7A`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="IEI: US Treasury vencimiento 5-7A", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`AGG: US Corporativos Grado Inversion`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="AGG: US Bonos Corporativos Grado Inversion", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`HYG: US Especulativos `, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="HYG: US Bonos Corporativos Especulativos", cex=1, bty = "n")
box(col="grey")

par(mar=c(0,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$`EMB: JPM EMBI`, 
     type="l", xaxt="n", xlab="", ylab="")
legend("topleft", legend="EMB: JPMorgan EMBI Index", cex=1, bty = "n")
box(col="grey")

par(mar=c(4,4,0,4))
plot(x = lubridate::ymd(rownames(bond.vol.df)), y = bond.vol.df$GOBIX.diario, 
     type="l", xlab="", ylab="")
legend("topleft", legend="GOBIXDR", cex=1, bty = "n")
box(col="grey")
par(mfrow=c(1,1), mar=c(4,4,4,4))
## 002. Performance relativo desde ATH 2021 (impacto cumulativo) ----
plot(df.indexed$GOBIXDR~df.indexed$date, type="l", ylab="", xlab="", main="GOBIXDR versus pares", ylim=c(0.65,1),
     cex.main=0.8, #change font size of title
     cex.sub=0.8, #change font size of subtitle
     cex.lab=0.8, #change font size of axis labels
     cex.axis=0.8)
lines(df.indexed[,'BND']~df.indexed$date, col="grey")
lines(df.indexed[,'AGG']~df.indexed$date, col="grey")
lines(df.indexed[,'HYG']~df.indexed$date, col="grey")
lines(df.indexed[,'EMB']~df.indexed$date, col="grey")
lines(df.indexed[,'LEMB']~df.indexed$date, col="grey")
lines(df.indexed[,'EMLC']~df.indexed$date, col="grey")
lines(df.indexed[,'EBND']~df.indexed$date, col="grey")
lines(df.indexed[,'EMHY']~df.indexed$date, col="grey")
lines(df.indexed[,'ELD']~df.indexed$date, col="grey")
lines(df.indexed[,'VWOB']~df.indexed$date, col="grey")
lines(df.indexed[,'EMTL']~df.indexed$date, col="grey")
lines(df.indexed[,'GOBIXDR']~df.indexed$date, col="red", lwd=2)
mtext("enero 2014 = 100", side=3, cex=0.75)
title(adj = 0.5, line = 0.85)
box(col="grey")
## 003. GARCH: distribución residuos -----
par(mfrow=c(1,2), mar=c(4,4,4,4))
plot(training.best_garch[[1]], which=8)
plot(training.best_garch[[1]], which=9)
par(mfrow=c(1,1))