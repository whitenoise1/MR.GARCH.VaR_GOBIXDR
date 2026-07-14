# Funciones extraídas de VR_y_RM_Aplicacion_VaR_en_GOBIXDR.R (copiadas textualmente).
# Versión modular: ver refactor/main.R. El script original sigue siendo la referencia canónica.
#
# NOTA: el script original REDEFINE plot.quadrants y rolling.var.backtest a mitad
# del análisis (con firmas/comportamiento distintos). Aquí se incluye la PRIMERA
# versión de cada una; las redefiniciones posteriores se mantienen en refactor/main.R
# en su posición original para preservar exactamente el comportamiento.

# ---- Métricas de riesgo ----
cVaR <- function(x, conf){
  threshold <- quantile(na.omit(x),conf)
  excesos.VaR.historico <- which(x < threshold)
  cVaR.historico <- mean(x[excesos.VaR.historico])
  return(cVaR.historico)
}

# ---- Utilidades ----
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

# aux function
quantInv <- function(distr, value) ecdf(distr)(value)

# ---- Pruebas estadísticas ----
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

# ---- Gráficos ----
# (primera versión; redefinida más adelante en main.R con argumento label.position)
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

# ---- Datos ----
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

# ---- GARCH: backtesting y simulación ----
# (primera versión; redefinida más adelante en main.R -- versión "solo downside"
#  y versión "upside + downside" con 4 bandas)
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

# ---- Resúmenes de riesgo ----
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
