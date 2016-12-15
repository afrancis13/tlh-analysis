library(quantmod)
library(forecast)
library(fGarch)
library(TSA)
library(tseries)
library(rugarch)
library(dplyr)
library(depmixS4)
library(ggplot2)
library(assertthat)

TRADING_DAYS_PER_YEAR = 256

#################### DATA STRUCTURES #################### 

# These are taken from reference documents publicly available
# on the Internet.

# Queue: http://www.exegetic.biz/blog/2013/11/implementing-a-queue-as-a-reference-class/

Queue <- setRefClass(
  Class = "Queue",
  fields = list(
    name = "character",
    data = "list"
  ),
  methods = list(
    size = function() {
      'Returns the number of items in the queue.'
      return(length(data))
    },
    #
    push = function(item) {
      'Inserts element at back of the queue.'
      data[[size()+1]] <<- item
    },
    #
    pop = function() {
      'Removes and returns head of queue (or raises error if queue is empty).'
      if (size() == 0) stop("queue is empty!")
      value <- data[[1]]
      data[[1]] <<- NULL
      value
    },
    #
    poll = function() {
      'Removes and returns head of queue (or NULL if queue is empty).'
      if (size() == 0) return(NULL)
      else pop()
    },
    #
    peek = function(pos = c(1)) {
      'Returns (but does not remove) specified positions in queue (or NULL if any one of them is not available).'
      if (size() < max(pos)) return(NULL)
      #
      if (length(pos) == 1) return(data[[pos]])
      else return(data[pos])
    },
    initialize=function(...) {
      callSuper(...)
      #
      # Initialise fields here (place holder)...
      #
      .self
    }
  )
)

# Priority Queue: https://www.r-bloggers.com/deriving-a-priority-queue-from-a-plain-vanilla-queue/

PriorityQueue <- setRefClass(
  "PriorityQueue",
  contains = "Queue",
  fields = list(
    priorities = "numeric"
  ),
  methods = list(
    push = function(item, priority) {
      'Inserts element into the queue, reordering according to priority.'
      callSuper(item)
      priorities <<- c(priorities, priority)
      order = order(priorities, decreasing = TRUE, partial = size():1)
      data <<- data[order]
      priorities <<- priorities[order]
    },
    pop = function() {
      'Removes and returns head of queue (or raises error if queue is empty).'
      if (size() == 0) stop("queue is empty!")
      priorities <<- priorities[-1]
      callSuper()
    })
)

#################### FUNCTIONS #################### 

# An algorithm to print the highest correlated holdings within an index fund.
# This would be used in something like Direct Indexing (DI), where individual/
# consistuent assets that are highely correlated are targeted for TLH.
getHighestCorrelations <- function(tickers) {
  closingPrices <- sapply(tickers, function(ticker) {
    tryCatch ({
      ticker_stock <- getSymbols(ticker, from = '01-01-2010', auto.assign = F)
      closingPriceColName <- sprintf("%s.Close", ticker)
      print(paste("Computing for ticker =", ticker))
      return(as.ts(ticker_stock[, closingPriceColName]))
    }, error = function(e) {
      return(NULL)
    })
  })
  
  # Use custom data structure.
  cc_pq <- PriorityQueue$new()
  
  # The dreaded for loop.
  for (i in 1:(length(tickers) - 1)) {
    stockOne <- closingPrices[[i]]
    for (j in (i+1):length(tickers)) {
      stockTwo <- closingPrices[[j]]
      if (!is.null(stockOne) && !is.null(stockTwo)) {
        trimmedStocks <- trimArrays(stockOne, stockTwo)
        stockOne <- trimmedStocks$vectorOne
        stockTwo <- trimmedStocks$vectorTwo
        corr <- cor(stockOne, stockTwo)
        name <- paste(tickers[i], tickers[j], sep = "-")
        cc_pq$push(name, corr)
      }
    }
  }
  
  for (k in 1:10) {
    print(paste(
      "The #", k, "most correlated stocks are (with correlation =",
      cc_pq$priorities[k], ")"
    ))
    print(cc_pq$pop())
  }
}

# Computes the best ARMA model from a grid of (p, q) possibilities
# by using Akaike's Information Criterion.
computeBestArma <- function(data, diff = 0, grid_length = 5) {
  minAic <- Inf
  best <- c(-1, diff, -1)
  for (i in 0:(grid_length - 1)) {
    for (j in 0:(grid_length - 1)) {
      tryCatch({
        arimaFit <- arima(data, order = c(i, diff, j))
        if (arimaFit$aic < minAic) {
          best <- c(i, diff, j)
          minAic <- arimaFit$aic
        }
      }, error = function(e){})
    }
  }
  print(best)
  return(best)
  sprintf("The minimum AIC for that parameter vector is %f", minAic)
}

# Computes the best GARCH model from a grid of (p, q) possibilities
# by using Akaike's Information Criterion.
computeBestGarch <- function(data, gridLength = 5) {
  minAic <- Inf
  best <- c(-1, -1)
  for (i in 1:(gridLength - 1)) {
    for (j in 1:(gridLength - 1)) {
      tryCatch({
        garchFit <- garch(data, order = c(i, j))
        if (AIC(garchFit) < minAic) {
          best <- c(i, j)
          minAic <- AIC(garchFit)
        }
      }, error = function(e){})
    }
  }
  print(best)
  return(best)
  sprintf("The minimum AIC for that parameter vector is %f", minAic)
}

# Performs a one-step ahead forecast for a time-series by computing the best ARIMA-GARCH Model
# using the previous windowLength data point, and then performing the prediction. Also provides
# a variance parameter. I return both in a list.
forecastArimaGarch <- function(data, oneDayAheadForecastLength = 500, windowLength = 1000) {
  forecasts <- rep(0, oneDayAheadForecastLength)
  forecastsSd <- rep(0, oneDayAheadForecastLength)
  
  for (i in 1:oneDayAheadForecastLength) {
    sprintf("Currently executing iteration %d", i)
    
    dayBeforeDayToPredict <- length(r.vti_ts_close) - (oneDayAheadForecastLength - i) - 1
    returnDataInWindow <- data[(dayBeforeDayToPredict - windowLength):dayBeforeDayToPredict]
    
    returnDataInWindow.bestArmaOrder <- computeBestArma(returnDataInWindow)
    returnDataInWindow.bestGarchOrder <- computeBestGarch(returnDataInWindow)
    
    spec = ugarchspec(
      variance.model = list(garchOrder = returnDataInWindow.bestGarchOrder),
      mean.model = list(
        armaOrder = c(returnDataInWindow.bestArmaOrder[1], returnDataInWindow.bestArmaOrder[2]),
        include.mean = T
      ),
      distribution.model = 'sged'
    )
    
    fit = tryCatch(
      ugarchfit(
        spec, returnDataInWindow, solver = 'hybrid'
      ), error=function(e) e, warning=function(w) w
    )
    
    if (is(fit, "warning")) {
      forecasts[i] <- NA
    } else {
      returnDataInWindow.forecast <- ugarchforecast(fit, n.ahead = 1)
      forecasts[i] <- returnDataInWindow.forecast@forecast$seriesFor[1]
      forecastsSd[i] <- returnDataInWindow.forecast@forecast$sigmaFor[1]  
    }
  }

  return(list(Forecast = forecasts, ForecastSd = forecastsSd))
}

# Recomputes the assets dataframe.
rebalanceAssets <- function(assetDf, valueToSell, priceToSell, priceToBuy, typeToSell, taxRate, precompute = F) {
  # Sort on loss potential.
  lossPotentialAssetDf <- assetDf %>% arrange(numShares * priceToSell - initialValue)

  totalTaxCredit = 0
  i = 1
  while (valueToSell > 0) {
    row <- lossPotentialAssetDf[i,]
    row.initialDeposit <- row$initialValue
    row.initialPrice <- row$initialPrice
    row.numShares <- row$numShares
    row.currentValue <- row.numShares * priceToSell
    row.type <- row$type

    # If current value of these assets exceed the number to share,
    # only sell number to share. This requires creating a new row
    # in the dataframe.
    if (row.type == typeToSell) {
      if (row.currentValue > valueToSell) {
        numSharesInValueToSell <- valueToSell / priceToSell
        originalPrice <- row.initialPrice

        loss <- originalPrice * numSharesInValueToSell - priceToSell * numSharesInValueToSell
        assert_that(loss >= 0)

        taxCredit <- taxRate * loss
        if (!precompute) {
          lossPotentialAssetDf[i, "initialValue"] <- row.currentValue - valueToSell
          lossPotentialAssetDf[i, "initialPrice"] <- priceToSell
          lossPotentialAssetDf[i, "currentValue"] <- row.currentValue - valueToSell
          lossPotentialAssetDf[i, "numShares"] <- (row.currentValue - valueToSell) / priceToSell
          
          lossPotentialAssetDf[(nrow(assetDf) + 1), "initialValue"] <- valueToSell + taxCredit
          lossPotentialAssetDf[(nrow(assetDf) + 1), "initialPrice"] <- priceToBuy
          lossPotentialAssetDf[(nrow(assetDf) + 1), "currentValue"] <- valueToSell + taxCredit
          lossPotentialAssetDf[(nrow(assetDf) + 1), "numShares"] <- (valueToSell + taxCredit) / priceToBuy

          if (row.type == 1) {
            lossPotentialAssetDf[(nrow(assetDf) + 1), "type"] <- 2          
          } else {
            lossPotentialAssetDf[(nrow(assetDf) + 1), "type"] <- 1
          }
        }

        valueToSell <- 0
      } else {
        loss <- row.initialDeposit - row.currentValue
        taxCredit <- taxRate * loss

        if (!precompute) {
          lossPotentialAssetDf[i, "initialValue"] <- row.currentValue + taxCredit
          lossPotentialAssetDf[i, "initialPrice"] <- priceToBuy
          lossPotentialAssetDf[i, "currentValue"] <- row.currentValue + taxCredit
          lossPotentialAssetDf[i, "numShares"] <- (row.currentValue + taxCredit) / priceToBuy          

          if (row.type == 1) {
            lossPotentialAssetDf[(nrow(assetDf) + 1), "type"] <- 2          
          } else {
            lossPotentialAssetDf[(nrow(assetDf) + 1), "type"] <- 1
          }
        }

        valueToSell <- valueToSell - row.currentValue
      }
      totalTaxCredit = totalTaxCredit + taxCredit
    }
    i = i + 1
  }
  
  if (precompute) {
    return(totalTaxCredit)
  } else {
    return(lossPotentialAssetDf)
  }
}

# Assumes the following: there are only two highly-correlated assets,
# that a wash sale occurs if the same individual sells and re-buys the same
# assets within 30 trading days.
# Recurring frequency is the number of days between each recurring deposit.
# I assume one quarter = 252 trading days / 4 = 63.
taxLossHarvestSim <- function(priceDataOne, priceDataTwo,
                              initialDeposit = 500000, recurringDeposit = 5000,
                              recurringFrequency = 30, taxRate = 0.40) {
  # Keeps track of recently sold items, so that we do not incur a "Wash Sale"
  washSaleDayAssetOne <- 0
  washSaleDayAssetTwo <- 0

  # Assumes you've already computed forecasts for the return data,
  # since this is heavily computationally intensive.
  returnForecasts <- read.csv(file = "returnForecasts.csv")
  returnForecastsSd <- read.csv(file = "returnForecastsSd.csv")
  
  minVolatility <- min(returnForecastsSd)
  maxVolatility <- max(returnForecastsSd)

  # Examine 3 years worth of data.
  start <- 1
  end <- 500

  # A list used to partition the assets into the identity of the security
  # (1 = VTI, 2 = SCHB), and the initial price upon first buying the security.
  # This is an important data structure to maintain throughout the simulation.
  initialPrice <- priceDataOne[start + (length(priceDataOne) - 500)]

  numShares <- initialDeposit / initialPrice
  assetsTLH <- data.frame(
    initialValue   = c(initialDeposit),
    initialPrice   = c(initialPrice),
    currentValue   = c(initialDeposit),
    numShares      = c(numShares),
    type           = c(1)
  )

  assetsNormal <- data.frame(
    initialValue   = c(initialDeposit),
    initialPrice   = c(initialPrice),
    currentValue   = c(initialDeposit),
    numShares      = c(numShares)
  )

  daysPerQuarter <- 63
  daysPassed <- 1
  
  totalValueWithHarvesting <- rep(0, 500)  
  totalValueWithoutHarvesting <- rep(0, 500)
  
  percentageAssetOne <- rep(0, 500)
  percentageAssetTwo <- rep(0, 500)

  # For each day in the last three years.
  for (i in start:end) {
    sprintf("Starting iteration %d", i)

    priceOne <- priceDataOne[i + (length(priceDataOne) - 500)]
    priceTwo <- priceDataTwo[i + (length(priceDataOne) - 500)]
    
    assetsTLH$currentValue <- (priceOne * assetsTLH$numShares * (assetsTLH$type %% 2) +
                               priceTwo * assetsTLH$numShares * ((assetsTLH$type %% 3) - 1))

    assetsNormal$currentValue <- priceOne * assetsNormal$numShares

    returnForecast <- returnForecasts$x[i]
    volatilityForecast <- returnForecastsSd$x[i]

    canBuyAssetOne <- washSaleDayAssetOne == 0
    canBuyAssetTwo <- washSaleDayAssetTwo == 0

    # Look for opportunities to harvest tax losses for asset one.
    # Only sell when the forecast is positive. If negative, wait
    # for it to go down. If it doesn't, the forecast should adjust,
    # and become positive, at which point we sell.
    if (returnForecast > 0) {
      assetsOne <- assetsTLH %>% filter(type == 1)
      assetsTwo <- assetsTLH %>% filter(type == 2)
      
      assetsToSellOne <- assetsOne %>% filter(initialValue > currentValue)
      assetsToSellTwo <- assetsTwo %>% filter(initialValue > currentValue)

      sumLossesAssetOne <- sum(assetsToSellOne$currentValue)
      sumLossesAssetTwo <- sum(assetsToSellTwo$currentValue)
      
      priceDifferenceOne <- assetsOne$initialPrice - priceOne
      priceDifferenceTwo <- assetsOne$initialPrice - priceTwo

      # Amount of assets sold proportional to standard deviations
      # above the average for volatility, by the formula,
      # Total Value To Sell = Total At Loss Value *
      # (volatility - min(volatility)) / (max(volatility) - min(volatility))
      
      # Only sell if you can gain at least $1000.

      # These are purely heuristics, and this should be tinkered with,
      # perhaps using a Machine Learning approach or local search approaches
      # (i.e. Simulated Annealing).
      
      # Also, sell the assets with larger losses first!

      threshold <- 200

      if (canBuyAssetOne && canBuyAssetTwo && (sumLossesAssetOne != 0 || sumLossesAssetTwo != 0)) {
        if (sumLossesAssetOne > sumLossesAssetTwo) {
          # Sell 1, Buy 2.
          valueToSell <- sumLossesAssetOne * ((volatilityForecast - minVolatility) / (maxVolatility - minVolatility))
          computedTaxGain <- rebalanceAssets(assetsTLH, valueToSell, priceOne, priceTwo, 1, taxRate, precompute = T)
          if (computedTaxGain > threshold) {
            assetsTLH <- rebalanceAssets(assetsTLH, valueToSell, priceOne, priceTwo, 1, taxRate)
            washSaleDayAssetTwo <- 30           
          }
        } else {
          # Sell 2, Buy 1.
          valueToSell <- sumLossesAssetTwo * ((volatilityForecast - minVolatility) / (maxVolatility - minVolatility))
          computedTaxGain <- rebalanceAssets(assetsTLH, valueToSell, priceTwo, priceOne, 2, taxRate, precompute = T)
          if (computedTaxGain > threshold) {          
            assetsTLH <- rebalanceAssets(assetsTLH, valueToSell, priceTwo, priceOne, 2, taxRate)
            washSaleDayAssetOne <- 30
          }
        }
      }
    }

    # Recurring deposit is made. You should always be able to do this (otherwise,
    # the algorithm is broken).

    if (daysPassed %% recurringFrequency == 0) {
      assetsNormal[(nrow(assetsNormal) + 1),] <- c(recurringDeposit, priceOne, recurringDeposit, recurringDeposit / priceOne)
      if (canBuyAssetOne) {
        assetsTLH[(nrow(assetsTLH) + 1),] <- c(recurringDeposit, priceOne, recurringDeposit, recurringDeposit / priceOne, 1)
        } else {
        assetsTLH[(nrow(assetsTLH) + 1),] <- c(recurringDeposit, priceOne, recurringDeposit, priceTwo, recurringDeposit / priceTwo, 2)        
      }
    }

    totalValueWithHarvesting[i] = sum(assetsTLH$currentValue)
    totalValueWithoutHarvesting[i] = sum(assetsNormal$currentValue)
    percentageAssetOne[i] = sum(assetsOne$currentValue) / sum(assetsTLH$currentValue)
    percentageAssetTwo[i] = sum(assetsTwo$currentValue) / sum(assetsTLH$currentValue)

    daysPassed = daysPassed + 1
    
    if (washSaleDayAssetTwo != 0) {
      washSaleDayAssetTwo = washSaleDayAssetTwo - 1
    }
    if (washSaleDayAssetOne != 0) {
      washSaleDayAssetOne = washSaleDayAssetOne - 1
    }
  }

  return(list(
    withTLH = totalValueWithHarvesting,
    withoutTLH = totalValueWithoutHarvesting,
    percentageOne = percentageAssetOne,
    percentageTwo = percentageAssetTwo
  ))
}

#################### DATA COLLECTION ####################

# Vanguard - an index fund with the underlying index as CRSP US Total Market.
vti_stock <- getSymbols("VTI", from = '2010-01-01', auto.assign = F)
vti_ts_close <- as.ts(vti_stock$VTI.Close)
r.vti_ts_close <- diff(log(vti_ts_close))

schb_stock <- getSymbols("SCHB", from = '2010-01-01', auto.assign = F)
schb_ts_close <- as.ts(schb_stock$SCHB.Close)
r.schb_ts_close <- diff(log(schb_ts_close))

#################### MODEL CREATION ####################

# Returns.

r.vti_ts_close.arima <- arima(r.vti_ts_close, order = c(3, 0, 3))
r.vti_ts_close.garch <- garch(r.vti_ts_close, order = c(1, 1))

# Intervention Analysis.

detectAO(r.vti_ts_close.arima)
detectIO(r.vti_ts_close.arima)

detectAO(r.vti_ts_close.garch)
detectIO(r.vti_ts_close.garch)

r.vti_ts_close.arimax <- arimax(r.vti_ts_close, order = c(3, 0, 3),
                                io = c(87, 95, 400, 402, 403, 404, 410, 468, 482, 1419))

# Prices.

vti_ts_close.arimaOrder <- computeBestArma(vti_ts_close, diff = 1)
vti_ts_close.garchOrder <- computeBestGarch(vti_ts_close)

vti_ts_close.arima <- arima(vti_ts_close, order = c(3, 1, 1))
vti_ts_close.garch <- garch(vti_ts_close, order = c(4, 1))

detectAO(vti_ts_close.arima)
detectIO(vti_ts_close.arima)

detectAO(vti_ts_close.garch)
detectIO(vti_ts_close.garch)

vti_ts_close.arimax <- arimax(vti_ts_close, order = c(3, 1, 1),
                              io = c(401, 403, 1419, 1420, 1422, 1631))

# These lines are extremely computationally-intensive with oneDayAheadForecastLength = 500.
# As a result, I have saved them to disk. Use the CSV.

r.vti_ts_close.forecastList <- forecastArimaGarch(r.vti_ts_close)

write.csv(r.vti_ts_close.forecastList$Forecast, file="returnForecasts.csv", row.names=FALSE)
write.csv(r.vti_ts_close.forecastList$ForecastSd, file="returnForecastsSd.csv", row.names=FALSE)

r.vti_ts_close.upperBound <- r.vti_ts_close.forecastList$Forecast + 1.96 * (r.vti_ts_close.forecastList$ForecastSd)
r.vti_ts_close.lowerBound <- r.vti_ts_close.forecastList$Forecast - 1.96 * (r.vti_ts_close.forecastList$ForecastSd)

vti_ts_close.forecastList <- forecastArimaGarch(vti_ts_close)

write.csv(vti_ts_close.forecastList$Forecast, file="PriceForecasts.csv", row.names=FALSE)
write.csv(vti_ts_close.forecastList$ForecastSd, file="PriceForecastsSd.csv", row.names=FALSE)

log.vti_ts_close <- log(vti_ts_close)
log.vti_ts_close.forecastList <- forecastArimaGarch(log.vti_ts_close)

write.csv(log.vti_ts_close.forecastList$Forecast, file="LogForecasts.csv", row.names=FALSE)
write.csv(log.vti_ts_close.forecastList$ForecastSd, file="LogForecastsSd.csv", row.names=FALSE)

log.vti_ts_close.upperBound <- log.vti_ts_close.forecastList$Forecast + 1.96 * (log.vti_ts_close.forecastList$ForecastSd)
log.vti_ts_close.lowerBound <- log.vti_ts_close.forecastList$Forecast - 1.96 * (log.vti_ts_close.forecastList$ForecastSd)

#################### PLOTS ####################

# Basics.
plot(vti_stock)
barChart(vti_stock)
lineChart(vti_stock, color.vol = F, TA = NULL)

# ACF, PACF, EACF.
acf(r.vti_ts_close)
pacf(r.vti_ts_close)
eacf(r.vti_ts_close)

# Periodogram.
spec(r.vti_ts_close, kernel = kernel("daniell", m = 20),
     taper = 0.05, ci.plot = T, main = "Smoothed Periodogram for VTI Returns",
     xlab = "Frequency", ylab = "Spectrum")

# Cross Correlation.

ccf(vti_ts_close, schb_ts_close, main = "CCF for VTI-SCHB")

# McLeod-Li Test.

McLeod.Li.test(y = r.vti_ts_close, main = "McLeod-Li Test on VTI Residuals")

# Residual Test.

r.vti_ts_close.arimagarch.resid <- r.vti_ts_close.arimax$residuals / r.vti_ts_close.garch$fitted.values[,"sigt"]
plot(r.vti_ts_close.arimagarch.resid, type = "p",
     main = "Residual Plot for ARIMA(3, 0, 3) - GARCH(1, 1)",
     ylab = "Standardized Residual Value")

qqnorm(r.vti_ts_close.arimagarch.resid,
       main = "Normal Q-Q Plot for ARIMA-GARCH Standardized Residuals")
qqline(r.vti_ts_close.arimagarch.resid)

vti_ts_close.arimagarch.resid <- vti_ts_close.arimax$residuals / vti_ts_close.garch$fitted.values[,"sigt"]

jarque.bera.test(as.vector(r.vti_ts_close.arimagarch.resid)[2:length(r.vti_ts_close.arimagarch.resid)])
shapiro.test(r.vti_ts_close.arimagarch.resid)

ggplot(data = data.frame(
    Time = 1:500,
    Prices = r.vti_ts_close[(length(r.vti_ts_close) - 499):length(r.vti_ts_close)],
    Forecasts = r.vti_ts_close.forecastList$Forecast,
    ForecastUb = r.vti_ts_close.upperBound,
    ForecastLb = r.vti_ts_close.lowerBound),
    aes(Time)
  ) + 
  geom_line(aes(y = Prices, colour = "Logged Price Data"), colour = "black") +
  geom_line(aes(y = Forecasts, colour = "Forecasts"), colour = "blue") +
  geom_line(aes(y = ForecastUb, colour = "95% Confidence Interval"), colour = "red", linetype = "dotted") + 
  geom_line(aes(y = ForecastLb), colour = "red", linetype = "dotted") +
  theme(plot.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(size = rel(2.0)),
        axis.title.y = element_text(size = rel(2.0)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        legend.position = c(1, 0)) +
  labs(title = "Daily Returns & One-Day-Ahead ARIMA-GARCH Predictions")  

ggplot(data = data.frame(
  Time = 1:500,
  Prices = log.vti_ts_close[(length(log.vti_ts_close) - 499):length(log.vti_ts_close)],
  Forecasts = log.vti_ts_close.forecastList$Forecast,
  ForecastUb = log.vti_ts_close.upperBound,
  ForecastLb = log.vti_ts_close.lowerBound),
  aes(Time)
) + 
  geom_line(aes(y = Prices, colour = "Logged Price Data"), colour = "black") +
  geom_line(aes(y = Forecasts, colour = "Forecasts"), colour = "blue") +
  geom_line(aes(y = ForecastUb, colour = "95% Confidence Interval"), colour = "red", linetype = "dotted") + 
  geom_line(aes(y = ForecastLb), colour = "red", linetype = "dotted") +
  theme(plot.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(size = rel(2.0)),
        axis.title.y = element_text(size = rel(2.0)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        legend.position = c(1, 0)) +
  labs(title = "Daily Prices & One-Day-Ahead ARIMA-GARCH Predictions")

ggplot(data = data.frame(Returns = r.vti_ts_close, Time = 1:500),
       aes(Time, Returns)) +
  geom_line(colour = "blue") +
  theme(plot.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(size = rel(2.0)),
        axis.title.y = element_text(size = rel(2.0)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5))) +
  labs(title = "Daily Returns, January 2010 through December 2016")

ggplot(data = data.frame(Prices = vti_ts_close, Time = 1:length(vti_ts_close)),
       aes(Time, Prices)) +
  geom_line(colour = "dark green") +
  theme(plot.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(size = rel(2.0)),
        axis.title.y = element_text(size = rel(2.0)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5))) +
  labs(title = "Daily Returns, January 2010 through December 2016")

#################### TLH EVALUATION ####################

taxAdjustedReturns <- taxLossHarvestSim(vti_ts_close, schb_ts_close)

ggplot(data = data.frame(
  Time = 1:500,
  withTLH = taxAdjustedReturns$withTLH,
  withoutTLH = taxAdjustedReturns$withoutTLH),
  aes(Time)
) + 
  geom_line(aes(y = withTLH, colour = "With TLH"), colour = "black") +
  geom_line(aes(y = withoutTLH, colour = "Without TLH"), colour = "green") +
  theme(plot.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(size = rel(2.0)),
        axis.title.y = element_text(size = rel(2.0)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        legend.position = c(1, 0)) +
  labs(y = "Tax-Adjusted Returns",
       title = "Simulated Tax-Adjusted Returns with TLH and Without TLH")

ggplot(data = data.frame(
  Time = 1:500,
  percentageVTI = taxAdjustedReturns$percentageOne,
  percentageSCHB = taxAdjustedReturns$percentageTwo),
  aes(Time)
) + 
  geom_line(aes(y = percentageVTI, colour = "VTI"), colour = "blue") +
  geom_line(aes(y = percentageSCHB, colour = "SCHB"), colour = "red") +
  theme(plot.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(size = rel(2.0)),
        axis.title.y = element_text(size = rel(2.0)),
        axis.text.x = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        legend.position = c(1, 0)) +
  labs(y = "Percentage of Portfolio Value",
       title = "Portfolio Composition Over Time With TLH")

#################### FURTHER RESEARCH  ####################

# Fifty top holdings.

schb_holdings <- read.csv("~/Downloads/SCHB_Fund_Holdings_ 11-04-2016.csv", 1)
schb_symbols <- sapply(schb_holdings$Symbol[1:100], as.character)

# This is rather computationally intensive.
# Output:
# [1] "The # 1 most correlated stocks are (with correlation = 0.992987990343391)"
# [2] "GOOGL-GOOG"
# [3] "The # 2 most correlated stocks are (with correlation = 0.980077639283554)"
# [4] "UNH-AVGO"
# [5] "The # 3 most correlated stocks are (with correlation = 0.977059518125652)"
# [6] "HON-TMO"
# [7] "The # 4 most correlated stocks are (with correlation = 0.976231044890752)"
# [8] "FB-LMT"
# [9] "The # 5 most correlated stocks are (with correlation = 0.973139579016713)"
# [10] "PFE-CMCSA"
# [11] "The # 6 most correlated stocks are (with correlation = 0.970285210113544)"
# [12] "CVS-LOW"
# [13] "The # 7 most correlated stocks are (with correlation = 0.969827537486666)"
# [14] "PEP-MMM"
# [15] "The # 8 most correlated stocks are (with correlation = 0.969480576346857)"
# [16] "BAC-MS"
# [17] "The # 9 most correlated stocks are (with correlation = 0.969292200183947)"
# [18] "ACN-COST"
# [19] "The # 10 most correlated stocks are (with correlation = 0.967620520087802)"
# [20] "HD-LOW"

getHighestCorrelations(schb_symbols)
