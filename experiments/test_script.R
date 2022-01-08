
pkgload::load_all()

# es_run <- function( y, alpha = 0.3 ) {
#   for( i in seq_len(length(y)-1) ) {
#     y[i+1] <- alpha * y[i+1] + (1-alpha) * y[i]
#   }
#   return(y)
# }

# test <- log(tsibbledata::pelt$Lynx + 10000)


get_symbols <- function( symbols, first_date = "2002-01-01")
{
  
  default <- future::plan()
  on.exit(future::plan(default), add = TRUE)
  future::plan("multicore");
  
  symbols <- BatchGetSymbols::BatchGetSymbols( symbols,
                                               first.date = first_date,
                                               do.complete.data = TRUE,
                                               do.parallel = TRUE)
  symbols <- split( as.data.frame(symbols$df.tickers),
                    f = symbols$df.tickers$ticker)
  
  return(symbols)
}

test <- get_symbols( "GOOG" )
test <- (log( test[["GOOG"]][["price.close"]]) - log( dplyr::lag(test[["GOOG"]][["price.close"]]))) 
test <- na.omit(test)
test <- (test - mean(test))/sd(test)


