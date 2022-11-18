# get_symbols <- function( symbols, first_date = "2002-01-01")
# {
#
#   default <- future::plan()
#   on.exit(future::plan(default), add = TRUE)
#   future::plan("multicore");
#
#   symbols <- yfR::yf_get( symbols,
#                           first_date = first_date,
#                           do_complete_data = TRUE,
#                           do_parallel = FALSE)
#   symbols <- split( symbols,
#                     f = symbols$ticker)
#
#   return(symbols)
# }
#
# test <- get_symbols( c("GOOG") )
#
# goog_test <- log( test[["GOOG"]][["price_close"]]/ dplyr::lag(test[["GOOG"]][["price_close"]]))
# goog_test <- na.omit(goog_test)
# goog_test <- (goog_test - mean(goog_test))/sd(goog_test)


pkgload::load_all()
vol_test <- rapcf( goog_test[2900:3750], burn_in = 400 )






