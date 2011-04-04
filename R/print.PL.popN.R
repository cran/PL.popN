print.PL.popN <-
function(x, ...)
  {
  cat("___________________________________________________________________________________\n")
  cat("\nPopulation size estimation from capture-recapture data in a closed population when \n")
  cat("using partial likelihood and a ")
  cat(x$model.type) 
  cat(". \n")
  cat("___________________________________________________________________________________\n")
  cat("\nNumber of individuals captured at least once D = ") 
  cat(x$D) 
  cat(" over tau = ")
  cat(x$tau)
  cat(" capture occasions.\n")
  cat("\n-----------------------------------------------------------------------------------\n")
  cat("Estimate of population size = ")
  cat(x$N.hat)
  cat(" with standard error = ")
  cat(x$N.hat.se)
  cat(".\n-----------------------------------------------------------------------------------\n")
}

