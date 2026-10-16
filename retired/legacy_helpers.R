# legacy_helpers.R -- CHECK_STEADY_STATE, retired from R/Polysome_functions_RAUC.R
# on 2026-10-09 (Claude Code session ef5348d2), moved verbatim (source lines
# 156-404 at commit 0399bd4). The plotting helpers and s2_approx that were
# moved with it were restored to R/Polysome_functions_RAUC.R afterwards.
# CHECK_STEADY_STATE is broken (tests/testthat/KNOWN_FAILURES.md #4, #5);
# the test suite's residual test replaces it. Kept for review; see README.md.

CHECK_STEADY_STATE <- function(unmarked, marked,parameters, imax=NA,post=F,verbose=T,tol=10){
  if(verbose==T){
    print("Steady State Check")
    print("Parameters:")
    print(parameters)
  }
  parameters <- unlist(parameters)
  if(is.na(imax)){
    imax <- length(unmarked)-1
  }
  tau <- rep(parameters[4],imax+1)*seq(0,imax,1)/parameters[4]
  kappa <- rep(parameters[1],imax+1)*seq(1,0,-1/imax)/parameters[4]
  mu <- rep(parameters[2]/parameters[4],imax+1)
  lambda <-  c(parameters[3]/parameters[4],rep(0,imax))
  delta = c(parameters[5]/parameters[4],rep(0,imax))
  marked_check <- numeric(imax+1)
  unmarked_check <- numeric(imax+1)
  for(i in 1:(imax+1)){
    if(i==1){
      unmarked_check[i] <- (-mu[i]-kappa[i])*unmarked[i] + tau[i+1]*unmarked[i+1] + lambda[i]
      marked_check[i] <- (mu[i])*unmarked[i] + tau[i+1]*marked[i+1] - delta[i]*marked[i]
    }else if(i<(imax+1)){
      unmarked_check[i] <- unmarked[i]*(-kappa[i]-tau[i]-mu[i]) + kappa[i-1]*unmarked[i-1]+tau[i+1]*unmarked[i+1]
      marked_check[i] <- (mu[i])*unmarked[i] + tau[i+1]*marked[i+1] - tau[i]*marked[i]
    }else {
      unmarked_check[i] <- unmarked[i]*(-tau[i]-mu[i]) + kappa[i-1]*unmarked[i-1]
      marked_check[i] <- (mu[i])*unmarked[i] - tau[i]*marked[i]
    }
    unmarked_check <- round(unmarked_check,tol)
    marked_check <- round(marked_check,tol)
    test <- numeric(2)
    if(verbose==T){
      print("Unmarked Class Check...")
      if(all(unmarked_check==0)){
        print("Passed")
        test[1] <- 0
      }else{
        print("Failed")
        test[1] <- 1
      }
      print("Marked Class Check...")
      if(all(marked_check==0)){
        print("Passed")
        test[2] <- 0
      }else{
        print("Failed")
        test[2] <- 1
      }
    } else {
      if(all(unmarked_check==0)){
        test[1] <- 0
      }else{
        test[1] <- 1
      }
      if(all(marked_check==0)){
        test[2] <- 0
      }else{
        test[2] <- 1
      }
    }
    if(post==T){
      out.tmp <- rbind(unmarked_check,marked_check)
      colnames(out.tmp) <- seq(0,imax,1)
      rownames(out.tmp) <- c("unmarked","marked")
      return(list(test,out.tmp))
    }else{
      return(test)
    }

  }
}
