##' # Load in an example data set
##' library(survRM2)
##' D <- rmst2.sample.data()
##'
##' # Calculate the RMST based on the area under the KM approach
##' rmstKM(Surv(time, status) ~ arm, data=D, trunc=5, alpha=0.05)


# Function required for rmstKM
rmstfunc <- function(dat, tau, alpha){
  
  indat <- dat
  
  # Taken from rmst2 package and embedded here for version control
  onearm <- function(time, status, tau, alpha, arm){
    ft         <- survival::survfit(Surv(time, status) ~ 1)
    idx        <- ft$time <= tau
    wk.time    <- sort(c(ft$time[idx], tau))
    wk.surv    <- ft$surv[idx]
    wk.n.risk  <- ft$n.risk[idx]
    wk.n.event <- ft$n.event[idx]
    time.diff  <- diff(c(0, wk.time))
    areas      <- time.diff * c(1, wk.surv)
    rmst       <- sum(areas)
    wk.var     <- ifelse((wk.n.risk - wk.n.event) == 0, 0, wk.n.event/(wk.n.risk * (wk.n.risk - wk.n.event)))
    wk.var     <- c(wk.var, 0)
    rmst.var   <- sum(cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
    rmst.se    <- sqrt(rmst.var)
    out        <- matrix(0, 1, 5)
    out[1, ]   <- c(paste("RMST ",arm,":",sep=""),
                    rmst, rmst.se, rmst - qnorm(1 - alpha/2) * rmst.se,
                    rmst + qnorm(1 - alpha/2) * rmst.se)
    return(list("RMST"=out,"var"=rmst.var))
  }
  
  multiarm <- function(val){
    onearm(indat$time  [indat$arm == val],
           indat$status[indat$arm == val],
           tau   = tau,
           alpha = alpha,
           arm   = val)
  }
  
  uniqarm <- as.character(unique(sort(indat$arm)))
  allarms <- lapply(uniqarm,multiarm)
  
  if (length(uniqarm) > 1) {
    
    combi <- expand.grid(seq(1,length(uniqarm)),seq(1,length(uniqarm)))
    combi <- t(apply(combi, 1, sort))
    combi <- combi[combi[,1] != combi[,2],]
    combi <- matrix(combi[!duplicated(combi),],ncol=2)
    
    difference <- function(comb){
      dat1 <- allarms[[comb[1]]]
      dat2 <- allarms[[comb[2]]]
      rmst.diff.10     <- as.numeric(dat1$RMST[2]) - as.numeric(dat2$RMST[2])
      rmst.diff.10.se  <- sqrt(dat1$var + dat2$var)
      rmst.diff.10.low <- rmst.diff.10 - qnorm(1 - alpha/2) * rmst.diff.10.se
      rmst.diff.10.upp <- rmst.diff.10 + qnorm(1 - alpha/2) * rmst.diff.10.se
      rmst.diff.pval   <- pnorm(-abs(rmst.diff.10)/rmst.diff.10.se) * 2
      rmst.diff.result <- c(paste("RMST Dif (",uniqarm[comb[1]]," - ",uniqarm[comb[2]],"):",sep=""),
                            rmst.diff.10,rmst.diff.10.se, rmst.diff.10.low,
                            rmst.diff.10.upp, rmst.diff.pval)
      out              <- matrix(0, 1, 6)
      out[1, ]         <- rmst.diff.result
      return(out)
    }
    
    out22          <- t(apply(combi, 1, difference))
    out2           <- matrix(out22[,-1],ncol=5)
    class(out2)    <- "numeric"
    rownames(out2) <- out22[,1]
    colnames(out2) <- c("Est.",
                        "se",
                        paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""),
                        paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""),
                        "p")
    
  } else {out2 <- NULL}
  
  
  out11 <- lapply(paste("allarms[[",seq(1,length(uniqarm)),"]]$RMST",sep=""),function(x) eval(parse(text=x)) )
  out11 <- as.matrix(do.call(rbind, out11))
  out1  <- matrix(out11[,-1],ncol=4)
  class(out1)    <- "numeric"
  rownames(out1) <- out11[,1]
  colnames(out1) <- c("Est.",
                      "se",
                      paste("lower .", round((1 - alpha) * 100, digits = 0), sep = ""),
                      paste("upper .", round((1 - alpha) * 100, digits = 0), sep = ""))
  
  return(list("RMST"=out1,"diff"=out2))
}

############################################################
# Calculate the RMST based on the area under the KM approach
############################################################
rmstKM <- function(formula, data, trunc, alpha = 0.05) {
  Call <- match.call()
  temp <- Call[c(1, match(c("formula", "data"), names(Call), nomatch=0))]
  temp[[1]] <- as.name("model.frame")
  
  m        <- eval.parent(temp)
  response <- model.extract(m, "response")
  Terms    <- terms(formula, c("strata", "cluster"))
  ord      <- attr(Terms, "order")
  ll       <- attr(Terms, "term.labels")
  
  if (!survival::is.Surv(response)) {
    stop("Response must be a survival object")
  }
  
  if (attr(response, "type") != "right") {
    stop("Require right-censored data")
  }
  
  if (length(ord) & any(ord != 1))
    stop("Interaction terms are not valid for this function")
  
  if (length(ll) > 1) {
    stop("The survival formula must only contain arm as a covariate")
  }
  
  if (!is.numeric(trunc) | trunc < 0){
    stop("The trunction time must be a positive number")
  }
  
  if (identical(ll,character(0))){
    m[,2] <- 1
  }
  
  indat <- data.frame("time" = as.numeric(m[,1])[1:dim(m)[1]],
                      "status" = as.numeric(m[,1])[-1:-dim(m)[1]],
                      "arm"  = m[,2])
  
  
  if (any(tapply(indat$time, indat$arm, max) < trunc)) {
    stop(paste("The truncation time must be shorter than the minimum of the largest observed time in each group: ",
               sprintf("%.3f",min(tapply(indat$time, indat$arm, max))),
               sep=""))
  }
  
  out <- rmstfunc(dat=indat, tau=trunc, alpha=alpha)
  
  result <- list("RMST"=out$RMST, "diff"=out$diff, call=Call)
  class(result) <- c("rmstKM", "rmst")
  return(result)
}



#######
library(survRM2)
D <- rmst2.sample.data()
##' # Calculate the RMST based on the area under the KM approach
rmstKM(formula = Surv(time, status) ~ arm, data=D, trunc=5, alpha=0.05)

