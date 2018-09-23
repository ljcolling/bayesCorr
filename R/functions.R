posteriorSummary <- function( paramSampleVec){
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  return( list( mean =  meanParam , median =  medianParam , mode = modeParam ))
}


cor.mcmc<-setClass("cor.mcmc",
                   slots = list(rho = "numeric",
                                hdi.l = "numeric",
                                hdi.u = "numeric",
                                posterior.dist = "numeric",
                                x.name = "character",
                                y.name = "character"))

setMethod("show",
          "cor.mcmc",
          function(object){
            cat("\n\tBayesian parameter estimate of a correlation coefficient\n\n")
            cat("data: ",object@x.name,"and",object@y.name,"\n")
            cat("95 percent highest density inveral:\n")
            cat(unname(object@hdi.l),unname(object@hdi.u),"\n")
            cat("estimate:\n")
            cat("\trho\n")
            cat(object@rho,"\n")
            #cat("Rho:", sprintf("%.2f",object@rho),"\n")
            #cat("Lower HDI:", sprintf("%.2f",object@hdi.l),"\n")
            #cat("Upper HDI:", sprintf("%.2f",object@hdi.u),"\n")
          })

#' @export
setGeneric("posterior", function(x)
  standardGeneric("posterior"))

setMethod("posterior",
          "cor.mcmc",
          function(x){
            x@posterior.dist
          })



setMethod("$",
          "cor.mcmc",
          function(x,name){
            if(name == "estimate"){
              returnvalue = x@rho
              names(returnvalue) = "rho"
              return(returnvalue)}

            if(name == "conf.int"){
              returnvalue =unname(c(x@hdi.l,x@hdi.l))
              returnvalue = `attributes<-`(returnvalue,list("conf.level" = .95))
              return(returnvalue)
            }

            if(name == "data.name"){
              return(paste(x@x.name,"and",x@y.name))
            }
          })

cor.mcmc <- function(x, y) {

  x.name = deparse(substitute(x))
  y.name = deparse(substitute(y))
  # make sure you have complete cases
  z = data.frame(x = x, y = y)
  z = z[complete.cases(z), ]

  x = z$x
  y = z$y

  model_string <- "
  model {
  for(i in 1:n) {
  x[i,1:2] ~ dmnorm(mu[], prec[ , ])
  }

  # Make the covariance matrix and invert it because jags uses precision parameterisation.
  prec[1:2,1:2] <- inverse(cov[,])
  cov[1,1] <- sigma[1] * sigma[1]
  cov[1,2] <- sigma[1] * sigma[2] * rho
  cov[2,1] <- sigma[1] * sigma[2] * rho
  cov[2,2] <- sigma[2] * sigma[2]

  # Priors
  sigma[1] ~ dunif(sigmaLow1, sigmaHigh1)
  sigma[2] ~ dunif(sigmaLow2, sigmaHigh2)
  rho ~ dunif(-1, 1)
  mu[1] ~ dnorm(muM1, muP1)
  mu[2] ~ dnorm(muM2, muP2)

  # Random samples from the estimated distribution
  # Useful for some plotting functions
  x_rand ~ dmnorm(mu[], prec[ , ])
  }
  "
  thisData = data.frame(x = x, y = y)
  data_list = list(
    x = thisData[, c("x", "y")],
    n = nrow(thisData),
    sigmaLow1 = sd(x) / 1000,
    sigmaHigh1 = sd(x) * 1000 ,
    sigmaLow2 = sd(y) / 1000,
    sigmaHigh2 = sd(y) * 1000,
    muM1 = mean(x),
    muP1 = 0.000001 * 1 / sd(x) ^ 2 ,
    muM2 = mean(y),
    muP2 = 0.000001 * 1 / sd(y) ^ 2
  )

  mad.x = mad(thisData$x)
  mad.y = mad(thisData$y)
  if (mad(thisData$x) == 0) {
    mad.x = sd(thisData$x)
  }
  if (mad(thisData$y) == 0) {
    mad.y = sd(thisData$y)
  }



  inits_list = list(
    mu = c(median(thisData$x), median(thisData$y)),
    rho = cor(thisData$x,
              thisData$y, method = "spearman"),
    sigma = c(mad.x, mad.y),
    .RNG.name = "base::Mersenne-Twister",
    .RNG.seed = 6
  )

  jags_model <-
    rjags::jags.model(
      textConnection(model_string),
      data = data_list,
      inits = inits_list,
      n.adapt = 500,
      n.chains = 3,
      quiet = TRUE
    )
  update(jags_model, 1000)

  mcmc_samples <-
    rjags::coda.samples(jags_model,
                 c("mu", "rho", "sigma", "x_rand"),
                 n.iter = 10000,
                 thin = 1)

  samples_mat <- as.matrix(mcmc_samples)
  rho = samples_mat[, "rho"]


  new(Class = "cor.mcmc",
      rho = posteriorSummary(paramSampleVec = rho)$mode,
      hdi.l = HDInterval::hdi(rho)[1],
      hdi.u = HDInterval::hdi(rho)[2],
      posterior.dist = rho,
      x.name = x.name,
      y.name = y.name)


}


#' Estimate the correlation parameter of a bi-variate normal
#'
#' This function estimates the rho parameter for bivariate normal using
#' MCMC. It use a very wide prior on the mean (normal) and sd (uniform) for the marginals.
#'
#' It is used in the same way at the base cor.test function except with settting
#' method to "mcmc"
#'
#' @usage cor.test(x,y,method = "mcmc")
#' @param x,y numeric vectors of data values. `x` and `y` must have the same length
#'
#' @examples
#' ## Hollander & Wolfe (1973), p. 187f.
#' ## Assessment of tuna quality.  We compare the Hunter L measure of
#' ##  lightness to the averages of consumer panel scores (recoded as
#'##  integer values from 1 to 6 and averaged over 80 such values) in
#' ##  9 lots of canned tuna.

#' x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
#' y <- c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)

#' ##  The alternative hypothesis of interest is that the
#' ##  Hunter L value is positively associated with the panel score.

#' cor.test(x, y, method = "mcmc")
#' @export
cor.test<-function(x,y,method = "pearson",...){
  if(method == "mcmc"){
    cor.mcmc(x = x,y = y)
  } else {
    stats::cor.test(x,y,...)
  }
}
