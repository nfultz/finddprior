#' Find Dirichlet Prior
#'
#' A package for finding dirichlet priors based on box constraints.
#'
#' @name finddprior
#' @docType package
NULL

#' Find a dirichlet prior for box constraints
#'
#' @param lower lower bounds on the posterior proportion
#' @param upper upper bounds on the posterior proportion
#' @param size_box bounds of psuedo-sample-size
#' @param x0 initial values of params
#' @param alpha tuning paramter
#' @param opts optimization params for nloptr
#' @param local_opts local optimization params for nloptr
#' @param n scaling parameter, sample size of the likelihood
#'
#' @export
#' @importFrom nloptr nloptr
find_dirichlet_prior <- function(lower, upper, ..., size_box=c(0,10*n), x0, alpha = 1, opts=OPTS, local_opts=LOCAL_OPTS, n=1) {
  call <- sys.call()
  k <- length(lower)


  # reparameterizing so lambda3 = 1 - lambda1 - lambda2
  # x3 in result will be garbage, be sure to recalculate

  eval_f <- function( theta ) {
      x <- theta
      s <- theta[k]
      x[k] <- 1 - sum(theta[-k])

    list( "objective" = (s**2 + (n*alpha)**2 * crossprod(x)[1,1]) / 2,
          "gradient"  = c( (n * alpha)**2 * x[-k], s)
    )
  }

  # eval_f <- function( theta ) {
  #   x <- theta
  #   s <- theta[k]
  #   x[k] <- 1 - sum(theta[-k])
  #
  #   loss <- s**2 + alpha**2 * crossprod(x)[1,1]
  #
  #   list( "objective" = log(loss),
  #       "gradient"  = 1/loss * c(s, alpha* x[-k])
  #   )
  # }

  skeleton <- -diag(k-1)
  skeleton <- rbind(skeleton,1)
  skeleton <- rbind(skeleton, -skeleton)

  # constraint functions
  # inequalities
  # s*lambda_i / (n + s) > dl => (n+s)*dl - s lambda_i < 0
  eval_g_ineq <- function( theta ) {
      x <- theta
      s <- theta[k]
      x[k] <- 1 - sum(theta[-k])

      constr <- c(
        (n + s) * lower - s * x ,
       -(n + s) * upper + s * x + n,
       -x[k]
      )

      grad <- c(
        rbind(s * skeleton, 1),
        lower - x, x - upper, 0
      )

      return( list( "constraints"=constr, "jacobian"=grad ) )
  }

  lb <- lower
  lb[k] <- size_box[1]

  ub <- upper
  ub[k] <- size_box[2]

  if(missing(x0))
    x0 <- lb



  res <- nloptr( x0=x0,
                 eval_f=eval_f,
                 lb=lb,
                 ub=ub,
                 eval_g_ineq=eval_g_ineq,
                 eval_g_eq=NULL,
                 opts=opts)

  theta <- res$solution

  structure(
    list(nloptr=res, s=theta[k], lambda=c(theta[-k], 1 - sum(theta[-k])), call=call, n=n),
    class='find_dirichlet_prior'
  )
}

LOCAL_OPTS <- list( "algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel"  = sqrt(.Machine$double.eps) )

OPTS <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                "xtol_rel"  = sqrt(.Machine$double.eps),
                "maxeval"   = 20000,
                "local_opts" = LOCAL_OPTS )


#' @export
autoplot.find_dirichlet_prior <- function(result) {
  stopifnot(isNamespaceLoaded('ggplot2'))

  k <- length(result$lambda)
  z <- with(result, data.frame(
    i=gl(k,1),
    lambda,
    lower,
    upper,
    a_lower = s*lambda / (n + s),
    a_upper=(n + s*lambda) / (n + s))
    )
 ggplot(z) + aes(x=i) +
   geom_line(aes(y=lower), color='red', linetype='dashed', group='') +
   geom_line(aes(y=upper), color='red', linetype='dashed', group='') +
   geom_pointrange(aes(y=lambda, ymin=a_lower, ymax=a_upper)) +
   geom_abline(slope=0,intercept=1/k, color='cyan', linetype='dotted')+
   ylab(expression(hat(p))) +
   xlab(expression(lambda)) +
   with(result,
     ggtitle(sprintf(
       'Feasible ranges of posterior proportions\nPsuedo-sample-size: %d (%4.2f%%)',
       as.integer(s), round(s * 100 / (s+n), 2)
     ))) +
   theme_classic() + theme(plot.title = element_text(hjust = 0.5))

}

#' @export
search_alpha <- function(lower, upper, log_alpha=seq(25,0, length.out=25), ...,
                         last_s=Inf, max_tries=40) {

  out <- list()

  result2 <- find_dirichlet_prior(lower, upper)
  for(a in log_alpha) {
    # message(alpha)
    r <- find_dirichlet_prior(lower,upper, alpha=exp(a), ...)
    for(i in 1:max_tries) {

      if(r$nloptr$status %in% c(0,3) && r$s < last_s) break;
      # message("xxx")
      a = a + rnorm(1,0,.01)
      if(length(out) == 0) {
        x0 <- r$nloptr$solution
      } else {
        # Choose a previous one at random
        j <- sample(length(out), 1, prob = seq_along(out))
        x0 <- out[[j]]$res$nlopt$solution
      }
      r <- find_dirichlet_prior(lower,upper, alpha=exp(a), x0=x0, ...)
    }
    if(i == max_tries) next;
    # message(a, "     /", i)
    last_s <- r$s

    out[[length(out) + 1]] <- list(alpha=a, res=r, tries=i)
  }



  pss <-do.call(rbind.data.frame,
           lapply(out, function(x, r=x$res, l=r$lambda, s=r$s, n=r$n, a=x$alpha)
             cbind(
               data.frame(a=a, w=s / (n + s) ),
               as.data.frame(t(l))
             )))
  structure(pss, class=c('find_dirichlet_prior_search', 'data.frame'))
}


#' @export
autoplot.find_dirichlet_prior_search <- function(pss) {
  stopifnot(isNamespaceLoaded('ggplot2'))

  #rewritten to not need reshape2::melt
  pss <- reshape(pss, direction='long', idvar = 'a', v.names='value', varying=2:ncol(result2), times =colnames(result2)[-1], timevar='variable')


  ggplot(pss) +aes(x=a,y=value, color=variable) + geom_line() + geom_point(shape=4)
}
