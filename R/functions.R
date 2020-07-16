#' Edgeworth expansion for Welch's t-statistic
#'
#' This function provides approximation for the cumulative distribution function
#' of the sampling distribution of the Welch's t-statistic using Normal distribution, first order or second order Edgeworth expansion.
#'
#' @param x a real number.
#' @param order the order of edgeworth expansion. Valid options are 0, 1, and 2. If set to 0,
#' it reduces to approximation based on the central limit theorem and returns the CDF of standard normal distribution
#' evaluated at x.
#' @param n1 sample size for the sample from the first population.
#' @param n2 sample size for the sample from the second population.
#' @param mu1 mean of the first population.
#' @param mu2 mean of the second population.
#' @param sigma1 standard deviation of the first population.
#' @param sigma2 standard deviation of the second population.
#' @param gamma1 skewness of the first population.
#' @param gamma2 skewness of the second population.
#' @param tau1 kurtosis of the first population.
#' @param tau2 kurtosis of the second population.
#'
#' @return Edgeworth expansion evaluated at x.
#' @export
#'
#' @examples
#' t_edgeworth(1.96, order=2,
#' n1=20, n2=30,
#' mu1=0, mu2=0,
#' sigma1=1, sigma2=0.5,
#' gamma1=1, gamma2=0,
#' tau1=6, tau2=0)

t_edgeworth <- function(x, order=2,
                      n1, n2,
                      mu1, mu2,
                      sigma1, sigma2,
                      gamma1, gamma2,
                      tau1, tau2){
    # Low-level parameters
    n <- n1+n2
    lambda1 <- n1/n
    lambda2 <- n2/n
    d1 <- sqrt(lambda2*sigma1^2/(lambda2*sigma1^2 + lambda1*sigma2^2))
    d2 <- sqrt(lambda1*sigma2^2/(lambda2*sigma1^2 + lambda1*sigma2^2))
    d3 <- d1^2/sqrt(lambda1)
    d4 <- d2^2/sqrt(lambda2)
    delta <- mu1 - mu2
    # High-level parameters
    w <- delta/sqrt(sigma1^2/n1 +sigma2^2/n2)
    a1 <- -gamma1*d1*d3/2 + gamma2*d2*d4/2
    a2 <- w*a1
    a3 <- a1*2/3
    b1 <- 3 * w/8 *( d3^2*(tau1+2) + d4^2*(tau2+2) )
    b2 <- ( gamma1*d1*d3 - gamma2*d2*d4 )^2 + (d3^2+d4^2) + w^2/8*( d3^2*(tau1+2) + d4^2*(tau2+2) )
    b3 <- w * ( 7/8*( gamma1*d1*d3 - gamma2*d2*d4 )^2 + (d3^2+d4^2)/2 )
    b4 <- -(tau1*d3^2+tau2*d4^2)/12 + 2/3*(gamma1*d1*d3 - gamma2*d2*d4)^2 + (d3^2+d4^2)/2 + a2^2/2
    b5 <- a2*a3
    b6 <- a3^2/2
    # Polynomials
    H1 <- x-w
    H2 <- (x-w)^2 - 1
    H3 <- (x-w)^3 - 3*(x-w)
    H4 <- (x-w)^4 - 6*(x-w)^2 + 3;
    H5 <- (x-w)^5-10*(x-w)^3+15*(x-w)
    p1 <- -( a1 + a2*H1 + a3*H2 )
    p2 <- -( b1 + b2*H1 + b3*H2 + b4*H3 + b5*H4 + b6*H5 )
    # return
    if (order == 0){
        out <- stats::pnorm(x-w)
    } else if (order == 1){
        out <- stats::pnorm(x-w) + p1*stats::dnorm(x-w)/sqrt(n)
    } else if (order == 2){
        out <- stats::pnorm(x-w) + p1*stats::dnorm(x-w)/sqrt(n) + p2*stats::dnorm(x-w)/n
    } else {
        stop("Unavailable order. The order must be 0, 1, or 2.")
    }
    names(out) <- NULL
    return(out)

}

#' Cornish-Fisher expansion for Welch's t-statistic
#'
#' This function provides approximation for the quantile function of the sampling distribution
#' of the Welch's t-statistic using Cornish-Fisher expansion (up to second order).
#'
#' @param p a probability value.
#' @param order the order of Cornish-Fisher expansion. Valid options are 0, 1, and 2. If set to 0,
#' it reduces to a normal approximation and it returns the p-th percentile of standard normal distribution.
#'
#' @param n1 sample size for the sample from the first population.
#' @param n2 sample size for the sample from the second population.
#' @param mu1 mean of the first population.
#' @param mu2 mean of the second population.
#' @param sigma1 standard deviation of the first population.
#' @param sigma2 standard deviation of the second population.
#' @param gamma1 skewness of the first population.
#' @param gamma2 skewness of the second population.
#' @param tau1 kurtosis of the first population.
#' @param tau2 kurtosis of the second population.
#'
#' @return Cornish-Fisher expansion value evaluated at p.
#' @export
#'
#' @examples
#' t_cornish_fisher(0.9, order=2,
#' n1=60, n2=30,
#' mu1=0, mu2=0,
#' sigma1=1, sigma2=0.5,
#' gamma1=1, gamma2=0,
#' tau1=6, tau2=0)
#'
#' t_cornish_fisher(0.3, order=1,
#' n1=60, n2=30,
#' mu1=0, mu2=0,
#' sigma1=1, sigma2=0.5,
#' gamma1=1, gamma2=0,
#' tau1=6, tau2=0)
#'
t_cornish_fisher <- function(p, order=2,
                           n1, n2,
                           mu1, mu2,
                           sigma1, sigma2,
                           gamma1, gamma2,
                           tau1, tau2){
    # Low-level parameters
    n <- n1+n2
    lambda1 <- n1/n
    lambda2 <- n2/n
    d1 <- sqrt(lambda2*sigma1^2/(lambda2*sigma1^2 + lambda1*sigma2^2))
    d2 <- sqrt(lambda1*sigma2^2/(lambda2*sigma1^2 + lambda1*sigma2^2))
    d3 <- d1^2/sqrt(lambda1)
    d4 <- d2^2/sqrt(lambda2)
    delta <- mu1 - mu2
    # High-level parameters
    w <- delta/sqrt(sigma1^2/n1 +sigma2^2/n2)
    a1 <- -gamma1*d1*d3/2 + gamma2*d2*d4/2
    a2 <- w*a1
    a3 <- a1*2/3
    b1 <- 3 * w/8 *( d3^2*(tau1+2) + d4^2*(tau2+2) )
    b2 <- ( gamma1*d1*d3 - gamma2*d2*d4 )^2 + (d3^2+d4^2) + w^2/8*( d3^2*(tau1+2) + d4^2*(tau2+2) )
    b3 <- w * ( 7/8*( gamma1*d1*d3 - gamma2*d2*d4 )^2 + (d3^2+d4^2)/2 )
    b4 <- -(tau1*d3^2+tau2*d4^2)/12 + 2/3*(gamma1*d1*d3 - gamma2*d2*d4)^2 + (d3^2+d4^2)/2 + a2^2/2
    b5 <- a2*a3
    b6 <- a3^2/2
    # # Polynomials
    z <- stats::qnorm(p)
    H1z <- z
    H2z <- z^2 - 1
    H3z <- z^3 - 3*z
    H4z <- z^4 - 6*z^2 + 3
    H5z <- z^5 - 10*z^3 + 15*z
    p1z <- -( a1 + a2*H1z + a3*H2z )
    p1_d <- -( a2 + a3*z*2 )
    p2z <- -( b1 + b2*H1z + b3*H2z + b4*H3z + b5*H4z + b6*H5z )
    p11 <- - p1z
    p21 <- p1z*p1_d - 1/2*z*p1z^2 - p2z
    # return
    if (order == 0){
        out <- z
    } else if (order == 1) {
        out <- z + p11/sqrt(n)
    } else if (order == 2) {
        out <- z + p11/sqrt(n) + p21/n
    } else {
        stop("Unavailable order. The order must be 0, 1, or 2.")
    }
    names(out) <- NULL
    return(out)
}

#' The TCFU test
#'
#' This test is suitable for testing the equality of
#' two-sample means for the populations having unequal variances.
#' When the populations are not normally distributed, this test
#' can provide better type I error control and more accurate power than a large-sample t-test using normal
#' approximation. The critical values of the test are computed based on the
#' Cornish-Fisher expansion of the Welch's t-statistic. The order of the
#' Cornish-Fisher expansion is allowed to be 0, 1, or 2. More details
#' please refer to Zhang and Wang (2020).
#'
#' @param x1 the first sample.
#' @param x2 the second sample.
#' @param effectSize the effect size of the test. The default value is 0.
#' @param alternative the alternative hypothesis: "greater" for upper-tailed, "less" for lower-tailed, and "two.sided" for two-sided alternative.
#' @param alpha the significance level. The default value is 0.05.
#' @param order the order of the Cornish-Fisher expansion.
#'
#' @return test statistic, critical value, p-value, reject decision at the given significance level.
#' @export
#'
#' @references
#' Zhang, H. and Wang, H. (2020). Transformation tests and their asymptotic power in two-sample comparisons. Manuscript in review.
#'
#' @examples
#' x1 <- rnorm(20, 1, 3)
#' x2 <- rnorm(21, 2, 3)
#' tcfu(x1, x2, alternative = 'two.sided')
tcfu <- function(x1, x2, effectSize = 0,
                 alternative = 'greater',
                 alpha = 0.05, order = 2){
    n1 <- length(x1)
    n2 <- length(x2)
    n <- n1+n2
    mean1 <- drop(matrix(rep(1, n1), nrow=1) %*% matrix(x1)/n1)
    mean2 <- drop(matrix(rep(1, n2), nrow=1) %*% matrix(x2)/n2)
    var1 <- drop(matrix(rep(1, n1), nrow=1) %*% matrix((x1-mean1)^2)/(n1-1))
    var2 <- drop(matrix(rep(1, n2), nrow=1) %*% matrix((x2-mean2)^2)/(n2-1))
    Sp <- sqrt( var1/n1+var2/n2 )
    y1 <- (x1-mean1)/sqrt(var1)
    y2 <- (x2-mean2)/sqrt(var2)
    gamma1 <- (n1/((n1-1)*(n1-2))) * sum(y1^3)
    gamma2 <- (n2/((n2-1)*(n2-2))) * sum(y2^3)
    tau1 <- n1* rep(1, n1) %*% ((x1-mean1)^4)/ ((n1-1)*var1)^2 +6/(n1+1)
    tau2 <- n2* rep(1, n2) %*% ((x2-mean2)^4)/ ((n2-1)*var2)^2 +6/(n2+1)
    t_stat <- (mean1 - mean2 - effectSize)/Sp

    cf_est <- function(p, order){
        t_cornish_fisher(p, order=order,
                       n1=n1, n2=n2, mu1=0, mu2=0, sigma1=sqrt(var1), sigma2=sqrt(var2),
                       gamma1=gamma1, gamma2= gamma2, tau1=tau1, tau2=tau2)
    }
    ee_est <- function(x, order){
        t_edgeworth(x, order=order,
                  n1=n1, n2=n2, mu1=0, mu2=0, sigma1=sqrt(var1), sigma2=sqrt(var2),
                  gamma1=gamma1, gamma2= gamma2, tau1=tau1, tau2=tau2)
    }
    # Return the testing results
    if (alternative == 'greater')
    {
        cutoff <- cf_est(1-alpha, order)
        reject <- (t_stat >= cutoff)
        pvalue <- 1-ee_est(t_stat, order)
    }
    else if (alternative=='less')
    {
        cutoff <- cf_est(alpha, order)
        reject <-  (t_stat <= cutoff)
        pvalue <- ee_est(t_stat, order)
    }
    else if (alternative=='two.sided')
    {
        cutoff <- c(cf_est(alpha/2, order), cf_est(1-alpha/2, order))
        reject <- (t_stat >= cutoff[2] || t_stat <= cutoff[1])
        pvalue <- 1-ee_est(abs(t_stat), order) + ee_est(-abs(t_stat), order)
    }
    else
    {
        stop("Undefined alternative.")
    }
    names(reject) <- NULL
    names(t_stat) <- NULL
    return(list(stat = t_stat, cutoff = cutoff, pvalue=pvalue, reject = reject))
}

#' The transformation based test
#'
#' This test is suitable for testing the equality of
#' two-sample means for the populations having unequal variances.
#' When the populations are not normally distributed, the sampling distribution of
#' the Welch's t-statistic may be skewed. This test conducts transformations
#' of the Welch's t-statistic to make the sampling distribution more symmetric.
#' For more details, please refer to Zhang and Wang (2020).
#'
#' @param x1 the first sample.
#' @param x2 the second sample.
#' @param effectSize the effect size of the test. The default value is 0.
#' @param alternative the alternative hypothesis: "greater" for upper-tailed, "less" for lower-tailed, and "two.sided" for two-sided alternative.
#' @param alpha the significance level. The default value is 0.05.
#' @param type the type of transformation to be used. Possible choices are 1 to 4. They correspond to the TT1 to TT4 in Zhang and Wang (2020).
#' Which type provides the best test depends on the relative skewness parameter A in Theorem 2.2 of Zhang and Wang (2020).
#' In general, if A is greater than 3, \code{type} =3 is recommended. Otherwise, \code{type}=1 or 4 is recommended. The \code{type}=2
#' transformation may be more conservative in some cases and more liberal in some other cases than the \code{type}=1 and 4 transformations.   For more details, please refer to Zhang and Wang (2020).
#'
#' @return test statistic, critical value, p-value, reject decision at the given significance level.
#' @export
#'
#' @references
#' Zhang, H. and Wang, H. (2020). Transformation tests and their asymptotic power in two-sample comparisons Manuscript in review.
#'
#'
#' @examples
#' x1 <- rnorm(20, 1, 3)
#' x2 <- rnorm(21, 2, 3)
#' tt(x1, x2, alternative = 'two.sided', type = 1)
#'
#' #Negative lognormal versus normal data
#'  n1=50;  n2=33
#'  x1 = -rlnorm(n1, meanlog = 0, sdlog = sqrt(1)) -0.3*sqrt((exp(1)-1)*exp(1))
#'  x2 = rnorm(n2, -exp(1/2), 0.5)
#'  tt(x1, x2, alternative = 'less', type = 1)
#'  tt(x1, x2, alternative = 'less', type = 2)
#'  tt(x1, x2, alternative = 'less', type = 3)
#'  tt(x1, x2, alternative = 'less', type = 4)
#'
#' #Lognormal versus normal data
#'  n1=50;  n2=33
#'  x1 = rlnorm(n1, meanlog = 0, sdlog = sqrt(1)) + 0.3*sqrt((exp(1)-1)*exp(1))
#'  x2 = rnorm(n2, exp(1/2), 0.5)
#'  tt(x1, x2, alternative = 'greater', type = 1)
#'  tt(x1, x2, alternative = 'greater', type = 2)
#'  tt(x1, x2, alternative = 'greater', type = 3)
#'  tt(x1, x2, alternative = 'greater', type = 4)


tt <- function(x1, x2, alternative='greater', effectSize = 0, alpha=0.05, type = 1)
{
    n1 <- length(x1)
    n2 <- length(x2)
    n <- n1+n2
    mean1 <- drop(matrix(rep(1, n1), nrow=1) %*% matrix(x1)/n1)
    mean2 <- drop(matrix(rep(1, n2), nrow=1) %*% matrix(x2)/n2)
    var1 <- drop(matrix(rep(1, n1), nrow=1) %*% matrix((x1-mean1)^2)/(n1-1))
    var2 <- drop(matrix(rep(1, n2), nrow=1) %*% matrix((x2-mean2)^2)/(n2-1))
    s1 <- sqrt(var1)
    s2 <- sqrt(var2)
    Sp <- sqrt( var1/n1+var2/n2 )
    y1 <- (x1-mean1)/sqrt(var1)
    y2 <- (x2-mean2)/sqrt(var2)
    gamma1 <- (n1/((n1-1)*(n1-2))) * sum(y1^3)
    gamma2 <- (n2/((n2-1)*(n2-2))) * sum(y2^3)
    t_stat <- (mean1 - mean2 - effectSize)/Sp
    lambda1 <- n1/n
    lambda2 <- n2/n

    A_hat <- ( s1^3*gamma1/(lambda1^2)- s2^3*gamma2/(lambda2^2) ) / ( s1^2/lambda1 + s2^2/lambda2 )^(3/2)
    u <- t_stat/sqrt(n)
    # Transformations
    if (type == 1){
        transformedStat <- sqrt(n)*( u + A_hat*u^2/3 + A_hat^2*u^3/27 + A_hat/(6*n) )
    } else if (type == 2){
        transformedStat <- sqrt(n)*( (2/3/sqrt(n)*A_hat )^(-1)*(exp(2/3/sqrt(n)*A_hat*u)-1) + A_hat/6/n )
    } else if (type == 3){
        transformedStat <- sqrt(n)*( u + u^2 + u^3/3 +  A_hat/6/n )
    } else if (type == 4){
        transformedStat <- sqrt(n)*( (2/3*A_hat)^(-1)*(exp(2/3*A_hat*u)-1) + A_hat/6/n )
    } else {
        stop("Undefined transformation type. The value of type should be 1, 2, 3, or 4.")
    }
    names(transformedStat) <- NULL
    # Testing result
    if (alternative == 'greater')
    {
        cutoff <- stats::qnorm(1-alpha)
        reject <- (transformedStat >= cutoff)
        pvalue <- 1 - stats::pnorm(transformedStat)
    }
    else if (alternative == 'less')
    {
        cutoff <- stats::qnorm(alpha)
        reject <- (transformedStat <= cutoff)
        pvalue <- stats::pnorm(transformedStat)
    }
    else if (alternative == 'two.sided')
    {
        cutoff <- c(stats::qnorm(1-alpha/2), stats::qnorm(alpha/2))
        reject <- (transformedStat >= cutoff[1] || transformedStat <= cutoff[2])
        pvalue <- 2*(1-stats::pnorm(abs(transformedStat)))
    }
    else
    {
        stop("Undefined alternative.")
    }
    return(list(stat=transformedStat, cutoff=cutoff, pvalue = 1-stats::pnorm(transformedStat), reject=reject))
}

#' Compute t-statistic
#'
#' This is a helper function for the bootstrap test. It computes the t-statistic in a fast way.
#'
#' @param x1 the first sample.
#' @param x2 the second sample.
#'
#' @return the Welch's t-statistic.
#' @keywords internal
compute_t <- function(x1,x2){
    n1 <- length(x1)
    n2 <- length(x2)
    mean1 <- drop(matrix(rep(1, n1), nrow=1) %*% matrix(x1)/n1)
    mean2 <- drop(matrix(rep(1, n2), nrow=1) %*% matrix(x2)/n2)
    var1 <- drop(matrix(rep(1, n1), nrow=1) %*% matrix((x1-mean1)^2)/(n1-1))
    var2 <- drop(matrix(rep(1, n2), nrow=1) %*% matrix((x2-mean2)^2)/(n2-1))
    Sp <- sqrt( var1/n1+var2/n2 )
    (mean1 - mean2)/Sp
}

#' Bootstrap_t test for two-sample comparisons
#'
#' This function provides bootstrap approximation to the sampling distribution of the the
#' Welch's t-statistic
#'
#' @param x1 the first sample.
#' @param x2 the second sample.
#' @param B number of resampling rounds. Default value is 1000.
#' @param alternative the alternative hypothesis: "greater" for upper-tailed, "less" for lower-tailed, and "two.sided" for two-sided alternative.
#'
#' @return the p-value of the bootstrap_t test.
#' @export
#'
#' @examples
#' x1 <- rnorm(100, 0, 1)
#' x2 <- rnorm(100, 0.5, 2)
#' boot_test(x1, x2)
boot_test <- function(x1, x2, B = 1000, alternative = 'greater'){
    computed_stat <- compute_t(x1, x2)
    t_vect <- rep(NA, B)
    data_combine <- c(x1, x2)
    n <- length(data_combine)
    n1 <- length(x1)
    n2 <- length(x2)
    for(i in 1:B){
        boot.1 <- sample(1:n, n1, replace=T)
        boot.2 <- setdiff(1:n, boot.1)
        t_vect[i] <- compute_t(data_combine[boot.1], data_combine[boot.2])
    }
    if(alternative == 'greater'){
        pvalue = mean(t_vect > computed_stat)
    }
    else if (alternative == 'less'){
        pvalue = mean(t_vect < computed_stat)
    }
    else{
        pvalue = mean((t_vect < -abs(computed_stat) | t_vect > abs(computed_stat) ))
    }
    return(pvalue)
}


#' Adjusting power to assure actual size is within significance level
#'
#' It is common to use Monte Carlo experiments to evaluate the performance of
#' hypothesis tests and compare the empirical power among competing
#' tests. High power is desirable but difficulty arises when the actual sizes of
#' competing tests are not comparable. A possible way of tackling this issue is
#' to adjust the empirical power according to the actual size. This function
#' incorporates three types of power adjustment methods.
#'
#' @param size the empirical size of a test.
#' @param power the empirical power of a test.
#' @param method the power adjustment method. 'ZW' is the method proposed by Zhang and Wang (2020), 'CYS' is the method proposed by Cavus et al. (2019),
#' and 'probit' is the "method 1: probit analysis" in Lloyd (2005).
#'
#' @return the power value after adjustment.
#' @export
#'
#' @references
#' Lloyd, C. J. (2005). Estimating test power adjusted for size. Journal of Statistical Computation and Simulation, 75(11):921-933.
#'
#' Cavus, M., Yazici, B., & Sezer, A. (2019). Penalized power approach to compare the power of the tests when Type I error probabilities are different. Communications in Statistics-Simulation and Computation, 1-15.
#'
#' Zhang, H. and Wang, H. (2020). Transformation tests and their asymptotic power in two-sample comparisons Manuscript in review.
#'
#'
#' @examples
#' adjust_power(size = 0.06, power = 0.8, method = 'ZW')
#' adjust_power(size = 0.06, power = 0.8, method = 'CYS')
#' adjust_power(size = 0.06, power = 0.8, method = 'probit')
adjust_power <- function(size, power, method = 'ZW'){
    if (method == 'probit'){
        return(stats::pnorm(stats::qnorm(power)-stats::qnorm(size)+stats::qnorm(0.05)))
    } else if(method=='ZW' || method=='CYS'){
        if (size < 0.05){
            return(power)
        } else {
            if (method == 'ZW'){
                stats::pnorm(stats::qnorm(power)-stats::qnorm(size)+stats::qnorm(0.05))
            } else if (method == 'CYS'){
                power/sqrt(1+abs(1-size/0.05))
            }
        }
    }
}


#' Power-adjustment based on non-parametric estimation of the ROC curve
#'
#' It is common to use Monte Carlo experiments to evaluate the performance of
#' hypothesis tests and compare the empirical power among competing
#' tests. High power is desirable but difficulty arises when the actual sizes of
#' competing tests are not comparable. A possible way of tackling this issue is
#' to adjust the empirical power according to the actual size. This function
#' implements the "method 2: non-parametric estimation of the ROC curve" in
#' Lloyd (2005). For more details, please refer to the paper.
#'
#' @param stat_h0 simulated test statistics under the null hypothesis.
#' @param stat_ha simulated test statistics under the alternative hypothesis.
#' @param target_range_lower the lower end of the size range.
#' @param target_range_upper the upper end of the size range.
#'
#' @return the adjusted power.
#' @export
#'
#' @references
#' Lloyd, C. J. (2005). Estimating test power adjusted for size. Journal of Statistical Computation and Simulation, 75(11):921-933.
#'
#'
#' @examples
#' stath0 <- rnorm(100)
#' statha <- rnorm(100, mean=1)
#' pauc(stath0, statha, 0.01, 0.1)
pauc <- function(stat_h0, stat_ha,
                 target_range_lower,
                 target_range_upper){
    flag <- stat_h0 > stats::quantile(stat_h0, probs = target_range_lower) &
        stat_h0 < stats::quantile(stat_h0, probs = target_range_upper)
    stat_h0 <- stat_h0[flag]
    rank_sum <- length(stat_h0) * length(stat_ha)
    adj_power <- stats::wilcox.test(stat_ha, stat_h0)$statistic/rank_sum
    names(adj_power) <- NULL
    return(adj_power)
}


