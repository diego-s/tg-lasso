#' @useDynLib tg.lasso fit_tg_lasso
fit.tg.lasso <- function(X.c, Y.c, W.c, lambda.1, lambda.2, 
  lambda.3, intervals, allow.missing.Y, min.difference) {
  M <- ncol(X.c)
  T <- ncol(Y.c)
  if (is.null(W.c)) W.c <- matrix(0, nrow = M, ncol = T)
  W.c <- .Call("fit_tg_lasso", X.c, Y.c, W.c, lambda.1, lambda.2, lambda.3, 
    rep(1, ncol(Y.c)), allow.missing.Y, min.difference)
  return(W.c)
}

#'
#' Temporal Group LASSO
#'
#' This function fits a regularized regression model with the following 
#' objective function:
#' 
#' L(W) = arg min _W 0.5 (||Y - XW||_F + 0.5 lambda_1 ||W||_F + 
#'   0.5 lambda_2 ||WR||_F + lambda_3 ||W||_2,1), 
#' 
#' where X is a feature matrix, Y is a matrix with one row for each patient and 
#' one column for each timepoint and R is a matrix where R_i,j is 1 if i = j, 
#' -1 if i = j + 1, and 0 otherwise. This is a regularized regression with two 
#' L2 penalties, the first one induces shrinkage of the coefficients. The second
#' one induces temporal smoothness by penalizing differences between adjacent 
#' coefficients corresponding to the same variable but different timepoints.
#'
#' @import lbfgs
#' @import methods
#' @import numDeriv
#'
#' @param X A feature matrix, with one row per patient and one column per 
#' feature.
#' @param Y A longitudinal outcome matrix, with one row per patient and one 
#' column per timepoint.
#' @param lambda.1 Tuning parameter for the coefficient shrinkage inducing 
#' penalty.
#' @param lambda.2 Tuning parameter for the temporal smoothness inducing 
#' penalty.
#' @param lambda.3.range Tuning parameter range for the group penalty for each 
#' feature across time points. It is encouraged to pass multiple values at a 
#' time in order to take advantage of warm starts.
#' @param allow.missing.X If TRUE, missing feature values are imputed using 
#' the average value of the missing feature. Default FALSE.
#' @param allow.missing.Y If TRUE, missing outcome values are allowed in the
#' Y longitudinal outcome matrix. In this case, missing Y values result in zero
#' errors, and therefore have no influence in the model. Default FALSE.
#' @param allow.zero.sd If TRUE, features with zero standard deviations are 
#' allowed in the X matrix. Default FALSE.
#' @param intevals An optional vector of time intervals between the time points 
#' in the columns of Y. By default, they are all set to 1.
#' 
#' @return A list having the following members:
#' \itemize{
#'   \item{"X"}{The original X feature matrix.}
#'   \item{"X.mean"}{The vector of column means of X.}
#'   \item{"X.sd"}{The vector of column standard deviations of X.}
#'   \item{"Y"}{The original Y longitudinal outcome matrix.}
#'   \item{"Y.mean"}{The vector of column means of Y.}
#'   \item{"Y.sd"}{The vector of column standard deviations of Y.}
#'   \item{"W"}{The matrix of fitted regression coefficients.}
#'   \item{"lambda.1"}{The value of the lambda.1 used to fit the model.}
#'   \item{"lambda.2"}{The value of the lambda.2 used to fit the model.}
#'   \item{"lambda.3"}{The value of the lambda.3 used to fit the model.}
#'   \item{"allow.missing.X"}{Whether missing X values were allowed when fitting
#'     the model.}
#'   \item{"allow.missing.Y"}{Whether missing Y values were allowed when fitting
#'     the model.}
#'   \item{"intervals"}{The vector of time intervals between the time points in
#'   Y.}
#' }
#'
#' @examples
#'
#' # generate data
#' n <- 1000
#' p <- 30
#' t <- 17
#' X <- matrix(rnorm(n * p, mean = 3, sd = 1.5), nrow = n, ncol = p)
#' W <- matrix(rnorm(p * t, mean = 0, sd = 1.2), nrow = p, ncol = t)
#' Y <- X %*% W
#' 
#' # test parameter estimation
#' X <- matrix(rnorm(n * p, mean = 3, sd = 1.5), nrow = n, ncol = p)
#' W <- matrix(rnorm((p + 1) * t, mean = 0, sd = 1.2), nrow = (p + 1), ncol = t)
#' Y <- cbind(rep(1, nrow(X)), X) %*% W
#' model <- tg.lasso(X, Y, 0.1, 0.1, 0.00001)
#' correlation <- cor(c(W), c(model$W.list[[1]][[1]][[1]]))
#' cat(sprintf("Parameter correlation: %s.\n", formatC(correlation, 
#'   digits = 3)))
#' stopifnot(correlation > 0.99)
#'
#' @export
#'
tg.lasso <- function(X, Y, lambda.1.range, lambda.2.range, lambda.3.range, 
  allow.missing.X = FALSE, allow.missing.Y = FALSE, allow.zero.sd = FALSE, 
  intervals = NULL, min.difference = 1e-7) {
  X <- as.matrix(X)
  X.mean <- apply(X, 2, function(x) mean(x, na.rm = TRUE))
  X.sd <- apply(X, 2, function(x) sd(x, na.rm = TRUE))
  if (allow.missing.X) {
    X.fill <- cbind(apply(X, 1, function(x) {
      x[is.na(x)] <- X.mean[is.na(x)]
      return(x)
    }))
    if (!ncol(X.fill) == ncol(X)) X.fill <- t(X.fill)
    X.c <- t((t(X.fill) - X.mean) / X.sd)
  } else {
    X.c <- t((t(X) - X.mean) / X.sd)
  }
  if (allow.zero.sd) X.c[is.nan(X.c)] <- 0
  Y <- as.matrix(Y)
  Y.mean <- apply(Y, 2, function(y) mean(y, na.rm = TRUE))
  Y.sd <- apply(Y, 2, function(y) sd(y, na.rm = TRUE))
  Y.c <- t((t(Y) - Y.mean) / Y.sd)
  W.c.list <- lapply(lambda.1.range, function(i) lapply(lambda.2.range, 
    function(j) list()))
  W.list <- lapply(lambda.1.range, function(i) lapply(lambda.2.range, 
    function(j) list()))
  X.ratio <- X.mean / X.sd
  if (allow.zero.sd) X.ratio[is.nan(X.ratio)] <- 0
  lambda.1.range <- sort(lambda.1.range)
  lambda.2.range <- sort(lambda.2.range)
  lambda.3.range <- sort(lambda.3.range)
  for (i in 1:length(lambda.1.range)) {
    for (j in 1:length(lambda.2.range)) {
      for (k in 1:length(lambda.3.range)) {
        if (k > 1) { W.c.last <- W.c.list[[i]][[j]][[k - 1]] }
        else if (j > 1) { W.c.last <- W.c.list[[i]][[j - 1]][[k]] }
        else if (i > 1) { W.c.last <- W.c.list[[i - 1]][[j]][[k]] }
        else { W.c.last <- NULL }
        W.c.list[[i]][[j]][[k]] <- fit.tg.lasso(X.c, Y.c, W.c.last, 
          lambda.1.range[i], lambda.2.range[j], lambda.3.range[k], intervals, 
          allow.missing.Y, min.difference)
        W.list[[i]][[j]][[k]] <- rbind(Y.mean - colSums(t(t(W.c.list[[i]][[j]][[
          k]]) * Y.sd) * X.ratio), t(t(W.c.list[[i]][[j]][[k]]) * Y.sd) / X.sd)
      }
    }
  }
  for (i in 1:length(lambda.1.range)) {
    for (j in 1:length(lambda.2.range)) {
      for (k in 1:length(lambda.3.range)) {
        if (allow.zero.sd) {
          W.c.list[[i]][[j]][[k]][!is.finite(W.c.list[[i]][[j]][[k]])] <- 0
          W.list[[i]][[j]][[k]][!is.finite(W.list[[i]][[j]][[k]])] <- 0
        }
        if (!is.null(dim(X)) && !is.null(colnames(X))) {
          rownames(W.list[[i]][[j]][[k]]) <- c("intercept", colnames(X))
        }
      }
    }
  }
  model <- list(X = X, X.mean = X.mean, X.sd = X.sd, Y = Y, Y.mean = Y.mean, 
    Y.sd = Y.sd, W.list = W.list, W.c.list = W.c.list, 
    lambda.1.range = lambda.1.range, lambda.2.range = lambda.2.range, 
    lambda.3.range = lambda.3.range, allow.missing.X = allow.missing.X, 
    allow.missing.Y = allow.missing.Y, allow.zero.sd = allow.zero.sd)
  class(model) <- append(class(model), "tg.lasso")
  return(model)
}

#'
#' Predict Smooth Temporal Model
#'
#' @export
#'
predict.tg.lasso <- function(model, X, lambda.1 = NULL, lambda.2 = NULL, 
  lambda.3 = NULL) {
  X <- as.matrix(X)
  if (model$allow.missing.X) {
    X.fill <- cbind(apply(X, 1, function(x) {
      x[is.na(x)] <- model$X.mean[is.na(x)]
      return(x)
    }))
    if (!ncol(X.fill) == ncol(X)) X.fill <- t(X.fill)
    X <- cbind(rep(1, nrow(X)), X.fill)
  } else {
    X <- cbind(rep(1, nrow(X)), X)
  }
  if (is.null(lambda.1)) {
    if (length(model$lambda.1.range) != 1) stop("Model fit with more than one \
      lambda 1 value, must specify lambda 1 parameter")
    lambda.1.index <- 1
  } else {
    lambda.1.index <- which(model$lambda.1.range == lambda.1)
  }
  if (is.null(lambda.2)) {
    if (length(model$lambda.2.range) != 1) stop("Model fit with more than one \
      lambda 2 value, must specify lambda 2 parameter")
    lambda.2.index <- 1
  } else {
    lambda.2.index <- which(model$lambda.2.range == lambda.2)
  }
  if (is.null(lambda.3)) {
    if (length(model$lambda.3.range) != 1) stop("Model fit with more than one \
      lambda 3 value, must specify lambda 3 parameter")
    lambda.3.index <- 1
  } else {
    lambda.3.index <- which(model$lambda.3.range == lambda.3)
  }
  if (min(length(lambda.1.index), length(lambda.2.index), 
    length(lambda.3.index)) == 0) 
    stop("Model wasn't fit with the provided lambda values.")
  Y.hat <- X %*% model$W.list[[lambda.1.index]][[lambda.2.index]][[
    lambda.3.index]]
  return(Y.hat)
}

#'
#' Cross Validate Smooth Temporal Model
#'
#' @export
#'
cross.validate.tg.lasso <- function(X, Y, lambda.1.range, lambda.2.range, 
  lambda.3.range,  k, allow.missing.X = FALSE, allow.missing.Y = FALSE, 
  allow.zero.sd = FALSE, intervals = NULL) {
  indices <- rep(1:k, length.out = nrow(X))
  results <- matrix(NA, nrow = length(lambda.1.range) * 
    length(lambda.2.range) * length(lambda.3.range), ncol = k)
  colnames(results) <- 1:k
  for (index in unique(indices)) {
    X.training <- cbind(X[indices != index,])
    Y.training <- cbind(Y[indices != index,])
    X.test <- cbind(X[indices == index,])
    Y.test <- cbind(Y[indices == index,])
    model <- tg.lasso(X.training, Y.training, lambda.1.range, lambda.2.range, 
      lambda.3.range, allow.missing.X = allow.missing.X, 
      allow.missing.Y = allow.missing.Y, allow.zero.sd = allow.zero.sd, 
      intervals = intervals)
    results[,index] <- c(sapply(lambda.1.range, function(lambda.1) 
        sapply(lambda.2.range, function(lambda.2) sapply(lambda.3.range, 
        function(lambda.3) {
        Y.hat <- predict(model, X.test, lambda.1, lambda.2, lambda.3)
        rmse <- sqrt(mean((Y.test - Y.hat) ^ 2, na.rm = TRUE))
        return(rmse)
    }))))
  }
  rmse.means <- rowMeans(results)
  rmse.sds <- apply(results, 1, sd)
  rmse.mean <- min(rmse.means)
  rmse.sd <- rmse.sds[which.min(rmse.means)]
  results <- cbind(rep(lambda.1.range, each = length(lambda.2.range) * 
    length(lambda.3.range)), rep(lambda.2.range, each = length(lambda.3.range), 
    times = length(lambda.1.range)), rep(lambda.3.range, 
    times = length(lambda.1.range) * length(lambda.2.range)), results)
  colnames(results)[1:3] <- c("lambda_1", "lambda_2", "lambda_3")
  results <- as.data.frame(results)
  lambda.1.min <- results$lambda_1[which.min(rmse.means)]
  lambda.2.min <- results$lambda_2[which.min(rmse.means)]
  lambda.3.min <- results$lambda_3[which.min(rmse.means)]
  cv <- list(X = X, Y = Y, results = results, rmse.means = rmse.means,
    rmse.sds = rmse.sds, rmse.mean = rmse.mean, rmse.sd = rmse.sd, 
    lambda.1.min = lambda.1.min, lambda.2.min = lambda.2.min, 
    lambda.3.min = lambda.3.min)
  return(cv)
}

#'
#' Tune Smooth Temporal Model
#'
#' @export
#'
tune.tg.lasso <- function(X, Y, lambda.1.range, lambda.2.range, 
  lambda.3.range, k, allow.missing.X = FALSE, allow.missing.Y = FALSE, 
  allow.zero.sd = FALSE, intervals = NULL) {
  size <- length(lambda.1.range) * length(lambda.2.range) * 
    length(lambda.3.range)
  results <- data.frame(lambda_1 = rep(NA, size), lambda_2 = rep(NA, size), 
    lambda_3 = rep(NA, size), mean = rep(NA, size), sd = rep(NA, size))
  for (i in 1:length(lambda.1.range)) {
    lambda.1 <- lambda.1.range[i]
    for(j in 1:length(lambda.2.range)) {
      lambda.2 <- lambda.2.range[j]
      cv <- cross.validate.tg.lasso(X, Y, lambda.1, lambda.2, 
        lambda.3.range, k, allow.missing.X = allow.missing.X, 
        allow.missing.Y = allow.missing.Y, allow.zero.sd = allow.zero.sd, 
        intervals = intervals)
      cat(sprintf("%s", apply(cbind(lambda.1, lambda.2, lambda.3.range, 
        cv$rmse.means, cv$rmse.sds), 1, function(row) sprintf(
        "lambda 1: %f - lambda 2: %f - lambda 3: %f - CV RMSE: %f (%f).\n", 
        row[1], row[2], row[3], row[4], row[5]))))
      l <- ((i - 1) * length(lambda.1.range) * length(lambda.2.range)) + 
        ((j - 1) * length(lambda.2.range)) + 1
      results[l:(l + length(lambda.3.range) - 1),] <- cbind(lambda.1, lambda.2, 
        lambda.3.range, cv$rmse.means, cv$rmse.sds)
    }
  }
  results <- results[order(results$lambda_1 + results$lambda_2 + 
    results$lambda_3),]
  results$upper <- results$mean + (results$sd / sqrt(k))
  index.min <- which.min(results$mean)
  index.1se <- rev(which(results$mean <= results$upper[index.min]))[1]
  lambda.1.min <- results$lambda_1[index.min]
  lambda.2.min <- results$lambda_2[index.min]
  lambda.3.min <- results$lambda_3[index.min]
  mean.min <- results$mean[index.min]
  sd.min <- results$sd[index.min]
  lambda.1.1se <- results$lambda_1[index.1se]
  lambda.2.1se <- results$lambda_2[index.1se]
  lambda.3.1se <- results$lambda_3[index.1se]
  mean.1se <- results$mean[index.1se]
  sd.1se <- results$sd[index.1se]
  tuning <- list(results = results, lambda.1.min = lambda.1.min, 
    lambda.2.min = lambda.2.min, lambda.3.min = lambda.3.min, 
    mean.min = mean.min, sd.min = sd.min, lambda.1.1se = lambda.1.1se, 
    lambda.2.1se = lambda.2.1se, lambda.3.1se = lambda.3.1se, 
    mean.1se = mean.1se, sd.1se = sd.1se)
  return(tuning)
}

#'
#' Backward Ranking
#'
#' Eliminates and replaces features one at a time, ranking features based on 
#' the performance obtained by eliminating them from the feature set. The 
#' minimum and maximum values of the mean CV RMSE, the features corresponding to
#' these values, as well as the mean and standard deviation of the CV RMSE.
#'
#' @param X A feature matrix, with one row per patient and one column per 
#' feature.
#' @param Y A longitudinal outcome matrix, with one row per patient and one 
#' column per timepoint.
#' @param lambda.1.range A vector of lambda.1 values to try.
#' @param lambda.2.range A vector of lambda.2 values to try.
#' @param lambda.3.range A vector of lambda.3 values to try.
#' @param allow.missing.X If TRUE, missing feature values are imputed using 
#' the average value of the missing feature. Default FALSE.
#' @param allow.missing.Y If TRUE, missing outcome values are allowed in the
#' Y longitudinal outcome matrix. In this case, missing Y values result in zero
#' errors, and therefore have no influence in the model. Default FALSE.
#' @param allow.zero.sd If TRUE, features with zero standard deviations are 
#' allowed in the X matrix. Default FALSE.
#' @param intevals An optional vector of time intervals between the time points 
#' in the columns of Y. By default, they are all set to 1.
#' 
#' @return A list having the following members:
#' \itemize{
#'   \item{"X"}{The original X feature matrix.}
#'   \item{"Y"}{The original Y longitudinal outcome matrix.}
#'   \item{"ranked.features"}{The set of features to rank (may be all 
#'     features, for example).}
#'   \item{"min.index"}{The column index of the feature having the smallest 
#'     mean CV RMSE.}
#'   \item{"min.feature"}{The feature having the smallest mean CV RMSE.}
#'   \item{"min.mean"}{The mean CV RMSE associated to the feature having the
#'     smallest mean CV RMSE.}
#'   \item{"min.sd"}{The sd of the CV RMSE associated to the feature having the
#'     smallest mean CV RMSE.}
#'   \item{"max.index"}{The column index of the feature having the largest 
#'     mean CV RMSE.}
#'   \item{"max.feature"}{The feature having the largest mean CV RMSE.}
#'   \item{"max.mean"}{The mean CV RMSE associated to the feature having the
#'     largest mean CV RMSE.}
#'   \item{"max.sd"}{The sd of the CV RMSE associated to the feature having the
#'     largest mean CV RMSE.}
#'   \item{"results"}{Table of results of the backward ranking."}
#' }
#'
#' @export
#'
backward.ranking <- function(X, Y, ranked.features, lambda.1.range, 
  lambda.2.range, lambda.3.range, k, allow.missing.X = FALSE, 
  allow.missing.Y = FALSE, allow.zero.sd = FALSE, intervals = NULL) {
  tuning <- tune.tg.lasso(X, Y, lambda.1.range, lambda.2.range, 
    lambda.3.range, k, allow.missing.X = allow.missing.X, 
    allow.missing.Y = allow.missing.Y, allow.zero.sd = allow.zero.sd, 
    intervals = intervals)
  results <- data.frame(feature = ranked.features, mean = rep(NA, 
    length(ranked.features)), sd = rep(NA, length(ranked.features)))
  for (i in 1:length(ranked.features)) {
    feature <- ranked.features[i]
    j <- which(colnames(X) == feature)
    cv <- cross.validate.tg.lasso(cbind(X[,-j]), Y, tuning$lambda.1.min, 
        tuning$lambda.2.min, tuning$lambda.3.min, k, 
        allow.missing.X = allow.missing.X, allow.missing.Y = allow.missing.Y, 
        allow.zero.sd = allow.zero.sd, intervals = intervals)
    results[i,2:3] <- c(cv$rmse.mean, cv$rmse.sd)
  }
  min.index <- which.min(results$mean)
  min.feature <- ranked.features[min.index]
  min.mean <- results$mean[min.index]
  min.sd <- results$sd[min.index]
  max.index <- which.max(results$mean)
  max.feature <- ranked.features[max.index]
  max.mean <- results$mean[max.index]
  max.sd <- results$sd[max.index]
  ranking <- list(X, Y, results = results, min.index = min.index, 
    min.feature = min.feature, min.mean = min.mean, min.sd = min.sd, 
    max.index = max.index, max.feature = max.feature, max.mean = max.mean, 
    max.sd = max.sd, results = results)
  return(ranking)
}

#'
#' Backward Selection
#'
#' @export
#'
backward.selection <- function(X, Y, lambda.1.range, lambda.2.range, 
  lambda.3.range, k, allow.missing.X = FALSE, allow.missing.Y = FALSE, 
  allow.zero.sd = FALSE, intervals = NULL) {
  initial.features <- colnames(X)
  results <- data.frame(feature = rep(NA, ncol(X)), mean = rep(NA, ncol(X)), 
    sd = rep(NA, ncol(X)))
  p <- ncol(X)
  for (i in 1:(p - 1)) {
    ranking <- backward.ranking(X, Y, colnames(X), lambda.1.range, 
      lambda.2.range, lambda.3.range, k, allow.missing.X = allow.missing.X, 
      allow.missing.Y = allow.missin.Y, allow.zero.sd = allow.zero.sd, 
      intervals = intervals)
    if (i < p) X <- cbind(X[,-ranking$min.index])
    results[i,1] <- ranking$min.feature
    results[i,2:ncol(results)] <- c(ranking$min.mean, ranking$min.sd)
    cat(sprintf("Selected feature: '%s' - Value: %s (+/- %s).\n", 
      ranking$min.feature, formatC(ranking$min.mean), formatC(ranking$min.sd)))
  }
  results[nrow(results),1] <- setdiff(initial.features, results$feature)
  min.index <- which.min(results$mean)
  min.feature.set <- results$feature[(min.index + 1):nrow(results)]
  min.mean <- results$mean[min.index]
  min.sd <- results$sd[min.index]
  selection <- list(min.index = min.index, min.feature.set = min.feature.set, 
    min.mean = min.mean, min.sd = min.sd, results = results)
  return(selection)
}
