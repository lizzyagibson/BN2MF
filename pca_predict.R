# PCA predict works

myprcomp <- function (x, ...) UseMethod("myprcomp")

myprcomp.default <- function(x, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL, ...) {
  chkDots(...)
  x <- as.matrix(x)
  x <- scale(x, center = center, scale = scale.)
  cen <- attr(x, "scaled:center")
  sc <- attr(x, "scaled:scale")
  if(any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance")
  s <- svd(x) # INSTEAD OF NU = 0
  s$d <- s$d / sqrt(max(1, nrow(x) - 1))
  if (!is.null(tol)) {
    ## we get rank at least one even for a 0 matrix.
    rank <- sum(s$d > (s$d[1L]*tol))
    if (rank < ncol(x)) {
      s$v <- s$v[, 1L:rank, drop = FALSE]
      s$d <- s$d[1L:rank]
    }
  }
  dimnames(s$v) <-
    list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
  r <- list(sdev = s$d, rotation = s$v,
            center = if(is.null(cen)) FALSE else cen,
            scale = if(is.null(sc)) FALSE else sc)
  if (retx) r$x <- x %*% s$v
  class(r) <- "prcomp"
  r
}


myprcomp.formula <- function (formula, data = NULL, subset, na.action, ...){
  mt <- terms(formula, data = data)
  if (attr(mt, "response") > 0L)
    stop("response not allowed in formula")
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$... <- NULL
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  ## this is not a `standard' model-fitting function,
  ## so no need to consider contrasts or levels
  if (.check_vars_numeric(mf))
    stop("PCA applies only to numerical variables")
  na.act <- attr(mf, "na.action")
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 0L
  x <- model.matrix(mt, mf)
  res <- prcomp.default(x, ...)
  ## fix up call to refer to the generic, but leave arg name as `formula'
  cl[[1L]] <- as.name("prcomp")
  res$call <- cl
  if (!is.null(na.act)) {
    res$na.action <- na.act
    if (!is.null(sc <- res$x))
      res$x <- napredict(na.act, sc)
  }
  res
}

plot.myprcomp <- function(x, main = deparse(substitute(x)), ...)
  screeplot.default(x, main = main, ...)

print.myprcomp <- function(x, print.x = FALSE, ...) {
  cat("Standard deviations:\n")
  print(x$sdev, ...)
  cat("\nRotation:\n")
  print(x$rotation, ...)
  if (print.x && length(x$x)) {
    cat("\nRotated variables:\n")
    print(x$x, ...)
  }
  invisible(x)
}

summary.myprcomp <- function(object, ...){
  chkDots(...)
  vars <- object$sdev^2
  vars <- vars/sum(vars)
  importance <- rbind("Standard deviation" = object$sdev,
                      "Proportion of Variance" = round(vars, 5),
                      "Cumulative Proportion" = round(cumsum(vars), 5))
  colnames(importance) <- colnames(object$rotation)
  object$importance <- importance
  class(object) <- "summary.myprcomp"
  object
}

print.summary.myprcomp <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("Importance of components:\n")
    print(x$importance, digits = digits, ...)
    invisible(x)
  }

predict.myprcomp <- function(object, newdata, ...) {
  chkDots(...)
  if (missing(newdata)) {
    if(!is.null(object$x)) return(object$x)
    else stop("no scores are available: refit with 'retx=TRUE'")
  }
  if(length(dim(newdata)) != 2L)
    stop("'newdata' must be a matrix or data frame")
  nm <- rownames(object$rotation)
  if(!is.null(nm)) {
    if(!all(nm %in% colnames(newdata)))
      stop("'newdata' does not have named columns matching one or more of the original columns")
    newdata <- newdata[, nm, drop = FALSE]
  } else {
    if(NCOL(newdata) != NROW(object$rotation) )
      stop("'newdata' does not have the correct number of columns")
  }
  ## next line does as.matrix
  scale(newdata, object$center, object$scale) %*% object$rotation
}

y <- rmatrix(20, 10)
y_scaled <- scale(y)
y_scaled_F <- scale(y, center = TRUE, scale = FALSE)
dim(y)

pca_y <- prcomp(y)
mypca_y <- myprcomp(y)

# predict scores
# y not scaled
pred_x1 <- y %*% pca_y$rotation
pred_x2 <- y %*% svd(y)$v  

# y scaled
pred_x1_s <- y_scaled %*% pca_y$rotation
pred_x2_s <- y_scaled %*% svd(y_scaled)$v 

# y centered
pred_x1_sF <- y_scaled_F %*% pca_y$rotation
pred_x2_sF <- y_scaled_F %*% svd(y_scaled_F)$v 

# y not scaled
mean((pca_y$x - pred_x1)^2)
mean((pca_y$x - pred_x2)^2)
mean((mypca_y$x - pred_x1)^2)
mean((mypca_y$x - pred_x2)^2)
# none work

# y scaled
mean((pca_y$x - pred_x1_s)^2)
mean((pca_y$x - pred_x2_s)^2)
mean((mypca_y$x - pred_x1_s)^2)
# THIS ONE WORKS
mean((mypca_y$x - pred_x2_s)^2)

# y centered
# FIRST TWO WORK
mean((pca_y$x - pred_x1_sF)^2)
mean((pca_y$x - pred_x2_sF)^2)
mean((mypca_y$x - pred_x1_sF)^2)
mean((mypca_y$x - pred_x2_sF)^2)

# Sanity check
mean((mypca_y$x - pca_y$x)^2)
