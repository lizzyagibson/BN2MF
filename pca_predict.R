# PCA predict works

myprcomp <- function (x, ...) UseMethod("myprcomp")

myprcomp.default <- function(x, retx = TRUE, center = TRUE, scale. = FALSE, tol = NULL, ...)
{
  chkDots(...)
  x <- as.matrix(x)
  x <- scale(x, center = center, scale = scale.)
  cen <- attr(x, "scaled:center")
  sc <- attr(x, "scaled:scale")
  if(any(sc == 0))
    stop("cannot rescale a constant/zero column to unit variance")
  s <- svd(x)
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


myprcomp.formula <- function (formula, data = NULL, subset, na.action, ...)
{
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

summary.myprcomp <- function(object, ...)
{
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

print.summary.myprcomp <-
  function(x, digits = max(3L, getOption("digits") - 3L), ...)
  {
    cat("Importance of components:\n")
    print(x$importance, digits = digits, ...)
    invisible(x)
  }

predict.myprcomp <- function(object, newdata, ...)
{
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

y <- scale(rmatrix(20, 10))
#y <- scale(matrix(c(1,0,0,1),nrow = 2, ncol = 2))
#y <- (matrix(c(1,0,0,1),nrow = 2, ncol = 2))

dim(y)
pca_y <- myprcomp(y)

# PCA scores
pca_y$x
dim(pca_y$x)

# predicted scores
predict(pca_y, y)
mean((pca_y$x - predict(pca_y, y))^2)

y_s <- scale(y, center = TRUE, scale = TRUE)

pred_x1 <- y %*% pca_y$rotation
pred_x2 <- y %*% svd(y )$v  # need an extra step to deal with rank deficient / nonsquare matrices

pca_y$x
pred_x1
pred_x2

mean((pca_y$x - pred_x1)^2)
cor(pca_y$x, pred_x1)

mean((pca_y$x - pred_x2)^2)
cor(pca_y$x, pred_x2)