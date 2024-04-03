# Reimplementation of the `smoothn` function in R
# https://www.biomecardio.com/matlab/smoothn.m
# https://github.com/SanParraguez/smoothn/blob/main/smoothn.py

# dctn --------------------------------------------------------------------

#' Multidimensional Discrete Cosine Transform along one or more axes
#'
#' This function uses the seperability property of n-D Discrete Consine
#' Transform (DCT) and apply the one dimensional DCT/IDCT along each specified
#' axis. See [gsignal::dct()] and [gsignal::dct2()] for 1-D DCT and 2-D DCT
#' respectively. See [gsignal::idct()] and [gsignal::idct2()] for
#' 1-D IDCT and 2-D IDCT respectively.
#'
#' @param x Array. A n-dimensional array.
#' @param d Integer. A vector of dimensions for applying the 1-D DCT/IDCT. If
#' `d = NULL`, all dimensions will be used.
#' @returns An array of the same shape as `x`. If `x` is a vector, the
#' output will also be a vector.
#' @examples
#'
#' # Apply DCT for rows and columns
#' dctn(cars)
#'
#' # Apply DCT for rows
#' dctn(cars, 2)
#'
#' # Reconstruct `cars`
#' idctn(dctn(cars))
#'
#' @name dctn
NULL

#' @describeIn dctn Discrete Cosine Transform
#' @export
dctn <- function(x, d = NULL) {

  if (is.vector(x)) return(gsignal::dct(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  x <- unname(x)

  if (length(dim(x)) == 1 && is.null(d)) return(gsignal::dct(x))
  if (length(dim(x)) == 2 && is.null(d)) return(gsignal::dct2(x))

  z <- x

  if (is.null(d)) d <- 1:length(dim(z))

  for (this_d in d) {
    z <- replace_d(z, this_d, gsignal::dct)
  }

  return(z)
}


# idctn -------------------------------------------------------------------

#' @describeIn dctn Inverse Discrete Cosine Transform
#' @export
idctn <- function(x, d = NULL) {

  if (is.vector(x)) return(gsignal::idct(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  x <- unname(x)

  if (length(dim(x)) == 1 && is.null(d)) return(gsignal::idct(x))
  if (length(dim(x)) == 2 && is.null(d)) return(gsignal::idct2(x))

  z <- x

  if (is.null(d)) {
    d <- 1:length(dim(z))
  } else {
    d <- unique(d)
  }

  for (this_d in d) {
    z <- replace_d(z, this_d, gsignal::idct)
  }

  return(z)
}


# slice_d -----------------------------------------------------------------

#' Get indices of vectors of array along an axis
#'
#' This function slices all vectors of an array along an axis.
#'
#' @param x Array. A n-dimensional array.
#' @param d Integer. The target dimension. e.g. 1 for slicing along row, 2 for
#' slicing along column.
#' @param get_data Boolean. Whether to retrieve data using the indices.
#' @returns A tibble with two or three columns. The first column is `vector_id`
#' providing a vector of unique IDs for each vector. The second column is
#' `index` providing a list of matrices for accessing each vector. The third
#' column is `data` if `get_data = TRUE` providing a list of vectors.
#'
#' @examples
#'
#' # Slice along row. For 2D array, returns all columns.
#' slice_d(as.matrix(cars), d = 1)
#'
#' # Slice along column. For 2D array, returns all rows.
#' slice_d(as.matrix(cars), d = 2)
#'
#' # Slice and access.
#' slice_d(as.matrix(cars), d = 1, get_data = TRUE)
#'
#' @export
slice_d <- function(x, d, get_data = FALSE) {

  if (is.data.frame(x)) x <- as.matrix(x)
  x <- unname(x)
  if (length(dim(x)) < d) stop("Dimension indexing out of bound!")
  if (length(dim(x)) == 1) return(tibble::tibble(set = 1,
                                                 index = list(t(matrix(1:dim(x)[1])))
                                                 )
                                  )

  # Shift the d dimension to the first position
  all_dim <- dim(x)
  all_dim <- c(all_dim[d], all_dim[-d])

  result <- do.call(expand.grid, lapply(all_dim, function(x) 1:x)) |>
    as.matrix() |>
    unname()

  result <- result[, order(c(d, (1:length(dim(x)))[-d]))]
  colnames(result) <- paste0("V", 1:length(dim(x)))

  result <- result |>
    tibble::as_tibble() |>
    dplyr::mutate(vector_id = rep(1:(dplyr::n()/dim(x)[d]), each = dim(x)[d])) |>
    tidyr::nest(.by = vector_id, .key = "index") |>
    dplyr::mutate(index = lapply(index, function(x) unname(as.matrix(x))))

  if (get_data) {
    result <- result |>
      dplyr::mutate(data =  lapply(index, function(idx) x[idx]))
  }

  return(result)
}


# replace_d ---------------------------------------------------------------

#' Replace elements in vectors of array along an axis
#'
#' This function replaces elements in vectors of an array along an axis.
#' See [slice_d()] for more details about slicing vectors along an axis.
#'
#' @param x Array. A n-dimensional array.
#' @param d Integer. The target dimension. e.g. 1 for slicing and replacing
#' along row, 2 for slicing and replacing along column.
#' @param fn Function. A function that takes an original vector as input and
#' outputs a new vector of the same length.
#' @param ... Additional arguments passed to `fn`.
#' @returns An array of the same shape as `x`.
#' @examples
#'
#' # Replace columns with column sums
#' replace_d(as.matrix(cars), d = 1, sum)
#'
#' # Replace rows with row sums
#' replace_d(as.matrix(cars), d = 2, sum)
#' @export
replace_d <- function(x, d, fn, ...) {

  if (is.data.frame(x)) x <- as.matrix(x)
  x <- unname(x)
  idx <- slice_d(x, d)

  z <- x
  for (i in 1:nrow(idx)) {
    z[idx$index[[i]]] <- fn(x[idx$index[[i]]], ...)
  }

  return(z)
}

# initial_guess -----------------------------------------------------------

initial_guess <- function(x) {

  z <- x

  # Replace missing values with nearest non-missing values
  if (any(is.na(x))) {
      na_index <- which(is.na(x), arr.ind = TRUE)
      avi_index <- which(!is.na(x), arr.ind = TRUE)
      for (i in 1:nrow(na_index)) {
        index_dist <- apply((na_index[i, ] - avi_index)^2, 1, sum)
        target_index <- avi_index[which.min(index_dist), ]
        target_value <- x[t(unname(target_index))]
        z[t(unname(na_index[i, ]))] <- target_value
      }
  }

  # TODO: coarse fast smoothing using one-tenth of the DCT coefficients

  return(z)
}


# smoothn -----------------------------------------------------------------

#' Robust spline smoothing for one and higher dimensions data
#'
#' This function implements an algorithm for fast and robust smoothing of data
#' in one and higher dimensions based on a penalized least squares method.
#' It utilizes the discrete cosine transform for efficient processing.
#' The algorithm automatically determines the amount of smoothing required
#' by minimizing the generalized cross-validation score.
#'
#' @param x Array. A gridded numerical dataset.
#' @param s Character/Numeric. A positive numerical value that controls the
#' degree of smoothing. A greater value will result in a smoother output. If
#' `s = "auto"`, the smoothing parameter will be chosen by minimizing the
#' generalized cross-validation score (GCV) score.
#' @param weight Data frame/Matrix/Array. A weight matrix for smoothing.
#' It should have the same size as `x`. Use equal weights by default.
#' @param tol Numeric. A positive tolerance value for determining if the
#' prediction converges.
#' @param initial_guess Data frame/Matrix/Array. Initial prediction for the
#' iterative process. Use a zero matrix by default. If `x` contains any
#' missing data or the provided weight matrix has unequal weights,
#' `[initial_guess()]` will be used to compute the initial prediction.
#' @param spacing Numeric. A vector of positive numeric values controls the
#' spacing between points in each dimension. Length of the vector should be the
#' same as the number of columns of `x`. Assume equal spacing in all dimensions
#' by default.
#' @param robust Boolean. Whether robust smoothing should be applied to
#' minimize or cancel the side effects of high leverage points and outliers.
#' @param robust_iter Integer. Number of iterations used in robust smoothing.
#' @param robust_inner_iter Integer. The maximum number of iterations allowed for
#' updating the prediction in one robust iteration.
#' @param robust_method Character. Weight function for robust smoothing.
#' One of "bisquare", "talworth" and "cauchy".
#' @returns test
#' @examples
#'
#' c("test")
#'
#' @export
smoothn <- function(x,
                    s = "auto",
                    weight = array(1, dim(x)),
                    tol = 1e-3,
                    initial_guess = array(0, dim(x)),
                    spacing = rep(1, dim(x)[2]),
                    robust = TRUE,
                    robust_iter = 3L,
                    robust_inner_iter = 1000L,
                    robust_method = "bisquare"
                    ) {

}
