pi <- 3.14159265358979323846
earth_radius <- 6371000


GPS_bound_to_XY <- function(bound, middle_lat) {
  middle_lat <- middle_lat * pi / 180
  bound[1] <- earth_radius * bound[1] * pi / 180 * cos(middle_lat)
  bound[2] <- earth_radius * bound[2] * pi / 180 * cos(middle_lat)
  bound[3] <- earth_radius * bound[3] * pi / 180
  bound[4] <- earth_radius * bound[4] * pi / 180
}

shift_bound_GPS <- function(bound, min_x, min_y) {
  bound[1] <- bound[1] + (-min_x)
  bound[2] <- bound[2] + (-min_x)
  bound[3] <- bound[3] + (-min_y)
  bound[4] <- bound[4] + (-min_y)
}

GPS_to_XY <- function(data, middle_lat) {
  middle_lat <- middle_lat * pi / 180
  data$x <- earth_radius * data$lon * pi / 180 * cos(middle_lat)
  data$y <- earth_radius * data$lat * pi / 180
  return(data)
}

XY_to_GPS <- function(result, middle_lat) {
  middle_lat <- middle_lat * pi / 180
  result$lon <- (result$x / (earth_radius * cos(middle_lat))) * (180 / pi)
  result$lat <- (result$y / earth_radius) * (180 / pi)
  return(result)
}

shift_GPS <- function(data) {
  min_x <- min(data$x, na.rm = TRUE)
  data$x <- data$x + (-min_x)
  min_y <- min(data$y, na.rm = TRUE)
  data$y <- data$y + (-min_y)
  return(list(data, min_x, min_y))
}

unshift_GPS <- function(data, min_x, min_y) {
  data$x <- data$x + min_x
  data$y <- data$y + min_y
  return(data)
}

shift_time <- function(data) {
  if (!("t" %in% names(data))) {
    return(0)
  }
  min_t <- min(data$t)
  data$t <- data$t + (-min_t)
  data$t <- data$t / 86400
  return(list(data, min_t))
}

unshift_time <- function(data, min_t) {
  if (!("t" %in% names(data))) {
    return(0)
  }
  data$t <- data$t * 86400
  data$t <- data$t + min_t
  return(data)
}
#' @importFrom utils capture.output
#' @importFrom utils write.csv
set_data <- function(data, GPS, middle_lat, KDV_type) {
  if (is.null(data)) {
    return()
  }
  if (GPS) {
    if (is.null(middle_lat)) {
      if (!("lat" %in% colnames(data))) {
        middle_lat <- 90
      } else {
        middle_lat <- (min(data$lat) + max(data$lat)) / 2
      }
    }
    if (!("lat" %in% colnames(data))) {
      data <- XY_to_GPS(data, middle_lat)
    } else {
      data <- GPS_to_XY(data, middle_lat)
    }
  }
  datas <- shift_GPS(data)
  data <- datas[[1]]
  min_x <- datas[[2]]
  min_y <- datas[[3]]
  if (KDV_type == "STKDV") {
    datas <- shift_time(data)
    data <- datas[[1]]
    min_t <- datas[[2]]
  }
  if (!("w" %in% colnames(data))) {
    data$w <- 1
  }
  if (KDV_type == "KDV") {
    data_str <- paste(capture.output(write.csv(data[, c("x", "y", "w")], row.names = FALSE)), collapse = "\n")
  } else {
    data_str <- paste(capture.output(write.csv(data[, c("x", "y", "t", "w")], row.names = FALSE)), collapse = "\n")
  }
  if (KDV_type == "STKDV") {
    return(list(data, data_str, middle_lat, min_x, min_y, min_t))
  }
  return(list(data, data_str, middle_lat, min_x, min_y))
}

set_bound <- function(data, GPS, bound, middle_lat, min_x, min_y) {
  tryCatch({
    if (length(bound) != 4) {
      bound <- NULL
    }
    if (all(bound == c(0, 0, 0, 0))) {
      bound <- NULL
    }
  }, error = function(e) {
    bound <- NULL
  })
  if (is.null(bound)) {
    bound <- c(min(data$x, na.rm = TRUE), max(data$x, na.rm = TRUE), min(data$y, na.rm = TRUE), max(data$y, na.rm = TRUE))
  } else {
    bound <- c(bound[3], bound[4], bound[1], bound[2])
    if (GPS) {
      GPS_bound_to_XY(bound, middle_lat)
      shift_bound_GPS(bound, min_x, min_y)
    }
  }
  return(bound)
}

set_t_bound <- function(data, t_bound, min_t) {
  tryCatch({
    if (length(t_bound) != 2) {
      t_bound <- NULL
    }
    if (all(t_bound == c(0, 0))) {
      t_bound <- NULL
    }
  }, error = function(e) {
    t_bound <- NULL
  })

  if (is.null(t_bound)) {
    t_bound <- c(min(data$t, na.rm = TRUE), max(data$t, na.rm = TRUE))
  } else {
    t_bound[1] <- t_bound[1] - min_t
    t_bound[1] <- t_bound[1] / 86400
    t_bound[2] <- t_bound[2] - min_t
    t_bound[2] <- t_bound[2] / 86400
  }
  return(t_bound)
}
