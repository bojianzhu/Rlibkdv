#' @title Use KDV
#' @description
#' Efficient and accurate kernel density visualization.
#'
#' @param longitude features' longitude
#' @param latitude features' latitude
#' @param bandwidth_s spatial bandwidth
#' @param row_pixels row pixels
#' @param col_pixels col pixels
#' @return kdv result
#' @importFrom utils read.csv
#' @examples
#' \donttest{
#' data(hk)
#' resKDV <- kdv(hk$lon, hk$lat, 1000, 800 ,640)
#' }
#' @export
kdv <- function(longitude, latitude, bandwidth_s =1000, row_pixels=800, col_pixels=640) {
  type <- "KDV"
  num_threads <- 8
  GPS <- TRUE
  kernel_s_type <- 1
  omit_zero <- 1
  middle_lat <- NULL
  bound <- vector("numeric", length = 4)
  data <- data.frame(latitude, longitude)
  colnames(data) <- c("lat", "lon")
  datas <- set_data(data, GPS, middle_lat, type)
  data <- datas[[1]]
  data_str <- datas[[2]]
  middle_lat <- datas[[3]]
  min_x <- datas[[4]]
  min_y <- datas[[5]]
  bound <- set_bound(data, GPS, bound, middle_lat, min_x, min_y)
  KDV_type <- ifelse(type == 'STKDV', 3, 1)
  args <- c(data_str, KDV_type, num_threads, bound[1], bound[2], bound[3], bound[4], row_pixels, col_pixels, kernel_s_type, bandwidth_s, 0, 0, 32, 1, 6, omit_zero)
  kdv_res <- kdvCpp(args)
  col_names <- c("x", "y", "val")
  result <- read.csv(textConnection(kdv_res), header = FALSE, col.names = col_names)
  result <- unshift_GPS(result, min_x, min_y)
  if (GPS) {
    result <- XY_to_GPS(result, middle_lat)
    result_cols <- c('lon', 'lat')
  } else {
    result_cols <- c('x', 'y')
  }
  result_cols <- c(result_cols, 'val')
  result <- result[, result_cols]
  return(result)
}
#' @title Use STKDV
#' @description
#' Efficient and accurate spatiotemporal kernel density visualization.
#'
#' @param longitude features' longitude
#' @param latitude features' latitude
#' @param time  features' time
#' @param bandwidth_s spatial bandwidth
#' @param bandwidth_t temporal bandwidth
#' @param row_pixels row pixels
#' @param col_pixels col pixels
#' @param t_pixels time pixels
#' @importFrom utils read.csv
#' @return stkdv result
#' @examples
#' \donttest{
#' data(hk)
#' resSTKDV <- stkdv(hk$lon, hk$lat, hk$t, 1000, 6, 800, 640, 32)
#' }
#' @export
stkdv <- function(longitude, latitude, time, bandwidth_s=1000, bandwidth_t=6, row_pixels=800, col_pixels=640,  t_pixels=32) {
  type <- "STKDV"
  num_threads <- 8
  GPS <- TRUE
  kernel_s_type <- 1
  kernel_t_type <- 1
  omit_zero <- 1
  middle_lat <- NULL
  bound <- vector("numeric", length = 4)
  t_bound <- vector("numeric", length = 2)
  data <- data.frame(latitude, longitude, time)
  colnames(data) <- c("lat", "lon", "t")
  datas <- set_data(data, GPS, middle_lat, type)
  data <- datas[[1]]
  data_str <- datas[[2]]
  middle_lat <- datas[[3]]
  min_x <- datas[[4]]
  min_y <- datas[[5]]
  bound <- set_bound(data, GPS, bound, middle_lat, min_x, min_y)
  min_t <- datas[[6]]
  t_bound <- set_t_bound(data, t_bound, min_t)
  KDV_type <- ifelse(type == 'STKDV', 3, 1)
  args <- c(data_str, KDV_type, num_threads, bound[1], bound[2], bound[3], bound[4], row_pixels, col_pixels, kernel_s_type, bandwidth_s, t_bound[1], t_bound[2], t_pixels, kernel_t_type, bandwidth_t, omit_zero)
  kdv_res <- kdvCpp(args)
  col_names <- c("x", "y", "t", "val")
  result <- read.csv(textConnection(kdv_res), header = FALSE, col.names = col_names)
  result <- unshift_GPS(result, min_x, min_y)
  result <- unshift_time(result, min_t)
  if (GPS) {
    result <- XY_to_GPS(result, middle_lat)
    result_cols <- c('lon', 'lat')
  } else {
    result_cols <- c('x', 'y')
  }
  result_cols <- c(result_cols, 'val')
  result_cols <- c(result_cols, 't')

  result <- result[, result_cols]
  result <- result[order(result$lon, result$t), ]
  rownames(result) <- NULL
  return(result)
}

#' @title Plot KDV
#'
#' @param data result of kdv
#' @return No return value, called to plot KDV heatmap
#' @export
#' @import leaflet
#' @import raster
#' @import sf
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' data(hk)
#' resKDV <- kdv(hk$lon, hk$lat, 1000, 800 ,640)
#' plotKDV(resKDV)
#' }

plotKDV <- function(data) {
  r <- raster::rasterFromXYZ(data)
  raster::crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
  pal <- leaflet::colorNumeric("Reds", raster::values(r), na.color = "transparent")
  leaflet::leaflet() %>%
    leaflet::addTiles() %>%
    leaflet::addRasterImage(r, colors = pal, opacity = 0.75)
}

#' @title Plot STKDV
#'
#' @param data result of stkdv
#' @return No return value, called to plot STKDV heatmap
#' @export
#' @import leaflet
#' @import raster
#' @import sf
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' data(hk)
#' resSTKDV <- stkdv(hk$lon, hk$lat, hk$t, 1000, 6, 800, 640, 32)
#' plotSTKDV(resSTKDV)
#' }

plotSTKDV <- function(data) {
  data$lon <- round(data$lon, 3)
  data$lat <- round(data$lat, 3)
  data_list <- split(data, data$t)
  names(data_list) <- format(as.POSIXct(as.numeric(names(data_list)), origin = "1970-01-01"), "%Y-%m-%d %H-%M-%S")
  r_list <- lapply(data_list, function(data) {
    r <- raster::rasterFromXYZ(data[, c("lon", "lat", "val")])
    raster::crs(r) <- "+proj=longlat +datum=WGS84 +no_defs"
    return(r)
  })
  domain = range(data$val, na.rm = TRUE)
  pal <- leaflet::colorNumeric("Reds", domain , na.color = "transparent")
  m <- leaflet::leaflet() %>% leaflet::addTiles()
  for (i in seq_along(r_list)) {
    m <- m %>% leaflet::addRasterImage(r_list[[i]], colors = pal, opacity = 0.75, group = names(data_list)[[i]])
  }
  map <- m %>% leaflet::addLayersControl(baseGroups = names(data_list))

  map
}

