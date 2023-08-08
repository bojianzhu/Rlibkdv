#' @title Use KDV
#' @description
#' Efficient and accurate kernel density visualization.
#'
#' @param longitude features' longitude
#' @param latitude features' latitude
#' @param bandwidth_s spatial bandwidth
#' @param row_pixels row pixels
#' @param col_pixels col pixels
#'
#' @return kdv result
#' @export
#'
#' @examples
#' data(hk)
#' resKDV <-kdv(hk$lon, hk$lat, 1000.0, 800, 640)
kdv <- function(longitude, latitude, bandwidth_s = 1000.0, row_pixels = 800, col_pixels = 640) {
  data <- data.frame(latitude, longitude)
  colnames(data) <- c("lat", "lon")
  file_path <- system.file("libkdv/kdv.py", package = "Rlibkdv")
  reticulate::source_python(file_path)
  result <- use_kdv(data, bandwidth_s, row_pixels, col_pixels)
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
#'
#' @return stkdv result
#' @export
#'
#' @examples
#' data(hk)
#' resSTKDV <- stkdv(hk$lon, hk$lat, hk$t, 1000.0, 6.0, 800, 640, 32)
stkdv <- function(longitude, latitude, time, bandwidth_s = 1000.0 , bandwidth_t = 6.0, row_pixels = 800, col_pixels = 640, t_pixels = 32) {
  data <- data.frame(latitude, longitude, time)
  colnames(data) <- c("lat", "lon", "t")
  file_path <- system.file("libkdv/kdv.py", package = "Rlibkdv")
  reticulate::source_python(file_path)
  result <- use_stkdv(data, bandwidth_s, bandwidth_t, row_pixels, col_pixels, t_pixels)
  return(result)
}

#' @title Plot KDV
#'
#' @param data result of kdv
#'
#' @export
#'
#' @examples
#' data(hk)
#' resKDV <-kdv(hk$lon, hk$lat, 1000.0, 800, 640)
#' plotKDV(resKDV)
plotKDV <- function(data) {
  r <- raster::rasterFromXYZ(data)
  raster::crs(r) <- "EPSG:4326"
  pal <- leaflet::colorNumeric("Reds", raster::values(r), na.color = "transparent")
  leaflet::leaflet() %>%
    leaflet::addTiles() %>%
    leaflet::addRasterImage(r, colors = pal, opacity = 0.75)
}

#' @title Plot STKDV
#'
#' @param data result of stkdv
#'
#' @export
#'
#' @examples
#' data(hk)
#' resSTKDV <- stkdv(hk$lon, hk$lat, hk$t, 1000.0, 6.0, 800, 640, 32)
#' plotSTKDV(resSTKDV)
plotSTKDV <- function(data) {
  data_list <- split(data, data$t)
  names(data_list) <- format(as.POSIXct(as.numeric(names(data_list)), origin = "1970-01-01"), "%Y-%m-%d %H-%M-%S")
  r_list <- lapply(data_list, function(data) {
    data$lon <- round(data$lon, 3)
    data$lat <- round(data$lat, 3)
    r <- raster::rasterFromXYZ(data[, c("lon", "lat", "val")])
    raster::crs(r) <- "EPSG:4326"
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

generate_heatmap <- function(data) {
# Split the data frame into a list based on time
  data_list <- split(data, data$t)
  names(data_list) <- format(as.POSIXct(as.numeric(names(data_list)), origin = "1970-01-01"), "%Y-%m-%d %H-%M-%S")
# Create a heatmap for each element in the list
  domain = range(data$val, na.rm = TRUE)
  pal <- leaflet::colorNumeric("Reds", domain , na.color = "transparent")
  r_list <- lapply(data_list, function(data) {
    data$lon <- round(data$lon, 3)
    data$lat <- round(data$lat, 3)
    r <- raster::rasterFromXYZ(data[, c("lon", "lat", "val")])
    raster::crs(r) <- "EPSG:4326"
    return(r)
  })
  m <- leaflet::leaflet() %>% leaflet::addTiles()
  for (i in seq_along(r_list)) {
    m <- m %>% leaflet::addRasterImage(r_list[[i]], colors = pal, opacity = 0.75, group = names(data_list)[[i]])
  }

  mapWithTimeSlide <- htmlwidgets::onRender(
    m %>%
      registerPlugin(
        htmltools::htmlDependency(
          "leaflet.timedimension", "1.1.1",
          src = c(href = "https://cdn.jsdelivr.net/npm/leaflet-timedimension@1/"),
          script = "dist/leaflet.timedimension.min.js",
          stylesheet = "dist/leaflet.timedimension.control.min.css"
        )
      ) %>%
      registerPlugin(
        htmltools::htmlDependency(
          "ISO8601-js-period", "0.2.1",
          src = c(href = "https://cdn.jsdelivr.net/npm/iso8601-js-period@0/"),
          script = "iso8601.min.js"
        )
      ),
    sprintf(
      "function(el, x) {
        var timeDimension = new L.TimeDimension({times: [%s]});
        this.timeDimension = timeDimension;
        this.timeDimensionControl = new L.Control.TimeDimension({playerOptions: {buffer: 5}});
        this.addControl(this.timeDimensionControl);
      }",
      paste0(unique(data$t), collapse = ",")
    )
  )
  mapWithTimeSlide
}

registerPlugin <- function(map, plugin) {
  map$dependencies <- c(map$dependencies, list(plugin))
  map
}
