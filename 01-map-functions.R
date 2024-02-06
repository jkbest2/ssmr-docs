read_wbd <- function(hu = "2", wbd_dir = "data/WBD_17_HU2_Shape") {
    fn <- paste0("WBDHU", hu, ".shp")
    path <- file.path(wbd_dir, "Shape", fn)
    read_sf(path)
}

check_nhd_feature <- function(feature) {
    nhd_features <- c(
            "Area",
            "AreaEventFC",
            "Flowline",
            "Line",
            "LineEventFC",
            "Point",
            "PointEventFC",
            "Waterbody")
    if (!(feature %in% nhd_features)) {
            stop(paste("feature must be one of",
             paste(nhd_features, collapse = ", ")))
    }
    invisible(TRUE)
}

read_nhd <- function(feature, nhd_dir = "data/NHD_H_1708_HU4_Shape") {
    check_nhd_feature(feature)
    fn <- paste0("NHD", feature, ".shp")
    path <- file.path(nhd_dir, "Shape", fn)
    read_sf(path)
}

read_dnr_hydro <- function(dnr_dir = "data/DNR_Hydrography_-_Water_Bodies") {
    read_sf(file.path(dnr_dir, "DNR_Hydrography_-_Water_Bodies_-_Forest_Practices_Regulation.shp"))
}

prep_geo <- function(geo_data, filter_poly = NULL) {
    updata <- st_zm(geo_data)
    if (!is.null(filter_poly)) {
        updata <- updata |>
            st_transform(st_crs(filter_poly)) |>
            st_filter(filter_poly, .predicate = st_within)
    }
    updata
}
