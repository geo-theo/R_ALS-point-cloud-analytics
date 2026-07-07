# Basic LiDAR point cloud cleaning and analytics prototype.
#
# Command line usage:
#   Rscript basic_lidar_analytics.R
#   Rscript basic_lidar_analytics.R --quick-test
#   Rscript basic_lidar_analytics.R path/to/input.laz path/to/output_dir
#   Rscript basic_lidar_analytics.R input.laz outputs/test --add-hag
#   Rscript basic_lidar_analytics.R input.laz outputs/test --read-filter="-keep_xy xmin ymin xmax ymax"
#
# RStudio Console usage:
#   source("basic_lidar_analytics.R")
#   run_quick_test()
#   run_default_lidar_analytics()
#
# Required packages:
#   install.packages(c("lidR", "data.table", "rlas"))

required_packages <- c("lidR", "data.table", "rlas")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them with: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(lidR))

quick_test_filter <- "-keep_xy 451350 4497900 451550 4498100"

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }

  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
}

usage <- function() {
  cat(
    paste(
      "Basic LiDAR point cloud cleaning and analytics",
      "",
      "Usage:",
      "  Rscript basic_lidar_analytics.R [input_las_or_laz] [output_dir] [options]",
      "",
      "Options:",
      "  --help                    Show this help message.",
      "  --no-write-cleaned        Skip writing the cleaned/preprocessed LAS/LAZ.",
      "  --add-hag                 Add height above ground (HAG) if ground points exist.",
      "  --drop-overlap            Drop LAS 1.4 overlap-flagged points when present.",
      "  --read-filter=FILTER      Pass a LASlib filter string to lidR::readLAS().",
      "  --z-trim=LOW,HIGH         Optional quantile trim on Z, e.g. 0.001,0.999.",
      "",
      "Default input:",
      "  USGS_LPC_NV_USFSR4_D23_1301_2064_EPSG-6341_LAS2023.laz",
      "",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  script_dir <- get_script_dir()

  config <- list(
    input_file = file.path(
      script_dir,
      "USGS_LPC_NV_USFSR4_D23_1301_2064_EPSG-6341_LAS2023.laz"
    ),
    output_dir = file.path(script_dir, "outputs", "lidar_basic_analytics"),
    read_filter = "",
    write_cleaned_laz = TRUE,
    add_hag = FALSE,
    drop_overlap = FALSE,
    z_trim = NULL
  )

  positional <- character()

  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      usage()
      quit(save = "no", status = 0)
    } else if (arg == "--no-write-cleaned") {
      config$write_cleaned_laz <- FALSE
    } else if (arg == "--add-hag") {
      config$add_hag <- TRUE
    } else if (arg == "--drop-overlap") {
      config$drop_overlap <- TRUE
    } else if (startsWith(arg, "--read-filter=")) {
      config$read_filter <- sub("^--read-filter=", "", arg)
    } else if (startsWith(arg, "--z-trim=")) {
      trim_value <- sub("^--z-trim=", "", arg)
      trim_parts <- as.numeric(strsplit(trim_value, ",", fixed = TRUE)[[1]])

      if (
        length(trim_parts) != 2 ||
          anyNA(trim_parts) ||
          trim_parts[[1]] < 0 ||
          trim_parts[[2]] > 1 ||
          trim_parts[[1]] >= trim_parts[[2]]
      ) {
        stop("--z-trim must look like --z-trim=0.001,0.999", call. = FALSE)
      }

      config$z_trim <- trim_parts
    } else if (startsWith(arg, "--")) {
      stop("Unknown option: ", arg, call. = FALSE)
    } else {
      positional <- c(positional, arg)
    }
  }

  if (length(positional) >= 1) {
    config$input_file <- positional[[1]]
  }

  if (length(positional) >= 2) {
    config$output_dir <- positional[[2]]
  }

  config$input_file <- normalizePath(config$input_file, mustWork = FALSE)
  config$output_dir <- normalizePath(config$output_dir, mustWork = FALSE)

  config
}

point_count <- function(las) {
  if (is.null(las) || is.null(las@data)) {
    return(0)
  }

  nrow(las@data)
}

append_preprocess_log <- function(log, step, before, after, note = "") {
  data.table::rbindlist(
    list(
      log,
      data.table::data.table(
        step = step,
        points_before = before,
        points_after = after,
        points_removed = before - after,
        note = note
      )
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

capture_las_check <- function(las) {
  tryCatch(
    utils::capture.output(lidR::las_check(las)),
    error = function(err) paste("las_check failed:", conditionMessage(err))
  )
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  max(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (length(x) < 2 || all(is.na(x))) {
    return(NA_real_)
  }

  stats::sd(x, na.rm = TRUE)
}

safe_quantile <- function(x, probs) {
  if (length(x) == 0 || all(is.na(x))) {
    return(rep(NA_real_, length(probs)))
  }

  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

numeric_summary <- function(x, attribute) {
  probs <- c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)
  quantiles <- safe_quantile(x, probs)

  data.table::data.table(
    attribute = attribute,
    statistic = c(
      "min",
      "p01",
      "p05",
      "p25",
      "median",
      "p75",
      "p95",
      "p99",
      "max",
      "mean",
      "sd"
    ),
    value = c(
      safe_min(x),
      quantiles[[1]],
      quantiles[[2]],
      quantiles[[3]],
      quantiles[[4]],
      quantiles[[5]],
      quantiles[[6]],
      quantiles[[7]],
      safe_max(x),
      safe_mean(x),
      safe_sd(x)
    )
  )
}

count_return_field <- function(dt, column_name) {
  counts <- dt[
    ,
    .(point_count = .N),
    by = .(value = as.character(get(column_name)))
  ][order(value)]

  counts[, field := column_name]
  data.table::setcolorder(counts, c("field", "value", "point_count"))
  counts
}

extract_last_epsg_from_wkt <- function(wkt_text) {
  if (is.null(wkt_text) || length(wkt_text) == 0 || is.na(wkt_text[[1]])) {
    return(NA_character_)
  }

  matches <- regmatches(
    wkt_text[[1]],
    gregexpr('ID\\["EPSG",[[:space:]]*[0-9]+\\]', wkt_text[[1]], perl = TRUE)
  )[[1]]

  if (length(matches) == 0 || matches[[1]] == "-1") {
    return(NA_character_)
  }

  sub('.*ID\\["EPSG",[[:space:]]*([0-9]+)\\].*', "\\1", matches[[length(matches)]])
}

class_lookup <- c(
  "0" = "Created, never classified",
  "1" = "Unclassified",
  "2" = "Ground",
  "3" = "Low vegetation",
  "4" = "Medium vegetation",
  "5" = "High vegetation",
  "6" = "Building",
  "7" = "Low point / noise",
  "8" = "Model key point",
  "9" = "Water",
  "10" = "Rail",
  "11" = "Road surface",
  "12" = "Overlap / reserved",
  "13" = "Wire guard",
  "14" = "Wire conductor",
  "15" = "Transmission tower",
  "16" = "Wire connector",
  "17" = "Bridge deck",
  "18" = "High noise"
)

write_table <- function(x, path) {
  data.table::fwrite(x, path, na = "")
  path
}

format_count <- function(x) {
  if (length(x) == 0 || is.na(x)) {
    return("")
  }

  format(round(as.numeric(x)), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_measure <- function(x, digits = 3) {
  if (length(x) == 0 || is.na(x)) {
    return("")
  }

  format(round(as.numeric(x), digits), big.mark = ",", scientific = FALSE, trim = TRUE)
}

summary_value <- function(overall_summary, metric) {
  idx <- match(metric, overall_summary$metric)

  if (is.na(idx)) {
    return(NA_character_)
  }

  overall_summary$value[[idx]]
}

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

format_report_cell <- function(x) {
  if (is.numeric(x)) {
    out <- ifelse(
      is.na(x),
      "",
      format(round(x, 4), big.mark = ",", scientific = FALSE, trim = TRUE)
    )
  } else {
    out <- as.character(x)
    out[is.na(out)] <- ""
  }

  out
}

table_to_html <- function(x, max_rows = 40) {
  if (is.null(x) || nrow(x) == 0) {
    return("<p>No data available.</p>")
  }

  x <- as.data.frame(x)
  total_rows <- nrow(x)

  if (total_rows > max_rows) {
    x <- utils::head(x, max_rows)
  }

  header <- paste0(
    "<tr>",
    paste0("<th>", html_escape(names(x)), "</th>", collapse = ""),
    "</tr>"
  )

  rows <- vapply(
    seq_len(nrow(x)),
    function(i) {
      cells <- vapply(x[i, , drop = FALSE], format_report_cell, character(1))
      paste0("<tr>", paste0("<td>", html_escape(cells), "</td>", collapse = ""), "</tr>")
    },
    character(1)
  )

  note <- if (total_rows > max_rows) {
    paste0("<p class=\"table-note\">Showing first ", max_rows, " of ", total_rows, " rows.</p>")
  } else {
    ""
  }

  paste0("<table>", header, paste(rows, collapse = ""), "</table>", note)
}

report_href <- function(path, output_dir) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  prefix <- paste0(normalized_dir, "/")

  if (startsWith(tolower(normalized_path), tolower(prefix))) {
    return(substr(normalized_path, nchar(prefix) + 1, nchar(normalized_path)))
  }

  basename(path)
}

sample_indices <- function(n, max_n, seed = 42L) {
  if (n <= max_n) {
    return(seq_len(n))
  }

  set.seed(seed)
  sort(sample.int(n, max_n))
}

sample_vector <- function(x, max_n, seed = 42L) {
  idx <- sample_indices(length(x), max_n, seed)
  x[idx]
}

plot_sample_note <- function(n, max_n) {
  if (n <= max_n) {
    return("")
  }

  paste0(
    " Plot rendering uses a deterministic sample of ",
    format_count(max_n),
    " points for speed; CSV summary statistics use all cleaned points."
  )
}

save_png_plot <- function(path, code, width = 1200, height = 800, res = 140) {
  code_expr <- substitute(code)

  tryCatch(
    {
      grDevices::png(path, width = width, height = height, res = res)
      on.exit(grDevices::dev.off(), add = TRUE)
      eval(code_expr, envir = parent.frame())
      path
    },
    error = function(err) {
      warning(
        "Could not create plot ",
        basename(path),
        ": ",
        conditionMessage(err),
        call. = FALSE
      )
      NULL
    }
  )
}

add_plot_manifest_row <- function(manifest, path, title, description) {
  if (is.null(path) || !file.exists(path)) {
    return(manifest)
  }

  data.table::rbindlist(
    list(
      manifest,
      data.table::data.table(
        title = title,
        description = description,
        path = path
      )
    ),
    use.names = TRUE
  )
}

write_report_plots <- function(las, summaries, preprocess_log, output_dir) {
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- data.table::data.table(
    title = character(),
    description = character(),
    path = character()
  )

  dt <- las@data
  plot_n <- nrow(dt)

  z_values <- sample_vector(dt$Z, max_n = 1000000)
  path <- save_png_plot(
    file.path(plot_dir, "elevation_histogram.png"),
    {
      graphics::par(mar = c(4.5, 4.5, 3.5, 1))
      graphics::hist(
        z_values,
        breaks = 80,
        col = "#4C78A8",
        border = "white",
        main = "Elevation Distribution",
        xlab = "Z elevation",
        ylab = "Point count"
      )
      graphics::grid(col = "#E6E8EC")
    }
  )
  manifest <- add_plot_manifest_row(
    manifest,
    path,
    "Elevation histogram",
    paste0("Distribution of cleaned Z elevation values.", plot_sample_note(plot_n, 1000000))
  )

  if ("Intensity" %in% names(dt)) {
    intensity_values <- sample_vector(dt$Intensity, max_n = 1000000, seed = 43L)
    path <- save_png_plot(
      file.path(plot_dir, "intensity_histogram.png"),
      {
        graphics::par(mar = c(4.5, 4.5, 3.5, 1))
        graphics::hist(
          intensity_values,
          breaks = 80,
          col = "#F28E2B",
          border = "white",
          main = "Intensity Distribution",
          xlab = "Intensity",
          ylab = "Point count"
        )
        graphics::grid(col = "#E6E8EC")
      }
    )
    manifest <- add_plot_manifest_row(
      manifest,
      path,
      "Intensity histogram",
      paste0("Distribution of cleaned intensity values.", plot_sample_note(plot_n, 1000000))
    )
  }

  if (!is.null(summaries$classification_summary)) {
    classes <- summaries$classification_summary
    labels <- paste0(classes$Classification, ": ", classes$class_description)
    path <- save_png_plot(
      file.path(plot_dir, "classification_counts.png"),
      {
        graphics::par(mar = c(4.5, 10, 3.5, 1))
        graphics::barplot(
          rev(classes$point_count),
          names.arg = rev(labels),
          horiz = TRUE,
          las = 1,
          col = "#59A14F",
          border = NA,
          main = "Point Counts by Classification",
          xlab = "Point count"
        )
        graphics::grid(nx = NULL, ny = NA, col = "#E6E8EC")
      }
    )
    manifest <- add_plot_manifest_row(
      manifest,
      path,
      "Classification counts",
      "Point counts by LAS classification after cleaning."
    )
  }

  if (!is.null(summaries$return_summary)) {
    return_groups <- split(summaries$return_summary, summaries$return_summary$field)
    path <- save_png_plot(
      file.path(plot_dir, "return_counts.png"),
      {
        graphics::par(mfrow = c(length(return_groups), 1), mar = c(4, 4.5, 3, 1))
        for (field_name in names(return_groups)) {
          group <- return_groups[[field_name]]
          graphics::barplot(
            group$point_count,
            names.arg = group$value,
            col = "#B07AA1",
            border = NA,
            main = paste("Counts by", field_name),
            xlab = field_name,
            ylab = "Point count"
          )
          graphics::grid(nx = NA, ny = NULL, col = "#E6E8EC")
        }
      },
      height = if (length(return_groups) > 1) 1000 else 700
    )
    manifest <- add_plot_manifest_row(
      manifest,
      path,
      "Return counts",
      "Counts by return number and number of returns."
    )
  }

  path <- save_png_plot(
    file.path(plot_dir, "preprocessing_removed_points.png"),
    {
      graphics::par(mar = c(4.5, 9, 3.5, 1))
      removed <- preprocess_log$points_removed
      labels <- gsub("_", " ", preprocess_log$step, fixed = TRUE)

      if (any(removed > 0, na.rm = TRUE)) {
        graphics::barplot(
          rev(removed),
          names.arg = rev(labels),
          horiz = TRUE,
          las = 1,
          col = "#E15759",
          border = NA,
          main = "Points Removed During Preprocessing",
          xlab = "Points removed"
        )
        graphics::grid(nx = NULL, ny = NA, col = "#E6E8EC")
      } else {
        graphics::plot.new()
        graphics::title("Points Removed During Preprocessing")
        graphics::text(0.5, 0.5, "No points were removed by preprocessing.")
      }
    }
  )
  manifest <- add_plot_manifest_row(
    manifest,
    path,
    "Preprocessing removals",
    "Point removals for each cleaning step."
  )

  if (all(c("X", "Y") %in% names(dt))) {
    idx <- sample_indices(nrow(dt), max_n = 250000, seed = 44L)
    x_values <- dt$X[idx]
    y_values <- dt$Y[idx]
    path <- save_png_plot(
      file.path(plot_dir, "xy_density_overview.png"),
      {
        graphics::par(mar = c(4.5, 4.5, 3.5, 1))
        graphics::smoothScatter(
          x_values,
          y_values,
          nrpoints = 0,
          colramp = grDevices::colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B")),
          xlab = "X",
          ylab = "Y",
          main = "XY Point Density Overview"
        )
        graphics::box()
      }
    )
    manifest <- add_plot_manifest_row(
      manifest,
      path,
      "XY density overview",
      paste0("Plan-view density overview of cleaned points.", plot_sample_note(plot_n, 250000))
    )
  }

  manifest
}

unique_output_path <- function(path) {
  if (!file.exists(path)) {
    return(path)
  }

  ext <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)
  suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")

  if (nzchar(ext)) {
    return(paste0(stem, "_", suffix, ".", ext))
  }

  paste0(path, "_", suffix)
}

get_las_epsg <- function(las) {
  epsg_value <- tryCatch(
    suppressWarnings(lidR::epsg(las)),
    error = function(err) NA_integer_
  )

  if (length(epsg_value) > 0 && !is.na(epsg_value) && epsg_value != 0) {
    return(as.character(epsg_value))
  }

  crs_epsg <- tryCatch(
    suppressWarnings(lidR::st_crs(las)$epsg),
    error = function(err) NA_integer_
  )

  if (length(crs_epsg) > 0 && !is.na(crs_epsg) && crs_epsg != 0) {
    return(as.character(crs_epsg))
  }

  wkt_text <- tryCatch(
    suppressWarnings(lidR::wkt(las)),
    error = function(err) NA_character_
  )

  extract_last_epsg_from_wkt(wkt_text)
}

get_las_crs <- function(las) {
  tryCatch(
    {
      crs_value <- suppressWarnings(lidR::st_crs(las))
      if (isTRUE(is.na(crs_value))) {
        return(NA_character_)
      }

      if (!is.null(crs_value$input) && nzchar(crs_value$input)) {
        return(crs_value$input)
      }

      as.character(crs_value)
    },
    error = function(err) NA_character_
  )
}

read_header_value <- function(header, name) {
  value <- header[[name]]

  if (is.null(value) || length(value) == 0) {
    return(NA_character_)
  }

  as.character(value[[1]])
}

preprocess_las <- function(las, config) {
  preprocess_log <- data.table::data.table(
    step = character(),
    points_before = numeric(),
    points_after = numeric(),
    points_removed = numeric(),
    note = character()
  )

  before <- point_count(las)
  las <- lidR::filter_poi(las, is.finite(X) & is.finite(Y) & is.finite(Z))
  after <- point_count(las)
  preprocess_log <- append_preprocess_log(
    preprocess_log,
    "drop_non_finite_xyz",
    before,
    after,
    "Removed records with non-finite X, Y, or Z values."
  )

  if ("Withheld" %in% names(las@data)) {
    before <- point_count(las)
    las <- lidR::filter_poi(las, !Withheld)
    after <- point_count(las)
    preprocess_log <- append_preprocess_log(
      preprocess_log,
      "drop_withheld",
      before,
      after,
      "Removed LAS points flagged as withheld."
    )
  } else if ("Withheld_flag" %in% names(las@data)) {
    before <- point_count(las)
    las <- lidR::filter_poi(las, !Withheld_flag)
    after <- point_count(las)
    preprocess_log <- append_preprocess_log(
      preprocess_log,
      "drop_withheld",
      before,
      after,
      "Removed LAS points flagged as withheld."
    )
  } else {
    preprocess_log <- append_preprocess_log(
      preprocess_log,
      "drop_withheld",
      point_count(las),
      point_count(las),
      "Skipped because no withheld flag column was found."
    )
  }

  if (isTRUE(config$drop_overlap)) {
    if ("Overlap" %in% names(las@data)) {
      before <- point_count(las)
      las <- lidR::filter_poi(las, !Overlap)
      after <- point_count(las)
      preprocess_log <- append_preprocess_log(
        preprocess_log,
        "drop_overlap",
        before,
        after,
        "Removed LAS 1.4 overlap-flagged points."
      )
    } else if ("Overlap_flag" %in% names(las@data)) {
      before <- point_count(las)
      las <- lidR::filter_poi(las, !Overlap_flag)
      after <- point_count(las)
      preprocess_log <- append_preprocess_log(
        preprocess_log,
        "drop_overlap",
        before,
        after,
        "Removed LAS 1.4 overlap-flagged points."
      )
    } else {
      preprocess_log <- append_preprocess_log(
        preprocess_log,
        "drop_overlap",
        point_count(las),
        point_count(las),
        "Skipped because no overlap flag column was found."
      )
    }
  }

  if ("Classification" %in% names(las@data)) {
    before <- point_count(las)
    las <- lidR::remove_noise(las)
    after <- point_count(las)
    preprocess_log <- append_preprocess_log(
      preprocess_log,
      "drop_noise_classes",
      before,
      after,
      "Removed points classified as low noise or high noise."
    )
  } else {
    preprocess_log <- append_preprocess_log(
      preprocess_log,
      "drop_noise_classes",
      point_count(las),
      point_count(las),
      "Skipped because no Classification column was found."
    )
  }

  before <- point_count(las)
  las <- lidR::filter_duplicates(las)
  after <- point_count(las)
  preprocess_log <- append_preprocess_log(
    preprocess_log,
    "drop_duplicate_xyz",
    before,
    after,
    "Removed points with duplicate XYZ coordinates."
  )

  if (!is.null(config$z_trim)) {
    z_limits <- safe_quantile(las@data$Z, config$z_trim)
    z_min <- z_limits[[1]]
    z_max <- z_limits[[2]]

    before <- point_count(las)
    las <- lidR::filter_poi(las, Z >= z_min & Z <= z_max)
    after <- point_count(las)
    preprocess_log <- append_preprocess_log(
      preprocess_log,
      "trim_z_quantiles",
      before,
      after,
      sprintf(
        "Kept Z values between quantiles %.6f and %.6f: %.3f to %.3f.",
        config$z_trim[[1]],
        config$z_trim[[2]],
        z_min,
        z_max
      )
    )
  }

  if (isTRUE(config$add_hag)) {
    ground_classes <- c(2L, 9L)
    has_ground <- "Classification" %in% names(las@data) &&
      any(las@data$Classification %in% ground_classes, na.rm = TRUE)

    if (has_ground) {
      before <- point_count(las)
      hag_result <- tryCatch(
        {
          lidR::height_above_ground(las, lidR::tin())
        },
        error = function(err) {
          attr(err, "hag_failed") <- TRUE
          err
        }
      )

      if (inherits(hag_result, "error")) {
        preprocess_log <- append_preprocess_log(
          preprocess_log,
          "add_height_above_ground",
          before,
          before,
          paste("Skipped because HAG interpolation failed:", conditionMessage(hag_result))
        )
      } else {
        las <- hag_result
        after <- point_count(las)
        preprocess_log <- append_preprocess_log(
          preprocess_log,
          "add_height_above_ground",
          before,
          after,
          "Added HAG attribute using a TIN surface from ground/water points."
        )
      }
    } else {
      preprocess_log <- append_preprocess_log(
        preprocess_log,
        "add_height_above_ground",
        point_count(las),
        point_count(las),
        "Skipped because no ground or water-classified points were found."
      )
    }
  }

  list(las = las, preprocess_log = preprocess_log)
}

build_summaries <- function(las, raw_header, raw_point_count, config) {
  dt <- las@data

  xmin <- safe_min(dt$X)
  xmax <- safe_max(dt$X)
  ymin <- safe_min(dt$Y)
  ymax <- safe_max(dt$Y)
  bbox_area <- (xmax - xmin) * (ymax - ymin)
  clean_point_count <- point_count(las)
  density <- if (is.finite(bbox_area) && bbox_area > 0) {
    clean_point_count / bbox_area
  } else {
    NA_real_
  }

  overall_summary <- data.table::data.table(
    metric = c(
      "input_file",
      "output_dir",
      "read_filter",
      "header_point_count",
      "loaded_point_count_raw",
      "cleaned_point_count",
      "points_removed",
      "epsg",
      "crs",
      "las_version",
      "point_data_format",
      "xmin",
      "xmax",
      "ymin",
      "ymax",
      "zmin_cleaned",
      "zmax_cleaned",
      "bbox_area_map_units2",
      "density_points_per_map_unit2"
    ),
    value = as.character(c(
      config$input_file,
      config$output_dir,
      config$read_filter,
      read_header_value(raw_header, "Number of point records"),
      raw_point_count,
      clean_point_count,
      raw_point_count - clean_point_count,
      get_las_epsg(las),
      get_las_crs(las),
      paste0(
        read_header_value(raw_header, "Version Major"),
        ".",
        read_header_value(raw_header, "Version Minor")
      ),
      read_header_value(raw_header, "Point Data Format ID"),
      xmin,
      xmax,
      ymin,
      ymax,
      safe_min(dt$Z),
      safe_max(dt$Z),
      bbox_area,
      density
    ))
  )

  numeric_tables <- list(numeric_summary(dt$Z, "Z_elevation"))

  if ("Intensity" %in% names(dt)) {
    numeric_tables <- c(numeric_tables, list(numeric_summary(dt$Intensity, "Intensity")))
  }

  if ("ScanAngle" %in% names(dt)) {
    numeric_tables <- c(numeric_tables, list(numeric_summary(dt$ScanAngle, "ScanAngle")))
  }

  if ("hag" %in% names(dt)) {
    numeric_tables <- c(numeric_tables, list(numeric_summary(dt$hag, "height_above_ground")))
  }

  numeric_attribute_summary <- data.table::rbindlist(numeric_tables, use.names = TRUE)

  classification_summary <- NULL

  if ("Classification" %in% names(dt)) {
    has_hag <- "hag" %in% names(dt)
    has_intensity <- "Intensity" %in% names(dt)

    classification_summary <- dt[
      ,
      .(
        point_count = .N,
        z_min = safe_min(Z),
        z_mean = safe_mean(Z),
        z_max = safe_max(Z),
        intensity_mean = if (has_intensity) safe_mean(Intensity) else NA_real_,
        hag_mean = if (has_hag) safe_mean(hag) else NA_real_
      ),
      by = Classification
    ][order(Classification)]

    classification_summary[
      ,
      percent := round(100 * point_count / sum(point_count), 4)
    ]
    classification_summary[
      ,
      class_description := unname(class_lookup[as.character(Classification)])
    ]
    classification_summary[
      is.na(class_description),
      class_description := "Unknown or user-defined"
    ]

    data.table::setcolorder(
      classification_summary,
      c(
        "Classification",
        "class_description",
        "point_count",
        "percent",
        "z_min",
        "z_mean",
        "z_max",
        "intensity_mean",
        "hag_mean"
      )
    )
  }

  return_summary <- list()

  if ("ReturnNumber" %in% names(dt)) {
    return_summary <- c(return_summary, list(count_return_field(dt, "ReturnNumber")))
  }

  if ("NumberOfReturns" %in% names(dt)) {
    return_summary <- c(
      return_summary,
      list(count_return_field(dt, "NumberOfReturns"))
    )
  }

  return_summary <- if (length(return_summary) > 0) {
    data.table::rbindlist(return_summary, use.names = TRUE)[order(field, value)]
  } else {
    NULL
  }

  list(
    overall_summary = overall_summary,
    numeric_attribute_summary = numeric_attribute_summary,
    classification_summary = classification_summary,
    return_summary = return_summary
  )
}

metric_card_html <- function(label, value, detail = "") {
  paste0(
    "<div class=\"metric-card\">",
    "<div class=\"metric-label\">", html_escape(label), "</div>",
    "<div class=\"metric-value\">", html_escape(value), "</div>",
    "<div class=\"metric-detail\">", html_escape(detail), "</div>",
    "</div>"
  )
}

plot_manifest_to_html <- function(plot_manifest, output_dir) {
  if (is.null(plot_manifest) || nrow(plot_manifest) == 0) {
    return("<p>No plots were created for this run.</p>")
  }

  paste(
    vapply(
      seq_len(nrow(plot_manifest)),
      function(i) {
        href <- report_href(plot_manifest$path[[i]], output_dir)
        title <- plot_manifest$title[[i]]
        description <- plot_manifest$description[[i]]

        paste0(
          "<figure>",
          "<a href=\"", html_escape(href), "\">",
          "<img src=\"", html_escape(href), "\" alt=\"", html_escape(title), "\">",
          "</a>",
          "<figcaption><strong>", html_escape(title), "</strong><br>",
          html_escape(description),
          "</figcaption>",
          "</figure>"
        )
      },
      character(1)
    ),
    collapse = "\n"
  )
}

output_file_list_to_html <- function(paths, output_dir) {
  if (length(paths) == 0) {
    return("<p>No output files were recorded.</p>")
  }

  items <- vapply(
    paths,
    function(path) {
      href <- report_href(path, output_dir)
      paste0(
        "<li><a href=\"",
        html_escape(href),
        "\">",
        html_escape(href),
        "</a></li>"
      )
    },
    character(1)
  )

  paste0("<ul class=\"file-list\">", paste(items, collapse = ""), "</ul>")
}

write_html_report <- function(
  report_path,
  plot_manifest,
  preprocess_log,
  summaries,
  config,
  generated_files
) {
  overall <- summaries$overall_summary

  loaded_points <- as.numeric(summary_value(overall, "loaded_point_count_raw"))
  cleaned_points <- as.numeric(summary_value(overall, "cleaned_point_count"))
  points_removed <- as.numeric(summary_value(overall, "points_removed"))
  density <- as.numeric(summary_value(overall, "density_points_per_map_unit2"))
  zmin <- as.numeric(summary_value(overall, "zmin_cleaned"))
  zmax <- as.numeric(summary_value(overall, "zmax_cleaned"))

  metric_cards <- paste(
    c(
      metric_card_html(
        "Loaded points",
        format_count(loaded_points),
        "Points read from the input file or read filter."
      ),
      metric_card_html(
        "Cleaned points",
        format_count(cleaned_points),
        "Points retained after preprocessing."
      ),
      metric_card_html(
        "Points removed",
        format_count(points_removed),
        "Total records removed by cleaning steps."
      ),
      metric_card_html(
        "EPSG",
        summary_value(overall, "epsg"),
        summary_value(overall, "crs")
      ),
      metric_card_html(
        "Elevation range",
        paste0(format_measure(zmin), " to ", format_measure(zmax)),
        "Cleaned Z values."
      ),
      metric_card_html(
        "Point density",
        format_measure(density),
        "Points per map-unit squared."
      )
    ),
    collapse = "\n"
  )

  read_filter <- summary_value(overall, "read_filter")
  if (!nzchar(read_filter)) {
    read_filter <- "None"
  }

  source_note <- data.table::data.table(
    field = c("Input file", "Output folder", "Read filter", "Generated at"),
    value = c(
      summary_value(overall, "input_file"),
      summary_value(overall, "output_dir"),
      read_filter,
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  )

  css <- paste(
    "body{font-family:Arial,sans-serif;line-height:1.45;margin:0;background:#f5f7fa;color:#17202a;}",
    "header{background:#1f2937;color:white;padding:32px 40px;}",
    "main{max-width:1180px;margin:0 auto;padding:28px 24px 48px;}",
    "section{background:white;border:1px solid #d9dee7;border-radius:8px;padding:22px;margin:0 0 22px;}",
    "h1,h2,h3{margin-top:0;} h1{font-size:30px;margin-bottom:8px;} h2{font-size:22px;}",
    "header .subtle{color:#D1D5DB;}",
    ".subtle{color:#596273;} .metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:14px;}",
    ".metric-card{border:1px solid #d9dee7;border-radius:8px;padding:14px;background:#fbfcfe;}",
    ".metric-label{font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:#596273;}",
    ".metric-value{font-size:24px;font-weight:700;margin-top:6px;word-break:break-word;}",
    ".metric-detail{font-size:13px;color:#596273;margin-top:5px;word-break:break-word;}",
    ".plot-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:18px;}",
    "figure{margin:0;border:1px solid #d9dee7;border-radius:8px;overflow:hidden;background:#fbfcfe;}",
    "figure img{display:block;width:100%;height:auto;background:white;}",
    "figcaption{padding:12px 14px;color:#384152;font-size:14px;}",
    "table{border-collapse:collapse;width:100%;font-size:14px;margin:8px 0 0;}",
    "th,td{border:1px solid #d9dee7;padding:8px 10px;text-align:left;vertical-align:top;}",
    "th{background:#eef2f7;} .table-wrap{overflow-x:auto;} .table-note{color:#596273;font-size:13px;}",
    ".file-list{columns:2;line-height:1.8;} code{background:#eef2f7;padding:2px 5px;border-radius:4px;}",
    "@media(max-width:760px){header{padding:24px;} .file-list{columns:1;}}",
    sep = "\n"
  )

  html <- paste0(
    "<!doctype html>\n",
    "<html lang=\"en\">\n",
    "<head>\n",
    "<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>LiDAR Point Cloud Analytics Report</title>\n",
    "<style>\n", css, "\n</style>\n",
    "</head>\n",
    "<body>\n",
    "<header>",
    "<h1>LiDAR Point Cloud Analytics Report</h1>",
    "<p class=\"subtle\">Automatic report for the cleaned point cloud and summary outputs.</p>",
    "</header>\n",
    "<main>\n",
    "<section>",
    "<h2>Run Overview</h2>",
    "<div class=\"metrics\">", metric_cards, "</div>",
    "<div class=\"table-wrap\">", table_to_html(source_note, max_rows = 10), "</div>",
    "</section>\n",
    "<section>",
    "<h2>Plots</h2>",
    "<p class=\"subtle\">Plots are generated from the cleaned point cloud. Large point clouds may be sampled for plotting speed only.</p>",
    "<div class=\"plot-grid\">", plot_manifest_to_html(plot_manifest, config$output_dir), "</div>",
    "</section>\n",
    "<section>",
    "<h2>Preprocessing</h2>",
    "<p>This table shows what each cleaning step did and how many points it removed.</p>",
    "<div class=\"table-wrap\">", table_to_html(preprocess_log), "</div>",
    "</section>\n",
    "<section>",
    "<h2>Classification Summary</h2>",
    "<div class=\"table-wrap\">", table_to_html(summaries$classification_summary), "</div>",
    "</section>\n",
    "<section>",
    "<h2>Numeric Attribute Summary</h2>",
    "<div class=\"table-wrap\">", table_to_html(summaries$numeric_attribute_summary), "</div>",
    "</section>\n",
    "<section>",
    "<h2>Return Summary</h2>",
    "<div class=\"table-wrap\">", table_to_html(summaries$return_summary), "</div>",
    "</section>\n",
    "<section>",
    "<h2>Generated Files</h2>",
    "<p>The CSV and TXT files remain available for detailed inspection, and the cleaned LAZ can be opened in GIS or LiDAR tools.</p>",
    output_file_list_to_html(generated_files, config$output_dir),
    "</section>\n",
    "</main>\n",
    "</body>\n",
    "</html>\n"
  )

  writeLines(html, report_path, useBytes = TRUE)
  report_path
}

write_outputs <- function(las, raw_check, clean_check, preprocess_log, summaries, config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  generated_files <- character()

  generated_files <- c(
    generated_files,
    write_table(
      preprocess_log,
      file.path(config$output_dir, "preprocessing_log.csv")
    )
  )

  generated_files <- c(
    generated_files,
    write_table(
      summaries$overall_summary,
      file.path(config$output_dir, "overall_summary.csv")
    )
  )

  generated_files <- c(
    generated_files,
    write_table(
      summaries$numeric_attribute_summary,
      file.path(config$output_dir, "numeric_attribute_summary.csv")
    )
  )

  if (!is.null(summaries$classification_summary)) {
    generated_files <- c(
      generated_files,
      write_table(
        summaries$classification_summary,
        file.path(config$output_dir, "classification_summary.csv")
      )
    )
  }

  if (!is.null(summaries$return_summary)) {
    generated_files <- c(
      generated_files,
      write_table(
        summaries$return_summary,
        file.path(config$output_dir, "return_summary.csv")
      )
    )
  }

  raw_check_path <- file.path(config$output_dir, "las_check_raw.txt")
  clean_check_path <- file.path(config$output_dir, "las_check_cleaned.txt")

  writeLines(raw_check, raw_check_path)
  writeLines(clean_check, clean_check_path)

  generated_files <- c(generated_files, raw_check_path, clean_check_path)

  plot_manifest <- write_report_plots(
    las,
    summaries,
    preprocess_log,
    config$output_dir
  )

  if (nrow(plot_manifest) > 0) {
    generated_files <- c(generated_files, plot_manifest$path)
  }

  if (isTRUE(config$write_cleaned_laz)) {
    input_ext <- tools::file_ext(config$input_file)
    output_ext <- if (nzchar(input_ext)) paste0(".", input_ext) else ".laz"
    output_name <- paste0(
      tools::file_path_sans_ext(basename(config$input_file)),
      "_cleaned",
      output_ext
    )
    clean_las_path <- unique_output_path(file.path(config$output_dir, output_name))

    lidR::writeLAS(las, clean_las_path, index = TRUE)
    generated_files <- c(generated_files, clean_las_path)

    index_path <- paste0(tools::file_path_sans_ext(clean_las_path), ".lax")
    if (file.exists(index_path)) {
      generated_files <- c(generated_files, index_path)
    }
  }

  report_path <- file.path(config$output_dir, "report.html")
  generated_files_path <- file.path(config$output_dir, "generated_files.txt")

  report_path <- write_html_report(
    report_path,
    plot_manifest,
    preprocess_log,
    summaries,
    config,
    c(generated_files, report_path, generated_files_path)
  )
  generated_files <- c(generated_files, report_path, generated_files_path)

  writeLines(generated_files, generated_files_path)

  generated_files
}

run_lidar_analytics <- function(config) {
  if (!file.exists(config$input_file)) {
    stop("Input file does not exist: ", config$input_file, call. = FALSE)
  }

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Reading LAS/LAZ: ", config$input_file)
  if (nzchar(config$read_filter)) {
    message("Using read filter: ", config$read_filter)
  }

  raw_header <- rlas::read.lasheader(config$input_file)
  las <- lidR::readLAS(config$input_file, filter = config$read_filter)

  if (is.null(las) || point_count(las) == 0) {
    stop("No points were loaded from the input file.", call. = FALSE)
  }

  raw_point_count <- point_count(las)

  message("Running raw LAS QA check.")
  raw_check <- capture_las_check(las)

  message("Cleaning/preprocessing point cloud.")
  preprocess_result <- preprocess_las(las, config)
  clean_las <- preprocess_result$las

  message("Running cleaned LAS QA check.")
  clean_check <- capture_las_check(clean_las)

  message("Building summary tables.")
  summaries <- build_summaries(clean_las, raw_header, raw_point_count, config)

  message("Writing outputs.")
  generated_files <- write_outputs(
    clean_las,
    raw_check,
    clean_check,
    preprocess_result$preprocess_log,
    summaries,
    config
  )

  message("Done. Generated files:")
  for (path in generated_files) {
    message("  - ", path)
  }

  invisible(generated_files)
}

run_default_lidar_analytics <- function(...) {
  config <- parse_args(character())

  overrides <- list(...)
  for (name in names(overrides)) {
    if (!name %in% names(config)) {
      stop("Unknown config option: ", name, call. = FALSE)
    }

    config[[name]] <- overrides[[name]]
  }

  run_lidar_analytics(config)
}

if (sys.nframe() == 0) {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  run_lidar_analytics(config)
}
