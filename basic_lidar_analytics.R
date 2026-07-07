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
options(lidR.progress = FALSE)

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
      "  Rscript basic_lidar_analytics.R --quick-test",
      "",
      "Options:",
      "  --help                    Show this help message.",
      "  --quick-test              Run a small sample crop and write outputs/quick_test.",
      "  --no-write-cleaned        Skip writing the cleaned/preprocessed LAS/LAZ.",
      "  --add-hag                 Add height above ground (HAG) if ground points exist.",
      "  --drop-overlap            Drop LAS 1.4 overlap-flagged points when present.",
      "  --read-filter=FILTER      Pass a LASlib filter string to lidR::readLAS().",
      "  --max-3d-points=N         Maximum sampled points in point_cloud_3d.html.",
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
    z_trim = NULL,
    quick_test = FALSE,
    max_3d_points = 75000L
  )

  positional <- character()

  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      usage()
      quit(save = "no", status = 0)
    } else if (arg == "--quick-test") {
      config$quick_test <- TRUE
    } else if (arg == "--no-write-cleaned") {
      config$write_cleaned_laz <- FALSE
    } else if (arg == "--add-hag") {
      config$add_hag <- TRUE
    } else if (arg == "--drop-overlap") {
      config$drop_overlap <- TRUE
    } else if (startsWith(arg, "--read-filter=")) {
      config$read_filter <- sub("^--read-filter=", "", arg)
    } else if (startsWith(arg, "--max-3d-points=")) {
      max_3d_points <- as.integer(sub("^--max-3d-points=", "", arg))

      if (is.na(max_3d_points) || max_3d_points < 1000) {
        stop("--max-3d-points must be an integer of at least 1000.", call. = FALSE)
      }

      config$max_3d_points <- max_3d_points
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

  if (isTRUE(config$quick_test)) {
    config$output_dir <- file.path(script_dir, "outputs", "quick_test")
    config$read_filter <- quick_test_filter
    config$write_cleaned_laz <- FALSE
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

format_duration <- function(start_time) {
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if (elapsed < 60) {
    return(paste0(round(elapsed, 1), " seconds"))
  }

  paste0(round(elapsed / 60, 1), " minutes")
}

format_bool <- function(x) {
  if (isTRUE(x)) {
    return("yes")
  }

  "no"
}

message_step <- function(text) {
  message("")
  message(text)
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

js_escape_string <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub('"', '\\"', x, fixed = TRUE)
  x <- gsub("\r", "", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  paste0('"', x, '"')
}

js_numeric_array <- function(x, digits = 5L, values_per_line = 16L) {
  values <- ifelse(
    is.na(x),
    "NaN",
    format(round(as.numeric(x), digits), scientific = FALSE, trim = TRUE)
  )

  groups <- split(values, ceiling(seq_along(values) / values_per_line))
  lines <- vapply(groups, paste, character(1), collapse = ",")
  paste0("[\n", paste(lines, collapse = ",\n"), "\n]")
}

js_integer_array <- function(x, values_per_line = 24L) {
  values <- ifelse(is.na(x), "0", as.character(as.integer(x)))
  groups <- split(values, ceiling(seq_along(values) / values_per_line))
  lines <- vapply(groups, paste, character(1), collapse = ",")
  paste0("[\n", paste(lines, collapse = ",\n"), "\n]")
}

classification_labels_js <- function(classes) {
  if (is.null(classes) || length(classes) == 0) {
    return("{}")
  }

  classes <- sort(unique(as.integer(classes)))
  entries <- vapply(
    classes,
    function(class_id) {
      label <- class_lookup[[as.character(class_id)]]
      if (is.null(label) || is.na(label)) {
        label <- "Unknown or user-defined"
      }

      paste0('"', class_id, '":', js_escape_string(paste0(class_id, ": ", label)))
    },
    character(1)
  )

  paste0("{", paste(entries, collapse = ","), "}")
}

write_interactive_3d_viewer <- function(las, config) {
  dt <- las@data

  if (!all(c("X", "Y", "Z") %in% names(dt)) || nrow(dt) == 0) {
    warning("Could not create 3D viewer because X, Y, or Z data is missing.", call. = FALSE)
    return(NULL)
  }

  viewer_path <- file.path(config$output_dir, "point_cloud_3d.html")
  sample_cap <- min(as.integer(config$max_3d_points), nrow(dt))
  idx <- sample_indices(nrow(dt), sample_cap, seed = 45L)
  sampled <- dt[idx]

  xmin <- safe_min(dt$X)
  xmax <- safe_max(dt$X)
  ymin <- safe_min(dt$Y)
  ymax <- safe_max(dt$Y)
  zmin <- safe_min(dt$Z)
  zmax <- safe_max(dt$Z)

  center_x <- mean(c(xmin, xmax))
  center_y <- mean(c(ymin, ymax))
  center_z <- mean(c(zmin, zmax))
  scale <- max(c(xmax - xmin, ymax - ymin, zmax - zmin), na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }

  positions <- as.vector(rbind(
    (sampled$X - center_x) / scale,
    (sampled$Z - center_z) / scale,
    -(sampled$Y - center_y) / scale
  ))

  has_intensity <- "Intensity" %in% names(sampled)
  has_classification <- "Classification" %in% names(sampled)
  has_hag <- "hag" %in% names(sampled)

  intensity_values <- if (has_intensity) js_numeric_array(sampled$Intensity, digits = 3L) else "null"
  classification_values <- if (has_classification) js_integer_array(sampled$Classification) else "null"
  hag_values <- if (has_hag) js_numeric_array(sampled$hag, digits = 3L) else "null"
  class_labels <- if (has_classification) classification_labels_js(sampled$Classification) else "{}"

  metadata <- data.table::data.table(
    field = c(
      "Cleaned points",
      "Sampled viewer points",
      "Sample cap",
      "Color modes",
      "X range",
      "Y range",
      "Z range"
    ),
    value = c(
      format_count(nrow(dt)),
      format_count(length(idx)),
      format_count(config$max_3d_points),
      paste(
        c(
          "Elevation",
          if (has_classification) "Classification" else NULL,
          if (has_intensity) "Intensity" else NULL,
          if (has_hag) "Height above ground" else NULL
        ),
        collapse = ", "
      ),
      paste(format_measure(xmin, 2), "to", format_measure(xmax, 2)),
      paste(format_measure(ymin, 2), "to", format_measure(ymax, 2)),
      paste(format_measure(zmin, 2), "to", format_measure(zmax, 2))
    )
  )

  metadata_items <- paste0(
    "<dl>",
    paste0(
      "<dt>", html_escape(metadata$field), "</dt><dd>", html_escape(metadata$value), "</dd>",
      collapse = ""
    ),
    "</dl>"
  )

  html <- paste0(
    "<!doctype html>\n",
    "<html lang=\"en\">\n",
    "<head>\n",
    "<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>Interactive 3D Point Cloud Viewer</title>\n",
    "<style>\n",
    "html,body{margin:0;width:100%;height:100%;overflow:hidden;background:#0f172a;color:#e5e7eb;font-family:Arial,sans-serif;}\n",
    "#viewer{display:block;width:100vw;height:100vh;background:#0f172a;}\n",
    ".panel{position:absolute;left:16px;top:16px;max-width:360px;background:rgba(15,23,42,.88);border:1px solid rgba(148,163,184,.35);border-radius:8px;padding:14px;box-shadow:0 12px 30px rgba(0,0,0,.28);}\n",
    "h1{font-size:18px;margin:0 0 8px;} p{margin:8px 0;color:#cbd5e1;font-size:13px;line-height:1.4;} button{border:1px solid rgba(148,163,184,.5);background:#1e293b;color:#f8fafc;border-radius:6px;padding:7px 10px;margin:3px;cursor:pointer;} button.active{background:#2563eb;border-color:#60a5fa;} button:disabled{opacity:.45;cursor:not-allowed;} label{display:block;margin:10px 0 4px;font-size:13px;color:#cbd5e1;} input[type=range]{width:100%;} dl{display:grid;grid-template-columns:auto 1fr;gap:4px 12px;margin:10px 0 0;font-size:12px;} dt{color:#94a3b8;} dd{margin:0;color:#e5e7eb;} .hint{position:absolute;right:16px;bottom:14px;color:#cbd5e1;background:rgba(15,23,42,.78);border-radius:8px;padding:8px 10px;font-size:12px;} .legend{font-size:12px;color:#cbd5e1;margin-top:8px;max-height:130px;overflow:auto;} .legend-row{display:flex;align-items:center;gap:7px;margin:3px 0;} .swatch{width:12px;height:12px;border-radius:50%;display:inline-block;}\n",
    "@media(max-width:720px){.panel{left:8px;right:8px;top:8px;max-width:none}.hint{display:none}}\n",
    "</style>\n",
    "</head>\n",
    "<body>\n",
    "<canvas id=\"viewer\"></canvas>\n",
    "<section class=\"panel\">\n",
    "<h1>Interactive 3D Point Cloud</h1>\n",
    "<p>Drag to rotate. Use the mouse wheel or trackpad to zoom. Double-click to reset the view.</p>\n",
    "<div id=\"modeButtons\"></div>\n",
    "<label for=\"pointSize\">Point size</label><input id=\"pointSize\" type=\"range\" min=\"1\" max=\"6\" step=\"0.5\" value=\"2\">\n",
    metadata_items,
    "<div id=\"legend\" class=\"legend\"></div>\n",
    "</section>\n",
    "<div class=\"hint\">X = easting, vertical = elevation, depth = northing</div>\n",
    "<script>\n",
    "'use strict';\n",
    "const pointData={\n",
    "positions:new Float32Array(", js_numeric_array(positions, digits = 6L), "),\n",
    "z:new Float32Array(", js_numeric_array(sampled$Z, digits = 3L), "),\n",
    "intensity:", if (has_intensity) paste0("new Float32Array(", intensity_values, ")") else "null", ",\n",
    "classification:", if (has_classification) paste0("new Int16Array(", classification_values, ")") else "null", ",\n",
    "hag:", if (has_hag) paste0("new Float32Array(", hag_values, ")") else "null", ",\n",
    "classLabels:", class_labels, ",\n",
    "pointCount:", length(idx), "\n",
    "};\n",
    "const canvas=document.getElementById('viewer');\n",
    "const gl=canvas.getContext('webgl',{antialias:true,preserveDrawingBuffer:false});\n",
    "if(!gl){document.body.innerHTML='<p style=\"padding:20px\">WebGL is not available in this browser.</p>';throw new Error('WebGL unavailable');}\n",
    "const vs='attribute vec3 aPosition;attribute vec3 aColor;uniform mat4 uMatrix;uniform float uPointSize;varying vec3 vColor;void main(){gl_Position=uMatrix*vec4(aPosition,1.0);gl_PointSize=uPointSize;vColor=aColor;}';\n",
    "const fs='precision mediump float;varying vec3 vColor;void main(){vec2 c=gl_PointCoord-vec2(0.5);if(dot(c,c)>0.25)discard;gl_FragColor=vec4(vColor,1.0);}';\n",
    "function shader(type,src){const s=gl.createShader(type);gl.shaderSource(s,src);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw new Error(gl.getShaderInfoLog(s));return s;}\n",
    "const program=gl.createProgram();gl.attachShader(program,shader(gl.VERTEX_SHADER,vs));gl.attachShader(program,shader(gl.FRAGMENT_SHADER,fs));gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw new Error(gl.getProgramInfoLog(program));gl.useProgram(program);\n",
    "const posBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,posBuffer);gl.bufferData(gl.ARRAY_BUFFER,pointData.positions,gl.STATIC_DRAW);const posLoc=gl.getAttribLocation(program,'aPosition');gl.enableVertexAttribArray(posLoc);gl.vertexAttribPointer(posLoc,3,gl.FLOAT,false,0,0);\n",
    "const colorBuffer=gl.createBuffer();const colorLoc=gl.getAttribLocation(program,'aColor');gl.enableVertexAttribArray(colorLoc);gl.vertexAttribPointer(colorLoc,3,gl.FLOAT,false,0,0);\n",
    "const matrixLoc=gl.getUniformLocation(program,'uMatrix');const sizeLoc=gl.getUniformLocation(program,'uPointSize');\n",
    "const axisPositions=new Float32Array([-0.55,-0.55,0,0.55,-0.55,0, -0.55,-0.55,0,-0.55,0.55,0, -0.55,-0.55,0,-0.55,-0.55,-0.55]);\n",
    "const axisColors=new Float32Array([0.95,0.25,0.25,0.95,0.25,0.25, 0.25,0.95,0.4,0.25,0.95,0.4, 0.35,0.6,1,0.35,0.6,1]);\n",
    "const axisPosBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,axisPosBuffer);gl.bufferData(gl.ARRAY_BUFFER,axisPositions,gl.STATIC_DRAW);const axisColorBuffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,axisColorBuffer);gl.bufferData(gl.ARRAY_BUFFER,axisColors,gl.STATIC_DRAW);\n",
    "const classPalette={0:[0.60,0.60,0.60],1:[0.62,0.65,0.70],2:[0.55,0.36,0.17],3:[0.49,0.70,0.26],4:[0.26,0.63,0.28],5:[0.10,0.37,0.13],6:[0.62,0.62,0.62],7:[0.90,0.24,0.23],8:[0.93,0.69,0.13],9:[0.12,0.47,0.71],10:[0.66,0.34,0.16],11:[0.36,0.36,0.36],12:[0.74,0.74,0.74],13:[0.58,0.40,0.74],14:[0.58,0.40,0.74],15:[0.58,0.40,0.74],16:[0.58,0.40,0.74],17:[0.50,0.50,0.50],18:[0.90,0.12,0.12]};\n",
    "const viridis=[[0,[0.267,0.005,0.329]],[0.25,[0.230,0.322,0.545]],[0.5,[0.128,0.567,0.551]],[0.75,[0.369,0.789,0.382]],[1,[0.993,0.906,0.144]]];\n",
    "function range(values){let lo=Infinity,hi=-Infinity;for(let i=0;i<values.length;i++){const v=values[i];if(Number.isFinite(v)){if(v<lo)lo=v;if(v>hi)hi=v;}}return [lo,hi];}\n",
    "function ramp(t){t=Math.max(0,Math.min(1,t));for(let i=1;i<viridis.length;i++){if(t<=viridis[i][0]){const a=viridis[i-1],b=viridis[i],f=(t-a[0])/(b[0]-a[0]);return [a[1][0]+(b[1][0]-a[1][0])*f,a[1][1]+(b[1][1]-a[1][1])*f,a[1][2]+(b[1][2]-a[1][2])*f];}}return viridis[viridis.length-1][1];}\n",
    "function fillColors(mode){const colors=new Float32Array(pointData.pointCount*3);let values=pointData.z;if(mode==='intensity')values=pointData.intensity;if(mode==='hag')values=pointData.hag;const r=values?range(values):[0,1];for(let i=0;i<pointData.pointCount;i++){let c;if(mode==='classification'&&pointData.classification){c=classPalette[pointData.classification[i]]||[0.86,0.86,0.86];}else{const denom=(r[1]-r[0])||1;c=ramp(((values?values[i]:pointData.z[i])-r[0])/denom);}colors[i*3]=c[0];colors[i*3+1]=c[1];colors[i*3+2]=c[2];}gl.bindBuffer(gl.ARRAY_BUFFER,colorBuffer);gl.bufferData(gl.ARRAY_BUFFER,colors,gl.DYNAMIC_DRAW);updateLegend(mode,r);render();}\n",
    "function updateLegend(mode,r){const el=document.getElementById('legend');if(mode==='classification'&&pointData.classification){const used=[...new Set(Array.from(pointData.classification))].sort((a,b)=>a-b);el.innerHTML=used.map(k=>{const c=classPalette[k]||[.86,.86,.86];const rgb=c.map(v=>Math.round(v*255)).join(',');return '<div class=\"legend-row\"><span class=\"swatch\" style=\"background:rgb('+rgb+')\"></span>'+ (pointData.classLabels[k]||k) +'</div>';}).join('');}else{const label=mode==='hag'?'Height above ground':mode.charAt(0).toUpperCase()+mode.slice(1);el.innerHTML='<div>'+label+' color range: '+r[0].toFixed(2)+' to '+r[1].toFixed(2)+'</div>';}}\n",
    "function matMul(a,b){const o=new Float32Array(16);for(let r=0;r<4;r++){for(let c=0;c<4;c++){o[c*4+r]=a[0*4+r]*b[c*4+0]+a[1*4+r]*b[c*4+1]+a[2*4+r]*b[c*4+2]+a[3*4+r]*b[c*4+3];}}return o;}\n",
    "function perspective(fov,aspect,near,far){const f=1/Math.tan(fov/2),nf=1/(near-far);return new Float32Array([f/aspect,0,0,0,0,f,0,0,0,0,(far+near)*nf,-1,0,0,2*far*near*nf,0]);}\n",
    "function translate(z){return new Float32Array([1,0,0,0,0,1,0,0,0,0,1,0,0,0,z,1]);}\n",
    "function rotX(a){const c=Math.cos(a),s=Math.sin(a);return new Float32Array([1,0,0,0,0,c,s,0,0,-s,c,0,0,0,0,1]);}\n",
    "function rotY(a){const c=Math.cos(a),s=Math.sin(a);return new Float32Array([c,0,-s,0,0,1,0,0,s,0,c,0,0,0,0,1]);}\n",
    "let yaw=0.65,pitch=-0.55,distance=2.4,pointSize=2,drag=false,lastX=0,lastY=0,currentMode='elevation';\n",
    "function matrix(){const aspect=canvas.width/canvas.height;return matMul(perspective(Math.PI/4,aspect,0.01,20),matMul(translate(-distance),matMul(rotX(pitch),rotY(yaw))));}\n",
    "function resize(){const dpr=Math.min(window.devicePixelRatio||1,2);const w=Math.floor(canvas.clientWidth*dpr),h=Math.floor(canvas.clientHeight*dpr);if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;gl.viewport(0,0,w,h);}render();}\n",
    "function render(){gl.enable(gl.DEPTH_TEST);gl.clearColor(0.06,0.09,0.16,1);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);gl.uniformMatrix4fv(matrixLoc,false,matrix());gl.uniform1f(sizeLoc,pointSize*(window.devicePixelRatio||1));gl.bindBuffer(gl.ARRAY_BUFFER,posBuffer);gl.vertexAttribPointer(posLoc,3,gl.FLOAT,false,0,0);gl.bindBuffer(gl.ARRAY_BUFFER,colorBuffer);gl.vertexAttribPointer(colorLoc,3,gl.FLOAT,false,0,0);gl.drawArrays(gl.POINTS,0,pointData.pointCount);gl.uniform1f(sizeLoc,1);gl.bindBuffer(gl.ARRAY_BUFFER,axisPosBuffer);gl.vertexAttribPointer(posLoc,3,gl.FLOAT,false,0,0);gl.bindBuffer(gl.ARRAY_BUFFER,axisColorBuffer);gl.vertexAttribPointer(colorLoc,3,gl.FLOAT,false,0,0);gl.drawArrays(gl.LINES,0,6);}\n",
    "canvas.addEventListener('pointerdown',e=>{drag=true;lastX=e.clientX;lastY=e.clientY;canvas.setPointerCapture(e.pointerId);});canvas.addEventListener('pointermove',e=>{if(!drag)return;const dx=e.clientX-lastX,dy=e.clientY-lastY;lastX=e.clientX;lastY=e.clientY;yaw+=dx*0.006;pitch=Math.max(-1.45,Math.min(1.45,pitch+dy*0.006));render();});canvas.addEventListener('pointerup',()=>{drag=false;});canvas.addEventListener('wheel',e=>{e.preventDefault();distance=Math.max(0.75,Math.min(8,distance*Math.exp(e.deltaY*0.001)));render();},{passive:false});canvas.addEventListener('dblclick',()=>{yaw=0.65;pitch=-0.55;distance=2.4;render();});\n",
    "document.getElementById('pointSize').addEventListener('input',e=>{pointSize=parseFloat(e.target.value);render();});\n",
    "const modes=[['elevation','Elevation'], pointData.classification?['classification','Classification']:null, pointData.intensity?['intensity','Intensity']:null, pointData.hag?['hag','Height above ground']:null].filter(Boolean);const modeBox=document.getElementById('modeButtons');modes.forEach(([mode,label])=>{const b=document.createElement('button');b.textContent=label;b.addEventListener('click',()=>{currentMode=mode;document.querySelectorAll('#modeButtons button').forEach(x=>x.classList.toggle('active',x===b));fillColors(mode);});modeBox.appendChild(b);if(mode==='elevation')b.classList.add('active');});\n",
    "window.addEventListener('resize',resize);fillColors(currentMode);resize();\n",
    "</script>\n",
    "</body>\n",
    "</html>\n"
  )

  writeLines(html, viewer_path, useBytes = TRUE)
  viewer_path
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
      "density_points_per_map_unit2",
      "max_3d_viewer_points"
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
      density,
      config$max_3d_points
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

interactive_viewer_section_html <- function(generated_files, output_dir) {
  viewer_path <- generated_files[basename(generated_files) == "point_cloud_3d.html"]

  if (length(viewer_path) == 0) {
    return("")
  }

  href <- report_href(viewer_path[[1]], output_dir)

  paste0(
    "<section>",
    "<h2>Interactive 3D Viewer</h2>",
    "<p>Open the sampled WebGL point cloud viewer to rotate, zoom, and color the cleaned point cloud by available attributes.</p>",
    "<p><a class=\"button-link\" href=\"", html_escape(href), "\">Open interactive 3D viewer</a></p>",
    "<p class=\"subtle\">The viewer is sampled for browser performance. CSV summaries and report metrics still use all cleaned points.</p>",
    "</section>\n"
  )
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
    ".button-link{display:inline-block;background:#2563eb;color:#fff;text-decoration:none;border-radius:6px;padding:10px 14px;font-weight:700;}",
    ".button-link:hover{background:#1d4ed8;}",
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
    interactive_viewer_section_html(generated_files, config$output_dir),
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

  message("Writing summary CSV files.")
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

  message("Writing LAS QA check reports.")
  writeLines(raw_check, raw_check_path)
  writeLines(clean_check, clean_check_path)

  generated_files <- c(generated_files, raw_check_path, clean_check_path)

  message("Creating PNG plots.")
  plot_manifest <- write_report_plots(
    las,
    summaries,
    preprocess_log,
    config$output_dir
  )

  if (nrow(plot_manifest) > 0) {
    generated_files <- c(generated_files, plot_manifest$path)
  }

  message("Writing sampled interactive 3D viewer.")
  viewer_path <- write_interactive_3d_viewer(las, config)
  if (!is.null(viewer_path) && file.exists(viewer_path)) {
    generated_files <- c(generated_files, viewer_path)
  }

  if (isTRUE(config$write_cleaned_laz)) {
    message("Writing cleaned LAZ file and spatial index.")
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
  } else {
    message("Skipping cleaned LAZ write because write_cleaned_laz is FALSE.")
  }

  report_path <- file.path(config$output_dir, "report.html")
  generated_files_path <- file.path(config$output_dir, "generated_files.txt")

  message("Writing HTML report.")
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
  start_time <- Sys.time()

  if (!file.exists(config$input_file)) {
    stop("Input file does not exist: ", config$input_file, call. = FALSE)
  }

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  message_step("Starting LiDAR point cloud analytics.")
  message("Input file: ", config$input_file)
  message("Output folder: ", config$output_dir)
  message("Quick test mode: ", format_bool(config$quick_test))
  message("Write cleaned LAZ: ", format_bool(config$write_cleaned_laz))
  message("3D viewer sample cap: ", format_count(config$max_3d_points))
  if (nzchar(config$read_filter)) {
    message("Read filter: ", config$read_filter)
  } else {
    message("Read filter: none")
  }

  message_step("Reading LAS/LAZ header and points.")
  raw_header <- rlas::read.lasheader(config$input_file)
  header_count <- read_header_value(raw_header, "Number of point records")
  message("Header point count: ", format_count(as.numeric(header_count)))

  utils::capture.output(
    las <- lidR::readLAS(config$input_file, filter = config$read_filter),
    type = "output"
  )

  if (is.null(las) || point_count(las) == 0) {
    stop("No points were loaded from the input file.", call. = FALSE)
  }

  raw_point_count <- point_count(las)
  message("Loaded point count: ", format_count(raw_point_count))

  message_step("Running raw LAS QA check.")
  raw_check <- capture_las_check(las)

  message_step("Cleaning/preprocessing point cloud.")
  preprocess_result <- preprocess_las(las, config)
  clean_las <- preprocess_result$las
  clean_point_count <- point_count(clean_las)
  total_removed <- raw_point_count - clean_point_count
  message("Cleaning complete.")
  message("Retained points: ", format_count(clean_point_count))
  message("Removed points: ", format_count(total_removed))

  removed_steps <- preprocess_result$preprocess_log[
    points_removed > 0,
    list(step, points_removed)
  ]
  if (nrow(removed_steps) > 0) {
    for (i in seq_len(nrow(removed_steps))) {
      message(
        "  - ",
        gsub("_", " ", removed_steps$step[[i]], fixed = TRUE),
        ": ",
        format_count(removed_steps$points_removed[[i]])
      )
    }
  } else {
    message("  - No points were removed by preprocessing.")
  }

  message_step("Running cleaned LAS QA check.")
  clean_check <- capture_las_check(clean_las)

  message_step("Building summary tables.")
  summaries <- build_summaries(clean_las, raw_header, raw_point_count, config)
  message("Summary tables ready.")

  message_step("Writing outputs.")
  generated_files <- write_outputs(
    clean_las,
    raw_check,
    clean_check,
    preprocess_result$preprocess_log,
    summaries,
    config
  )

  report_path <- generated_files[basename(generated_files) == "report.html"]
  if (length(report_path) == 0) {
    report_path <- file.path(config$output_dir, "report.html")
  } else {
    report_path <- report_path[[1]]
  }

  generated_files_path <- generated_files[basename(generated_files) == "generated_files.txt"]
  if (length(generated_files_path) == 0) {
    generated_files_path <- file.path(config$output_dir, "generated_files.txt")
  } else {
    generated_files_path <- generated_files_path[[1]]
  }

  message_step("Done.")
  message("Runtime: ", format_duration(start_time))
  message("Generated file count: ", length(generated_files))
  message("Open this first: ", report_path)
  message("Full output list: ", generated_files_path)

  invisible(generated_files)
}

run_default_lidar_analytics <- function(...) {
  config <- parse_args(character())

  config <- apply_config_overrides(config, list(...))

  run_lidar_analytics(config)
}

run_quick_test <- function(...) {
  config <- parse_args("--quick-test")

  config <- apply_config_overrides(config, list(...))

  run_lidar_analytics(config)
}

apply_config_overrides <- function(config, overrides) {
  for (name in names(overrides)) {
    if (!name %in% names(config)) {
      stop("Unknown config option: ", name, call. = FALSE)
    }

    config[[name]] <- overrides[[name]]
  }

  if ("input_file" %in% names(overrides)) {
    config$input_file <- normalizePath(config$input_file, mustWork = FALSE)
  }

  if ("output_dir" %in% names(overrides)) {
    config$output_dir <- normalizePath(config$output_dir, mustWork = FALSE)
  }

  config
}

if (sys.nframe() == 0) {
  config <- parse_args(commandArgs(trailingOnly = TRUE))
  run_lidar_analytics(config)
}
