# R ALS Point Cloud Analytics

Basic R workflow for cleaning/preprocessing an airborne LiDAR point cloud and creating simple analytical summaries, plots, and an HTML report.

The prototype uses the included sample file:

```text
USGS_LPC_NV_USFSR4_D23_1301_2064_EPSG-6341_LAS2023.laz
```

## What To Open First

After a run finishes, open:

```text
outputs/.../report.html
```

The HTML report brings the key metrics, plots, summary tables, and file links into one place so you do not have to inspect every CSV manually.

## Install R Packages

The script requires:

```r
install.packages(c("lidR", "data.table", "rlas"))
```

## Quick Test In RStudio

Use this first to confirm the workflow is working. It runs a small crop of the sample LAZ and skips writing the cleaned LAZ file, so it is much faster than a full run.

```r
setwd("C:/Users/theoj/Files/Github/R_ALS-point-cloud-analytics")
source("basic_lidar_analytics.R")
run_quick_test()
```

Outputs are written to:

```text
outputs/quick_test/
```

Open:

```text
outputs/quick_test/report.html
```

## Full Run In RStudio

This processes the full sample LAZ using the default settings.

```r
setwd("C:/Users/theoj/Files/Github/R_ALS-point-cloud-analytics")
source("basic_lidar_analytics.R")
run_default_lidar_analytics()
```

Outputs are written to:

```text
outputs/lidar_basic_analytics/
```

## Run From PowerShell

If `Rscript` is on your PATH:

```powershell
cd C:\Users\theoj\Files\Github\R_ALS-point-cloud-analytics
Rscript basic_lidar_analytics.R --quick-test
```

On this machine, R is installed here:

```powershell
cd C:\Users\theoj\Files\Github\R_ALS-point-cloud-analytics
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' basic_lidar_analytics.R --quick-test
```

For the full run:

```powershell
& 'C:\Program Files\R\R-4.5.2\bin\Rscript.exe' basic_lidar_analytics.R
```

## Main Outputs

- `report.html`: automatic summary report with plots and tables.
- `plots/`: PNG charts used by the report.
- `overall_summary.csv`: run settings, CRS, point counts, extent, elevation range, and density.
- `preprocessing_log.csv`: cleaning steps and how many points each step removed.
- `classification_summary.csv`: point counts and elevation summaries by LAS classification.
- `numeric_attribute_summary.csv`: min, max, mean, standard deviation, and percentiles for numeric attributes.
- `return_summary.csv`: counts by return number and number of returns.
- `las_check_raw.txt`: `lidR::las_check()` report before cleaning.
- `las_check_cleaned.txt`: `lidR::las_check()` report after cleaning.
- `*_cleaned.laz`: cleaned point cloud, if writing cleaned LAZ is enabled.
- `*_cleaned.lax`: spatial index for the cleaned LAZ.
- `generated_files.txt`: list of files created by the run.

## Useful Options

Run a small test crop:

```powershell
Rscript basic_lidar_analytics.R --quick-test
```

Skip writing the cleaned LAZ:

```powershell
Rscript basic_lidar_analytics.R --no-write-cleaned
```

Run a custom spatial crop:

```powershell
Rscript basic_lidar_analytics.R input.laz outputs/custom_crop --read-filter="-keep_xy xmin ymin xmax ymax"
```

Add height above ground if ground points are available:

```powershell
Rscript basic_lidar_analytics.R input.laz outputs/with_hag --add-hag
```

## Notes

The full sample LAZ has about 26.7 million points, so a full run can take time and produce a large cleaned LAZ. Use `run_quick_test()` first, then run the full analysis once the report looks good.
