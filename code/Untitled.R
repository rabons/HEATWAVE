# Read in .RData and export to GeoTIFF
## Read all .RData files in subfolder "RData"
# Use terra library for doing so
library(terra)
library(tools)
# Get the path to this file and set it as WD
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Set work directory to current file location

rdata_files <- list.files(path = "RData", pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)
for (file in rdata_files) {
  load(file) # Load .RData file
  rdata_name <- file_path_sans_ext(basename(file)) # Get base name without extension
  
  if (inherits(raster_data, "Raster")) {
    output_filename <- file.path("GeoTIFF", paste0(rdata_name, ".tif")) # Define output filename
    dir.create("GeoTIFF", showWarnings = FALSE) # Create output directory if it doesn't exist
    writeRaster(raster_data, filename = output_filename, format = "GTiff", overwrite = TRUE) # Export to GeoTIFF
    cat("Exported:", output_filename, "\n")
    }
}
cat("All done!\n")