
# ============================================================
# Helper function to read and label each sample
# ============================================================

read_sample <- function(path, sample_id, condition) {
  
  message("Reading sample: ", sample_id, " | ", condition)
  
  obj <- readRDS(path)
  
  # Add a sample-specific prefix to avoid duplicated cell barcodes.
  # This is safer than replacing characters with gsub().
  obj <- RenameCells(
    object = obj,
    add.cell.id = sample_id
  )
  
  obj$sample_id <- sample_id
  obj$condition <- condition
  obj$orig.ident <- sample_id
  
  return(obj)
}



read_sample_remove_doublets <- function(sample_id, condition, raw_file, tag) {
  
  message("Reading ", sample_id, " / ", condition)
  
  obj <- readRDS(
    file.path(datadir, "raw", condition, raw_file)
  )
  
  doublet_file <- file.path(
    datadir,
    "barcodes.doublets",
    condition,
    paste0("barcodes.", sample_id, ".rds")
  )
  
  if (!file.exists(doublet_file)) {
    stop("Missing doublet barcode file: ", doublet_file)
  }
  
  doublets <- readRDS(doublet_file)
  
  obj <- subset(
    x = obj,
    cells = doublets,
    invert = TRUE
  )
  
  obj <- rename_like_old(
    obj = obj,
    sample_id = sample_id
  )
  
  obj$orig.ident <- condition
  obj$tag <- tag
  
  return(obj)
}



rename_like_old <- function(obj, sample_id) {
  
  if (sample_id == "m1") {
    return(obj)
  }
  
  suffix <- sub("^m", "", sample_id)
  new_names <- gsub("1", suffix, colnames(obj))
  
  obj <- RenameCells(
    object = obj,
    new.names = new_names
  )
  
  return(obj)
}




collapse_to_single_leading_underscore <- function(obj) {
  
  new_names <- colnames(obj)
  new_names <- sub("^_+", "", new_names)
  new_names <- paste0("_", new_names)
  
  if (anyDuplicated(new_names)) {
    stop("Duplicated cell names after canonical renaming.")
  }
  
  if (!identical(colnames(obj), new_names)) {
    obj <- RenameCells(
      object = obj,
      new.names = new_names
    )
  }
  
  return(obj)
}


scale <- function(data){

all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

data <- RunPCA(data, features = VariableFeatures(object = data))

return(data)

}




PCA <- function(data){
## Selecting PCAs

# First metric 
#The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.

# Determine percent of variation associated with each PC
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Second metric
# The point where the percent change in variation between the consequtive PCs is less than 0.1%.

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.

co2
list <- list(co1,co2)
return(list)

}