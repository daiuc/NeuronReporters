## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
# if (requireNamespace("workflowr", quietly = TRUE)) {
#   message("Loading .Rprofile for the current workflowr project")
#   library("workflowr")
# } else {
#   message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
# }

myopen <- function(file) {
  if (!tools::file_ext(file) %in% c("Rmd", "qmd")) {
    stop(cat("\nWrong file type. Only open/create .Rmd or .qmd files\n"))
  }

  current_path = dirname(".")
  doc_path = dirname(file)

  if (!file.exists(file)) {
    if (!dir.exists(doc_path)) dir.create(doc_path)
    file.create(file)
  }

  rstudioapi::navigateToFile(file)

  if (current_path != doc_path) {
    setwd(doc_path)
    cat(paste0("Set working dir to: \n", doc_path, "\n"))
  }
}
