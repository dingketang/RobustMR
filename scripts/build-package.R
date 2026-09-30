# Run from the package root: Rscript scripts/build-package.R
needed <- c("devtools", "roxygen2", "testthat")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install these development packages first: ",
       paste(missing, collapse = ", "), call. = FALSE)
}
if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the RobustMR package root.", call. = FALSE)
}

devtools::document()
devtools::test()
result <- devtools::check(document = FALSE, manual = FALSE, error_on = "never")
if (length(result$errors) || length(result$warnings)) {
  stop("Resolve R CMD check errors and warnings before building a release.",
       call. = FALSE)
}
archive <- devtools::build(path = "..")
message("Built source archive: ", archive)
