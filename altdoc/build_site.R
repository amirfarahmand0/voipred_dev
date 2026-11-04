# Builds the package website (https://amirfarahmand0.github.io/voipred/) with
# altdoc and Quarto, from the package root:  Rscript altdoc/build_site.R
# The site is built in a temporary copy of the package, in which the tutorial
# (vignettes/articles/voipred.qmd) is the only vignette (the collaborator guides
# in vignettes/ are not on the website) and altdoc/home.md is the home page.
# The result is written to docs/, to be copied to the website repository.
pkg <- normalizePath(".")
tmp <- file.path(tempdir(), "voipred_site")
unlink(tmp, recursive = TRUE)
dir.create(tmp)
keep <- setdiff(list.files(pkg, all.files = TRUE, no.. = TRUE),
                c(".git", "docs", "reports", "_quarto"))
file.copy(file.path(pkg, keep), tmp, recursive = TRUE)
vig <- file.path(tmp, "vignettes")
unlink(list.files(vig, pattern = "\\.(Rmd|qmd)$", full.names = TRUE))
file.copy(file.path(pkg, "vignettes", "articles", c("voipred.qmd", "references.bib")), vig,
          overwrite = TRUE)
unlink(file.path(vig, "articles"), recursive = TRUE)
unlink(file.path(tmp, "NEWS.md"))
file.copy(file.path(pkg, "altdoc", "home.md"), file.path(tmp, "README.md"), overwrite = TRUE)
altdoc::render_docs(path = tmp, verbose = FALSE)
unlink(file.path(pkg, "docs"), recursive = TRUE)
file.copy(file.path(tmp, "docs"), pkg, recursive = TRUE)
# redirects from the addresses of the former (pkgdown) website
redirect <- function(from, to) {
  dir.create(dirname(file.path(pkg, "docs", from)), recursive = TRUE, showWarnings = FALSE)
  writeLines(sprintf(paste0('<!DOCTYPE html><html><head><meta charset="utf-8">',
                            '<meta http-equiv="refresh" content="0; url=%s">',
                            '<link rel="canonical" href="%s"></head><body>',
                            '<a href="%s">%s</a></body></html>'), to, to, to, to),
             file.path(pkg, "docs", from))
}
redirect("articles/voipred.html", "../vignettes/voipred.html")
redirect("articles/index.html", "../vignettes/voipred.html")
redirect("reference/index.html", "../man/nb_val.html")
message("Website built in docs/")
