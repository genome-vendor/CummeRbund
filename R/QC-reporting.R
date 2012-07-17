

############
#Utility functions
############
#include_graph <- function(width = 1, filename) {
#	paste("\includegraphics[width=", width, "\linewidth]{",
#	filename, "}", sep = "")
#}
#include_tbl <- function(width = 1, filename) {
#	print(xtable(filename), table.placement = "",
#	latex.environments = "", include.rownames = FALSE,
#	floating = FALSE)
#}
#subfloat_graph <- function(width, filename, caption = "") {
#	paste("\subfloat[", caption, "]{", "\begin{minipage}[h]{",
#	width, "\linewidth}\centering", include_graph(width = 1,
#	filename), "\end{minipage}}", sep = "")
#}
#subfloat_tbl <- function(width, filename, caption) {
#	paste("\subfloat[", caption, "]{", "\begin{minipage}[h]{",
#	width, "\linewidth}\centering", print(xtable(filename),
#	file = stderr(), table.placement = "",
#	latex.environments = "", include.rownames = FALSE,
#	floating = FALSE), "\end{minipage}}",
#	sep = "")
#}