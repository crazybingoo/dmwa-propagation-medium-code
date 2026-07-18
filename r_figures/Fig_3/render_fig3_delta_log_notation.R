#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}
source_script <- file.path(script_dir, "make_fig3_eta_decomposition_font10p5_label12p5_no_colorbar_b_legend_inside.R")
source_lines <- readLines(source_script, warn = FALSE, encoding = "UTF-8")

replacements <- c(
  'out_base <- file.path(base_dir, "Fig3_eta_decomposition_font10p5_label12p5_no_colorbar_b_legend_inside")' =
    'out_base <- file.path(base_dir, "Fig3_eta_decomposition_delta_log_notation")',
  'x = 5.64, y = 14.25, label = "+dlog R",' =
    'x = 5.64, y = 14.25, label = "\\u0394log R",',
  'x = 5.64, y = -13.55, label = paste0("-dlog ", lambda1_sym),' =
    'x = 5.64, y = -13.55, label = paste0("\\u2212\\u0394log ", lambda1_sym),',
  'y = paste0("Contribution to relative change,\\n100 ", times_sym, " dlog ", eta_italic_sym)' =
    'y = paste0("Contribution to relative change,\\n100 ", times_sym, " \\u0394log ", eta_italic_sym)',
  'y = paste0("100 ", times_sym, " dlog(", R_italic_sym, "/", lambda_italic_sym, subscript_1, ")")' =
    'y = paste0("100 ", times_sym, " \\u0394log(", R_italic_sym, "/", lambda_italic_sym, subscript_1, ")")',
  'y = "Change from pre, 100 x dlog"' =
    'y = "Change from pre, 100 \\u00d7 \\u0394log"'
)

for (old in names(replacements)) {
  matched <- trimws(source_lines) == old
  if (sum(matched) != 1L) {
    stop("Expected one exact match for: ", old, "; found ", sum(matched), ".")
  }
  source_lines[matched] <- unname(replacements[[old]])
}

export_line <- grep("^save_pub_r\\(fig, out_base\\)$", source_lines)
if (length(export_line) != 1L) {
  stop("Could not identify the unique Fig. 3 export line.")
}

eval(parse(text = source_lines[seq_len(export_line)]), envir = globalenv())
