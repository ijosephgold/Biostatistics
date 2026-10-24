export_plot <- function(filename, plot = last_plot(), width = 7, height = 5) {
  if (!dir.exists("output")) dir.create("output")
  filepath <- file.path("output", filename)
  
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(filepath, plot = plot, width = width, height = height, dpi = 300)
    message(paste("Saved ggplot to:", filepath))
  } else {
    dev.copy(jpeg, filename = filepath, width = width * 100, height = height * 100)
    dev.off()
    message(paste("Saved base R plot to:", filepath))
  }
}