# Run with source("project/create-logo.R") in Positron/the R console,
# or Rscript project/create-logo.R. Saves the image used by the GitHub README.

library(ggplot2)
library(dplyr)
library(sf)
library(terra)
library(hexSticker)
library(png)

# Design parameters -----------------------------------------------------------
word <- "spseg"
point_layout <- "hexagonal" # "random" or "hexagonal"
point_count <- 5000L # Exact for random; approximate for hexagonal
point_size <- 0.6
seed <- 357
noise <- 0.3
letter_strength <- 2
center_strength <- 3 # Weight on background log-distance, outside letters
radius <- 0.94
colour_quantile <- 0.9
viridis_option <- "D"
word_bounds <- c(left = -0.84, right = 0.84, bottom = -0.42, top = 0.58)
font_family <- "sans"
background <- "#101519"
border_colour <- "#303940"
border_size <- 0.65
output_dpi <- 300
output_file <- file.path("man", "figures", "logo.png")
show_logo <- TRUE # Display in interactive Positron/R sessions
shape1 <- 0.2
shape2 <- 0.8

# Find the output location -----------------------------------------------------
script <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(script)) {
  file.path(dirname(sub("^--file=", "", script[1L])), "..")
} else if (file.exists("DESCRIPTION")) {
  "."
} else {
  ".."
}
output <- file.path(normalizePath(root), output_file)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

# Only command-line runs need a temporary device for font measurements.
# Interactive runs retain the session's normal plot device for print(logo).
if (!interactive()) {
  drawing_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(drawing_file)
  drawing_device <- grDevices::dev.cur()
}

# Turn the lettering into polygon shapes --------------------------------------
mask_file <- tempfile(fileext = ".png")
grDevices::png(mask_file, width = 1200, height = 500, res = 144, bg = "white")
graphics::par(mar = rep(0, 4), family = font_family)
graphics::plot.new()
graphics::plot.window(c(-0.5, 0.5), c(-0.5, 0.5), xaxs = "i", yaxs = "i")
font_size <- 0.94 / graphics::strwidth(word, font = 2)
graphics::text(0, 0, word, font = 2, cex = font_size, adj = c(0.5, 0.5))
invisible(grDevices::dev.off())
mask_image <- png::readPNG(mask_file)
unlink(mask_file)
gray <- if (length(dim(mask_image)) == 2L) mask_image else mask_image[,, 1]
letter_raster <- terra::rast(ifelse(gray < 0.5, 1, NA_real_))
terra::ext(letter_raster) <- terra::ext(
  word_bounds["left"],
  word_bounds["right"],
  word_bounds["bottom"],
  word_bounds["top"]
)
terra::crs(letter_raster) <- ""
letter_shapes <- st_as_sf(terra::as.polygons(letter_raster, aggregate = TRUE))

# Choose point locations -------------------------------------------------------
set.seed(seed)
point_layout <- match.arg(point_layout, c("random", "hexagonal"))
if (point_layout == "random") {
  dots <- tibble(
    x = runif(
      2 * point_count,
      -sqrt(3) * radius / 2,
      sqrt(3) * radius / 2
    ),
    y = runif(2 * point_count, -radius, radius)
  ) |>
    filter(abs(y) <= radius - abs(x) / sqrt(3)) |>
    slice_head(n = point_count)
} else {
  domain <- st_as_sfc(st_bbox(c(
    xmin = -sqrt(3) * radius / 2,
    ymin = -radius,
    xmax = sqrt(3) * radius / 2,
    ymax = radius
  )))
  grid_spacing <- radius * sqrt(3 / point_count)
  grid <- st_make_grid(
    domain,
    cellsize = grid_spacing,
    square = FALSE,
    what = "centers"
  )
  xy <- st_coordinates(grid)
  dots <- tibble(x = xy[, 1], y = xy[, 2]) |>
    filter(
      abs(x) <= sqrt(3) * radius / 2,
      abs(y) <= radius - abs(x) / sqrt(3)
    )
}


# Assign continuous resources by noisy spatial utility -------------------------
n_points <- nrow(dots)
point_shapes <- st_as_sf(dots, coords = c("x", "y"), remove = FALSE)
in_word <- lengths(st_intersects(point_shapes, letter_shapes)) > 0L

dots <- dots |>
  mutate(
    in_word = lengths(st_intersects(point_shapes, letter_shapes)) > 0L,
    log_distance = log(sqrt(x^2 + y^2 + 0.1) / radius),
    center_score = ifelse(in_word, 0, log_distance),
    letter_score = as.numeric(in_word),
    score = letter_strength *
      scale(letter_score) +
      center_strength * scale(center_score),
    score = as.numeric(scale(score)) + noise * rnorm(n_points)
  ) |>
  arrange(score) |>
  mutate(z = sort(rbeta(n_points, beta_a, beta_b)))

# Cap the display scale only; retain the original resource values.
colour_cap <- unname(quantile(dots$z, colour_quantile))
dots <- dots |> mutate(colour_value = pmin(z, colour_cap))
dot_plot <- ggplot(dots, aes(x, y, colour = z)) +
  geom_point(size = point_size, stroke = 0) +
  scale_colour_viridis_c() +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

# Save, then show the logo in an interactive plot pane --------------------------
# hexSticker places the zero-centered subplot in its own frame at (1, 1).
sticker_file <- tempfile(fileext = ".png")
logo <- hexSticker::sticker(
  dot_plot,
  package = "",
  p_family = font_family,
  u_family = font_family,
  s_x = 1,
  s_y = 1,
  s_width = 2,
  s_height = 2,
  h_fill = background,
  h_color = border_colour,
  h_size = border_size,
  filename = sticker_file,
  dpi = output_dpi
)
logo <- logo + theme(plot.margin = margin(2, 2, 2, 2, unit = "mm"))
hexSticker::save_sticker(output, logo, dpi = output_dpi)
unlink(sticker_file)
if (!interactive()) {
  invisible(grDevices::dev.off(drawing_device))
  unlink(drawing_file)
}
if (interactive() && show_logo) {
  print(logo)
}
message("Saved ", n_points, " ", point_layout, " points to: ", output)
