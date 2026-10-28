# Regenerate the GitHub README logo with Rscript project/create-logo.R,
# or source("project/create-logo.R") from the repository root in RStudio.
# Edit the dot-resources parameters below to change the design.

local({
  script <- grep("^--file=", commandArgs(), value = TRUE)
  root <- if (length(script)) {
    file.path(dirname(sub("^--file=", "", script[1L])), "..")
  } else if (file.exists("DESCRIPTION")) {
    "."
  } else {
    ".."
  }
  description <- file.path(root, "DESCRIPTION")
  if (!file.exists(description) || read.dcf(description)[1L, "Package"] != "spseg") {
    stop("Run from the spseg repository root or project/ directory.", call. = FALSE)
  }
  root <- normalizePath(root, mustWork = TRUE)
  required <- c("ggplot2", "dplyr", "sf", "terra", "hexSticker", "png")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install the logo drawing packages: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  output <- file.path(root, "man", "figures", "logo.png")
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

  draw_logo <- function() {
    # Give font measurement a temporary device rather than creating Rplots.pdf.
    drawing_file <- tempfile(fileext = ".pdf")
    grDevices::pdf(drawing_file)
    drawing_device <- grDevices::dev.cur()
    on.exit({
      if (drawing_device %in% grDevices::dev.list()) grDevices::dev.off(drawing_device)
      unlink(drawing_file)
    }, add = TRUE)

    word_mask <- function(word = "spseg", width = 1200, height = 500) {
      file <- tempfile(fileext = ".png")
      on.exit(unlink(file), add = TRUE)
      grDevices::png(file, width = width, height = height, res = 144, bg = "white")
      device <- grDevices::dev.cur()
      on.exit(
        {
          if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
        },
        add = TRUE
      )
      graphics::par(mar = rep(0, 4), family = "sans")
      graphics::plot.new()
      graphics::plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
      size <- 0.94 / graphics::strwidth(word, font = 2)
      graphics::text(0.5, 0.5, word, font = 2, cex = size, adj = c(0.5, 0.5))
      grDevices::dev.off(device)
      image <- png::readPNG(file)
      gray <- if (length(dim(image)) == 2L) image else image[,, 1]
      gray < 0.5
    }

    mask <- word_mask()

    word_bounds <- c(left = 0.19, right = 1.81, bottom = 0.58, top = 1.58)
    letter_raster <- terra::rast(ifelse(mask, 1, NA_real_))
    terra::ext(letter_raster) <- terra::ext(
      word_bounds["left"],
      word_bounds["right"],
      word_bounds["bottom"],
      word_bounds["top"]
    )
    terra::crs(letter_raster) <- ""
    letter_shapes <- sf::st_as_sf(terra::as.polygons(
      letter_raster,
      aggregate = TRUE
    ))
    stopifnot(all(sf::st_is_valid(letter_shapes)))

    set.seed(2026)
    point_count <- 10000L
    noise <- 0.4
    letter_strength <- 2
    center_strength <- 2
    radius <- 0.94

    dots <- dplyr::tibble(
      x = runif(
        2 * point_count,
        1 - sqrt(3) * radius / 2,
        1 + sqrt(3) * radius / 2
      ),
      y = runif(2 * point_count, 1 - radius, 1 + radius)
    ) |>
      dplyr::filter(abs(y - 1) <= radius - abs(x - 1) / sqrt(3)) |>
      dplyr::slice_head(n = point_count)

    point_shapes <- sf::st_as_sf(dots, coords = c("x", "y"), remove = FALSE)
    in_word <- lengths(sf::st_intersects(point_shapes, letter_shapes)) > 0L

    resource_mean <- 10
    resource_gini <- 0.3
    resource_sd <- sqrt(2) * qnorm((1 + resource_gini) / 2)
    resources <- sort(rlnorm(
      point_count,
      log(resource_mean) - resource_sd^2 / 2,
      resource_sd
    ))
    distance <- with(dots, sqrt((x - 1)^2 + (y - 1)^2)) / radius
    center_score <- (!in_word) * log(distance - 0.99 * min(distance))
    letter_score <- in_word
    score <- letter_strength *
      scale(letter_score) +
      center_strength * scale(center_score)

    score <- scale(score) + noise * rnorm(point_count)

    assigned <- numeric(point_count)
    assigned[order(score)] <- resources
    dots <- dots |>
      dplyr::mutate(in_word = in_word, utility = score, resources = assigned)

    # Bound the display scale, not the actual resource distribution.
    colour_cap <- unname(quantile(resources, 0.9))
    dots <- dots |> dplyr::mutate(colour_value = pmin(resources, colour_cap))

    dot_plot <- ggplot2::ggplot(dots, ggplot2::aes(x, y, colour = colour_value)) +
      ggplot2::geom_point(size = 0.38, stroke = 0) +
      ggplot2::scale_colour_viridis_c(option = "D", limits = c(0, colour_cap)) +
      ggplot2::coord_fixed(xlim = c(0, 2), ylim = c(0, 2), expand = FALSE) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(0, 0, 0, 0))

    sticker_file <- tempfile(fileext = ".png")
    logo <- hexSticker::sticker(
      dot_plot,
      package = "",
      p_family = "sans",
      u_family = "sans",
      s_x = 1,
      s_y = 1,
      s_width = 2,
      s_height = 2,
      h_fill = "#101519",
      h_color = "#303940",
      h_size = 0.65,
      filename = sticker_file,
      dpi = 300
    )
    logo <- logo + ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "mm"))
    hexSticker::save_sticker(output, logo, dpi = 300)
    unlink(sticker_file)
    message("Saved logo: ", output)
    invisible(output)
  }
  draw_logo()
})
