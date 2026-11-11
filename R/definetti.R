# ---------------------------------------------------------------------------
# theme_void() base (the simplex has no x/y axes) plus the text/legend/caption
# elements.  Every argument is a plain value; "..." is passed on to theme() so
# anything not listed can still be overridden.
#   legend_position: "upper-right" (default), "upper-left", "lower-right",
#     "lower-left" (inset, single box), or any ggplot2 position
#     ("right", "bottom", "none", ...).
#' @import ggplot2
#' @importFrom utils head
NULL

utils::globalVariables(c("x", "y", "xend", "yend", ".a", ".lty", ".x", ".y", ".colour",
                         ".shape", ".alpha", ".size", ".label", ".g", "hj", "vj", ".lab",
                         "angle", ".view", ".ga", ".glty"))

#' ggplot2 theme for de Finetti plots
#'
#' A `theme_void()`-based theme (the simplex has no x/y axes) with text, legend,
#' caption and background elements. Used by [definetti_plot()]; can also be added
#' to any ggplot.
#'
#' @param base_size,base_family,base_face Base font size (pt), family and face.
#' @param title_size,subtitle_size,caption_size,legend_text_size,legend_title_size,strip_size
#'   Font sizes in points; `NULL` derives them from `base_size`.
#' @param title_face Face of the title.
#' @param title_hjust Title/subtitle justification.
#' @param caption_hjust Caption justification (0 = left).
#' @param legend_position One of `"upper-right"`, `"upper-left"`, `"lower-right"`,
#'   `"lower-left"` (single inset box) or any ggplot2 position such as `"right"`,
#'   `"bottom"`, `"none"`.
#' @param background Plot background colour.
#' @param ... Further arguments passed to [ggplot2::theme()].
#' @return A ggplot2 theme.
#' @export
theme_definetti <- function(base_size = 12, base_family = "", base_face = "plain",
                            title_size = NULL, title_face = "bold", title_hjust = 0.5,
                            subtitle_size = NULL, caption_size = NULL,
                            caption_hjust = 0,
                            legend_position = "upper-right",
                            legend_text_size = NULL, legend_title_size = NULL,
                            strip_size = NULL, background = "white", ...) {
  nz <- function(x, d) if (is.null(x)) d else x
  th <- theme_void(base_size = base_size, base_family = base_family) +
    theme(
      text             = element_text(face = base_face),
      plot.title       = element_text(size = nz(title_size, base_size * 1.2),
                                      face = title_face, hjust = title_hjust,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = nz(subtitle_size, base_size),
                                      hjust = title_hjust, margin = margin(b = 4)),
      plot.caption     = element_text(size = nz(caption_size, base_size * 0.9),
                                      hjust = caption_hjust, margin = margin(t = 8)),
      plot.caption.position = "plot",
      plot.background  = element_rect(fill = background, colour = NA),
      plot.margin      = margin(8, 8, 8, 8),
      strip.text       = element_text(size = nz(strip_size, base_size), face = "bold",
                                      margin = margin(b = 4)),
      legend.text      = element_text(size = nz(legend_text_size, base_size * 0.85)),
      legend.title     = element_text(size = nz(legend_title_size, base_size * 0.9)),
      legend.key.size  = unit(0.9, "lines"),
      legend.box       = "vertical",
      legend.background = element_rect(fill = background, colour = "grey70",
                                       linewidth = 0.3),
      legend.margin    = margin(3, 4, 3, 4)
    )
  corners <- list("upper-right" = c(1, 1), "upper-left" = c(0, 1),
                  "lower-right" = c(1, 0), "lower-left" = c(0, 0))
  if (legend_position %in% names(corners)) {
    th <- th + theme(legend.position = "inside",
                     legend.position.inside = corners[[legend_position]],
                     legend.justification.inside = corners[[legend_position]])
  } else {
    th <- th + theme(legend.position = legend_position)
  }
  th + theme(...)
}

# ---------------------------------------------------------------------------
# Geometry helpers (internal)
# ---------------------------------------------------------------------------
# Vertex coordinates (rows = components, cols = x,y,z) with unit edge length.
# K = 4: apex up (component 1), base triangle below.
.df_vertices <- function(K) {
  if (K == 2) return(rbind(c(0, 0, 0), c(1, 0, 0)))
  if (K == 3) return(rbind(c(0, 0, 0), c(1, 0, 0), c(0.5, sqrt(3) / 2, 0)))
  if (K == 4) {
    R  <- sqrt(6) / 4                      # circumradius, unit edge
    r  <- 1 / sqrt(3)                      # base-triangle circumradius
    th <- c(-90, 30, 150) * pi / 180       # first base vertex faces the viewer
    return(rbind(c(0, 0, R), cbind(r * cos(th), r * sin(th), -R / 3)))
  }
  stop("definetti: K must be 2, 3 or 4 (got ", K, ")", call. = FALSE)
}

# Orthographic projection.  view = c(azimuth, elevation) in degrees.
# depth: larger = farther from the viewer.
.df_project <- function(XYZ, view) {
  a <- view[1] * pi / 180; e <- view[2] * pi / 180
  x1 <- XYZ[, 1] * cos(a) - XYZ[, 2] * sin(a)
  y1 <- XYZ[, 1] * sin(a) + XYZ[, 2] * cos(a)
  z1 <- XYZ[, 3]
  data.frame(x = x1, y = z1 * cos(e) + y1 * sin(e), depth = y1 * cos(e) - z1 * sin(e))
}

# Composition matrix (n x K) -> data.frame(x, y, depth, near); near in [0,1],
# 1 = nearest the viewer (only meaningful for K = 4).
.df_xy <- function(P, K, view = c(20, 20), V = .df_vertices(K)) {
  XYZ <- P %*% V
  if (K < 4) return(data.frame(x = XYZ[, 1], y = XYZ[, 2], depth = 0, near = 1))
  pr <- .df_project(XYZ, view)
  R  <- sqrt(6) / 4                        # every point lies within this radius
  pr$near <- pmin(1, pmax(0, (R - pr$depth) / (2 * R)))
  pr
}

# Text justification so a label sits outside the anchor, along direction (dx,dy).
.df_just <- function(dx, dy) {
  m <- pmax(abs(dx), abs(dy)); m[m == 0] <- 1
  data.frame(hj = 0.5 - 0.5 * dx / m, vj = 0.5 - 0.5 * dy / m)
}

# Validate (and optionally normalise) the frequency matrix.  Fails loudly.
.df_check_freqs <- function(P, cols, normalize, tol) {
  if (anyNA(P))
    stop("definetti: NA in frequency columns, rows: ",
         paste(head(which(!stats::complete.cases(P)), 10), collapse = ", "), call. = FALSE)
  if (any(P < -tol))
    stop("definetti: negative values in frequency columns, rows: ",
         paste(head(which(apply(P < -tol, 1, any)), 10), collapse = ", "), call. = FALSE)
  s <- rowSums(P)
  if (normalize) {
    if (any(s <= 0)) stop("definetti: normalize=TRUE but some rows sum to 0", call. = FALSE)
    return(P / s)
  }
  bad <- which(abs(s - 1) > tol)
  if (length(bad))
    stop("definetti: ", length(bad), " row(s) do not sum to 1 (sums range ",
         format(min(s[bad]), digits = 4), " to ", format(max(s[bad]), digits = 4),
         "; first rows: ", paste(head(bad, 10), collapse = ", "),
         "). Use normalize = TRUE to rescale.", call. = FALSE)
  P
}

# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------
# data, freq_cols : data.frame and the names of its K composition columns
# Mappings        : colour_by, shape_by, label_by (column names);
#                   connect_by (column naming a trajectory id; points joined in
#                   row order within id, optional arrow)
# Text            : title, subtitle, caption (+ figure_number -> "Figure N: "),
#                   vertex_labels (default = freq_cols), axis_titles (one per
#                   component; K = 2 and 3 only), legend_title, shape_title
# Ticks/grid      : tick_labels, n_ticks (intervals per axis), show_grid
# 4 components    : view = c(azimuth, elevation) or list of such pairs
#                   (-> facets), depth_cue, depth_alpha_far
# Points          : point_size, point_alpha, point_colour, point_shape, palette
# Font sizes (pt) : base_size, vertex_label_size, tick_label_size,
#                   axis_title_size, label_size; faces/families likewise
# Extras          : show_hw (K = 3: columns are homozygote 1, heterozygote,
#                   homozygote 2; draws the Hardy-Weinberg curve)
# Input handling  : normalize (rescale rows to sum 1), tol (sum tolerance)
#' de Finetti (simplex) plots for 2, 3 or 4 components
#'
#' Plots compositions (allele frequencies, or the ancestry fractions of an
#' individual) whose K components sum to 1. K = 2 gives a line, K = 3 a triangle,
#' K = 4 a tetrahedron projected to the page (or, with `colour_component`, a
#' triangle with the fourth component as colour). Returns a ggplot object.
#'
#' @param data A data.frame.
#' @param freq_cols Names of the K (2 to 4) composition columns.
#' @param colour_by,shape_by,label_by Columns mapped to colour, shape or text labels.
#' @param connect_by Column giving a trajectory id; points are joined in row order within id.
#' @param title,subtitle,caption Plot text.
#' @param figure_number If given, the caption
#' is prefixed `"Figure N: "`.
#' @param caption_wrap Characters per caption line;
#' `NULL` derives it from `fig_width` and the caption size.
#' @param fig_width Intended
#'   figure width in inches (used only for caption wrapping).
#' @param vertex_labels Corner labels (default `freq_cols`).
#' @param axis_titles Optional titles, one per component (K = 2, 3).
#' @param legend_title,shape_title Legend titles.
#' @param tick_labels Draw tick labels?
#' @param n_ticks Tick intervals per axis.
#' @param label_at Levels that get a label.
#' @param show_grid,grid_at Light gridlines
#' at the given levels.
#' @param edge_ticks Tick marks on the tetrahedron edges (K = 4).
#' @param face_grid K = 4: draw the `grid_at` levels as lines on the faces of the tetrahedron (faces turned away from the viewer are dashed and faint).
#' @param colour_component Name of one of `freq_cols` to show as colour instead of
#'   geometry (K = 3, 4). The component's corner is placed at the centre of the others.
#' @param colour_scale A ggplot2 colour scale for `colour_component` (overrides the next two).
#' @param colour_option Viridis palette name ("viridis", "magma", "cividis", ...),
#'   uniform in perceived lightness.
#' @param colour_range Two colours for a ramp interpolated in CIELAB, not RGB.
#' @param view `c(azimuth, elevation)` in degrees, or a list of them (facets), for K = 4.
#' @param depth_cue,depth_alpha_far Fade and shrink far points (K = 4).
#' @param point_size,point_alpha,point_colour,point_shape,palette Point appearance;
#'   `palette` is a named or unnamed vector of discrete group colours.
#' @param label_size Size (pt) of `label_by` text.
#' @param path_colour,path_linewidth,path_arrow Trajectory appearance.
#' @param show_hw,hw_colour,hw_linetype K = 3 only: the Hardy-Weinberg curve, reading the
#'   columns as homozygote 1, heterozygote, homozygote 2.
#' @param frame_colour,frame_linewidth,grid_colour,grid_linetype,grid_linewidth Line styling.
#' @param vertex_label_size,vertex_label_face,tick_label_size,axis_title_size,axis_title_face
#'   Font sizes (pt) and faces.
#' @param vertex_offset,tick_offset,axis_title_offset Label distances (simplex units).
#' @param jitter_height,seed K = 2 vertical jitter and its seed (the global RNG is untouched).
#' @param margin Padding around the plot in simplex units.
#' @param normalize Rescale rows to sum to 1 instead of failing.
#' @param tol Tolerance
#'   for the row-sum check.
#' @param base_size,base_family,base_face See [theme_definetti()].
#' @param legend_position As in [theme_definetti()]; default `"upper-right"`, or `"bottom"` for K = 2.
#' @param ... Further arguments to [theme_definetti()].
#' @return A ggplot object.
#' @examples
#' f <- read.csv(system.file("extdata", "individuals.ancestry_k4.csv", package = "definetti"))
#' definetti_plot(f, paste0("Q", 1:4), colour_by = "group")
#' definetti_plot(f, paste0("Q", 1:4), colour_component = "Q4")
#' @export
definetti_plot <- function(data, freq_cols,
                           colour_by = NULL, shape_by = NULL, label_by = NULL,
                           connect_by = NULL,
                           title = NULL, subtitle = NULL, caption = NULL,
                           figure_number = NULL, caption_wrap = NULL, fig_width = 6.5,
                           vertex_labels = freq_cols, axis_titles = NULL,
                           legend_title = NULL, shape_title = NULL,
                           tick_labels = TRUE, n_ticks = 10, label_at = c(0, 0.5, 1),
                           show_grid = TRUE, grid_at = 0.5, edge_ticks = TRUE, face_grid = FALSE,
                           colour_component = NULL, colour_scale = NULL,
                           colour_option = "viridis", colour_range = NULL,
                           view = c(20, 20), depth_cue = TRUE, depth_alpha_far = 0.35,
                           point_size = 2.5, point_alpha = 0.8,
                           point_colour = "black", point_shape = 16, palette = NULL,
                           label_size = NULL,
                           path_colour = "grey40", path_linewidth = 0.4, path_arrow = TRUE,
                           show_hw = FALSE, hw_colour = "grey30", hw_linetype = "dashed",
                           frame_colour = "black", frame_linewidth = 0.6,
                           grid_colour = "grey80", grid_linetype = "solid",
                           grid_linewidth = 0.3,
                           vertex_label_size = NULL, vertex_label_face = "bold",
                           tick_label_size = NULL,
                           axis_title_size = NULL, axis_title_face = "plain",
                           vertex_offset = 0.06, tick_offset = 0.035,
                           axis_title_offset = 0.11, jitter_height = 0.04,
                           seed = 1, margin = 0.08,
                           normalize = FALSE, tol = 1e-6,
                           base_size = 12, base_family = "", base_face = "plain",
                           legend_position = NULL, ...) {
  if (is.null(legend_position))
    legend_position <- if (length(freq_cols) == 2) "bottom" else "upper-right"
  # ---- input checks (all failures are errors, never silent) ----------------
  K <- length(freq_cols)
  if (!K %in% 2:4) stop("definetti: need 2, 3 or 4 frequency columns (got ", K, ")", call. = FALSE)
  data <- as.data.frame(data)
  need <- c(freq_cols, colour_by, shape_by, label_by, connect_by)
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("definetti: column(s) not in data: ", paste(miss, collapse = ", "), call. = FALSE)
  if (length(vertex_labels) != K) stop("definetti: vertex_labels needs length ", K, call. = FALSE)
  if (!is.null(axis_titles) && length(axis_titles) != K)
    stop("definetti: axis_titles needs length ", K, call. = FALSE)
  if (!is.null(axis_titles) && K == 4)
    message("definetti: axis_titles are not used for K = 4 (vertex labels only)")
  if (show_hw && K != 3) stop("definetti: show_hw is only defined for K = 3", call. = FALSE)
  P <- .df_check_freqs(as.matrix(data[, freq_cols, drop = FALSE]), freq_cols, normalize, tol)
  n <- nrow(P)

  # ---- colour_component: view the simplex along one component's vertex -----
  # The chosen component is dropped from the geometry (its vertex is placed at
  # the centroid of the remaining ones, i.e. the view looking down that vertex)
  # and its frequency is mapped to a continuous colour scale.
  V_geo <- NULL; cc_val <- NULL; cc_label <- NULL
  if (!is.null(colour_component)) {
    if (K < 3) stop("definetti: colour_component needs K = 3 or 4", call. = FALSE)
    if (!is.null(colour_by)) stop("definetti: use either colour_by or colour_component", call. = FALSE)
    cc <- match(colour_component, freq_cols)
    if (is.na(cc)) stop("definetti: colour_component must be one of freq_cols", call. = FALSE)
    keep <- setdiff(seq_len(K), cc)
    cc_val <- P[, cc]; cc_label <- vertex_labels[cc]
    P <- P[, c(keep, cc), drop = FALSE]
    vertex_labels <- vertex_labels[keep]
    if (!is.null(axis_titles)) axis_titles <- axis_titles[keep]
    K <- K - 1
    Vk <- .df_vertices(K); V_geo <- rbind(Vk, colMeans(Vk))
  }

  # ---- font sizes in points -> geom_text mm --------------------------------
  nz <- function(x, d) if (is.null(x)) d else x
  vl_pt <- nz(vertex_label_size, base_size * 1.1)
  tk_pt <- nz(tick_label_size,   base_size * 0.8)
  at_pt <- nz(axis_title_size,   base_size)
  lb_pt <- nz(label_size,        base_size * 0.75)
  mm <- function(pt) pt / ggplot2::.pt

  # ---- aesthetic columns ---------------------------------------------------
  aes_col <- function(v, discrete) {
    if (is.null(v)) return(rep("all", n))
    x <- data[[v]]
    if (discrete || !is.numeric(x)) factor(x) else x
  }
  colour_v <- if (is.null(cc_val)) aes_col(colour_by, FALSE) else cc_val
  shape_v  <- aes_col(shape_by, TRUE)
  label_v  <- if (is.null(label_by)) NA_character_ else as.character(data[[label_by]])
  path_v   <- if (is.null(connect_by)) NA else data[[connect_by]]

  # ---- views (only K = 4 has more than one) --------------------------------
  views <- if (K == 4) { if (is.list(view)) view else list(view) } else list(c(0, 0))
  vnames <- if (K == 4) vapply(views, function(v) sprintf("Azimuth %g, elevation %g", v[1], v[2]), "")
            else "view"
  vnames <- factor(vnames, levels = unique(vnames))

  pts <- list(); paths <- list(); frame_seg <- list(); vtx <- list(); anch <- list(); fseg <- list()
  grid_seg <- NULL; poly <- NULL; ticks <- NULL; ttl <- NULL; hw <- NULL

  # K = 2 vertical jitter, local seed so the user's RNG stream is untouched
  yj <- rep(0, n)
  if (K == 2 && jitter_height > 0) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seed); yj <- stats::runif(n, -jitter_height, jitter_height)
    if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
  }

  tt <- seq(0, 1, length.out = n_ticks + 1)           # tick levels incl. ends
  ti <- tt[-c(1, length(tt))]                         # interior levels
  lab_lv <- function(t) vapply(t, function(x) any(abs(x - label_at) < 1e-9), TRUE)
  tick_seg <- NULL
  fmt <- function(x) formatC(x, format = "g", digits = 2)

  for (iv in seq_along(views)) {
    vw <- views[[iv]]; vn <- vnames[iv]
    xy <- .df_xy(P, K, vw, if (is.null(V_geo)) .df_vertices(K) else V_geo)
    cue <- if (K == 4 && depth_cue) depth_alpha_far + (1 - depth_alpha_far) * xy$near else 1
    szf <- if (K == 4 && depth_cue) 0.6 + 0.4 * xy$near else 1
    d <- data.frame(.x = xy$x, .y = xy$y + yj, .depth = xy$depth,
                    .alpha = point_alpha * cue, .size = point_size * szf,
                    .colour = colour_v, .shape = shape_v, .label = label_v,
                    .path = path_v, .view = vn)
    paths[[iv]] <- d
    pts[[iv]] <- if (K == 4) d[order(-d$.depth), ] else d
    anch[[iv]] <- data.frame(x = d$.x, y = d$.y)

    V <- .df_vertices(K)
    if (K == 2) {
      frame_seg[[iv]] <- data.frame(x = 0, y = 0, xend = 1, yend = 0, .a = 1, .lty = "solid", .view = vn)
      tl <- tick_offset * 0.5
      tick_seg <- rbind(
        data.frame(x = tt, y = 0, xend = tt, yend = -tl, .view = vn),
        data.frame(x = tt, y = 0, xend = tt, yend =  tl, .view = vn))
      lt <- tt[lab_lv(tt)]
      if (tick_labels) ticks <- rbind(
        data.frame(x = lt, y = -tick_offset, .lab = fmt(lt), hj = 0.5, vj = 1, .view = vn),
        data.frame(x = lt, y =  tick_offset, .lab = fmt(1 - lt), hj = 0.5, vj = 0, .view = vn))
      vx <- c(-vertex_offset, 1 + vertex_offset)
      vtx[[iv]] <- data.frame(x = vx, y = 0, .lab = vertex_labels, hj = c(1, 0), vj = 0.5, .view = vn)
      if (!is.null(axis_titles)) ttl <- data.frame(
        x = 0.5, y = c(axis_title_offset, -axis_title_offset), .lab = axis_titles,
        angle = 0, hj = 0.5, vj = c(0, 1), .view = vn)
      anch[[length(anch) + 1]] <- data.frame(
        x = c(vx, 0, 1), y = c(0, 0, axis_title_offset, -axis_title_offset))
    }
    if (K == 3) {
      poly <- data.frame(x = V[, 1], y = V[, 2], .view = vn)
      cen <- colMeans(V[, 1:2])
      mk <- function(i, t, j) { b <- matrix(0, length(t), 3); b[, i] <- t; b[, j] <- 1 - t
                                b %*% V[, 1:2, drop = FALSE] }
      gl <- list(); tl <- list(); tt_ <- list(); tk <- list()
      gi <- grid_at[grid_at > 0 & grid_at < 1]
      lt <- tt[lab_lv(tt)]
      for (i in 1:3) {
        ni <- i %% 3 + 1; pi_ <- (i + 1) %% 3 + 1
        A <- mk(i, gi, ni); B <- mk(i, gi, pi_)
        gl[[i]] <- data.frame(x = A[, 1], y = A[, 2], xend = B[, 1], yend = B[, 2], .view = vn)
        Tm <- mk(i, ti, ni)
        dd <- V[ni, 1:2] - V[i, 1:2]; nn <- c(dd[2], -dd[1]); nn <- nn / sqrt(sum(nn^2))
        mid <- (V[ni, 1:2] + V[i, 1:2]) / 2
        if (sum(nn * (mid - cen)) < 0) nn <- -nn
        js <- .df_just(nn[1], nn[2])
        tk[[i]] <- data.frame(x = Tm[, 1], y = Tm[, 2], xend = Tm[, 1] + nn[1] * tick_offset * 0.5,
                              yend = Tm[, 2] + nn[2] * tick_offset * 0.5, .view = vn)
        L <- mk(i, lt, ni)
        tl[[i]] <- data.frame(x = L[, 1] + nn[1] * tick_offset, y = L[, 2] + nn[2] * tick_offset,
                              .lab = fmt(lt), hj = js$hj, vj = js$vj, .view = vn)
        ang <- atan2(dd[2], dd[1]) * 180 / pi
        if (ang > 90) ang <- ang - 180; if (ang <= -90) ang <- ang + 180
        tt_[[i]] <- data.frame(x = mid[1] + nn[1] * axis_title_offset,
                               y = mid[2] + nn[2] * axis_title_offset,
                               .lab = if (is.null(axis_titles)) NA_character_ else axis_titles[i],
                               angle = ang, hj = 0.5, vj = 0.5, .view = vn)
      }
      if (show_grid && length(gi)) grid_seg <- do.call(rbind, gl)
      tick_seg <- do.call(rbind, tk)
      if (tick_labels) ticks <- do.call(rbind, tl)
      if (!is.null(axis_titles)) ttl <- do.call(rbind, tt_)
      dv <- cbind(V[, 1] - cen[1], V[, 2] - cen[2]); dv <- dv / sqrt(rowSums(dv^2))
      js <- .df_just(dv[, 1], dv[, 2])
      vtx[[iv]] <- data.frame(x = V[, 1] + dv[, 1] * vertex_offset, y = V[, 2] + dv[, 2] * vertex_offset,
                              .lab = vertex_labels, hj = js$hj, vj = js$vj, .view = vn)
      anch[[length(anch) + 1]] <- rbind(
        vtx[[iv]][, c("x", "y")], if (!is.null(ticks)) ticks[, c("x", "y")],
        if (!is.null(ttl)) ttl[, c("x", "y")])
      if (show_hw) {
        p <- seq(0, 1, length.out = 201)
        h <- cbind(p^2, 2 * p * (1 - p), (1 - p)^2) %*% V[, 1:2]
        hw <- data.frame(x = h[, 1], y = h[, 2])
      }
    }
    if (K == 4) {
      pv <- .df_project(V, vw)
      pr <- t(utils::combn(4, 2))
      hid <- apply(pr, 1, function(e) { o <- setdiff(1:4, e); all(pv$depth[o] <= 0) })
      mid_d <- (pv$depth[pr[, 1]] + pv$depth[pr[, 2]]) / 2
      near <- pmin(1, pmax(0, (sqrt(6) / 4 - mid_d) / (sqrt(6) / 2)))
      frame_seg[[iv]] <- data.frame(
        x = pv$x[pr[, 1]], y = pv$y[pr[, 1]], xend = pv$x[pr[, 2]], yend = pv$y[pr[, 2]],
        .a = if (depth_cue) depth_alpha_far + (1 - depth_alpha_far) * near else 1,
        .lty = ifelse(hid, "dashed", "solid"), .view = vn)
      if (face_grid) {
        # lines of constant component frequency drawn on each of the 4 faces;
        # a face is seen from outside iff the vertex opposite it is farther than the centre
        gj <- grid_at[grid_at > 0 & grid_at < 1]
        for (m in 1:4) {
          fv <- setdiff(1:4, m); vis <- pv$depth[m] > 0
          for (i in fv) for (t in gj) {
            o <- setdiff(fv, i)
            A <- t * V[i, ] + (1 - t) * V[o[1], ]; B <- t * V[i, ] + (1 - t) * V[o[2], ]
            q <- .df_project(rbind(A, B), vw)
            fseg[[length(fseg) + 1]] <- data.frame(
              x = q$x[1], y = q$y[1], xend = q$x[2], yend = q$y[2],
              .ga = if (vis) 1 else 0.45, .glty = if (vis) grid_linetype else "dashed", .view = vn)
          }
        }
      }
      if (edge_ticks) {
        tk <- list(); tl <- list(); lt <- ti[lab_lv(ti)]
        for (e in seq_len(nrow(pr))) {
          a <- c(pv$x[pr[e, 1]], pv$y[pr[e, 1]]); b <- c(pv$x[pr[e, 2]], pv$y[pr[e, 2]])
          dd <- b - a; nn <- c(dd[2], -dd[1]); nl <- sqrt(sum(nn^2))
          if (nl < 1e-8) next
          nn <- nn / nl; if (sum(nn * (a + b) / 2) < 0) nn <- -nn   # outward from centre
          Q <- cbind(a[1] + ti * dd[1], a[2] + ti * dd[2])
          tk[[e]] <- data.frame(x = Q[, 1], y = Q[, 2], xend = Q[, 1] + nn[1] * tick_offset * 0.5,
                                yend = Q[, 2] + nn[2] * tick_offset * 0.5, .view = vn)
          if (tick_labels && length(lt)) {
            Lq <- cbind(a[1] + lt * dd[1], a[2] + lt * dd[2]); js <- .df_just(nn[1], nn[2])
            tl[[e]] <- data.frame(x = Lq[, 1] + nn[1] * tick_offset, y = Lq[, 2] + nn[2] * tick_offset,
                                  .lab = fmt(lt), hj = js$hj, vj = js$vj, .view = vn)
          }
        }
        tick_seg <- rbind(tick_seg, do.call(rbind, tk))
        if (length(tl)) ticks <- rbind(ticks, do.call(rbind, tl))
      }
      dv <- cbind(pv$x, pv$y); nr <- sqrt(rowSums(dv^2))
      dv[nr < 1e-8, ] <- cbind(0, 1)[rep(1, sum(nr < 1e-8)), ]; nr[nr < 1e-8] <- 1
      dv <- dv / nr
      js <- .df_just(dv[, 1], dv[, 2])
      vtx[[iv]] <- data.frame(x = pv$x + dv[, 1] * vertex_offset, y = pv$y + dv[, 2] * vertex_offset,
                              .lab = vertex_labels, hj = js$hj, vj = js$vj, .view = vn)
      anch[[length(anch) + 1]] <- rbind(vtx[[iv]][, c("x", "y")], if (!is.null(ticks)) ticks[, c("x", "y")])
    }
  }
  pts <- do.call(rbind, pts); frame_seg <- do.call(rbind, frame_seg); vtx <- do.call(rbind, vtx)
  anch <- do.call(rbind, anch)
  fseg <- if (length(fseg)) do.call(rbind, fseg) else NULL
  # make room for the horizontal extent of the vertex labels (heuristic width)
  w <- nchar(vtx$.lab) * 0.022 * vl_pt / 12
  anch <- rbind(anch, data.frame(x = c(vtx$x - vtx$hj * w, vtx$x + (1 - vtx$hj) * w),
                                 y = rep(vtx$y, 2)))

  # ---- assemble ggplot -----------------------------------------------------
  p <- ggplot()
  if (!is.null(poly))
    p <- p + geom_polygon(data = poly, aes(x = x, y = y), fill = NA,
                          colour = frame_colour, linewidth = frame_linewidth)
  if (!is.null(grid_seg))
    p <- p + geom_segment(data = grid_seg, aes(x = x, y = y, xend = xend, yend = yend),
                          colour = grid_colour, linetype = grid_linetype, linewidth = grid_linewidth)
  if (!is.null(fseg))
    p <- p + geom_segment(data = fseg, aes(x = x, y = y, xend = xend, yend = yend,
                                           alpha = .ga, linetype = .glty),
                          colour = grid_colour, linewidth = grid_linewidth)
  if (!is.null(tick_seg))
    p <- p + geom_segment(data = tick_seg, aes(x = x, y = y, xend = xend, yend = yend),
                          colour = frame_colour, linewidth = frame_linewidth * 0.7)
  if (!is.null(hw))
    p <- p + geom_path(data = hw, aes(x = x, y = y), colour = hw_colour, linetype = hw_linetype,
                       linewidth = 0.5)
  if (K != 3)
    p <- p + geom_segment(data = frame_seg,
                          aes(x = x, y = y, xend = xend, yend = yend, alpha = .a, linetype = .lty),
                          colour = frame_colour, linewidth = frame_linewidth)
  if (!is.null(connect_by)) {
    pa <- do.call(rbind, paths)
    pa$.g <- interaction(pa$.path, pa$.view, drop = TRUE)
    p <- p + geom_path(data = pa, aes(x = .x, y = .y, group = .g), colour = path_colour,
                       linewidth = path_linewidth,
                       arrow = if (path_arrow) grid::arrow(length = grid::unit(0.08, "inches")) else NULL)
  }
  p <- p + geom_point(data = pts, aes(x = .x, y = .y, colour = .colour, shape = .shape,
                                      alpha = .alpha, size = .size), stroke = 0.6)
  if (!is.null(label_by))
    p <- p + geom_text(data = pts, aes(x = .x, y = .y, label = .label), hjust = 0, nudge_x = 0.012,
                       size = mm(lb_pt), family = base_family, check_overlap = TRUE)
  if (!is.null(ticks))
    p <- p + geom_text(data = ticks, aes(x = x, y = y, label = .lab, hjust = hj, vjust = vj),
                       size = mm(tk_pt), family = base_family, colour = "grey30")
  if (!is.null(ttl))
    p <- p + geom_text(data = ttl, aes(x = x, y = y, label = .lab, angle = angle, hjust = hj, vjust = vj),
                       size = mm(at_pt), family = base_family, fontface = axis_title_face)
  p <- p + geom_text(data = vtx, aes(x = x, y = y, label = .lab, hjust = hj, vjust = vj),
                     size = mm(vl_pt), family = base_family, fontface = vertex_label_face)

  # scales: identity for computed alpha/size/linetype; colour/shape as requested
  p <- p + scale_alpha_identity() + scale_size_identity() + scale_linetype_identity()
  if (!is.null(cc_val)) {
    # Perceptual colour scaling: viridis-family palettes are uniform in CIELAB
    # lightness; a custom colour_range = c(low, high) is interpolated in CIELAB
    # (space = "Lab"), never in RGB.
    p <- p + (if (!is.null(colour_scale)) colour_scale
              else if (!is.null(colour_range))
                scale_colour_gradient(low = colour_range[1], high = colour_range[2],
                                      limits = c(0, 1), space = "Lab")
              else scale_colour_viridis_c(option = colour_option, limits = c(0, 1)))
  } else if (is.null(colour_by)) {
    p <- p + scale_colour_manual(values = c(all = point_colour), guide = "none")
  } else if (!is.null(palette) && is.factor(colour_v)) {
    p <- p + scale_colour_manual(values = palette)
  }
  if (is.null(shape_by)) {
    p <- p + scale_shape_manual(values = c(all = point_shape), guide = "none")
  }
  p <- p + labs(colour = nz(legend_title, nz(colour_by, cc_label)), shape = nz(shape_title, shape_by))

  if (!is.null(caption)) {
    dots <- list(...)
    cap_pt <- if (!is.null(dots$caption_size)) dots$caption_size else base_size * 0.9
    cw <- if (is.null(caption_wrap)) floor(fig_width * 72 / (0.6 * cap_pt)) else caption_wrap
    caption <- paste(strwrap(caption, cw), collapse = "\n")
    if (!is.null(figure_number)) caption <- paste0("Figure ", figure_number, ": ", caption)
  }
  p <- p + labs(title = title, subtitle = subtitle, caption = caption)

  xl <- range(anch$x, na.rm = TRUE) + c(-1, 1) * margin
  yl <- range(anch$y, na.rm = TRUE) + c(-1, 1) * margin
  p <- p + coord_fixed(xlim = xl, ylim = yl, expand = FALSE, clip = "off")
  if (K == 4 && length(views) > 1) p <- p + facet_wrap(vars(.view))
  p + theme_definetti(base_size = base_size, base_family = base_family, base_face = base_face,
                      legend_position = legend_position, ...)
}

.df_need <- function(K, freq_cols) {
  if (length(freq_cols) != K)
    stop("definetti: this function needs exactly ", K, " frequency columns (got ",
         length(freq_cols), ")", call. = FALSE)
}

#' @rdname definetti_plot
#' @export
definetti_plot_2 <- function(data, freq_cols, ...) { .df_need(2, freq_cols); definetti_plot(data, freq_cols, ...) }
#' @rdname definetti_plot
#' @export
definetti_plot_3 <- function(data, freq_cols, ...) { .df_need(3, freq_cols); definetti_plot(data, freq_cols, ...) }
#' @rdname definetti_plot
#' @export
definetti_plot_4 <- function(data, freq_cols, ...) { .df_need(4, freq_cols); definetti_plot(data, freq_cols, ...) }
