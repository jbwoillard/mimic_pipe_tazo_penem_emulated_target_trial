# =============================================================================
# 82_dag_figure.R — DAG redessiné (C1, docs/v3.3_plan.md) : coordonnées
#   explicites, boîtes rectangulaires, segments coupés au bord des boîtes,
#   une pointe de flèche par arête ; QA : arêtes tracées == arêtes dagitty,
#   extrémités sur les bords, aucun segment ne traverse un autre nœud.
#   Le modèle (nœuds, arêtes) est inchangé : DAG_NODES_V32 / DAG_EDGES_V32.
# =============================================================================

DAG_COORDS_V33 <- data.table::data.table(
  node = c("Treatment", "AKI", "Comorbidity", "BaselineRenal", "CreatTrajectory", "Severity", "InitContext",
           "CareUnit", "Infection", "Nephrotoxins", "Vancomycin", "Period", "PrescriberPreference", "Allergy", "Colonization"),
  x = 1.4 * c(0, 14, 7, 11, 8.8, 2.7, -1.5, 5, 1.2, 11.8, 12.5, 7, -4.5, -4.5, -4.5),
  y = 1.4 * c(0, 0, 7.5, 5.3, 1.2, 4.8, 4, -5.2, -3.6, -2.8, -5, -1.0, 2, 0, -2.4))

DAG_CHAR_W_V33 <- 0.155; DAG_LINE_H_V33 <- 0.38; DAG_PAD_V33 <- 0.45

#' Géométrie : boîtes (centre, demi-largeur, demi-hauteur) et segments coupés
dag_geometry_v33 <- function(g = build_dag_v32(), coords = DAG_COORDS_V33) {
  nd <- merge(coords, DAG_NODES_V32[, .(node, label, measured)], by = "node")
  nd[, `:=`(hw = (vapply(strsplit(label, "\n"), function(s) max(nchar(s)), numeric(1)) * DAG_CHAR_W_V33 + DAG_PAD_V33) / 2,
            hh = (lengths(strsplit(label, "\n")) * DAG_LINE_H_V33 + DAG_PAD_V33) / 2)]
  e <- data.table::as.data.table(dagitty::edges(g))[e == "->", .(from = as.character(v), to = as.character(w))]
  # point de sortie du segment centre->centre sur le bord du rectangle (+ marge)
  clip <- function(cx, cy, hw, hh, dx, dy, m = 0.05) { t <- min(if (dx != 0) (hw + m) / abs(dx) else Inf, if (dy != 0) (hh + m) / abs(dy) else Inf)
    c(cx + t * dx, cy + t * dy) }
  seg <- e[, { a <- nd[node == from]; b <- nd[node == to]; dx <- b$x - a$x; dy <- b$y - a$y
    p0 <- clip(a$x, a$y, a$hw, a$hh, dx, dy); p1 <- clip(b$x, b$y, b$hw, b$hh, -dx, -dy)
    .(x = p0[1], y = p0[2], xend = p1[1], yend = p1[2], latent = !a$measured) }, by = .(from, to)]
  list(nodes = nd, edges = seg)
}

.seg_hits_box_v33 <- function(x0, y0, x1, y1, cx, cy, hw, hh) {
  # échantillonnage fin du segment : un point strictement dans le rectangle = traversée
  t <- seq(0, 1, length.out = 400); px <- x0 + t * (x1 - x0); py <- y0 + t * (y1 - y0)
  any(abs(px - cx) < hw & abs(py - cy) < hh)
}

#' QA du DAG : liste d'arêtes, extrémités, pointes de flèche, traversées
qa_dag_v33 <- function(p, geo, g = build_dag_v32()) {
  b <- ggplot2::ggplot_build(p)
  li <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1)))
  drawn <- data.table::rbindlist(lapply(li, function(i) data.table::as.data.table(b$data[[i]])[, .(x, y, xend, yend)]))
  arrows_ok <- all(vapply(p$layers[li], function(l) !is.null(l$geom_params$arrow), logical(1)))
  model <- data.table::as.data.table(dagitty::edges(g))[e == "->", .(key = paste(v, w, sep = "->"))]$key
  # correspondance segment tracé -> arête : chaque segment doit partir du bord de son nœud d'origine et arriver au bord du nœud de destination
  nd <- geo$nodes
  near <- function(px, py) { d <- nd[, pmax(abs(px - x) - hw, abs(py - y) - hh)]; i <- which.min(d); list(node = nd$node[i], dist = d[i]) }
  m <- drawn[, { a <- near(x, y); z <- near(xend, yend); .(from = a$node, to = z$node, d0 = a$dist, d1 = z$dist) }, by = seq_len(nrow(drawn))]
  drawn_keys <- paste(m$from, m$to, sep = "->")
  ge <- geo$edges
  hits <- vapply(seq_len(nrow(ge)), function(k) { s <- ge[k]; o <- nd[!node %in% c(s$from, s$to)]
    paste(o$node[vapply(seq_len(nrow(o)), function(j) .seg_hits_box_v33(s$x, s$y, s$xend, s$yend, o$x[j], o$y[j], o$hw[j], o$hh[j]), logical(1))], collapse = ",") }, "")
  cross <- ge[, .(from, to)][, hits := hits]
  list(n_model = length(model), n_drawn = nrow(drawn), same_edges = setequal(model, drawn_keys) && !anyDuplicated(drawn_keys),
       endpoints_on_borders = all(m$d0 > 0 & m$d0 < 0.15 & m$d1 > 0 & m$d1 < 0.15), arrows = arrows_ok,
       crossings = if (any(cross$hits != "")) cross[hits != "", paste0(from, "->", to, " crosses ", hits)] else character(0))
}

plot_dag_v33 <- function(geo) {
  nd <- geo$nodes; ed <- geo$edges
  nd[, role := fifelse(node %in% c("Treatment", "AKI"), "key", fifelse(measured, "measured", "latent"))]
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = ed[latent == FALSE], ggplot2::aes(x, y, xend = xend, yend = yend), linewidth = 0.3, colour = "grey25",
                          arrow = grid::arrow(length = grid::unit(4, "pt"), type = "closed")) +
    ggplot2::geom_segment(data = ed[latent == TRUE], ggplot2::aes(x, y, xend = xend, yend = yend), linewidth = 0.3, colour = "grey45", linetype = "dashed",
                          arrow = grid::arrow(length = grid::unit(4, "pt"), type = "closed")) +
    ggplot2::geom_rect(data = nd, ggplot2::aes(xmin = x - hw, xmax = x + hw, ymin = y - hh, ymax = y + hh, fill = role, linetype = role),
                       colour = "grey20", linewidth = 0.3) +
    ggplot2::geom_text(data = nd, ggplot2::aes(x, y, label = label, colour = role), size = 2.6, lineheight = 0.9) +
    ggplot2::scale_fill_manual(values = c(key = "#d9d9d9", measured = "white", latent = "white")) +
    ggplot2::scale_colour_manual(values = c(key = "black", measured = "black", latent = "grey40")) +
    ggplot2::scale_linetype_manual(values = c(key = "solid", measured = "solid", latent = "dashed")) +
    ggplot2::coord_fixed(clip = "off") + ggplot2::theme_void(base_size = 9) +
    ggplot2::labs(caption = "Grey boxes: exposure and outcome. * Unmeasured in MIMIC-IV (dashed boxes and arrows). Every arrow is an edge of the dagitty model.") +
    ggplot2::theme(legend.position = "none", plot.caption = ggplot2::element_text(size = 7, hjust = 0), plot.margin = ggplot2::margin(6, 6, 6, 6))
}
