#' Describe a population network
#'
#' Defines shortest-path proximity from a supplied edge list. Population
#' allocation, geographic adjacency, snapping, and road-network preparation
#' take place outside spseg. Use [spseg()] with `network` and `node_id` to
#' attach observations to the graph. Multiple observations may share a node;
#' nodes appearing only in the edge list connect observations without adding
#' population. Observation nodes absent from the edge list are isolated.
#'
#' @param edges A data frame or matrix containing one edge per row.
#' @param from,to Single column names or positions containing endpoint IDs.
#'   IDs must be nonmissing character, factor, or finite numeric values.
#' @param cost Optional single column name or position containing finite,
#'   nonnegative traversal costs. Omit for shortest-path hop counts. These
#'   costs are distinct from the distance kernels used in [spseg()].
#' @param directed Logical. If `FALSE`, each edge can be traversed both ways.
#'   Otherwise traversal follows `from` to `to`, outward from the focal node.
#' @param units Distance-unit label for results, defaulting to `"hops"` or
#'   `"cost"`. It documents the supplied costs without converting them.
#' @return A `spseg_network` specification. Self-loops and parallel edges are
#'   allowed; they do not multiply a node's population.
#' @examples
#' edges <- data.frame(from=c("a","b"), to=c("b","c"), minutes=c(2,3))
#' network <- spseg_network(edges, cost="minutes", units="minutes")
#' nodes <- data.frame(node=c("a","b","c"), group=c("A","B","B"))
#' spseg(nodes, list(group=spseg_categorical(values="group")), network=network,
#'       node_id="node", bands=3, weighting="unweighted")
#' @export
spseg_network <- function(edges, from = "from", to = "to", cost = NULL,
                          directed = FALSE, units = NULL) {
  if (!is.data.frame(edges) && !is.matrix(edges))
    stop("'edges' must be a data frame or matrix", call. = FALSE)
  edges <- as.data.frame(edges)
  endpoint <- function(column) {
    if (length(column) != 1L) stop("Select one edge column", call. = FALSE)
    .spseg_node_ids(.spseg_select(edges, column)[[1L]])
  }
  source <- endpoint(from)
  destination <- endpoint(to)
  weighted <- !is.null(cost)
  if (weighted) {
    if (length(cost) != 1L) stop("Select one cost column", call. = FALSE)
    costs <- .spseg_select(edges, cost)[[1L]]
    if (!is.numeric(costs) || any(!is.finite(costs)) || any(costs < 0))
      stop("Edge costs must be finite nonnegative numbers", call. = FALSE)
  } else costs <- rep(1, nrow(edges))
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed))
    stop("'directed' must be TRUE or FALSE", call. = FALSE)
  if (is.null(units)) units <- if (weighted) "cost" else "hops"
  if (!is.character(units) || length(units) != 1L || is.na(units) || !nzchar(units))
    stop("'units' must be a nonempty distance-unit label", call. = FALSE)
  structure(list(from = source, to = destination, cost = as.double(costs),
    directed = directed, weighted = weighted, units = units), class = "spseg_network")
}

.spseg_node_ids <- function(ids) {
  if (!(is.character(ids) || is.factor(ids) || is.numeric(ids)) ||
      anyNA(ids) || (is.numeric(ids) && any(!is.finite(ids))))
    stop("Node IDs must be nonmissing character or finite numeric values", call. = FALSE)
  ids <- as.character(ids)
  if (any(!nzchar(ids))) stop("Node IDs must not be empty", call. = FALSE)
  ids
}

.spseg_network_input <- function(x, data, coords, network, node_id) {
  if (!inherits(network, "spseg_network"))
    stop("Supply a network created by spseg_network()", call. = FALSE)
  # Revalidate public specifications before passing integer node IDs to C++.
  network <- spseg_network(data.frame(from = network$from, to = network$to,
    cost = network$cost), cost = if (network$weighted) "cost" else NULL,
    directed = network$directed, units = network$units)
  if (length(node_id) != 1L) stop("Select a 'node_id' column for network input", call. = FALSE)
  if (!is.null(coords)) checked <- .spseg_geometry(x, data, coords) else {
    spatial <- inherits(x, "sf")
    checked <- list(coords = NULL,
      data = if (!is.null(data)) as.data.frame(data) else
        if (spatial) sf::st_drop_geometry(x) else as.data.frame(x),
      geometry = if (spatial) sf::st_geometry(x) else NULL,
      crs = if (spatial) sf::st_crs(x) else sf::st_crs(NA_character_))
    if (spatial && nrow(checked$data) != nrow(x))
      stop("Supply one attribute row per feature", call. = FALSE)
  }
  if (nrow(checked$data) < 1L) stop("Supply at least one observation", call. = FALSE)
  ids <- .spseg_node_ids(.spseg_select(checked$data, node_id)[[1L]])
  labels <- unique(c(ids, network$from, network$to))
  graph <- list(from = match(network$from, labels) - 1L,
    to = match(network$to, labels) - 1L, cost = network$cost,
    node = match(ids, labels) - 1L, labels = labels,
    directed = network$directed, weighted = network$weighted)
  list(checked = checked, graph = graph, ids = ids,
       summary = list(directed = network$directed, weighted = network$weighted,
         units = network$units, nodes = length(labels), edges = length(network$from)))
}
