# spseg
R package 'spseg' for measuring and visualising population segregation

Pre-1.0 development: the interface may change before the first stable release.

## Installing Local Changes

Install from the repository root in a fresh R session, before loading spseg or
opening its installed help:

```r
pak::pak(".")
```

If spseg was loaded in an analysis session during installation, restart that
session before continuing. In RStudio, use **Session > Restart R**. Calling
`library(spseg)` again does not reload the existing namespace or compiled code.
The restart is also needed when installation ran in a separate terminal or R
process while the analysis session still had spseg loaded.

Reinstalling a loaded package can leave the session's lazy-load state or cached
help inconsistent with the replacement files. A later help request or function
call can then report `internal error 1 in R_decompress1 with libdeflate`, even
when the new installation is healthy. Restart R and retry before reinstalling
packages or changing compression settings. If the warning persists in a fresh
session, record the warning and the installation path from `find.package("spseg")`
so that the affected installation can be investigated.

Start with `vignette("spseg-getting-started")` for a reproducible synthetic-population
example covering individual data, regular and Voronoi aggregation, and
count-preserving population surfaces. The older `spseg-overview` vignette is
retained as a legacy guide to the bundled population patterns.

The development interface measures several attributes with shared spatial
neighborhood searches. Euclidean distance is the default:

```r
library(spseg)

spseg(
  people_sf,
  variables = list(
    ethnicity = spseg_categorical(values = "ethnicity"),
    education = spseg_ordinal(
      values = "education",
      levels = c("primary", "secondary", "tertiary")
    ),
    income = spseg_continuous(values = "income"),
    age = spseg_continuous(values = "age")
  ),
  measures = "all",
  bands = c(250, 500, 1000),
  scope = "multigroup"
)
```

For a plain data frame, add `coords = c("x", "y")`. For aggregate tables, use
`counts` selectors. The three variable constructors are `spseg_categorical()`,
`spseg_ordinal()`, and `spseg_continuous()`. Continuous input defaults to
`representation = "exact"`; select `representation = "intervals"` for interval
membership or ordered interval counts. Exact-value count tables also require
`scores`; interval input uses CDF interpolation without imputing individual values.
With `neighbors = "knn"`, each variable selects its own eligible population:
`bands = 100` includes 100 eligible people before distance weighting, even
when the variables cover different populations. Neighbor discovery is shared.

Continuous variables default to synchronized shares of four global income
groups, using interior quartile probabilities `.25`, `.5`, and `.75`. Set
`summary = list(type = "quantiles")` for local quantiles of global income
ranks and their focal-group averages; optionally add `mean_rank = TRUE` for
the distinct mean-rank statistic. The same probability grid defines global
focal groups and local quantile positions. Shares and quantile summaries can
be extracted with `as.data.frame(result, what = "exposure")`; inspect its
`summary`, `statistic`, `from_group`, `prob`, `exposure`, and `rank` columns
to distinguish share rows from rank profiles. Local summaries are retained
with `output = "full"` or `"localenv"` and extracted with
`as.data.frame(result, what = "env")` or `sf::st_as_sf(result)`.

For exact income, `summary = list(type = "shares", cuts = c(...))` uses
absolute value thresholds instead of probabilities. For interval input,
the requested probabilities are retained. Categories crossing a requested boundary
are proportionally split in global-percentile space, with the same allocation
applied to focal populations and local target shares. This assumes locations share
the same within-category rank distribution; no monetary maximum is needed.
A brief warning flags interpolation. Retained neighborhood output includes
`bottom_category_share`, `top_category_share`, and `interpolated_share` (the local
share in split categories). Allocation fractions are recorded in
`result$analyses$income$definition$summary$groups$allocation`. Probabilities aligned
with actual interval boundaries avoid apportionment. Interval labels alone do not
identify absolute income cuts. Changing summaries does not change the
full-resolution continuous scalar indices. The explicit legacy continuous
`localenv` and `exposure` arguments remain temporarily supported with a
deprecation warning; do not combine them with `summary`.

Use `scope = "pairwise"` to obtain binary segregation indices between resolved
income groups. These group-pair indices use the original variable neighborhoods
and do not replace continuous scalar indices.

See `vignette("spseg-multiple-variables")` for definitions, sources, deterministic
tie handling, population-count kNN neighborhoods, and exposure cuts. Polygon count specifications also support shared grid and pycnophylactic
surface construction through the same interface.

For network neighborhoods, supply an edge list and map observation rows to nodes:

```r
network <- spseg_network(edges, from="source", to="destination",
                         cost="minutes", units="minutes", directed=FALSE)
spseg(nodes, variables, network=network, node_id="node",
      bands=c(5,10), weighting="unweighted")
```

Omit `cost` for hop distances. Geometry is optional and never determines network
relationships. Users prepare adjacency, road networks, and population allocation
outside spseg. See `vignette("spseg-networks")` for examples, including preparation
with sfnetworks. Boost.Geometry and Boost.Graph are supplied through BH; the
measurement engine remains single-threaded.

spseg originated as a fork of the [seg package](https://github.com/syunhong/seg)
by Seong-Yun Hong and David O'Sullivan. It now provides its own interface and
implementation while retaining attribution for reused code, example data, and
other materials.
