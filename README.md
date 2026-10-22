# spseg
R package 'spseg' for measuring and visualising population segregation

Pre-1.0 development: the interface may change before the first stable release.

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

For neighborhood construction, use `output = "localenv"` (or `"full"` to also
calculate indices). Categorical and ordinal variables retain category shares by
default. Continuous variables retain local deciles of global ranks for either
input representation. For quartiles, specify
`localenv = list(type = "quantiles", probs = c(.25, .5, .75))` in
`spseg_ordinal()` or `spseg_continuous()` (either representation). Use
`localenv = "shares"` for local category shares. These presentation settings
do not coarsen the index calculations. Extract summaries with
`as.data.frame(result, what = "env")` or `sf::st_as_sf(result)`.

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
