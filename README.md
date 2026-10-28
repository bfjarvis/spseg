# spseg

<img src="man/figures/logo.png" align="right" width="240" alt="spseg hex logo: yellow dot-pattern lettering on a viridis population surface" />

**spseg** is an R package for measuring segregation in individual records,
small-area population counts, and networks. It implements categorical, ordinal,
and continuous segregation measures based on the work of Reardon and colleagues,
with a common interface for geographic and network neighborhoods.

Use it to study how segregation changes with neighborhood scale, describe the
population surrounding each location, and calculate groups' exposure to one
another. Examples include segregation by ethnicity, education, and income.

spseg provides:

- Dissimilarity, information theory, relative diversity, and square-root indices,
  plus exposure and isolation summaries.
- Fixed-distance and population nearest-neighbor neighborhoods, with spatial
  weighting options.
- Analysis of several variables in one call, sharing neighborhood searches.
- Local population compositions and income-rank summaries, with tabular and
  spatial outputs for plotting.
- Polygon population surfaces and neighborhoods defined by network hops or
  travel costs.

<br clear="right" />

## Installation

Install the development version from GitHub:

```r
install.packages("pak")
pak::pak("bfjarvis/spseg")
```

spseg requires R 4.2 or later. The GitHub installation compiles C++ code and
installs dependencies, including `sf` and `terra`. When updating an installation,
restart R before loading the updated package.

The package is in pre-1.0 development. Once a CRAN release is available, its
installation command will be:

```r
install.packages("spseg")
```

The examples below use tidyverse for data preparation and extraction:

```r
install.packages("tidyverse")
```

## A first analysis

Each row in this small example represents a person, with coordinates, a group,
an ordered education category, and an income. The coordinates are synthetic
planar units. For real geographic data, use a suitable projected coordinate system
so distances have a meaningful unit, such as meters.

```r
library(spseg)
library(tidyverse)

people <- tibble(
  x = c(0, 0, 1, 1, 2, 2, 3, 3),
  y = 0,
  group = c("A", "A", "A", "B", "B", "B", "A", "B"),
  education = ordered(
    c("basic", "basic", "middle", "middle", "higher", "higher", "middle", "higher"),
    levels = c("basic", "middle", "higher")
  ),
  income = c(10, 15, 20, 30, 40, 60, 80, 100)
)

variables <- list(
  group = spseg_categorical(values = "group"),
  education = spseg_ordinal(values = "education"),
  income = spseg_continuous(values = "income")
)

result <- spseg(
  people,
  coords = c("x", "y"),
  variables = variables,
  neighbors = "radius",
  bands = c(1, 2),
  weighting = "biweight",
  measures = "all",
  output = "full"
)

result
```

`variables` names the analyses and describes their inputs. Each `values` argument
selects one column. `bands` sets the neighborhood radii; the biweight kernel gives
nearer residents more weight. `measures = "all"` requests the scalar indices and
exposure summaries. `output = "full"` also retains local environments for inspection
and mapping.

For an `sf` object, pass it directly as the first argument; its active geometry
supplies the locations. Points use their coordinates and raw polygons use their
centroids.

### Explore the results

Extract segregation indices as a tidy table:

```r
indices <- result |>
  as.data.frame(what = "indices") |>
  as_tibble()

indices |>
  select(variable, band, measure, value)
```

Compare values across scales within each variable and measure. The categorical,
ordinal, and continuous indices describe their respective kinds of variation.

Extract exposure to each group, averaged over members of each focal group:

```r
result |>
  as.data.frame(what = "exposure", variables = "group", bands = 2) |>
  select(from_group, to_group, value)
```

The diagonal entries are isolation: a group's average exposure to itself.
Result tables use one numeric `value` column. Calculation metadata remains in
`result$analyses[[name]]$definition`; inspect status and kNN information separately
with `as.data.frame(result, what = "diagnostics")`.
For local environments, extract each location's neighborhood composition:

```r
result |>
  as.data.frame(what = "env", variables = "group", bands = 2) |>
  as_tibble()
```

Use `sf::st_as_sf(result, variables = "group", bands = 2)` to obtain the same
local summaries with geometry, ready for mapping. `what = "pairwise"` extracts
binary segregation indices between categorical groups or resolved income groups.

For larger analyses, the default `output = "indices"` retains indices and grouped
exposures. Use `output = "localenv"` when the local environments are the desired
result.

### Choose a neighborhood

Use `neighbors = "radius"` for physical distance bands. For adaptive neighborhoods,
use `neighbors = "knn"` and give population thresholds:

```r
nearest <- spseg(
  people,
  coords = c("x", "y"),
  variables = variables,
  neighbors = "knn",
  bands = c(4, 8),
  weighting = "unweighted"
)
```

These thresholds count eligible people, including the focal population. With
aggregate count data, they count residents represented by the rows. Each variable
uses its own eligible population.

## Income and other continuous variables

Continuous indices use the underlying income distribution. The `summary` option
chooses how to describe local environments and groups' experiences. The default
reports shares of four global quartile groups. For a different grouping, supply
interior probabilities:

```r
spseg_continuous(values = "income", summary = list(probs = c(.2, .4, .6, .8)))
```

To describe the income ranks found at local quantiles, select quantile summaries:

```r
spseg_continuous(
  values = "income",
  summary = list(type = "quantiles", probs = c(.25, .5, .75), mean_rank = TRUE)
)
```

This reports global income ranks at each neighborhood's quartiles, plus its mean
rank. The same probabilities define the focal groups used to average those local
summaries. Extract the grouped profiles with `what = "exposure"`; the `statistic`,
`prob`, and `value` columns identify their meaning. These are averages of local
rank statistics, not exposure probabilities of encountering a target group.

For absolute income thresholds, use
`summary = list(type = "shares", cuts = c(30000, 60000))`.

## Aggregated population counts

For small-area tables, select category-count columns with `counts`. This example
has two areas and three ordered income intervals:

```r
areas <- tibble(
  x = c(0, 1), y = 0,
  group_A = c(40, 10), group_B = c(10, 40),
  income_low = c(30, 10), income_middle = c(15, 20), income_high = c(5, 20)
)

area_result <- spseg(
  areas,
  coords = c("x", "y"),
  variables = list(
    group = spseg_categorical(counts = c("group_A", "group_B")),
    income = spseg_continuous(
      counts = c("income_low", "income_middle", "income_high"),
      representation = "intervals"
    )
  ),
  bands = 0,
  weighting = "unweighted"
)
```

With these distinct area locations, `bands = 0` treats each area's own residents
as its local environment, giving non-spatial segregation across areas. Positive
distance bands combine nearby areas into spatial neighborhoods. Ordered-category
counts use `spseg_ordinal(counts = c(...))` in the same way.

For interval income, column order defines the income ordering. spseg interpolates
in global-rank space and proportionally allocates interval populations to the
requested summary groups. A warning identifies when groups split supplied
intervals. The model assumes the same within-interval rank distribution at each
location; allocation fractions and local diagnostics are available in the result.
Extract the neighborhood diagnostics with `as.data.frame(result, what = "interpolation")`.

For polygon `sf` count data, `surface = "grid"` distributes population across a
fine grid before measuring neighborhoods. `surface = "pycno"` adds
pycnophylactic smoothing, preserving source-zone category totals through
rescaling. See the [surface guide](vignettes/spseg-surface-construction.Rmd) for
resolution choices and count-conservation audits.

## Network neighborhoods

Supply an edge list to define proximity by shortest-path hops or travel costs.
Map each population record to a node with `node_id`. Here the people from the
first example occupy four nodes connected by travel-time edges:

```r
edges <- tibble(
  from = c("0", "1", "2"),
  to = c("1", "2", "3"),
  minutes = c(2, 3, 1)
)

network <- spseg_network(edges, cost = "minutes", units = "minutes")
nodes <- people |> mutate(node = as.character(x))

network_result <- spseg(
  nodes,
  variables = variables,
  network = network,
  node_id = "node",
  bands = c(2, 5),
  weighting = "unweighted"
)
```

The bands here are minutes of shortest-path travel. Omit `cost` for hop-count
distance, or use `directed = TRUE` for outward travel along directed edges.
The [network guide](vignettes/spseg-networks.Rmd) includes adjacency and road-network
preparation with `sfnetworks`.

## Learn more

- [Getting started](vignettes/spseg-getting-started.Rmd): a complete workflow from
  simulated people to polygon counts and parallel comparisons across patterns.
- [Multiple variables and methods](vignettes/spseg-multiple-variables.Rmd): measure
  definitions, sources, input choices, and interpretation.
- Installed help: `?spseg`, `?spseg_variables`, and `?spseg_result`.
- [GitHub Issues](https://github.com/bfjarvis/spseg/issues): bug reports, feature
  requests, and questions about methods or workflows.

spseg builds on methods developed by Reardon and colleagues and originated as a
fork of [seg](https://github.com/syunhong/seg) by Seong-Yun Hong and David
O'Sullivan. Attribution for reused code, example data, and other materials is
retained in the package.
