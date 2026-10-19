# spseg
R package 'spseg' for measuring and visualising population segregation

Pre-1.0 development: the interface may change before the first stable release.

The development interface measures several attributes with shared spatial
neighborhood searches:

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
`counts` selectors; `spseg_grouped()` describes continuous intervals, while
`spseg_continuous(counts = ..., scores = ...)` describes exact observed values.
With `neighbors = "knn"`, each variable selects its own eligible population:
`bands = 100` includes 100 eligible people before distance weighting, even
when the variables cover different populations. Neighbor discovery is shared.

See `vignette("spseg-multiple-variables")` for definitions, sources, deterministic
tie handling, population-count kNN neighborhoods, and exposure cuts. Polygon count specifications also support shared grid and pycnophylactic
surface construction through the same interface.

spseg originated as a fork of the [seg package](https://github.com/syunhong/seg)
by Seong-Yun Hong and David O'Sullivan. It now provides its own interface and
implementation while retaining attribution for reused code, example data, and
other materials.
