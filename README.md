# seg
R package 'seg' for measuring and visualising population segregation

The development interface measures several attributes with shared spatial
neighborhood searches:

```r
spseg(
  people_sf,
  variables = list(
    ethnicity = seg_categorical(values = "ethnicity"),
    education = seg_ordinal(
      values = "education",
      levels = c("primary", "secondary", "tertiary")
    ),
    income = seg_continuous(values = "income"),
    age = seg_continuous(values = "age")
  ),
  measures = "all",
  bands = c(250, 500, 1000),
  scope = "multigroup"
)
```

For a plain data frame, add `coords = c("x", "y")`. For aggregate tables, use
`counts` selectors; `seg_grouped()` describes continuous intervals, while
`seg_continuous(counts = ..., scores = ...)` describes exact observed values.
With `neighbors = "knn"`, each variable selects its own eligible population:
`bands = 100` includes 100 eligible people before distance weighting, even
when the variables cover different populations. Neighbor discovery is shared.

See `vignette("seg-multiple-variables")` for definitions, sources, deterministic
tie handling, population-count kNN neighborhoods, and exposure cuts. The
multiple-variable interface currently supports raw spatial support; existing
categorical grid and pycnophylactic workflows remain available.
