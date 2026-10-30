test_that("categorical column shorthand agrees with explicit specifications", {
  d <- data.frame(x = c(0, 1, 3, 4), y = 0,
    race = c("a", "b", "a", "b"), residents = c(2, 3, 4, 1),
    a = c(2, 0, 4, 0), b = c(0, 3, 0, 1))
  expect_message(
    micro <- spseg(d, "race", coords = c("x", "y"),
      population = "residents", bands = c(0, 2), output = "full"),
    "Assuming categorical observed memberships.*observed category memberships"
  )
  explicit <- spseg(d, list(race = spseg_categorical(vars = "race")),
    coords = c("x", "y"), population = "residents", bands = c(0, 2),
    output = "full")
  expect_identical(names(micro$analyses), "race")
  expect_equal(micro$analyses, explicit$analyses)
  expect_message(
    counts <- spseg(d, c("a", "b"), coords = c("x", "y"),
      bands = c(0, 2), output = "full"),
    "Assuming aggregated categorical data.*category counts"
  )
  explicit_counts <- spseg(d,
    list(categorical = spseg_categorical(vars = c("a", "b"))),
    coords = c("x", "y"), bands = c(0, 2), output = "full")
  expect_identical(names(counts$analyses), "categorical")
  expect_equal(counts$analyses, explicit_counts$analyses)
  expect_equal(micro$analyses$race$indices, counts$analyses$categorical$indices)

  points <- sf::st_as_sf(d, coords = c("x", "y"))
  expect_message(spatial <- spseg(points, c("a", "b"),
    bands = c(0, 2), output = "full"), "Assuming aggregated")
  expect_equal(spatial$analyses, counts$analyses)
  expect_equal(spatial$geometry, sf::st_geometry(points))
  expect_message(separate <- spseg(as.matrix(d[c("x", "y")]), c("a", "b"),
    data = d, bands = c(0, 2), output = "full"), "Assuming aggregated")
  expect_equal(separate$analyses, counts$analyses)
  expect_error(expect_message(
    spseg(d, "absent", coords = c("x", "y"), bands = 0),
    "Assuming categorical"), "Unknown column")
  for (columns in list(character(), NA_character_, "", c("a", "a"))) {
    expect_error(spseg(d, columns, coords = c("x", "y"), bands = 0),
      "non-empty, non-missing, and unique")
  }
  expect_message(spseg(d, list(race = spseg_categorical(vars = "race")),
    coords = c("x", "y"), bands = 0), NA)
})

test_that("vars selects observed columns or ordered count columns", {
  for (constructor in list(spseg_categorical, spseg_ordinal, spseg_continuous)) {
    observed <- constructor("race")
    expect_identical(observed$values, "race")
    expect_null(observed$counts)
    counts <- constructor(c("b", "a"))
    expect_identical(counts$counts, c("b", "a"))
    expect_null(counts$values)
    expect_identical(constructor(2L)$values, 2L)
    expect_identical(constructor(c(3L, 1L))$counts, c(3L, 1L))
    expect_identical(constructor("race", weights = "n")$weights, "n")
    expect_error(constructor(c("a", "b"), weights = "n"),
      "weights.*one column")
    for (invalid in list(NULL, character(), NA_character_, "", c("a", "a"),
                         0, -1, 1.5, Inf, NA_real_, TRUE, list("a"),
                         matrix("a"), 1i)) {
      expect_error(constructor(invalid), "'vars' must select")
    }
    expect_error(constructor(values = "race"), "unused argument")
    expect_error(constructor(counts = c("a", "b")), "unused argument")
  }
})

test_that("weighted long-form categories agree with wide neighborhood counts", {
  wide <- data.frame(x = c(0, 1, 3, 4), y = 0,
    a = c(8, 0, 4, 1), b = c(2, 5, 3, 0), c = c(0, 5, 3, 9))
  groups <- c("a", "b", "c")
  rows <- rep(seq_len(nrow(wide)), each = length(groups))
  long <- data.frame(wide[rows, c("x", "y")],
    race = rep(groups, nrow(wide)),
    n = as.vector(t(as.matrix(wide[groups]))))
  for (weighting in c("unweighted", "biweight")) {
    counts <- spseg(wide, list(race = spseg_categorical(vars = groups)),
      coords = c("x", "y"), bands = c(0, 2, 10),
      weighting = weighting, output = "full")
    members <- spseg(long,
      list(race = spseg_categorical(vars = "race", weights = "n")),
      coords = c("x", "y"), bands = c(0, 2, 10),
      weighting = weighting, output = "full")
    expect_equal(members$analyses$race$indices, counts$analyses$race$indices,
      tolerance = 1e-12)
    for (b in seq_along(counts$bands)) {
      expect_equal(unname(members$analyses$race$env[[b]]),
        unname(counts$analyses$race$env[[b]][rows, , drop = FALSE]),
        tolerance = 1e-12)
    }
  }
})

test_that("the variable interface is required and returns one result class", {
  d <- data.frame(x=1:4,y=0,a=c(1,2,3,4),b=c(4,3,2,1))
  vars <- list(groups=spseg_categorical(vars=c("a","b")))
  expect_error(spseg(d,coords=c("x","y"),bands=1), "variables")
  expect_error(spseg(as.matrix(d[c("x","y")]),as.matrix(d[c("a","b")]),bands=1),
               "named list|named.*list")
  result <- spseg(d,vars,coords=c("x","y"),bands=1,output="full")
  expect_identical(class(result), "spseg_result")
  expect_setequal(names(result$analyses$groups$indices$overall),c("d","h","r","s","p"))
  expect_null(result$indices)
  expect_null(result$env)
  expect_error(as.data.frame(result,shape="long"),"Unused extraction")
  expect_error(sf::st_as_sf(result,shape="wide"),"Unused extraction")
  expect_identical(as.list(result), unclass(result))
})

test_that("three constructors distinguish exact and interval continuous input", {
  exports <- getNamespaceExports("spseg")
  expect_true(all(c("spseg_categorical","spseg_continuous","spseg_ordinal") %in% exports))
  expect_false("spseg_grouped" %in% exports)
  expect_false(exists("spseg_grouped",envir=asNamespace("spseg"),inherits=FALSE))
  exact <- spseg_continuous(vars="income")
  intervals <- spseg_continuous(vars="band",representation="intervals")
  expect_identical(exact$type,"continuous")
  expect_identical(intervals$type,"continuous")
  expect_identical(exact$representation,"exact")
  expect_identical(intervals$representation,"intervals")
  expect_null(exact$interpolation)
  expect_identical(intervals$interpolation,"linear")
})

test_that("interval membership and counts give the same continuous calculation", {
  d <- data.frame(x=c(0,1,3,4),y=0,w=c(2,4,3,1),
    band=ordered(c("low","high","middle","high"),levels=c("low","middle","high")))
  columns <- c("low","middle","high")
  for(label in columns) d[[label]] <- ifelse(d$band==label,d$w,0)
  variables <- list(
    members=spseg_continuous(vars="band",weights="w",representation="intervals",
      summary=list(probs=c(.2,.5))),
    counts=spseg_continuous(vars=columns,representation="intervals",
      summary=list(probs=c(.2,.5))))
  for(mode in c("radius","knn")) for(summary in c("quantiles","shares")) {
    specs <- lapply(variables,function(s) {
      s$summary$type <- summary
      s
    })
    r <- spseg(d,coords=c("x","y"),variables=specs,neighbors=mode,
      bands=if(mode=="radius") c(0,2,10) else c(1,5,10),
      weighting="unweighted",output="full")
    for(part in c("indices","env","normalizers","status","knn"))
      expect_equal(r$analyses$members[[part]],r$analyses$counts[[part]],tolerance=1e-12)
    expect_equal(unname(r$analyses$members$definition$shares),
                 unname(r$analyses$counts$definition$shares))
    expect_identical(r$analyses$members$definition$type,"continuous")
    expect_identical(r$analyses$members$definition$representation,"intervals")
    expect_identical(r$analyses$members$definition$method,"linear_cdf")
    expect_equal(r$analyses$members$normalizers$information,.5)
    expect_equal(r$diagnostics$focal_searches,nrow(d))
    for (what in c("indices", "exposure", "env")) {
      expect_false(any(c("representation", "type", "method") %in%
        names(as.data.frame(r, what = what))))
    }
  }
  d$band <- as.character(d$band)
  expect_error(spseg(d,coords=c("x","y"),variables=list(v=variables$members),bands=1),
               "ordered")
  explicit <- spseg(d,coords=c("x","y"),variables=list(
    v=spseg_continuous(vars="band",weights="w",representation="intervals",levels=columns,
      summary=list(probs=c(.2,.5)))),
    bands=1)
  expect_identical(explicit$analyses$v$definition$representation,"intervals")
})

test_that("exact defaults and explicit representation agree", {
  d <- data.frame(x=1:4,y=0,v=c(0,0,2,5))
  r <- spseg(d,coords=c("x","y"),variables=list(
    default=spseg_continuous(vars="v"),
    exact=spseg_continuous(vars="v",representation="exact"),
    ordinal=spseg_ordinal(vars="v",levels=c(0,2,5))),bands=c(0,2),output="full")
  expect_equal(r$analyses$default,r$analyses$exact)
  table <- as.data.frame(r)
  expect_false(any(c("representation", "type", "method") %in% names(table)))
  expect_identical(r$analyses$exact$definition$representation, "exact")
  expect_null(r$analyses$ordinal$definition$representation)
})

test_that("continuous representation options reject incompatible interpretations", {
  expect_error(spseg_continuous(vars="v",representation="bins"),"exact.*intervals")
  expect_error(spseg_continuous(vars="a",representation="intervals",scores=1),"scores")
  expect_error(spseg_continuous(vars="a",representation="intervals",
                               exposure=list(probs=c(0,1))),"unused argument")
  expect_error(spseg_continuous(vars="a",representation="intervals",
                               interpolation="polynomial"),"linear")
  expect_error(spseg_continuous(vars="v",levels=c("a","b")),"levels")
  expect_error(spseg_continuous(vars="v",interpolation="linear"),"interpolation")
  d <- data.frame(x=1:2,y=0,v=ordered(c("a","b")))
  expect_error(spseg(d,coords=c("x","y"),variables=list(v=spseg_continuous(vars="v")),
                     bands=1),"finite numeric")
})

test_that("all count variable types share grid surfaces and neighborhood discovery", {
  geom <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin=0,ymin=0,xmax=2,ymax=2))),n=c(2,2))
  areas <- sf::st_sf(a=c(8,2,4,6),b=c(2,8,6,4),residents=rep(10,4),
    low=c(2,3,4,5),middle=c(3,1,2,1),high=c(1,2,2,2),
    inc0=c(0,2,3,1),inc1=c(3,2,2,4),
    age1=c(1,2,1,3),age2=c(2,1,2,1),age3=c(2,3,4,2),geometry=geom)
  vars <- list(race=spseg_categorical(vars=c("a","b")),
                education=spseg_ordinal(vars=c("low","middle","high")),
                income=spseg_continuous(representation = "intervals", vars=c("inc0","inc1"),
                  summary=list(probs=6/17)),
                age=spseg_continuous(vars=c("age1","age2","age3"),scores=c(20,20,40)))
  totals <- vapply(vars, function(z) sum(as.matrix(sf::st_drop_geometry(areas)[z$counts])),numeric(1))
  for (surface in c("grid","pycno")) for (neighbors in c("radius","knn")) {
    band <- if(neighbors=="radius") .9 else 3
    result <- spseg(areas,vars,population="residents",bands=band,surface=surface,
                     neighbors=neighbors,nrow=4,ncol=4,converge=1,
                     surface_audit=TRUE,output="full",scope="overall")
    expect_equal(result$diagnostics$focal_searches,nrow(result$coords))
    expect_named(result$surface_info$variables,names(vars))
    for (name in names(vars)) {
      alone <- spseg(areas,vars[name],population="residents",bands=band,surface=surface,
                       neighbors=neighbors,nrow=4,ncol=4,converge=1,
                       output="full",scope="overall")
      expect_equal(result$analyses[[name]],alone$analyses[[name]],tolerance=1e-10)
      expect_equal(result$analyses[[name]]$definition$eligible_population,unname(totals[name]),
                   tolerance=1e-10)
      if(neighbors=="knn") expect_equal(result$analyses[[name]]$knn$min_population,band)
    }
    expect_true(all(vapply(result$analyses,function(z) all(z$status=="ok"),logical(1))))
  }
})
