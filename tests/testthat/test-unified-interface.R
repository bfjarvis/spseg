test_that("the variable interface is required and returns one result class", {
  d <- data.frame(x=1:4,y=0,a=c(1,2,3,4),b=c(4,3,2,1))
  vars <- list(groups=spseg_categorical(counts=c("a","b")))
  expect_error(spseg(d,coords=c("x","y"),bands=1), "variables")
  expect_error(spseg(as.matrix(d[c("x","y")]),as.matrix(d[c("a","b")]),bands=1),
               "named list|named.*list")
  result <- spseg(d,vars,coords=c("x","y"),bands=1,output="full")
  expect_identical(class(result), "spseg_result")
  expect_setequal(names(result$analyses$groups$indices$multigroup),c("d","h","r","s","p"))
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
  exact <- spseg_continuous(values="income")
  intervals <- spseg_continuous(values="band",representation="intervals")
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
    members=spseg_continuous(values="band",weights="w",representation="intervals",
      summary=list(probs=c(.2,.5))),
    counts=spseg_continuous(counts=columns,representation="intervals",
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
    expect_true(all(as.data.frame(r)$representation=="intervals"))
    expect_true(all(as.data.frame(r,what="exposure")$representation=="intervals"))
    expect_true(all(as.data.frame(r,what="env")$representation=="intervals"))
  }
  d$band <- as.character(d$band)
  expect_error(spseg(d,coords=c("x","y"),variables=list(v=variables$members),bands=1),
               "ordered")
  explicit <- spseg(d,coords=c("x","y"),variables=list(
    v=spseg_continuous(values="band",weights="w",representation="intervals",levels=columns,
      summary=list(probs=c(.2,.5)))),
    bands=1)
  expect_identical(explicit$analyses$v$definition$representation,"intervals")
})

test_that("exact defaults and explicit representation agree", {
  d <- data.frame(x=1:4,y=0,v=c(0,0,2,5))
  r <- spseg(d,coords=c("x","y"),variables=list(
    default=spseg_continuous(values="v"),
    exact=spseg_continuous(values="v",representation="exact"),
    ordinal=spseg_ordinal(values="v",levels=c(0,2,5))),bands=c(0,2),output="full")
  expect_equal(r$analyses$default,r$analyses$exact)
  table <- as.data.frame(r)
  expect_true(all(table$representation[table$type=="continuous"]=="exact"))
  expect_true(all(is.na(table$representation[table$type=="ordinal"])))
})

test_that("continuous representation options reject incompatible interpretations", {
  expect_error(spseg_continuous(values="v",representation="bins"),"exact.*intervals")
  expect_error(spseg_continuous(counts="a",representation="intervals",scores=1),"scores")
  expect_error(spseg_continuous(counts="a",representation="intervals",
                               exposure=list(probs=c(0,1))),"Interval exposure")
  expect_error(spseg_continuous(counts="a",representation="intervals",
                               interpolation="polynomial"),"linear")
  expect_error(spseg_continuous(values="v",levels=c("a","b")),"levels")
  expect_error(spseg_continuous(values="v",interpolation="linear"),"interpolation")
  d <- data.frame(x=1:2,y=0,v=ordered(c("a","b")))
  expect_error(spseg(d,coords=c("x","y"),variables=list(v=spseg_continuous(values="v")),
                     bands=1),"finite numeric")
})

test_that("all count variable types share grid surfaces and neighborhood discovery", {
  geom <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin=0,ymin=0,xmax=2,ymax=2))),n=c(2,2))
  areas <- sf::st_sf(a=c(8,2,4,6),b=c(2,8,6,4),residents=rep(10,4),
    low=c(2,3,4,5),middle=c(3,1,2,1),high=c(1,2,2,2),
    inc0=c(0,2,3,1),inc1=c(3,2,2,4),
    age1=c(1,2,1,3),age2=c(2,1,2,1),age3=c(2,3,4,2),geometry=geom)
  vars <- list(race=spseg_categorical(counts=c("a","b")),
                education=spseg_ordinal(counts=c("low","middle","high")),
                income=spseg_continuous(representation = "intervals", counts=c("inc0","inc1"),
                  summary=list(probs=6/17)),
                age=spseg_continuous(counts=c("age1","age2","age3"),scores=c(20,20,40)))
  totals <- vapply(vars, function(z) sum(as.matrix(sf::st_drop_geometry(areas)[z$counts])),numeric(1))
  for (surface in c("grid","pycno")) for (neighbors in c("radius","knn")) {
    band <- if(neighbors=="radius") .9 else 3
    result <- spseg(areas,vars,population="residents",bands=band,surface=surface,
                     neighbors=neighbors,nrow=4,ncol=4,converge=1,
                     surface_audit=TRUE,output="full",scope="multigroup")
    expect_equal(result$diagnostics$focal_searches,nrow(result$coords))
    expect_named(result$surface_info$variables,names(vars))
    for (name in names(vars)) {
      alone <- spseg(areas,vars[name],population="residents",bands=band,surface=surface,
                       neighbors=neighbors,nrow=4,ncol=4,converge=1,
                       output="full",scope="multigroup")
      expect_equal(result$analyses[[name]],alone$analyses[[name]],tolerance=1e-10)
      expect_equal(result$analyses[[name]]$definition$eligible_population,unname(totals[name]),
                   tolerance=1e-10)
      if(neighbors=="knn") expect_equal(result$analyses[[name]]$knn$min_population,band)
    }
    expect_true(all(vapply(result$analyses,function(z) all(z$status=="ok"),logical(1))))
  }
})
