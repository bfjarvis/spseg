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

test_that("all count variable types share grid surfaces and neighborhood discovery", {
  geom <- sf::st_make_grid(sf::st_as_sfc(sf::st_bbox(
    c(xmin=0,ymin=0,xmax=2,ymax=2))),n=c(2,2))
  areas <- sf::st_sf(a=c(8,2,4,6),b=c(2,8,6,4),residents=rep(10,4),
    low=c(2,3,4,5),middle=c(3,1,2,1),high=c(1,2,2,2),
    inc0=c(0,2,3,1),inc1=c(3,2,2,4),
    age1=c(1,2,1,3),age2=c(2,1,2,1),age3=c(2,3,4,2),geometry=geom)
  vars <- list(race=spseg_categorical(counts=c("a","b")),
                education=spseg_ordinal(counts=c("low","middle","high")),
                income=spseg_grouped(counts=c("inc0","inc1")),
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
