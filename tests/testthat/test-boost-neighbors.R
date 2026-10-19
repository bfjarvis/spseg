test_that("radius searches include closed boundaries and coincident points", {
  d <- data.frame(x=c(0,9,0,10),y=c(0,4,0,0),a=c(2,0,0,0),b=c(0,2,2,2))
  vars <- list(g=spseg_categorical(counts=c("a","b")))
  for (band in c(0,sqrt(97))) {
    indexed <- spseg(d,vars,coords=c("x","y"),bands=band,
                     weighting="unweighted",output="full")
    brute <- spseg(d,vars,coords=c("x","y"),bands=band,
                   weighting="unweighted",output="full",search="brute")
    expect_equal(indexed$analyses,brute$analyses,tolerance=1e-12)
    expect_equal(unname(indexed$analyses$g$env[[1]][1,]),
                 if(band==0)c(.5,.5) else c(1/3,2/3))
  }
})

test_that("KNN includes tied boundary populations proportionately", {
  d <- data.frame(x=c(0,-1,1,0),y=c(0,0,0,1),a=c(10,20,0,0),b=c(0,0,10,30))
  vars <- list(g=spseg_categorical(counts=c("a","b")))
  run <- function(d,search="boost") spseg(d,vars,coords=c("x","y"),
    neighbors="knn",bands=25,weighting="unweighted",output="full",search=search)
  indexed <- run(d)
  # Ten focal people, then 15/60 of each population at distance one.
  expect_equal(unname(indexed$analyses$g$env[[1]][1,]),c(.6,.4))
  expect_equal(indexed$analyses,run(d,"brute")$analyses,tolerance=1e-12)
  permutation <- c(1,4,2,3)
  reordered <- run(d[permutation,])
  expect_equal(reordered$analyses$g$env[[1]],
               indexed$analyses$g$env[[1]][permutation,],tolerance=1e-12)
  expect_equal(reordered$analyses$g$indices,indexed$analyses$g$indices,tolerance=1e-12)
})

test_that("coincident observations share a fractional boundary after the focal", {
  d <- data.frame(x=0,y=0,g=c("a","b","b","b"))
  r <- spseg(d,list(g=spseg_categorical(values="g")),coords=c("x","y"),
             neighbors="knn",bands=2,weighting="unweighted",output="full")
  expect_equal(unname(r$analyses$g$env[[1]][1,]),c(.5,.5))
  expect_equal(unname(r$analyses$g$env[[1]][2,]),c(1/6,5/6))
  expect_equal(r$analyses$g$knn$min_population,2)
  expect_equal(r$analyses$g$knn$max_population,2)
})
