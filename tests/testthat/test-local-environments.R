test_that("weighted local quantiles use global midranks and preserve ties", {
  d <- data.frame(x=c(0,1,3,4),y=0,v=c(0,0,20,100),w=c(1,3,2,4))
  options <- list(type="quantiles",probs=c(0,.1,.5,.9,1))
  r <- spseg(d,coords=c("x","y"),variables=list(
    v=spseg_continuous(values="v",weights="w",localenv=options)),
    bands=c(1,4),weighting="unweighted",output="full")
  expect_equal(unname(r$analyses$v$env[[1]][1,]),rep(.2,5))
  expect_equal(unname(r$analyses$v$env[[1]][3,]),c(.5,.5,.8,.8,.8))
  expect_equal(unname(r$analyses$v$env[[2]][1,]),c(.2,.2,.5,.8,.8))
  expect_identical(r$analyses$v$definition$localenv$scale,"global_rank")

  counts <- transform(d,low=ifelse(v==0,w,0),mid=ifelse(v==20,w,0),
                       high=ifelse(v==100,w,0))
  aggregate <- spseg(counts,coords=c("x","y"),variables=list(
    v=spseg_continuous(counts=c("low","mid","high"),scores=c(0,20,100),
                       localenv=options)),bands=c(1,4),weighting="unweighted",output="full")
  expect_equal(r$analyses$v$env,aggregate$analyses$v$env,tolerance=1e-12)
  expect_equal(r$analyses$v$indices,aggregate$analyses$v$indices,tolerance=1e-12)
  d$v <- exp(d$v/100)
  transformed <- spseg(d,coords=c("x","y"),variables=list(
    v=spseg_continuous(values="v",weights="w",localenv=options)),
    bands=c(1,4),weighting="unweighted",output="full")
  expect_equal(r$analyses$v$env,transformed$analyses$v$env)
})

test_that("grouped quantiles and shares use the same linear CDF interpolant", {
  d <- data.frame(x=1:3,y=0,a=c(2,0,0),b=c(0,1,1),c=c(0,1,5))
  r <- spseg(d,coords=c("x","y"),variables=list(
    q=spseg_grouped(counts=c("a","b","c"),
      localenv=list(type="quantiles",probs=c(0,.25,.5,.75,1))),
    s=spseg_grouped(counts=c("a","b","c"),
      localenv=list(type="shares",probs=c(0,.1,.5,1)))),
    bands=0,weighting="unweighted",output="full")
  expected <- rbind(c(0,.05,.1,.15,.2),c(.2,.3,.4,.7,1),c(.2,.46,.64,.82,1))
  expect_equal(unname(r$analyses$q$env[[1]]),expected,tolerance=1e-12)
  expect_equal(unname(r$analyses$s$env[[1]]),
               rbind(c(.5,.5,0),c(0,7/12,5/12),c(0,11/36,25/36)),tolerance=1e-12)
  expect_equal(r$analyses$q$indices,r$analyses$s$indices,tolerance=1e-12)
  d$b1 <- d$b*.3; d$b2 <- d$b*.7
  split <- spseg(d,coords=c("x","y"),variables=list(
    q=spseg_grouped(counts=c("a","b1","b2","c"),
      localenv=list(type="quantiles",probs=c(0,.25,.5,.75,1)))),
    bands=0,weighting="unweighted",output="full")
  expect_equal(r$analyses$q$env,split$analyses$q$env,tolerance=1e-12)
})

test_that("local summaries follow the kernel and proportional kNN boundaries", {
  d <- data.frame(x=c(0,1,-1,3),y=0,v=1:4,w=c(1,2,4,3))
  ranks <- (cumsum(d$w)-d$w/2)/sum(d$w)
  probabilities <- c(0,.25,.5,.75,1)
  inverse <- function(weights) {
    present <- which(weights>0)
    cumulative <- cumsum(weights[present])/sum(weights)
    vapply(probabilities,function(p) {
      if(p==1) return(tail(ranks[present],1))
      ranks[present[which(cumulative>=p)[1]]]
    },numeric(1))
  }
  spec <- list(v=spseg_continuous(values="v",weights="w",
    localenv=list(type="quantiles",probs=probabilities)))
  radius <- spseg(d,coords=c("x","y"),variables=spec,bands=4,
                   weighting="exponential",output="full")
  expect_equal(unname(radius$analyses$v$env[[1]][1,]),
               inverse(d$w*exp(-abs(d$x)/4)))
  nearest <- spseg(d,coords=c("x","y"),variables=spec,neighbors="knn",bands=4,
                    weighting="unweighted",output="full")
  expect_equal(unname(nearest$analyses$v$env[[1]][1,]),inverse(c(1,1,2,0)))
  expect_equal(nearest$analyses$v$knn$min_population,4)
})

test_that("presentation settings leave all indices and exposure unchanged", {
  d <- data.frame(x=1:6,y=0,v=c(1,1,2,3,5,8),
    o=ordered(c("a","a","b","c","c","c")))
  specs <- list(
    default=spseg_continuous(values="v"),
    quartiles=spseg_continuous(values="v",
      localenv=list(type="quantiles",probs=c(.25,.5,.75))),
    shares=spseg_continuous(values="v",localenv="shares"),
    custom=spseg_continuous(values="v",
      localenv=list(type="shares",breaks=c(-Inf,2,Inf))),
    ordinal=spseg_ordinal(values="o"),
    categories=spseg_ordinal(values="o",localenv="shares"))
  for(mode in c("radius","knn")) {
    r <- spseg(d,coords=c("x","y"),variables=specs,bands=c(1,3),
                 neighbors=mode,output="full")
    indices <- spseg(d,coords=c("x","y"),variables=specs,bands=c(1,3),neighbors=mode)
    expect_equal(as.data.frame(r),as.data.frame(indices),tolerance=1e-12)
    for(name in c("quartiles","shares","custom"))
      expect_equal(r$analyses$default$indices,r$analyses[[name]]$indices)
    expect_equal(r$analyses$ordinal$indices,r$analyses$categories$indices)
    for(name in c("shares","custom","categories"))
      expect_equal(rowSums(r$analyses[[name]]$env[[2]]),rep(1,nrow(d)),ignore_attr=TRUE)
    expect_equal(ncol(r$analyses$custom$env[[1]]),2)
    expect_equal(ncol(r$analyses$default$env[[1]]),9)
  }
})

test_that("local shares work independently of exposure calculation", {
  d <- data.frame(x=1:10,y=0,v=1:10,o=ordered(rep(letters[1:2],each=5)))
  r <- spseg(d,coords=c("x","y"),variables=list(
    v=spseg_continuous(values="v",localenv="shares",exposure=list(probs=c(0,.5,1))),
    o=spseg_ordinal(values="o",localenv=list(type="shares",probs=c(0,.5,1)))),
    measures="information",bands=20,weighting="unweighted",output="full")
  for(a in r$analyses) {
    expect_equal(unname(a$env[[1]]),matrix(.5,10,2))
    expect_null(a$indices$multigroup$p)
  }
})

test_that("empty local summaries and mixed SF extraction retain focal rows", {
  d <- data.frame(x=c(0,10,20),y=0,g=c("a","b","b"),v=c(NA,1,2))
  x <- sf::st_as_sf(d,coords=c("x","y"),crs=3857,sf_column_name="location")
  r <- spseg(x,variables=list(g=spseg_categorical(values="g"),
    v=spseg_continuous(values="v",na.action="omit")),bands=0,output="localenv")
  expect_true(all(is.na(r$analyses$v$env[[1]][1,])))
  table <- as.data.frame(r,what="env")
  expect_equal(nrow(table),3*(2+9))
  quantiles <- subset(table,summary=="quantiles")
  expect_equal(sort(unique(quantiles$prob)),seq(.1,.9,.1))
  expect_true(all(is.na(quantiles$share)))
  expect_true(all(is.na(quantiles$category)))
  expect_s3_class(sf::st_as_sf(r,variables="v"),"sf")
  expect_equal(nrow(sf::st_as_sf(r,variables="v")),27L)
})

test_that("continuous environment storage stays compact at large bandwidths", {
  d <- data.frame(x=seq(0,1,length.out=500),y=0,v=seq_len(500))
  r <- spseg(d,coords=c("x","y"),variables=list(v=spseg_continuous(values="v")),
    bands=c(1,2,3),weighting="unweighted",output="full")
  expect_equal(dim(r$analyses$v$env[[3]]),c(500L,9L))
  expect_lt(as.numeric(object.size(r)),1024^2)
  # The regional reference is stored once, not repeated per neighborhood.
  expect_equal(r$analyses$v$definition$scores,as.numeric(d$v))
})

test_that("local summary options reject ambiguous or invalid specifications", {
  for(p in list(numeric(),c(.5,.5),c(.8,.2),c(-.1,.5),NA_real_,"median"))
    expect_error(spseg_continuous(values="v",localenv=list(type="quantiles",probs=p)),
                 "probabilities")
  expect_error(spseg_ordinal(values="o",localenv="unknown"),"localenv")
  expect_error(spseg_continuous(values="v",localenv=list(type="quantiles",breaks=0:1)),
               "not 'breaks'")
  expect_error(spseg_grouped(counts="a",localenv=list(type="shares",probs=c(.1,.9))),
               "from 0 to 1")
  expect_error(spseg_continuous(values="v",localenv=list(type="shares",probs=0:1,breaks=0:1)),
               "only one")
  expect_error(spseg_continuous(values="v",localenv=list(type="shares",breaks=c(2,1))),
               "breaks")
  d <- data.frame(x=1:2,y=0,a=1,b=2)
  expect_error(spseg(d,coords=c("x","y"),variables=list(g=spseg_grouped(counts=c("a","b"),
    localenv=list(type="shares",breaks=0:1))),bands=0),"exact continuous")
})
