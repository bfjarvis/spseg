ordered_reference <- function(counts, local, type, scores = seq_len(ncol(counts))) {
  mass <- rowSums(counts)
  alpha <- mass / sum(mass)
  global <- colSums(counts) / sum(counts)
  q <- cumsum(global)
  q[length(q)] <- 1
  z <- cumsum(global) - global / 2
  p <- local / rowSums(local)
  functions <- list(h = function(x) {
    y <- x*0; ok <- x > 0 & x < 1
    y[ok] <- -x[ok]*log(x[ok])-(1-x[ok])*log1p(-x[ok]); y
  }, r = function(x) x*(1-x), s = function(x) 2*sqrt(pmax(x*(1-x),0)))
  if (type == "ordinal") {
    cq <- q[-length(q)]
    cp <- t(apply(p, 1, cumsum))[, -length(q), drop = FALSE]
    out <- vapply(functions, function(f) 1-sum(alpha*rowSums(f(cp)))/sum(f(cq)), numeric(1))
    out["d"] <- sum(alpha*rowSums(abs(sweep(cp, 2, cq)))) / (2*sum(cq*(1-cq)))
    return(out)
  }
  if (type == "continuous") {
    width <- diff(z)
    cp <- t(apply(p, 1, cumsum))[, -length(q), drop = FALSE]
    cq <- q[-length(q)]
    out <- vapply(functions, function(f) {
      1-sum(alpha*as.vector(f(cp) %*% width))/sum(width*f(cq))
    }, numeric(1))
    out["d"] <- sum(alpha*as.vector(abs(sweep(cp,2,cq)) %*% width)) /
      (2*sum(width*cq*(1-cq)))
    return(out)
  }
  knots <- c(0, q)
  # Numeric integration supplies an independent reference for the analytical
  # linear-CDF formula, including nonuniform global bin shares.
  out <- vapply(names(functions), function(key) {
    f <- functions[[key]]
    values <- vapply(seq_len(nrow(p)), function(i) {
      cumulative <- c(0, cumsum(p[i, ]))
      cumulative[length(cumulative)] <- 1
      sum(vapply(seq_len(ncol(p)), function(k) {
        g <- function(x) cumulative[k] + (x-knots[k]) / diff(knots[k:(k+1)]) *
          diff(cumulative[k:(k+1)])
        integrate(function(x) f(g(x)),knots[k],knots[k+1],rel.tol=1e-10)$value
      }, numeric(1)))
    }, numeric(1))
    1-sum(alpha*values)/integrate(f,0,1,rel.tol=1e-10)$value
  }, numeric(1))
  distance <- vapply(seq_len(nrow(p)), function(i) {
    cumulative <- c(0,cumsum(p[i,]))
    cumulative[length(cumulative)] <- 1
    sum(vapply(seq_len(ncol(p)), function(k) {
      g <- function(x) cumulative[k] + (x-knots[k])/diff(knots[k:(k+1)]) *
        diff(cumulative[k:(k+1)])
      a <- cumulative[k]-knots[k]; b <- cumulative[k+1]-knots[k+1]
      breaks <- knots[k:(k+1)]
      if (a*b<0) breaks <- sort(c(breaks,knots[k]+diff(breaks)*abs(a)/(abs(a)+abs(b))))
      sum(vapply(seq_len(length(breaks)-1), function(j)
        integrate(function(x) abs(g(x)-x),breaks[j],breaks[j+1],
                  abs.tol=1e-10,rel.tol=1e-9)$value,numeric(1)))
    },numeric(1)))
  },numeric(1))
  out["d"] <- 3*sum(alpha*distance)
  out
}

test_that("categorical measures match independent formulas across spatial options", {
  set.seed(103)
  d <- data.frame(x=runif(18),y=runif(18),a=rpois(18,4),b=rpois(18,6),c=rpois(18,3))
  variables <- list(first=seg_categorical(counts=c("a","b","c")),
                    second=seg_categorical(counts=c("a","b","c")))
  for (neighborhood in c("radius","knn")) for (kernel in c("unweighted","biweight","inverse","exponential")) {
    bands <- if (neighborhood=="radius") c(0,.35,.8) else c(10,40,100)
    result <- spseg(d, variables, coords=c("x","y"), bands=bands,
                    neighbors=neighborhood, weighting=kernel, output="full")
    for (band in seq_along(bands)) {
      reference <- categorical_reference(d[c("a","b","c")],
                                           result$analyses$first$env[[band]])
      for (key in c("d","h","r","s")) {
        expect_lte(abs(unname(result$analyses$first$indices$multigroup[[key]][band])-
                         reference$multigroup[[key]]), 1e-12)
        expect_equal(result$analyses$first$indices$pairwise[[key]][[band]],
                     reference$pairwise[[key]][[1]], tolerance=1e-12)
      }
      expect_equal(result$analyses$first$indices$multigroup$p[[band]],
                   reference$multigroup$p[[1]], tolerance=1e-12)
    }
    expect_equal(result$analyses$first,result$analyses$second)
    expect_equal(result$diagnostics$focal_searches,nrow(d))
  }
})

test_that("ordinal, grouped, and empirical calculators match independent references", {
  set.seed(67)
  counts <- matrix(sample(1:20,24,TRUE),6,dimnames=list(NULL,c("a","b","c","d")))
  d <- data.frame(x=1:6,y=0,counts)
  distance <- abs(outer(d$x,d$x,"-"))
  w <- (1-(distance/2.5)^3)^3; w[distance>2.5] <- 0
  local <- w %*% counts
  variables <- list(o=seg_ordinal(counts=colnames(counts)),
                    g=seg_grouped(counts=colnames(counts)),
                    c=seg_continuous(counts=colnames(counts),scores=c(-5,0,1,100)))
  result <- spseg(d,coords=c("x","y"),variables=variables,bands=2.5)
  for (v in names(variables)) {
    reference <- ordered_reference(counts,local,c(o="ordinal",g="grouped",c="continuous")[[v]])
    for (key in names(reference)) expect_equal(
      unname(result$analyses[[v]]$indices$multigroup[[key]]),unname(reference[key]),tolerance=2e-8)
  }
})

test_that("raw tied and weighted observations equal exact-value count input", {
  d <- data.frame(x=c(0,1,3,4),y=0,v=c(-3,-3,2,9),w=c(2,5,3,4))
  expanded <- d[rep(seq_len(nrow(d)),d$w),]
  counts <- cbind(a=ifelse(d$v==-3,d$w,0),b=ifelse(d$v==2,d$w,0),c=ifelse(d$v==9,d$w,0))
  a <- spseg(d,coords=c("x","y"),population="w",
             variables=list(v=seg_continuous(values="v")),bands=c(0,2,10),output="full")
  b <- spseg(expanded,coords=c("x","y"),variables=list(v=seg_continuous(values="v")),
             bands=c(0,2,10))
  c <- spseg(cbind(d,counts),coords=c("x","y"),population="w",
             variables=list(v=seg_continuous(counts=c("a","b","c"),scores=c(-3, -3+5,9))),
             bands=c(0,2,10))
  expect_equal(a$analyses$v$indices,b$analyses$v$indices,tolerance=1e-12)
  expect_equal(a$analyses$v$indices,c$analyses$v$indices,tolerance=1e-12)
  expect_equal(a$analyses$v$definition$exposure$shares,c(7,3,4)/14)
  expect_equal(sum(a$analyses$v$env[[1]][[1]]$share),1)
  expect_lte(length(a$analyses$v$env[[1]][[1]]$share),length(unique(d$v)))
  permutation <- c(4,2,1,3)
  changed <- d[permutation,]; changed$v <- exp(changed$v)
  r <- spseg(changed,coords=c("x","y"),population="w",
             variables=list(v=seg_continuous(values="v")),bands=c(0,2,10))
  expect_equal(a$analyses$v$indices,r$analyses$v$indices,tolerance=1e-12,
               ignore_attr=TRUE)
})

test_that("categorical root is Hutchens for two groups and Atkinson for several", {
  counts <- rbind(c(50,50,0),c(0,50,50),c(50,0,50))
  colnames(counts) <- c("a","b","c")
  d <- data.frame(x=1:3,y=0,counts)
  result <- spseg(d,coords=c("x","y"),variables=list(g=seg_categorical(counts=c("a","b","c"))),
                  bands=0,measures="square_root")
  expect_equal(unname(result$analyses$g$indices$multigroup$s),1)
  expect_equal(result$analyses$g$indices$pairwise$s[[1]][upper.tri(counts)],rep(.5,3))
  binary <- spseg(as.matrix(d[c("x","y")]), variables = list(groups = seg_categorical(counts = seq_len(ncol(counts[,1:2])))), data = counts[,1:2],bands=c(0,2),
                  measures=c("square_root","information"),output="full")
  totals <- colSums(counts[,1:2])
  direct <- 1-sum(sqrt(counts[,1]/totals[1]*counts[,2]/totals[2]))
  expect_equal(unname(binary$analyses$groups$indices$multigroup$s[1]),direct)
  expect_equal(binary$analyses$groups$indices$multigroup$s[1],binary$analyses$groups$indices$pairwise$s[[1]][1,2],ignore_attr=TRUE)
  expect_equal(binary$analyses$groups$env[[1]],counts[,1:2]/rowSums(counts[,1:2]),ignore_attr=TRUE)
  expect_true("square_root" %in% as.data.frame(binary)$measure)
  recomputed <- categorical_reference(counts[,1:2],binary$analyses$groups$env[[1]])
  expect_equal(recomputed$multigroup$s,unname(binary$analyses$groups$indices$multigroup$s[1]))
  expect_equal(recomputed$pairwise$s[[1]],binary$analyses$groups$indices$pairwise$s[[1]])
})

test_that("two-category ordinal measures reduce to binary categorical measures", {
  d <- data.frame(x=1:3,y=0,a=c(80,60,10),b=c(20,40,90))
  result <- spseg(d,coords=c("x","y"),variables=list(
    cat=seg_categorical(counts=c("a","b")),ord=seg_ordinal(counts=c("a","b"))),bands=0)
  expect_equal(result$analyses$cat$indices$multigroup,result$analyses$ord$indices$multigroup,
               tolerance=1e-12)
  expect_equal(unname(result$analyses$ord$indices$multigroup$d),8/15)
})

test_that("proportional subdivision of a continuous interval preserves the interpolant", {
  d <- data.frame(x=1:4,y=0,a=c(1,4,6,2),b=c(2,6,1,3),c=c(8,2,3,6))
  d$b1 <- d$b*.3; d$b2 <- d$b*.7
  result <- spseg(d,coords=c("x","y"),variables=list(
    original=seg_grouped(counts=c("a","b","c")),
    split=seg_grouped(counts=c("a","b1","b2","c"))),
    measures=c("information","diversity","dissimilarity","square_root"),bands=c(0,2))
  expect_equal(result$analyses$original$indices,result$analyses$split$indices,tolerance=1e-10)
})

test_that("exposure cut selection affects only exposure, and preserves ties and RNG", {
  d <- data.frame(x=1:20,y=0,v=rep(c(0,5,10,20),c(12,3,3,2)))
  set.seed(701); seed <- .Random.seed
  result <- spseg(d,coords=c("x","y"),variables=list(
    deciles=seg_continuous(values="v"),
    values=seg_continuous(values="v",exposure=list(breaks=c(-Inf,6,Inf)))),
    bands=c(0,3),output="full")
  expect_identical(.Random.seed,seed)
  expect_equal(nrow(result$analyses$deciles$indices$multigroup$p[[1]]),4L)
  expect_equal(nrow(result$analyses$values$indices$multigroup$p[[1]]),2L)
  for (key in c("d","h","r","s")) expect_equal(
    result$analyses$deciles$indices$multigroup[[key]],
    result$analyses$values$indices$multigroup[[key]])
  expect_equal(rowSums(result$analyses$values$indices$multigroup$p[[2]]),c(1,1),ignore_attr=TRUE)
  unique <- transform(d, v=seq_len(nrow(d)))
  ten <- spseg(unique[1:10,],coords=c("x","y"),
                variables=list(v=seg_continuous(values="v")),bands=1)
  expect_equal(ten$analyses$v$definition$exposure$shares,rep(.1,10))
})

test_that("deterministic new-interface calls do not create an RNG state", {
  existed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if (existed) old <- get(".Random.seed",envir=.GlobalEnv)
  on.exit({
    if (existed) assign(".Random.seed",old,envir=.GlobalEnv) else
      if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
        rm(".Random.seed",envir=.GlobalEnv)
  })
  if (existed) rm(".Random.seed",envir=.GlobalEnv)
  d <- data.frame(x=1:4,y=0,v=c(0,0,1,2))
  spseg(d,coords=c("x","y"),variables=list(v=seg_continuous(values="v")),bands=2)
  expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
})

test_that("single populated intervals and near-endpoint grouped integrals are stable", {
  d <- data.frame(x=1:3,y=0,a=c(1,4,2),b=0)
  r <- spseg(d,coords=c("x","y"),variables=list(v=seg_grouped(counts=c("a","b"))),bands=2)
  expect_equal(as.data.frame(r)$value,rep(0,4),tolerance=1e-12)
  d$a <- c(1e-8,2e-8,3e-8); d$b <- c(1,2,3)
  d$c <- c(1e8,2e8,3e8)
  r <- spseg(d,coords=c("x","y"),variables=list(v=seg_grouped(counts=c("a","b","c"))),bands=2)
  expect_equal(as.data.frame(r)$value,rep(0,4),tolerance=1e-10)
})

test_that("mixed SF results, sparse output, and search backends work together", {
  d <- data.frame(x=c(0,0,1,3,4),y=0,g=c("a","b","a","c","c"),
                  o=ordered(c("low","low","mid","high","mid"),levels=c("low","mid","high")),
                  v=c(0,0,2,3,4))
  sfdata <- sf::st_as_sf(d,coords=c("x","y"),sf_column_name="location_geometry",crs=3857)
  variables <- list(g=seg_categorical(values="g"),o=seg_ordinal(values="o"),v=seg_continuous(values="v"))
  for (neighborhood in c("radius","knn")) {
    a <- spseg(sfdata,variables=variables,bands=c(1,3),neighbors=neighborhood,output="full")
    b <- spseg(d,coords=c("x","y"),variables=variables,bands=c(1,3),neighbors=neighborhood,
               search="brute",output="full")
    expect_equal(a$analyses,b$analyses,tolerance=1e-12)
    expect_equal(a$diagnostics$focal_searches,nrow(d))
    expect_s3_class(sf::st_as_sf(a,variables="v",bands=3),"sf")
    expect_true(all(c("variable","method","location","rank","share") %in%
                      names(as.data.frame(a,what="env"))))
    expect_true(all(c("from_group","to_group") %in% names(as.data.frame(a,what="exposure"))))
    expect_true(all(c("group_a","group_b") %in% names(as.data.frame(a,what="pairwise"))))
  }
})

test_that("populations, omissions, constants, and invalid specifications are explicit", {
  d <- data.frame(x=1:4,y=0,a=c(2,0,4,1),b=c(0,5,2,4),v=c(1,NA,1,1),pop=10)
  variables <- list(g=seg_categorical(counts=c("a","b")),v=seg_continuous(values="v",na.action="omit"))
  r <- spseg(d,coords=c("x","y"),population="pop",variables=variables,bands=2,output="full")
  expect_true(all(is.na(as.data.frame(r,variables="v")$value)))
  expect_equal(r$analyses$v$status, c("2"="constant"))
  expect_equal(r$population,d$pop)
  expect_equal(r$analyses$v$indices$multigroup$p[[1]][1,1],1)
  expect_error(spseg(d,coords=c("x","y"),variables=list(v=seg_continuous(values="v")),bands=2),"Missing")
  d$c <- 1
  expect_error(spseg(d,coords=c("x","y"),variables=list(
    a=seg_categorical(counts=c("a","b")),b=seg_categorical(counts=c("a","c"))),bands=2),
    "populations differ")
  expect_error(seg_continuous(values="v",counts="a"),"exactly one")
  expect_error(spseg(d,coords=c("x","y"),variables=variables,bands=2,surface="grid"),"counts")
  expect_error(spseg(d,coords=c("x","y"),variables=list(v=seg_ordinal(values="a")),bands=2),"ordered")
  expect_error(spseg(d,coords=c("x","y"),variables=list(v=seg_continuous(values="a")),bands=2,
                     measures="unknown"),"Unknown")
  env <- spseg(d,coords=c("x","y"),population="pop",variables=variables,bands=2,output="localenv")
  expect_error(as.data.frame(env),"not calculated")
})

test_that("kNN reaches 100 eligible people rather than stopping at 100 residents", {
  d <- data.frame(x=seq_len(140), y=0, race=rep(c("a","b"),70),
                  income=as.numeric(seq_len(140)))
  d$income[89:100] <- NA_real_
  variables <- list(race=seg_categorical(values="race"),
                    income=seg_continuous(values="income",na.action="omit"))
  r <- spseg(d,coords=c("x","y"),variables=variables,neighbors="knn",
             bands=100,weighting="unweighted",scope="multigroup",output="full")
  expected <- c(1:88,101:112)
  expect_equal(r$analyses$income$env[[1]][[1]]$value,expected)
  expect_equal(r$analyses$income$env[[1]][[1]]$share,rep(.01,100))
  for (name in names(variables)) {
    expect_equal(r$analyses[[name]]$knn$min_population,100)
    expect_equal(r$analyses[[name]]$knn$max_population,100)
  }
  adults <- d[!is.na(d$income),]
  income_only <- spseg(adults,coords=c("x","y"),neighbors="knn",bands=100,
                       variables=list(income=seg_continuous(values="income")),
                       weighting="unweighted",scope="multigroup")
  expect_equal(r$analyses$income$indices,income_only$analyses$income$indices,
               tolerance=1e-12)
  expect_equal(r$diagnostics$focal_searches,nrow(d))
  expect_identical(r$neighbors$population,"variable")
})

test_that("aggregate kNN populations and boundary fractions are variable-specific", {
  d <- data.frame(x=c(0,1,3,5,9),y=0,
                  race_a=c(4,10,9,12,15),race_b=c(6,10,31,8,25),
                  income_low=c(0,4,2,8,6),income_high=c(0,6,2,16,4))
  vars <- list(race=seg_categorical(counts=c("race_a","race_b")),
                income=seg_grouped(counts=c("income_low","income_high")))
  for (kernel in c("unweighted","biweight","inverse","exponential")) {
    r <- spseg(d,coords=c("x","y"),variables=vars,neighbors="knn",bands=c(7,20),
               weighting=kernel,output="full")
    for (name in names(vars)) {
      alone <- spseg(d,coords=c("x","y"),variables=vars[name],neighbors="knn",
                     bands=c(7,20),weighting=kernel,output="full")
      expect_equal(r$analyses[[name]],alone$analyses[[name]],tolerance=1e-12)
      expect_equal(r$analyses[[name]]$knn$min_population,c(7,20),tolerance=1e-12)
      expect_equal(r$analyses[[name]]$knn$max_population,c(7,20),tolerance=1e-12)
    }
    brute <- spseg(d,coords=c("x","y"),variables=vars,neighbors="knn",bands=c(7,20),
                   weighting=kernel,output="full",search="brute")
    expect_equal(r$analyses,brute$analyses,tolerance=1e-12)
    if (kernel=="unweighted") {
      # At location 1, race selects its 10 residents plus half of location 2;
      # income skips location 1 and selects 10 + 4 + 6 eligible people.
      expect_equal(unname(r$analyses$race$env[[2]][1,]),c(9,11)/20)
      expect_equal(unname(r$analyses$income$env[[2]][1,]),c(8,12)/20)
    }
    expect_equal(r$diagnostics$focal_searches,nrow(d))
  }
})

test_that("kNN reports unattainable eligible thresholds and audits exact focal truncation", {
  d <- data.frame(x=1:3,y=0,race_a=c(40,20,30),race_b=c(60,30,20),
                  income=c(1,NA,2))
  expect_error(spseg(d,coords=c("x","y"),neighbors="knn",bands=100,
    variables=list(race=seg_categorical(counts=c("race_a","race_b")),
                   income=seg_continuous(values="income",na.action="omit"))),
    "Variable 'income'.*population 2.*threshold 100")
  r <- spseg(d,coords=c("x","y"),neighbors="knn",bands=10,output="full",
    variables=list(race=seg_categorical(counts=c("race_a","race_b"))))
  expect_equal(r$analyses$race$knn$min_population,10)
  expect_equal(r$analyses$race$knn$max_population,10)
  expect_equal(r$analyses$race$knn$max_radius,0)
  expect_equal(unname(r$analyses$race$env[[1]][1,]),c(.4,.6))
  fractional <- data.frame(x=1:2000,y=0,a=.04,b=.06)
  r <- spseg(fractional,coords=c("x","y"),neighbors="knn",bands=200,
              variables=list(g=seg_categorical(counts=c("a","b"))))
  expect_equal(r$analyses$g$knn$min_population,200,tolerance=1e-12)
  expect_equal(r$analyses$g$knn$max_population,200,tolerance=1e-12)
})

test_that("shared candidate discovery includes distance ties across query boundaries", {
  d <- data.frame(x=rep(0,80),y=0,race=rep(c("a","b"),40),income=as.numeric(1:80))
  d$income[seq(3,80,3)] <- NA_real_
  vars <- list(race=seg_categorical(values="race"),
                income=seg_continuous(values="income",na.action="omit"))
  kd <- spseg(d,coords=c("x","y"),variables=vars,neighbors="knn",bands=5,output="full")
  brute <- spseg(d,coords=c("x","y"),variables=vars,neighbors="knn",bands=5,
                 output="full",search="brute")
  expect_equal(kd$analyses,brute$analyses,tolerance=1e-12)
})

test_that("fractional kNN boundaries do not expand the kernel radius through rounding", {
  xy <- cbind(x=0:2,y=0)
  counts <- cbind(a=c(1,0,2),b=c(0,3,4))
  coordinate_input <- spseg(xy, variables = list(groups = seg_categorical(counts = seq_len(ncol(counts)))), data = counts,neighbors="knn",bands=3,
                   weighting="exponential",output="full")
  # Location 1 selects itself and 2/3 of location 2, with radius 1.
  expected <- c(1,2*exp(-1/(1+sqrt(.Machine$double.eps))))
  expect_equal(unname(coordinate_input$analyses$groups$env[[1]][1,]),expected/sum(expected),tolerance=1e-12)
  d <- data.frame(xy,counts)
  multi <- spseg(d,coords=c("x","y"),neighbors="knn",bands=3,
                  weighting="exponential",output="full",
                  variables=list(g=seg_categorical(counts=c("a","b"))))
  expect_equal(unname(multi$analyses$g$env[[1]][1,]),expected/sum(expected),
               tolerance=1e-12)
})

test_that("shared kNN radius expansion retains nonzero-distance boundary candidates", {
  set.seed(707)
  n <- 1000L
  d <- data.frame(x=runif(n),y=runif(n),race=sample(letters[1:3],n,TRUE),
                  income=rnorm(n))
  d$income[seq(1,n,4)] <- NA_real_
  variables <- list(race=seg_categorical(values="race"),
                    income=seg_continuous(values="income",na.action="omit"))
  r <- spseg(d,coords=c("x","y"),variables=variables,neighbors="knn",
             bands=100,measures="information",scope="multigroup")
  brute <- spseg(d,coords=c("x","y"),variables=variables,neighbors="knn",
                 bands=100,measures="information",scope="multigroup",search="brute")
  expect_equal(r$analyses,brute$analyses,tolerance=1e-12)
  for (analysis in r$analyses) {
    expect_equal(analysis$knn$min_population,100)
    expect_equal(analysis$knn$max_population,100)
  }
})
