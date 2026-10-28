# Small all-pairs reference, independent of Boost traversal and stopping rules.
network_distances <- function(nodes, edges, weighted=FALSE, directed=FALSE) {
  d <- matrix(Inf,length(nodes),length(nodes),dimnames=list(nodes,nodes))
  diag(d) <- 0
  for (e in seq_len(nrow(edges))) {
    a <- as.character(edges$from[e]); b <- as.character(edges$to[e])
    cost <- if(weighted) edges$cost[e] else 1
    d[a,b] <- min(d[a,b],cost)
    if(!directed) d[b,a] <- min(d[b,a],cost)
  }
  for (k in seq_along(nodes)) d <- pmin(d,outer(d[,k],d[k,],"+"))
  d
}

network_environment <- function(counts, distances, target, knn=FALSE,
                                weighting="unweighted") {
  mass <- rowSums(counts)
  result <- counts * 0
  for (i in seq_len(nrow(counts))) {
    distance <- distances[i,]
    inclusion <- numeric(length(mass))
    radius <- target
    if(knn) {
      inclusion[i] <- if(mass[i]>0) min(1,target/mass[i]) else 0
      remaining <- target-inclusion[i]*mass[i]
      radius <- 0
      for (r in sort(unique(distance[is.finite(distance)]))) {
        if(remaining<1e-12) break
        shell <- which(distance==r & seq_along(mass)!=i)
        total <- sum(mass[shell])
        if(total<=0) next
        inclusion[shell] <- min(1,remaining/total)
        remaining <- remaining-sum(inclusion[shell]*mass[shell])
        radius <- r
      }
      stopifnot(remaining<1e-10)
    } else inclusion[distance<=target] <- 1
    weight <- inclusion
    if(weighting=="exponential") {
      bw <- if(radius>0) radius * if(knn) (1+sqrt(.Machine$double.eps)) else 1 else 1
      weight <- weight*exp(-distance/bw)
    }
    composition <- colSums(counts*weight)
    result[i,] <- composition/sum(composition)
  }
  result
}

test_that("network paths agree with Euclidean paths for every variable type", {
  d <- data.frame(node=letters[1:5],x=c(0,1,3,5,8),y=0,
    a=c(2,0,3,1,2),b=c(0,3,1,2,1),v=c(10,20,20,40,60),
    o=ordered(c(1,2,2,3,1)))
  edges <- data.frame(from=letters[1:4],to=letters[2:5],cost=diff(d$x))
  net <- spseg_network(edges,cost="cost",units="minutes")
  variables <- list(g=spseg_categorical(counts=c("a","b")),
    o=spseg_ordinal(values="o"),b=spseg_continuous(representation = "intervals", counts=c("a","b"),
      summary=list(probs=sum(d$a)/sum(d$a+d$b))),
    v=spseg_continuous(values="v"))
  for (mode in c("radius","knn")) for (kernel in c("unweighted","biweight","inverse","exponential")) {
    bands <- if(mode=="radius") c(0,2,5) else c(2,4)
    spatial <- spseg(d,variables,coords=c("x","y"),population=rep(1,5),
      neighbors=mode,bands=bands,weighting=kernel,output="full")
    relational <- spseg(d,variables,network=net,node_id="node",population=rep(1,5),
      neighbors=mode,bands=bands,weighting=kernel,output="full")
    expect_equal(relational$analyses,spatial$analyses,tolerance=1e-12)
    expect_equal(relational$neighbors$distance_units,"minutes")
    expect_equal(relational$neighbors$space,"network")
    expect_equal(relational$neighbors$units,if(mode=="knn")"population" else "minutes")
  }
})

test_that("weighted and unweighted neighborhoods match an independent reference", {
  set.seed(490)
  n <- 12
  edges <- rbind(data.frame(from=1:n,to=c(2:n,1)),
                 data.frame(from=sample(n,30,TRUE),to=sample(n,30,TRUE)))
  edges$cost <- sample(0:4,nrow(edges),TRUE)
  d <- data.frame(node=1:n,a=sample(0:4,n,TRUE),b=sample(1:3,n,TRUE),
    c=sample(0:2,n,TRUE),e=sample(1:2,n,TRUE))
  vars <- list(g=spseg_categorical(counts=c("a","b")),
               h=spseg_continuous(representation = "intervals", counts=c("c","e"),
                 summary=list(probs=sum(d$c)/sum(d$c+d$e))))
  for (directed in c(FALSE,TRUE)) for (weighted in c(FALSE,TRUE)) {
    net <- spseg_network(edges,cost=if(weighted)"cost" else NULL,directed=directed)
    distance <- network_distances(as.character(1:n),edges,weighted,directed)
    for (knn in c(FALSE,TRUE)) for (kernel in c("unweighted","exponential")) {
      bands <- if(knn)c(2.5,15) else c(0,1,3)
      r <- spseg(d,vars,network=net,node_id="node",bands=bands,
        neighbors=if(knn)"knn" else "radius",weighting=kernel,output="full")
      for (v in names(vars)) for (b in seq_along(bands)) {
        counts <- as.matrix(d[vars[[v]]$counts])
        expected <- network_environment(counts,distance,bands[b],knn,kernel)
        # Continuous summaries label rank bins rather than original categories.
        expect_equal(unname(r$analyses[[v]]$env[[b]]),unname(expected),tolerance=1e-12)
      }
      if(knn) for (v in names(vars)) {
        expect_equal(r$analyses[[v]]$knn$min_population,bands)
        expect_equal(r$analyses[[v]]$knn$max_population,bands)
      }
    }
  }
})

test_that("KNN completes hop and zero-cost boundary groups", {
  # Focal plus 10 first-, 70 second-, and 500 third-order neighbors.
  edges <- data.frame(from=c(rep(0,10),rep(1,70),rep(11,500)),to=1:580)
  d <- data.frame(node=0:580,a=c(rep(1,81),rep(0,500)),b=c(rep(0,81),rep(1,500)))
  r <- spseg(d,list(g=spseg_categorical(counts=c("a","b"))),
    network=spseg_network(edges),node_id="node",neighbors="knn",bands=100,
    weighting="unweighted",output="full")
  expect_equal(unname(r$analyses$g$env[[1]][1,]),c(.81,.19))
  expect_equal(r$analyses$g$knn$min_population,100)
  expect_equal(r$analyses$g$knn$max_population,100)
  edges <- data.frame(from=c("a","b","c"),to=c("b","c","d"),cost=c(1,0,0))
  d <- data.frame(node=letters[1:4],a=c(1,1,0,0),b=c(0,0,1,1))
  r <- spseg(d,list(g=spseg_categorical(counts=c("a","b"))),
    network=spseg_network(edges,cost="cost"),node_id="node",neighbors="knn",
    bands=2,weighting="unweighted",output="full")
  expect_equal(unname(r$analyses$g$env[[1]][1,]),c(2/3,1/3))
})

test_that("transit nodes, parallel edges and shared observation nodes conserve mass", {
  edges <- data.frame(from=c("a","t","a","a"),to=c("t","b","t","a"),cost=c(1,1,3,0))
  d <- data.frame(node=c("a","a","b"),g=c("A","B","B"))
  vars <- list(g=spseg_categorical(values="g"))
  r <- spseg(d,vars,network=spseg_network(edges,cost="cost"),node_id="node",
             bands=c(0,2),weighting="unweighted",output="full")
  expect_equal(unname(r$analyses$g$env[[1]][1,]),c(.5,.5))
  expect_equal(unname(r$analyses$g$env[[2]][1,]),c(1/3,2/3))
  expect_equal(r$network$nodes,3)
  k <- spseg(d,vars,network=spseg_network(edges,cost="cost"),node_id="node",
    neighbors="knn",bands=1.5,weighting="unweighted",output="full")
  expect_equal(unname(k$analyses$g$env[[1]][1,]),c(2/3,1/3))
  expect_equal(k$analyses$g$knn$min_population,1.5)
  expect_setequal(as.data.frame(k,what="env")$node,c("a","b"))
})

test_that("network neighborhoods preserve row and edge order invariance", {
  d <- data.frame(node=letters[1:6],a=c(10,0,5,0,2,1),b=c(0,20,0,10,3,4))
  edges <- data.frame(from=letters[1:5],to=letters[2:6])
  vars <- list(g=spseg_categorical(counts=c("a","b")))
  run <- function(d,edges) spseg(d,vars,network=spseg_network(edges),node_id="node",
    neighbors="knn",bands=25,weighting="exponential",output="full")
  r <- run(d,edges)
  p <- c(4,1,6,2,5,3)
  reordered <- run(d[p,],edges[5:1,])
  expect_equal(reordered$analyses$g$env[[1]],r$analyses$g$env[[1]][p,],tolerance=1e-12)
  expect_equal(reordered$analyses$g$indices,r$analyses$g$indices,tolerance=1e-12)
})

test_that("unreachable eligible populations produce informative errors", {
  edges <- data.frame(from="a",to="b")
  d <- data.frame(node=c("a","b","c"),g=c("A","B","A"),v=c(1,NA,2))
  net <- spseg_network(edges)
  expect_error(spseg(d,list(g=spseg_categorical(values="g")),network=net,
    node_id="node",neighbors="knn",bands=2),"Variable 'g'.*focal row 3.*node 'c'.*population is 1")
  expect_error(spseg(d,list(v=spseg_continuous(values="v",na.action="omit")),
    network=net,node_id="node",neighbors="knn",bands=2),
    "Variable 'v'.*node 'a'.*population is 1")
  d <- d[1:2,]
  expect_error(spseg(d,list(g=spseg_categorical(values="g")),
    network=spseg_network(edges,directed=TRUE),node_id="node",neighbors="knn",bands=2),
    "focal row 2.*node 'b'.*population is 1")
  isolated <- spseg_network(data.frame(from=character(),to=character()))
  r <- spseg(d,list(g=spseg_categorical(values="g")),network=isolated,
              node_id="node",bands=0,output="full")
  expect_equal(unname(r$analyses$g$env[[1]]),diag(2))
})

test_that("geometry is optional and never determines network distances", {
  d <- data.frame(node=c("a","b"),x=c(0,100),y=0,g=c("A","B"))
  net <- spseg_network(data.frame(from="a",to="b"))
  vars <- list(g=spseg_categorical(values="g"))
  r <- spseg(d,vars,network=net,node_id="node",bands=1,
             weighting="unweighted",output="full")
  expect_null(r$coords)
  expect_error(sf::st_as_sf(r),"network-only")
  spatial <- sf::st_as_sf(d,coords=c("x","y"),sf_column_name="shape",crs=3857)
  mapped <- spseg(spatial,vars,network=net,node_id="node",bands=1,
                  weighting="unweighted",output="full")
  expect_equal(mapped$analyses,r$analyses)
  expect_s3_class(sf::st_as_sf(mapped),"sf")
  coordinates <- spseg(d,vars,coords=c("x","y"),network=net,node_id="node",
    bands=1,weighting="unweighted",output="full")
  expect_s3_class(sf::st_as_sf(coordinates),"sf")
  small <- spseg(d,vars,network=net,node_id="node",bands=1)
  expect_null(small$data)
  expect_null(small$node_id)
  expect_equal(small$network$edges,1)
  expect_equal(small$network$units,"hops")
})

test_that("network inputs reject invalid or ambiguous specifications", {
  edges <- data.frame(from="a",to="b",cost=1)
  expect_error(spseg_network(transform(edges,cost=-1),cost="cost"),"nonnegative")
  expect_error(spseg_network(transform(edges,cost=NA_real_),cost="cost"),"finite")
  expect_error(spseg_network(transform(edges,from=NA_character_)),"nonmissing")
  expect_error(spseg_network(transform(edges,to="")),"empty")
  expect_error(spseg_network(edges,directed=NA),"directed")
  expect_error(spseg_network(edges,units=1),"units")
  d <- data.frame(node=c("a","b"),g=c("A","B"))
  vars <- list(g=spseg_categorical(values="g"))
  net <- spseg_network(edges)
  expect_error(spseg(d,vars,network=edges,node_id="node",bands=1),"spseg_network")
  expect_error(spseg(d,vars,network=net,bands=1),"node_id")
  expect_error(spseg(d,vars,network=net,node_id="node"),"explicit.*bands")
  expect_error(spseg(d,vars,network=net,node_id="node",bands=1,surface="grid"),"raw")
  expect_error(spseg(d,vars,network=net,node_id="node",bands=1,search="brute"),"boost")
  expect_error(spseg(d,vars,node_id="node",bands=1),"requires a network")
})
