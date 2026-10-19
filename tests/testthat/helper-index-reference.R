# Independent nominal formulas applied to stored local compositions.
categorical_reference <- function(counts, local) {
  counts <- as.matrix(counts)
  local <- as.matrix(local)
  scalar <- function(n, p) {
    mass <- rowSums(n)
    keep <- mass > 0
    n <- n[keep, , drop = FALSE]; p <- p[keep, , drop = FALSE]
    p <- p / rowSums(p)
    mass <- mass[keep]; alpha <- mass / sum(mass)
    shares <- colSums(n) / sum(n)
    entropy <- function(v) -sum(v[v > 0] * log(v[v > 0]))
    e <- entropy(shares); diversity <- sum(shares*(1-shares))
    active <- shares > 0
    list(d = sum(alpha*rowSums(abs(sweep(p, 2, shares))))/(2*diversity),
         r = 1-sum(alpha*rowSums(p*(1-p)))/diversity,
         h = 1-sum(alpha*apply(p, 1, entropy))/e,
         s = 1-sum(alpha*exp(rowMeans(log(p[,active,drop=FALSE])))) /
           exp(mean(log(shares[active]))))
  }
  multigroup <- scalar(counts, local)
  p <- local/rowSums(local)
  exposure <- t(counts) %*% p / colSums(counts)
  dimnames(exposure) <- list(colnames(counts), colnames(counts))
  multigroup$p <- list(exposure)
  matrices <- lapply(c("d","r","h","s"), function(key) {
    out <- matrix(NA_real_,ncol(counts),ncol(counts),
                   dimnames=list(colnames(counts),colnames(counts)))
    for (a in seq_len(ncol(counts)-1L)) for (b in seq.int(a+1L,ncol(counts))) {
      if (all(colSums(counts[,c(a,b),drop=FALSE]) > 0))
        out[a,b] <- out[b,a] <- scalar(counts[,c(a,b),drop=FALSE],
                                      local[,c(a,b),drop=FALSE])[[key]]
    }
    list(out)
  })
  names(matrices) <- c("d","r","h","s")
  list(multigroup=multigroup,pairwise=matrices)
}
