# Test/fixture references, independent of package rank tables and grouping.
# Callers provide the actual frozen metric distances; geometry has separate
# cross-engine tests. These ordinary-scale references use direct binomials,
# Beta CDF differences, explicit averaging and complete deletion reconstruction.
.distance_reference_table <- function(n, k, s, method) {
  if (method == "exact") return(vapply(seq_len(n), function(r) {
    sum(vapply(seq_len(k), function(j) {
      choose(r-1L,j-1L)*choose(n-r,s-j)
    }, numeric(1L))) / (k*choose(n,s))
  }, numeric(1L)))
  rowMeans(vapply(seq_len(k), function(j) {
    diff(pbeta((0:n)/n,j,s-j+1L))
  }, numeric(n)))
}

.distance_reference_weights <- function(dx, k = 1L, s = NULL, method = "exact",
                                         s1 = NULL, s2 = NULL, dim = NULL) {
  dx <- as.matrix(dx)
  n <- nrow(dx)
  if (is.null(s1)) table <- .distance_reference_table(n,k,s,method)
  else {
    coefficient <- -1/((s2/s1)^(2/dim)-1)
    table <- coefficient*.distance_reference_table(n,1L,s1,method) +
      (1-coefficient)*.distance_reference_table(n,1L,s2,method)
  }
  out <- matrix(0,n,ncol(dx))
  for (j in seq_len(ncol(dx))) {
    assigned <- numeric(n)
    assigned[order(dx[,j],seq_len(n))] <- table
    for (i in seq_len(n)) out[i,j] <- mean(assigned[dx[,j]==dx[i,j]])
  }
  out
}

.distance_reference_fit <- function(dx, signal, k = 1L, s = NULL, method = "exact",
                                    s1 = NULL, s2 = NULL, dim = NULL, ratio = FALSE) {
  dx <- as.matrix(dx); signal <- as.matrix(signal)
  n <- nrow(dx); m <- ncol(dx); q <- ncol(signal)
  weight <- function(d) .distance_reference_weights(d,k,s,method,s1,s2,dim)
  W <- weight(dx)
  components <- crossprod(signal,W)
  theta <- if (ratio) components[1L,]/components[2L,] else NULL
  size <- if (ratio) m else q*m
  # Explicit dimensions also preserve n x 1 when both q and m equal one.
  loo <- t(matrix(vapply(seq_len(n), function(i) {
    value <- crossprod(signal[-i,,drop=FALSE],weight(dx[-i,,drop=FALSE]))
    if (ratio) as.numeric((value[1L,]-theta*value[2L,])/components[2L,]) else as.numeric(value)
  }, numeric(size)), nrow=size, ncol=n))
  covariance <- crossprod(sweep(loo,2L,colMeans(loo)))*(n-1L)/n
  list(W=W, point=if(ratio)matrix(theta,nrow=1L)else components,
       covariance=covariance, se=matrix(sqrt(diag(covariance)),nrow=if(ratio)1L else q), loo=loo)
}

.distance_reference_check <- function(actual, expected, label, tolerance = 2e-12) {
  actual <- as.matrix(actual); expected <- as.matrix(expected)
  if (length(actual)!=length(expected) || any(!is.finite(actual)) ||
      any(!is.finite(expected)) ||
      max(abs(as.numeric(actual)-as.numeric(expected)))>tolerance*max(1,abs(expected))) {
    stop("Independent distance-averaging reference mismatch: ",label,call.=FALSE)
  }
}

.distance_reference_l2 <- function(X,X0) {
  vapply(seq_len(nrow(X0)),function(j) {
    sqrt(rowSums(sweep(X,2L,X0[j,],"-")^2))
  },numeric(nrow(X)))
}

.distance_reference_source_metadata <- function(root) {
  rfiles <- list.files(file.path(root,"R"),"[.]R$",full.names=FALSE)
  cppfiles <- list.files(file.path(root,"src"),"[.](cpp|h|hpp)$",recursive=TRUE,full.names=FALSE)
  headers <- list.files(file.path(root,"inst/include"),"[.](h|hpp)$",recursive=TRUE,full.names=FALSE)
  relative <- sort(c("DESCRIPTION",file.path("R",rfiles),file.path("src",cppfiles),
    file.path("inst/include",headers)))
  inventory <- tempfile("distance-reference-source-")
  on.exit(unlink(inventory),add=TRUE)
  writeLines(paste(relative,tools::md5sum(file.path(root,relative))),inventory)
  c(paste0("// R base commit (working source may differ): ",
      system2("git",c("-C",shQuote(root),"rev-parse","HEAD"),stdout=TRUE)),
    paste0("// R source version: ",read.dcf(file.path(root,"DESCRIPTION"))[1L,"Version"]),
    paste0("// R source fingerprint (sorted source path/MD5 inventory): ",tools::md5sum(inventory)))
}
