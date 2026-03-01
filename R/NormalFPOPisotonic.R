NormalFPOPisotonic <- function
(data.vec,
 penalty=NULL
){
  n.data <- length(data.vec)
  stopifnot(2 <= n.data)
  stopifnot(is.numeric(data.vec))
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty) == 1)
  stopifnot(0 <= penalty)
  cost.vec <- double(n.data)
  end.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.vec <- integer(n.data)
  result.list <- .C(
    "NormalFPOPisotonic_interface",
    data.vec=as.double(data.vec),
    n.data=as.integer(n.data),
    penalty=as.numeric(penalty),
    cost.vec=as.double(cost.vec),
    end.vec=as.integer(end.vec),
    mean.vec=as.double(mean.vec),
    intervals.vec=as.integer(intervals.vec),
    PACKAGE="PeakSegOptimal")
  n.segments <- sum(result.list$end.vec != -2L)
  if(n.segments == 0) n.segments <- 1L
  seg.mean <- rev(result.list$mean.vec[1:n.segments])
  seg.end <- rev(result.list$end.vec[1:n.segments])
  result.list$seg.mean <- seg.mean
  result.list$seg.end <- seg.end + 1L
  result.list$n.segments <- n.segments
  result.list
}
