#' Understand Phred quality
#'
#' Convert between Phred quality index and error probablity
#' @author Elves H Duarte
#' @param val A error (range 0 to 1) or quality (greather than 0) value to be converted whether to Phred quality or error probability.
#' @param get Specify to what you want to convert the value: 'error' convert a Phred value to error and 'quality' convert a desired error probability to Phred quality
#' @return A error probability of a Phred quality or the Phred quality associated with a error.
#' @export
phred=function(val, get="error"){
  if(get=="error")
  {
    if(val < 0)
      stop("Quality value have to be greather than 0.")
    return(10^(val/-10))
  }
  if(get=="quality")
  {
    if(val < 0 | val > 1)
      stop("0 < error < 1")
    return(-10*log(val)/log(10))
  }
  else
    stop("get option must be 'error' or 'quality'")
}
