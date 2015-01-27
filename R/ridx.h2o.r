#' Parameters: Refraction index of water
#' 
#' @param l The wavelength in micrometers
#' @seealso \code{\link{ridx.h2o}}
#' @examples
#' l<- .532
#' ridx.h2o(l)
ridx.h2o<-function(l){
   return(1.3253+.0027553/l^2+.00003779/l^4)
}