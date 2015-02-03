#' Parameters: Refraction index of water
#' 
#' @param l The wavelength in micrometers
#' @return The index of refraction of the water as function of the wavelength
#' @export
#' @seealso \code{\link{ridx.h2o}}
#' @examples
#' l<-seq(.38,.75,.01)
#' u<-ridx.h2o(l)
#' plot(l,u,type='l')
ridx.h2o<-function(l){
   return(1.3253+.0027553/l^2+.00003779/l^4)
}