#' Parameters: Refraction index of poliestirene.
#' 
#' @param l The wavelength in micrometers
#' @return The index of refraction of the water as function of the wavelength
#' @export
#' @seealso \code{\link{ridx.h2o}}
#' @examples
#' l<-seq(.38,.75,.01)
#' u<-ridx.pol(l)
#' plot(l,u,type='l')
ridx.pol<-function(l){
   return(1.5725+.0031080/l^2+.00034779/l^4)
}