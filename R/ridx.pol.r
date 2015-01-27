#' Parameters: Refraction index of poliestirene.
#' 
#' @param l The wavelength in micrometers
#' @seealso \code{\link{ridx.h2o}}
#' @examples
#' l<- .532
#' ridx.pol(l)
# INDEX OF REFRACTION POLIESTIRENE
ridx.pol<-function(l){
   return(1.5725+.0031080/l^2+.00034779/l^4)
}