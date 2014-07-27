\name{lcfe.cbl}
\alias{lcfe.cbl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lcfe.cbl(n, x, NMAX = 2000, code = "C")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{NMAX}{
%%     ~~Describe \code{NMAX} here~~
}
  \item{code}{
%%     ~~Describe \code{code} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, x, NMAX = 2000, code = "C") 
{
    nmaxo <- NMAX
    fn <- 0
    if (!code \%in\% c("C", "R")) {
        stop("Code must be \"C\" or \"R\"")
    }
    if (code == "C") {
        if (is.complex(x)) {
            u <- .C("lcfc_cbl", n = as.integer(n), x = as.complex(x), 
                NMAX = as.integer(NMAX), fn = as.complex(fn))
        }
        else {
            u <- .C("lcfe_cbl", n = as.integer(n), x = as.double(x), 
                NMAX = as.integer(NMAX), fn = as.double(fn))
        }
        if (u$NMAX > nmaxo) {
            stop(paste("The accuracy was not reached in", NMAX, 
                "iterations!"))
        }
        return(u$fn)
    }
    else {
        n <- as.integer(n)
        eo <- .Machine$double.xmin
        ACC <- 10^-50
        fn <- lcfe.afs(n, x)
        if (abs(fn) < eo) {
            fn <- eo
        }
        Pn <- fn
        Qn <- 0
        j <- 0
        Dn <- 10
        while (abs(Dn - 1) > ACC) {
            if (j > NMAX) {
                NMAX <- j
                break
            }
            j <- j + 1
            an <- -1
            bn <- lcfe.afs(2 * (n + j), x)
            Pn <- bn + an/Pn
            if (abs(Pn) < eo) {
                Pn <- eo
            }
            Qn <- bn + an * Qn
            if (abs(Qn) < eo) {
                Qn <- eo
            }
            Qn <- 1/Qn
            Dn <- Pn * Qn
            fn <- fn * Dn
        }
        if (NMAX > nmaxo) {
            stop(paste("The accuracy was not reached in", NMAX, 
                "iterations!"))
        }
        return(fn)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
