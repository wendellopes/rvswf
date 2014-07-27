\name{lcfa.cyl}
\alias{lcfa.cyl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lcfa.cyl(nmax, x, code = "C", NMAX = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nmax}{
%%     ~~Describe \code{nmax} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{code}{
%%     ~~Describe \code{code} here~~
}
  \item{NMAX}{
%%     ~~Describe \code{NMAX} here~~
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
function (nmax, x, code = "C", NMAX = 200) 
{
    if (!code \%in\% c("C", "R")) {
        stop("Code must be \"C\" or \"R\"")
    }
    if (code == "C") {
        dummy <- rep(0, nmax + 1)
        u <- .C("lcfa_cyl", nmax = as.integer(nmax), x = as.double(x), 
            gn = as.double(dummy), Dn = as.double(dummy), NMAX = as.integer(NMAX))
        return(data.frame(gn = u$gn, Dn = u$Dn))
    }
    else {
        S <- function(n, x) {
            return(n/x)
        }
        Dn <- rep(0, nmax + 1)
        gn <- rep(1, nmax + 1)
        Dn[nmax + 1] <- lcfe.cbl(nmax, x)
        gn[nmax + 1] <- lcfe.cbd(nmax, x)
        for (n in nmax:1) {
            gn[n] <- S(n, x) + Dn[n + 1]
            Dn[n] <- S(n - 1, x) - 1/gn[n]
        }
        return(data.frame(gn, Dn))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
