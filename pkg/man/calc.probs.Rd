\name{calc.probs}
\alias{calc.probs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate line origin probabilities
}
\description{
Uses triM to calculate line origin probabilities at specified locations.
}
\usage{
calc.probs(data, interval = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
dataset in the internal MAPfastR format
}
  \item{interval}{
Either the interval, in cM, at which you want the probabilities estimated, or "markers", in which case probabilities will be estimated at marker positions.
}
}
\details{...
%%  ~~ If necessary, more details than the description above ~~
}
\value{...
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
L. Crooks, C. Nettelblad, and O. Carlborg, "An improved method for estimating chromosomal line origin in QTL analysis of crosses between outbred lines," G3: Genes Genomes Genetics, iss. 1: 57-64, 2011.}
\author{
L. Crooks
}
\note{...
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{...
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## ---- See Package Vignettes ----
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
