\name{autosome.genome.scan}
\alias{autosome.genome.scan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function performs a QTL scan on autosomes
}
\description{
Function performs a scan for QTL on autosomes. 
}
\usage{
autosome.genome.scan(pheno, data, factors = NULL, covariates = NULL, 
                     bk.QTL = NULL, error = normal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
  phenotype name
}
  \item{data}{
data in the internal MAPfastR format (with calculated line origin probabilities)
}
  \item{factors}{
vector of factors characterising the trait (if any)
}
  \item{covariates}{
vector of covariate names
}
  \item{bk.QTL}{...
%%     ~~Describe \code{bk.QTL} here~~
}
  \item{error}{...
%%     ~~Describe \code{error} here~~
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
L. Crooks, C. Nettelblad, and O. Carlborg, "An improved method for estimating chromosomal line origin in QTL analysis of crosses between outbred lines," G3: Genes Genomes Genetics, iss. 1: 57-64, 2011.
}
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
##---- See Package Vignettes ----
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
