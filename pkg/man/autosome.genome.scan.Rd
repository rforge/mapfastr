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
  \item{bk.QTL}{
Pre identified QTL to be used as covariate. Default = NULL
}
  \item{error}{
Distribution of the residuals in the model. Normal or Binomial
}
}
\details{...
Performs a QTL scan, using the Haley Knott regression on the line origin in /data. F values, the test statistic from the Haley Knott regression, for each  interval is generated.
}
\value{...
List containing: max.F.value, position.max.F, estimated.effects.max.F, variables.in.starting.model, all.F.values, chrom.boundaries, results.by.chromosome
}
\references{
L. Crooks, C. Nettelblad, and O. Carlborg, "An improved method for estimating chromosomal line origin in QTL analysis of crosses between outbred lines," G3: Genes Genomes Genetics, iss. 1: 57-64, 2011.
}
\author{
Lucy Crooks
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{...
\code{\link{calc.probs}}
}
\examples{
##---- See Package Vignettes ----
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
