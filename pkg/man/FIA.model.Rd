\name{FIA.model}
\alias{FIA.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitting FIA model for estimating within-line correlation 
}
\description{
\code{FIA.model} estimates within-line correlation based on the flexible intercross analysis (FIA) model (Ronnegard et al. 2008 Genetics)
}
\usage{
FIA.model(data = NULL, phenotype = NULL, fixed.effects = NULL, 
          chr = 1, IBD.index = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{data}{
The "MAPfastR" data object.
}
  \item{phenotype}{
A vector of phenotype to be analyzed.
}
  \item{fixed.effects}{
A vector or matrix with columns of fixed effects to be included in the linear mixed model.
}
  \item{chr}{
An index of the chromosome that the fitted locus is located on.
}
  \item{IBD.index}{
An index of the IBD matrix of the fitted locus.
}
}
\details{
Before fitting a FIA model, one needs to run \code{MCIBD.chro} or \code{MCIBD.genome} twice with different settings of \code{segregation} argument. The first run should be default setting, whereas the second should assume within-line fixation (see \bold{Details} in \code{MCIBD} for details).
}
\value{
\item{model}{
A FIA random effect model object fitted by \code{hglm}.
}
  \item{rho}{
An estimate of within-line correlation (see the reference paper for details).
}
}
\references{
Ronnegard, L., Besnier, F. and Carlborg, O (2008). An Improved Method for Quantitative Trait Loci Detection and Identification of Within-Line Segregation in F2 Intercross Designs. \emph{Genetics}, \bold{178}(4), 2315-2326
}
\author{
Xia Shen
}
\seealso{
\code{\link{MCIBD.chro}}, \code{\link{MCIBD.genome}}
}
\examples{
##---- See Package Vignettes ----
}