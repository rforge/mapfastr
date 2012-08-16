\name{FIA.scan}
\alias{FIA.scan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance component QTL scan using FIA 
}
\description{
FIA (flexible intercross analysis, Ronnegard et al. 2008 Genetics) does genome-wide QTL scan using fast score test based on linear mixed models.
}
\usage{
FIA.scan(data = NULL, phenotype = NULL, fixed.effects = NULL, 
         chr = 1:2, estimate.ro = FALSE, print.score = TRUE, 
         pdf.figure = TRUE, figure.file = "FIA_scan.pdf")
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
A vector of indices of chromosomes to be scanned.
}
  \item{estimate.ro}{
A logical value indicating whether the within-line segregation would be estimated (See Ronnegard et al. 2008 Genetics).
}
  \item{print.score}{
A logical value indicating whether the score will be printed during calculation.
}
  \item{pdf.figure}{
A logical value indicating whether pdf figures will be generated.
}
  \item{figure.file}{
A string of pdf file name for figures.
}
}
\details{
This function requires IBD matrices calculated beforehand. See, e.g., \code{MCIBD.chro} or \code{MCIBD.geno}.
}
\value{
A list of data frames with score statistics, test for negativeness of gradient, and their positions.
}
\references{
Ronnegard, L., Besnier, F. and Carlborg, O (2008). An Improved Method for Quantitative Trait Loci Detection and Identification of Within-Line Segregation in F2 Intercross Designs. \emph{Genetics}, \bold{178}(4), 2315-2326
}
\author{
Xia Shen
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MCIBD}}, \code{\link{MCIBD.chro}}, \code{\link{MCIBD.genome}}
}
\examples{
##---- See Package Vignettes ----
}
