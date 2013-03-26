\name{MCIBD}
\alias{MCIBD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimating Identity-By-Descent (IBD) Matrices using Monte Carlo Sampling
}
\description{
\code{MCIBD} estimates the IBD matrix at a given locus or the epistatic IBD matrix for two linked loci. Monte Carlo sampling is used to approach to the matrix estimator, where the output file from triM is required. Segregation of the founder alleles can be set up in the estimation. Parallelzation is available using \code{snowfall} package if executing on a multi-core computer or cluster.
}
\usage{
MCIBD(loci, data = NULL, trim.output = NULL, chr = 1, 
      IBD.type = "genotypic", output.Z = "none", segregation = NULL, 
      mc.size = 99, hpc = FALSE, n.cpus = 2) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loci}{typically an integer specifying a test locus, or a vector containing two linked loci.}
  \item{data}{data object.}
  \item{trim.output}{the output probabilities from triM.}
  \item{chr}{chromosome index.}
  \item{IBD.type}{a string specifying the output type of IBD matrix, \code{"genotypic"} or \code{"gametic"}.}
  \item{output.Z}{a string indicates the type of output. If \code{"av"}, average incidence matrix Z is estimated instead of the corresponding IBD matrix; If \code{"all"}, the IBD matrix is still estimated as default but with all the incidence matrix imputes output in folder "Zall"; If \code{"pc"}, an incidence matrix containing the principle components of the corresponding IBD matrix is output.}
  \item{segregation}{a vector specifying the segregation of founder alleles. (See Details.)}
  \item{mc.size}{an integer setting the sample size of the Monte Carlo sampling for one CPU.}
  \item{hpc}{logical. If \code{TRUE}, high performance computing is carried out by parallelization.}
  \item{n.cpus}{an integer telling the number of cores that parallelization is executed on. Only useful when \code{hpc} is \code{TRUE}. Not recommended when mc.size is small.}
}
\details{
To set up the \code{segregation} of the founder alleles, suppose that we have \emph{m} founder(s) from one line and \emph{n} from another line. Then in total there are \emph{2 * (m + n)} alleles needed to be set up. According to the order in the pedigree information file, these alleles for the founders should be filled in \code{segregation} as a vector with length \emph{2 * (m + n)}. For instance, if there is one male in one line and three females in another line, the vector might be created like \code{c(1,1,2,2,3,3,2,3)}, which means that the male has two identical alleles whereas the genotypes of the females are (2,2), (3,3) and (2,3), respectively. If \code{NULL}, all the founder alleles are assumed to be different from each other.
}
\value{
An output IBD matrix is saved in the working directory as a file named by the loci with extension \code{.ibd}.
}
\references{
Shen, X., Nettelblad, C., Ronnegard, L. and Carlborg, O. (2009). \emph{Flexible Identity-By-Descent Matrix Estimation with Given Base Generation Structures in F2 Intercross Designs}. \bold{Package 'MCIBD' Vignette}. \cr \cr
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