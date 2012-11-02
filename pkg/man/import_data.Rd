\name{import_data}
\alias{import_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Import data into the internal MAPfastR format.
}
\description{
Function used to import data stored in cnf2freq or cri-map format into the internal MAPfastR format.
}
\usage{
import_data(pedigree_file = "pedigree.txt", phenotype_file = "phenotypes.txt", 
            genotype_file = "marker_data.txt", marker_file = "marker_info.txt", 
            backcross = FALSE, backcross.line = NA, backcross.parent = NA, 
            heterogam = 1, sex.chrom = "last", sex.restrict = FALSE, 
            n_gen_map = 2, ref.map = "average", 
            format = "cnf2freq", crimap_data_dir = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedigree_file}{
file containing pedigree information in Cri-Map or cnf2freq format.
}
  \item{phenotype_file}{
file containing phenotype information in Cri-Map or cnf2freq format.
}
  \item{genotype_file}{
file containing genotyping information in Cri-Map or cnf2freq format.
}
  \item{marker_file}{
file containing information about markers in cnf2freq format.
}
  \item{backcross}{
parameter determining whether backcross has been used (1) or not (0)
}
  \item{backcross.line}{
backcross line id
}
  \item{backcross.parent}{
backcross parent id
}
  \item{heterogam}{
heterogametic sex (number of X chromosomes in females)
}
  \item{sex.chrom}{
number of the sex chromosome, by defoult it is the last one
}
  \item{sex.restrict}{
if backcross has been sex restricted, set to 
}
  \item{n_gen_map}{
number of input genetic maps used (2 or 1) VIP!!!
}
  \item{ref.map}{
parameter determining whether genetic map is based on average distance, map from dames or map from sires
}
  \item{format}{
parameter specifying whether cnf2freq or Cri-Map format has been used
}
  \item{crimap_data_dir}{
parameter specifying directory with Cri-Map data
}
}
\details{
The function reads a set of input files, in either triM or crimap format, and stores the data in a MAPfastR data object. Downstream analyses will extend the object with their results, while leaving the input data intact.
}
\value{...
A list, constituting a MAPfastR data object, containing data from the input files as well as the specified values for the input parameters.
}
\references{
Green P, Fallsand K, Crooks S: CRI-MAP Documentation version 2.4. Washington University School of Medicine St. Louis, MO. (1990).

Nettelblad, C., Holmgren, S. & Crooks, L. cnF2freq: Efficient Determination of Genotype and Haplotype Probabilities in Outbred Populations Using Markov Models. Bioinformatics and Computational Biology 307-319 (2009).
}
\author{
Mats Pettersson
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{...
\code{\link{calc.probs}}
}
\examples{
## ---- See Package Vignettes ----
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ import }
\keyword{ cnf2freq }
\keyword{ crimap }% __ONLY ONE__ keyword per line
