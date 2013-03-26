\name{run_NOIA}
\alias{run_NOIA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run NOIA model within MAPfastR.
}
\description{
A wrapper function that allows running of the NOIA model using internal MAPfastR data format}
\usage{
run_NOIA(data, phenotype, markers, noia_func = "linear")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
a data object in the internal MAPfastR format, result returned by import_data function
}
\item{phenotype}{
name of the phenotype of interest
}
\item{markers}{
a vector of markers to be used
}
\item{noia_func}{
a function to be used in the NOIA model: linear or multilinear
}
 
}
\details{
For details see the documentation of NOIA package.
}
\value{...
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Alvarez-Castro JM, Carlborg O. (2007). A unified model for functional and statistical epistasis and its application in quantitative trait loci analysis. Genetics 176(2):1151-1167.

Alvarez-Castro JM, Le Rouzic A, Carlborg O. (2008). How to perform meaningful estimates of genetic effects. PLoS Genetics 4(5):e1000062.

Cheverud JM, Routman, EJ. (1995). Epistasis and its contribution to genetic variance components. Genetics 139:1455-1461.

Hansen TF, Wagner G. (2001) Modeling genetic architecture: A multilinear theory of gene interactions. Theoretical Population Biology 59:61-86.

Le Rouzic A, Alvarez-Castro JM. (2008). Estimation of genetic effects and genotype-phenotype maps. Evolutionary Bioinformatics 4.

Zeng ZB, Wang T, Zou W. (2005). Modelling quantitative trait loci and interpretation of models. Genetics 169: 1711-1725.
}
\author{
Mats Pettersson
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
%\keyword{ import }
%\keyword{ cnf2freq }
%\keyword{ crimap }% __ONLY ONE__ keyword per line
