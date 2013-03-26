\name{MAPfastR-package}
\alias{MAPfastR-package}
\alias{MAPfastR}
\docType{package}
\title{
QTL Mapping in Outbred Line Crosses
}
\description{
MAPfastR is a fast and comprehensive software package for analyzing QTL data from outbred line-crosses that is developed for flexible analyses of large datasets. It includes an online developer and community based support system. MAPfastR is designed to efficiently handle outbred line cross data and has been shown to make efficient use of all available data from dense SNP-chips while remaining computationally efficient (Crooks et al., 2011). MAPfastR provide functionality for both standard line-cross analyses assuming within-line QTL fixation (Crooks et al., 2011), line-cross analyses allowing for within-line segregation (FIA; Ronnegard et al., 2008) and tests for epistatic interactions (Carlborg and Andersson, 2002). In addition to the standard functionality, the software comes with additional add-on packages that allow more experienced users take advantage of modules for analyzes of deep (Advanced Intercross Line) pedigrees. MAPfastR is implemented in the R language (with optimization of the more computationally intensive algorithms in C++), accepts several standard input formats and is available for Windows, Linux and MacOS.
}
\details{
\tabular{ll}{
Package: \tab MAPfastR\cr
Type: \tab Package\cr
Version: \tab 1.0-3\cr
Date: \tab 2013-01-24\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: \tab hglm, noia, cnf2freq, cnf2freqibdtracer, sfsmisc, snow, snowfall\cr
}
}
\author{
MAPfastR Development Team 

Maintainer: MAPfastR Developers <mapfastr@googlegroups.com>
}
\references{
J Alvarez-Castro and O Carlborg. A general model for functional and statistical epistasis and its application in QTL analysis. Genetics, (176), 2007. \cr

J Alvarez-Castro, A Le Rouzic, and O Carlborg. How to perform meaningful estimates of genetic effects. PLoS Genetics, (4), 2008. \cr

O Carlborg and L Andersson. The use of randomization testing for detection of multiple epistatic QTL. Genetical Research, (79), 2002. \cr

L Crooks, C Nettelblad, and O Carlborg. An improved method for estimating chromosomal line origin in QTL analysis of crosses between outbred lines. G3: Genes Genomes Genetics, (1: 57-64), 2011. \cr

RM Nelson, C Nettelblad, ME Pettersson, X Shen, L Crooks, F Besnier, J Alvarez Castro, L Ronnegard, W Ek, Z Sheng, M Kierczak, S Holmgren, and O Carlborg. MAPfastR: QTL mapping in outbred line crosses. Bioinformatics, 2012.\cr

C Nettelblad, S Holmgren, L Crooks, and O Carlborg. cnf2freq: Efficient determination of genotype and haplotype probabilities in outbred populations using markov models. In Bioinformatics and Computational Biology, number 5462 in Lecture notes in Bioinformatics. Springer, 2009.\cr

L Ronnegard, F Besnier, and O Carlborg. An improved method for QTL detection and identification of within-line segregation in F2 intercross designs. Genetics, (178), 2008.\cr

X Shen, L Ronnegard and O Carlborg. How to deal with genotype uncertainty in variance component QTL analyses. Genetics Research, Cambridge, (93), 2011.\cr
}
\keyword{ MAPfastR, outbred, QTL, MCIBD, FIA, NOIA }
