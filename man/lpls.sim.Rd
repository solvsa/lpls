\name{lpls.sim}
\alias{lpls.sim}
\title{
Function for simulating data for testing the lpls-code
}
\description{
This function simulates three matrices suitable as input to lpls-regression.
}
\usage{
lpls.sim(q = 5, n = 30, p = 20, m = 6, comps = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{
Number of columns of X1
}
  \item{n}{
Number of rows of X1 and X2
}
  \item{p}{
Number of columns of X2 and rows of X3
}
  \item{m}{
Number of columns of X3
}
  \item{comps}{
Number of relevant components.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing X1, X2 and X3.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Solve Saebo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
testdata <- lpls.sim()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ lpls-data }
\keyword{ simulation }% __ONLY ONE__ keyword per line
