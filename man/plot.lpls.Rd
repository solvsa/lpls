\name{plot.lpls}
\alias{plot.lpls}
\title{
Function for making correlation loadings plot from lpls-regression
}
\description{
This function enables plotting of correlation loadings from an endo- or exo-lpls object.
}
\usage{
plot.lpls(fit, comps = c(1, 2), doplot = c(TRUE, TRUE, TRUE), level = c(2, 2, 2), arrow = c(1, 0, 1), xlim = c(-1, 1), ylim = c(-1, 1), samplecol = 4, pathcol = 2, varcol = "grey70", varsize = 1, sampleindex = 1:dim(fit$R22)[1], pathindex = 1:dim(fit$R3)[1], varindex = 1:dim(fit$R21)[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
A fitted lpls-object
}
  \item{comps}{
A vector of two elements indicating which components to plot.
}
  \item{doplot}{
A logical vector of three elements corresponding to X1, X2 and X3, respectively. A value \code{TRUE} will plot the loadings for the corresponding matrix.
}
  \item{level}{
A vector of length three giving the plot level for each matrix; \code{level}=1 gives dots, whereas \code{level}=2 uses dimnames as plotting symbols.
}
  \item{arrow}{
A vector of length three indicating the use of arrows from the origin to the loading value (1), or no arrow (0) for each matrix.
}
  \item{xlim}{
The plot limits of the x-axis
}
  \item{ylim}{
The plot limits of the y-axis
}
  \item{samplecol}{
The color for the R22-values (samples of X2) 
}
  \item{pathcol}{
The color for R3
}
  \item{varcol}{
The color for the R21-values (variables of X2)
}
  \item{varsize}{
Size of plotting symbol for the R21-values (variables of X2)
}
  \item{sampleindex}{
Which samples to plot
}
  \item{pathindex}{
Which columns of X3 to plot
}
  \item{varindex}{
Which variables of X2 to plot
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Solve Sæbø
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    simdata <- lpls.sim()
    X1 <- simdata$X1
    X2 <- simdata$X2
    X3 <- simdata$X3
    
    #To run exo-LPLS with orthogonal scores:
    fit.exo <- lpls(X1,X2,t(X3), npc=2, type="exo_ort")
    #Correlation plots:
    plot(fit.exo)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }
\keyword{ lpls }% __ONLY ONE__ keyword per line
