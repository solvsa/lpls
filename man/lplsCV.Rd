\name{lplsCV}
\alias{lplsCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for Cross-Validation in exo or endo-lpls.
}
\description{
The cross-validation is based on a fitted lpls object of type exo or endo (not exo_ort).
The predictive ability of the model is measured in terms of root mean sums of squares of prediction (RMSEP)
over the set of components from 1 to \code{npc} as defined in the call to the lpls-object.
Cross-validation is only implmented across segments defined as rows of \code{X1} and \code{X2} (horizontal prediction)
or across columns of \code{X2} (rows of \code{X3}) (vertical prediction). Cross-validation requires that all missing values have
been imputed in the model fit, that is, option \code{impute=TRUE} must be used in the call to \code{lpls}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lplsCV(fit, segments1 = NULL, segments2 = NULL, trace=TRUE)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{fit}{
A fitted lpls object of type endo or exo.
}
  \item{segments1}{
A list of cross-validation segments for horizontal prediction. For leave-one-out-CV across \code{n} objects 
use \code{segments1=as.list(1:n)}.
}
  \item{segments2}{
A list of cross-validation segments for vertical prediction. For leave-one-out-CV across \code{p} variables 
use \code{segments1=as.list(1:p)}.
}
  \item{trace}{
Logical. If \code{TRUE} the iteration number of the cross-validation loop is printed.
}
}
\details{
If no segment list is supplied a horizontal leave-one-out CV is performed across the components \code{1:npc}.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{rmsep}{The \code{rmsep} values across the components \code{1:npc}}
  \item{pred}{An array of predicted values. The last dimension of pred is the number of components used in the prediction.}
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
\code{\link{lpls}}
}
\examples{

    simdata <- lpls.sim()
    X1 <- simdata$X1
    X2 <- simdata$X2
    X3 <- simdata$X3

    #To run endo-LPSL:
    fit.endo <- lpls(X1,X2,t(X3), npc=2, type="endo")
    
    #To cross-validate horizontally
    cv <- lplsCV(fit.endo, segments1=as.list(1:dim(X1)[1]))
    #To cross-validate vertically
    cv <- lplsCV(fit.endo, segments2=as.list(1:dim(X2)[2]))
    #Three-fold CV, horizontal
    cv <- lplsCV(fit.endo, segments1=as.list(1:10, 11:20, 21:30))
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ validation }
\keyword{ lpls }% __ONLY ONE__ keyword per line
