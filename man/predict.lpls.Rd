\name{predict.lpls}
\alias{predict.lpls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lpls prediction
}
\description{
Uses a fitted lpls-model for predictions. If the model is an endo-lpls model, the prediction is inadvert towards
X2. Otherwise, in case of exo-lpls, the prediction is extrovert in the direction (X1 or X2) as specified by the user.
}
\usage{
predict.lpls(obj, X1new = NULL, X2new = NULL, X3new = NULL, exo.direction = c("X1", "X3"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A fitted model object of class ``lpls''. This could be either an endo- or an exo-lpls model object.
}
  \item{X1new}{
New predictor values for X1 in case of endo-lpls. Also a new (or old) X3 must be specified.
}
  \item{X2new}{
New predictor values for X2 in case of exo-lpls. X2 must have the same number of columns as X2 used in model fit in case of
prediction of X1, and the same number of rows as X2 in case of X3 prediction.
}
  \item{X3new}{
New predictor values for X3 in case of endo-lpls. Also a new (or old) X1 must be specified.
}
  \item{exo.direction}{
Applicable for exo-lpls prediction. If X1 is to be predicted then \code{exo.direction = "X1"}, or if X3 is to be predicted, then
\code{exp.direction = "X3"}
}
}
\details{
If the corner matrix X2 has been double centered in the endo-LPLS-model fit, then the predictions may have a bias since the row- and/or 
column means of the new observations are not known. In such cases it may be a better option to use grand-mean centering of X2 only.
}
\value{
A list object \code{pred} with element
%%  If it is a LIST, use
  \item{pred}{A vector or a matrix of predicted values}
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
    
    #Simulating data
    simdata <- lpls.sim()
    X1 <- simdata$X1
    X2 <- simdata$X2
    X3 <- simdata$X3

    #To run exo-LPLS with non-orthogonal scores:
    fit.exo <- lpls(X1,X2,t(X3), npc=2, type="exo")
    #Predict X1
    predict(fit.exo, X2new=X2, exo.direction="X1")
    #Predict X3
    predict(fit.exo, X2new=X2, exo.direction="X3")

    #To run endo-LPSL:
    fit.endo <- lpls(X1,X2,t(X3), npc=2, type="endo")
    #Predict X2 from X1 and X3 (in this case fitted values):
    predict(fit.endo,X1new=X1,X3new=t(X3))
}
    
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ prediction }
\keyword{ lpls }% __ONLY ONE__ keyword per line
