\name{lpls}
\alias{lpls}
\title{ Function for lpls data exploration and regression}
\description{
This function allows statistical analysis of three data matrices arranged in an L-shape
to be modelled by lpls. Two versions of lpls is implemented, exo- and endo-lpls. See
references for details.
}
\usage{
lpls(X1, X2, X3, npc = 2, doublecenter = TRUE, scaledata = c(FALSE, FALSE, FALSE), 
type = c("exo"), impute = FALSE, niter = 25, subsetX1 = NULL, subsetX3 = NULL)
}

\arguments{
  \item{X1}{
A data matrix with \code{n} rows and \code{q} columns.
}
  \item{X2}{
A data matrix with \code{n} rows and \code{p} columns.
}
  \item{X3}{
A data matrix with \code{p} rows and \code{m} columns.
}
  \item{npc}{
The number of latent components to be extracted from each matrix
}
  \item{doublecenter}{
Logical, whether \code{X2} should be both row- and column-centered. If \code{FALSE} the matrix is only centered using the overall mean. 
}
  \item{scaledata}{
A logical vector with three elements indicating whether the columns of the three matrices should be standardized with the standard deviation.
}
  \item{type}{
Type of lpls. Character, either \code{"endo"}, \code{"exo"} or \code{"exo_ort"}. The latter is the exo-lpls with orthogonal scores. Prediction is only
implmented for the non-orthogonal exo-lpls.
}
  \item{impute}{Logical. Should missing values be imputed by SVD before analysis? If \code{FALSE} the missing values
  will be ignored in the NIPALS projections. If cross-validation is to be performed \code{impute=TRUE} must be used.
}
  \item{niter}{
The latent vectors are found using the NIPALS algorithm. This parameter defines the number of iterations for the NIPALS.
}
  \item{subsetX1}{
an optional vector specifying a subset of observations (rows of X1 and X2) to be used in the fitting process
}
  \item{subsetX3}{
an optional vector specifying a subset of variables (columns of X2 and rows X3) to be used in the fitting process
}
}
\details{
See chapter Three-block data modeling by endo- and exo-LPLS regression by S?b?, S., Martens, M and Martens, M. for details
on the endo- and exo-LPSL algorithms. Endo-LPLS is an inward regression of X1 and X3 on X2, whereas exo-LPSL is an outward
regression of X2 on X1 and X3. The \code{exo_ort} algorithm returns orthogonal scores and should be chosen for visual exploration
in correlation loading plots. If exo-LPSL with prediction is the main purpose of the model then the non-orthogonal \code{exo} 
type LPSL should be chosen for which the predict function has prediction implemented.
}
\value{
 \item{call }{The function call}
 \item{ncomp}{The number of components extracted}
 \item{coefficients}{Regression coefficient matrices: \code{B1} = The matrix of regression coefficients for predicting X1 from X2 if \code{type} is \code{exo_ort}
 \code{B3} = The matrix of regression coefficients for predicting X3 from X2 if \code{type} is \code{exo_ort}.
 \code{C} = A matrix of regression coeficients for predicting X2 from X1 and X3 in case of \code{type} is \code{endo}.
 }
 \item{scores}{The latent (normalized) scorevectors as extracted from NIPALS. One set for each extracted component of the LPLS model. 
 Thus, all scores matrices have \code{npc} number of columns. Score-matrices \code{T11} and \code{T12} correspond to matrix \code{X1},
 \code{T21} and \code{T22} to matrix \code{X2}, and \code{T31} and \code{T32} to matrix \code{X3}.}
 \item{loadings}{Loading matrices: \code{P1}=Loading(s) for \code{X1} (endo-LPLS), \code{P21} and \code{P22} =Loading(s) for \code{X2} (exo-LPSL),
 \code{P3}=Loading(s) for \code{X3} (endo-LPLS), and \code{D} = Kernel loadings for \code{X2} (endo-LPLS).}
 \item{corloadings}{Correlation loadings parallel to the regular P-loadings but without dimension. \code{R1} = Correlation loadings for \code{X1},
 \code{R21} and \code{R22} = Correlation loadings for \code{X2}, and \code{R3} = Correlation loadings for \code{X3}, \code{R2rmean}=Correlation loadings for the row means of \code{X2},
 \code{R2cmean}=Correlation loadings for the col means of \code{X2}.}
 \item{means}{\code{mX1}=The column means of matrix \code{X1},\code{mX3}=The column means of matrix \code{X3}, \code{rowmX2}=The row means of matrix \code{X2},
 \code{colmX2}=The column means of matrix \code{X2},\code{grandmX2}=The overall mean of matrix \code{X2}.}
 \item{data}{The original data matrices.}
 \item{residuals}{The residual matrices after subtracting the contribution from the \code{npc} components.}
 \item{options}{\code{Doublecenter} = Logical. \code{TRUE} if doublecentering of \code{X2} was performed. \code{scaledata} = Logical vector. \code{TRUE} if 
 respective data matrix was column scaled. \code{type} = The chosen type of lpls. Either \code{"endo"}, \code{"exo"} or \code{"exo_ort"}}
 \item{vars}{Proportions of explained un-corrected sums of squares for each component for each data matrix.}
}


\references{
Sæbø, S., Martens, M. and Martens H. (2010) Three-block data modeling by endo- and
exo-LPLS regression. In Handbook of Partial Least Squares: Concepts, Methods and
Applications. Esposito Vinzi, V.; Chin, W.W.; Henseler, J.; Wang, H. (Eds.). Springer.
}

\author{
Solve Sæbø
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{plot.lpls}},\code{\link{predict.lpls}},\code{\link{lplsCV}},
}

\examples{
    
    #Simulating data
    simdata <- lpls.sim()
    X1 <- simdata$X1
    X2 <- simdata$X2
    X3 <- simdata$X3
    
    #To run exo-LPLS with orthogonal scores:
    fit.exo <- lpls(X1,X2,t(X3), npc=2, type="exo_ort")
    #Correlation plots:
    plot(fit.exo)

    #To run exo-LPLS with non-orthogonal scores:
    fit.exo <- lpls(X1,X2,t(X3), npc=2, type="exo")
    #Predict X1
    predict(fit.exo, X2new=X2,exo.direction="X1")
    #Predict X3
    predict(fit.exo, X2new=X2,exo.direction="X3")

    #To run endo-LPSL:
    fit.endo <- lpls(X1,X2,t(X3), npc=2, type="endo")

    #Correlation loadings plots
    plot(fit.endo)

    #Predict X2 from X1 and X3:
    predict(fit.endo,X1new=X1,X3new=t(X3))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ model }
\keyword{ regression }
