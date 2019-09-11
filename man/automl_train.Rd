\name{automl_train}
\alias{automl_train}
\title{automl_train}
\description{
The multi deep neural network automatic train function (several deep neural networks are trained with automatic hyperparameters tuning, best model is kept)\cr
This function launches the \link{automl_train_manual} function by passing it parameters
for each particle at each converging step
}
\usage{
automl_train(Xref, Yref, autopar = list(), hpar = list(), mdlref = NULL)
}
\arguments{
\item{Xref}{ inputs matrix or data.frame (containing numerical values only)
}

\item{Yref}{ target matrix or data.frame (containing numerical values only)
}

\item{autopar}{ list of parameters for hyperparameters optimization, see \link{autopar} section\cr
Not mandatory (the list is preset and all arguments are initialized with default value) but it is advisable to adjust some important arguments for performance reasons (including processing time)
}

\item{hpar}{ list of parameters and hyperparameters for Deep Neural Network, see \link{hpar} section\cr
Not mandatory (the list is preset and all arguments are initialized with default value) but it is advisable to adjust some important arguments for performance reasons (including processing time)
}

\item{mdlref}{ model trained with \link{automl_train} to start training with saved \link{hpar} and \link{autopar}
(not the model)\cr
nb: manually entered parameters above override loaded ones}
}

\examples{
\dontrun{
##REGRESSION (predict Sepal.Length given other Iris parameters)
data(iris)
xmat <- cbind(iris[,2:4], as.numeric(iris$Species))
ymat <- iris[,1]
amlmodel <- automl_train(Xref = xmat, Yref = ymat)
}
##CLASSIFICATION (predict Species given other Iris parameters)
data(iris)
xmat = iris[,1:4]
lab2pred <- levels(iris$Species)
lghlab <- length(lab2pred)
iris$Species <- as.numeric(iris$Species)
ymat <- matrix(seq(from = 1, to = lghlab, by = 1), nrow(xmat), lghlab, byrow = TRUE)
ymat <- (ymat == as.numeric(iris$Species)) + 0
#with gradient descent and random hyperparameters sets
amlmodel <- automl_train(Xref = xmat, Yref = ymat,
                          autopar = list(numiterations = 1, psopartpopsize = 1, seed = 11),
                          hpar = list(numiterations = 10))
}
