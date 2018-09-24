\name{automl_train_manual}
\alias{automl_train_manual}
\title{automl_train_manual}
\description{
The base deep neural network train function (one deep neural network trained without automatic hyperparameters tuning)
}
\usage{
automl_train_manual(Xref, Yref, hpar = list())
}
\arguments{
\item{Xref}{ inputs matrix or data.frame (containing numerical values only)}

\item{Yref}{ target matrix or data.frame (containing numerical values only)}

\item{hpar}{ list of parameters and hyperparameters for Deep Neural Network, see \link{hpar} section\cr
Not mandatory (the list is preset and all arguments are initialized with default value) but it is advisable to adjust some important arguments for performance reasons (including processing time)}
}

\examples{
##REGRESSION (predict Sepal.Length given other Iris parameters)
data(iris)
xmat <- cbind(iris[,2:4], as.numeric(iris$Species))
ymat <- iris[,1]
#with gradient descent
amlmodel <- automl_train_manual(Xref = xmat, Yref = ymat,
                                hpar = list(learningrate = 0.01,
                                            minibatchsize = 2^2))
\dontrun{
#with PSO
amlmodel <- automl_train_manual(Xref = xmat, Yref = ymat,
                                hpar = list(modexec = 'trainwpso',
                                            numiterations = 20,
                                            psopartpopsize = 30))
#with PSO and custom cost function
f <- 'J=abs((y-yhat)/y)'
f <- c(f, 'J=sum(J[!is.infinite(J)],na.rm=TRUE)')
f <- c(f, 'J=(J/length(y))')
f <- paste(f, collapse = ';')
amlmodel <- automl_train_manual(Xref = xmat, Yref = ymat,
                                hpar = list(modexec = 'trainwpso',
                                            numiterations = 20,
                                            psopartpopsize = 30,
                                            costcustformul = f))

##CLASSIFICATION (predict Species given other Iris parameters)
data(iris)
xmat = iris[,1:4]
lab2pred <- levels(iris$Species)
lghlab <- length(lab2pred)
iris$Species <- as.numeric(iris$Species)
ymat <- matrix(seq(from = 1, to = lghlab, by = 1), nrow(xmat), lghlab, byrow = TRUE)
ymat <- (ymat == as.numeric(iris$Species)) + 0
#with gradient descent and 2 hidden layers
amlmodel <- automl_train_manual(Xref = xmat, Yref = ymat,
                                hpar = list(layersshape = c(10, 10, 0),
                                            layersacttype = c('tanh', 'relu', 'sigmoid'),
                                            layersdropoprob = c(0, 0, 0)))
#with gradient descent and no hidden layer (logistic regression)
amlmodel <- automl_train_manual(Xref = xmat, Yref = ymat,
                                hpar = list(layersshape = c(0),
                                            layersacttype = c('sigmoid'),
                                            layersdropoprob = c(0)))
#with PSO and softmax
amlmodel <- automl_train_manual(Xref = xmat, Yref = ymat,
                                hpar = list(modexec = 'trainwpso',
                                            layersshape = c(10, 0),
                                            layersacttype = c('relu', 'softmax'),
                                            layersdropoprob = c(0, 0),
                                            numiterations = 50,
                                            psopartpopsize = 50))
}
}
