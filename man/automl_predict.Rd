\name{automl_predict}
\alias{automl_predict}
\title{automl_predict}
\description{
Predictions function, to apply a trained model on datas
}
\usage{
automl_predict(model, X, layoutputnum)
}
\arguments{
\item{model}{ model trained previously with \link{automl_train} or \link{automl_train_manual}}

\item{X}{ inputs matrix or data.frame (containing numerical values only)}

\item{layoutputnum}{ which layer number to output (default 0 for last layer)}
}
\examples{
##REGRESSION (predict Sepal.Length given other parameters)
data(iris)
xmat <- as.matrix(cbind(iris[,2:4], as.numeric(iris$Species)))
ymat <- iris[,1]
amlmodel <- automl_train_manual(Xref = xmat, Yref = ymat,
 hpar = list(modexec = 'trainwpso', verbose = FALSE))
res <- cbind(ymat, automl_predict(model = amlmodel, X = xmat))
colnames(res) <- c('actual', 'predict')
head(res)
#
\dontrun{
##CLASSIFICATION (predict Species given other Iris parameters)
data(iris)
xmat = iris[,1:4]
lab2pred <- levels(iris$Species)
lghlab <- length(lab2pred)
iris$Species <- as.numeric(iris$Species)
ymat <- matrix(seq(from = 1, to = lghlab, by = 1), nrow(xmat),
 lghlab, byrow = TRUE)
ymat <- (ymat == as.numeric(iris$Species)) + 0
amlmodel <- automl_train_manual(Xref = xmat, Yref = ymat,
 hpar = list(modexec = 'trainwpso', verbose = FALSE))
res <- cbind(ymat, round(automl_predict(model = amlmodel, X = xmat)))
colnames(res) <- c(paste('act',lab2pred, sep = '_'),
 paste('pred',lab2pred, sep = '_'))
head(res)
}
}
