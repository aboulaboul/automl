\name{hpar}
\alias{hpar}
\title{Deep Neural Net parameters and hyperparameters}
\description{
List of Neural Network parameters and hyperparameters to train with gradient descent or particle swarm optimization\cr
Not mandatory (the list is preset and all arguments are initialized with default value) but it is advisable to adjust some important arguments for performance reasons (including processing time)
}
\arguments{

\item{modexec}{ \sQuote{trainwgrad} (the default value) to train with gradient descent (suitable for all volume of data)\cr
\sQuote{trainwpso} to train using Particle Swarm Optimization, each particle represents a set of neural network weights (CAUTION: suitable for low volume of data, time consuming for medium to large volume of data)}\cr

\emph{Below specific arguments to \sQuote{trainwgrad} execution mode}

\item{learningrate}{ learningrate alpha (default value 0.001)\cr
#tuning priority 1}

\item{beta1}{ see below}
\item{beta2}{ \sQuote{Momentum} if beta1 different from 0 and beta2 equal 0 )\cr
\sQuote{RMSprop} if beta1 equal 0 and beta2 different from 0\cr
\sQuote{adam optimization} if beta1 different from 0 and beta2 different from 0 (default)\cr
(default value beta1 equal 0.9 and beta2 equal 0.999)\cr
#tuning priority 2}

\item{lrdecayrate}{ learning rate decay value (default value 0, no learning rate decay, 1e-6 should be a good value to start with)\cr
#tuning priority 4}

\item{chkgradevery}{ epoch interval to run gradient check function (default value 0, for debug only)}

\item{chkgradepsilon}{ epsilon value for derivative calculations and threshold test in gradient check function (default 0.0000001)}\cr

\emph{Below specific arguments to \sQuote{trainwpso} execution mode}

\item{psoxxx}{ see \link{pso} for PSO specific arguments details}

\item{costcustformul}{ custom cost formula (default \sQuote{}, no custom cost function)\cr
standard input variables: yhat (prediction), y (target actual value)\cr
custom input variables: any variable declared in hpar may be used via alias mydl (ie: hpar(list = (foo = 1.5)) will be used in custom cost formula as mydl$foo))\cr
result: J\cr
see \sQuote{automl_train_manual} example using Mean Average Percentage Error cost function\cr
nb: X and Y matrices used as input into automl_train_manual or automl_train_manual functions are transposed (features in rows and cases in columns)}\cr

\emph{Below arguments for both execution modes}

\item{numiterations}{ number of training epochs (default value 50))\cr
#tuning priority 1}

\item{seed}{ seed for reproductibility (default 4)}

\item{minibatchsize}{ mini batch size, 2 to the power 0 for stochastic gradient descent (default 2 to the power 5)
#tuning priority 3}

\item{layersshape}{ number of nodes per layer, each nodes number initialize a hidden layer\cr
output layer nodes number, may be left to 0 it will be automatically set by Y matrix shape\cr
default value one hidden layer with 10 nodes: c(10, 0)\cr
#tuning priority 4}

\item{layersacttype}{ activation function for each layer; \sQuote{linear} for no activation or \sQuote{sigmoid}, \sQuote{relu} or \sQuote{reluleaky} or \sQuote{tanh} or \sQuote{softmax} (softmax for output layer only supported in trainwpso exec mode)\cr
output layer activation function may be left to \sQuote{}, default value \sQuote{linear} for regression, \sQuote{sigmoid} for classification\cr
nb: layersacttype parameter vector must have same length as layersshape parameter vector\cr
default value c(\sQuote{relu}, \sQuote{})\cr
#tuning priority 4}

\item{layersdropoprob}{ drop out probability for each layer, continuous value from 0 to less than 1 (give the percentage of matrix weight values to drop out randomly)\cr
nb: layersdropoprob parameter vector must have same length as layersshape parameter vector\cr
default value no drop out: c(0, 0)\cr
#tuning priority for regularization}

\item{printcostevery}{ epoch interval to test and print costs (train and cross validation cost: default value 10, for 1 test every 10 epochs)}

\item{testcvsize}{ size of cross validation sample, 0 for no cross validation sample (default 10, for 10 percent)}

\item{testgainunder}{ threshold to stop the training if the gain between last train or cross validation cost is smaller than the threshold, 0 for no stop test (default 0.000001)}

\item{costtype}{ cost type function name \sQuote{mse} or \sQuote{crossentropy} or \sQuote{custom}\cr
\sQuote{mse} for Mean Squared Error, set automatically for continuous target type\cr
\sQuote{crossentropy} set automatically for binary target type\cr
\sQuote{custom} set automatically if \sQuote{costcustformul} different from \sQuote{}
}

\item{lambda}{ regularization term added to cost function (default value 0, no regularization)
}

\item{batchnor_mom}{ batch normalization momentum for j and B (default value 0.9, no batch normalization if equal to 0)
}

\item{epsil}{ epsilon the low value to avoid dividing by 0 or log(0) in cost function, etc ... (default value 1e-12)\cr
}

\item{verbose}{ to display or not the costs and the shapes (default TRUE)}\cr

\emph{back to  \link{automl_train}, \link{automl_train_manual}}
}

\seealso{
Deep Learning specialization from Andrew NG on Coursera
}
