\name{autopar}
\alias{autopar}
\title{parameters for automatic hyperparameters optimization}
\description{
List of parameters to allow multi deep neural network automatic hyperparameters tuning with Particle Swarm Optimization\cr
Not mandatory (the list is preset and all arguments are initialized with default value) but it is advisable to adjust some important arguments for performance reasons (including processing time)
}
\arguments{
\item{psopartpopsize}{ number of particles in swarm, the main argument that should be tuned (default value 8, which is quite low)\cr
#tuning priority 1}

\item{psoxxx}{ see \link{pso} for other PSO specific arguments details}

\item{numiterations}{ number of convergence steps between particles (hyperparameters), default value 3)\cr
#tuning priority 1}

\item{auto_modexec}{ if \sQuote{TRUE} the type of Neural Net optimization will be randomly choosen between \sQuote{trainwgrad} and \sQuote{trainwpso} for each particle\cr
default value is \sQuote{FALSE} (so default value of argument \sQuote{modexec} in \link{automl_train_manual} function, actually \sQuote{trainwgrad} as default is more suited to large data volume)\cr
the value can be forced if defined in \link{hpar} list}

\item{auto_runtype}{ if \sQuote{2steps} the 2 following steps will be run automatically (default value is \sQuote{normal}):\cr
1st overfitting, the goal is performance\cr
2nd regularization, the goal is generalization \cr}

\item{auto_minibatchsize}{ see below}
\item{auto_minibatchsize_min}{ see below}
\item{auto_minibatchsize_max}{ \sQuote{auto_minibatch} default value \sQuote{TRUE} for automatic adjustment of \sQuote{minibatchsize} argument in \link{automl_train_manual} function\cr
the minimum and maximum value for \sQuote{minibatchsize} corespond to 2 to the power value (default 0 for \sQuote{auto_minibatchsize_min} and 9 for \sQuote{auto_minibatchsize_max})}

\item{auto_learningrate}{ see below}
\item{auto_learningrate_min}{ see below}
\item{auto_learningrate_max}{ \sQuote{auto_learningrate} default value \sQuote{TRUE} for automatic adjustment of \sQuote{learningrate} argument in \link{automl_train_manual} function\cr
the minimum and maximum value for \sQuote{learningrate} correspond to 10 to the power negative value (default -5 for \sQuote{auto_learningrate_min} and -2 for \sQuote{auto_learningrate_max})}

\item{auto_beta1}{ see below}
\item{auto_beta2}{ \sQuote{auto_beta1} and \sQuote{auto_beta2} default value \sQuote{TRUE} for automatic adjustment of \sQuote{beta1} and \sQuote{beta2} argument in \link{automl_train_manual} function}

\item{auto_psopartpopsize}{ see below}
\item{auto_psopartpopsize_min}{ see below}
\item{auto_psopartpopsize_max}{ \sQuote{auto_psopartpopsize} default value \sQuote{TRUE} for automatic adjustment of \sQuote{psopartpopsize} argument in \link{automl_train_manual} function (concern only \sQuote{modexec} set to \sQuote{trainwpso})\cr
the minimum and maximum value for \sQuote{psopartpopsize} ; default 2 for \sQuote{auto_psopartpopsize_min} and 50 for \sQuote{auto_psopartpopsize_max})}

\item{auto_lambda}{ see below}
\item{auto_lambda_min}{ see below}
\item{auto_lambda_max}{ \sQuote{auto_lambda} default value \sQuote{FALSE} for automatic adjustment of \sQuote{lambda} regularization argument in \link{automl_train_manual} function\cr
the minimum and maximum value for \sQuote{lambda} correspond to 10 to the power value  (default -2) for \sQuote{auto_lambda_min} and (default 4) for \sQuote{auto_lambda_max})}

\item{auto_psovelocitymaxratio}{ see below}
\item{auto_psovelocitymaxratio_min}{ see below}
\item{auto_psovelocitymaxratio_max}{ \sQuote{auto_psovelocitymaxratio} default value \sQuote{TRUE} for automatic adjustment of \sQuote{psovelocitymaxratio} PSO velocity max ratio argument in \link{automl_train_manual} function\cr
the minimum and maximum value for \sQuote{psovelocitymaxratio}; default 0.01 for \sQuote{auto_psovelocitymaxratio_min} and 0.5 for \sQuote{auto_psovelocitymaxratio_max}}

\item{auto_layers}{ see below (\sQuote{auto_layers} default value \sQuote{TRUE} for automatic adjustment of layers shape in \link{automl_train_manual} function)}
\item{auto_layers_min}{ (linked to \sQuote{auto_layers} above, set \link{hpar} \sQuote{layersshape} and \sQuote{layersacttype}) the minimum number of hidden layers (default 1 no hidden layer)}
\item{auto_layers_max}{ (linked to \sQuote{auto_layers} above, set \link{hpar} \sQuote{layersshape} and \sQuote{layersacttype}) the maximum number of hidden layers (default 2)}
\item{auto_layersnodes_min}{ (linked to \sQuote{auto_layers} above, set \link{hpar} \sQuote{layersshape} and \sQuote{layersacttype}) the minimum number of nodes per layer (default 3)}
\item{auto_layersnodes_max}{ (linked to \sQuote{auto_layers} above, set \link{hpar} \sQuote{layersshape} and \sQuote{layersacttype}) the maximum number of nodes per layer (default 33)}

\item{auto_layersdropo}{ see below}
\item{auto_layersdropoprob_min}{ see below}
\item{auto_layersdropoprob_max}{ \sQuote{auto_layersdropo} default value \sQuote{FALSE} for automatic adjustment of \link{hpar} \sQuote{layersdropoprob} in \link{automl_train_manual} function)\cr
the minimum and maximum value for \sQuote{layersdropoprob}; default 0.05 for \sQuote{auto_layersdropoprob_min} and 0.75 for \sQuote{auto_layersdropoprob_max}}

\item{seed}{ seed for reproductibility (default 4)}

\item{nbcores}{ number of cores used to parallelize particles optimization, not available on Windows (default 1, automatically reduced if not enough cores)}

\item{verbose}{ to display or not the costs at each iteration for each particle (default TRUE)}\cr

\item{subtimelimit}{ time limit in seconds for sub modelizations to avoid waiting too long for a specific particle to finish its modelization (default 3600)}\cr

\emph{back to \link{automl_train}}
}
