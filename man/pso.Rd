\name{pso}
\alias{pso}
\title{PSO parameters and hyperparameters}
\description{
List of parameters and hyperparameters for Particle Swarm Optimization
}
\arguments{
All PSO parameters and hyperparameters are preset with default value\cr

\item{psopartpopsize}{ number of particles in swarm (discrete value)\cr
(\sQuote{autopar} context: default value 8, which means that 8 different neural net hyperparameters sets will be tested\cr
(\sQuote{hpar} context: default value 50, which means that 50 neural net weights sets will be tested\cr
#tuning priority 1}\cr

\emph{CAUTION: you should only change the values below if you know what you are doing}

\item{psovarvalmin}{ Minimum value for particles positions (default value -10)}

\item{psovarvalmax}{ maximum value for particles positions (default value 10)}

\item{psovelocitymaxratio}{ ratio applied to limit velocities (continuous value between 0 and 1, default value 0.2)}

\item{psoinertiadampratio}{ inertia damp ratio (continuous value between 0 and 1, default value 1 equivalent to OFF)}

\item{psokappa}{ kappa (default value 1)}

\item{psophi1}{ Phi 1 (default value 2.05)}

\item{psophi2}{ Phi 2 (default value 2.05)}\cr

\emph{back to \link{autopar}, \link{hpar}}
}

\seealso{
PSO video tutorial from Yarpiz
}
