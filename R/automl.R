sbamlmodcp <- function(mymodelref, mymodel)
{
 mymodel$error <- mymodelref$error
 mymodel$hpar$layersshape <- mymodelref$hpar$layersshape
 mymodel$hpar$layersacttype <- mymodelref$hpar$layersacttype
 for (l in 1:mymodelref$hpar$nblayers)
 {
  mymodel[[paste("mydl_W", l, sep = '')]] <- mymodelref[[paste("mydl_W", l, sep = '')]]
  mymodel[[paste("mydl_vdW", l, sep = '')]] <- mymodelref[[paste("mydl_vdW", l, sep = '')]]
  mymodel[[paste("mydl_sdW", l, sep = '')]] <- mymodelref[[paste("mydl_sdW", l, sep = '')]]
  mymodel[[paste("mydl_bB", l, sep = '')]] <- mymodelref[[paste("mydl_bB", l, sep = '')]]
  mymodel[[paste("mydl_vdbB", l, sep = '')]] <- mymodelref[[paste("mydl_vdbB", l, sep = '')]]
  mymodel[[paste("mydl_sdbB", l, sep = '')]] <- mymodelref[[paste("mydl_sdbB", l, sep = '')]]
  mymodel[[paste("mydl_j", l, sep = '')]] <- mymodelref[[paste("mydl_j", l, sep = '')]]
  mymodel[[paste("mydl_vdj", l, sep = '')]] <- mymodelref[[paste("mydl_vdj", l, sep = '')]]
  mymodel[[paste("mydl_sdj", l, sep = '')]] <- mymodelref[[paste("mydl_sdj", l, sep = '')]]
  mymodel[[paste("mydl_M", l, sep = '')]] <- mymodelref[[paste("mydl_M", l, sep = '')]]
  mymodel[[paste("mydl_vM", l, sep = '')]] <- mymodelref[[paste("mydl_vM", l, sep = '')]]
  mymodel[[paste("mydl_V", l, sep = '')]] <- mymodelref[[paste("mydl_V", l, sep = '')]]
  mymodel[[paste("mydl_vV", l, sep = '')]] <- mymodelref[[paste("mydl_vV", l, sep = '')]]
 }
 return(mymodel)
}

sbamlmodhparcp <- function(mdlref, myhpartoignore, autopar = NULL, hpar)
{
 valmodimpflagok <- 1
 if (is.null(autopar))
 {
  if (any(c("layersshape", "layersacttype") %in% names(hpar)))
  {
   mylastlog <- ("wrng: layersshape or layersacttype defined (mdlref won't be used)")
   warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
   valmodimpflagok <- 0
  }
  if (valmodimpflagok == 1 & "layersdropoprob" %in% names(hpar))
  {
   if (length(hpar$layersdropoprob) != length(mdlref$hpar$layersdropoprob))
   {
    mylastlog <- ("wrng: layersdropoprob defined with incompatible shape (mdlref won't be used)")
    warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
    valmodimpflagok <- 0
   }
  }
  if (valmodimpflagok == 1 & "modexec" %in% names(hpar))
  {
   if (hpar$modexec == "trainwpso" & mdlref$hpar$modexec == "trainwgrad")
   {
    myhpartoignore = c(myhpartoignore, 'psopartpopsize', 'psonbvar2optim', 'psonumiteration',
                       'psovarvalmin', 'psovarvalmax', 'psovelocitymaxratio',
                       'psoinertiadampratio', 'psokappa', 'psophi1', 'psophi2', 'psoseed')
   } else if (hpar$modexec == "trainwgrad" & mdlref$hpar$modexec == "trainwpso")
   {
    myhpartoignore = c(myhpartoignore, 'learningrate', 'lrdecayrate', 'chkgradevery',
                       'chkgradepsilon', 'beta1', 'beta2')
    if (!"minibatchsize" %in% names(hpar)) {myhpartoignore = c(myhpartoignore, "minibatchsize")}
   }
  }
  if (valmodimpflagok == 1 & all(c('autopar', 'mylspsoparamlvauto') %in% names(hpar)))
  {
   myhpartoignore = c(myhpartoignore, c('useautopar', 'overfitautopar', 'verbose'))
  }
 }
 # forcing hpar with mdlref$hpar in agreement with hpar and myhpartoignore
 if (valmodimpflagok == 1)
 {
  for (i in 1:length(mdlref$hpar))
  {
   if (!names(mdlref$hpar)[i] %in% myhpartoignore &
       !names(mdlref$hpar)[i] %in% names(hpar))
   {
    hpar[[names(mdlref$hpar)[i]]] <- mdlref$hpar[[names(mdlref$hpar)[i]]]
   }
  }
 }
 # autopar validated before but may be changed
 if (!is.null(autopar))
 {
  if (autopar$auto_modexec == FALSE) {hpar$modexec <- mdlref$hpar$modexec}
  if (autopar$auto_minibatchsize == FALSE) {hpar$minibatchsize <- mdlref$hpar$minibatchsize}
  if (autopar$auto_learningrate == FALSE) {hpar$learningrate <- mdlref$hpar$learningrate}
  if (autopar$auto_beta1 == FALSE) {hpar$beta1 <- mdlref$hpar$beta1}
  if (autopar$auto_beta2 == FALSE) {hpar$beta2 <- mdlref$hpar$beta2}
  if (autopar$auto_psopartpopsize == FALSE) {hpar$psopartpopsize <- mdlref$hpar$psopartpopsize}
  if (autopar$auto_lambda == FALSE) {hpar$lambda <- mdlref$hpar$lambda}
  if (autopar$auto_psovelocitymaxratio == FALSE) {hpar$psovelocitymaxratio <- mdlref$hpar$psovelocitymaxratio}
  if (autopar$auto_layers == FALSE)
  {
   if (!"layersshape" %in% names(hpar)) {hpar$layersshape <- mdlref$hpar$layersshape}
   if (!"layersacttype" %in% names(hpar)) {hpar$layersacttype <- mdlref$hpar$layersacttype}
   #if layersdropoprob not in hpar, will be defined
  }
  if (any(c("layersshape", "layersacttype", "layersdropoprob") %in% names(hpar)))
  {
   autopar$auto_layers_min <- autopar$auto_layers_max <- mdlref$hpar$nblayers
  }
 }
 return(list(autopar = autopar, hpar = hpar, valmodimpflagok = valmodimpflagok))
}

sbamlmodpsoparam <- function(mylspsoparam, mylspsopart)
{
  mylspsoparam$phi = mylspsoparam$phi1 + mylspsoparam$phi2
  mylspsoparam$chi = 2 * mylspsoparam$kappa / abs(2 - mylspsoparam$phi - sqrt(mylspsoparam$phi^2 - 4 * mylspsoparam$phi))
  mylspsoparam$w = mylspsoparam$chi;# inertia coef (no damping ratio with constriction coefs)
  mylspsoparam$c1 = mylspsoparam$chi * mylspsoparam$phi1;# personal acceleration coef for cognitive composant
  mylspsoparam$c2 = mylspsoparam$chi * mylspsoparam$phi2;# social acceleration coef for social composant
  mylspsoparam$maxvel = (mylspsoparam$varvalmax - mylspsoparam$varvalmin) * mylspsoparam$velocitymaxratio
  mylspsoparam$minvel = -mylspsoparam$maxvel
  mylspsoparam$activpart = 0

  if (is.null(mylspsoparam$modeinit))
  {
   mylastlog <- ("error: no model trained")
   stop(mylastlog);
  } else {
   if (mylspsoparam$modeinit %in% c('angdlpso', 'autopar'))
   {
    mylspsopart <- list(
     Position = NULL,
     Velocity = 0,
     Cost = NA,
     Best = list(
      Position = NULL,
      Cost = NA,
      idpart = NA
     )
    )
    mylspsopart <- rep(list(mylspsopart), mylspsoparam$partpopsize)

    #pop members initialization
    for (i in 1:mylspsoparam$partpopsize)
    {
     if (mylspsoparam$modeinit == 'angdlpso')
     {
      #random solution
      mydl <- mylspsoparam$angdl$mydl
      mydl <- sbamlmoddlparaminit(mydl, mydl$hpar$layersshape, mydl$hpar$layersacttype,
                                  mylspsoparam$seed + i, mydl$hpar$batchnor_mom)
      for (l in 1:mydl$hpar$nblayers)
      {
       # n = mydl$hpar$layersshape[l + 1]
       # m = mydl$hpar$layersshape[l]
       mylspsopart[[i]]$Position <- c(mylspsopart[[i]]$Position, c(mydl[[paste("mydl_W", l, sep = '')]]))
       mylspsopart[[i]]$Position <- c(mylspsopart[[i]]$Position, c(mydl[[paste("mydl_bB", l, sep = '')]]))
       if (mydl$hpar$batchnor_mom != 0)
       {
        mylspsopart[[i]]$Position <- c(mylspsopart[[i]]$Position, c(mydl[[paste("mydl_j", l, sep = '')]]))
       }
       #to avoid stagnation for bB & j
       myprov <- length(mylspsopart[[i]]$Position[mylspsopart[[i]]$Position == 0])
       set.seed(mylspsoparam$seed + i)
       mylspsopart[[i]]$Position[mylspsopart[[i]]$Position == 0] <- stats::runif(n = myprov, min = 0, max = 0.05)
       myprov <- length(mylspsopart[[i]]$Position[mylspsopart[[i]]$Position == 1])
       set.seed(mylspsoparam$seed + i)
       mylspsopart[[i]]$Position[mylspsopart[[i]]$Position == 1] <- stats::runif(n = myprov, min = 0.95, max = 1)
       rm(myprov)
      }
      mylspsoparam$nbvar2optim <- length(mylspsopart[[i]]$Position)
     } else if (mylspsoparam$modeinit == 'autopar')
     {
      #random solution
      set.seed(mylspsoparam$seed + i)
      mylspsopart[[i]]$Position <- stats::runif(n = mylspsoparam$nbvar2optim,
                                                min = mylspsoparam$varvalmin,
                                                max = mylspsoparam$varvalmax)
     }

     if (i == 1)
     {
      #global best
      mylspsoparam$globalbest = list(
       Position = rep(0, mylspsoparam$nbvar2optim),
       Cost = Inf,
       idpart = NA); #the worst possible cost value
      #logs of best each cost operation on each iteration
      mylspsoparam$globbestcostlog = matrix(data = rep(NA, 2 * mylspsoparam$numiteration),
                                            nrow = mylspsoparam$numiteration, ncol = 2,
                                            byrow = FALSE)
      colnames(mylspsoparam$globbestcostlog) = c('idpart', 'cost')
     }
     #initialize velocity
     mylspsopart[[i]]$Velocity <- rep(0, mylspsoparam$nbvar2optim)
     #evaluation
     mylspsopart[[i]]$Cost <- Inf
     #update personal best
     mylspsopart[[i]]$Best$Position <- mylspsopart[[i]]$Position
     mylspsopart[[i]]$Best$Cost <- mylspsopart[[i]]$Cost
     mylspsopart[[i]]$Best$idpart <- i
     if (mylspsoparam$modeinit == 'angdlpso') {mylspsopart[[i]]$angdl$mydl <- mydl}
    }
   }
   return(list(mylspsoparam = mylspsoparam, mylspsopart = mylspsopart))
  }
}

sbamlmodpsooptkern <- function(mylspsoparam, mylspsoactivpart)
{
  # nb: mylspsoparam$activpart must correspond 2 mylspsoactivpart (4 seed)
  #
  # pso kernel sub for active particle
  set.seed(mylspsoparam$seed + mylspsoparam$activpart)
  # update velocity & position
  mylspsoactivpart$Velocity <- mylspsoparam$w * mylspsoactivpart$Velocity +
    (mylspsoparam$c1 * stats::runif(n = mylspsoparam$nbvar2optim, min = 0, max = 1) * (mylspsoactivpart$Best$Position - mylspsoactivpart$Position)) +
    (mylspsoparam$c2 * stats::runif(n = mylspsoparam$nbvar2optim, min = 0, max = 1) * (mylspsoparam$globalbest$Position - mylspsoactivpart$Position))
  mylspsoactivpart$Position <- mylspsoactivpart$Position + mylspsoactivpart$Velocity

  # apply velocity limits
  mylspsoactivpart$Velocity[mylspsoactivpart$Velocity < mylspsoparam$minvel] <- mylspsoparam$minvel
  mylspsoactivpart$Velocity[mylspsoactivpart$Velocity > mylspsoparam$maxvel] <- mylspsoparam$maxvel

  # apply lower and upper bounds limits
  mylspsoactivpart$Position[mylspsoactivpart$Position < mylspsoparam$varvalmin] <- mylspsoparam$varvalmin
  mylspsoactivpart$Position[mylspsoactivpart$Position > mylspsoparam$varvalmax] <- mylspsoparam$varvalmax
  if (mylspsoparam$modecost %in% c('angdlpso'))
  {
    # cost
    mylspsoactivpart <- sbamlmodpsocost(mylspsoparam = mylspsoparam, mylspsoactivpart = mylspsoactivpart)
    if (mylspsoactivpart$Cost < mylspsoactivpart$Best$Cost)
    {
      mylspsoactivpart$Best$Position <- mylspsoactivpart$Position
      mylspsoactivpart$Best$Cost <- mylspsoactivpart$Cost
      mylspsoactivpart$Best$idpart <- mylspsoparam$activpart
    }
    if (mylspsoactivpart$Cost < mylspsoparam$globalbest$Cost)
    {
      mylspsoparam$globalbest <- mylspsoactivpart$Best
      mylspsoparam$angdl$mydl <- mylspsoactivpart$angdl$mydl
    }
  }
  return(list(mylspsoparam = mylspsoparam, mylspsoactivpart = mylspsoactivpart))
}

sbamlmodpsocost <- function(mylspsoparam, mylspsoactivpart)
{
  if (mylspsoparam$modecost == 'angdlpso')
  {
    mydl <- sbamlpart2mydl(mylspsoactivpart = mylspsoactivpart)
    myprovres <- sbamlmoddlfw(mydl, X = mylspsoparam$angdl$X,
                              modexec = mydl$hpar$modexec,
                              nblayers = mydl$hpar$nblayers,
                              layersacttype = mydl$hpar$layersacttype,
                              layersdropoprob = mydl$hpar$layersdropoprob,
                              seed = mydl$hpar$seed,
                              batchnor_mom = mydl$hpar$batchnor_mom,
                              epsil = mydl$hpar$epsil)
    Yhat <- myprovres[["A"]]; mydl <- myprovres[["mydl"]]; rm(myprovres)
    mycost <- sbamlmoddlcost(mydl, y = mylspsoparam$angdl$Y,
                             yhat = Yhat,
                             costtype = mydl$hpar$costtype,
                             lambda = mydl$hpar$lambda,
                             nblayers = mydl$hpar$nblayers,
                             epsil = mydl$hpar$epsil)
    mylspsoactivpart$Cost <- mycost
    mylspsoactivpart$angdl$mydl <- mydl
    return(mylspsoactivpart)
  }
}

sbamlmoddlcost <- function(mydl, y, yhat, costtype, lambda, nblayers, epsil)
{
  if (costtype == 'crossentropy')
  {
    m = dim(y)[2]
    J <- y * log(yhat + epsil) + (1 - y) * log(1 - yhat + epsil)
    J[J > 0] <- 0;# if yhat + epsil > 1
    J <- (-1 / m) * sum(J)
  } else if (costtype == 'mse')
  {
    m = dim(y)[2]
    J <- (1 / (2 * m)) * sum((y - yhat) ^ 2)
  } else if (costtype == 'custom')
  {
    mytest <- try(eval(parse(text = mydl$hpar$costcustformul)), TRUE)
    if (grepl('error',class(mytest)[1])) {J = Inf}
  }
  if (lambda != 0)
  {
    m = dim(y)[2]
    l2regul <- 0
    for (l in 1:nblayers)
    {
      W <- mydl[[paste("mydl_W", l, sep = '')]]
      l2regul <- l2regul + (lambda / (2 * m)) * sum(W ^ 2)
    }
    J <- J + l2regul
  }
  return(J)
}

sbamlmoddlcostbk <- function(y, yhat, costtype)
{
  if (costtype == 'crossentropy')
  {
    dJ <- sbamlmatopebroadcst(sbamlmatopebroadcst(y, yhat, 'divid') * (-1),
                              sbamlmatopebroadcst((1 - y), (1 - yhat), 'divid'),
                              'add')
  } else if (costtype == 'mse')
  {
    dJ <- sbamlmatopebroadcst(yhat, y  * (-1), 'add')
  }
  return(dJ)
}

sbamlmoddldropo <- function(mymatA, mymatdD, mydropo)
{
  #inverted dropout
  mymin <- min(mymatA)
  mymax <- max(mymatA)
  mymatA <- sbamlmatopebroadcst(mymatA, mymatdD, 'mult') / (1 - mydropo)
  mymatA[mymatA > mymax] <- mymax
  mymatA[mymatA < mymin] <- mymin
  #why not to avoid result > 1 in cost function ...
  return(mymatA)
}

sbamlmoddlgradchk <- function(mydl, X, Y, nblayers, chkgradepsilon, epsil)
{
  mybigthetasize <- 0
  myvectrefpos <- myvectrefposname <- mygradapprox <- mygrad <- c()
  for (l in 1:nblayers)
  {
    myprov <- dim(mydl[[paste("mydl_W", l, sep = '')]])
    myvectrefpos <- c(myvectrefpos, mybigthetasize + myprov[1] * myprov[2])
    mybigthetasize <- mybigthetasize + myprov[1] * myprov[2]
    myvectrefposname <- c(myvectrefposname, paste('W',l,sep = ''))
    myprov <- dim(mydl[[paste("mydl_bB", l, sep = '')]])
    myvectrefpos <- c(myvectrefpos, mybigthetasize + myprov[1] * myprov[2])
    mybigthetasize <- mybigthetasize + myprov[1] * myprov[2]
    myvectrefposname <- c(myvectrefposname, paste('bB',l,sep = ''))
    if (mydl$hpar$batchnor_mom != 0)
    {
      myprov <- dim(mydl[[paste("mydl_j", l, sep = '')]])
      myvectrefpos <- c(myvectrefpos, mybigthetasize + myprov[1] * myprov[2])
      mybigthetasize <- mybigthetasize + myprov[1] * myprov[2]
      myvectrefposname <- c(myvectrefposname, paste('j',l,sep = ''))
      mygrad <- c(mygrad,
                  as.vector(mydl[[paste("mydl_dW", l, sep = '')]]),
                  as.vector(mydl[[paste("mydl_dbB", l, sep = '')]]),
                  as.vector(mydl[[paste("mydl_dj", l, sep = '')]]))
    } else {
      mygrad <- c(mygrad,
                  as.vector(mydl[[paste("mydl_dW", l, sep = '')]]),
                  as.vector(mydl[[paste("mydl_dbB", l, sep = '')]]))
    }
  }
  for (i in 1:mybigthetasize)
  {
    myrefpos <- 1; myflagcontinue <- 1
    while (myflagcontinue == 1)
    {
      if (i <= myvectrefpos[myrefpos])
      {
        myflagcontinue <- 0
      } else {
        myrefpos <- myrefpos + 1
      }
    }
    mymat <- mydl[[paste("mydl_",myvectrefposname[myrefpos] , sep = '')]]
    if (myrefpos == 1)
    {
      myposinmatrix <- i
    } else {
      myposinmatrix <- i - (myvectrefpos[myrefpos - 1])
    }
    myvalref <- mymat[myposinmatrix]
    mymat[myposinmatrix] <- myvalref + chkgradepsilon
    mydl[[paste("mydl_",myvectrefposname[myrefpos] , sep = '')]] <- mymat
    myprovres <- sbamlmoddlfw(mydl, X,
                              mydl$hpar$modexec,
                              nblayers,
                              layersacttype = mydl$hpar$layersacttype,
                              layersdropoprob = mydl$hpar$layersdropoprob,
                              seed = mydl$hpar$seed,
                              batchnor_mom = mydl$hpar$batchnor_mom, epsil = mydl$hpar$epsil)
    Yhat <- myprovres[["A"]]; mydl <- myprovres[["mydl"]]; rm(myprovres)
    Jplusepsilon <- sbamlmoddlcost(mydl, Y, Yhat,
                                   costtype = mydl$hpar$costtype,
                                   lambda = mydl$hpar$lambda,
                                   nblayers = nblayers,
                                   epsil = epsil)
    mymat[myposinmatrix] <- myvalref - chkgradepsilon
    mydl[[paste("mydl_",myvectrefposname[myrefpos] , sep = '')]] <- mymat
    myprovres <- sbamlmoddlfw(mydl, X,
                              mydl$hpar$modexec,
                              nblayers,
                              layersacttype = mydl$hpar$layersacttype,
                              layersdropoprob = mydl$hpar$layersdropoprob,
                              seed = mydl$hpar$seed,
                              batchnor_mom = mydl$hpar$batchnor_mom, epsil = mydl$hpar$epsil)
    Yhat <- myprovres[["A"]]; mydl <- myprovres[["mydl"]]; rm(myprovres)
    Jminusepsilon = sbamlmoddlcost(mydl, Y, Yhat,
                                   costtype = mydl$hpar$costtype,
                                   lambda = mydl$hpar$lambda,
                                   nblayers = nblayers,
                                   epsil = epsil)
    mygradapprox[i] = (Jplusepsilon - Jminusepsilon) / (2 * chkgradepsilon)
    mymat[myposinmatrix] <- myvalref
    mydl[[paste("mydl_",myvectrefposname[myrefpos] , sep = '')]] <- mymat
  }
  numerator = sbamlllnvectnorm((mygradapprox - mygrad), 2)
  denominator = sbamlllnvectnorm(mygradapprox, 2) + sbamlllnvectnorm(mygrad, 2)
  difference = numerator / denominator
  if (difference > chkgradepsilon)
  {
    cat(paste("Backward propagation to monitor! Difference =", difference, '\n'))
    cat(paste('Gradapprox:',
              paste(utils::head(mygradapprox), collapse = ', '), ' ... ', paste(utils::tail(mygradapprox), collapse = ', '),
              'Grad:',
              paste(utils::head(mygrad), collapse = ', '), ' ... ', paste(utils::tail(mygrad), collapse = ', '), '\n', sep = '\n'))
  } else {
    cat(paste("Backward propagation works fine! Difference =", difference, '\n'))
  }
}

automl_train <- function(Xref, Yref, autopar = list(), hpar = list(), mdlref = NULL)
{
  myautomatablehpar <- c('modexec','minibatchsize','learningrate',
                         'beta1','beta2','psopartpopsize',
                         'lambda','psovelocitymaxratio',
                         'layersshape','layersacttype','layersdropoprob')
  sbamlpart2parlln1 <- function(autopar = list(), hpar = list(), mylspsoparam = list())
  {
   sbamlpart2parlln2 <- function(mylspsoactivpart)
   {
    myactivpartnbr <- mylspsoactivpart$Best$idpart
    mylspsoparam$activpart <- myactivpartnbr
    myprovres <- sbamlmodpsooptkern(mylspsoparam = mylspsoparam, mylspsoactivpart = mylspsoactivpart)
    mylspsoparam <- myprovres[["mylspsoparam"]]
    mylspsoactivpart <- myprovres[["mylspsoactivpart"]]
    rm(myprovres)
    hpar$verbose <- FALSE
    hpar$useautopar <- TRUE
    hpar$seedautopar <- autopar$seed;#for train/test constant sampling
    hpar$seed <- autopar$seed + myactivpartnbr
    if (autopar$auto_modexec == TRUE & !"modexec" %in% names(hpar))
    {
     set.seed(autopar$seed + myactivpartnbr)
     hpar$modexec <- c('trainwgrad', 'trainwpso')[sample(2, 1, replace = TRUE)]
    } else if (autopar$auto_modexec == FALSE & !"modexec" %in% names(hpar))
    {
     hpar$modexec <- 'trainwgrad'
    }
    mydep = 1
    if (autopar$auto_minibatchsize == TRUE & !"minibatchsize" %in% names(hpar))
    {
     if (hpar$modexec != 'trainwpso')
     {
      hpar$minibatchsize <- eval(parse(text =
                                        paste0("2^round(",
                                               sbamlprepscalem11201(mylspsoactivpart$Position[mydep]),
                                               " * (",
                                               autopar$auto_minibatchsize_max,
                                               " - (",
                                               autopar$auto_minibatchsize_min,
                                               ")) + (",
                                               autopar$auto_minibatchsize_min,
                                               "), digits = 0)")
      ))
     }
    }
    mydep = mydep + 1
    if (autopar$auto_learningrate == TRUE & !"learningrate" %in% names(hpar))
    {
     hpar$learningrate <- eval(parse(text =
                                      paste0("10^(",
                                             sbamlprepscalem11201(mylspsoactivpart$Position[mydep]),
                                             " * (",
                                             autopar$auto_learningrate_max,
                                             " - (",
                                             autopar$auto_learningrate_min,
                                             ")) + (",
                                             autopar$auto_learningrate_min,
                                             "))")
     ))
    }
    mydep = mydep + 1
    if (autopar$auto_beta1 == TRUE & !"beta1" %in% names(hpar))
    {
     hpar$beta1 <- eval(parse(text =
                               paste0("1 - 10^(-",
                                      sbamlprepscalem11201(mylspsoactivpart$Position[mydep]),
                                      " - 1)")
     ))
    }
    mydep = mydep + 1
    if (autopar$auto_beta2 == TRUE & !"beta2" %in% names(hpar))
    {
     hpar$beta2 <- eval(parse(text =
                               paste0("1 - 10^(-",
                                      sbamlprepscalem11201(mylspsoactivpart$Position[mydep]),
                                      " - 2)")
     ))
    }
    mydep = mydep + 1
    if (autopar$auto_psopartpopsize == TRUE & !"psopartpopsize" %in% names(hpar))
    {
     hpar$psopartpopsize <- eval(parse(text =
                                        paste0("round(",
                                               sbamlprepscalem11201(mylspsoactivpart$Position[mydep]),
                                               " * (",
                                               autopar$auto_psopartpopsize_max,
                                               " - (",
                                               autopar$auto_psopartpopsize_min,
                                               ")) + (",
                                               autopar$auto_psopartpopsize_min,
                                               "), digits = 0)")
     ))
    }
    mydep = mydep + 1
    if (autopar$auto_lambda == TRUE & !"lambda" %in% names(hpar))
    {
     hpar$lambda <- eval(parse(text =
                                paste0("10^(",
                                       sbamlprepscalem11201(mylspsoactivpart$Position[mydep]),
                                       " * (",
                                       autopar$auto_lambda_max,
                                       " - (",
                                       autopar$auto_lambda_min,
                                       ")) + (",
                                       autopar$auto_lambda_min,
                                       "))")
     ))
    }
    mydep = mydep + 1
    if (autopar$auto_psovelocitymaxratio == TRUE & !"psovelocitymaxratio" %in% names(hpar))
    {
     hpar$psovelocitymaxratio <- eval(parse(text =
                                             paste0(sbamlprepscalem11201(mylspsoactivpart$Position[mydep]),
                                                    " * (",
                                                    autopar$auto_psovelocitymaxratio_max,
                                                    " - (",
                                                    autopar$auto_psovelocitymaxratio_min,
                                                    ")) + (",
                                                    autopar$auto_psovelocitymaxratio_min,
                                                    ")")
     ))
    }
    mydep <- mydep + 1
    if (any(c("layersshape", "layersacttype", "layersdropoprob") %in% names(hpar)))
    {
     if ("layersshape" %in% names(hpar))
     {
      mylayernbr <- length(hpar$layersshape)
     } else if ("layersacttype" %in% names(hpar))
     {
      mylayernbr <- length(hpar$layersacttype)
     } else if ("layersdropoprob" %in% names(hpar))
     {
      mylayernbr <- length(hpar$layersdropoprob)
     }
    } else {
     mylayernbr <- eval(parse(text =
                               paste0("round(",
                                      sbamlprepscalem11201(mylspsoactivpart$Position[mydep]),
                                      " * (",
                                      autopar$auto_layers_max,
                                      " - (",
                                      autopar$auto_layers_min,
                                      ")) + (",
                                      autopar$auto_layers_min,
                                      "), digits = 0)")
     ))
    }
    mydep <- mydep + 1
    if (autopar$auto_layers == TRUE & !"layersshape" %in% names(hpar))
    {
     myformula <- "c("
     for (i in mydep:(mydep + mylayernbr - 1))
     {
      if (i != (mydep + mylayernbr - 1))
      {
       myformula <- c(myformula, eval(parse(text =
                                             paste0("round(",
                                                    sbamlprepscalem11201(mylspsoactivpart$Position[i]),
                                                    " * (",
                                                    autopar$auto_layersnodes_max,
                                                    " - (",
                                                    autopar$auto_layersnodes_min,
                                                    ")) + (",
                                                    autopar$auto_layersnodes_min,
                                                    "), digits = 0)")
       )))
       myformula <- c(myformula, ',')
      } else {
       myformula <- c(myformula, '0)')
      }
     }
     hpar$layersshape <-  eval(parse(text = paste(myformula, collapse = '')))
    }
    mydep <- mydep + autopar$auto_layers_max
    if (autopar$auto_layers == TRUE & !"layersacttype" %in% names(hpar))
    {
     myformula <- "c("
     myhiddenacttype <- c('sigmoid', 'relu', 'reluleaky', 'tanh')
     for (i in mydep:(mydep + mylayernbr - 1))
     {
      if (i != (mydep + mylayernbr - 1))
      {
       myformula <- c(myformula, "'", myhiddenacttype[eval(parse(text =
                                                                  paste0("round(",
                                                                         sbamlprepscalem11201(mylspsoactivpart$Position[i]),
                                                                         " * (",
                                                                         length(myhiddenacttype),
                                                                         " - (",
                                                                         1,
                                                                         ")) + (",
                                                                         1,
                                                                         "), digits = 0)")
       ))], "'")
       myformula <- c(myformula, ',')
      } else {
       myformula <- c(myformula, "'')")
      }
     }
     hpar$layersacttype <-  eval(parse(text = paste(myformula, collapse = '')))
    }
    mydep <- mydep + autopar$auto_layers_max
    if (!"layersdropoprob" %in% names(hpar))
    {
     myformula <- "c("
     for (i in mydep:(mydep + mylayernbr - 1))
     {
      if (i != (mydep + mylayernbr - 1))
      {
       if (autopar$auto_layersdropo == FALSE)
       {
        myformula <- c(myformula, 0)
       } else {
        myformula <- c(myformula, eval(parse(text =
                                              paste0(sbamlprepscalem11201(mylspsoactivpart$Position[i]),
                                                     " * (",
                                                     autopar$auto_layersdropoprob_max,
                                                     " - (",
                                                     autopar$auto_layersdropoprob_min,
                                                     ")) + (",
                                                     autopar$auto_layersdropoprob_min,
                                                     ")")
        )))
       }
       myformula <- c(myformula, ',')
      } else {
       myformula <- c(myformula, '0)')
      }
     }
     hpar$layersdropoprob <-  eval(parse(text = paste(myformula, collapse = '')))
    }
    myprovres <- sbamlvalidhpar(hpar = hpar, myruntype = 'auto')
    hpar <- myprovres[["hpar"]]
    valhparflagok <- myprovres[["tstflagok"]]
    rm(myprovres)
    #
    if (valhparflagok == 1)
    {
     setTimeLimit(cpu = Inf, elapsed = autopar$subtimelimit, transient = FALSE)
     mymodel <- try(automl_train_manual(Xref = Xref, Yref = Yref, hpar = hpar), TRUE)
     setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    } else {
     mymodel <- NULL
    }
    if (grepl('error', class(mymodel)[1]) | is.null(mymodel))
    {
     mymodel <- list(error = list(tr = Inf, cv = Inf))
     mylspsoactivpart$Cost <- Inf
    } else
    {
     mylspsoactivpart$Cost <- sbamlmodresult(trerr = mymodel$error$tr,
                                             cverr = mymodel$error$cv,
                                             cvflag = mymodel$hpar$testcvsize != 0 & mymodel$hpar$overfitautopar == FALSE)
     if (mylspsoactivpart$Cost < mylspsoactivpart$Best$Cost)
     {
      mylspsoactivpart$Best$Position <- mylspsoactivpart$Position
      mylspsoactivpart$Best$Cost <- mylspsoactivpart$Cost
      mylspsoactivpart$Best$idpart <- myactivpartnbr
      mylspsoactivpart$Best$model <- mymodel
     }
    }
    return(mylspsoactivpart)
   }
  }
  #
  hparref <- sbamlvalidhpar(hpar = hpar, myruntype = 'rawmanual')[["hpar"]]
  autoparref <- sbamlvalidautopar(autopar = autopar, myruntype = 'rawmanual')[["autopar"]]
  if (!is.null(mdlref))
  {
   myprovres <- sbamltestmodel(mymodel = mdlref, mymodeltype = 'auto', guess1stlayer = TRUE)
   mdlref <- myprovres[["mymodel"]]
   valmodautoflagok <- myprovres[["tstflagok"]]
   rm(myprovres)
   if (valmodautoflagok == 1)
   {
    #initializing autopar with mdlref$autopar if not defined
    for (i in 1:length(names(mdlref$autopar)))
    {
     if (!names(mdlref$autopar)[i] %in% names(autopar))
     {
      autopar[[names(mdlref$autopar)[i]]] <- mdlref$autopar[[names(mdlref$autopar)[i]]]
     }
    }
    #initializing hpar with mdlref$hpar
    autopar <- sbamlvalidautopar(autopar = autopar, myruntype = 'auto')[["autopar"]]
    myprovres <- sbamlmodhparcp(mdlref, myhpartoignore = myautomatablehpar, autopar, hpar)
    autopar <- myprovres[["autopar"]]
    hpar <- myprovres[["hpar"]]
    rm(myprovres)
   }
  } else {valmodautoflagok <- 0}
  #
  Xref <- t(as.matrix(Xref))
  Yref <- t(as.matrix(Yref))
  sbamlmatstest(Xref, Yref)
  Xref <- sbamlmattest(Xref, "Xref")
  Yref <- sbamlmattest(Yref, "Yref")
  #
  if (valmodautoflagok == 0)
  {
   autopar <- sbamlvalidautopar(autopar = autoparref, myruntype = 'manual')[["autopar"]]
  }
  #
  if (autopar$auto_runtype == "2steps")
  {
   myautonbloop <- 2
  } else {
   myautonbloop <- 1
  }
  #
  for (autoloopnbr in 1:myautonbloop)
  {
   if (autopar$auto_runtype == "2steps" & autoloopnbr == 1)
   {
    hpar$overfitautopar <- TRUE;#TRUE to overfit (cv but no cv test)
    autopar[["auto_lambda"]] <- FALSE
    autopar[["auto_layersdropo"]] <- FALSE
    myprovres <- sbamlvalidautopar(autopar = autopar, myruntype = 'auto')
    autopar <- myprovres[["autopar"]]
    if (myprovres[["tstflagok"]] == 0)
    {
     mylastlog <- ("wrng: so mdlref autopar won't be used")
     warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
     autopar <- autoparref
     autopar[["auto_lambda"]] <- FALSE
     autopar[["auto_layersdropo"]] <- FALSE
     autopar <- sbamlvalidautopar(autopar = autopar, myruntype = 'manual')[["autopar"]]
    }
    rm(myprovres)

   } else if (autopar$auto_runtype == "2steps" & autoloopnbr == 2)
   {
    hpar$overfitautopar <- FALSE;#FALSE to generalize (cv w cv test)
    autopar[["auto_modexec"]] <- FALSE
    autopar[["auto_minibatchsize"]] <- FALSE
    autopar[["auto_learningrate"]] <- FALSE
    autopar[["auto_beta1"]] <- FALSE
    autopar[["auto_beta2"]] <- FALSE
    autopar[["auto_psopartpopsize"]] <- FALSE
    autopar[["auto_psovelocitymaxratio"]] <- FALSE
    autopar[["auto_layers"]] <- FALSE
    #
    autopar[["auto_lambda"]] <- TRUE
    autopar[["auto_layersdropo"]] <- TRUE
    #
    myprovres <- sbamltestmodel(mymodel = mdlref, mymodeltype = 'auto', guess1stlayer = TRUE)
    mdlref <- myprovres[["mymodel"]]
    valmodautoflagok <- myprovres[["tstflagok"]]
    rm(myprovres)
    if (valmodautoflagok == 1)
    {
     #initializing autopar with mdlref$autopar if not defined
     for (i in 1:length(names(mdlref$autopar)))
     {
      if (!names(mdlref$autopar)[i] %in% names(autopar))
      {
       autopar[[names(mdlref$autopar)[i]]] <- mdlref$autopar[[names(mdlref$autopar)[i]]]
      }
     }
     #initializing hpar with mdlref$hpar
     autopar <- sbamlvalidautopar(autopar = autopar, myruntype = 'manual')[["autopar"]]
     myprovres <- sbamlmodhparcp(mdlref, myhpartoignore = myautomatablehpar, autopar, hpar)
     autopar <- myprovres[["autopar"]]
     hpar <- myprovres[["hpar"]]
     rm(myprovres)
    }
   }
   #
   mylspsoparam <- list(
   modeinit = autopar$psomodeinit,
   modecost = autopar$psomodecost,
   nbvar2optim = autopar$psonbvar2optim,
   partpopsize = autopar$psopartpopsize,
   numiteration = autopar$numiterations,
   varvalmin = -1, varvalmax = 1,
   velocitymaxratio = autopar$psovelocitymaxratio,
   inertiadampratio = autopar$psoinertiadampratio,
   kappa = autopar$psokappa,
   phi1 = autopar$psophi1,
   phi2 = autopar$psophi2,
   seed = autopar$seed)
   myprovres <- sbamlmodpsoparam(mylspsoparam = mylspsoparam, mylspsopart = NULL)
   mylspsoparam <- myprovres[["mylspsoparam"]]
   mylspsopart <- myprovres[["mylspsopart"]]
   rm(myprovres)
   #
   mymc <- sbamlparllmc(autopar$nbcores)
   sbamlpart2parll <- sbamlpart2parlln1(autopar = autopar, hpar = hpar, mylspsoparam = mylspsoparam)
   #
   for (iternbr in 1:autopar$numiterations)
   {
    if (mymc != 1)
    {
     mylsres <- parallel::mclapply(X = mylspsopart, FUN = sbamlpart2parll,
                                   mc.preschedule = FALSE, mc.set.seed = FALSE,
                                   mc.silent = FALSE, mc.cores = getOption("mc.cores", mymc),
                                   mc.cleanup = TRUE, mc.allow.recursive = FALSE)
    } else {
     mylsres <- lapply(X = mylspsopart, FUN = sbamlpart2parll)
    }
    mylspsopart <- mylsres
    for (i in 1:autopar$psopartpopsize)
    {
     myflagok <- 1
     if (class(mylspsopart[[i]]) != "list")
     {
      myflagok <- 0
     } else
     {
      if (any(!c('Cost', 'Best') %in% names(mylspsopart[[i]])))
      {
       myflagok <- 0
      } else
      {
       if (any(!c('Cost', 'model', 'Position') %in% names(mylspsopart[[i]]$Best)))
       {
        myflagok <- 0
       } else
       {
        if (any(!c('error', 'hpar') %in% names(mylspsopart[[i]]$Best$model)))
        {
         myflagok <- 0
        } else
        {
         if (mylspsopart[[i]]$Best$model$error$tr == Inf)
         {
          myflagok <- 0
         }
        }
       }
      }
     }
     if (iternbr == 1 & i == 1 & autopar$verbose == TRUE)
     {
      if (autopar$auto_runtype == "2steps")
      {
       cat(paste('STEP: ', autoloopnbr,' (',ifelse(autoloopnbr == 1, 'overfitting', 'regularization'),')\n', sep = ''))
      }
      if (myflagok == 1) {cat(paste('(cost: ', mylspsopart[[i]]$Best$model$hpar$costtype,')\n', sep = ''))}
     }
     mylastlog <- paste("iteration", iternbr,
                        "particle", i, sep = " ")
     if (myflagok == 1)
     {
      if (is.na(mylspsopart[[i]]$Best$model$error$cv))
      {
       mylastlog <- paste(mylastlog, "weighted err:", round(mylspsopart[[i]]$Best$Cost, digits = 5), sep = " ")
      } else
      {
       mylastlog <- paste(mylastlog, "weighted err:", round(mylspsopart[[i]]$Best$Cost, digits = 5),
                          "(train:", round(mylspsopart[[i]]$Best$model$error$tr, digits = 5),
                          "cvalid:", round(mylspsopart[[i]]$Best$model$error$cv, digits = 5),
                          ")", sep = " ")
      }
      if (mylspsopart[[i]]$Best$Cost < mylspsoparam$globalbest$Cost)
      {
       mylspsoparam$globalbest <- mylspsopart[[i]]$Best
       mylastlog <- paste(mylastlog, "BEST MODEL KEPT", sep = " ")
      }
     } else
     {
      mylastlog <- paste(mylastlog, "no exploitable results or timeout", sep = " ")
     }
     if (autopar$verbose == TRUE) {cat(paste0(mylastlog, "\n"))}
    }
    mylspsoparam$globbestcostlog[iternbr, 1] <- mylspsoparam$globalbest$idpart
    mylspsoparam$globbestcostlog[iternbr, 2] <- mylspsoparam$globalbest$Cost
   }
   mymodel <- mylspsoparam$globalbest$model
   mylspsoparam$globalbest$model <- NULL
   mymodel$mylspsoparamlvauto <- mylspsoparam
   mymodel$autopar <- autopar
   myprovres <- sbamltestmodel(mymodel = mymodel, mymodeltype = 'auto')
   mymodel <- myprovres[["mymodel"]]
   rm(myprovres)
   if (autopar$auto_runtype == "2steps" & autoloopnbr == 1)
   {
    hpar <- hparref; autopar <- autoparref
    mdlref <- mymodel
   } else {
    return(mymodel)
   }
  }
}

automl_train_manual <- function(Xref, Yref, hpar = list(), mdlref = NULL)
{
  hpar <- sbamlvalidhpar(hpar = hpar, myruntype = 'rawmanual')[["hpar"]]
  if (!is.null(mdlref))
  {
   myprovres <- sbamltestmodel(mymodel = mdlref, mymodeltype = 'manual', guess1stlayer = TRUE)
   mdlref <- myprovres[["mymodel"]]
   valmodimpflagok <- myprovres[["tstflagok"]]
   rm(myprovres)
   if (valmodimpflagok == 1)
   {
    myprovres <- sbamlmodhparcp(mdlref, myhpartoignore = NULL, autopar = NULL, hpar)
    hpar <- myprovres[["hpar"]]
    valmodimpflagok <- myprovres[["valmodimpflagok"]]
    rm(myprovres)
    if ('autopar' %in% names(mdlref))
    {
     hpar$verbose <- TRUE
     hpar$useautopar <- FALSE
    }
   }
  } else {valmodimpflagok <- 0}
  #
  hpar <- sbamlvalidhpar(hpar = hpar, myruntype = 'manual')[["hpar"]]
  if (hpar$useautopar == FALSE)
  {
    Xref <- t(as.matrix(Xref))
    Yref <- t(as.matrix(Yref))
    sbamlmatstest(Xref, Yref)
    Xref <- sbamlmattest(Xref, "Xref")
    Yref <- sbamlmattest(Yref, "Yref")
  }
  mydl <- list()
  mydl$error$cv <- mydl$error$tr <- Inf
  if (valmodimpflagok == 1)
  {
   mydl <- sbamlmodcp(mdlref, mydl)
  }
  mydl[['hpar']] <- hpar
  mydl$hpar$layersshape <- c(dim(Xref)[1], mydl$hpar$layersshape)
  mydl$hpar$layersshape[length(mydl$hpar$layersshape)] <- dim(Yref)[1]
  mydl$hpar$nblayers <- length(mydl$hpar$layersshape) - 1
  if (mydl$hpar$layersacttype[mydl$hpar$nblayers] == '' | mydl$hpar$modexec == 'trainwgrad')
  {
    if (any(!Yref %in% c(0, 1)))
    {
      mydl$hpar$layersacttype[mydl$hpar$nblayers] <- 'linear'
    } else {
      mydl$hpar$layersacttype[mydl$hpar$nblayers] <- 'sigmoid'
    }
  }
  if (mydl$hpar$costcustformul != '' & mydl$hpar$modexec == 'trainwpso')
  {
    mydl$hpar$costtype <- 'custom'
  }
  if (mydl$hpar$costtype == '' | mydl$hpar$modexec == 'trainwgrad')
  {
    if (any(!Yref %in% c(0, 1)))
    {
      mydl$hpar$costtype <- 'mse'
    } else {
      mydl$hpar$costtype <- 'crossentropy'
    }
  }
  mynbtestsoktostop <- 1;#legacy hyperparam dropped to simplify
  mytestokcpt <- 0;#to avoid error with test while not defined
  mytrlastcost <- mycvlastcost <- mydl$hpar$epsil ^ -1;#why not to get a huge number
  if (mydl$hpar$testcvsize != 0)
  {
    #shuffling handled upstream
    mycvpctg <- mydl$hpar$testcvsize / 100
    mycvgainunder <- mydl$hpar$testgainunder / 100
    mythld <- floor(dim(Xref)[2] * mycvpctg)
    if (mydl$hpar$useautopar == FALSE)
    {
     set.seed(mydl$hpar$seed)
    } else {
     set.seed(mydl$hpar$seedautopar)
    }
    mysample <- sample(x = dim(Xref)[2], size = mythld, replace = FALSE)
    Xcv <- matrix(Xref[,mysample], nrow = dim(Xref)[1], byrow = F)
    Ycv <- matrix(Yref[,mysample], nrow = dim(Yref)[1], byrow = F)
    Xref <- matrix(Xref[,-mysample], nrow = dim(Xref)[1], byrow = F)
    Yref <- matrix(Yref[,-mysample], nrow = dim(Yref)[1], byrow = F)
  } else {
    mycvpctg <- mycvgainunder <- 0
  }
  if (mydl$hpar$minibatchsize != 0 & dim(Xref)[2] > mydl$hpar$minibatchsize)
  {
    nbbatch <- ceiling(dim(Xref)[2] / mydl$hpar$minibatchsize)
  } else {
    mydl$hpar$minibatchsize <- dim(Xref)[2]
    nbbatch <- 1
  }
  if (mydl$hpar$modexec == 'trainwgrad')
  {
    if (valmodimpflagok == 0)
    {
     mydl <- sbamlmoddlparaminit(mydl, mydl$hpar$layersshape, mydl$hpar$layersacttype, mydl$hpar$seed, mydl$hpar$batchnor_mom)
    }
  } else if (mydl$hpar$modexec == 'trainwpso')
  {
    mylspsoparam <- list(
      modeinit = mydl$hpar$psomodeinit,
      modecost = mydl$hpar$psomodecost,
      partpopsize = mydl$hpar$psopartpopsize,
      nbvar2optim = NULL, numiteration = mydl$hpar$numiterations,
      varvalmin = mydl$hpar$psovarvalmin, varvalmax = mydl$hpar$psovarvalmax,
      velocitymaxratio = mydl$hpar$psovelocitymaxratio,
      inertiadampratio = mydl$hpar$psoinertiadampratio,
      kappa = mydl$hpar$psokappa, phi1 = mydl$hpar$psophi1, phi2 = mydl$hpar$psophi2,
      seed = mydl$hpar$seed)
    mylspsoparam$angdl$mydl <- mydl
    myprovres <- sbamlmodpsoparam(mylspsoparam = mylspsoparam, mylspsopart = NULL)
    mylspsoparam <- myprovres[["mylspsoparam"]]
    mylspsopart <- myprovres[["mylspsopart"]]
    rm(myprovres)
    if (valmodimpflagok == 1 & "mylspsoparam" %in% names(mdlref))
    {
     mylspsoparam$globalbest <- mdlref$mylspsoparam$globalbest
     mylspsoparam$globalbest$idpart <- 1
     mylspsopart[[1]]$Position <- mylspsopart[[1]]$Best$Position <- mdlref$mylspsoparam$globalbest$Position
     mylspsopart[[1]]$Cost <- mylspsopart[[1]]$Best$Cost <- mdlref$mylspsoparam$globalbest$Cost
    }
  }

  myflagcontinue <- 1
  t <- 1
  if (mydl$hpar$verbose == TRUE)
  {
    cat(paste('(cost: ', mydl$hpar$costtype,')\n', sep = ''))
  }
  for (epochnum in 1:mydl$hpar$numiterations)
  {
    minibatchposstart <- 1
    for (batchnum in 1:nbbatch)
    {
      #ds shuffled upstream ;-)
      minibatchposend <- batchnum * mydl$hpar$minibatchsize
      if (minibatchposend > dim(Xref)[2]) {minibatchposend <- dim(Xref)[2]}
      if (minibatchposstart > dim(Xref)[2]) {minibatchposstart <- dim(Xref)[2]}
      if (nbbatch != 1 | (epochnum == 1 & nbbatch == 1))
      {
        X <- matrix(Xref[,minibatchposstart:minibatchposend], nrow = dim(Xref)[1], byrow = F)
        mydl[["mydl_A0"]] <- X
        Y <- matrix(Yref[,minibatchposstart:minibatchposend], nrow = dim(Yref)[1], byrow = F)
        if (mydl$hpar$modexec == 'trainwpso')
        {
          mylspsoparam$angdl$X <- X
          mylspsoparam$angdl$Y <- Y
        }
      }
      if (mydl$hpar$modexec == 'trainwpso')
      {
        for (i in 1:mylspsoparam$partpopsize)
        {
          mylspsoparam$activpart <- i
          myprovres <- sbamlmodpsooptkern(mylspsoparam = mylspsoparam, mylspsoactivpart = mylspsopart[[i]])
          mylspsoparam <- myprovres[["mylspsoparam"]]
          mylspsopart[[i]] <- myprovres[["mylspsoactivpart"]]
          rm(myprovres)
        }
        mylspsoparam$globbestcostlog[epochnum, 1] <- mylspsoparam$globalbest$idpart
        mylspsoparam$globbestcostlog[epochnum, 2] <- mylspsoparam$globalbest$Cost
        mydl <- mylspsoparam$angdl$mydl;#to set best angdlpso param
      }
      myprovres <- sbamlmoddlfw(mydl, X, mydl$hpar$modexec, mydl$hpar$nblayers, mydl$hpar$layersacttype,
                                mydl$hpar$layersdropoprob, mydl$hpar$seed, mydl$hpar$batchnor_mom,
                                mydl$hpar$epsil, epochnum, batchnum)
      Yhat <- myprovres[["A"]]; mydl <- myprovres[["mydl"]]; rm(myprovres)
      if (mydl$hpar$printcostevery != 0)
      {
        if ((epochnum %% mydl$hpar$printcostevery == 0 | epochnum == mydl$hpar$numiterations) & batchnum %% nbbatch == 0)
        {
          mytrcost <- sbamlmoddlcost(mydl, Y, Yhat, mydl$hpar$costtype,
                                     mydl$hpar$lambda, mydl$hpar$nblayers,
                                     mydl$hpar$epsil)
          mydl$error$tr <- mytrcost
          if (mydl$hpar$verbose == TRUE) {cat(paste('cost epoch',epochnum , ': ', mytrcost, sep = ''))}
          if (mydl$hpar$testcvsize != 0 & mydl$hpar$overfitautopar == FALSE)
          {
            mylastlog <- " "; mytestokcpt <- 0
            myprovres <- sbamlmoddlfw(mydl, Xcv, 'predict', mydl$hpar$nblayers, mydl$hpar$layersacttype,
                                      mydl$hpar$layersdropoprob, mydl$hpar$seed, mydl$hpar$batchnor_mom,
                                      mydl$hpar$epsil)
            Ycvhat <- myprovres[["A"]]; mydl <- myprovres[["mydl"]]; rm(myprovres)
            mycvcost <- sbamlmoddlcost(mydl, Ycv, Ycvhat, mydl$hpar$costtype, mydl$hpar$lambda,
                                       mydl$hpar$nblayers, mydl$hpar$epsil)
            mydl$error$cv <- mycvcost
            if (mydl$hpar$verbose == TRUE) {cat(paste(' (cv cost: ', mycvcost, ')', sep = ''))}
            if (abs(mytrlastcost - mytrcost) < mydl$hpar$testgainunder)
            {
              mylastlog <- paste(mylastlog, "[test trgainunder OK]", sep = " ")
              mytestokcpt <- mytestokcpt + 1
            }
            mytrlastcost <- mytrcost
            if (abs(mycvlastcost - mycvcost) < mycvgainunder)
            {
              mylastlog <- paste(mylastlog, "[test cvgainunder OK]", sep = " ")
              mytestokcpt <- mytestokcpt + 1
            }
            if (mylastlog != ' ')
            {
              if (mydl$hpar$verbose == TRUE) {cat(paste(mylastlog), '\n', sep = '')}
            }
            mycvlastcost <- mycvcost
          } else {
            mydl$error$cv <- NA
          }
          if (mydl$hpar$verbose == TRUE) {cat(paste(' (LR: ', mydl$hpar$learningrate, ')',
                    '\n', sep = ' '))}
          if (mytestokcpt >= mynbtestsoktostop)
          {
            if (mydl$hpar$verbose == TRUE) {cat(mylastlog)}
            myflagcontinue <- 0
          }
        }
      }
      if (mydl$hpar$modexec == 'trainwgrad' & (epochnum != mydl$hpar$numiterations & batchnum != nbbatch))
      {
        mydl <- sbamlmoddlbk(mydl, Y, Yhat, mydl$hpar$costtype, mydl$hpar$nblayers,
                             mydl$hpar$layersacttype, mydl$hpar$lambda, mydl$hpar$layersdropoprob,
                             mydl$hpar$batchnor_mom, mydl$hpar$epsil, epochnum, batchnum)
        if (mydl$hpar$chkgradevery != 0)
        {
          if (epochnum %% mydl$hpar$chkgradevery == 0 & batchnum %% nbbatch == 0)
          {
            sbamlmoddlgradchk(mydl, X, Y, mydl$hpar$nblayers, mydl$hpar$chkgradepsilon, mydl$hpar$epsil)
          }
        }
        mydl <- sbamlmoddlparamupdt(mydl, mydl$hpar$nblayers, mydl$hpar$learningrate, mydl$hpar$beta1,
                                    mydl$hpar$beta2, t, mydl$hpar$epsil, mydl$hpar$batchnor_mom)
      }
      t <- t + 1
      minibatchposstart <- minibatchposend + 1
    }
    if (myflagcontinue == 0)
    {
      break
    }
    if (mydl$hpar$lrdecayrate != 0)
    {
      mydl$hpar$learningrate <- mydl$hpar$learningrate / (1 + mydl$hpar$lrdecayrate * epochnum)
    }
  }
  if (mydl$hpar$modexec == 'trainwpso')
  {
    mylspsoparam$angdl <- NULL
    mydl$mylspsoparam <- mylspsoparam
  }
  if (mydl$hpar$verbose == TRUE)
  {
    cat(paste('   dim X', ': [', paste(dim(Xref), collapse = ','), ']\n', sep = ''))
    for (l in 1:mydl$hpar$nblayers)
    {
      W <- mydl[[paste("mydl_W", l, sep = '')]]
      bB <- mydl[[paste("mydl_bB", l, sep = '')]]
      cat(paste('   dim W', l, ': [', paste(dim(W), collapse = ','),
                '] (min|max: ',min(W),', ',max(W),
                ')\n   dim bB', l, ': [',  paste(dim(bB), collapse = ','),
                '] (min|max: ',min(bB),', ',max(bB),
                ')\n', sep = ''))
      if (mydl$hpar$batchnor_mom != 0)
      {
        j <- mydl[[paste("mydl_j", l, sep = '')]]
        cat(paste('   dim j', l, ': [', paste(dim(j), collapse = ','),
                  '] (min|max: ',min(j),', ',max(j),
                  ')\n', sep = ''))
      }
    }
    cat(paste('   dim Y', ': [', paste(dim(Yref), collapse = ','), ']\n', sep = ''))
  }
  myprovres <- sbamltestmodel(mymodel = mydl, mymodeltype = 'manual', guess1stlayer = FALSE)
  mydl <- myprovres[["mymodel"]]
  valmodflagok <- myprovres[["tstflagok"]]
  rm(myprovres)
  if (valmodflagok == 0 & mydl$hpar$useautopar == TRUE)
  {
   mylastlog <- ("error: model training error")
   stop(mylastlog);#to return error to automl_train
  }
  for (l in 0:mydl$hpar$nblayers)
  {
   # lightening model
   try(mydl[[paste0("mydl_Z",l)]] <- NULL)
   try(mydl[[paste0("mydl_A",l)]] <- NULL)
   try(mydl[[paste0("mydl_Zn",l)]] <- NULL)
   try(mydl[[paste0("mydl_Zt",l)]] <- NULL)
  }
  return(mydl)
}

automl_predict <- function(model, X, layoutputnum = 0)
{
  Yhat <- sbamlmoddlfw(model, t(X),
                       'predict',
                       ifelse(layoutputnum == 0, model$hpar$nblayer, layoutputnum),
                       model$hpar$layersacttype,
                       c(rep(0, length(model$hpar$layersacttype))),
                       model$hpar$seed,
                       batchnor_mom = model$hpar$batchnor_mom, epsil = model$hpar$epsil)[["A"]]
  if (dim(Yhat)[1] == 1)
  {
    Yhat <- as.vector(Yhat)
  } else {
    Yhat <- as.data.frame(t(Yhat), row.names = NULL, stringsAsFactors = FALSE)
  }
  return(Yhat)
}

sbamlmoddlbk <- function(mydl, Y, Yhat, costtype, nblayers, layersacttype, lambda, layersdropoprob, batchnor_mom, epsil,
                         epochnum = 0, batchnum = 0)
{
  m <- dim(Y)[2]
  for (l in nblayers:1)
  {
    if (l == nblayers)
    {
      dA <- sbamlmoddlcostbk(Y,
                             Yhat,
                             costtype)
      mydl[[paste("mydl_dA", l, sep = '')]] <- dA
    } else
    {
     dA <- mydl[[paste("mydl_dA", l, sep = '')]]
     if (layersdropoprob[l] != 0 & epochnum != 0 & batchnum != 0)
     {
      D <- mydl[[paste("mydl_D", l, sep = '')]]
      dA <- sbamlmoddldropo(dA, D, layersdropoprob[l])
     }
    }
    if (batchnor_mom == 0)
    {
      dZ <- sbamlmoddlsubbkact(mydl, dA, l, layersacttype[l], batchnor_mom)
    } else
    {
      dZt <- sbamlmoddlsubbkact(mydl, dA, l, layersacttype[l], batchnor_mom)
      # BN
      Z <- mydl[[paste("mydl_Z", l, sep = '')]]
      M <- mydl[[paste("mydl_M", l, sep = '')]]
      V <- mydl[[paste("mydl_V", l, sep = '')]]
      Zn <- mydl[[paste("mydl_Zn", l, sep = '')]]

      dbB <- matrix(apply(dZt, 1, sum), ncol = 1, byrow = FALSE) * (1 / m)
      mydl[[paste("mydl_dbB", l, sep = '')]] <- dbB
      dj <- matrix(apply(sbamlmatopebroadcst(Zn,
                                             dZt, 'mult')
                         , 1, sum), ncol = 1, byrow = FALSE) * (1 / m)
      mydl[[paste("mydl_dj", l, sep = '')]] <- dj
      #dZ
      j <- mydl[[paste("mydl_j", l, sep = '')]]
      Z_M <- sbamlmatopebroadcst(Z, M * (-1), 'add')
      Z_Ms <- Z_M ^ 2
      Ves <- (V + epsil) ^ (1 / 2)
      Vesinv <- Ves ^ (-1)
      dZn <- sbamlmatopebroadcst(dZt, j, 'mult')
      dZ_M1 <- sbamlmatopebroadcst(dZn, Vesinv, 'mult')
      dVesinv <- matrix(apply(sbamlmatopebroadcst(dZn, Z_M, 'mult'), 1, sum),
                        ncol = 1, byrow = FALSE);# * (1 / m)
      dVes <- sbamlmatopebroadcst(dVesinv, (Ves ^ -2) * (-1), 'mult')
      dV <- sbamlmatopebroadcst(dVes, (Ves ^ -1) * 0.5, 'mult')
      dZ_Ms <- dV * (1 / m)
      dZ_M2 <- sbamlmatopebroadcst(Z_M * 2, dZ_Ms, 'mult')
      dZ_M <- sbamlmatopebroadcst(dZ_M1, dZ_M2, 'add')
      dM <- matrix(apply(dZ_M * (-1), 1, sum), ncol = 1, byrow = FALSE)
      dZ <- sbamlmatopebroadcst(dZ_M, dM * (1 / m), 'add')
    }
    Aprev <- mydl[[paste("mydl_A",(l - 1), sep = '')]]
    W <- mydl[[paste("mydl_W", l, sep = '')]]
    mygrads <- sbamlmoddlsubbklin(mydl, dZ, Aprev, W, l, lambda,
                                  batchnor_mom)
    mydl[[paste("mydl_dW", l, sep = '')]] <- mygrads$dW
    if (batchnor_mom == 0)
    {
      mydl[[paste("mydl_dbB", l, sep = '')]] <- mygrads$dbB
    }
    if (l != 1)
    {
      mydl[[paste("mydl_dA", (l - 1), sep = '')]] <- mygrads$dAprev
    }
  }
  return(mydl)
}

sbamlmoddlfw <- function(mydl, X, modexec, nblayers, layersacttype, layersdropoprob, seed, batchnor_mom, epsil,
                         epochnum = 0, batchnum = 0)
{
  Aprev <- X
  for (l in 1:nblayers)
  {
    W <- mydl[[paste("mydl_W", l, sep = '')]]
    bB <- mydl[[paste("mydl_bB", l, sep = '')]]
    if (batchnor_mom == 0)
    {
      Z <- sbamlmatopebroadcst(sbamlmatopebroadcst(W, Aprev, 'dot'), bB, 'add')
      A <- sbamlmoddlsubfwact(Z, layersacttype[l])
    } else {
      Z <- sbamlmatopebroadcst(W, Aprev, 'dot')
      if (modexec %in% c('trainwgrad', 'trainwpso'))
      {
        m = dim(Z)[2]
        M <- (1 / m) * matrix(apply(Z, 1, sum), nrow = dim(Z)[1], ncol = 1, byrow = FALSE)
        Z_M <- sbamlmatopebroadcst(Z, M * (-1), 'add')
        V <- (1 / m) * matrix(apply(Z_M ^ 2, 1, sum), nrow = dim(Z)[1], ncol = 1, byrow = FALSE)
        Vesinv <- (V + epsil) ^ (-1 / 2)
        Zn <- sbamlmatopebroadcst(Z_M, Vesinv, 'mult')
        mydl[[paste("mydl_M", l, sep = '')]] <- M
        mydl[[paste("mydl_V", l, sep = '')]] <- V
        if (modexec == 'trainwpso')
        {
          vM <- mydl[[paste("mydl_vM", l, sep = '')]]
          mydl[[paste("mydl_vM", l, sep = '')]] <- sbamlexponweigtdavg(M, vM, batchnor_mom, tbiascorrection = 1)
          vV <- mydl[[paste("mydl_vV", l, sep = '')]]
          mydl[[paste("mydl_vV", l, sep = '')]] <- sbamlexponweigtdavg(V, vV, batchnor_mom, tbiascorrection = 1)
        }
      } else {
        M <- mydl[[paste("mydl_vM", l, sep = '')]]
        Zn <- sbamlmatopebroadcst(Z, M * (-1), 'add')
        V <- mydl[[paste("mydl_vV", l, sep = '')]]
        Zn <- sbamlmatopebroadcst(Zn, (V + epsil) ^ (-1 / 2), 'mult')
      }
      j <- mydl[[paste("mydl_j", l, sep = '')]]
      Zt <- sbamlmatopebroadcst(sbamlmatopebroadcst(Zn, j, 'mult'), bB, 'add')
      A <- sbamlmoddlsubfwact(Zt, layersacttype[l])
    }
    if (modexec %in% c('trainwgrad', 'trainwpso'))
    {
      if (layersdropoprob[l] != 0 & epochnum != 0 & batchnum != 0)
      {
        set.seed(ifelse(epochnum + batchnum > .Machine$integer.max,
                        epochnum + batchnum - .Machine$integer.max,
                        epochnum + batchnum))
        D <- matrix(stats::runif(n = dim(A)[1] * dim(A)[2], min = 0, max = 1),
                    nrow = dim(A)[1], ncol = dim(A)[2], byrow = FALSE)
        D <- D > layersdropoprob[l]
        A <- sbamlmoddldropo(A, D, layersdropoprob[l])
      }
      mydl[[paste("mydl_Z", l, sep = '')]] <- Z
      mydl[[paste("mydl_A", l, sep = '')]] <- A
      if (layersdropoprob[l] != 0 & epochnum != 0 & batchnum != 0)
      {
        mydl[[paste("mydl_D", l, sep = '')]] <- D
      }
      if (batchnor_mom != 0)
      {
        mydl[[paste("mydl_Zn", l, sep = '')]] <- Zn
        mydl[[paste("mydl_Zt", l, sep = '')]] <- Zt
      }
    }
    if (l != nblayers) {Aprev <- A}
  }
  return(list(A = A, mydl = mydl))
}

sbamlmatopebroadcst <- function(mat1, mat2, optype = 'dot')
{
  # optype
  #  'add' add matrices
  #  'mult' mutiply matrices element wise
  #  'divide' divide matrices element wise
  #  'dot' mutiply matrices
  #
  if (is.null(dim(mat1)) | is.null(dim(mat2)))
  {
    mylastlog <- ("fctmatopebroadcst mat1 or/n 2 is not a matrix")
    stop(mylastlog)
  } else {
    n1 = dim(mat1)[1]
    m1 = dim(mat1)[2]
    n2 = dim(mat2)[1]
    m2 = dim(mat2)[2]
    mat1[is.na(mat1) | is.infinite(mat1)] <- 0
    mat2[is.na(mat2) | is.infinite(mat2)] <- 0
  }
  if ((optype %in% c('add', 'mult', 'divid') & n1 != n2) |
      (optype %in% c('dot') & m1 != n2))
  {
    mylastlog <- paste("fctmatopebroadcst mat1 n mat2 w wrong shape for", optype, sep = " ")
    stop(mylastlog)
  }
  if (optype %in% c('add', 'mult', 'divid'))
  {
    if (m2 == 1 & m2 != m1)
    {
      #broadcast
      mat2 <- matrix(rep(mat2[,1], m1), ncol = m1, byrow = FALSE)
    }
    if (optype == 'add')
    {
      matres <- mat1 + mat2
    } else if (optype == 'mult')
    {
      matres <- mat1 * mat2
    } else if (optype == 'divid')
    {
      mat2[mat2 == 0] <- 1e-12
      matres <- mat1 / mat2
    }
  } else if (optype == 'dot')
  {
    matres <- mat1 %*% mat2
  }
  return(matres)
}

sbamlmoddlsubbkact <- function(mydl, dA, l, activtype, batchnor_mom)
{
  if (activtype == 'sigmoid')
  {
    A <- mydl[[paste("mydl_A", l, sep = '')]]
    dZ <- sbamlmatopebroadcst(dA, sbamlmodsigmoidbk(A), 'mult')
  } else  if (activtype == 'relu')
  {
    if (batchnor_mom == 0)
    {
      Z <- mydl[[paste("mydl_Z", l, sep = '')]]
    } else {
      Z <- mydl[[paste("mydl_Zt", l, sep = '')]]
    }
    dZ <- sbamlmatopebroadcst(dA, sbamlmodrelubk(Z), 'mult')
  } else if (activtype == 'reluleaky')
  {
    if (batchnor_mom == 0)
    {
      Z <- mydl[[paste("mydl_Z", l, sep = '')]]
    } else {
      Z <- mydl[[paste("mydl_Zt", l, sep = '')]]
    }
    dZ <- sbamlmatopebroadcst(dA, sbamlmodreluleakybk(Z), 'mult')
  } else if (activtype == 'tanh')
  {
    A <- mydl[[paste("mydl_A", l, sep = '')]]
    dZ <- sbamlmatopebroadcst(dA, sbamlmodtanhbk(A), 'mult')
  } else if (activtype == 'linear')
  {
    dZ <- dA
  }
  return(dZ)
}

sbamlmoddlsubbklin <- function(mydl, dZ, Aprev, W, l, lambda, batchnor_mom)
{
  m <- dim(Aprev)[2]
  dW <- sbamlmatopebroadcst(dZ, t(Aprev), 'dot') * (1 / m)
  if (lambda != 0)
  {
    dW <- sbamlmatopebroadcst(dW , W * (lambda / m), 'add')
  }
  if (batchnor_mom == 0)
  {
    dbB <- matrix(apply(dZ, 1, sum), ncol = 1, byrow = FALSE) * (1 / m)
    res <- list(
      'dW' = dW,
      'dbB' = dbB
    )
  } else {
    res <- list(
      'dW' = dW
    )
  }

  if (l != 1)
  {
    dAprev <- sbamlmatopebroadcst(t(W), dZ, 'dot')
    res[['dAprev']] = dAprev
  }
  return(res)
}

sbamlmoddlsubfwact <- function(Z, activtype)
{
  if (activtype == 'sigmoid')
  {
    A <- sbamlmodsigmoid(Z)
  } else  if (activtype == 'relu')
  {
    A <- sbamlmodrelu(Z)
  } else if (activtype == 'reluleaky')
  {
    A <- sbamlmodreluleaky(Z)
  } else if (activtype == 'tanh')
  {
    A <- tanh(Z)
  } else if (activtype == 'linear')
  {
    A <- Z
  } else if (activtype == 'softmax')
  {
    A <- sbamlmodsoftmax(Z)
  }
  return(A)
}

sbamlmoddlparaminit <- function(mydl, layersshape, layersacttype, seed, batchnor_mom)
{
  for (l in 1:(length(layersshape) - 1))
  {
    n = layersshape[l + 1]
    m = layersshape[l]
    set.seed(seed + l)
    mymat <- matrix(stats::rnorm(m * n, mean = 0, sd = 1), nrow = n, ncol = m, byrow = FALSE)
    if (layersacttype[l] %in% c('relu', 'reluleaky'))
    {
      mymat <- mymat * sqrt(2 / m)
    } else if (layersacttype[l] == 'tanh')
    {
      mymat <- mymat * sqrt(1 / m)
    } else
    {
      mymat <- mymat * 0.01
    }

    mydl[[paste("mydl_W", l, sep = '')]] <- mymat
    mymat <- matrix(rep(0, m * n), nrow = n, ncol = m)
    mydl[[paste("mydl_vdW", l, sep = '')]] <- mymat
    mydl[[paste("mydl_sdW", l, sep = '')]] <- mymat
    mymat <- matrix(rep(0, n),nrow = n, ncol = 1)
    mydl[[paste("mydl_bB", l, sep = '')]] <- mymat
    mydl[[paste("mydl_vdbB", l, sep = '')]] <- mymat
    mydl[[paste("mydl_sdbB", l, sep = '')]] <- mymat
    #if (batchnor_mom != 0);#legacy: same architecture 4 all since mdlref
    mymat <- matrix(rep(1, n),nrow = n, ncol = 1)
    mydl[[paste("mydl_j", l, sep = '')]] <- mymat
    mydl[[paste("mydl_vdj", l, sep = '')]] <- mymat
    mydl[[paste("mydl_sdj", l, sep = '')]] <- mymat
    mymat <- matrix(rep(0, n),nrow = n, ncol = 1)
    mydl[[paste("mydl_M", l, sep = '')]] <- mymat
    mydl[[paste("mydl_vM", l, sep = '')]] <- mymat
    mydl[[paste("mydl_V", l, sep = '')]] <- mymat
    mydl[[paste("mydl_vV", l, sep = '')]] <- mymat
  }
  return(mydl)
}

sbamlmoddlparamupdt <- function(mydl, nblayers, dllearngrate, beta1, beta2, t, epsil, batchnor_mom)
{
  myflagok <- 1
  for (l in 1:nblayers)
  {
    W <- mydl[[paste("mydl_W", l, sep = '')]]
    dW <- mydl[[paste("mydl_dW", l, sep = '')]]
    bB <- mydl[[paste("mydl_bB", l, sep = '')]]
    dbB <- mydl[[paste("mydl_dbB", l, sep = '')]]
    if (beta1 != 0)
    {
      vdW <- mydl[[paste("mydl_vdW", l, sep = '')]]
      vdW <- sbamlexponweigtdavg(dW, vdW, beta1, t)
      mydl[[paste("mydl_vdW", l, sep = '')]] <- vdW
      vdbB <- mydl[[paste("mydl_vdbB", l, sep = '')]]
      vdbB <- sbamlexponweigtdavg(dbB, vdbB, beta1, t)
      mydl[[paste("mydl_vdbB", l, sep = '')]] <- vdbB
    }
    if (beta2 != 0)
    {
      sdW <- mydl[[paste("mydl_sdW", l, sep = '')]]
      sdW <- sbamlexponweigtdavg(dW ^ 2, sdW, beta2, t)
      mydl[[paste("mydl_sdW", l, sep = '')]] <- sdW
      sdbB <- mydl[[paste("mydl_sdbB", l , sep = '')]]
      sdbB <- sbamlexponweigtdavg(dbB ^ 2, sdbB, beta2, t)
      mydl[[paste("mydl_sdbB", l, sep = '')]] <- sdbB
    }
    if (beta1 != 0 & beta2 == 0 )
    {
      #gradient with Momentum
      dW <- vdW
      dbB <- vdbB
    } else if (beta1 == 0 & beta2 != 0 )
    {
      #gradient with RMSprop
      dW <- sbamlmatopebroadcst(dW, sqrt(sdW) + epsil, 'divid')
      dbB <- sbamlmatopebroadcst(dbB, sqrt(sdbB) + epsil, 'divid')
    } else if (beta1 != 0 & beta2 != 0 )
    {
      #gradient with adam optimization
      dW <- sbamlmatopebroadcst(vdW, sqrt(sdW) + epsil, 'divid')
      dbB <- sbamlmatopebroadcst(vdbB, sqrt(sdbB) + epsil, 'divid')
    }
    W <- sbamlmatopebroadcst(W, (-dllearngrate * dW), 'add')
    mydl[[paste("mydl_W", l, sep = '')]] <- W
    bB <- sbamlmatopebroadcst(bB, (-dllearngrate * dbB), 'add')
    mydl[[paste("mydl_bB", l, sep = '')]] <- bB
    if (any(is.na(W)) | any(is.na(bB)))
    {
      mylastlog <- "error: NaN produced in W or bB \n"
      myflagok <- 0
    }
    # BN
    if (batchnor_mom != 0)
    {
      j <- mydl[[paste("mydl_j", l, sep = '')]]
      dj <- mydl[[paste("mydl_dj", l, sep = '')]]
      if (beta1 != 0)
      {
        vdj <- mydl[[paste("mydl_vdj", l, sep = '')]]
        vdj <- sbamlexponweigtdavg(dj, vdj, beta1, t)
        mydl[[paste("mydl_vdj", l, sep = '')]] <- vdj
      }
      if (beta2 != 0)
      {
        sdj <- mydl[[paste("mydl_sdj", l, sep = '')]]
        sdj <- sbamlexponweigtdavg(dj ^ 2, sdj, beta2, t)
        mydl[[paste("mydl_sdj", l, sep = '')]] <- sdj
      }
      if (beta1 != 0 & beta2 == 0 )
      {
        #gradient with Momentum
        dj <- vdj
      } else if (beta1 == 0 & beta2 != 0 )
      {
        #gradient with RMSprop
        dj <- sbamlmatopebroadcst(dj, sqrt(sdj) + epsil, 'divid')
      } else if (beta1 != 0 & beta2 != 0 )
      {
        #gradient with adam optimization
        dj <- sbamlmatopebroadcst(vdj, sqrt(sdj) + epsil, 'divid')
      }
      j <- sbamlmatopebroadcst(j, (-dllearngrate * dj), 'add')
      if (any(is.na(j)))
      {
        mylastlog <- "error: NaN produced in j \n"
        myflagok <- 0
      }
      mydl[[paste("mydl_j", l, sep = '')]] <- j

      M <- mydl[[paste("mydl_M", l, sep = '')]]
      vM <- mydl[[paste("mydl_vM", l, sep = '')]]
      mydl[[paste("mydl_vM", l, sep = '')]] <- sbamlexponweigtdavg(M, vM, batchnor_mom, t)
      V <- mydl[[paste("mydl_V", l, sep = '')]]
      vV <- mydl[[paste("mydl_vV", l, sep = '')]]
      mydl[[paste("mydl_vV", l, sep = '')]] <- sbamlexponweigtdavg(V, vV, batchnor_mom, t)
    }
    if (myflagok == 0)
    {
      cat(mylastlog)
      stop(mylastlog, call. = FALSE)
    }
  }
  return(mydl)
}

sbamlexponweigtdavg <- function(valtosmooth, prevval, Beta, tbiascorrection = 1)
{
  if (tbiascorrection != 1)
  {
    if (is.null(dim(valtosmooth)))
    {
      newval <- Beta * prevval + valtosmooth * (1 - Beta)
    } else {
      newval <- sbamlmatopebroadcst(Beta * prevval,
                                    valtosmooth * (1 - Beta),
                                    'add')
    }
  } else {
    newval <- valtosmooth
    #newval <- newval * 1 / (1 - Beta ^ tbiascorrection)
  }
  return(newval)
}

sbamlmodrelu <- function(x)
{
  res = ifelse(x < 0, 0, x)
  return(res)
}

sbamlmodrelubk <- function(x)
{
  res = ifelse(x < 0, 0, 1)
  return(res)
}

sbamlmodreluleaky <- function(x)
{
  res = ifelse(x < 0, 0.01 * x, x)
  return(res)
}

sbamlmodreluleakybk <- function(x)
{
  res = ifelse(x < 0, 0.01, 1)
  return(res)
}

sbamlmodsigmoid <- function(z)
{
  g = (1 + exp(-z)) ^ (-1)
  return(g)
}

sbamlmodsigmoidbk <- function(x)
{
  res = x * (1 - x)
  return(res)
}

sbamlmodsoftmax <- function(z)
{
  g = exp(z)
  gs = matrix(apply(g, 2, sum), ncol = ncol(g), nrow = nrow(g), byrow = TRUE)
  g = sbamlmatopebroadcst(g, gs, 'divid')
  return(g)
}

sbamlmodtanhbk <- function(x)
{
  res = 1 - x^2
  return(res)
}

sbamlllnvectnorm <- function(x, n = 2)
{
  res = (sum(x^n))^(1 / n)
  return(res)
}

sbamlgetintdec <- function(x, part='int', nbedec = 0)
{
  #transform negatives val to positives ones
  x = abs(x)
  nbedec2 = nbedec + 1
  pint = as.integer(floor(x + (1/(10 ^ nbedec2))))
  pdec = as.integer(floor(10 ^ nbedec*((x + (1/(10 ^ nbedec2)) - pint))))
  res <- ifelse(part == 'int', pint,pdec)
  return(res)
}

sbamlmattest <- function(mymat, mymatname)
{
  if (is.numeric(mymat))
  {
    if (any(is.na(mymat)))
    {
      mymat[is.na(mymat)] <- 0
      mylastlog <- paste("wrng: NAs replaced by 0s in", mymatname, sep = " ")
      warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
    }
    if (any(is.infinite(mymat)))
    {
      mymat[is.infinite(mymat)] <- 0
      mylastlog <- paste("wrng: infinite values replaced by 0s in", mymatname, sep = " ")
      warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
    }
  } else {
    mylastlog <- paste("error:", mymatname, "type is not numeric", sep = " ")
    stop(mylastlog)
  }
 return(mymat)
}

sbamlmatstest <- function(Xref, Yref)
{
  #on transposed matrices
  if (dim(Xref)[2] != dim(Yref)[2])
  {
    mylastlog <- ("error: Xref and Yref have different number of rows")
    stop(mylastlog)
  }
}

sbamlpart2mydl <- function(mylspsoactivpart)
{
  mydl <- mylspsoactivpart$angdl$mydl
  myvectrefpos1 <- 1
  for (l in 1:mydl$hpar$nblayers)
  {
    myvectrefposname <- paste("mydl_W", l, sep = '')
    myprov <- dim(mydl[[myvectrefposname]])
    myvectrefpos2 <- c(myvectrefpos1 + myprov[1] * myprov[2] - 1)
    mymat <- matrix(mylspsoactivpart$Position[myvectrefpos1:myvectrefpos2], nrow = myprov[1], ncol = myprov[2], byrow = F)
    mydl[[myvectrefposname]] <- mymat
    myvectrefpos1 <- myvectrefpos2 + 1
    myvectrefposname <- paste("mydl_bB", l, sep = '')
    myprov <- dim(mydl[[myvectrefposname]])
    myvectrefpos2 <- c(myvectrefpos1 + myprov[1] * myprov[2] - 1)
    mymat <- matrix(mylspsoactivpart$Position[myvectrefpos1:myvectrefpos2], nrow = myprov[1], ncol = myprov[2], byrow = F)
    mydl[[myvectrefposname]] <- mymat
    myvectrefpos1 <- myvectrefpos2 + 1
    if (mydl$hpar$batchnor_mom != 0)
    {
      myvectrefposname <- paste("mydl_j", l, sep = '')
      myprov <- dim(mydl[[myvectrefposname]])
      myvectrefpos2 <- c(myvectrefpos1 + myprov[1] * myprov[2] - 1)
      mymat <- matrix(mylspsoactivpart$Position[myvectrefpos1:myvectrefpos2], nrow = myprov[1], ncol = myprov[2], byrow = F)
      mydl[[myvectrefposname]] <- mymat
      myvectrefpos1 <- myvectrefpos2 + 1
    }
  }
  return(mydl)
}

sbamlvalidautopar <- function(autopar, myruntype = 'manual')
{
 #myruntype:
 # = 'manual' produce error on error
 # = 'auto' produce warning on error
 # = 'rawmanual' only raw test produce error on error
 # = 'rawauto' only raw test produce warning on error
 myflagok <- 1
 if (class(autopar) != 'list')
 {
  if (myruntype %in% c('manual', 'rawmanual'))
  {
   mylastlog <- ("error: autopar is not a list")
   stop(mylastlog)
  } else if (myruntype %in% c('auto', 'rawauto')) {
   myflagok <- 0
   mylastlog <- ("wrng: autopar is not a list")
   warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
  }
 }
 if (myflagok == 1 & !myruntype %in% c('rawmanual', 'rawauto'))
 {
  #
  myprov <- names(autopar)
  if (!"seed" %in% myprov) {autopar[["seed"]] <- 4}
  if (!"nbcores" %in% myprov) {autopar[["nbcores"]] <- 1}
  if (!"verbose" %in% myprov) {autopar[["verbose"]] <- TRUE}
  if (!"subtimelimit" %in% myprov) {autopar[["subtimelimit"]] <- 3600}
  if (!"psomodeinit" %in% myprov) {autopar[["psomodeinit"]] <- 'autopar'}
  if (!"psomodecost" %in% myprov) {autopar[["psomodecost"]] <- 'autopar'}
  if (!"psopartpopsize" %in% myprov) {autopar[["psopartpopsize"]] <- 8}
  if (!"numiterations" %in% myprov) {autopar[["numiterations"]] <- 3}
  if (!"psovelocitymaxratio" %in% myprov) {autopar[["psovelocitymaxratio"]] <- 0.2}
  if (!"psoinertiadampratio" %in% myprov) {autopar[["psoinertiadampratio"]] <- 1}
  if (!"psokappa" %in% myprov) {autopar[["psokappa"]] <- 1}
  if (!"psophi1" %in% myprov) {autopar[["psophi1"]] <- 2.05}
  if (!"psophi2" %in% myprov) {autopar[["psophi2"]] <- 2.05}
  if (!"auto_modexec" %in% myprov) {autopar[["auto_modexec"]] <- FALSE}
  if (!"auto_runtype" %in% myprov) {autopar[["auto_runtype"]] <- 'normal'}
  if (!"auto_minibatchsize" %in% myprov) {autopar[["auto_minibatchsize"]] <- TRUE}
  if (!"auto_minibatchsize_min" %in% myprov) {autopar[["auto_minibatchsize_min"]] <- 0}
  if (!"auto_minibatchsize_max" %in% myprov) {autopar[["auto_minibatchsize_max"]] <- 9}
  if (!"auto_learningrate" %in% myprov) {autopar[["auto_learningrate"]] <- TRUE}
  if (!"auto_learningrate_min" %in% myprov) {autopar[["auto_learningrate_min"]] <- -5}
  if (!"auto_learningrate_max" %in% myprov) {autopar[["auto_learningrate_max"]] <- -2}
  if (!"auto_beta1" %in% myprov) {autopar[["auto_beta1"]] <- TRUE}
  if (!"auto_beta2" %in% myprov) {autopar[["auto_beta2"]] <- TRUE}
  if (!"auto_psopartpopsize" %in% myprov) {autopar[["auto_psopartpopsize"]] <- TRUE}
  if (!"auto_psopartpopsize_min" %in% myprov) {autopar[["auto_psopartpopsize_min"]] <- 2}
  if (!"auto_psopartpopsize_max" %in% myprov) {autopar[["auto_psopartpopsize_max"]] <- 50}
  if (!"auto_lambda" %in% myprov) {autopar[["auto_lambda"]] <- FALSE}
  if (!"auto_lambda_min" %in% myprov) {autopar[["auto_lambda_min"]] <- -2}
  if (!"auto_lambda_max" %in% myprov) {autopar[["auto_lambda_max"]] <- 4}
  if (!"auto_psovelocitymaxratio" %in% myprov) {autopar[["auto_psovelocitymaxratio"]] <- TRUE}
  if (!"auto_psovelocitymaxratio_min" %in% myprov) {autopar[["auto_psovelocitymaxratio_min"]] <- 0.01}
  if (!"auto_psovelocitymaxratio_max" %in% myprov) {autopar[["auto_psovelocitymaxratio_max"]] <- 0.5}
  if (!"auto_layers" %in% myprov) {autopar[["auto_layers"]] <- TRUE}
  if (!"auto_layers_min" %in% myprov) {autopar[["auto_layers_min"]] <- 1}
  if (!"auto_layers_max" %in% myprov) {autopar[["auto_layers_max"]] <- 2}
  if (!"auto_layersnodes_min" %in% myprov) {autopar[["auto_layersnodes_min"]] <- 3}
  if (!"auto_layersnodes_max" %in% myprov) {autopar[["auto_layersnodes_max"]] <- 33}
  if (!"auto_layersdropo" %in% myprov) {autopar[["auto_layersdropo"]] <- FALSE}
  if (!"auto_layersdropoprob_min" %in% myprov) {autopar[["auto_layersdropoprob_min"]] <- 0.05}
  if (!"auto_layersdropoprob_max" %in% myprov) {autopar[["auto_layersdropoprob_max"]] <- 0.75}
  rm(myprov)
  #
  myprov <- c(autopar$auto_minibatchsize,
              autopar$auto_learningrate,
              autopar$auto_beta1,
              autopar$auto_beta2,
              autopar$auto_psopartpopsize,
              autopar$auto_lambda,
              autopar$auto_psovelocitymaxratio)
  if (autopar$auto_layers == TRUE | autopar$auto_layersdropo == TRUE)
  {
   myprov <- c(myprov, rep(1, (1 + 3 * autopar$auto_layers_max)));#because layers nb, nodes nbs, act types, dropout
  }
  #100 reservations for positions initialization reproductibility
  if (length(myprov) > 100)
  {
   mylastlog <- ("wrng: positions initialization reproductibility")
   warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
  }
  autopar$psonbvar2optim <- ifelse(length(myprov) > 100, length(myprov), 100);

  if (sum(myprov) == 0 & autopar$numiterations > 1)
  {
   if (myruntype %in% c('manual'))
   {
    mylastlog <- ("error: no auto hpar set to TRUE in autopar")
    stop(mylastlog)
   } else if (myruntype %in% c('auto')) {
    myflagok <- 0
    mylastlog <- ("wrng: no auto hpar set to TRUE in autopar")
    warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
   }
  }
  rm(myprov)
 }
 return(list(autopar = autopar, tstflagok = myflagok))
}

sbamlvalidhpar <- function(hpar, myruntype = 'manual')
{
 #myruntype:
 # = 'manual' produce error on error
 # = 'auto' produce warning on error
 # = 'rawmanual' only raw test produce error on error
 # = 'rawauto' only raw test produce warning on error
 myflagok <- 1
 if (class(hpar) != 'list')
 {
  if (myruntype %in% c('manual', 'rawmanual'))
  {
   mylastlog <- ("error: hpar is not a list")
   stop(mylastlog)
  } else if (myruntype %in% c('auto', 'rawauto')) {
   myflagok <- 0
   mylastlog <- ("wrng: hpar is not a list")
   warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
  }
 }
 if (myflagok == 1 & !myruntype %in% c('rawmanual', 'rawauto'))
 {
  if (!"modexec" %in% names(hpar)) {hpar[["modexec"]] <- 'trainwgrad'}
  if (hpar$modexec == 'trainwgrad')
  {
   hpar$psopartpopsize <- 0
   hpar$psonbvar2optim <- 0
   hpar$psonumiteration <- 0
   hpar$psovarvalmin <- 0
   hpar$psovarvalmax <- 0
   hpar$psovelocitymaxratio <- 0
   hpar$psoinertiadampratio <- 0
   hpar$psokappa <- 0
   hpar$psophi1 <- 0
   hpar$psophi2 <- 0
   hpar$psoseed <- 0
  } else if (hpar$modexec == 'trainwpso')
  {
   if (!"minibatchsize" %in% names(hpar)) {hpar[["minibatchsize"]] <- 0}
   hpar$learningrate <- 0
   hpar$lrdecayrate <- 0
   hpar$chkgradevery <- 0
   hpar$chkgradepsilon <- 0
   hpar$beta1 <- 0
   hpar$beta2 <- 0
  }
  myprov <- names(hpar)
  if (!"useautopar" %in% myprov) {hpar[["useautopar"]] <- FALSE}
  if (!"overfitautopar" %in% myprov) {hpar[["overfitautopar"]] <- FALSE}
  if (!"seed" %in% myprov) {hpar[["seed"]] <- 4}
  if (!"minibatchsize" %in% myprov) {hpar[["minibatchsize"]] <- 2^5}
  if (!"layersshape" %in% myprov) {hpar[["layersshape"]] <- c(10, 0)}
  if (!"layersacttype" %in% myprov) {hpar[["layersacttype"]] <- c('relu', '')}
  if (!"layersdropoprob" %in% myprov) {hpar[["layersdropoprob"]] <- c(0, 0)}
  if (!"costtype" %in% myprov) {hpar[["costtype"]] <- ''}
  if (!"costcustformul" %in% myprov) {hpar[["costcustformul"]] <- ''}
  if (!"printcostevery" %in% myprov) {hpar[["printcostevery"]] <- 10}
  if (!"learningrate" %in% myprov) {hpar[["learningrate"]] <- 0.001}
  if (!"lrdecayrate" %in% myprov) {hpar[["lrdecayrate"]] <- 0}
  if (!"numiterations" %in% myprov) {hpar[["numiterations"]] <- 50}
  if (!"chkgradevery" %in% myprov) {hpar[["chkgradevery"]] <- 0}
  if (!"chkgradepsilon" %in% myprov) {hpar[["chkgradepsilon"]] <- 1e-7}
  if (!"lambda" %in% myprov) {hpar[["lambda"]] <- 0}
  if (!"beta1" %in% myprov) {hpar[["beta1"]] <- 0.9}
  if (!"beta2" %in% myprov) {hpar[["beta2"]] <- 0.999}
  if (!"epsil" %in% myprov) {hpar[["epsil"]] <- 1e-12}
  if (!"batchnor_mom" %in% myprov) {hpar[["batchnor_mom"]] <- 0.9}
  if (!"testcvsize" %in% myprov) {hpar[["testcvsize"]] <- 10}
  if (!"testgainunder" %in% myprov) {hpar[["testgainunder"]] <- 0.000001}
  if (!"psomodeinit" %in% myprov) {hpar[["psomodeinit"]] <- 'angdlpso'}
  if (!"psomodecost" %in% myprov) {hpar[["psomodecost"]] <- 'angdlpso'}
  if (!"psopartpopsize" %in% myprov) {hpar[["psopartpopsize"]] <- 50}
  if (!"psovarvalmin" %in% myprov) {hpar[["psovarvalmin"]] <- -10}
  if (!"psovarvalmax" %in% myprov) {hpar[["psovarvalmax"]] <- 10}
  if (!"psovelocitymaxratio" %in% myprov) {hpar[["psovelocitymaxratio"]] <- 0.2}
  if (!"psoinertiadampratio" %in% myprov) {hpar[["psoinertiadampratio"]] <- 1}
  if (!"psokappa" %in% myprov) {hpar[["psokappa"]] <- 1}
  if (!"psophi1" %in% myprov) {hpar[["psophi1"]] <- 2.05}
  if (!"psophi2" %in% myprov) {hpar[["psophi2"]] <- 2.05}
  if (!"verbose" %in% myprov) {hpar[["verbose"]] <- TRUE}
  if (length(hpar$layersshape) != length(hpar$layersacttype) |
      length(hpar$layersshape) != length(hpar$layersdropoprob))
  {
   mylastlog <- paste("layersshape, layersacttype and layersdropoprob must have same length (",
                      length(hpar$layersshape), length(hpar$layersacttype), length(hpar$layersdropoprob),
                      ")", sep = " ")
   if (myruntype == 'manual')
   {
    mylastlog <- paste("error: ", mylastlog, sep = " ")
    stop(mylastlog)
   } else {
    myflagok <- 0
    mylastlog <- paste("wrng: ", mylastlog, sep = " ")
    warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
   }
  }
  rm(myprov)
 }
 return(list(hpar = hpar, tstflagok = myflagok))
}

sbamlmodresult <- function(trerr, cverr, cvflag = TRUE)
{
  if (is.na(trerr)) {trerr <- Inf}
  if (is.na(cverr)) {cverr <- Inf}
  if (cvflag == TRUE)
  {
    err  <- abs(trerr - cverr) + ((1.5 * cverr + trerr) / 2.5)
  } else {
    err <- trerr
  }
  if (is.na(err)) {err <- Inf}
  return(err)
}

sbamlparllmc <- function(nbcores)
{
  myprov <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
  mymc <- ifelse(nbcores <= myprov,
                 nbcores,
                 myprov)
  if (gregexpr("dows", Sys.info()["sysname"]) > 0)
  {
    mymc <- 1
  }
  return(mymc)
}

sbamlprepscalem11201 <- function(x)
{
  z <- (x + 1)/2
  return(z)
}

sbamltestmodel <- function(mymodel, mymodeltype = 'manual', guess1stlayer = FALSE)
{
 #myruntype = 'manual' or 'auto' produce only warning on error
 myflagok <- 1
 if (class(mymodel) != 'list')
 {
  mymodel <- list()
  mylastlog <- ("wrng: no exploitable model")
  warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
  myflagok <- 0
 } else
 {
  if (any(!c('error', 'hpar') %in% names(mymodel)))
  {
   mylastlog <- ("wrng: no exploitable model")
   warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
   myflagok <- 0
  } else
  {
   if (myflagok == 1 & any(!c('useautopar',
                              'layersshape') %in% names(mymodel$hpar)))
   {
    mylastlog <- ("wrng: missing hpar parameters in model")
    warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
    myflagok <- 0
   }
   if (myflagok == 1 & guess1stlayer == TRUE)
   {
    mymodel$hpar[["layersshape"]] <- mymodel$hpar$layersshape[2:length(mymodel$hpar$layersshape)]
   }
   if (myflagok == 1 & (is.infinite(mymodel$error$tr) | is.infinite(mymodel$error$cv)))
   {
    mylastlog <- ("wrng: model create infinite error")
    warning(mylastlog, immediate. = FALSE, noBreaks. = TRUE)
    myflagok <- 0
   }
  }
  if (myflagok == 1 & mymodeltype == 'auto' & any(!c('autopar',
                                                     'mylspsoparamlvauto') %in% names(mymodel)))
  {
   myflagok <- 0
   mylastlog <- "wrng: model not trained with automl_train (ignored)"
  }
 }
 if (myflagok == 0)
 {
  mymodel <- NULL
 }
 return(list(mymodel = mymodel, tstflagok = myflagok))
}
