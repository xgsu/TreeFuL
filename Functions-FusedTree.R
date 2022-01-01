
# #########################################
# #########################################
# FUNCTIONS FOR BOOTSTRAP CALIBRATION
# #########################################
# #########################################


library(rpart)
# install.packages("rpart.plot")
library(rpart.plot)
library(tidyverse)
# install.packages("genlasso")
library(genlasso)
# install.packages("party")
library(party)
# install.packages("partykit")
library("partykit")
# install.packages("binaryLogic")
library(binaryLogic)
library("randomcoloR")  # COLORING
library("viridis")
library(mlbench)   # DATASETS
library("MASS")


# ================================
# FUNCTIONS THAT GENERATE DATA
# ================================
# MARS & TREE MODEL
rdat <- function(n, p=5, model=c("MARS1", "MARS2", "Tree1", "Tree2", "UpDown"), 
                      sigma=1, digits=3,
                      b=5, ngroup=5)  # ONE OPTION IN UpDown MODEL
{
	X <- NULL
	for (j in 1:p) {
		x <- round(runif(n), digits = digits)
		assign(paste("x", j, sep=""), x)
		X <- cbind(X, x)
	}
	if (model=="MARS1") {mu <- 0.1*exp(4*x1) + 4/(1+exp(-20*(x2-0.5))) + 3*x3 + 2*x4 + x5}
	else if (model=="MARS2") {mu <- 10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + x5}
	else if (model=="linear") {mu <- 2 + b*x1 + b*x2 + b*x3} 
	else if (model=="Tree1") {mu <- 2 + b* sign(0.25<= x1 & x1 <= 0.75)*sign(x2>= 0.25 & x2 <=0.75)}
	else if (model=="Tree2") {mu <- 2 + b*sign(x1 <= 0.5) + b*sign(x2 <=0.5)}
	else if (model=="UpDown") {
	  percentage <-  seq(0, 1, 1/ngroup)
	  # xcut <- cut(x1, breaks=quantile(x, probs=percentage), include.lowest=TRUE)
	  xcut <- cut(x1, breaks=percentage, include.lowest=TRUE)
	  X0 <- model.matrix(~ factor(xcut) - 1)
	  mu <- X0 %*% (rep(c(b, 0), ngroup)[1:ngroup])
	}
	else stop("The arugment model= needs to be either A or B.")
	y <- mu + rnorm(n, mean = 0, sd=sigma)
	dat <- data.frame(cbind(y, X, mu))
	names(dat) <- c("y", paste("x", 1:p, sep=""), "mu")
	return(dat)
}



# ===============================================
# FUNCTION cart0() GROWS A LARGE INITIAL TREE 
# ===============================================
# n0 = "minsplit" IN rpart

cart0 <- function(formula, data, method="anova", n0=NROW(data)/15, 
                  max.leaves = NULL)  # SOME PRUNING
{
	control0 <- rpart.control(minsplit=n0, minbucket=round(n0/3), 
	          maxdepth=8, cp=0, maxcompete=0,   # NUMBER OF COMPETITIVE SPLITS
            maxsurrogate=0, usesurrogate=2, surrogatestyle=0,  	# SURROGATE SPLITS FOR MISSING DATA
            xval=0)  # NO V-FOLD CV NEEDED   
	tre0 <- rpart(formula=formula, data=data, method=method, control=control0, 
	              model=TRUE);   # IMPORTANT TO SET model=TRUE TO HAVE TRAINING DATA AVAILABLE FOR FUTURE USE
	# PRUNE TRE0 A LITTLE
	if (!is.null(max.leaves)) {
	  CP.table <- tre0$cptable %>% as.data.frame
	  if (max(CP.table$nsplit) > max.leaves) {  
	    cp0 <- CP.table$CP[which(CP.table$nsplit <= max.leaves)] %>% last 
	    tre0 <- prune(tre0, cp=cp0)
	  }
	}
	return(tre0)
}


# ==================================
# SEND A TREE DOWN A DATASET 
# ==================================
senddown <- function(tree, data){
  data$node <- rpart:::pred.rpart(tree, x=rpart:::rpart.matrix(data)); 
	return(data)
}

# USING party - (BETTER AT HANDLING CATEGORICAL PREDICTORS?)
send.down <- function(tree, data){
  tre <- as.party(tree)
  data$node <- predict(tre, newdata = data, type = "node")
  return(data)
}




# ================
# FUSED LASSO FIT
# ================
Flasso <- function(y, node, 
                   adaptive=TRUE, m0=0.01, # M0=1,
                   modify.D=TRUE, alpha0=0.001, 
                   return.X=FALSE) {
  if (length(y)!= length(node)) stop("Error: y and node should have the same length.")
  n.node <- length(unique(node))
  fit0.anova <- lm(y~factor(node)-1, x =TRUE)
  X <- fit0.anova$x
  # X <- model.matrix(~factor(node)-1)
  
  if (is.ordered(node)) {
    Nd <- levels(node) 
    D <- getD1d(n.node)
    if (modify.D) {    # MODIFY GRAPH ADJACENCY MATRIX
      for (j in 1:(n.node-1)) {
        pvalue.t <- t.test(x=y[node==Nd[j]], y = y[node==Nd[j+1]], 
                           alternative="two.sided", var.equal = FALSE)$p.value 
        # print(pvalue.t)
        if (pvalue.t <= alpha0) {D[j, j] <- D[j, j+1] <- 0}
      }
    }
  } else {
    Nd <- unique(node)
    graf <- make_full_graph(n.node, directed = FALSE)
    D <- getDg(graf)
    k0 <- 0
    if (modify.D) { # MODIFY GRAPH ADJACENCY MATRIX
      for (j in 1:(n.node-1)) {
        for (j1 in (j+1):n.node) {
          k0 <- k0+1
          # print(c(j, j1, Nd[j], Nd[j1], length(y[node==Nd[j]]), length(y[node==Nd[j1]])))
          pvalue.t <- t.test(x=y[node==Nd[j]], y = y[node==Nd[j1]], 
                             alternative="two.sided", var.equal = FALSE)$p.value 
          if (pvalue.t <= alpha0) {D[k0, j] <- D[k0, j1] <- 0}
        }
      }
    }
  }
  row.rm <- which(rowSums(abs(D))==0); # print(length(row.rm))
  if (length(row.rm)>0) D <- D[-row.rm, ]  # REMOVE ROWS OF ALL 0s
  if (adaptive) {
    beta.hat <- coef(fit0.anova) %>% as.numeric()
    if (length(row.rm)>0) beta.hat <- beta.hat[-row.rm]
    # d1 <- as.vector(abs(D%*%beta.hat))
    d1 <- abs(diff(beta.hat))   # BETTER ALTERNATIVE
    if (max(d1) == min(d1)) stop("All beta-hats are the same already.")
    d1 <- (d1-min(d1)) / (max(d1)-min(d1))  # SCALING
    # THRESHOLD TO AVOID NUMERICAL DIFFICULTY WITH FUZED LASSO
    d1[d1 <= m0] <- m0; # d1[d1>= M0] <- M0
    d0 <- 1/d1
    D <- diag(d0) %*% D
  }
  
  fit <-  fusedlasso(y, X, D=D, gamma=0, approx=TRUE, minlam=0, 
                     rtol = 1e-15, btol = 1e-15, eps = 1e-10)
  if (return.X) return(list(fit=fit, X=X))
  else return(fit)
}


# ==============================================================
# ADD 0 COLUMNS TO A DATA FRAME OR MATRIX TO SPECIFIC POSITIONS
# ==============================================================

add0.columns <- function(x, columns){
  if (!is.data.frame(x) && !is.matrix(x)) stop("Error: x should be a matrix or data frame.")
  if (is.matrix(x)) x <- as.data.frame(x)
  if (min(columns) <1) stop("Error: columns must be a vector of integers >= 1.")
  p <- NCOL(x)
  for (j in 1:length(columns)) {
    pos.col <- columns[j]
    # print(cbind(j, pos.col, p))
    if (pos.col <= p) {x %>% add_column(0, .before=pos.col) -> x}
    p <- NCOL(x)
  }  
  for (j in rev(columns)) {if (j > p)  x<- x %>% add_column(0, .after=p)}  
  x <- data.matrix(x)  # RETURN A MATRIX
  return(x)
}
# add0.columns(x=matrix(1:10, 2, 5), columns=c(1, 3, 5, 8, 9))




# =======================================================
# WRAP-UP FUNCTION FOR TEST-SAMPLE BASED CART VIA rpart
# =======================================================
cart.test <- function(formula, data, test, n0=NROW(data)/20, 
                       plot.it=TRUE, save.subtrees=FALSE, details=FALSE){
  yname <- all.vars(formula)[1]
  y <- dat[, yname]
  tre0 <- cart0(formula, data=data, n0=n0)
  CP <- tre0$cptable %>% as.data.frame %>% select(CP) %>% unlist() %>% as.numeric()
  nsplit <- tre0$cptable %>% as.data.frame %>% select(nsplit) %>% unlist() %>% as.numeric()
  nsubtree <- length(CP)
  SSE <- rep(0, nsubtree)
  if (save.subtrees) TREE <- list(1:nsubtree)
  y.test <- test[, yname]
  for (i in 1:nsubtree) {
    tree.i <- rpart::prune(tre0, cp=CP[i])
    if (save.subtrees) TREE[[i]] <- tree.i
    yhat.i <- predict(tree.i, newdata=test)
    SSE[i] <- sum((y.test-yhat.i)^2)
  }
  n1 <- NROW(test)
  AIC <- n1*log(SSE) + 2 * nsplit
  BIC <- n1*log(SSE) + log(n1) *nsplit
  if (plot.it) {
    par(mfrow=c(1,3), mar=rep(4,4))
    plot(nsplit, SSE, type="b", col="blue", xlab="# splits")
    plot(nsplit, AIC, type="b", col="blue", xlab="# splits")
    plot(nsplit, BIC, type="b", col="blue", xlab="# splits")
  }
  i.star <- c(which.min(SSE), which.min(AIC), which.min(BIC))
  Best.cp <- CP[i.star]
  n.splits <- nsplit[i.star];n.splits
  # TO HANDLE ERROR WITH LIST OF rpart TREES
  btree.SSE <- rpart::prune(tre0, cp=Best.cp[1])  
  btree.AIC <- rpart::prune(tre0, cp=Best.cp[2])  
  btree.BIC <- rpart::prune(tre0, cp=Best.cp[3])  
  
  names(Best.cp) <- names(n.splits) <- c("SSE", "AIC", "BIC")
  if (save.subtrees)  {
    Btree <- TREE[i.star]
    names(Btree) <- c("SSE", "AIC", "BIC")
    return(list(tre0=tre0, CP=CP, nsplit=nsplit, TREE=TREE, 
              Best.cp=Best.cp, n.splits=n.splits, Btree=Btree,
              btree.SSE=btree.SSE, btree.AIC=btree.AIC, btree.BIC=btree.BIC))
  } else {
    return(list(tre0=tre0, CP=CP, nsplit=nsplit, 
                Best.cp=Best.cp, n.splits=n.splits,
                btree.SSE=btree.SSE, btree.AIC=btree.AIC, btree.BIC=btree.BIC))
  }
    
}




##################################
# V-FOLD CROSS-VALIDATION
##################################

# ===========================================================
# WRAP-UP FUNCTION FOR 1-SE V-FOLD CV BASED CART VIA rpart
# ===========================================================

cart.CV <- function(formula, data, method="anova", control=NULL, V=10,
                 size.selection=c("0SE", "1SE"), plot.it=FALSE){
  if (is.null(control)) control <- rpart.control(minsplit=20, minbucket=10, 
                                                 maxdepth=8, cp=0, maxcompete=0,  # NUMBER OF COMPETITIVE SPLITS
                                                 maxsurrogate=0, usesurrogate=2, surrogatestyle=0,  	# SURROGATE SPLITS FOR MISSING DATA
                                                 xval=V)  
  tre0 <- rpart(formula=formula, data=data, method=method, control=control); 
  if (plot.it) plotcp(tre0, minline = TRUE) # 1SE
  if (size.selection=="0SE") {
    opt <- which.min(tre0$cptable[,"xerror"])
    best.cp <- tre0$cptable[opt, "CP"]; # print(cp.best) 
    best.tree <- prune(tre0, cp = best.cp)
  } else if (size.selection=="1SE") {
    cv.error <- (tre0$cptable)[,4]
    SE1 <- min(cv.error) + ((tre0$cptable)[,5])[which.min(cv.error)]      # 1SE; CAN BE EASILY MODIFIED AS aSE FOR SOME a
    position <- min((1:length(cv.error))[cv.error <= SE1]); # print(position)
    # n.size  <- (tre0$cptable)[,2] + 1  #TREE SIZE IS ONE PLUS NUMBER OF SPLITS. 
    # best.size <- n.size[position]; # best.size
    best.cp <-  sqrt(tre0$cptable[position,1]*tre0$cptable[(position-1),1]); # print(best.cp)
    # best.cp <- tre0$cptable[position,1]; print(best.cp)
    best.tree <- prune(tre0, cp=best.cp)
  }
  else stop("The values of size.selection= must be either 0SE or 1SE")
  btree.size <- sum(best.tree$frame=="<leaf>")
  return(list(btree=best.tree, cp=best.cp, btree.size=btree.size, tree0=tre0))
}







# =============================================
# FUNCTION Grow.fusedTree() GROWS A FUSEDTREE
# =============================================

Grow.fusedTree <- function(formula, dat1, 
                           Lambda=seq(0, 50, length.out = 200),
                                n0=NROW(data)/20, max.leaves = NULL,  # PRUNING
                                adaptive=TRUE,
                                ordinalize=TRUE,  # ORDINALIZE NODE?
                                modify.D=TRUE, alpha0=0.0001, # MODIFY GRAPH ADJACENCY MATRIX D?
                                add0Lambda=TRUE, digits=4, 
                                details=FALSE)
{
  # CONSTRUCT FUSED TREES WITH dat1
  # ---------------------------------
  yname <- all.vars(formula)[1]
  y <- dat1[, yname]
  tre0 <- cart0(formula, data=dat1, n0=n0, max.leaves =max.leaves)
  dat <- send.down(tre0, data=dat1)
  node.train <- dat$node
  node.levels <- unique(node.train)
  Node.levels <- n.nodes <- length(node.levels) 
  
  if (ordinalize) {
    # USE INFO IN y
    # aggregate(y, by=list(node.train), FUN=mean) %>% 
    #   arrange(x) %>% 
    #   select(Group.1) %>% 
    #   t() %>% as.character -> node.levels
    tre0$frame %>% 
      mutate(nd=row.names(.), row=1:n()) %>% 
      subset(var=="<leaf>") %>% 
      arrange(yval) %>%   select(nd, row) -> Node.levels # MORE DETAILS
    Node.levels %>%
      select(row) %>% t() %>% as.character -> node.levels 
    if (details) { print(node.levels)
                print(tre0$frame[tre0$frame$var == "<leaf>", ])
                print(node.levels)}
    node.train <- ordered(node.train, levels=node.levels)
  }
  if (details) 
    print(aggregate(y, list(node.train), FUN=mean))
  out.Flasso <- Flasso(y=y, node=node.train, adaptive=adaptive,
                       modify.D=modify.D, alpha0=alpha0, return.X=TRUE)
  fit.Flasso <- out.Flasso$fit
  if (details) {
    print(" ===================================")
    print(fit.Flasso$lambda)
    print(fit.Flasso$beta)
    plot(fit.Flasso, type = "primal") 
  }
  
  if (is.null(Lambda)) {
    Y.fitted <- fit.Flasso$fit
    Lambda0 <- fit.Flasso$lambda
    if (add0Lambda) {
      Lambda0 <- c(Lambda0, 0)
      Y.fitted <- cbind(Y.fitted, 0)   ###
    }
    Lambda <- Lambda0 
  } else {
    Y.fitted <- predict.genlasso(fit.Flasso, Xnew=out.Flasso$X, 
                                 lambda = Lambda)$fit
  }
  # n.lambda <- length(Lambda)
  dimnames(Y.fitted)[[2]] <- Lambda %>% round(digits=digits)  
  # print(dimnames(Y.fitted))
  if (is.element(0, Lambda)) {
    Y.fitted[, "0"] <- predict(tre0)   
  }
  Y.fitted <- Y.fitted %>% round(digits=digits)
  # apply(Y.fitted, 2, FUN = function(x) length(unique(x)))
  return(list(tre0=tre0, fit.Flasso=fit.Flasso, ordinalize=ordinalize,
              node.train=node.train, node.levels=node.levels,  Node.levels=Node.levels,
              Lambda=Lambda, Y.fitted=Y.fitted, 
              digits=digits, formula=formula, y=y))
    
}
  
# =======================================================  
# FUNCTION Apply.fusedTree() APPLY FUSED TREES TO dat2
# =======================================================

Apply.fusedTree <- function(fit.Grow.fusedTree, dat2, Lambda=NULL, details=FALSE)  
{
  # EXTRACT RESULTS FROM fit.Grow.fusedTree
  tre0 <- fit.Grow.fusedTree$tre0;
  fit.Flasso <- fit.Grow.fusedTree$fit.Flasso
  ordinalize <- fit.Grow.fusedTree$ordinalize
  node.levels <- fit.Grow.fusedTree$node.levels
  n.nodes <- length(node.levels)
  Node.levels <- fit.Grow.fusedTree$Node.levels
  Y.fitted <- fit.Grow.fusedTree$Y.fitted
  digits <- fit.Grow.fusedTree$digits
  form <- fit.Grow.fusedTree$formula
  yname <- all.vars(form)[1]
  y <- fit.Grow.fusedTree$y
  Lambda0 <- fit.Grow.fusedTree$Lambda
  if (is.null(Lambda)) Lambda <- Lambda0
  n.lambda <- length(Lambda)
  
  # APPLY TO dat2
  test0 <- send.down(tre0, data=dat2)
  n.nodes.test <- length(unique(test0$node))
  if (ordinalize) test0$node <- ordered(test0$node, levels=node.levels)
  else test0$node <- factor(test0$node, levels=node.levels)
  X0.test <- model.matrix(~factor(node)-1, data=test0)
  X.test <- X0.test
  # print(c(n.nodes, n.nodes.test))
  # print(node.levels)
  if (n.nodes != n.nodes.test) { 
    # A KNOTTY PROBLEM - NEEDS TO ADD 0 COLUMNS TO X.test IN THIS CASE
    X0.test <- as.data.frame(X0.test)
    if (details) print(c(n.nodes, n.nodes.test))
    tbl0 <- table(test0$node)
    cols.add0 <- which(tbl0==0)
    X.test <- add0.columns(X0.test, cols.add0)
  } 
  Y.hat <- predict.genlasso(fit.Flasso, Xnew=X.test, 
                            lambda = Lambda)$fit %>% round(digits=digits)
  dimnames(Y.hat)[[2]] <- Lambda %>% round(digits=digits)  
  
  
  # COMPUTE PREDICTED VALUES AND PREDICTION SSE WITH Lambda
  # --------------------------------------------------------
  y.test <- dat2[, yname]
  n1 <- NROW(dat2); 
  SSE <- N.grp <- rep(0, n.lambda)
  Fits.anova <- as.list(1:n.lambda)
  YHAT <- matrix(0, n1, n.lambda)
  for (k in 1:n.lambda){
    lambda.k <- Lambda[k]
    # print(cbind(k, Lambda[k]))
    grp <- Y.fitted[,k]
    if (lambda.k==0) {
      grp.test <- predict(tre0, newdata=dat2, type="vector") %>% round(digits=digits)
    } else grp.test <- Y.hat[,k]
    # print(unique(grp)); print(unique(grp.test)); 
    if (sum(!is.element(unique(grp.test), unique(grp))) > 0) 
      stop("Hmmm! Check groups of predicted/fitted values.")
    n.grp <- length(unique(grp))
    N.grp[k] <- n.grp
    if (n.grp<=1) {fit.k <- lm(y~1)} else {fit.k <- lm(y~factor(grp)-1)}
    Fits.anova[[k]] <- fit.k
    yhat <- predict.lm(fit.k, newdata=data.frame(grp=grp.test))  
    YHAT[, k] <- t(yhat)
    # if (n.grp<=1) {fit1.k <- lm(y.test~1)} 
    # else {fit1.k <- lm(y.test~factor(Y.hat[,k])-1)}
    # yhat <- fitted(fit1.k)     # USING TEST DATA SOLELY
    SSE[k] <- sum((y.test-yhat)^2)
  }
  colnames(YHAT) <- Lambda
  out <- list(formula=form, digits=digits, tre0=tre0, fit.Flasso=fit.Flasso, 
              Lambda=Lambda, N.grp=N.grp, ordinalize=ordinalize, node.levels=node.levels, 
              Node.levels=Node.levels, Fits.anova=Fits.anova, 
              PreSSE=SSE, YHAT=YHAT)  
  return(out)
}




# ============================================================
# WRAPPER FOR fusedTree WITH BOTH TEST SAMPLE AND V-FOLD CV
# ============================================================

fusedTree <- function(formula, data, test=NULL, selection=c("TestSample", "CV"),
    n0=NROW(data)/20, max.leaves = NULL,  # ANY PRUNING?
    ordinalize=TRUE,  # ORDINALIZE NODE?
    adaptive=TRUE,
    modify.D=TRUE, alpha0=0.0001, # MODIFY GRAPH ADJACENCY MATRIX D?
    add0Lambda=TRUE, digits=4,
    nfold=10,                  # V-FOLD CV; V=nfold
    plot.it=TRUE, x.plot="#groups", 
    details=FALSE)
{
  Grow0 <- Grow.fusedTree(formula=formula, dat1=data, Lambda=NULL,
                          n0=n0, max.leaves = max.leaves,  # PRUNING
                          adaptive=adaptive, # ADAPTIVE FUSED LASSO
                          ordinalize=ordinalize,  # ORDINALIZE NODE?
                          modify.D=modify.D, alpha0=alpha0, # MODIFY GRAPH ADJACENCY MATRIX D?
                          add0Lambda=add0Lambda, digits=digits, 
                          details=details)
  # yname <- all.vars(formula)[1]
  tre0 <- Grow0$tre0;
  fit.Flasso <- Grow0$fit.Flasso
  Y.fitted <- Grow0$Y.fitted
  N.grp <- apply(Y.fitted, 2, n_distinct)
  Lambda0 <- Grow0$Lambda
  n.lambda <- length(Lambda0)
  node.levels <- Grow0$node.levels
  n.nodes <- length(node.levels) 
  Node.levels <- Grow0$Node.levels
  
  if (selection=="TestSample") {
    # ===================
    # TEST SAMPLE METHOD
    # ===================
    if (is.null(test)) stop("Error: No test sample for the test sample method.")
    out0 <- Apply.fusedTree(fit.Grow.fusedTree=Grow0, dat2=test, Lambda=NULL, 
                            details=details)
    SSE <- out0$PreSSE
    Fits.anova <- out0$Fits.anova
    n1 <- NROW(test)}
  else if (selection=="CV") {
    # ===================
    # V-FOLD CV METHOD
    # ===================    
    if (!is.null(test)) print("You have got a test sample, but want to do V-fold CV. Are you sure?")
    y <- Grow0$y
    # OBTAIN Fits.anova
    Fits.anova <- as.list(1:n.lambda)
    for (k in 1:n.lambda){
      grp <- Y.fitted[,k]
      if (N.grp[k] <=1) {fit.k <- lm(y~1)} else {fit.k <- lm(y~factor(grp)-1)}
      Fits.anova[[k]] <- fit.k
    }
    
    # THE LOOP FOR V-FOLD CV
    n <- n1 <- NROW(data)
    id.cv <- sample(1:nfold, size=n, replace=TRUE) 
    SSE <- 0
    for (v in 1:nfold){
      train.v <- dat %>% subset(id.cv!=v)
      test.v <- dat %>% subset(id.cv==v)
      grow.v <- Grow.fusedTree(formula=formula, dat1=train.v, Lambda=Lambda0,
                              n0=n0, max.leaves = max.leaves,  # PRUNING
                              adaptive=adaptive,
                              ordinalize=ordinalize,  # ORDINALIZE NODE?
                              modify.D=modify.D, alpha0=alpha0, # MODIFY GRAPH ADJACENCY MATRIX D?
                              add0Lambda=add0Lambda, digits=digits, 
                              details=details)
      out.v <- Apply.fusedTree(fit.Grow.fusedTree=grow.v, dat2=test.v, 
                               Lambda=Lambda0, details=details)
      SSE <- SSE + out.v$PreSSE
    }
  } else stop("Error: selection must be either TestSample or CV.")

  
  # FIND OUT THE BEST LAMBDAS
  AIC <- n1*log(SSE) + 2 * N.grp
  BIC <- n1*log(SSE) + log(n1) *N.grp
  if (plot.it) {
    if (x.plot == "lambda") {
      par(mfrow=c(1,3))
      plot(Lambda0, SSE, type="b", col="blue", xlab="lambda")
      plot(Lambda0, AIC, type="b", col="blue", xlab="lambda")
      plot(Lambda0, BIC, type="b", col="blue", xlab="lambda")
    }
    else if (x.plot=="#groups") {
      par(mfrow=c(1,3))
      plot(N.grp, SSE, type="b", col="blue", xlab="#groups")
      plot(N.grp, AIC, type="b", col="blue", xlab="#groups")
      plot(N.grp, BIC, type="b", col="blue", xlab="#groups")
    } 
    else stop("Error: Wrong Specification of x.plot.")
  }
  K.star <- c(which.min(SSE), which.min(AIC), which.min(BIC))
  Best.lambda <- Lambda0[K.star]
  n.groups <- N.grp[K.star]
  groups.best <- Y.fitted[, K.star] 
  Fits.anova.best <- Fits.anova[K.star] # ANOVA FITS

  # LEAF-COLORING FOR TREE
  node.train <- Grow0$node.train
  data.frame(node=node.train, groups.best) %>% 
      distinct() %>% 
      arrange(node) %>% 
      setNames(c("node", "SSE", "AIC", "BIC")) %>% 
      select(SSE, AIC, BIC) %>%
      mutate(across(everything(), ~ as.numeric(factor(.x)))) %>% 
      as.matrix -> Color.grp
    Colors <- matrix(NA, nrow =NROW(tre0$frame), ncol=3) %>% 
      as.data.frame() %>% 
      setNames(c("SSE", "AIC", "BIC"))
    Colors[which(tre0$frame=="<leaf>"), ] <- Color.grp
    
    out <- list(formula=formula, digits=digits, tre0=tre0, fit.Flasso=fit.Flasso, 
                Best.lambda=Best.lambda, ordinalize=ordinalize, digits=digits,
                node.levels=node.levels, Node.levels=Node.levels,
                Fits.anova=Fits.anova, Fits.anova.best=Fits.anova.best,
                n.groups=n.groups, Colors=Colors+1)  # TO AVOID BLACK COLOR (=1)
    return(out)
}
  



# =========================================================================
# shear() TRUNCATES TREEFUL INTO THE FINAL TREE AND PREPARES LEAVE COLORS
# =========================================================================

shear <- function(tre, paint, details=FALSE) {
  nd <- row.names(tre$frame) %>% as.integer
  if (details) print(cbind(nd, paint))
  which.desc <- rpart:::descendants(nodes=as.numeric(nd), include = FALSE)
  Nd <- t(replicate(n=length(nd), nd)) 
  # Descendant <- Nd * which.desc
  # which(nd[rowSums(Descendant) !=0])   # INTERNAL NODES
  Paint <- t(replicate(n=length(paint), paint)) 
  leaf.color <- Paint*which.desc
  leaf.color[leaf.color==0] <- NA
  ncolor.leaf <- apply(leaf.color, 1, FUN=function(x) length(unique(na.omit(x))))
  n.nodes.trunc <- sum(ncolor.leaf==1)
  if (details) {print(ncolor.leaf); print(n.nodes.trunc)}
  if (!is.na(n.nodes.trunc) & n.nodes.trunc >= 1) { 
    nd.truncate <- nd[ncolor.leaf==1]
    color.replace <- apply(leaf.color[ncolor.leaf==1, ], 1, FUN=function(x) unique(na.omit(x)))
    which.nd.2b.trunc <- which(colSums(which.desc[ncolor.leaf==1, ncolor.leaf==1])  == 0)
    nd.truncate <- nd.truncate[which.nd.2b.trunc]
    color.replace <- color.replace[which.nd.2b.trunc]
    if(details) {print("===== NODES TO BE TRUNCATED: ====== "); print(nd.truncate)}
    #
    tre <- snip.rpart(tre, toss = nd.truncate)  # THE TRUNCATED TREE
    color.new <- paint
    color.new[is.element(nd, nd.truncate)] <- color.replace
    nd.new <- row.names(tre$frame) %>% as.integer  
    paint <- color.new[is.element(nd, nd.new)]  # THE NEW COLOR
  }
  # PREPARE THE COLOR PALLETTE
  n.clrs <- length(unique(na.omit(paint)))
  palet <- viridis_pal(alpha = 1, begin = 0.2, end = 1, option="D")(n.clrs)
  color.tree <- palet[paint-1] %>% replace_na("gray50")
  return(list(tre=tre, paint=paint, color.tree=color.tree))
}
 


# ==============================
# PLOT TREE WITH COLORED LEAVES
# ==============================
plot.tree <- function(tree, color) {
  prp(tree, extra = 1, faclen=0, nn = T, roundint=FALSE,
      fallen.leaves=FALSE, type=0, 
      branch.col="darkgray",
      split.col="gray30", nn.col="gray45",
      # split.box.col = "lightgray", 
      # split.border.col = "black",
      # shadow.col = "gray", col="red",
      box.col=color) 
}



# ===============================
# PREDICT WITH A NEW DATA SET
# ===============================

predict.fusedTree <- function(fit.fusedTree, newdata){
  tre0 <- fit.fusedTree$tre0
  fit.Flasso <- fit.fusedTree$fit.Flasso
  Fits.anova.best <- fit.fusedTree$Fits.anova.best
  n.nodes <- length(fit.fusedTree$node.levels) 
  Best.lambda <- fit.fusedTree$Best.lambda
  digits <- fit.fusedTree$digits
  
  test <- send.down(tre0, data=newdata)
  if (fit.fusedTree$ordinalize) test$node <- ordered(test$node, levels=fit.fusedTree$node.levels)
  node.pred <- test$node  #################
  X0.test <- model.matrix(~factor(node)-1, data=test)
  X.test <- X0.test
  n.nodes.test <- length(unique(test$node))
  if (n.nodes != n.nodes.test) { 
    tbl0 <- table(test$node)
    cols.add0 <- which(tbl0==0)
    X.test <- add0.columns(X0.test, cols.add0)
  } 
  # print(Best.lambda)
  
  yhat.fusedtree <- yhat.flasso <- matrix(0, NROW(newdata), length(Best.lambda)) 
  for (j in 1:length(Best.lambda)) {
    fit.j <- Fits.anova.best[[j]]
    lamb <- Best.lambda[j]
    if (lamb==0) yhat.flasso[,j] <- predict(tre0, newdata, type="vector") %>% round(digits=digits)
    else yhat.flasso[,j] <- predict.genlasso(fit.Flasso, Xnew=X.test, 
                                      lambda = lamb)$fit %>% round(digits=digits)
    yhat.fusedtree[, j] <- predict.lm(fit.j, newdata=data.frame(grp=yhat.flasso[,j]))
  }
  # PREPARE THE OUTPUT
  group.pred <- yhat.flasso %>% 
    as.data.frame() %>% 
    setNames(c("SSE", "AIC", "BIC")) %>% 
    mutate(across(everything(), ~ as.numeric(factor(.x)))) 
  
  yhat.fusedtree <- yhat.fusedtree %>% 
    as.data.frame() %>% 
    setNames(c("SSE", "AIC", "BIC")) 
  return(list(node.pred=node.pred, group.pred=group.pred, 
              yhat.flasso=yhat.flasso, yhat.fusedtree=yhat.fusedtree))
}





#