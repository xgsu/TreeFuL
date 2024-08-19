
rm(list=ls(all=TRUE))
source("Functions-TreeFuL.R")

#######################################
# EXAMPLE I:  OBESITY DATA
#######################################

dat <- read.csv("dat-obesity.csv", header = T)
head(dat); dim(dat); names(dat);  

# SOME EXPLOATAION
dat %>%
  group_by(Fam.Hist) %>%
  summarise_at(vars(BMI), list(mean = mean, sd=sd))

# MODEL
form <- BMI~ Gender + Age + Fam.Hist + FAVC  + FCVC + NCP + CAEC + SMOKE + 
  CH2O + SCC + FAF + TUE + CALC + MTRANS


# ===========================
# FIT CART & TREEFUL MODELS
# ===========================
set.seed(123)
fit.TreeFuL <- TreeFuL(formula0=form, dat=dat, n.folds=10, plot.it=TRUE,
              start.rpart="CART.0SE", minsplit=50, minbucket=20, maxdepth=10)
names(fit.TreeFuL)
fit.TreeFuL$cp.btrees    # TUNING PARAMETERS
fit.TreeFuL$size.btrees
TREES <- fit.TreeFuL$TREES; names(TREES)

# OBTAIN THE CART TREES
CART.initial <- TREES$CART.initial   # THE LARGE INITIAL TREE
CART.0SE <- TREES$CART.0SE  # 0SE CART TREE
# CART.1SE <- TREES$CART.1SE  # 1SE CART TREE
plot.tree(CART.0SE)

# OBTAIN TREEFUL TREES
TreeFuL.0SE <- TREES$TreeFuL.0SE  # 0SE TREEFUL TREE
# TreeFuL.1SE <- TREES$TreeFuL.1SE  # 1SE TREEFUL TREE
plot.tree(TreeFuL.0SE)

# FIGURE WITH BOTH CART AND TREEFUL (1SE) TREES
postscript(file="fig03-obisity-tree.eps", horizontal = FALSE)
par(mfrow=c(2, 1), mar=rep(4,4))
plot.tree(tree=CART.1SE)
title(main="(a) CART", line =3)

plot.tree(tree=TreeFuL.1SE)
title(main="(b) TreeFuL", line =3)
dev.off()

names(TreeFuL.1SE)
# plot.tree(tree=TreeFuL.0SE$tre, color = TreeFuL.0SE$paint)


# ============================================
# REPEATED TEST-SAMPLES TO INSPECT INSTABILITY
# ============================================

set.seed(123)
nrun <- 100
n <- NROW(dat)
Bootstrap <- TRUE
OUT <- NULL
VARS <- as.list(1:nrun)
yname <- form[[2]] |> as.character()
for (i in 1:nrun) {
  print(i)
  if (Bootstrap) id.train <- sample(1:n, size=n, replace=TRUE)
  else id.train <- sample(1:n, size=ceiling(n*2/3), replace=FALSE)
  train <- dat[id.train, ]; test <- dat[-id.train, ]
  y.test <- test[, yname]
  fit.i <- TreeFuL(formula0=form, dat=train, n.folds=10, plot.it=FALSE,
              start.rpart="CART.0SE", minsplit=50, minbucket=20, maxdepth=10)
  sizes.i <- fit.i$size.btrees
  TREES.i <- fit.i$TREES;
  btre.CART <- TREES.i$CART.0SE 
  
  yhat.CART <- predict(btre.CART, newdata=test, type = "vector")
  mse.CART <- mean((yhat.CART- y.test)^2) 
  btre.TreeFuL <- TREES.i$TreeFuL.0SE 
  # EXTRACT SPLITTING VARIABLES
  Xs.TreeFuL <- btre.TreeFuL$tre$frame$var  
  Xs.TreeFuL <- Xs.TreeFuL[Xs.TreeFuL!="<leaf>"]
  Xs.TreeFuL <- unique(Xs.TreeFuL)
  VARS[[i]] <- Xs.TreeFuL
  # COMPUTE MSE AND RELATIVE AIC
  yhat.TreeFuL <- predict.TreeFuL(fit.TreeFuL=fit.i, newdata = test)$yhat
  mse.TreeFuL <- mean((yhat.TreeFuL-y.test)^2) 
  AIC.CART <- NROW(test)*log(mse.CART) + 2 * sizes.i[3]
  AIC.TreeFuL <- NROW(test)*log(mse.TreeFuL) + 2 * sizes.i[5]
  AIC.relative  <- (AIC.CART - AIC.TreeFuL) / AIC.CART
  OUT <- rbind(OUT, c(size.cart=sizes.i[3], size.treeful=sizes.i[5], 
                      mse.cart=mse.CART, mse.treeful=mse.TreeFuL, 
                      AIC.relative=AIC.relative))
}
OUT <- OUT %>% as.data.frame() %>% 
  setNames(c("size.cart", "size.treeful", "mse.cart", "mse.treeful", "AIC.relative"))
OUT

apply(OUT, 2, mean)
apply(OUT, 2, sd)
apply(OUT, 2, median)

# OVERLAYED BAR (DENSITY) PLOT  OF TREE SIZES 
dat.size <- OUT %>% 
  select("size.cart", "size.treeful") %>% 
  setNames(c("CART", "TreeFuL")) %>% 
  pivot_longer(cols= everything(), names_to = "Method", values_to="Size")
fig.density <- ggplot(dat.size, aes(x = Size, fill = Method)) + 
  geom_histogram(position = "dodge") +
  scale_fill_manual(values=c("#FF8C69", "#69b3a2")) +
  # geom_density(alpha = 0.5) +
  ylab("Density") + xlab("Final Tree Size") +
  ggtitle('(a) Tree Size') + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

# PARALLEL BOXPLOT OF MSE
dat.mse <- OUT %>% 
  select("mse.treeful", "mse.cart") %>% 
  setNames(c("TreeFuL", "CART")) %>%
  pivot_longer(cols= everything(), names_to = "Method", values_to="MSE") %>% 
  mutate(Method=factor(Method, levels=c("TreeFuL", "CART"), 
                       ordered =TRUE))

fit.boxplot <- ggplot(dat.mse, aes(x=Method, y=MSE, fill=Method)) +
  geom_boxplot(notch = TRUE) +
  scale_fill_manual(values=c("#69b3a2", "#FF8C69")) +
  # coord_flip() + 
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, 
               color="blue", fill="blue") +
  theme(legend.position="none") +
  # xlab("MSE") +
  ggtitle('(b) MSE')  + 
  theme(plot.title = element_text(hjust = 0.42, face = "bold")) 


# DENSITY PLOT of RELATIVE AIC
dat.mse <- OUT %>% 
  select("AIC.relative") %>% 
  ggplot(aes(x = AIC.relative)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1, linetype = 1,
               colour = 2) +
  geom_vline(xintercept = 0, colour = "blue", 
             linetype = 2,)-> fig.AIC


# BAR PLOT OF SPLITTING VARIABLES
VARS0 <- unlist(VARS)  %>%  as.vector() %>% as.data.frame() %>% 
  setNames(c("split.var")) %>% 
  group_by(split.var) %>%   
  mutate(n.count = n()) %>% 
  ggplot(aes(x=reorder(split.var,-n.count))) +
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  xlab("Splitting Variables") + 
  ylab("Frequency") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5)) -> fig.splitvar


# library(ggpubr)
ggarrange(fig.density, fit.boxplot, 
          fig.AIC, fig.splitvar, # labels = c("(a)", "(b)"),
          ncol = 2, nrow = 2)
ggsave(file="fig4-obesity.eps", width = 8, height = 11)






################################
################################
# EXAMPLE II: SIMULATED DATA
################################
################################

rm(list=ls(all=TRUE))
source("Functions-TreeFuL.R")

set.seed(321)
nrun <- 10
N <- c(100, 500, 1000, 5000, 10000)
OUT.time <- matrix(0, 2, length(N))
form <- y~ . - mu 
for (k in 1:length(N)) {
  n <- N[k]
  for (i in 1:nrun) {
    print(cbind(n=n, run=i))
    dat <- rdat(n=n, p=5, model="Tree1", b=5)
    
    # CART
    time.start <- proc.time() 
    control0 <- rpart.control(minsplit=20, minbucket=5, maxdepth=10, 
                cp=0, maxcompete=0, maxsurrogate=0, usesurrogate=2, surrogatestyle=0, xval=10)  
    tre0 <- rpart(formula=form, data=dat, method="anova", control=control0); 
    size.tre0  <- sum(tre0$frame$var=="<leaf>")
    opt <- which.min(tre0$cptable[,"xerror"])
    cp.0SE <- tre0$cptable[opt, "CP"]; # print(cp.0SE) 
    best.tree.0SE <- prune(tre0, cp = cp.0SE)
    size.0SE <- sum(best.tree.0SE$frame$var=="<leaf>")
    cv.error <- (tre0$cptable)[,4]
    SE1 <- min(cv.error) + ((tre0$cptable)[,5])[which.min(cv.error)]      # 1SE; CAN BE EASILY MODIFIED AS aSE FOR SOME a
    position <- min((1:length(cv.error))[cv.error <= SE1]); # print(position)
    # n.size  <- (tre0$cptable)[,2] + 1  #TREE SIZE IS ONE PLUS NUMBER OF SPLITS. 
    # best.size <- n.size[position]; # best.size
    cp.1SE <-  sqrt(tre0$cptable[position,1]*tre0$cptable[(position-1),1]); # print(best.cp)
    # cp.1SE <- tre0$cptable[position,1]; print(cp.1SE)
    best.tree.1SE <- prune(tre0, cp=cp.1SE)
    size.1SE <- sum(best.tree.1SE$frame$var=="<leaf>")
    size.btrees <- c(size.tre0=size.tre0, size.0SE=size.0SE, size.1SE=size.1SE)
    cp.btrees <- c(0, cp.0SE, cp.1SE)
    time.stop <- proc.time() 
    OUT.time[1, k] <- OUT.time[1, k] + sum((time.stop-time.start)[1:2]) 
    
    # TREEFUL
    time.start <- proc.time() 
    fit.TreeFuL <- TreeFuL(formula0=form, dat=dat, n.folds=10, plot.it=FALSE,
                           start.rpart="CART.0SE", minsplit=20, minbucket=5, maxdepth=10)
    time.stop <- proc.time() 
    OUT.time[2, k] <- OUT.time[1, k] + sum((time.stop-time.start)[1:2]) 
  }
}
OUT.time <- OUT.time/nrun
OUT.time
  
write.csv(OUT.time, file="time.csv", row.names = FALSE)  
  


# =========================================== END ===============================================  #

#
