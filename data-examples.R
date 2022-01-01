
# rm(list=ls(all=TRUE))
source("Functions-FusedTree.R")


####################################
# DATA EXAMPLES
####################################

# BOSTON HOUSING DATA
data(BostonHousing2, package="mlbench")
dat <- BostonHousing2
head(dat); dim(dat); names(dat);  anyNA(dat); 
form <- cmedv~ crim + zn + indus + chas + nox + rm + age + dis + 
      rad + tax + ptratio + b + lstat

# CPUS DATA
data(cpus, package="MASS")
dat <- cpus
head(dat); dim(dat); names(dat);  anyNA(dat); 
form <- log10(perf) ~ syct+mmin+mmax+cach+chmin+chmax


# -------------
# CART 1SE TREE
# -------------

n0 <- 20
cntrl <- rpart.control(minsplit=n0, minbucket=round(n0/3), 
                         maxdepth=8, cp=0, maxcompete=0,  # NUMBER OF COMPETITIVE SPLITS
                         maxsurrogate=0, usesurrogate=2, surrogatestyle=0,  	# SURROGATE SPLITS FOR MISSING DATA
                         xval=10)  

set.seed(123)
fit.cart <- cart.CV(formula=form, data=dat, method="anova", 
                    control=cntrl, V=10, 
                    size.selection="1SE", plot.it=FALSE)
size.cart <- fit.cart$btree.size; size.cart
btre.cart <- fit.cart$btree; # btre.cart

# USING plot.rpart
# plot(btre.cart, compress = TRUE); text(btre.cart, use.n = TRUE)

par(mfrow=c(1, 1), mar=rep(4,4))
plot.tree(tree=btre.cart, color="gray95")


# -------------
# TREEFUL
# -------------

fit.fusedtree <-  fusedTree(formula=form, data=dat, selection="CV",
                            nfold=10, n0=n0, max.leaves = 20)
size.fusedtree <- fit.fusedtree$n.groups[3]; size.fusedtree
c(size.cart, size.fusedtree)

tre0 <- fit.fusedtree$tre0
color.bic <- fit.fusedtree$Colors$BIC

n.clrs <- length(unique(na.omit(color.bic)))
palet.viridis <- viridis_pal(alpha = 1, begin = 0.2, end = 1, option="D")(n.clrs)  # VIRIDIS COLORS
palet.random <- distinctColorPalette(n.clrs)  # RANDOM COLORS
col.tre0 <- palet.viridis[color.bic-1] %>% replace_na("gray50")
col1.tre0 <- palet.random[color.bic-1] %>% replace_na("gray50")

# TRUNCATE TREEFUL 
Truc <- shear(tre=tre0, paint=color.bic, details = FALSE)
tre.final <- Truc$tre
color.final <- Truc$color.tree
par(mfrow=c(2,1))
plot.tree(tree=tre0, color=col.tre0)
plot.tree(tree=tre.final, color=color.final)

# FIGURE WITH BOTH CART AND TREEFUL TREES
postscript(file="fig3-cpus.eps", horizontal = FALSE)
par(mfrow=c(2, 1), mar=rep(4,4))
plot.tree(tree=btre.cart, color="gray95")
title(main="(a) CART", line =3)

plot.tree(tree=tre.final, color=color.final)
title(main="(b) TreeFuL", line =3)
dev.off()



# =========================================== END ===============================================  #

#
