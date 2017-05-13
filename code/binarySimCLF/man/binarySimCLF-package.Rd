\name{binarySimCLF-package}
\alias{binarySimCLF-package}
\alias{binarySimCLF}
\docType{package}
\title{
    A Package For Generating Correlated Binary Data
}
\description{
This software is based on Qaqish (2003).
}
\details{
\tabular{ll}{
Package: \tab binarySimCLF\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-04-20\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
The main routines are \code{blrchk1}, \code{blrchk}, \code{mbsclf}, \code{mbsclf1}, \code{ranXch},
and \code{ranBin2}.  Each contains examples that are pretty sufficient in explaining
how to use the routines.  We warn the user to pay special attention to the documentation
for the functions \code{mbsclf} and \code{mbsclf1} -- especially the part about
checking for CLF compatibility.  It is crucial that this step is taken for
every mean vector and correlation matrix.  Otherwise, the output is nonsense.
}
\author{
Kunthel By and Bahjat F. Qaqish

Maintainer: Kunthel By <kby@bios.unc.edu>
}
\references{
Qaqish, B. F., A family of multivariate binary distributions for
simulating correlated binary variables with specified marginal means and
correlations. \emph{Biometrika} \bold{92}:455-463, 2003.
}
\examples{
#------------------------------------------------------------------------------#

# Simulation under exchangeable correlation matrix.

# range of marginal means is lo:hi
lo = 0.2;
hi = 0.3;

rho = 0.3;      # correlation parameter
n = 11;         # cluster size

m = 1000;        # number of clusters.

# mean vector
mu =(0:(n-1)) * (hi-lo)/(n-1) + lo ;

# exchangeable correlation matrix
R = xch(n,rho);

# Converts R to a covariance matrix.
V = cor2var(R,mu); V

# Intermediate matrix needed for the CLF simulation.
B = allReg(V); B

# checks for CLF compatibility
clf.compat = blrchk1(mu,B); clf.compat

# Generates 1000 correlated response vectors each with mean vector mu and
# correlation matrix R.
if (clf.compat){
  y = mbsclf(m,mu,B)
}

# Sample mean. This should be compared to mu - looks similar.
colMeans(y$y)

# Sample correlation.  This should be compared to R - looks similar.
cor(y$y)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Instead of using the mbsclf function, we can also use the ranXch function.
# This is more efficient for exchangeable correlation.  Consider:
mu = c(0.25, 0.24, 0.26, 0.35);
R = xch(4,0.25);
V = cor2var(R,mu);
# Checks CLF compatibility.
clf.compat = blrchk(mu,V); clf.compat

# Generates 1000 clusters each with mean vector mu and  exchangeable correlation
# parameter rho = 0.25.
if (clf.compat)
{
    y = matrix(rep(NA,4*1000),1000,byrow=TRUE);
    for (i in 1:1000)
    {
        y[i,] <- t( ranXch(mu, R)$y );
    }
}
colMeans(y);
cor(y)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# An example where the mean and the correlation matrix are not CLF compatible.
# This means that there are no multivariate binary distributions in the conditional linear
# family that satisfies the specified mean vector and the correlation matrix.
mu = c(0.25, 0.24, 0.26, 0.35);
R2 = xch(4,0.8); R2
V2 = cor2var(R2,mu); V2
B2 = allReg(V2); B2
clf.compat = blrchk1(mu,B2) #Not CLF compatible
if (!clf.compat){print("Not CLF compatible")} else {
    y2 = mbsclf(200,mu,B2)
}

# Suppose that you didn't check for CLF compatibility and you just went ahead
# and simulated some data based on mu and R2.  If mbsclf() ran without problems
# and returned a data matrix, that data matrix is NONSENSE; i.e., the rows are not iid
# with mean mu and correlation matrix R2.  The most likely scenario however is
# that the simulation routine mbsclf() will fail.
y2 = mbsclf(10,mu,B2); y2

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Unstructured correlation matrix.
r1 = c(1,0.1,0.2);
r2 = c(0.1, 1, 0.3);
r3 = c(0.2, 0.3, 1);
R = rbind(r1,r2,r3);

mu = c(0.4,0.5,0.6);

V = cor2var(R,mu); V
B = allReg(V);

# Checks CLF compatibility.
clf.compat = blrchk1(mu,B)
if (clf.compat)
{
    print("CLF compatible");
    y = mbsclf(100,mu,B);
} else { print("Not CLF compatible"); }


colMeans(y$y);
cor(y$y);
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
}

