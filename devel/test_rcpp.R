rm(list = ls())
library("Rcpp")
library("inline")

objfct_sv_bonds_grid <- '
	RcppVector<double> betac(beta);
        RcppVector<double> tauc(tau);
        RcppMatrix<double> mc(m);
        RcppMatrix<double> cfc(cf);
        RcppVector<double> wc(w);
        RcppVector<double> pc(p);
        double schas;
        //schas = exp(5.0);
        double mym = 3;
        int i = 0;
        int j = 0;
        double s;

        s =  (betac(1) + betac(2) * ((1 - exp(-mc(i,j)/tauc(0)))/(mc(i,j)/tauc(0))) +
          betac(2) * (((1 - exp(-mc(i,j)/tauc(0)))/(mc(i,j)/tauc(0))) - exp(-mc(i,j)/tauc(0))) +
         betac(3) * (((1 - exp(-mc(i,j)/tauc(1)))/(mc(i,j)/tauc(1))) - exp(-mc(i,j)/tauc(1))));

        double d;
        d = cfc(i,j)*exp(-s*mc(i,j));
        cfc(1,1) = 11;
        Rprintf( "Your value is: %f\\n", cfc(1,1));
	Rprintf( "RcppMatrix: n = %d times k= %d\\n", mc.rows(), mc.cols() ) ;
	Rprintf( "RcppMatrix: elem(1,2) = %f\\n", mc(0,2) ) ;
        return Rcpp::wrap(d);
        '
funx <- cfunction(signature(beta="numeric", tau="numeric", m="numeric", cf="numeric", w="numeric", p="numeric"), objfct_sv_bonds_grid, Rcpp=TRUE, verbose=FALSE)

beta <- c(4,3,2,1)
tau <- c(5,15)
m <- matrix(1:9,ncol= 3)
cf <- matrix(rnorm(100), 10, 10)
w <- rnorm(10)
p <- rnorm(10)

funx(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p)


## objfct_sv_bonds_grid <- function(beta, tau, m, cf, w, p) {
##       bsv <- c(beta[1:3],tau[1],beta[4],tau[2])
##       phat <- bond_prices("sv",bsv,m,cf)$bond_prices
##       loss_function(p, phat, w)
##     }
