rm(list = ls())
library("Rcpp")
library("inline")

library(termstrc)
data(govbonds)
bdata <- prepro_bond("GERMANY",govbonds,c(0,10))
beta <- c(1,2,3,4)
tau <- c(1,10)
m <- bdata$m[[1]]
cf <- bdata$cf[[1]]
w <- bdata$duration[[1]][,3]
p <- bdata$p[[1]]

objfct_sv_bonds_gridCpp <- '
	RcppVector<double> betac(beta);
	RcppVector<double> tauc(tau);
	RcppMatrix<double> mc(m);
	RcppMatrix<double> cfc(cf);
	RcppVector<double> wc(w);
	RcppVector<double> pc(p);
	RcppVector<double> phat(pc.size());
	
	int i = 0;
	int j = 0;
	double s;
        double mse = 0;
	
	for (j = 0; j<mc.cols(); j++) {
         // Rprintf("j is %d:\\n", j);
                i = 0;
		while (i<mc.rows() && mc(i,j)>0) {
                //Rprintf("i is %d and mc(i,j) is %f:\\n", i, mc(i,j));
			s =  (betac(0) + betac(1) * ((1 - exp(-mc(i,j)/tauc(0)))/(mc(i,j)/tauc(0))) +
				  betac(2) * (((1 - exp(-mc(i,j)/tauc(0)))/(mc(i,j)/tauc(0))) - exp(-mc(i,j)/tauc(0))) +
				  betac(3) * (((1 - exp(-mc(i,j)/tauc(1)))/(mc(i,j)/tauc(1))) - exp(-mc(i,j)/tauc(1))))/100;

			phat(j) += cfc(i,j)*exp(-s*mc(i,j));
			i++;
		}
	mse += pow(pc(j) - phat(j),2)*wc(j);
	}
	
	return Rcpp::wrap(mse);
        '
objfct_sv_bonds_gridC <- cfunction(signature(beta="numeric", tau="numeric", m="numeric", cf="numeric", w="numeric", p="numeric"), objfct_sv_bonds_gridCpp, Rcpp=TRUE, verbose=FALSE)

objfct_sv_bonds_gridC(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p)


objfct_sv_bonds_grid <- function(beta, tau, m, cf, w, p) {
      bsv <- c(beta[1:3],tau[1],beta[4],tau[2])
      phat <- bond_prices("sv",bsv,m,cf)$bond_prices
      loss_function(p, phat, w)
    }

objfct_sv_bonds_grid(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p)
objfct_sv_bonds_gridC(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p)

system.time(for (i in 1:10000) objfct_sv_bonds_gridC(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p))
system.time(for (i in 1:10000) objfct_sv_bonds_grid(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p))
