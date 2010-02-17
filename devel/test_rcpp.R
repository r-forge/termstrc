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

grad_sv_bonds_gridCpp <- '
	RcppVector<double> betac(beta);
	RcppVector<double> tauc(tau);
	RcppMatrix<double> mc(m);
	RcppMatrix<double> cfc(cf);
	RcppVector<double> wc(w);
	RcppVector<double> pc(p);
        RcppVector<double> phat(pc.size()); // delete later

        RcppMatrix<double> emt1(mc.rows(), mc.cols());
        RcppMatrix<double> emt2(mc.rows(), mc.cols());
        RcppMatrix<double> t1emt1(mc.rows(), mc.cols());
        RcppMatrix<double> emt1tm(mc.rows(), mc.cols());
        RcppMatrix<double> emt2tm(mc.rows(), mc.cols());

        double gbeta1 = 0;
        double gbeta2 = 0;
	double gbeta3 = 0;
        double gbeta5 = 0;

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
        Rcpp::NumericVector xa(4);
        xa[0] = 4;
        xa[3] = 99;
	return xa;
	//return Rcpp::wrap(mse);
        '
grad_sv_bonds_gridC <- cfunction(signature(beta="numeric", tau="numeric", m="numeric", cf="numeric", w="numeric", p="numeric"), grad_sv_bonds_gridCpp, Rcpp=TRUE, verbose=FALSE)

grad_sv_bonds_gridC(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p)

grad_sv_bonds_grid <- function(beta, tau, m, cf, w, p){

  emt1 <- exp(-m/tau[1])
  emt2 <- exp(-m/tau[2])
  t1emt1 <- tau[1]*(1 - emt1)
  emt1tm <- (-emt1 + t1emt1/m)
  emt2tm <- (-emt2 + tau[2]*(1 - emt2)/m)
  

  a <- exp((-beta[1] - beta[3]*emt1tm - beta[4]*emt2tm - (beta[2]*t1emt1)/m)*m/100)

 
  acf <- a*cf
  b <- -2*w*(p-cSums(acf,na.rm=TRUE))
  d <- acf/100
  dm <- d*m
  
  gbeta1 <- sum(b*(-cSums(dm,na.rm=TRUE)))
  gbeta2 <- sum(b*(-cSums(d*t1emt1,na.rm=TRUE)))
  gbeta3 <- sum(b*(-cSums(dm*emt1tm, na.rm=TRUE)))
  gbeta5 <- sum(b*(-cSums(dm*emt2tm, na.rm=TRUE)))
  
  c(gbeta1,gbeta2,gbeta3,gbeta5)
}

cSums <- function (x, na.rm = FALSE, dims = 1L) {
    dn <- dim(x)
    n <- prod(dn[1L:dims])
    dn <- dn[-(1L:dims)]
    z <-  .Internal(colSums(x, n, prod(dn), na.rm))
    z
}

objfct_sv_bonds_grid(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p)
objfct_sv_bonds_gridC(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p)

system.time(for (i in 1:10000) objfct_sv_bonds_gridC(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p))
system.time(for (i in 1:10000) objfct_sv_bonds_grid(beta=beta, tau=tau, m=m, cf=cf, w=w, p=p))
