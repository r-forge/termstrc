rm(list = ls())
library("Rcpp")
library("inline")

library(termstrc)
data(govbonds)
bdata <- prepro_bond("GERMANY",govbonds,c(0,10))
beta <- rep(1,4)
tau <- c(1,10)
m <- bdata$m[[1]]
cf <- bdata$cf[[1]]
w <- bdata$duration[[1]][,3]
p <- bdata$p[[1]]

objfct_sv_bonds_grid <- function(beta, tau, m, cf, w, p) {
      bsv <- c(beta[1:3],tau[1],beta[4],tau[2])
      phat <- bond_prices("sv",bsv,m,cf)$bond_prices
      loss_function(p, phat, w)
    }

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
                i = 0;
		while (i<mc.rows() && mc(i,j)>0) {
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

grad_sv_bonds_gridCpp <- '
	RcppVector<double> betac(beta);
	RcppVector<double> tauc(tau);
	RcppMatrix<double> mc(m);
	RcppMatrix<double> cfc(cf);
	RcppVector<double> wc(w);
	RcppVector<double> pc(p);
	
	const int N = mc.rows();
	const int M = mc.cols();
	
	RcppMatrix<double> emt1(N, M);
	RcppMatrix<double> emt2(N, M);
	RcppMatrix<double> t1emt1(N, M);
	RcppMatrix<double> emt1tm(N, M);
	RcppMatrix<double> emt2tm(N, M);
	RcppMatrix<double> acf(N, M);
	RcppMatrix<double> dm(N, M);
	RcppVector<double> csacf(M);
	RcppVector<double> csdm(M);
	RcppVector<double> csdt1emt1(M);
	RcppVector<double> csdmem1tm(M);
	RcppVector<double> csdmem2tm(M);
	RcppVector<double> b(M);
	Rcpp::NumericVector gbeta(4);
	
	int i = 0;
	int j = 0;
		
	for (j = 0; j<M; j++) {
		i = 0;
		while (i<N && mc(i,j)>0) {
			emt1(i,j) = exp(-mc(i,j)/tauc(0));
			emt2(i,j) = exp(-mc(i,j)/tauc(1));
			t1emt1(i,j) = tauc(0)*(1 - emt1(i,j));
			emt1tm(i,j) = - emt1(i,j) + t1emt1(i,j)/mc(i,j);
			emt2tm(i,j) = - emt2(i,j) + tauc(1)*(1 - emt2(i,j))/mc(i,j);					   
			acf(i,j) = cfc(i,j)*exp((-betac(0) - betac(2)*emt1tm(i,j) - betac(3)*emt2tm(i,j) - (betac(1)*t1emt1(i,j))/mc(i,j))*mc(i,j)/100);
			dm(i,j) = acf(i,j)/100*mc(i,j);
			csacf(j) += acf(i,j);
			csdm(j) += dm(i,j) ;
			csdt1emt1(j) += acf(i,j)/100*t1emt1(i,j);
			csdmem1tm(j) += dm(i,j)*emt1tm(i,j);
			csdmem2tm(j) += dm(i,j)*emt2tm(i,j);
			i++;
		}
		b(j) = -2*wc(j)*(pc(j) - csacf(j));
		gbeta[0] += b(j)*(-csdm(j));
		gbeta[1] += b(j)*(-csdt1emt1(j));
		gbeta[2] += b(j)*(-csdmem1tm(j));
		gbeta[3] += b(j)*(-csdmem2tm(j));
	}
	
	return gbeta;
        '
grad_sv_bonds_gridC <- cfunction(signature(beta="numeric", tau="numeric", m="numeric", cf="numeric", w="numeric", p="numeric"), grad_sv_bonds_gridCpp, Rcpp=TRUE, verbose=FALSE)


objfct_sv_bonds_grid(beta, tau, m, cf, w, p)
objfct_sv_bonds_gridC(beta, tau, m, cf, w, p)
grad_sv_bonds_grid(beta, tau, m, cf, w, p)
grad_sv_bonds_gridC(beta, tau, m, cf, w, p)

control = list()
outer.iterations = 30
outer.eps = 1e-04

ui <- rbind(c(1,0,0,0),                 # beta0 > 0
            c(1,1,0,0))                 # beta0 + beta1 > 0
ci <- c(0,0)


lsparam <- constrOptim(theta = rep(0.01,4),
                       f = objfct_sv_bonds_grid,
                       grad = NULL,
                       ui = ui,
                       ci = ci,
                       mu = 1e-04,
                       control = control,
                       method = "Nelder-Mead",
                       outer.iterations = outer.iterations,
                       outer.eps = outer.eps,
                       tau, m, cf, w, p) ## additional inputs for f and grad


lsparam2 <- constrOptim(theta = rep(0.01,4),
                        f = objfct_sv_bonds_grid,
                        grad = grad_sv_bonds_grid,
                        ui = ui,
                        ci = ci,
                        mu = 1e-04,
                        control = control,
                        method = "BFGS",
                        outer.iterations = outer.iterations,
                        outer.eps = outer.eps,
                        tau, m, cf, w, p) ## additional inputs for f and grad

lsparam3 <- constrOptim(theta = rep(0.01,4),
                       f = objfct_sv_bonds_gridC,
                       grad = grad_sv_bonds_gridC,
                       ui = ui,
                       ci = ci,
                       mu = 1e-04,
                       control = control,
                       method = "BFGS",
                       outer.iterations = outer.iterations,
                       outer.eps = outer.eps,
                       tau, m, cf, w, p) ## additional inputs for f and grad

## TESTING

N <- 20

test1 <- function() {
  for (i in 1:N) {
    lsparam <- constrOptim(theta = rep(0.01,4),
                       f = objfct_sv_bonds_grid,
                       grad = NULL,
                       ui = ui,
                       ci = ci,
                       mu = 1e-04,
                       control = control,
                       method = "Nelder-Mead",
                       outer.iterations = outer.iterations,
                       outer.eps = outer.eps,
                       tau, m, cf, w, p) ## additional inputs for f and grad
  }
}

test2 <- function() {
  for (i in 1:N) {
    lsparam2 <- constrOptim(theta = rep(0.01,4),
                        f = objfct_sv_bonds_grid,
                        grad = grad_sv_bonds_grid,
                        ui = ui,
                        ci = ci,
                        mu = 1e-04,
                        control = control,
                        method = "BFGS",
                        outer.iterations = outer.iterations,
                        outer.eps = outer.eps,
                        tau, m, cf, w, p) ## additional inputs for f and grad
  }
}

test3 <- function() {
  for (i in 1:N) {
    lsparam <- constrOptim(theta = rep(0.01,4),
                       f = objfct_sv_bonds_gridC,
                       grad = NULL,
                       ui = ui,
                       ci = ci,
                       mu = 1e-04,
                       control = control,
                       method = "Nelder-Mead",
                       outer.iterations = outer.iterations,
                       outer.eps = outer.eps,
                       tau, m, cf, w, p) ## additional inputs for f and grad
  }
}

test4 <- function() {
  for (i in 1:N) {
    lsparam2 <- constrOptim(theta = rep(0.01,4),
                        f = objfct_sv_bonds_gridC,
                        grad = grad_sv_bonds_gridC,
                        ui = ui,
                        ci = ci,
                        mu = 1e-04,
                        control = control,
                        method = "BFGS",
                        outer.iterations = outer.iterations,
                        outer.eps = outer.eps,
                        tau, m, cf, w, p) ## additional inputs for f and grad
  }
}

test5 <- function() {
  for (i in 1:N) {
    lsparam2 <- constrOptim(theta = rep(0.01,4),
                        f = objfct_sv_bonds_gridC,
                        grad = grad_sv_bonds_grid,
                        ui = ui,
                        ci = ci,
                        mu = 1e-04,
                        control = control,
                        method = "BFGS",
                        outer.iterations = outer.iterations,
                        outer.eps = outer.eps,
                        tau, m, cf, w, p) ## additional inputs for f and grad
  }
}

system.time(test1())
system.time(test2())
system.time(test3())
system.time(test4())
system.time(test5())


objfct_sv_bonds_grid2 <- function(beta, tau, m, cf, w, p) {
      bsv <- c(beta[1:3],tau[1],beta[4],tau[2])
      spot_rates <- spr_sv(bsv, m)/100
      spot_rates[is.nan(spot_rates)] <- 0        
      
      phat <- apply(cf*exp(-m*spot_rates), 2, sum)
      sum(w*((p - phat)^2))
    }





Rprof()
test1()
Rprof(NULL)
summaryRprof()

Rprof()
test2()
Rprof(NULL)
summaryRprof()

test3 <- function() {
  for (i in 1:5) {
 grad_sv_bonds_grid(beta, tau, m, cf, w, p)
}
}

Rprof()
test3()
Rprof(NULL)
summaryRprof()

Rprof()
 for(i in 1:10000) grad_sv_bonds_gridC(beta,tau,m,cf,w,p)
 Rprof(NULL)
summaryRprof()

Rprof()
 for(i in 1:10000) grad_sv_bonds_grid(beta,tau,m,cf,w,p)
 Rprof(NULL)
summaryRprof()
