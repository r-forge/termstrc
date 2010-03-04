#include <Rcpp.h>

RcppExport SEXP objfct_sv_bonds_gridCpp(SEXP beta, SEXP tau, SEXP m, SEXP cf, SEXP w, SEXP p){

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
}

