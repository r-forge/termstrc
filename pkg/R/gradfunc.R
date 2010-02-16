
grad_sv_bonds_grid <- function(beta, tau, m, cf, w, p){

a <- exp((-beta[1] - beta[3]*(-exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m) - beta[4]*(-exp(-m/tau[2]) + (tau[2]*(1 - exp(-m/tau[2])))/m) - (beta[2]*tau[1]*(1 - exp(-m/tau[1])))/m)*m)


# general term : sum(-2(p-A*cf)*(-apply(A*Cf*m*st,2,sum))*w
#st = specific term 

}



                                                                                                                                    # beta2 specific term
                                                                                                                                    #   beta4
#   1 - exp(-m/beta4)
# beta3 specific term
#  -exp(-m/tau[1]) + (tau[1]*(1 - exp(-m/tau[1])))/m

# beta 4 specific term
#          -(beta[2]/(tau[1]*exp(m/tau[1]))) + (beta[2]*(1 - exp(-m/tau[1])))/m + 
#   beta[3]*(-(1/(tau[1]*exp(m/tau[1]))) + (1 - exp(-m/tau[1]))/m - m/(Power(tau[1],2)*exp(m/tau[1])))

# beta5 specific term
#          -exp(-m/beta6) + (beta6*(1 - exp(-m/beta6)))/m
# beta6 specific term
#beta6
#          -(1/(beta6*exp(m/beta6))) + (1 - exp(-m/beta6))/m - m/(Power(beta6,2)*exp(m/beta6))
