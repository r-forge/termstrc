 load("govbonds.RData")
  ns_res <- estim_nss(govbonds, c("GERMANY", "FRANCE", "BELGIUM", "SPAIN"),matrange = c(0,30), method = "ns", deltatau = 0.2)
  print(ns_res)
  plot(ns_res)
  summary(ns_res)

