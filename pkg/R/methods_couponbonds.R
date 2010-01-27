print.dyncouponbonds <- function(x, ...) {
  cat("This is a dynamic dataset of coupon bonds.\n")
  cat(paste("There are",length(x), "observations between",x[[1]][[1]]$TODAY, "and",x[[length(x)]][[1]]$TODAY,".\n"))
}

print.couponbonds <- function(x, ...) {
  cat("This is a dataset of coupon bonds for:\n")
  cat(names(x),",","\n")
  cat(paste("observed at ", x[[1]]$TODAY,".","\n",sep=""))
 
  }

