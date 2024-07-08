#' Plot survival curves
#' @param S1
#' @param S0
#' @param Time
#' @param freq.time Map time interval to coarser intervals
#' @param Label Treatment labels (comparator, treatment)
#' @param Main Title of the plot
#' @export

plotS <- function(S1, S0, Time=seq(0, 70, by=1), freq.time=90,
                  Label=c("Thiazide-like diuretics", "ACE inhibitor"), Main="Cross-fitted TMLE, S(t)"){

  plot(Time, c(1, S0), main=Main, xlab=paste0("Time (", freq.time, " days)"), ylab="S(t)", type="l", col="blue", ylim=c(0.94, 1))
  lines(Time, c(1, S1), type="l", col="red")
  legend("topright", Label,fill=c("blue","red"), box.lty=0)

}










