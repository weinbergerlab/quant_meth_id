plot.wave2<-function(cwt.obj) {
  split.screen(rbind(c(0, 1, 0.3, 1), c(0, 1, 0, 0.3)  )) #lr bt
  
  screen(1)
  par(mar=c(0.5, 2, 1, .5))
  image(cwt.obj$time, log2(cwt.obj$periods), t(cwt.obj$local), col=tim.64, breaks=quantile(cwt.obj$local, probs=seq(0, 1, by=1/64)), axes=F, xlab="", ylab="", xlim=(range(cwt.obj$time)), ylim=range(log2(cwt.obj$periods)))
  axis(side=1, at=round(cwt.obj$time)+.5, labels=F)
  axis(side=1, at=seq(ceiling(min(cwt.obj$time)), floor(max(cwt.obj$time)), by=floor(length(cwt.obj$time)*cwt.obj$dt/5))+.5, labels=seq(ceiling(min(cwt.obj$time)), floor(max(cwt.obj$time)), by=floor(length(cwt.obj$time)*cwt.obj$dt/5)), tcl=-.8)
  axis(side=2, at=(-1:floor(log2(max(cwt.obj$periods)))), labels=2^(-1:floor(log2(max(cwt.obj$periods)))))
  contour(cwt.obj$time, log2(cwt.obj$periods), t(cwt.obj$local/cwt.obj$local.sig), add=T, levels=c(1), drawlabels=F, col="black", lwd=2, lty=1)	
  polygon(c(min(cwt.obj$time), min(cwt.obj$time), cwt.obj$time, max(cwt.obj$time), max(cwt.obj$time)), c(max(cwt.obj$periods), min(log2(cwt.obj$coi)),  log2(cwt.obj$coi), min(log2(cwt.obj$coi)), max(cwt.obj$periods)), border=NA, col=hsv(0, 0, .5, alpha=.5))
  box()
  
  screen(2)
  par(mar=c(0.5, 2, 0.25, .5))	
  plot( cwt.obj$x^2, bty='l', type='l', ylab='', xaxt='n')
  
  close.screen(all = TRUE)
}