biennial.func<-function(ds, period.range.short=c(0.8,1.2),period.range.long=c(1.8,2.2)){
  annual.indices<-which(ds$periods>=period.range.short[1] & ds$periods<=period.range.short[2])
  biennial.indices<-which(ds$periods>=period.range.long[1] & ds$periods<=period.range.long[2])
  annual.power<-mean(ds$global[annual.indices])
  biennial.power<-mean(ds$global[biennial.indices])
  biennial.annual.ratio<-biennial.power/annual.power
}