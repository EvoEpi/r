#testing whether the distribution of x and y are identical

x<-sleep$extra[sleep$group==1]
y<-sleep$extra[sleep$group==2]

pooled<-c(x,y) #pooled sample, assumes equal variance
boot.t<-0
for (i in 1:1000) {
  sample.index<-sample(c(1:length(pooled)),replace=TRUE)
  sample.x<-pooled[sample.index][1:length(x)]
  sample.y<-pooled[sample.index][-c(1:length(y))]
  boot.t[i]<-t.test(sample.x,sample.y)$statistic
}
p.pooled<-(1 + sum(abs(boot.t) > abs(t.test(x,y)$statistic))) / (1000+1)
p.pooled
