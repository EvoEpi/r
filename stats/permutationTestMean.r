#test whether population means are equal, without making any assumptions about their variance
#sample from H0 separately, no assumption about equal variance

x<-sleep$extra[sleep$group==1]
y<-sleep$extra[sleep$group==2]
xt<-x - mean(x) + mean(sleep$extra)
yt<-y - mean(y) + mean(sleep$extra)

boot.t<-c(1:1000)
for (i in 1:1000){
  sample.x<-sample(xt,replace=TRUE)
  sample.y<-sample(yt,replace=TRUE)
  boot.t[i]<-t.test(sample.x,sample.y)$statistic
}
p.h0<-(1 + sum(abs(boot.t) > abs(t.test(x,y)$statistic))) / (1000+1)
p.h0
