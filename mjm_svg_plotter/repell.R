
# df<-read.table("test2.dat",header=F)
dfpos<-read.table("xm",header=F)
xt=dfpos$V1;yt=dfpos$V2
x=dfpos$V3;y=dfpos$V4
# also have dvx, dvy, v
sv=.1
# grad v, v, theta
dvx=dfpos$V7*sv
dvy=dfpos$V8*sv
# this needs to find the second smallest ... doh these were zero... 
#xt=x
# yt=y
xmax=-4
print(max(x,xt))
xmin=-8
ymax=-4
ymin=-8
plot(x,y,col=df$V3,pch=3,cex=3,xlim=c(xmin,xmax),ylim=c(ymin,ymax))
plot(x,y,col=df$V3,pch=3,cex=3)
#plot(x,y,col=df$V3,pch=3,cex=3,xlim=c(-4,-1.5),ylim=c(-4,-1.5))
lines(c(-20,20),c(-20,20))
#text(x,y,labels=paste(df$V9,"\n",df$V10))
#text(xt,yt,labels=paste(df$V10))
text(xt,yt,labels=paste(dfpos$V5))
#text(xt,yt,labels=(df$V10))
grid(col="black")
for(i in 1:length(x))
{
lines(c(xt[i],x[i]),c(yt[i],y[i]))
#lines(c(x[i]+dvx[i],x[i]),c(y[i]+dvy[i],y[i]))


}

