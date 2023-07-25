
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
bxx=order(x,decreasing=FALSE)
byy=order(y,decreasing=FALSE)
i=1
while ( x[bxx[i]]==x[bxx[1]]) { i=i+1;}
ffx=x[bxx[i]]
print (i)
bl=min(ffx)+log(2/3)
i=1
while ( y[byy[i]]==y[byy[1]]) { i=i+1;}
ffy=y[byy[i]]
print (i)


bly=min(ffy)+log(2/3)
#xt=x
# yt=y
xmax=6
xmin=-6
ymax=6
ymin=-6
plot(x,y,col=df$V3,pch=3,cex=3,xlim=c(xmin,xmax),ylim=c(ymin,ymax))
#plot(x,y,col=df$V3,pch=3,cex=3,xlim=c(-4,-1.5),ylim=c(-4,-1.5))
lines(c(-20,20),c(-20,20))
lines(c(bl,bl),c(-20,20),col="red")
lines(c(-20,20),c(bly,bly),col="green")
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

