
x=(0:50)/50
lin=1-(1-x)
two=1-(1-x)^2
three=1-(1-x)^3
four=1-(1-x)^4
plot(x,lin,pch="1",xlim=c(0,.2))
lines(x,lin,pch="1")
points(x,two,col="blue",pch="2")
lines(x,two,col="blue",pch="2")
points(x,three,col="red",pch="3")
lines(x,three,col="red",pch="3")
points(x,four,col="green",pch="4")
lines(x,four,col="green",pch="4")
grid(col="black")


