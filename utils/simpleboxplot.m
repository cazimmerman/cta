function simpleboxplot(x,y,c)

yprime = prctile(y,[10 25 50 75 90]);

plot([x x],[yprime(1) yprime(5)],'color',c,'linewidth',0.5)
plot([x x],[yprime(2) yprime(4)],'color',c,'linewidth',4)
plot(x+[-.25 +.25],[yprime(3) yprime(3)],'color',c,'linewidth',0.5)