#»­Í¼±à¼­¿ò

library(igraph)
termi=50
ter=50

xzhou=(1:termi)/termi

par(mar = c(5,7,1,1)) # Í¼µÄÉÏÏÂ×óÓÒ¾àÀë


#ER or SFÍ¼   MC¶Ô±ÈMMCA 
#RouA
plot(Xzhou,BETA1,ylim=c(0.5,0.9),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta),ylab=expression(rho^A),tck=0.02,pch=20,col="steelblue3",type="l",lwd=4,lty=1,)
points(xzhou,BETA_A_AVER,pch=1,col="red",bg="yellow",cex=1.3,lty=1)
legend(0, 0.9,c("MMCA","MC"),pch=c(NA,1),col=c("steelblue3","red"),lty=c(1,0),lwd=c(4,NA),cex=1.5)
#RouI
plot(Xzhou,BETA2,ylim=c(0,1),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta),ylab=expression(rho^I),tck=0.02,pch=20,col="steelblue3",type="l",lwd=4,lty=1,)
points(xzhou,BETA_I_AVER,pch=1,col="red",bg="yellow",cex=1.3,lty=1)
legend(0, 1,c("MMCA","MC"),pch=c(NA,1),col=c("steelblue3","red"),lty=c(1,0),lwd=c(4,NA),cex=1.5)

#ER vs SFÍ¼ MC,RouI
plot(Xzhou,BETA_I_AVER,ylim=c(0,1),las=1,cex.lab=1.4,xlab = expression(beta),ylab=expression(rho^I),tck=0.02,pch=4,col="violet",cex=1.3,type="p",lwd=2,lty=0)
points(xzhou,BETA_I_AVER,pch=6,col="skyblue3",cex=1.3,lwd=2,lty=0)
legend(0, 1,c("ER network","SF network"),pch=c(4,6),col=c("violet","skyblue3"),lty=c(0,0),lwd=c(2,2))

#gamma2 = 0 ,0.5 ,1
#RouA

plot(Xzhou,BETA1,ylim=c(0.4,1),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta),ylab=expression(rho^A),tck=0.02,pch=20,col="gold",cex=1.3,type="p",lwd=2,lty=1)
points(xzhou,BETA1,pch=20,col="orange",cex=1.3,lwd=2,lty=1)
points(xzhou,BETA1,pch=20,col="red",cex=1.3,lwd=2,lty=1)
legend(0.75, 0.65,cex=1.5,c(expression(gamma[2]==0.0),expression(gamma[2]==0.5),expression(gamma[2]==1.0)),pch=c(20,20,20),col=c("gold","orange","red"),lty=c(0,0,0),lwd=c(2,2,2))

#RouI

plot(Xzhou,BETA2,ylim=c(0,1),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta),ylab=expression(rho^I),tck=0.02,pch=1,col="gold",cex=1.2,type="p",lwd=1,lty=1)
points(xzhou,BETA2,pch=4,col="orangered",cex=1.2,lwd=1,lty=1)
points(xzhou,BETA2,pch=20,col="red",cex=1.2,lwd=1,lty=1)
legend(0, 1,cex=1.5,c(expression(gamma[2]==0.0),expression(gamma[2]==0.5),expression(gamma[2]==1.0)),pch=c(1,4,20),col=c("gold","orangered","red"),lty=c(0,0,0),lwd=c(1,1,1))

#alfa= 0.1 0.5 0.9
#RouA
plot(Xzhou,BETA1,ylim=c(0,1),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta),ylab=expression(rho^A),tck=0.02,pch=20,col="gold",cex=1.3,type="p",lwd=2,lty=1)
points(xzhou,BETA1,pch=20,col="orange",cex=1.3,lwd=2,lty=1)
points(xzhou,BETA1,pch=20,col="red",cex=1.3,lwd=2,lty=1)
legend(0.75, 0.37,cex=1.5,c(expression(alpha==0.1),expression(alpha==0.5),expression(alpha==0.9)),pch=c(20,20,20),col=c("gold","orange","red"),lty=c(0,0,0),lwd=c(2,2,2), adj = c(0.2, 0.25))

#RouI

plot(Xzhou,BETA2,ylim=c(0,1),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta),ylab=expression(rho^I),tck=0.02,pch=1,col="gold",cex=1.2,type="p",lwd=1,lty=1)
points(xzhou,BETA2,pch=4,col="orangered",cex=1.2,lwd=1,lty=1)
points(xzhou,BETA2,pch=20,col="red",cex=1.2,lwd=1,lty=1)
legend(0, 1,cex=1.5,c(expression(alpha==0.1),expression(alpha==0.5),expression(alpha==0.9)),pch=c(1,4,20),col=c("gold","orangered","red"),lty=c(0,0,0),lwd=c(1,1,1),adj = c(0.2, 0.25))

#gamma1 = 0 ,0.2,0.6 ,1

#RouA
plot(Xzhou,BETA1,ylim=c(0.5,0.9),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta),ylab=expression(rho^A),tck=0.02,pch=1,col="gold",cex=1.3,type="p",lwd=2,lty=1)
points(xzhou,BETA1,pch=4,col="orange",cex=1.3,lwd=2,lty=1)
points(xzhou,BETA1,pch=20,col="red",cex=1.3,lwd=2,lty=1)
points(xzhou,BETA1,pch=2,col="maroon3",cex=1.3,lwd=2,lty=1)
legend(0.75, 0.7,cex=1.5,c(expression(gamma[1]==0),expression(gamma[1]==0.2),expression(gamma[1]==0.6),expression(gamma[1]==1)),pch=c(1,4,20,2),col=c("gold","orange","red","maroon3"),lty=c(0,0,0,0),lwd=c(2,2,2,2), adj = c(0.2, 0.25))

#RouI

plot(Xzhou,BETA2,ylim=c(0,1),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta),ylab=expression(rho^I),tck=0.02,pch=1,col="gold",cex=1.2,type="p",lwd=1,lty=1)
points(xzhou,BETA2,pch=4,col="orangered",cex=1.2,lwd=1,lty=1)
points(xzhou,BETA2,pch=20,col="red",cex=1.2,lwd=1,lty=1)
points(xzhou,BETA2,pch=2,col="maroon3",cex=1.2,lwd=1,lty=1)
legend(0, 1,cex=1.5,c(expression(gamma[1]==0),expression(gamma[1]==0.2),expression(gamma[1]==0.6),expression(gamma[1]==1)),pch=c(1,4,20,2),col=c("gold","orange","red","maroon3"),lty=c(0,0,0,0),lwd=c(2,2,2,2), adj = c(0.2, 0.25))

#µ÷lamda
plot(beta_C,Lam,xlim=c(0.02,0.075),ylim=c(0,1),las=1,cex.lab=1.4,xlab = expression(beta[c]),ylab=expression(lambda),tck=0.02,pch=20,col="gold",cex=1.2,type="p",lwd=1,lty=1)
points(beta_C,Lam,pch=20,col="orangered",cex=1.2,lwd=1,lty=1)
points(beta_C,Lam,pch=20,col="red",cex=1.2,lwd=1,lty=1)
points(beta_C,Lam,pch=20,col="maroon3",cex=1.2,lwd=1,lty=1)
legend(0.018, 1,c(expression(paste(delta[1]==0.8,",",mu==0.4)),expression(paste(delta[1]==0.8,",",mu==0.6)),expression(paste(delta[1]==0.2,",",mu==0.6))),pch=c(20,20,20,20),col=c("gold","orange","red","maroon3"),lty=c(0,0,0,0),lwd=c(2,2,2,2), adj = c(0.2, 0.25))

#D-BC
plot(beta_C,Lam,xlim=c(0.048,0.103),ylim=c(0,1),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta[c]),ylab=expression(delta[1]),tck=0.02,pch=20,col="red",cex=1.2,type="p",lwd=1,lty=1)
legend(0.087, 1,cex=1.5,expression(mu==0.6),col="red",pch=20,lty=0,lwd=2, adj = c(0.2, 0.25))
#M-BC
plot(beta_C,Lam,xlim=c(0,0.082),ylim=c(0,1),las=1,cex.lab=2.3,cex.axis=2,xlab = expression(beta[c]),ylab=expression(mu),tck=0.02,pch=20,col="red",cex=1.2,type="p",lwd=1,lty=1)
legend(0, 1,cex=1.5,expression(delta[1]==0.8),col="red",pch=20,lty=0,lwd=2, adj = c(0.2, 0.25))