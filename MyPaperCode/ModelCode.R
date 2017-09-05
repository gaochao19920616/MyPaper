N=1000
csmd_wl=0.02
csmd_xx=0.51
stp=20

lamda=0.3
delta1=0.6
delta2=0.3
beta_U=0.25
gama=0.2
miu=0.7
alfa=-0.5


library(igraph)

h=barabasi.game(N,m=5,directed = FALSE)   #BA网络
h=erdos.renyi.game(N,p=0.01)      #ER网络

hh=erdos.renyi.game(N,p=0.002)     #随机加2000条边


B=get.adjacency(h,sparse=FALSE)
C=get.adjacency(hh,sparse=FALSE)
A=B+C                                   #信息层比物理层多2000边

k=1:N
for(i in 1:N)
{
  k[i]=sum(B[,i])
}                            #物理层节点度







# 循环10次

BETA_A_NEW=matrix(0,nr=10,nc=50)
BETA_I_NEW=matrix(0,nr=10,nc=50)


for(rep in 1:10) {    


#  MC调beta的循环

termi=50
BETA_A=1:termi 
BETA_I=1:termi

for(l in 1:termi) {
  
  beta_U=l/termi
  
  
  
  
  #MC
  
  x=runif(N,min=0,max=1)
  for (i in 1:N)
  {if (x[i]>csmd_wl)
  {x[i]=0}
    else
    {x[i]=1}
  }                    #设置物理层初始状态
  
  m=runif(N,min=0,max=1)
  for (i in 1:N)
  {if (m[i]>csmd_xx)
  {m[i]=0}
    else
    {m[i]=1}
  }                    #设置信息层初始状态
  
 
  
  
  
  
  
  #循环
  
  for(t in 1:stp){
    
    
    n=1:N ;n[]=0
    y=1:N ;y[]=0       # 设置N维空向量
    
    
    
    #第一步，更新信息层
    
    for(i in 1:N)
    {if(m[i]==1)   #信息层为知道时
    {if(x[i]==1)  # 节点被感染
    {if(runif(1,min=0,max=1)>delta2)  #不忘记
    {n[i]=1}
      else                             #忘记
      {n[i]=0}
    }
      else             #节点未感染
      {if(runif(1,min=0,max=1)>delta1)
      {n[i]=1}
        else                              #不忘记
        {n[i]=0}
      }                                 #忘记
    }
      else           #信息层不知道
      {n[i]=0
      for (j in 1:N)
      {       
        if ((A[j,i]==1)&(m[j]==1)&(runif(1,min=0,max=1)<lamda))
        {n[i]=1}           #信息层变为知道
      }
      }
      
      
      
      
      
      #第二步，更新物理层
      
      if(x[i]==1)
      { if(runif(1,min=0,max=1)<miu)
      {y[i]=0}
        else
        {y[i]=1}
      }  
      else
      {if(n[i]==0)
      {y[i]=0
      
      for(j in 1:N)
      {
        if((B[j,i]==1)&(x[j]==1)&(runif(1,min=0,max=1)<beta_U))
        {y[i]=1}
      }      
      }
        else
        {
          y[i]=0
          for(j in 1:N)
          {
            if ((B[j,i]==1)&(x[j]==1)&(runif(1,min=0,max=1)<(gama*beta_U)))
            {y[i]=1}
          }
        }
      }
      
      
      
      
      
      
      
      #第三步，上传信息  
      
      
      
      if((y[i]==1)&(n[i]==0))
      { if(runif(1,min=0,max=1)<k[i]^alfa)
      {n[i]=0}
        
        else
        {n[i]=1}
      }
      
    }
    
    
    m=n            #将n更新给m      
    x=y            #将y更新给x
    
  }    
  
  
  BETA_A[l]=sum(m)/N 
  BETA_I[l]=sum(x)/N 
}  

BETA_A_NEW[rep,]=BETA_A
BETA_I_NEW[rep,]=BETA_I

}


for (av in 1:50){
BETA_A_AVER[av]=sum(BETA_A_NEW[,av])/10
BETA_I_AVER[av]=sum(BETA_I_NEW[,av])/10
}


#   MC     ROU VS BETA作图
xzhou=(1:termi)/termi

plot(xzhou,BETA_A_AVER,ylim=c(0,1)) 
plot(xzhou,BETA_I_AVER,ylim=c(0,1)) 

points(xzhou,BETA_A_AVER)
points(xzhou,BETA_I_AVER) 








#MMCA调beta

ter=50
BETA1=1:ter 
BETA2=1:ter

for(l in 1:ter) {
  
  beta_U=l/ter
  
  
  # MMCA
  
  
  MMCA=20
  
  PAI=1:N;PAI[]=0.01
  PAS=1:N;PAS[]=0.5
  PUI=1:N;PUI[]=0.01
  PUS=1:N;PUS[]=0.48
  
  PAI_UPDATE=1:N;PAI_UPDATE[]=0
  PAS_UPDATE=1:N;PAS_UPDATE[]=0
  PUI_UPDATE=1:N;PUI_UPDATE[]=0
  PUS_UPDATE=1:N;PUS_UPDATE[]=0
  
  r=1:N ;r[]=0
  qu=1:N ;qu[]=0
  qa=1:N ;qa[]=0
  
  
  
  
  R=matrix(0,nr=N,nc=N)
  QU=matrix(0,nr=N,nc=N)
  QA=matrix(0,nr=N,nc=N)
  
  
  for(t in 1:MMCA){
    
    for (i in 1:N){
      
      
      
      for(j in 1:N){
        R[j,i]=1-A[j,i]*(PAI[j]+PAS[j])*lamda
        QU[j,i]=1-B[j,i]*(PAI[j]+PUI[j])*beta_U
        QA[j,i]=1-B[j,i]*(PAI[j]+PUI[j])*gama*beta_U
        
      }  
      r[i]=prod(R[,i])
      qu[i]=prod(QU[,i])
      qa[i]=prod(QA[,i])
      
      
      
      
      PUS_UPDATE[i]=PUI[i]*r[i]*miu+PAI[i]*delta2*miu+PUS[i]*r[i]*qu[i]+PAS[i]*qu[i]*delta1 
      
      PUI_UPDATE[i]=PUI[i]*r[i]*(1-miu)*(k[i]^alfa)+PAI[i]*delta2*(1-miu)*(k[i]^alfa)+PUS[i]*r[i]*(1-qu[i])*(k[i]^alfa)+PAS[i]*delta1*(1-qu[i])*(k[i]^alfa)        
      
      PAS_UPDATE[i]=PUI[i]*(1-r[i])*miu+PAI[i]*(1-delta2)*miu+PUS[i]*(1-r[i])*qa[i]+PAS[i]*(1-delta1)*qa[i]
      
      PAI_UPDATE[i]=PUI[i]*(r[i]*(1-miu)*(1-(k[i]^alfa))+(1-r[i])*(1-miu))+PAI[i]*(delta2*(1-miu)*(1-(k[i]^alfa))+(1-delta2)*(1-miu))+PUS[i]*(r[i]*(1-qu[i])*(1-(k[i]^alfa))+(1-r[i])*(1-qa[i]))+PAS[i]*(delta1*(1-qu[i])*(1-(k[i]^alfa))+(1-delta1)*(1-qa[i]))                         
    }
    
    
    PUI=PUI_UPDATE
    PAI=PAI_UPDATE
    PUS=PUS_UPDATE
    PAS=PAS_UPDATE
  }
  
  PA=PAS+PAI
  
  PI=PAI+PUI
  
  
  
  
  BETA1[l]=sum(PA)/N
  BETA2[l]=sum(PI)/N
  
  
  
}


#roua vs beta ,rou vs beta 图
Xzhou=(1:ter)/ter

plot(Xzhou,BETA1,ylim=c(0,1),xlab="beta_U",ylab="rou_I",pch=20,col="blue",bg="blue",type="l",lwd="4")
plot(Xzhou,BETA2,ylim=c(0,1),xlab="beta_U",ylab="rou_I",pch=20,col="blue",bg="blue",type="l",lwd="4")


#MMCA算阈值

H=matrix(0,nr=1000,nc=1000)
for(i in 1:N){
  for(j in 1:N){
    H[i,j]=(1-(1-gama)*PA[i])*B[i,j]
  }
}

beta_MMCA=miu/max(eigen(H)$values)





#MMCA调 lamda

ter=50
Lam=(1:ter)/50 
beta_C=1:ter

for(l in 1:ter) {
  
  lamda=l/ter
  
  
  # MMCA
  
  
  MMCA=20
  
  PAI=1:N;PAI[]=0.01
  PAS=1:N;PAS[]=0.5
  PUI=1:N;PUI[]=0.01
  PUS=1:N;PUS[]=0.48
  
  PAI_UPDATE=1:N;PAI_UPDATE[]=0
  PAS_UPDATE=1:N;PAS_UPDATE[]=0
  PUI_UPDATE=1:N;PUI_UPDATE[]=0
  PUS_UPDATE=1:N;PUS_UPDATE[]=0
  
  r=1:N ;r[]=0
  qu=1:N ;qu[]=0
  qa=1:N ;qa[]=0
  
  
  
  
  R=matrix(0,nr=N,nc=N)
  QU=matrix(0,nr=N,nc=N)
  QA=matrix(0,nr=N,nc=N)
  
  
  for(t in 1:MMCA){
    
    for (i in 1:N){
      
      
      
      for(j in 1:N){
        R[j,i]=1-A[j,i]*(PAI[j]+PAS[j])*lamda
        QU[j,i]=1-B[j,i]*(PAI[j]+PUI[j])*beta_U
        QA[j,i]=1-B[j,i]*(PAI[j]+PUI[j])*gama*beta_U
        
      }  
      r[i]=prod(R[,i])
      qu[i]=prod(QU[,i])
      qa[i]=prod(QA[,i])
      
      
      
      
      PUS_UPDATE[i]=PUI[i]*r[i]*miu+PAI[i]*delta2*miu+PUS[i]*r[i]*qu[i]+PAS[i]*qu[i]*delta1 
      
      PUI_UPDATE[i]=PUI[i]*r[i]*(1-miu)*(k[i]^alfa)+PAI[i]*delta2*(1-miu)*(k[i]^alfa)+PUS[i]*r[i]*(1-qu[i])*(k[i]^alfa)+PAS[i]*delta1*(1-qu[i])*(k[i]^alfa)        
      
      PAS_UPDATE[i]=PUI[i]*(1-r[i])*miu+PAI[i]*(1-delta2)*miu+PUS[i]*(1-r[i])*qa[i]+PAS[i]*(1-delta1)*qa[i]
      
      PAI_UPDATE[i]=PUI[i]*(r[i]*(1-miu)*(1-(k[i]^alfa))+(1-r[i])*(1-miu))+PAI[i]*(delta2*(1-miu)*(1-(k[i]^alfa))+(1-delta2)*(1-miu))+PUS[i]*(r[i]*(1-qu[i])*(1-(k[i]^alfa))+(1-r[i])*(1-qa[i]))+PAS[i]*(delta1*(1-qu[i])*(1-(k[i]^alfa))+(1-delta1)*(1-qa[i]))                         
    }
    
    
    PUI=PUI_UPDATE
    PAI=PAI_UPDATE
    PUS=PUS_UPDATE
    PAS=PAS_UPDATE
  }
  
  PA=PAS+PAI
  
  PI=PAI+PUI
  
  
  H=matrix(0,nr=1000,nc=1000)
  for(i in 1:N){
    for(j in 1:N){
      H[i,j]=(1-(1-gama)*PA[i])*B[i,j]
    }
  }
  
  beta_C[l]=miu/max(eigen(H)$values)  
  
  
  
  
}









# 调lamda  dt=0.2 mu=0.6



ter=50
Lam=((1:ter)-1)/50 
beta_C=1:ter

for(l in 1:ter) {
  
  lamda=(l-1)/ter
  
  
  #迭代法解方程组
  DD=20
  PA=1:N;PA[]=0.02
  
  PA_UPDATE=1:N;PA_UPDATE[]=0
  
  
  r=1:N ;r[]=0
  R=matrix(0,nr=N,nc=N)
  
  for(t in 1:DD){
    
    for (i in 1:N){
      
      
      
      for(j in 1:N){
        R[j,i]=1-A[j,i]*PA[j]*lamda
        
        
      }  
      r[i]=prod(R[,i])
      
      PA_UPDATE[i]=(1-PA[i])*(1-r[i])+PA[i]*(1-delta1) 
      
      
    }
    
    
    PA=PA_UPDATE
    
  }
  H=matrix(0,nr=1000,nc=1000)
  for(i in 1:N){
    for(j in 1:N){
      H[i,j]=(1-(1-gama)*PA[i])*B[i,j]
    }
  }
  
  beta_C[l]=miu/max(eigen(H)$values)  
  
  
  
  
}








