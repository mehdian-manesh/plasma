rm(list=ls())
#---------- dimensions ---------
x0=1        #[cm]
E0=100      #[V cm-1]
mu0=3e5     #[cm2 v-1 s-1]
n0=5.53e7   #[cm-3] @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
t0=3.33e-8  #[s]
phi0=100    #[V]
D0=3e7      #[cm2 s-1]
u0=100      #[ev]
A0=3e7      #[cm s-1]
k_1_ion0=0.542   #[cm3 s-1]
Gama0=1.659e15   #[cm-2 s-1] @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
GamaE0=0.02655     #[ev cm-2 s-1]

#--------- constants, Bundary and initial values --------
V_applide=300/phi0
length_sys=1/x0
#pressure=1     #[torr]
n_gas=3.22e16/n0 #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
D_ele=3e5/D0
#mu_ele   ## alwase must be factored
gama=0.06 ## secondary electron yield

U_ion_ev=15.8/u0

# tf=1e-7/t0 ## time duration of simulation @@@@@@@@@@@@@@@@@@@@@@@@@
# N_t=10000  ## number of time steps
# dt=tf/N_t        #(tf/N_t)
# 
# #d=x0     ## Length of system
# N_x=100    ## Number of spatial steps
# dx=length_sys/N_x        #(d/N_x)
xf=1/x0
tf=1e-7/t0#1e-4/t0
dt=1e-11/t0 #اگر لازم شد اوردش رو هنوز هم کمتر کن
dx=0.001/x0
N_t=ceiling(tf/dt)
N_x=ceiling(xf/dx)

n_ele=matrix(0,nrow =N_x ,ncol =N_t ) #each time step is stored in a culomn of this matrix
n_ion=matrix(0,nrow = N_x,ncol = N_t)
phi=matrix(0,nrow = N_x,ncol = N_t+1)
u_ele=matrix(0,nrow = N_x,ncol = N_t)

#n_ele[1,]=0 ????????????????????????????
n_ele[N_x,]=0
n_ele[,1]=1  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

n_ion[,1]=1  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

phi[1,]=0
phi[N_x,]=V_applide

u_ele[,1]=1/u0
u_ele[1,]=5/u0
u_ele[N_x,]=0

#------------- functions ----------------
#mu_ionf=function(ix,it) 2.66e4/n_ion[ix,ifelse(it==1,1,it-1)]/(1 + (7.721e-3*abs(Ef(ix,ifelse(it==1,1,it-1)))/n_ion[ix,ifelse(it==1,1,it-1)])^1.5)^0.33   #(mu/mu_ele)
#D_ionf=function(ix,it) mu_ionf(ix,it)*8.33e-10  ##  T_ion is assumed 0.025 in [ev]  #(D/D0)
mu_ionf=function(ix,it) 4.411e19/mu0/n_gas/n0/(1+(7.721e+14*abs(Ef(ix,ifelse(it==1,1,it-1)))/n_gas*E0/n0)^1.5)^0.33 #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
D_ionf=function(ix,it) mu_ionf(ix,it)*mu0*0.025/D0 #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#M1_elef=function(ix,it) phi[ix,it+1]-phi[ifelse(ix==N_x,N_x,ix+1),it+1]
#M2_elef=function(ix,it) phi[ifelse(ix==1,1,ix-1),it+1]-phi[ix,it+1]
#M1_ionf=function(ix,it) 40*(phi[ix,it+1]-phi[ifelse(ix==N_x,N_x,ix+1),it+1])
#M2_ionf=function(ix,it) 40*(phi[ifelse(ix==1,1,ix-1),it+1]-phi[ix,it+1])

M1_elef=function(ix,it) (phi[ix+1,it+1]-phi[ix,it+1])*1/D_ele
M2_elef=function(ix,it) (phi[ix,it+1]-phi[ix-1,it+1])*1/D_ele
M1_ionf=function(ix,it) -(phi[ix+1,it+1]-phi[ix,it+1])*(mu_ionf(ix+1,it)+mu_ionf(ix,it))/(D_ionf(ix+1,it)+D_ionf(ix,it))
M2_ionf=function(ix,it) -(phi[ix,it+1]-phi[ix-1,it+1])*(mu_ionf(ix+1,it)+mu_ionf(ix,it))/(D_ionf(ix+1,it)+D_ionf(ix,it))

A1f=function(ix,it) -M1_elef(ix,it)*D_ele/dx/(exp(M1_elef(ix,it))-1)
A2f=function(ix,it) M1_elef(ix,it)*D_ele/dx/(exp(M1_elef(ix,it))-1)*exp(M1_elef(ix,it)) + M2_elef(ix,it)*D_ele/dx/(exp(M2_elef(ix,it))-1)
A3f=function(ix,it) -M2_elef(ix,it)*D_ele/dx/(exp(M2_elef(ix,it))-1)*exp(M2_elef(ix,it))

B1f=function(ix,it) -M1_ionf(ix,it)*D_ionf(ix+1,it+1)/dx/(exp(M1_ionf(ix,it))-1)
B2f=function(ix,it) M1_ionf(ix,it)*D_ionf(ix,it+1)/dx/(exp(M1_ionf(ix,it))-1)*exp(M1_ionf(ix,it)) + M2_ionf(ix,it)*D_ionf(ix,it+1)/dx/(exp(M2_ionf(ix,it))-1)
B3f=function(ix,it) -M2_ionf(ix,it)*D_ionf(ix-1,it+1)/dx/(exp(M2_ionf(ix,it))-1)*exp(M2_ionf(ix,it))

k_1_ionf=function(ix,it) 2.34e-8*(2/3*u_ele[ix,it]*u0)^0.59*exp(-17.44*3/2/u_ele[ix,it]/u0)/k_1_ion0
#k_1_ionf=function(ix,it) 2.34e-8/k_1_ion0*(u_ele[ix,it]*u0)^0.59*exp(-17.44/u_ele[ix,it]/u0)

Ef=function(ix,it) if(ix<N_x) return((phi[ix,it]-phi[ix+1,it])/dx) else return((phi[ix-1,it]-phi[ix,it])/dx)
#Gama_elef=function(ix,it) -n_ele[ix,it]*Ef(ix,it) - D_ele*(n_ele[ix,it]-n_ele[ifelse(ix==1,1,ix-1),it])/dx
Gama_elef=function(ix,it) {
  iix=floor(ix)
  return((n_ele[iix,it]*D_ele*exp(M1_elef(iix,it)) - n_ele[iix+1,it]*D_ele)/dx*M1_elef(iix,it)/(exp(M1_elef(iix,it))-1))
}
thomasf=function(zryb) { # zryb is a n*4 matrix that 4th column is for RHS of equations
  nrz=nrow(zryb)
  for(i in 2:nrz) {
    zryb[i,2]=zryb[i,2]-(zryb[i,1]/zryb[i-1,2])*zryb[i-1,3]
    zryb[i,4]=zryb[i,4]-(zryb[i,1]/zryb[i-1,2])*zryb[i-1,4]
  }
  x=matrix(0,nrow = nrz)
  x[nrz]=zryb[nrz,4]/zryb[nrz,2]
  for(i in (nrz-1):1){
    x[i]=(zryb[i,4]-zryb[i,3]*x[i+1])/zryb[i,2]
  }
  return(x)
}

#-------------------- main -----------------------

# First, we have to calculate phi[,1] and phi[,2] (phi in time step 1 & 2)
# As we can't calculate phi[,1] (becuse the n_ele[,0] is not exist), we assume that phi[,1]=phi[,2]
zrphi=matrix(1,nrow = N_x-2,ncol = 4)
zrphi[,2]=-2
zrphi[1,4]=-dx^2 * (n_ion[2,1]-n_ele[2,1])-phi[1,1]
for(i in 3:(N_x-2)){
  zrphi[i-1,4]=-dx^2 * (n_ion[i,1]-n_ele[i,1])
}
zrphi[N_x-2,4]=-dx^2 * (n_ion[N_x-1,1]-n_ele[N_x-1,1])-phi[N_x,1]
phi[2:(N_x-1),2]=thomasf(zrphi)
phi[,1]=phi[,2]
rm(zrphi)
cat("00.00 %")
for(i_t in 1:(N_t-1)){
  #calculate n_ele[,i_t+1] _____________________________________
  zrnel=matrix(0,nrow = N_x-2,ncol = 4)
  i_x=2
  zrnel[i_x-1,2]= dx + dt*A2f(i_x,i_t)
  zrnel[i_x-1,3]= dt*A1f(i_x,i_t)
  zrnel[i_x-1,4]= dx*n_ele[i_x,i_t]*(1 + dt*k_1_ionf(i_x,i_t)*n_gas) -dt*A3f(i_x,i_t)*n_ele[i_x-1,i_t+1]
  for(i_x in 3:(N_x-2)){
    zrnel[i_x-1,1]= dt*A3f(i_x,i_t)
    zrnel[i_x-1,2]= dx + dt*A2f(i_x,i_t)
    zrnel[i_x-1,3]= dt*A1f(i_x,i_t)
    zrnel[i_x-1,4]= dx*n_ele[i_x,i_t]*(1 + dt*k_1_ionf(i_x,i_t)*n_gas)
  }
  i_x=N_x-1
  zrnel[i_x-1,1]= dt*A3f(i_x,i_t)
  zrnel[i_x-1,2]= dx + dt*A2f(i_x,i_t)
  zrnel[i_x-1,4]= dx*n_ele[i_x,i_t]*(1 + dt*k_1_ionf(i_x,i_t)*n_gas) - dt*A1f(i_x,i_t)*n_ele[i_x+1,i_t+1]
  n_ele[2:(N_x-1),i_t+1]=thomasf(zrnel)
  rm(zrnel)
  
  #calculate n_ion[,i_t+1] _______________________________________
  #Neumann bundary condition ?????????????????????????????
  zrnio=matrix(0,nrow = N_x-2,ncol = 4)
  i_x=2
  zrnio[i_x-1,2]= dx + dt*B2f(i_x,i_t)+dt*B3f(i_x,i_t) #@@@@@@@@@@@@@@@@@@@@@@@@@@@
  zrnio[i_x-1,3]= dt*B1f(i_x,i_t)
  zrnio[i_x-1,4]= dx*n_ion[i_x,i_t]*(1 + dt*k_1_ionf(i_x,i_t)*n_gas) #@@@@@@@@@@@@@@@@@@@@@@@@@@@
  for(i_x in 3:(N_x-2)){
    zrnio[i_x-1,1]= dt*B3f(i_x,i_t)
    zrnio[i_x-1,2]= dx + dt*B2f(i_x,i_t)
    zrnio[i_x-1,3]= dt*B1f(i_x,i_t)
    zrnio[i_x-1,4]= dx*n_ion[i_x,i_t]*(1 + dt*k_1_ionf(i_x,i_t)*n_gas)
  }
  i_x=N_x-1
  zrnio[i_x-1,1]= dt*B3f(i_x,i_t)
  zrnio[i_x-1,2]= dx + dt*B2f(i_x,i_t)
  zrnio[i_x-1,4]= dx*n_ion[i_x,i_t]*(1 + dt*k_1_ionf(i_x,i_t)*n_gas) + dt*B1f(i_x,i_t) #@@@@@@@@@@@@@@@@@
  n_ion[2:(N_x-1),i_t+1]=thomasf(zrnio)
  rm(zrnio)
  n_ion[1,i_t+1]=n_ion[2,i_t+1] #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  n_ion[N_x,i_t+1]=n_ion[N_x-1,i_t+1] #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  #calculate u_ele[,i_t+1] _________________________________________
  zruel=matrix(0,nrow = N_x-2,ncol = 4)
  i_x=2
  zruel[i_x-1,2]= (dx + 1.667*dt*A2f(i_x,i_t))*n_ele[i_x,i_t+1]   # 1.667=5/3
  zruel[i_x-1,3]= 1.667*dt*A1f(i_x,i_t)*n_ele[i_x+1,i_t+1]
  zruel[i_x-1,4]= dx*u_ele[i_x,i_t]*n_ele[i_x,i_t] - dx*dt*Gama_elef(i_x,i_t)*Ef(i_x,i_t) - dx*dt*U_ion_ev* k_1_ionf(i_x,i_t)*n_ele[i_x,i_t]*n_gas - 1.667*dt*A3f(i_x,i_t)*n_ele[i_x-1,i_t+1]*u_ele[i_x-1,i_t+1]
  for(i_x in 3:(N_x-2)){
    zruel[i_x-1,1]= 1.667*dt*A3f(i_x,i_t)*n_ele[i_x-1,i_t+1]
    zruel[i_x-1,2]= (dx + 1.667*dt*A2f(i_x,i_t))*n_ele[i_x,i_t+1]
    zruel[i_x-1,3]= 1.667*dt*A1f(i_x,i_t)*n_ele[i_x+1,i_t+1]
    zruel[i_x-1,4]= dx*u_ele[i_x,i_t]*n_ele[i_x,i_t] - dx*dt*Gama_elef(i_x,i_t)*Ef(i_x,i_t) - dx*dt*U_ion_ev* k_1_ionf(i_x,i_t)*n_ele[i_x,i_t]*n_gas
  }
  i_x=N_x-1
  zruel[i_x-1,1]= 1.667*dt*A3f(i_x,i_t)*n_ele[i_x-1,i_t+1]
  zruel[i_x-1,2]= (dx + 1.667*dt*A2f(i_x,i_t))*n_ele[i_x,i_t+1]
  zruel[i_x-1,4]= dx*u_ele[i_x,i_t]*n_ele[i_x,i_t] - dx*dt*Gama_elef(i_x,i_t)*Ef(i_x,i_t) - dx*dt*U_ion_ev* k_1_ionf(i_x,i_t)*n_ele[i_x,i_t]*n_gas - 1.667*dt*A1f(i_x,i_t)*n_ele[i_x+1,i_t+1]*u_ele[i_x+1,i_t+1]
  u_ele[2:(N_x-1),i_t+1]=thomasf(zruel)
  rm(zruel)
  
  #calculate phi[,i_t+2] _______________________________________________
  zrphi=matrix(1,nrow = N_x-2,ncol = 4)
  zrphi[,2]=-2
  zrphi[1,4]=-dx^2*(n_ion[2,i_t+1]-n_ele[2,i_t+1]) - phi[1,i_t+2]
  for(i_x in 3:(N_x-2)){
    zrphi[i_x-1,4]=-dx^2*(n_ion[i_x,i_t+1]-n_ele[i_x,i_t+1])
  }
  zrphi[N_x-2,4]=-dx^2*(n_ion[N_x-1,i_t+1]-n_ele[N_x-1,i_t+1]) - phi[N_x,i_t+2]
  phi[2:(N_x-1),i_t+2]=thomasf(zrphi)
  rm(zrphi)
  #______________________________________________________________
  
  #n_ele[1,ifelse(i_t+2>N_t,i_t+1,i_t+2)]=n_ion[1,i_t]*gama*mu_ionf(1,i_t)
  cat(paste("\b\b\b\b\b\b\b\b",i_t/N_t*100,"%"))
  stopifnot(!is.nan(n_ion[2,i_t]))
}
####proc.time()[3]