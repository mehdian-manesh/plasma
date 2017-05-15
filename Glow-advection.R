rm(list = ls())
starttime=proc.time()[3]
#===================================== initial values =======================================
##-------------- Dimentions ---------------
x0=1      #[cm]
E0=100    #[V cm-1]
mu0=3e5   #[cm2 V-1 s-1]
n0=5.53e7 #[cm-3]
t0=3.33e-8#[s]
phi0=100  #[V]
D0=3e7    #[cm2 s-1]
u0=100    #[ev]
K_1_io_0=0.543  #[cm3 s-1]
Gama_0=1.659e15 #[cm-2 s-1]
Gama_E_0=0.0265 #[ev cm-2 s-1]
P=1#[torr]
##--------------- Dimentionless values ----------------
V_applied=300/phi0
xf=1/x0
u_ion=15.8/u0
gama=0.06
n_gas=3.22e16/n0
#_______changables_________
tf=1e-7/t0#1e-4/t0
dt=1e-10/t0 #اگر لازم شد اوردش رو هنوز هم کمتر کن
dx=0.001/x0
N_t=ceiling(tf/dt)
N_x=ceiling(xf/dx)
##-------------- Defenition & boundary condition ----------------
n_ele=matrix(nrow = N_x,ncol = N_t)
u_ele=matrix(nrow = N_x,ncol = N_t)
n_ion=matrix(nrow = N_x,ncol = N_t)
phi=matrix(0,nrow = N_x,ncol = N_t)

phi[1,]=0   #i_x=1 is the cathod
phi[N_x,]=V_applied #i_x=N_x is the anode

n_ele[N_x,]=0
n_ele[,1]=1e3/n0

u_ele[1,]=5/u0
u_ele[N_x,]=0
u_ele[2:(N_x-1),1]=1/u0

n_ion[,1]=1e3/n0

mu_ele=3e5/P/mu0
D_ele=3e5/P/D0
#=================================== functions =============================================
mu_ionf=function(ix,it) 4.411e19/mu0/n_gas/n0/(1+(7.721e+14*abs(Ef(ix,it))/n_gas*E0/n0)^1.5)^0.33
D_ionf=function(ix,it) mu_ionf(ix,it)*mu0*0.025/D0

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
M_elef=function(ix,it) (phi[ix+1,it]-phi[ix,it])*(mu_ele/D_ele)
M_ionf=function(ix,it) -(phi[ix+1,it]-phi[ix,it])*(mu_ionf(ix+1,it)+mu_ionf(ix,it))/(D_ionf(ix+1,it)+D_ionf(ix,it))

Gama_elef=function(ix,it) {
  iix=floor(ix)
  return((n_ele[iix,it]*D_ele*exp(M_elef(iix,it)) - n_ele[iix+1,it]*D_ele)/dx*M_elef(iix,it)/(exp(M_elef(iix,it))-1))
}
Gama_ionf=function(ix,it) {
  iix=floor(ix)
  return((n_ion[iix,it]*D_ionf(iix,it)*exp(M_ionf(iix,it)) - n_ion[iix+1,it]*D_ionf(iix+1,it))/dx*M_ionf(iix,it)/(exp(M_ionf(iix,it))-1))
}
Gama_E_elef=function(ix,it) {
  iix=floor(ix)
  return((u_ele[iix,it]*n_ele[iix,it]*D_ele*exp(M_elef(iix,it)) - u_ele[iix+1,it]*n_ele[iix+1,it]*D_ele)/dx*M_elef(iix,it)/(exp(M_elef(iix,it))-1))
}
Ef=function(ix,it) {
  if(ix<N_x) {
    return((phi[ix,it]-phi[ix+1,it])/dx) 
  } else {
    return((phi[ix-1,it]-phi[ix,it])/dx)
  }
}

k_1_ionf=function(ix,it) 2.34e-8*(2/3*u_ele[ix,it]*u0)^0.59*exp(-17.44*3/2/u_ele[ix,it]/u0)/K_1_io_0


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

#====================================== Main =================================================
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

phi[,1]=phi[,2] #???????????????????????????آیا این لازم است؟

cat("00.00 %")

for(i_t in 1:(N_t-1)){
  #------------------------------ cathod -------------------------------------
  i_x=1
  #_______n_ele_________
  n_ele[i_x,i_t+1]=n_ele[i_x,i_t] - dt/dx*(Gama_elef(i_x+1,i_t)-Gama_elef(i_x,i_t)) + k_1_ionf(i_x,i_t)*n_ele[i_x,i_t]*n_gas*dt
  #n_ion has neumann bondary condition.***
  #u_ele is fixed.
  #---------------------------------------------------------------------------
  for(i_x in 2:(N_x-1)){
    n_ele[i_x,i_t+1]=n_ele[i_x,i_t] - dt/dx*(Gama_elef(i_x+0.5,i_t)-Gama_elef(i_x-0.5,i_t)) + k_1_ionf(i_x,i_t)*n_ele[i_x,i_t]*n_gas*dt
    n_ion[i_x,i_t+1]=n_ion[i_x,i_t] - dt/dx*(Gama_ionf(i_x+0.5,i_t)-Gama_ionf(i_x-0.5,i_t)) + k_1_ionf(i_x,i_t)*n_ion[i_x,i_t]*n_gas*dt
    u_ele[i_x,i_t+1]=(u_ele[i_x,i_t]*n_ele[i_x,i_t] - 5/3*dt/dx*(Gama_E_elef(i_x+0.5,i_t)-Gama_E_elef(i_x-0.5,i_t)) - dt*Gama_elef(i_x,i_t)*Ef(i_x,i_t) - u_ion*k_1_ionf(i_x,i_t)*n_ele[i_x,i_t]*n_gas*dt )/n_ele[i_x,i_t+1]
  }
  #_______Neumann B.C. for n_ion_________
  n_ion[1,i_t+1]=n_ion[2,i_t+1]
  n_ion[N_x,i_t+1]=n_ion[N_x-1,i_t+1]
  #----------------------------- anode ----------------------------------------
  #n_ele is fixed
  #n_ele[N_x,i_t+1]=n_ele[N_x,i_t] - dt/dx*(Gama_elef(i_x,i_t)-Gama_elef(i_x-1,i_t)) + k_1_ionf(i_x,i_t)*n_ele[i_x,i_t]*n_gas*dt
  #n_ion has neumann bondary condition.***
  #u_ele is fixed.  
  #------------------------ calculate phi[,i_t+2] -----------------------------
  zrphi[,1:3]=c(1,-2,1)
  zrphi[1,4]=-dx^2*(n_ion[2,i_t+1]-n_ele[2,i_t+1]) - phi[1,i_t+2]
  for(i_x in 3:(N_x-2)){
    zrphi[i_x-1,4]=-dx^2*(n_ion[i_x,i_t+1]-n_ele[i_x,i_t+1])
  }
  zrphi[N_x-2,4]=-dx^2*(n_ion[N_x-1,i_t+1]-n_ele[N_x-1,i_t+1]) - phi[N_x,i_t+2]
  phi[2:(N_x-1),i_t+2]=thomasf(zrphi)
  #----------------------------------------------------------------------------
  
  #n_ele[1,ifelse(i_t+2>N_t,i_t+1,i_t+2)]=n_ion[1,i_t]*gama*mu_ionf(1,i_t)
  cat(paste("\b\b\b\b\b\b\b\b",i_t/N_t*100,"%"))
  stopifnot(!is.nan(n_ion[2,i_t]))
}
print(proc.time()[3]-starttime)