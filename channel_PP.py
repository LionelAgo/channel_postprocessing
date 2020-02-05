

import numpy as np
import matplotlib.pyplot as plt
import h5py
import Mesh_hex
import Plots
import importlib
importlib.reload(Mesh_hex)
importlib.reload(Plots)

from Mesh_hex import extract_grid
from Mesh_hex import Read_solution
from Mesh_hex import Read_stats
from Mesh_hex import MeanXZ
from Mesh_hex import MeanZ
from Mesh_hex import wallunit
from Plots import plot_2d

rhom = 0.014
tauw = 5.92464e-05
alp=1.025
alp=1
utau=alp*(tauw/rhom)**0.5
mu = 5.05967e-06
Retau=(rhom/mu)*utau
#T=20000
#Tp=(rhom/mu)*T*utau**2
#100/((rhom/mu)*utau**2)


nx=62
ny=19
nz=60
order=5
ff = '/home/cfd/Desktop/Kara/channel180/Channel-2U2U.pyfrm'
npart=12
Mesh,Lx,Ly,Lz =extract_grid(npart,ff,nx,ny,nz,order)

#%%
X,Y,Z=np.meshgrid(Lx,Ly,Lz, indexing='ij')

yp=(Ly[:int(1+np.floor(len(Ly)/2))]+1)*Retau

#%%                   Snapshots

timestep=2050
#U=Read_solution(2050)
U=Read_solution(2050,1)

#plt.pcolormesh(np.squeeze(X[:,0,:]), np.squeeze(Z[:,0,:]), np.squeeze(U[:,8,:]), cmap='jet')
plt.imshow(np.squeeze(U[:,8,:]).T,extent=[np.squeeze(X[0,0,0]), np.squeeze(X[-1,0,0]) , np.squeeze(Z[0,0,0]), np.squeeze(Z[0,0,-1])], cmap='jet', origin='lower', interpolation='spline16')


#plt.pcolormesh(np.squeeze(X[:,:,150]), np.squeeze(Y[:,:,150]), np.squeeze(u[:,:,150]), cmap='jet')
plt.imshow(np.squeeze(U[:,:,150]).T,extent=[np.squeeze(X[0,0,150]), np.squeeze(X[-1,0,150]) , np.squeeze(Y[0,0,150]), np.squeeze(Y[0,-1,150])], cmap='jet', origin='lower', interpolation='spline16')

plt.pcolormesh(np.squeeze(Z[150,:,:]), np.squeeze(Y[150,:,:]), np.squeeze(u[150,:,:]), cmap='jet')



#%%




#%%                   Statitics
U=Read_stats('u',part=1)
UU=Read_stats('uu',part=1)
UU=UU-U*U

UUU=Read_stats('uuu',part=1)
UUU=UUU-3*UU*U-U*U*U

V=Read_stats('v',part=1)
VV=Read_stats('vv',part=1)
VV=VV-V*V

W=Read_stats('w',part=1)
WW=Read_stats('ww',part=1)
WW=WW-W*W


UV=Read_stats('uv',part=1)
UV=UV-U*V


k=0.5*(UU+VV+WW)

#%%


f, ax = plt.subplots(1, 1, sharey=True,figsize=(8, 5),  facecolor='w', edgecolor='k') 

Var=k

# plot  <var>_x,z,t vs y
_, Varp=MeanXZ(Var,utau,normi=2)
ax.semilogx(yp,Varp,'r-',linewidth=3, label=r'$\left< k\right>^+_{x,z,t}$')

# plot  <var(xi,zi)>_t  vs y with xi and zi random
tmp= np.squeeze(Var[np.random.randint(0,372),:,np.random.randint(0,360)])
tmp=wallunit(tmp,utau,2)   
ax.semilogx(yp,tmp,'b--',linewidth=3, label=r'$\left<k \right>^+_{t}$')
ax.set_xlabel(r'$y^+$', fontsize=20)
ax.set_ylabel(r'$\left<k \right>^+$', fontsize=20)
ax.set_ylim(ymin=0)
ax.set_xlim(yp[0],yp[-1])
ax.legend()
plt.rcParams.update({'font.size': 18})

#%%

f, ax = plt.subplots(1, 1, sharey=True,figsize=(8, 5),  facecolor='w', edgecolor='k') 
Var=(UU+VV+WW)/2

_, Varp=MeanZ(Var,utau,normi=2)
#ax.semilogx(yp,up,'r-',linewidth=3, label=r'$\left< k\right>^+_{x,z,t}$')

xx,yy=np.meshgrid(np.squeeze(X[:,0,0]), yp)

plot_2d(Varp,xx,yy,ylog=True,ylab='y^+',xlab='x^+')

#plt.savefig('/home/cfd/Desktop/Kara/channel180/plots/k.pdf',  bbox_inches='tight')
