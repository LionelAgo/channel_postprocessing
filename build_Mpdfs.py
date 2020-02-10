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



B0=[]
B1=[]
B2=[]
B3=[]

nbins = 60
Pdf=np.zeros((int(len(Ly)/2+1),nbins,nbins+2,nbins+4,nbins+6))

m=0
for num in np.arange(2050,6300,50):
    m+=1
    suff=f'{num}.0000.pyfrs'
    print(suff)
   
    rho=Read_solution(num,0)
    rt=np.transpose(rho,(1,0,2))
    rt=rt.reshape(rt.shape[0],-1)

    u=Read_solution(num,1)
    u=u/rho
    ut=np.transpose(u,(1,0,2))
    ut=ut.reshape(ut.shape[0],-1)
    
    v=Read_solution(num,2)
    v=v/rho
    vt=np.transpose(v,(1,0,2))
    vt=vt.reshape(vt.shape[0],-1)

    w=Read_solution(num,3)
    w=w/rho
    wt=np.transpose(w,(1,0,2))
    wt=wt.reshape(wt.shape[0],-1)
    
    

    
    for yi in range(int(len(Ly)/2+1)):

        r_=rt[(yi,113-yi),:]
        r_=np.concatenate((r_[0,:],r_[1,:]),0)
        u_=ut[(yi,113-yi),:]
        u_=np.concatenate((u_[0,:],u_[1,:]),0)
    
        v_=vt[(yi,113-yi),:]
        v_=np.concatenate((v_[0,:],-1*v_[1,:]),0)

        w_=wt[(yi,113-yi),:]
        w_=np.concatenate((w_[0,:],w_[1,:]),0)
        

        
        if m==1:
            Mpdf, [b0, b1, b2,b3]=np.histogramdd((r_,u_,v_,w_),bins=(nbins,nbins+2,nbins+4,nbins+6),density=True)
            B0.append(b0)
            B1.append(b1)
            B2.append(b2)
            B3.append(b3)
        else:
            b0=B0[yi]
            b1=B1[yi]
            b2=B2[yi]
            b3=B3[yi]
            Mpdf,[_,_,_,_]=np.histogramdd((r_,u_,v_,w_),bins=(b0,b1,b2,b3),density=True)
            
        Pdf[yi,:,:,:,:]+=Mpdf
    
    #x1 = (b1[1:]+b1[:-1])*0.5
    #x2 = (b2[1:]+b2[:-1])*0.5
    #x3 = (b3[1:]+b3[:-1])*0.5
    #x0 = (b0[1:]+b0[:-1])*0.5
    

    
hf = h5py.File('MPdf_Full_nbins60.h5', 'w')  
hf.create_dataset('pdf', data=Pdf, compression="gzip", compression_opts=9)
hf.create_dataset('yp', data=yp)
#hf.create_dataset('utau', data=utau)
g1 = hf.create_group('limites')    
g1.create_dataset('B0',data=B0, compression="gzip", compression_opts=9)
g1.create_dataset('B1',data=B1, compression="gzip", compression_opts=9)
g1.create_dataset('B2',data=B2, compression="gzip", compression_opts=9)
g1.create_dataset('B3',data=B3, compression="gzip", compression_opts=9)

hf.close()

