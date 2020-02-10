import numpy as np
import h5py


rep=''

print('Load data (r,u,v,w)-Jpdf')

rhom = 0.014
tauw = 5.92464e-05
utau=(tauw/rhom)**0.5
mu = 5.05967e-06
Retau=(rhom/mu)*utau



hf = h5py.File(rep+'MPdf_Full_nbins60.h5', 'r')
Pdf = hf.get('pdf')
yp = hf.get('yp')
#utau = hf.get('utau')
g1 = hf.get('limites')
B0 = g1.get('B0')
B1 = g1.get('B1')
B2 = g1.get('B2')
B3 = g1.get('B3')

Pdf=np.array(Pdf)
yp=np.array(yp)
B0=np.array(B0)
B1=np.array(B1)
B2=np.array(B2)
B3=np.array(B3)
hf.close()
print(f'Size of the M-pdf : {Pdf.shape}')




def Xmeans(n0,n1,n2,n3,yi,Pdf=Pdf,B0=B0,B1=B1,B2=B2,B3=B3):
    b0=B0[yi]
    b1=B1[yi]
    b2=B2[yi]
    b3=B3[yi]
    
    x1 = (b1[1:]+b1[:-1])*0.5
    x2 = (b2[1:]+b2[:-1])*0.5
    x3 = (b3[1:]+b3[:-1])*0.5
    x0 = (b0[1:]+b0[:-1])*0.5
    
    Mpdf=np.squeeze(Pdf[yi,:,:,:,:])
    A=np.trapz(np.trapz(np.trapz(np.trapz(Mpdf,x0,axis=0),x1,axis=0),x2,axis=0),x3,axis=0)
    Mpdf=Mpdf/A
    
    
    X0,X1,X2,X3=np.meshgrid(x0,x1,x2,x3, indexing='ij')
    
    t=np.trapz(np.trapz(np.trapz(np.trapz((X0**n0)*(X1**n1)*(X2**n2)*(X3**n3)*Mpdf,x0,axis=0),x1,axis=0),x2,axis=0),x3,axis=0)
    return t

def Xproducts(n0,n1,n2,n3,yi,Favre=False,Pdf=Pdf,B0=B0,B1=B1,B2=B2,B3=B3):
    
    b0=B0[yi]
    b1=B1[yi]
    b2=B2[yi]
    b3=B3[yi]
    
    x1 = (b1[1:]+b1[:-1])*0.5
    x2 = (b2[1:]+b2[:-1])*0.5
    x3 = (b3[1:]+b3[:-1])*0.5
    x0 = (b0[1:]+b0[:-1])*0.5
    
    Mpdf=np.squeeze(Pdf[yi,:,:,:,:])
    A=np.trapz(np.trapz(np.trapz(np.trapz(Mpdf,x0,axis=0),x1,axis=0),x2,axis=0),x3,axis=0)
    Mpdf=Mpdf/A
    
    X0,X1,X2,X3=np.meshgrid(x0,x1,x2,x3, indexing='ij')
    
    
    Rm=0
    Um=0
    Vm=0
    Wm=0
    if Favre==True:
        
       if n0>0:
            Rm=Xmeans(1,0,0,0,yi)
       if n1>0:    
            Um=Xmeans(1,1,0,0,yi)
       if n2>0:    
            Vm=Xmeans(1,0,1,0,yi)
       if n3>0:    
            Wm=Xmeans(1,0,0,1,yi)
        
        
       t=np.trapz(np.trapz(np.trapz(np.trapz(((X0*X1-Um)**n1)*((X0*X2-Vm)**n2)*((X0*X3-Wm)**n3)*Mpdf,x0,axis=0),x1,axis=0),x2,axis=0),x3,axis=0)
       t=t/(Rm**(n1+n2+n3))
    else:    
        n0=0
        if n0>0:
            Rm=Xmeans(1,0,0,0,yi)
        if n1>0:    
            Um=Xmeans(0,1,0,0,yi)
        if n2>0:    
            Vm=Xmeans(0,0,1,0,yi)
        if n3>0:    
            Wm=Xmeans(0,0,0,1,yi)
    
        t=np.trapz(np.trapz(np.trapz(np.trapz(((X0-Rm)**n0)*((X1-Um)**n1)*((X2-Vm)**n2)*((X3-Wm)**n3)*Mpdf,x0,axis=0),x1,axis=0),x2,axis=0),x3,axis=0)
    
    
    
    
    return t    
