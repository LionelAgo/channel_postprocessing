import matplotlib.pyplot as plt
import numpy as np


def plot_2d(tp,x1,x0,a=1,bd=11,cmap='jet', xlog=False,ylog=False, cmin=None,cmax=None,xmin=None,xmax=None,ymin=None,ymax=None,ylab='',xlab='' ):
    
    
    if cmin is None:
        cmin=tp.min()
    if cmax is None:    
        cmax=tp.max()
    
    
    bd=np.linspace(a*cmin, a*cmax, np.int(bd))
    
    if xmax is None:
        xmax=x1[-1]
    if xmin is None:
        xmin=x1[0]
    if ymax is None:
        ymax=x0[-1]
    if ymin is None:  
        ymin=x0[0]
    fig, ax = plt.subplots(1, 1, sharey=True,figsize=(12, 8), dpi=80, facecolor='w', edgecolor='k')
    ax.contourf(x1,x0,tp,cmap=cmap, origin='lower',extent=(xmin,xmax,ymin,ymax), levels=bd) 
    ax.contour(x1,x0,tp, origin='lower',extent=(xmin,xmax,ymin,ymax),colors='black', levels=bd,linewidth=1) 
    
    if ylog==True:
        ax.set_yscale('log')
    if xlog==True:
        ax.set_xscale('log')
    xl=f'${xlab}$'
    ax.set_xlabel(xl, fontsize=20)
    yl=f'${ylab}$'
    ax.set_ylabel(yl, fontsize=20)
  
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    #return fig    
    #plt.show()