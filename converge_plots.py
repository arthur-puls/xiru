# This module has some plotting-related stuff
import matplotlib.pyplot as plt
import numpy as np
from parsemoog import diff_ovec

def make_dif_feplot(fe1ref, fe2ref, mh, refmh, stdout='a.out'):
    "Plots differential Boltzmann diagrams and saves to disc"
    epsl, gpar, mpar, rwsl, dfe1, dfe2, ep1, rw1, rw2, epint, rwint = diff_ovec(stdout, fe1ref, fe2ref, mh, refmh, mkplot=True)
    
    eplim = [0.0, 5.5]
    rwlim = [-5.90, -4.44]
    ylim = [np.concatenate((dfe1, dfe2)).min() - 0.05, np.concatenate((dfe1, dfe2)).max() + 0.05]
    ymid = np.around(0.50*(ylim[1] - ylim[0]) + ylim[0], decimals=2)
    ytickdown = np.around(0.20*(ylim[1] - ylim[0]) + ylim[0], decimals=2)
    ytickup = np.around(0.80*(ylim[1] - ylim[0]) + ylim[0], decimals=2)
    yticks = (ytickdown, ymid, ytickup)

    fe1h = np.mean(dfe1)

    epx = np.linspace(eplim[0], eplim[1], 100)
    rwx = np.linspace(rwlim[0], rwlim[1], 100)
    
    xoffset = 0.02
    
    figfe, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(6.4,4.7))
    
    ax1.hlines(fe1h+np.std(dfe1, ddof=1), eplim[0], eplim[1], color='g', linestyles='-.', linewidth=0.8)
    ax1.hlines(fe1h-np.std(dfe1, ddof=1), eplim[0], eplim[1], color='g', linestyles='-.', linewidth=0.8)
    ax2.hlines(fe1h+np.std(dfe1, ddof=1), rwlim[0], rwlim[1], color='g', linestyles='-.', linewidth=0.8)
    ax2.hlines(fe1h-np.std(dfe1, ddof=1), rwlim[0], rwlim[1], color='g', linestyles='-.', linewidth=0.8)
    ax3.hlines(fe1h+np.std(dfe1, ddof=1), rwlim[0], rwlim[1], color='g', linestyles='-.', linewidth=0.8)
    ax3.hlines(fe1h-np.std(dfe1, ddof=1), rwlim[0], rwlim[1], color='g', linestyles='-.', linewidth=0.8)
    
    ax1.plot(epx, epx*epsl + epint, 'k--')
    ax2.plot(rwx, rwx*rwsl + rwint, 'k-.')
    ax3.plot(rwx, fe1h*np.ones(100), 'k:')
    
    ax1.scatter(ep1, dfe1, marker='o', c='r')
    ax2.scatter(rw1, dfe1, marker='o', c='r')
    ax3.scatter(rw1, dfe1, marker='o', c='r')
    ax3.scatter(rw2, dfe2, marker='s', c='b')
    
    ax1.set_xlim(eplim)
    ax2.set_xlim(rwlim)
    ax3.set_xlim(rwlim)
    ax1.set_ylim(ylim)
    ax2.set_ylim(ylim)
    ax3.set_ylim(ylim)
    
    ax1.set_xlabel('E.P. (eV)')
    ax3.set_xlabel('R.W.')
    
    figfe.text(xoffset, 0.5, '$\delta$A(Fe)', rotation=90, va='center', ha='center')
    
    ax1.set_yticks(yticks)
    ax2.set_yticks(yticks)
    ax3.set_yticks(yticks)
    
    ax1.tick_params(bottom=True,
    	top=False,
    	direction='in',
    	labelbottom=True,
    	labeltop=False)

    ax2.tick_params(bottom=True,
    	top=False,
    	direction='in',
    	labelbottom=False,
    	labeltop=False)

    ax3.tick_params(bottom=True,
    	top=False,
    	direction='in',
    	labelbottom=True,
    	labeltop=False)
    
    plt.savefig('feplot.png', dpi=200)